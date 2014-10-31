(in-package :elprep)

;; filters

(defun replace-reference-sequence-dictionary (reference-sequence-dictionary)
  "A filter for replacing the reference sequence dictionary in a sam-header."
  (lambda (header)
    ;; changing the reference-sequence-dictionary may destroy the sorting order
    (when-let (sorting-order (getf (sam-header-hd header) :so))
      (when (string= sorting-order "coordinate")
        (unless (loop with previous-pos = -1
                      with old-dictionary = (sam-header-sq header)
                      for entry in reference-sequence-dictionary
                      for sname = (getf entry :sn)
                      for pos = (position sname old-dictionary :test 'string= :key (lambda (entry) (getf entry :sn)))
                      do (cond ((null pos))
                               ((> pos previous-pos)
                                (setq previous-pos pos))
                               (t (return nil)))
                      finally (return t))
          ; "unknown", not "unsorted", because the offending snames may not occur in the reads
          (sam-header-ensure-hd header :so "unknown"))))
    (let ((reference-sequence-table
           (make-hash-table :test 'equal :single-thread t)))
      (declare (hash-table reference-sequence-table))
      (loop for entry in reference-sequence-dictionary
            do (setf (gethash (getf entry :sn) reference-sequence-table) t))
      (setf (sam-header-sq header) reference-sequence-dictionary)
      (lambda ()
        (lambda (alignment)
          (declare (sam-alignment alignment) #.*optimization*)
          (gethash (sam-alignment-rname alignment) reference-sequence-table))))))

(defun replace-reference-sequence-dictionary-from-sam-file (sam-file)
  "A filter for replacing the reference sequence dictionary in a sam-header with one parsed from the given SAM/DICT file."
  (with-open-stream (in (open-sam sam-file :direction :input :header-only t))
    (let* ((header (parse-sam-header in))
           (dictionary (sam-header-sq header))) 
      (replace-reference-sequence-dictionary dictionary))))

(defun filter-unmapped-reads (header)
  "A filter for removing unmapped sam-alignment instances, based on FLAG."
  (declare (ignore header))
  (lambda ()
    (lambda (alignment)
      (declare (sam-alignment alignment) #.*optimization*)
      (= 0 (logand (sam-alignment-flag alignment) +unmapped+)))))

(defun filter-unmapped-reads-strict (header)
  "A filter for removing unmapped sam-alignment instances, based on FLAG, or POS=0, or RNAME=*."
  (declare (ignore header))
  (lambda ()
    (lambda (alignment)
      (declare (sam-alignment alignment) #.*optimization*)
      (or (= 0 (logand (sam-alignment-flag alignment) +unmapped+))
          (= (sam-alignment-pos alignment) 0)
          (string= (sam-alignment-rname alignment) "*")))))

(defun filter-duplicate-reads (header)
  "A filter for removing duplicate sam-alignment instances, based on FLAG."
  (declare (ignore header))
  (lambda ()
    (lambda (alignment)
      (declare (sam-alignment alignment) #.*optimization*)
      (= 0 (logand (sam-alignment-flag alignment) +duplicate+)))))

(defun filter-optional-reads (header)
  "A filter for removing sam-alignment instances that represent optional information in elPrep."
  (when (sam-header-user-tag header :|@sr|)
    (remf (sam-header-user-tags header) :|@sr|)
    (lambda ()
      (lambda (alignment)
        (declare (sam-alignment alignment) #.*optimization*)
        (not (sam-alignment-tag alignment :|sr|))))))
           
(defun add-or-replace-read-group (read-group)
  "A filter for adding or replacing the read group both in sam-header and each sam-alignment."
  (lambda (header)
    (setf (sam-header-rg header) (list read-group))
    (let ((id (getf read-group :ID))) 
      (lambda ()
        (lambda (alignment)
          (declare (sam-alignment alignment) #.*optimization*)
          (setf (sam-alignment-tag alignment :rg) id)
          t)))))

(defun parse-read-group-from-string (string)
  "Parse an @RG line in a string.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3,
   except that entries are separated by white space instead of tabulators."
  (declare (string string) #.*fixnum-optimization*)
  (let ((record '()) 
        (string-buffer (make-array 80 
                                   :element-type 'base-char 
                                   :adjustable t 
                                   :fill-pointer 0 
                                   :single-thread t)))
    (with-input-from-string (stream string)
      (labels ((read-char-or-error ()
                 (or (read-char stream nil nil)
                     (error "Unexpected end of read group string ~S." string)))
               (scan-tag ()
                 (unless (peek-char t stream nil nil)
                   (assert (getf record :ID))
                   (return-from parse-read-group-from-string record))
                 (let ((tag (make-array 2 :element-type 'base-char)))
                   (setf (sbchar tag 0) (char-upcase (read-char-or-error))
                         (sbchar tag 1) (char-upcase (read-char-or-error)))
                   (when (char/= (read-char-or-error) #\:)
                     (error "Expected tag separator for tag ~A in read group string ~S." tag string))
                   tag))
               (scan-string (tag)
                 (setf (fill-pointer string-buffer) 0)
                 (when (char= (peek-char nil stream nil #\Space) #\Space)
                   (error "Expected tag value for tag ~S in read group string ~S." tag string))
                 (loop do (vector-push-extend (read-char-or-error) string-buffer)
                       until (char= (peek-char nil stream nil #\Space) #\Space))
                 (subseq string-buffer 0 (fill-pointer string-buffer))))
        (loop for tag = (scan-tag)
              for val = (scan-string tag) do
              (nconcf record (string-case (tag :default (if (sam-header-user-tag-p tag)
                                                            (list (unique (intern-key tag) record) val)
                                                          (error "Unknown tag ~A in read group string ~S." tag string)))
                               ("ID" (list (unique :ID record) val))
                               ("CN" (list (unique :CN record) val))
                               ("DS" (list (unique :DS record) val))
                               ("DT" (list (unique :DT record) (parse-date-time val)))
                               ("FO" (list (unique :FO record) val))
                               ("KS" (list (unique :KS record) val))
                               ("LB" (list (unique :LB record) val))
                               ("PG" (list (unique :PG record) val))
                               ("PI" (list (unique :PI record) (parse-integer val)))
                               ("PL" (list (unique :PL record) val))
                               ("PU" (list (unique :PU record) val))
                               ("SM" (list (unique :SM record) val)))))))))

(defun add-pg-line (id &rest args &key pn cl ds vn)
  "A filter for adding a @PG tag to a sam-header, and ensuring that it is the first one in the chain."
  (declare (dynamic-extent args) (ignore pn cl ds vn))
  (let ((line (copy-list args)))
    (lambda (header)
      (let ((pgs (sam-header-pg header))
            (new-line line))
      ; ensure id is unique
        (loop while (find id pgs :key (lambda (pg) (getf pg :id)) :test 'equal)
              do (setq id (format nil "~A ~4,'0X" id (random #x10000))))
      ; determine PP tag
        (when-let (next-id (loop for pg in pgs
                                 for id = (getf pg :id)
                                 unless (find id pgs :key (lambda (pg) (getf pg :pp)) :test 'equal)
                                 return id))
          (setq new-line (list* :pp next-id new-line)))
      ; add ID tag
        (setq new-line (list* :id id new-line))
      ; add @PG line
        (push new-line (sam-header-pg header))
      ; no next-level filters necessary
        nil))))

(defun rename-chromosomes (header)
  "A filter for prepending \"chr\" to the reference sequence names in a sam-header, and in RNAME and RNEXT in each sam-alignment."
  (loop for plist in (sam-header-sq header)
        for sn = (getf plist :SN)
        when sn do (setf (getf plist :SN) (string-append "chr" sn)))
  (lambda ()
    (lambda (alignment)
      (declare (sam-alignment alignment) #.*optimization*)
      (let ((rname (sam-alignment-rname alignment)))
        (cond ((string= rname "="))
              ((string= rname "*"))
              (t (setf (sam-alignment-rname alignment) (string-append "chr" rname)))))
      (let ((rnext (sam-alignment-rnext alignment)))
        (cond ((string= rnext "="))
              ((string= rnext "*"))
              (t (setf (sam-alignment-rnext alignment) (string-append "chr" rnext)))))
      t)))

(defun add-refid (header)
  "A filter for adding the refid (index in the reference sequence dictionary) to sam-alignment instances."
  (let ((reference-sequence-table (make-hash-table :test 'equal :single-thread t)))
    (loop with ctr = -1
          for sn-form in (sam-header-sq header)
          do (setf (gethash (getf sn-form :SN) reference-sequence-table) (incf ctr)))
      (lambda ()
        (lambda (alignment)
          ; fill in refid
          (setf (sam-alignment-refid alignment)
                (gethash (sam-alignment-rname alignment) reference-sequence-table -1))))))