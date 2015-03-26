(in-package :elprep)
(in-simple-base-string-syntax)

(defun explain-flag (flag)
  "Return a symbolic representation of the FLAG entry in a SAM file alignment line.
   This is primarily for debugging purposes.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4."
  (let ((result '()))
    (macrolet ((test (&rest bits)
                 `(progn ,@(loop for bit in bits
                                 for bitn = (symbol-name bit)
                                 for bitk = (intern (subseq bitn 1 (1- (length bitn))) :keyword)
                                 collect `(when (/= (logand flag ,bit) 0)
                                            (push ,bitk result))))))
      (test +supplementary+
            +duplicate+
            +qc-failed+
            +secondary+
            +last+
            +first+
            +next-reversed+
            +reversed+
            +next-unmapped+
            +unmapped+
            +proper+
            +multiple+))
    result))
 
(defun sam-alignment-differ (aln1 aln2)
  "Return false if the two sam-alignments have the same mandatory fields.
   Otherwise, return a symbol indicating which field differs."
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  ; check that all mandatory fields are =
  (or (when (string/= (the base-string (sam-alignment-qname aln1)) (the base-string (sam-alignment-qname aln2))) 'qname)
      (when (/= (the fixnum (sam-alignment-flag aln1)) (the fixnum (sam-alignment-flag aln2))) 'flag)
      (when (string/= (the base-string (sam-alignment-rname aln1)) (the base-string (sam-alignment-rname aln2))) 'rname)
      (when (/= (the int32 (sam-alignment-pos aln1)) (the int32 (sam-alignment-pos aln2))) 'pos)
      (when (/= (sam-alignment-mapq aln1) (sam-alignment-mapq aln2)) 'mapq)
      (when (string/= (the base-string (sam-alignment-cigar aln1)) (the base-string (sam-alignment-cigar aln2))) 'cigar)
      (when (string/= (the base-string (sam-alignment-rnext aln1)) (the base-string (sam-alignment-rnext aln2))) 'rnext)
      (when (string/= (the base-string (sam-alignment-qual aln1)) (the base-string (sam-alignment-qual aln2))) 'qual)))

(defun sam-alignment-same (aln1 aln2)
  "Return true if the two sam-alignments have the same mandatory fields, false otherwise."
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  (and (string= (the base-string (sam-alignment-qname aln1)) (the base-string (sam-alignment-qname aln2)))
       (= (the fixnum (sam-alignment-flag aln1)) (the fixnum (sam-alignment-flag aln2)))
       (string= (the base-string (sam-alignment-rname aln1)) (the base-string (sam-alignment-rname aln2)))
       (= (the int32 (sam-alignment-pos aln1)) (the int32 (sam-alignment-pos aln2)))
       (= (sam-alignment-mapq aln1) (sam-alignment-mapq aln2))
       (string= (the base-string (sam-alignment-cigar aln1)) (the base-string (sam-alignment-cigar aln2)))
       (string= (the base-string (sam-alignment-rnext aln1)) (the base-string (sam-alignment-rnext aln2)))
       (string= (the base-string (sam-alignment-qual aln1)) (the base-string (sam-alignment-qual aln2)))))

(defun real-diffs (alns1 alns2)
  "Return a list of sam-alignments in alns1 for which no alignments in alns2 exist that have the same mandatory fields."
  (loop for aln1 in alns1
        unless (find aln1 alns2 :test #'sam-alignment-same)
        collect aln1))

(defun compare-sams (sam1-file sam2-file)
  "Parse both SAM files, then compare the mandatory fields of all alignments one by one."
  (let ((sam1 (make-sam))
        (sam2 (make-sam))
        (working-directory (get-working-directory)))
    (run-pipeline (merge-pathnames sam1-file working-directory) sam1)
    (run-pipeline (merge-pathnames sam2-file working-directory) sam2)
    ; sort the sams by qname
    (setf (sam-alignments sam1) (stable-sort (sam-alignments sam1) #'string< :key #'sam-alignment-qname))
    (setf (sam-alignments sam2) (stable-sort (sam-alignments sam2) #'string< :key #'sam-alignment-qname))
    (format t "sam1:~s alns sam2:~s alns ~%" (length (sam-alignments sam1)) (length (sam-alignments sam2)))
    (let ((differences1 nil)
          (differences2 nil))
      (loop for aln1 in (sam-alignments sam1) ; filter diffs
            for aln2 in (sam-alignments sam2)
            do (let ((d (sam-alignment-differ aln1 aln2))) 
                 (when d 
                   (push aln1 differences1) 
                   (push aln2 differences2))))
      (real-diffs differences1 differences2)))) ; sort slightly different order in elprep so get out real diffs

(defun verify-order-kept (sam-file)
  "Assume the SAM file is sorted by coordinate order. Verify that this is stil the case."
  (format t "verifying order kept ~%")
  (let ((sam (make-sam))
        (working-directory (get-working-directory)))
    (run-pipeline (merge-pathnames sam-file working-directory) sam)
    (let ((pos (sam-alignment-pos (first (sam-alignments sam))))
          (rname (sam-alignment-rname (first (sam-alignments sam))))
          (ctr 1))
      (loop for aln in (rest (sam-alignments sam))
            do (let ((new-pos (sam-alignment-pos aln))
                     (new-rname (sam-alignment-rname aln)))
                 (cond ((and (< new-pos pos) (string= rname new-rname ))
                        (format t "Not sorted: previous pos: ~s,~s current pos: ~s,~s. ~s reads were in the right order. ~%" rname pos new-rname new-pos ctr) 
                        (return nil))
                       (t 
                        (incf ctr)
                        (setf rname new-rname)
                        (setf pos new-pos))))
            finally (return t)))))

(defun count-duplicates (sam-file)
  "Return the number of alignments in the SAM file that are marked as duplicates."
  (let ((sam (make-sam)))  
    (run-pipeline (merge-pathnames sam-file (get-working-directory)) sam)
    (loop for aln in (sam-alignments sam)
          count (sam-alignment-duplicate-p aln))))

; code for splitting up sam files into chromosomes

(define-symbol-macro optional-data-tag "sr:i:1")

(defun split-file-per-chromosome (input output-path output-prefix output-extension)
  "A function for splitting a sam file into : a file containing all unmapped reads, a file containing all pairs where reads map to different chromosomes, a file per chromosome containing all pairs where the reads map to that chromosome. There are no requirements on the input file for splitting."
  (let ((files (directory input)))
    (multiple-value-bind
        (first-in first-program)
        (open-sam (first files) :direction :input)
      (let ((header (parse-sam-header first-in))
            (header-file (merge-pathnames output-path (make-pathname :name (format nil "~a-header" output-prefix) :type "dict")))
            (splits-path (merge-pathnames (make-pathname :directory '(:relative "splits")) output-path))
            (chroms-encountered (make-single-thread-hash-table :test #'buffer= :hash-function #'buffer-hash))
            (chroms-encountered-minmax (make-single-thread-hash-table :test #'buffer= :hash-function #'buffer-hash))
            (buf-unmapped (make-buffer "*")))
        ; create a directory for split files
        (ensure-directories-exist splits-path)
        ; tag the header as one created with elPrep split
        (setf (sam-header-user-tag header :|@sr|) (list "co:This file was created using elprep split."))
        ; fill in a file for unmapped reads
        (setf (gethash buf-unmapped chroms-encountered)
              (multiple-value-bind
                  (file program) ; create a new headerless file
                  (open-sam (merge-pathnames splits-path (make-pathname :name (format nil "~a-unmapped" output-prefix) :type output-extension)) :direction :output)
                (format-sam-header file header)
                (cons file program)))
        ; fill in a minmax for unmapped reads
        (setf (gethash buf-unmapped chroms-encountered-minmax) (cons nil (cons 0 0)))
        (loop for sn-form in (sam-header-sq header)
              do (assert (not (null sn-form)))
              (let* ((chrom (getf sn-form :SN))
                     (buf-chrom (make-buffer chrom)))
                (setf (gethash buf-chrom chroms-encountered)
                      (multiple-value-bind
                          (file program)
                          (open-sam (merge-pathnames splits-path (make-pathname :name (format nil "~a-~a" output-prefix chrom) :type output-extension)) :direction :output)
                        (format-sam-header file header)
                        (cons file program)))
                   ; minmax table
                (setf (gethash buf-chrom chroms-encountered-minmax) (cons sn-form (cons nil nil)))))
        (unwind-protect
            (with-open-sam (spread-reads-stream (merge-pathnames output-path (make-pathname :name (format nil "~a-spread" output-prefix) :type output-extension)) :direction :output)
              (format-sam-header spread-reads-stream header)
              (let ((buf-= (make-buffer "=")))
                (let ((rname (make-buffer))
                      (rnext (make-buffer))
                      (aln-string (make-buffer))
                      (pos (make-buffer)))
                  (flet ((process-file (in)
                           (loop do (read-line-into-buffer in aln-string)
                                 until (buffer-emptyp aln-string)
                                 do (progn (buffer-partition aln-string #\Tab 2 rname 3 pos 6 rnext)
                                      (let ((file (car (gethash rname chroms-encountered))))
                                        (cond ((or (buffer= buf-= rnext) (buffer= buf-unmapped rname) (buffer= rname rnext))
                                               (write-buffer aln-string file)
                                               (write-newline file)
                                                 ; fill in minmax info
                                               (let ((ppos (buffer-parse-integer pos))
                                                     (minmax (cdr (gethash rname chroms-encountered-minmax))))
                                                 (cond ((not (car minmax)) (setf (car minmax) ppos) (setf (cdr minmax) ppos))
                                                       ((< ppos (car minmax)) (setf (car minmax) ppos))
                                                       ((> ppos (cdr minmax)) (setf (cdr minmax) ppos)))))
                                              (t ; the read is part of a pair mapping to two different chromosomes
                                               (write-buffer aln-string spread-reads-stream)
                                               (write-newline spread-reads-stream)
                                               ; duplicate the info in the chromosome file so it can be used; mark the read as duplicate info
                                               (write-buffer aln-string file)
                                               (write-tab file)
                                               (writestr file optional-data-tag)
                                               (write-newline file))))))
                           (reinitialize-buffer rname)
                           (reinitialize-buffer rnext)
                           (reinitialize-buffer aln-string)
                           (reinitialize-buffer pos)))
                    (process-file first-in)
                    (close-sam first-in first-program)
                    (loop for in-file in (rest files)
                          do (with-open-sam (in in-file :direction :input)
                               (skip-sam-header in)
                               (process-file in)))))))
          (loop for (file . program) being each hash-value of chroms-encountered
                do (close-sam file program))
          ; write out min/max positions in new header
          (loop for (sn-form . minmax) being each hash-value of chroms-encountered-minmax
                do (when sn-form
                     (setf (cdr (last sn-form)) (list :|mn| (format nil "~D" (car minmax)) :|mx| (format nil "~D" (cdr minmax))))))
          ; write the header file for the split files
          (with-open-sam (header-out header-file :direction :output)
            (format-sam-header header-out header)))))))

(defun merge-sorted-files-split-per-chromosome (input-path output input-prefix input-extension header)
  "A function for merging files that were split with elPrep and sorted in coordinate order."
  ; Extract the header to identify the files names. 
  ; Assume that all files are sorted per cooordinate order, i.e. first sorted on refid entry according to sequence dictionary, then sorted on position entry. 
  ; There is a file per chromosome in the sequence dictionary. These contain all reads that map to that chromosome. 
  ; On top of that, there is a file that contains the unmapped (or *) reads and a file that contains the reads that map to different chromosomes.
  ; Merge these files in the order of the sequence dictionary. Put the unmapped reads as the last entries.
  ; When merging a particular chromosome file into the merged file, make sure that reads that map to different chromosomes are merged in correctly.
  ; So while mergin a particular chromosome file, pop and compare against reads in the file for reads that map to different chromosomes until the next chromosome
  ; is encountered on the refid position.
  ; when a file is empty, close it and remove it from the list of files to merge
  ; loop for identifying and opening the files to merge
  (with-open-sam (spread-reads-file (merge-pathnames input-path (make-pathname :name (format nil "~a-spread" input-prefix) :type input-extension)) :direction :input)
    (skip-sam-header spread-reads-file)
    ; merge loop
    (with-open-sam (out output :direction :output)
      (format-sam-header out header)
      (let ((spread-read (make-buffer)) ; for storing entries from the spread-read file
            (spread-read-refid (make-buffer))
            (spread-read-pos (make-buffer))
            (chromosome-read (make-buffer)) ; for storing reads from the chromsome file we are currently merging
            (chromosome-read-refid (make-buffer))
            (chromosome-read-pos (make-buffer))
            (common-read-refid (make-buffer)))
        ; then merge the rest of the files
        (loop for sn-form in (sam-header-sq header)
              for chrom = (getf sn-form :SN)
              for file-name = (merge-pathnames input-path (make-pathname :name (format nil "~a-~a" input-prefix chrom) :type input-extension))
              when (probe-file file-name) do
              (reinitialize-buffer common-read-refid)
              (buffer-extend common-read-refid chrom)
              (block chromosome-loop
                (with-open-sam (file file-name :direction :input)
                  (skip-sam-header file)
                  (read-line-into-buffer file chromosome-read)
                  (when (buffer-emptyp chromosome-read) (return-from chromosome-loop)) 
                  (buffer-partition chromosome-read #\Tab 2 chromosome-read-refid 3 chromosome-read-pos)
                  (assert (buffer= chromosome-read-refid common-read-refid))
                  (when (buffer-emptyp spread-read) ; if the buffer is not empty, the current entry is potentially an entry for this file and it should not be overwritten
                    (read-line-into-buffer spread-reads-file spread-read)
                    (buffer-partition spread-read #\Tab 2 spread-read-refid 3 spread-read-pos))
                  (unless (buffer-emptyp spread-read)
                    (when (buffer= spread-read-refid chromosome-read-refid)
                      (let ((pos1 (buffer-parse-integer spread-read-pos))
                            (pos2 (buffer-parse-integer chromosome-read-pos)))
                        (loop do (cond ((< pos1 pos2)
                                        (write-buffer spread-read out)
                                        (write-newline out)
                                        (read-line-into-buffer spread-reads-file spread-read)
                                        (cond ((buffer-emptyp spread-read)
                                               (loop-finish))
                                              (t (buffer-partition spread-read #\Tab 2 spread-read-refid 3 spread-read-pos)
                                                 (setq pos1 (buffer-parse-integer spread-read-pos)))))
                                       (t (write-buffer chromosome-read out)
                                          (write-newline out)
                                          (read-line-into-buffer file chromosome-read)
                                          (cond ((buffer-emptyp chromosome-read)
                                                 (loop-finish))
                                                (t (buffer-partition chromosome-read #\Tab 3 chromosome-read-pos)
                                                   (setq pos2 (buffer-parse-integer chromosome-read-pos))))))
                              while (buffer= chromosome-read-refid spread-read-refid)))))
                  ; copy remaining reads in the file, if any
                  (when (not (buffer-emptyp chromosome-read))
                    (write-buffer chromosome-read out)
                    (write-newline out))
                  (copy-stream file out)))
              ; copy remaining reads in the spread file, if any, that are one the same chromosome as the file was
              (when (not (buffer-emptyp spread-read))
                (loop while (buffer= spread-read-refid common-read-refid)
                      do
                      (write-buffer spread-read out)
                      (write-newline out)
                      (read-line-into-buffer spread-reads-file spread-read)
                      (when (buffer-emptyp spread-read) (return))
                      (buffer-partition spread-read #\Tab 2 spread-read-refid)
                      finally
                      (buffer-partition spread-read #\Tab 3 spread-read-pos))))
        ; merge the remaining reads in the spread-reads file
        (when (not (buffer-emptyp spread-read))
          (write-buffer spread-read out)
          (write-newline out))
        (copy-stream spread-reads-file out))
        ; merge the unmapped reads
      (with-open-sam (unmapped-file (merge-pathnames input-path (make-pathname :name (format nil "~a-unmapped" input-prefix) :type input-extension)) :direction :input)
                     (skip-sam-header unmapped-file)
                     (copy-stream unmapped-file out)))))

(defun merge-unsorted-files-split-per-chromosome (input-path output input-prefix input-extension header)
  "A function for merging files that were split with elPrep and are unsorted"
  ; merge loop
  (with-open-sam (out output :direction :output)
      (format-sam-header out header)
      ; first merge the unmapped reads
      (with-open-sam (unmapped-file (merge-pathnames input-path (make-pathname :name (format nil "~a-unmapped" input-prefix) :type input-extension)) :direction :input)
                     (skip-sam-header unmapped-file)
                     (copy-stream unmapped-file out))
      ; merge spread reads
      (with-open-sam (spread-reads-file (merge-pathnames input-path (make-pathname :name (format nil "~a-spread" input-prefix) :type input-extension)) :direction :input)
                     (skip-sam-header spread-reads-file)
                     (copy-stream spread-reads-file out))
      ; merge the rest of the files
      (loop for sn-form in (sam-header-sq header)
            for chrom = (getf sn-form :SN)
            for file-name = (merge-pathnames input-path (make-pathname :name (format nil "~a-~a" input-prefix chrom) :type input-extension))
            when (probe-file file-name) do
            (with-open-sam (file file-name :direction :input)
                           (skip-sam-header file)
                           (copy-stream file out)))))

(declaim (inline parse-sam-alignment-from-stream))

(defun parse-sam-alignment-from-stream (stream)
  "Read a line from a stream and parse it as a SAM alignment line.
   Return NIL if stream is at end of file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4."
  (let ((line (read-line stream nil)))
    (when line (parse-sam-alignment line))))

(defun compare-sam-files (file1 file2 &optional (output "/dev/stdout"))
  "A function for comparing two sam files. The input files must be sorted by coordinate order."
  (labels ((get-alns (stream next-aln)
             (loop with group-aln = (or next-aln (parse-sam-alignment-from-stream stream))
                   with alns = (list group-aln)
                   for aln = (parse-sam-alignment-from-stream stream)
                   while (and aln
                              (= (sam-alignment-pos group-aln)
                                 (sam-alignment-pos aln))
                              (string= (sam-alignment-rname group-aln)
                                       (sam-alignment-rname aln)))
                   do (push aln alns)
                   finally (return (values (sort alns (lambda (aln1 aln2)
                                                        (or (string< (sam-alignment-qname aln1)
                                                                     (sam-alignment-qname aln2))
                                                            (when (string= (sam-alignment-qname aln1)
                                                                           (sam-alignment-qname aln2))
                                                              (< (sam-alignment-flag aln1)
                                                                 (sam-alignment-flag aln2))))))
                                           aln))))
           (plist-to-sorted-alist (plist)
             (sort (loop for (key value) on plist by #'cddr collect (cons key value))
                   #'string< :key (lambda (object) (string (car object)))))
           (compare-alns (out alns1 alns2)
             (loop for aln1 in alns1
                   for aln2 in alns2
                   for difference = (cond ((string/= (sam-alignment-qname aln1)
                                                     (sam-alignment-qname aln2)) "qname (1)")
                                          ((/=       (sam-alignment-flag  aln1)
                                                     (sam-alignment-flag  aln2)) "flag (2)")
                                          ((string/= (sam-alignment-rname aln1)
                                                     (sam-alignment-rname aln2)) "rname (3)")
                                          ((/=       (sam-alignment-pos   aln1)
                                                     (sam-alignment-pos   aln2)) "pos (4)")
                                          ((/=       (sam-alignment-mapq  aln1)
                                                     (sam-alignment-mapq  aln2)) "mapq (5)")
                                          ((string/= (sam-alignment-cigar aln1)
                                                     (sam-alignment-cigar aln2)) "cigar (6)")
                                          ((string/= (sam-alignment-rnext aln1)
                                                     (sam-alignment-rnext aln2)) "rnext (7)")
                                          ((/=       (sam-alignment-pnext aln1)
                                                     (sam-alignment-pnext aln2)) "pnext (8)")
                                          ((/=       (sam-alignment-tlen  aln1)
                                                     (sam-alignment-tlen  aln2)) "tlen (9)")
                                          ((string/= (sam-alignment-seq   aln1)
                                                     (sam-alignment-seq   aln2)) "seq (10)")
                                          ((string/= (sam-alignment-qual  aln1)
                                                     (sam-alignment-qual  aln2)) "qual (11)")
                                          (t (let ((tags1 (plist-to-sorted-alist (sam-alignment-tags aln1)))
                                                   (tags2 (plist-to-sorted-alist (sam-alignment-tags aln2))))
                                               (when (or (/= (length tags1) (length tags2))
                                                         (loop for (nil . val1) in tags1
                                                               for (nil . val2) in tags2
                                                               thereis (or (not (eq (type-of val1) (type-of val2)))
                                                                           (etypecase val1
                                                                             (character (char/= val1 val2))
                                                                             (number    (/= val1 val2))
                                                                             (string    (string/= val1 val2))
                                                                             (array     (not (equalp val1 val2)))))))
                                                 "optional tags"))))
                   when difference do
                   (format t "alignments differ for ~a entry: ~%" difference)
                   (format-sam-alignment *standard-output* aln1)
                   (format-sam-alignment *standard-output* aln2)
                   (format-sam-alignment out aln1)
                   (format-sam-alignment out aln2))))
    (with-open-sam (in1 file1 :direction :input)
      (skip-sam-header in1)
      (with-open-sam (in2 file2 :direction :input)
        (skip-sam-header in2)
        (with-open-sam (out output :direction :output)
          (loop 
           for prev-aln1 = nil
           for prev-aln2 = nil
           for alns1 = (multiple-value-bind (alns next) (get-alns in1 prev-aln1) (setf prev-aln1 next) alns)
           for alns2 = (multiple-value-bind (alns next) (get-alns in2 prev-aln2) (setf prev-aln2 next) alns)
           for l1 = (length alns1)
           for l2 = (length alns2)
           for index from 1
           sum l1 into nr-of-reads-matched
           while (or alns1 alns2) do
           (if (= l1 l2)
               (compare-alns out alns1 alns2)
             (let ((pos (sam-alignment-pos (or (first alns1) (first alns2))))
                   (rname (sam-alignment-rname (or (first alns1) (first alns2)))))
               (format t "Files contain an unequal number of read entries at the same position.~%")
               (format t "File ~a has ~a reads at position ~a ~a.~%" file1 l1 rname pos) 
               (format t "File ~a has ~a reads at position ~a ~a.~%" file2 l2 rname pos)))
           (when (zerop (mod index 1000000))
             (format t "~a reads compared and matched.~%" nr-of-reads-matched))
           finally (format t "~a reads compared and matched.~%" nr-of-reads-matched)))))))
