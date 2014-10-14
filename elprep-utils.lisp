(in-package :elprep)

(defun sam-alignment-differ (aln1 aln2)
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  ; check that all mandatory fields are =
  (or (when (string/= (the simple-base-string (sam-alignment-qname aln1)) (the simple-base-string (sam-alignment-qname aln2))) 'qname)
      (when (/= (the fixnum (sam-alignment-flag aln1)) (the fixnum (sam-alignment-flag aln2))) 'flag)
      (when (string/= (the simple-base-string (sam-alignment-rname aln1)) (the simple-base-string (sam-alignment-rname aln2))) 'rname)
      (when (/= (the int32 (sam-alignment-pos aln1)) (the int32 (sam-alignment-pos aln2))) 'pos)
      (when (/= (sam-alignment-mapq aln1) (sam-alignment-mapq aln2)) 'mapq)
      (when (string/= (the simple-base-string (sam-alignment-cigar aln1)) (the simple-base-string (sam-alignment-cigar aln2))) 'cigar)
      (when (string/= (the simple-base-string (sam-alignment-rnext aln1)) (the simple-base-string (sam-alignment-rnext aln2))) 'rnext)
      (when (string/= (the simple-base-string (sam-alignment-qual aln1)) (the simple-base-string (sam-alignment-qual aln2))) 'qual)))

(defun sam-alignment-same (aln1 aln2)
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  (and (string= (the simple-base-string (sam-alignment-qname aln1)) (the simple-base-string (sam-alignment-qname aln2)))
       (= (the fixnum (sam-alignment-flag aln1)) (the fixnum (sam-alignment-flag aln2)))
       (string= (the simple-base-string (sam-alignment-rname aln1)) (the simple-base-string (sam-alignment-rname aln2)))
       (= (the int32 (sam-alignment-pos aln1)) (the int32 (sam-alignment-pos aln2)))
       (= (sam-alignment-mapq aln1) (sam-alignment-mapq aln2))
       (string= (the simple-base-string (sam-alignment-cigar aln1)) (the simple-base-string (sam-alignment-cigar aln2)))
       (string= (the simple-base-string (sam-alignment-rnext aln1)) (the simple-base-string (sam-alignment-rnext aln2)))
       (string= (the simple-base-string (sam-alignment-qual aln1)) (the simple-base-string (sam-alignment-qual aln2)))))

(defun real-diffs (alns1 alns2)
  (loop for aln1 in alns1
        unless (find aln1 alns2 :test 'sam-alignment-same)
        collect aln1))

(defun compare-sams (sam1-file sam2-file)
  ; parse both sams to memory, then do a 1 by 1 comparison on the alignments for all obligatory fields
  (let ((sam1 (make-sam))
        (sam2 (make-sam)))
    (run-pipeline (current-pathname sam1-file) sam1)
    (run-pipeline (current-pathname sam2-file) sam2)
    ; sort the sams by qname
    (setf (sam-alignments sam1) (stable-sort (sam-alignments sam1) 'string< :key 'sam-alignment-qname))
    (setf (sam-alignments sam2) (stable-sort (sam-alignments sam2) 'string< :key 'sam-alignment-qname))
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
  ; assume the input is coordinate sorted; verify if this is still the case
  (format t "verifying order kept ~%")
  (let ((sam (make-sam)))  
    (run-pipeline (current-pathname sam-file) sam)
    (let ((pos (sam-alignment-pos (first (sam-alignments sam))))
          (rname (sam-alignment-rname (first (sam-alignments sam))))
          (ctr 1))
      (loop for aln in (rest (sam-alignments sam))
            do (let ((new-pos (sam-alignment-pos aln))
                     (new-rname (sam-alignment-rname aln)))
                 (cond ((and (< new-pos pos) (string= rname new-rname )) (format t "Not sorted: previous pos: ~s,~s current pos: ~s,~s. ~s reads were in the right order. ~%" rname pos new-rname new-pos ctr) 
                        (return nil))
                       (t 
                        (incf ctr)
                        (setf rname new-rname)
                        (setf pos new-pos))))
            finally (return t)))))

(defun count-duplicates (sam-file)
  (let ((sam (make-sam)))  
    (run-pipeline (current-pathname sam-file) sam)
    (loop for aln in (sam-alignments sam)
          count (sam-alignment-duplicate-p aln))))

; code for splitting up sam files into chromosomes

(defun split-file-per-chromosome (input output-path output-prefix output-extension)
  (with-open-stream (in (open-sam input :direction :input))
    (let* ((header (parse-sam-header in))
           (chroms-encountered (make-hash-table :test 'buffer= :hash-function 'buffer-hash :single-thread t))
           (buf-unmapped (make-buffer)))
      ; fill in a file for unmapped reads
      (reinitialize-buffer buf-unmapped)
      (buffer-extend buf-unmapped "*")
      (setf (gethash buf-unmapped chroms-encountered)
            (let ((file (open-sam (format nil "~a~a-unmapped~a" output-path output-prefix output-extension) :direction :output)))
              (format-sam-header file header) ; fill in the header
              file))
      (loop for sn-form in (sam-header-sq header)
            do (let ((chrom (getf sn-form :SN))
                     (buf-chrom (make-buffer)))
                 (reinitialize-buffer buf-chrom)
                 (buffer-extend buf-chrom chrom)
                 (setf (gethash buf-chrom chroms-encountered)
                       (let ((file (open-sam (format nil "~a~a-~a~a" output-path output-prefix chrom output-extension) :direction :output)))
                         (format-sam-header file header)
                         file))))
      (let ((file-spread-reads (open-sam (format nil "~a~a-spread~a" output-path output-prefix output-extension) :direction :output))
            (buf-= (make-buffer)))
        (format-sam-header file-spread-reads header)
        (buffer-extend buf-= "=") 
        (let ((rname (make-buffer))
              (rnext (make-buffer))
              (aln-string (make-buffer)))
          (loop until (end-of-file-p (peekc in))
                do (progn (reinitialize-buffer aln-string)
                     (reinitialize-buffer rname)
                     (reinitialize-buffer rnext)
                     (read-line-into-buffer in aln-string)
                     (buffer-partition aln-string #\Tab 2 rname 6 rnext)
                     (let ((file (cond ((or (buffer= buf-= rnext) (buffer= buf-unmapped rname) (buffer= rname rnext)) (gethash rname chroms-encountered))
                                       (t file-spread-reads))))
                       (write-buffer aln-string file)))))
        ; close files
        (close file-spread-reads)
        (loop for file being each hash-value of chroms-encountered do (close file))))))

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
  (let ((files (loop for sn-form in (sam-header-sq header)
                     collect (let* ((chrom (getf sn-form :SN))
                                    (file (open-sam (format nil "~a~a-~a~a" input-path input-prefix chrom input-extension) :direction :input)))
                               (skip-sam-header file)
                               file)))
        (unmapped-file (open-sam (format nil "~a~a-unmapped~a" input-path input-prefix input-extension) :direction :input))
        (spread-reads-file (open-sam (format nil "~a~a-spread~a" input-path input-prefix input-extension) :direction :input)))
      ; merge loop
      (with-open-stream (out (open-sam output :direction :output))
        (format-sam-header out header)
        (let ((spread-read (make-buffer)) ; for storing entries from the spread-read file
              (spread-read-refid (make-buffer))
              (spread-read-pos (make-buffer))
              (chromosome-read (make-buffer)) ; for storing reads from the chromsome file we are currently merging
              (chromosome-read-refid (make-buffer))
              (chromosome-read-pos (make-buffer)))
          (flet ((reset-spread-read () (reinitialize-buffer spread-read) (reinitialize-buffer spread-read-pos) (reinitialize-buffer spread-read-refid))
                 (reset-chromosome-read () (reinitialize-buffer chromosome-read) (reinitialize-buffer chromosome-read-pos)))
            (loop for file in files
                  do 
                  (when (buffer-emptyp spread-read) ; if the buffer is not empty, the current entry is potentially an entry for this file and it should not be overwritten
                    (read-line-into-buffer file spread-read)
                    (buffer-partition spread-read #\Tab 2 spread-read-refid 4 spread-read-pos))
                  (read-line-into-buffer file chromosome-read)
                  (buffer-partition chromosome-read #\Tab 2 chromosome-read-refid 4 chromosome-read-pos)
                  (when (buffer= spread-read-refid chromosome-read-refid)
                    (loop do 
                          (let ((pos1 (buffer-parse-integer spread-read-pos))
                                (pos2 (buffer-parse-integer chromosome-read-pos)))
                            (cond ((< pos1 pos2) 
                                   (write-buffer spread-read out) 
                                   (reset-spread-read)
                                   (read-line-into-buffer spread-reads-file spread-read)
                                   (buffer-partition spread-read #\Tab 2 spread-read-refid 4 spread-read-pos))
                                  (t
                                   (write-buffer chromosome-read out)
                                   (reset-chromosome-read)
                                   (read-line-into-buffer file chromosome-read)
                                   (buffer-partition chromosome-read #\Tab 4 chromosome-read-pos))))
                          until (or (not (buffer= chromosome-read-refid spread-read-refid)) (end-of-file-p (peekc file)))
                          finally do (reinitialize-buffer chromosome-read)))
                  ; copy remaining reads in the file; is there a faster way to concatenate files?
                  (loop until (end-of-file-p (peekc file))
                        do 
                        (read-line-into-buffer file chromosome-read)
                        (write-buffer chromosome-read out)
                        (reinitialize-buffer chromosome-read)))
            ; merge the unmapped reads
            (reinitialize-buffer chromosome-read)
            (loop until (end-of-file-p (peekc unmapped-file))
                  do 
                  (read-line-into-buffer unmapped-file chromosome-read)
                  (write-buffer chromosome-read out)
                  (reinitialize-buffer chromosome-read)))))
      ; close the files
      (close unmapped-file)
      (close spread-reads-file)
      (loop for file in files do (close file))))

(defun compare-sam-files (file1 file2 &optional (output "/dev/stdout"))
  (labels ((get-alns (stream next-aln)
             (loop for aln = (or next-aln (elprep:parse-sam-alignment stream))
                   then (elprep:parse-sam-alignment stream)
                   while (and aln (or (null alns)
                                      (= (elprep:sam-alignment-pos (first alns))
                                         (elprep:sam-alignment-pos aln))))
                   collect aln into alns
                   finally (return (list (sort alns (lambda (aln1 aln2)
                                                      (or (string< (elprep:sam-alignment-qname aln1)
                                                                   (elprep:sam-alignment-qname aln2))
                                                          (when (string= (elprep:sam-alignment-qname aln1)
                                                                         (elprep:sam-alignment-qname aln2))
                                                            (< (elprep:sam-alignment-flag aln1)
                                                               (elprep:sam-alignment-flag aln2))))))
                                         aln))))
           (plist-to-sorted-alist (plist)
             (sort (loop for (key value) on plist by 'cddr collect (cons key value))
                   'string< :key (lambda (object) (string (car object)))))
           (compare-alns (out alns1 alns2)
             (loop for aln1 in alns1
                   for aln2 in alns2
                   always
                   (let ((difference nil))
                     (or (and (or (string= (elprep:sam-alignment-qname aln1)
                                           (elprep:sam-alignment-qname aln2))
                                  (setf difference "qname (1)") nil)
                              (or (= (elprep:sam-alignment-flag aln1)
                                     (elprep:sam-alignment-flag aln2))
                                  (prog1 nil (setf difference "flag (2)")))
                              (or (string= (elprep:sam-alignment-rname aln1)
                                           (elprep:sam-alignment-rname aln2))
                                  (prog1 nil (setf difference "rname (3)")))
                              (or (= (elprep:sam-alignment-pos aln1)
                                     (elprep:sam-alignment-pos aln2))
                                  (prog1 nil (setf difference "pos (4)")))
                              (or (= (elprep:sam-alignment-mapq aln1)
                                     (elprep:sam-alignment-mapq aln2))
                                  (prog1 nil (setf difference "mapq (5)")))
                              (or (string= (elprep:sam-alignment-cigar aln1)
                                           (elprep:sam-alignment-cigar aln2))
                                  (prog1 nil (setf difference "cigar (6)")))
                              (or (string= (elprep:sam-alignment-rnext aln1)
                                           (elprep:sam-alignment-rnext aln2))
                                  (prog1 nil (setf difference "rnext (7)")))
                              (or (= (elprep:sam-alignment-pnext aln1)
                                     (elprep:sam-alignment-pnext aln2))
                                  (prog1 nil (setf difference "pnext (8)")))
                              (or (= (elprep:sam-alignment-tlen aln1)
                                     (elprep:sam-alignment-tlen aln2))
                                  (prog1 nil (setf difference "tlen (9)")))
                              (or (string= (elprep:sam-alignment-seq aln1)
                                           (elprep:sam-alignment-seq aln2))
                                  (prog1 nil (setf difference "seq (10)")))
                              (or (string= (elprep:sam-alignment-qual aln1)
                                           (elprep:sam-alignment-qual aln2))
                                  (prog1 nil (setf difference "qual (11)")))
                              (let ((tags1 (plist-to-sorted-alist (elprep:sam-alignment-tags aln1)))
                                    (tags2 (plist-to-sorted-alist (elprep:sam-alignment-tags aln2))))
                                (or (and (= (length tags1) (length tags2))
                                         (loop for (nil . val1) in tags1
                                               for (nil . val2) in tags2
                                               always (and (eq (type-of val1) (type-of val2))
                                                           (etypecase val1
                                                             (character (char= val1 val2))
                                                             (number    (= val1 val2))
                                                             (string    (string= val1 val2))
                                                             (array     (equalp val1 val2))))))
                                    (prog1 nil (setf difference "optional tags")))))
                         (prog1 nil
                           (format t "alignments differ for ~a entry: ~%" difference)
                           (format-sam-alignment out aln1)
                           (format-sam-alignment out aln2))))))
           (make-alns-process (file)
             (let ((mailbox (mp:make-mailbox)))
               (mp:process-run-function 
                "get-alns" '(:internal-server T)
                (lambda (file mailbox)
                  (with-open-stream (stream (open-sam file :direction :input))
                    (elprep:parse-sam-header stream)
                    (loop for (alns next-aln) = (get-alns stream next-aln)
                          while alns do (mp:mailbox-send mailbox alns)
                          finally (mp:mailbox-send mailbox nil))))
                file mailbox)
               mailbox)))
    (with-open-stream (out (open-sam output :direction :output))
      (loop with mailbox1 = (make-alns-process file1)
            with mailbox2 = (make-alns-process file2)
            for alns1 = (mp:mailbox-read mailbox1)
            for alns2 = (mp:mailbox-read mailbox2)
            for index from 1
            while (or alns1 alns2)
            always (let ((l1 (length alns1))
                         (l2 (length alns2)))
                     (and (or (= l1 l2)
                              (prog1 nil
                                (let ((pos (if (> l1 0) (sam-alignment-pos (first alns1)) (sam-alignment-pos (first alns2)))))
                                  (format t "Files contain an unequal number of read entries at the same position.~%")
                                  (format t "File ~a has ~a reads at position ~a.~%" file1 l1 pos) 
                                  (format t "File ~a has ~a reads at position ~a.~%" file2 l2 pos))))
                          (compare-alns out alns1 alns2)))
            do (when (zerop (mod index 1000000)) (format t "~a reads compared and matched.~%" index)) 
            finally (format t "~a reads compared and matched.~%" index)))))
