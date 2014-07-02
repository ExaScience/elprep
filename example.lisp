(in-package :cl-user)

(defun best-practices-pipeline-no-filters (sam-file sam-out-file)
  (let ((hdr nil))
    (format t "~%Feed SAM file into Memory.~%")
    (time (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file)
                                    :filters (list (lambda (header)
                                                     (setf hdr header)
                                                     (lambda ()
                                                       (lambda (aln) 
                                                         aln))))))
    hdr))

(defun best-practices-pipeline-simple-filters (sam-file dict-file sam-out-file)
  (format t "~%Feed SAM file into Memory.~%")
  (time (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file) 
                                  :filters (list (sam-filters:replace-reference-sequence-dictionary-from-sam-file dict-file)
                                                 'sam-filters:filter-unmapped-reads
                                                 (sam-filters:add-or-replace-read-group
                                                  '(:ID "g"
                                                    :LB "lib1"
                                                    :PL "illumina"
                                                    :PU "unit1"
                                                    :SM "sample1"))))))

(defun best-practices-pipeline-rename-chromosome (sam-file sam-out-file)
  (format t "~%Feed SAM file into Memory.~%")
  (time (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file) 
                                  :filters (list 'sam-filters::rename-chromosomes
                                                 ))))

(defun best-practices-remove* (sam-file sam-out-file)
  (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file) 
                            :filters (list (lambda (header) (lambda () (lambda (alignment) (let ((rname (sam-alignment-rname alignment)))
                                                                                             (not (or (string= rname "=") (string= rname "*"))))))))
                                                 ))

(defun best-practices-pipeline-one-simple-filter (sam-file sam-out-file)
  (format t "~%Feed SAM file into Memory.~%")
  (time (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file) 
                                  :filters (list (sf::rename-reference-sequence-dictionary (sf::na12878_mapping))))))


(defun best-practices-pipeline-only-simple-filters (sam-file dict-file)
  (format t "~%Feed SAM file into Memory.~%")
  (time (sam-filters:run-pipeline (current-pathname sam-file)  nil
                                  :filters (list (sam-filters:replace-reference-sequence-dictionary-from-sam-file dict-file)
                                                 'sam-filters:filter-unmapped-reads
                                                 (sam-filters:add-or-replace-read-group
                                                  '(:ID "g"
                                                    :LB "lib1"
                                                    :PL "illumina"
                                                    :PU "unit1"
                                                    :SM "sample1"))))))

#|(defun best-practices-pipeline-all-filters (sam-file dict-file sam-out-file nof-reads)
  (format t "~%Feed SAM file into Memory.~%")
  (time (sam-filters:run-pipeline (current-pathname sam-file) (current-pathname sam-out-file) 
                                  :filters (list (sam-filters:replace-reference-sequence-dictionary-from-sam-file dict-file)
                                                 ;'sam-filters:clean-sam
                                                 'sam-filters:filter-unmapped-reads
                                                 (sam-filters:add-or-replace-read-group
                                                  '(:ID "g"
                                                    :LB "lib1"
                                                    :PL "illumina"
                                                    :PU "unit1"
                                                    :SM "sample1"))
                                                 (sam-filters:mark-duplicates-filter nof-reads)))))|#

(defun best-practices-pipeline-all-filters (sam-file dict-file sam-out-file &key (sorting-order :keep sorting-order-p))
  (format t "~%Feed SAM file into Memory.~%")
  (let ((deduplicated-reads (sam-filters:make-sam)))
    (time (sam-filters:run-pipeline (current-pathname sam-file) deduplicated-reads 
                                    :filters (list (sam-filters:replace-reference-sequence-dictionary-from-sam-file dict-file)
                                                 ;'sam-filters:clean-sam
                                                   'sam-filters:filter-unmapped-reads
                                                   (sam-filters:add-or-replace-read-group
                                                    '(:ID "g"
                                                      :LB "lib1"
                                                      :PL "illumina"
                                                      :PU "unit1"
                                                      :SM "sample1"))
                                                   'sam-filters:mark-duplicates-filter)))
    ; sort
    (let* ((effective-sorting-order :keep)
           (list-of-reads (sam-filters:sam-alignments deduplicated-reads))
           (hd (sam-filters:sam-header-hd (sam-filters:sam-header deduplicated-reads)))
           (vn (or (getf hd :vn) "1.4"))
           (so (or (getf hd :so) "unknown")))
      (ecase sorting-order
        (:keep)
        (:unknown    (setq so "unknown"))
        (:unsorted   (setq so "unsorted"))
        (:queryname  (unless (string= so "queryname")
                       (setq effective-sorting-order :queryname so "queryname")))
        (:coordinate (unless (string= so "coordinate")
                       (setq effective-sorting-order :coordinate so "coordinate"))))
      (unless (eq effective-sorting-order :keep) 
        (setf (sam-filters:sam-alignments deduplicated-reads) ; sort the list
              (ecase effective-sorting-order
                (:queryname  (stable-sort list-of-reads 'string< :key 'alignment-qname))
                (:coordinate (stable-sort list-of-reads 'sam-filters:coordinate<)))))
      ; fill in sorting order in header
      (cond (hd (setf (getf hd :vn) vn) (setf (getf hd :so) so))
            (t
             (setf (sam-filters:sam-header-hd (sam-filters:sam-header deduplicated-reads))
                   (list :vn vn :so so))))
      ;(break deduplicated-reads)
      ; write to file
      (sam-filters:run-pipeline deduplicated-reads (current-pathname sam-out-file)))))
  
#|(defstruct sam-alignment
  (qname   "" :type simple-base-string)
  (flag    0  :type fixnum)
  (rname   "" :type simple-base-string)
  (pos     0  :type int32)
  (mapq    0  :type fixnum)
  (cigar   "" :type simple-base-string)
  (rnext   "" :type simple-base-string)
  (pnext   0  :type int32)
  (tlen    0  :type int32)
  (seq     "" :type simple-base-string)
  (qual    "" :type simple-base-string)
  (tags   '() :type list)  ; additional tags in SAM files
  (xtags  '() :type list)  ; additional tags for other storage formats
  (temps  '() :type list))
|#

(defun sam-alignment-differ (aln1 aln2)
  (declare (sf:sam-alignment aln1 aln2) 
           (optimize (speed 3) (space 0) (safety 0) (debug 0)
                     (compilation-speed 0)))
  ; check that all mandatory fields are =
  (or (when (not (string= (the simple-base-string (sf:sam-alignment-qname aln1)) (the simple-base-string (sf:sam-alignment-qname aln2)))) 'qname)
      (when (not (= (the fixnum (sf:sam-alignment-flag aln1)) (the fixnum (sf:sam-alignment-flag aln2)))) 'flag)
      (when (not (string= (the simple-base-string (sf:sam-alignment-rname aln1)) (the simple-base-string (sf:sam-alignment-rname aln2)))) 'rname)
      (when (not (= (the sf:int32 (sf:sam-alignment-pos aln1)) (the sf:int32 (sf:sam-alignment-pos aln2)))) 'pos)
      (when (not (= (sf:sam-alignment-mapq aln1) (sf:sam-alignment-mapq aln2))) 'mapq)
      (when (not (string= (the simple-base-string (sf:sam-alignment-cigar aln1)) (the simple-base-string (sf:sam-alignment-cigar aln2)))) 'cigar)
      (when (not (string= (the simple-base-string (sf:sam-alignment-rnext aln1)) (the simple-base-string (sf:sam-alignment-rnext aln2)))) 'rnext)
      (when (not (string= (the simple-base-string (sf:sam-alignment-qual aln1)) (the simple-base-string (sf:sam-alignment-qual aln2)))) 'qual)))

(defun compare-sams (sam1-file sam2-file)
  ; parse both sams to memory, then do a 1 by 1 comparison on the alignments for all obligatory fields
  (let ((sam1 (sf:make-sam))
        (sam2 (sf:make-sam)))
    (sf:run-pipeline (current-pathname sam1-file) sam1 :filters (list (lambda (header) (declare (ignore header)) (lambda () (lambda (aln) aln)))) :ordered t)
    (sf:run-pipeline (current-pathname sam2-file) sam2 :filters (list (lambda (header) (declare (ignore header)) (lambda () (lambda (aln) aln)))) :ordered t)
    (format t "sam1:~s alns sam2:~s alns ~%" (length (sf:sam-alignments sam1)) (length (sf:sam-alignments sam2)))
    (let ((differences
           (loop for aln1 in (sf:sam-alignments sam1)
                 for aln2 in (sf:sam-alignments sam2)
                 append (let ((d (sam-alignment-differ aln1 aln2))) (when d (list (list d aln1 aln2)))))))
      (values (when differences t) differences))))

(defun spit-out-alns (sam-file)
  (let ((sam1 (sf:make-sam)))
    (sf:run-pipeline (current-pathname sam-file) sam1 :filters (list (lambda (header) (declare (ignore header)) (lambda () (lambda (aln) aln)))) :ordered t)
    sam1))

;/scratch/caherzee/new-iprep