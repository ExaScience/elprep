;;; Copyright (c) 2013, 2014 by Charlotte Herzeel, IMEC vzw.
;;; IMEC Confidential.
;;; Copyright (c) 2013, 2014 by Pascal Costanza, Intel Corporation.
;;; Intel Confidential.

(in-package :sam-filters)

(defun NA12878_mapping ()
  (let ((mapping (list "1" "chr1" "2" "chr2" "3" "chr3" "4" "chr4" "5" "chr5" "6" "chr6" 
                       "7" "chr7" "8" "chr8" "9" "chr9" "10" "chr10" "11" "chr11" "12" "chr12"
                       "13" "chr13" "14" "chr14" "15" "chr15" "16" "chr16" "17" "chr17" 
                       "18" "chr18" "19" "chr19" "20" "chr20" "21" "chr21" "22" "chr22"
                       "X" "chrX" "Y" "chrY" "MT" "chrMT")))
    (let ((table (make-hash-table :test 'equal :single-thread t)))
      (loop for old of-type string in mapping by (function cddr)
            for new of-type string in (rest mapping) by (function cddr)
            do (setf (gethash old table) new))
      table)))

(defun rename-reference-sequence-dictionary (mapping)
  (declare (hash-table mapping))
  (lambda (header)
    (let ((sq (sam-header-sq header)))
      (loop for sn in sq
            do (let* ((old (getf sn :sn))
                      (new (gethash old mapping))) 
                 (setf (getf sn :sn) new)))
      ;(break)
      (lambda ()
        (lambda (alignment)
          (declare (sam-alignment alignment)
                   (optimize (speed 3) (space 0) (debug 0) (safety 0)
                             (compilation-speed 0)))
          (let ((new (gethash (sam-alignment-rname alignment) header)))
            (when new
              (setf (sam-alignment-rname alignment) (the string new))
              alignment)))))))
