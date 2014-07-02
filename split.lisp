;;; Copyright (c) 2013, 2014 by Charlotte Herzeel, IMEC vzw.
;;; IMEC Confidential.
;;; Copyright (c) 2013, 2014 by Pascal Costanza, Intel Corporation.
;;; Intel Confidential.

(defun splitter-process (mbox output) 
  (mp:process-run-function "splitter" '() (lambda ()
                                            (do ((msg (mp:mailbox-read mbox) (mp:mailbox-read mbox)))
                                                ((eql msg :done))
                                              (with-open-file (out (format nil "~s-~s" output (sam-alignment-rname msg)) :direction :output :if-exists :append)
                                                (format-sam-alignment out aln))))))

(defmethod split-file-per-chromosome ((input pathname))
  (let ((nr-of-threads *nr-of-threads*))
    (with-open-file (in input :direction :input)
      (let ((header (parse-sam-header in)))
        ; each thread will be responsible for processing a chromosome (or a bunch of chromosomes)
        (let ((workers (make-array number-of-threads)))
          (dotimes (i nr-of-threads)
            (let ((mbox (mp:make-mailbox)))
              (setf (aref workers i) mbox)
              (make-splitter-worker mbox input)))
          (do ((aln (parse-sam-alignment in) (parse-sam-alignment in)))
              ((not aln) (dotimes (i nr-of-threads) (mp:mailbox-send (aref workers i) :done)))
            (let ((mbox (aref workers (mod (sam-alignment-refid aln) nr-of-threads))))
              (mp:mailbox-send mbox aln))))))))  
