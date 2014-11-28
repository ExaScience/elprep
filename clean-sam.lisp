(in-package :elprep)
(in-simple-base-string-syntax)

(eval-when (#+sbcl :compile-toplevel :load-toplevel :execute)
  (defun make-cigar-operations-table-consumes-reference-bases ()
    "Map CIGAR operations to boolean flags indicating whether they consume reference bases."
    (let ((table (make-array (1+ (- +max-cigar-operation+
                                    +min-cigar-operation+))
                             :initial-element nil
                             #+lispworks :allocation #+lispworks :long-lived
                             #+lispworks :single-thread #+lispworks t)))
      (loop for char across "mdn=x" do (setf (svref table (- (char-code char) +min-cigar-operation+)) t))
      table)))

(defglobal *cigar-consumes-reference-bases-table*
  (make-cigar-operations-table-consumes-reference-bases)
  "Map CIGAR operations to boolean flags indicating whether they consume reference bases.")

(defun sam-alignment-end (aln cigar-list)
  "Sums the lengths of all CIGAR operations that consume reference bases."
  (declare (sam-alignment aln) (list cigar-list) #.*optimization*)
  (let ((cigar-reference-length
         (loop with table = *cigar-consumes-reference-bases-table*
               for (key . val) of-type (base-char . int32) in cigar-list
               when (svref table (- (char-code key) +min-cigar-operation+))
               ;when (or (eq key :M) (eq key :D) (eq key :N) (eq key :=) (eq key :X))
               sum val)))
    (declare (int32 cigar-reference-length))
    (+ (sam-alignment-pos aln) cigar-reference-length -1)))

(declaim (inline reference-sequence-length))

(defun reference-sequence-length (aln reference-sequence-table)
  "Get reference sequence length for the given sam-alignment."
  (declare (sam-alignment aln) (hash-table reference-sequence-table) #.*optimization*)
  (gethash (sam-alignment-rname aln) reference-sequence-table))

(declaim (inline operator-consumes-read-bases-p))

(defun operator-consumes-read-bases-p (operator)
  "Does the CIGAR operation consume read bases?"
  (declare (base-char operator) #.*optimization*)
  (or (char= operator #\M) (char= operator #\I) (char= operator #\S) (char= operator #\=) (char= operator #\X)))

(declaim (inline operator-consumes-reference-bases-p))

(defun operator-consumes-reference-bases-p (operator)
  "Does the CIGAR operation consume reference bases?"
  (declare (base-char operator) #.*optimization*)
  (or (char= operator #\M) (char= operator #\D) (char= operator #\N) (char= operator #\=) (char= operator #\X)))

(defun get-read-length-from-cigar (cigar-list)
  "Sums the lengths of all CIGAR operations that consume read bases."
  (declare (list cigar-list) #.*optimization*)
  (loop for (operator . val) of-type (base-char . int32) in cigar-list
        when (operator-consumes-read-bases-p operator) 
        sum val))

(declaim (inline element-stradless-clipped-read))

(defun element-stradless-clipped-read (new-cigar operator relative-clipping-position clipped-bases)
  (declare (stream new-cigar) (base-char operator) (int32 relative-clipping-position clipped-bases) #.*optimization*)
  (if (operator-consumes-read-bases-p operator)
    (if (operator-consumes-reference-bases-p operator)
      (when (> relative-clipping-position 0)
        (format new-cigar "~C~D" operator relative-clipping-position))
      (setf clipped-bases (+ clipped-bases relative-clipping-position)))
    (when (/= relative-clipping-position 0)
      (error "Unexpected non-0 relative clipping position: ~s" relative-clipping-position)))
  (format new-cigar "S~D" clipped-bases))
  
(defun soft-clip-end-of-read (clip-from cigar-list)
  (declare (int32 clip-from) (list cigar-list) #.*optimization*)
  (let ((pos 0))
    (declare (int32 pos))
    (decf clip-from)
    (with-output-to-string (new-cigar nil :element-type 'base-char)
      (loop named cigar-loop for (operator . val) of-type (base-char . int32) in cigar-list
            do (let ((end-pos (+ pos (if (operator-consumes-read-bases-p operator) val 0))))
                 (declare (int32 end-pos))
                 (if (< end-pos clip-from)
                   (format new-cigar "~C~D" operator val)
                   (let ((clipped-bases (+ (get-read-length-from-cigar cigar-list) clip-from)))
                     (declare (int32 clipped-bases))
                     (element-stradless-clipped-read new-cigar operator (- clip-from pos) clipped-bases)
                     (return-from cigar-loop)))
                 (setf pos (+ end-pos pos)))))))

(defun clean-sam (header)
  "A filter for soft-clipping a sam-alignment at the end of a reference sequence, and set MAPQ to 0 if unmapped."
  (let ((reference-sequence-table ; create a reference dictionary that maps reference name onto the length of that reference
         (make-single-thread-hash-table :test #'equal)))
    (loop for sn-form in (sam-header-sq header)
          do (setf (gethash (getf sn-form :SN) reference-sequence-table) (getf sn-form :LN)))
    (lambda ()
      (lambda (alignment)
        (declare (sam-alignment alignment) #.*optimization*)
        (if (sam-alignment-unmapped-p alignment)
          (setf (sam-alignment-mapq alignment) 0)
          (let ((ref-seq-length (reference-sequence-length alignment reference-sequence-table))
                (scanned-cigar (scan-cigar-string 'list (sam-alignment-cigar alignment))))
            (declare (int32 ref-seq-length) (list scanned-cigar))
            (when (> (sam-alignment-end alignment scanned-cigar) ref-seq-length)
              (let ((clip-from (+ (- ref-seq-length (sam-alignment-pos alignment)) 1)))
                (setf (sam-alignment-cigar alignment) (soft-clip-end-of-read clip-from scanned-cigar))))))
        t))))

  
     