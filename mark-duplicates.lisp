(in-package :elprep)
(in-simple-base-string-syntax)

(declaim (inline sam-alignment-rg (setf sam-alignment-rg)))

(defun sam-alignment-rg (aln)
  "Access the read group optional field of a sam-alignment."
  (sam-alignment-tag aln :rg))

(defun (setf sam-alignment-rg) (new-value aln)
  "Access the read group optional field of a sam-alignment."
  (setf (sam-alignment-tag aln :rg) new-value))

(defun make-phred-score-table ()
  "Map Phred qualities to a reasonable range and an error flag indicating if it is outside a valid range."
  (let ((score-table (make-array 512 :element-type 'octet
                                 #+lispworks :allocation #+lispworks :long-lived
                                 #+lispworks :single-thread #+lispworks t)))
    (loop for char below 256
          for pos = (ash char 1) do
          (if (or (< char 33) (> char 126))
            (setf (aref score-table pos)      0
                  (aref score-table (1+ pos)) 1)
            (let ((qual (- char 33)))
              (setf (aref score-table pos)      (if (>= qual 15) qual 0)
                    (aref score-table (1+ pos)) 0))))
    score-table))

(declaim (notinline compute-phred-score))

(defun compute-phred-score (aln)
  "Sum the adapted Phred qualities of a sam-alignment."
  (declare (sam-alignment aln) #.*optimization*)
  (let ((string (sam-alignment-qual aln))
        (score-table (load-time-value (make-phred-score-table) t))
        (score 0) (error 0))
    (declare (base-string string) ((simple-array octet (512)) score-table) (fixnum score error))
    (multiple-value-bind (string* offset) (unwrap-displaced-array string)
      (declare (simple-base-string string*) (fixnum offset))
      (let ((end (the fixnum (+ (length string) offset))))
        (declare (fixnum end))
        (loop for i of-type fixnum from offset below end
              for pos of-type fixnum = (ash (char-code (schar string* i)) 1) do
              (setq score (+ score (aref score-table pos))
                    error (logior error (aref score-table (the fixnum (1+ pos))))))
        (assert (= error 0))))
    score))

(define-symbol-macro upcase-cigar-operations "MIDNSHPX=")

(defconstant +min-upcase-cigar-operation+ (reduce #'min upcase-cigar-operations :key #'char-code)
  "The smallest CIGAR operation, upper case only.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(defconstant +max-upcase-cigar-operation+ (reduce #'max upcase-cigar-operations :key #'char-code)
  "The largest CIGAR operation, upper case only.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(declaim (inline cigar-aux-pos))

(defun cigar-aux-pos (char)
  "Position in a vector of CIGAR operations, starting with the smallest CIGAR operation, upper case only."
  (declare (base-char char) #.*optimization*)
  (- (char-code char) +min-upcase-cigar-operation+))

(defun make-unclipped-aux-tables ()
  "Map CIGAR operations to flags indicating whether they are clipped and/or reference operations, upper case only."
  (let* ((tablesize (1+ (- +max-upcase-cigar-operation+
                           +min-upcase-cigar-operation+)))
         (clipped   (make-array tablesize
                                :initial-element 0
                                :element-type 'octet
                                #+lispworks :allocation #+lispworks :long-lived
                                #+lispworks :single-thread #+lispworks t))
         (reference (make-array tablesize
                                :initial-element 0
                                :element-type 'octet
                                #+lispworks :allocation #+lispworks :long-lived
                                #+lispworks :single-thread #+lispworks t)))
    (flet ((set-clipped-flag (char)
             (setf (aref clipped (cigar-aux-pos char)) 1))
           (set-reference-flag (char)
             (setf (aref reference (cigar-aux-pos char)) 1)))
      (set-clipped-flag #\S)
      (set-clipped-flag #\H)
      (set-reference-flag #\M)
      (set-reference-flag #\D)
      (set-reference-flag #\N)
      (set-reference-flag #\=)
      (set-reference-flag #\X))
    (cons clipped reference)))

(declaim (notinline %compute-unclipped-position))

(defun %compute-unclipped-position (flag pos cigar)
  "Compute unclipped position of a sam-alignment, based on its FLAG, POS, and CIGAR string."
  (declare (uint16 flag) (int32 pos) (simple-vector cigar) #.*fixnum-optimization*)
  (if (= 0 (length cigar)) pos
    (let* ((tables (load-time-value (make-unclipped-aux-tables) t))
           (clipped-table (car tables))
           (reference-table (cdr tables)))
      (declare (cons tables) ((simple-array octet) clipped-table reference-table))
      (the int32
           (if (/= (logand flag +reversed+) )
             (1- (+ pos (loop for i of-type fixnum = (the fixnum (1- (length cigar))) then (the fixnum (1- i))
                              for (key . value) of-type (base-char . fixnum) = (svref cigar i)
                              for p of-type fixnum = (cigar-aux-pos key)
                              for c of-type fixnum = (aref clipped-table p)
                              for r of-type fixnum = (aref reference-table p)
                              for clipped of-type fixnum = c then (the fixnum (* c clipped))
                              sum (the fixnum (* (logior r clipped) value)) fixnum
                              until (= i 0))))
             (- pos (loop for i of-type fixnum below (length cigar)
                          for (key . value) of-type (base-char . fixnum) = (svref cigar i)
                          for p of-type fixnum = (cigar-aux-pos key)
                          until (= 0 (aref clipped-table p))
                          sum value fixnum)))))))

(declaim (inline compute-unclipped-position))

(defun compute-unclipped-position (aln)
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (let ((flag  (sam-alignment-flag aln))
        (pos   (sam-alignment-pos aln))
        (cigar (scan-cigar-string 'vector (sam-alignment-cigar aln))))
    (declare (uint16 flag) (int32 pos) (simple-vector cigar))
    (%compute-unclipped-position flag pos cigar)))

(declaim (inline sam-alignment-adapted-pos (setf sam-alignment-adapted-pos)
                 sam-alignment-adapted-score (setf sam-alignment-adapted-score)))

(defun sam-alignment-adapted-pos (aln)
  "Access the unclipped position temporary field in the sam-alignment."
  (sam-alignment-temp aln :pos))

(defun (setf sam-alignment-adapted-pos) (val aln)
  "Access the unclipped position temporary field in the sam-alignment."
  (setf (sam-alignment-temp aln :pos) val))

(defun sam-alignment-adapted-score (aln)
  "Access adapted Phred score temporary field in the sam-alignment."
  (sam-alignment-temp aln :score))

(defun (setf sam-alignment-adapted-score) (val aln)
  "Access adapted Phred score temporary field in the sam-alignment."
  (setf (sam-alignment-temp aln :score) val))

(declaim (inline adapt-alignment))

(defun adapt-alignment (rg-table alignment)
  "Adapt the sam-alignment: Make read group unique; fill in unclipped position; fill in Phred score."
  ; make rg unique
  (let ((rg (sam-alignment-rg alignment)))
    (when rg
      (setf (sam-alignment-rg alignment)
            (or (gethash rg rg-table)
                (with-modify-hash (key value found) (rg-table rg)
                  (if found value rg))))))
  ; adapt position
  (setf (sam-alignment-adapted-pos alignment)
        (compute-unclipped-position alignment))
  ; adapt score
  (setf (sam-alignment-adapted-score alignment) 
        (compute-phred-score alignment))
  alignment)

(declaim (inline mark-as-duplicate))

(defun mark-as-duplicate (aln)
  "Set the PCR/optical duplicate FLAG in the sam-alignment."
  (declare (sam-alignment aln) #.*optimization*)
  (setf (sam-alignment-flag aln) (logior (sam-alignment-flag aln) +duplicate+))
  t)

(declaim (inline make-handle))

(defstruct (handle (:constructor make-handle (object hash)))
  "An indirection for compare-and-swap operations.
   The struct handle has a constructor make-handle that takes an object and a hash value as parameters.
   Read-only accessor handle-hash of type fixnum refers to the hash value.
   Accessor handle-object referes to the object.
   Primary use of this struct is to enable compare-and-swap operations in classify-fragment and classify-pair."
  (hash 0 :type fixnum :read-only t)
  (object nil))

(setf (documentation 'make-handle 'function)
      "Constructor for struct handle that takes an object and a hash value as parameters."
      (documentation 'handle-p 'function)
      "Default predicate for struct handle."
      (documentation 'copy-handle 'function)
      "Default copier function for struct handle."
      (documentation 'handle-hash 'function)
      "Read the handle hash value of type fixnum."
      (documentation 'handle-object 'function)
      "Access the handle object.")

(defun handle-fragment= (f1 f2)
  "Do the two handles refer to sam-alignment instances at the same unclipped position and with the same direction?"
  (declare (handle f1 f2) #.*optimization*)
  (let ((f1 (handle-object f1))
        (f2 (handle-object f2)))
    (declare (sam-alignment f1 f2))
    (and (eq (sam-alignment-rg f1) (sam-alignment-rg f2))
         (= (the int32 (sam-alignment-refid f1)) (the int32 (sam-alignment-refid f2)))
         (= (the int32 (sam-alignment-adapted-pos f1)) (the int32 (sam-alignment-adapted-pos f2)))
         (eq (sam-alignment-reversed-p  f1) (sam-alignment-reversed-p  f2)))))

(defun fragment-hash (f)
  "Hash function that corresponds to handle-fragment=, but operates an a sam-alignment directly."
  (declare (sam-alignment f) #.*optimization*)
  (logxor (sxhash (the base-string (sam-alignment-rg f)))
          (sxhash (the int32 (sam-alignment-refid f)))
          (sxhash (the int32 (sam-alignment-adapted-pos f)))
          (sxhash (the boolean (sam-alignment-reversed-p f)))))

(declaim (inline true-fragment-p true-pair-p))

(defun true-fragment-p (aln)
  "Is this sam-alignment definitely not part of a pair?"
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) #.(logior +multiple+ +next-unmapped+)) +multiple+))

(defun true-pair-p (aln)
  "Is this sam-alignment definitely part of a pair?"
  (declare (sam-alignment aln) #.*optimization*)
  (= (logand (sam-alignment-flag aln) #.(logior +multiple+ +next-unmapped+)) +multiple+))

(defun classify-fragment (aln fragments deterministic) ; fragments is the hash table to store intermediate "best" alns, passed as a parameter
  "For each list of sam-alignment instances with the same unclipped position and direction, all except the one with the highest score are marked as duplicates.
   If there are fragments in such a list that are actually part of pairs, all the true fragments are marked as duplicates and the pairs are left untouched."
  (let* ((hash (fragment-hash aln))
         (key (make-handle aln hash))
         (split (hash-table-split key fragments))
         (best (or (gethash key split)
                   (let ((entry (make-handle aln hash)))
                     (with-hash-table-locked split
                       (let ((value (gethash key split)))
                         (cond (value value)
                               (t (setf (gethash entry split) entry)
                                  (return-from classify-fragment)))))))))
    (declare (dynamic-extent key))
    (cond ((true-fragment-p aln)
           (loop with aln-score = (sam-alignment-adapted-score aln)
                 for best-aln = (handle-object best)
                 until (if (true-pair-p best-aln)
                         (mark-as-duplicate aln)
                         (let ((best-aln-score (sam-alignment-adapted-score best-aln)))
                           (cond ((> best-aln-score aln-score)
                                  (mark-as-duplicate aln))
                                 ((= best-aln-score aln-score)
                                  (if deterministic
                                    (if (string> (sam-alignment-qname aln) (sam-alignment-qname best-aln))
                                      (mark-as-duplicate aln)
                                      (when (compare-and-swap (handle-object best) best-aln aln)
                                        (mark-as-duplicate best-aln)))
                                    (mark-as-duplicate aln)))
                                 ((compare-and-swap (handle-object best) best-aln aln)
                                  (mark-as-duplicate best-aln)))))))
          (t ; the aln is a true pair object, there may be a true fragment stored in the hash table which we then need to mark and swap out
           (loop for best-aln = (handle-object best)
                 until (cond ((true-pair-p best-aln) t) ; stop, the best in the hash tab is a pair, this is marked via mark-duplicates
                             ((compare-and-swap (handle-object best) best-aln aln)
                              (mark-as-duplicate best-aln))))))))

(defun sam-alignment-pair= (aln1 aln2)
  "Are the two sam-alignment fragments part of the same pair?"
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  (and (eq      (sam-alignment-rg    aln1) (sam-alignment-rg    aln2))
       (string= (sam-alignment-qname aln1) (sam-alignment-qname aln2))))

(defun sam-alignment-pair-hash (aln)
  "Hash function that corresponds to sam-alignment-pair=."
  (declare (sam-alignment aln) #.*optimization*)
  (or (sam-alignment-temp aln :pair-hash)
      (let ((hash (logxor (sxhash (the base-string (sam-alignment-rg aln)))
                          (sxhash (the base-string (sam-alignment-qname aln))))))
        (setf (sam-alignment-temps aln)
              (list* :pair-hash hash (sam-alignment-temps aln)))
        hash)))

(declaim (inline make-pair))

(defstruct (pair (:constructor make-pair (score aln1 aln2)))
  "A pair consisting of two fragment sam-alignment instances.
   The struct pair has a constructor make-pair that takes a score and two sam-alignment instances as parameters.
   Read-only accessor pair-score of type fixnum refers to the score.
   Read-only accessor pair-aln1 of type sam-alignment refers to the first fragment.
   Read-only accessor pair-aln2 of type sam-alignment refers to the second fragment."
  (score  0 :type fixnum :read-only t)
  (aln1 nil :type (or null sam-alignment) :read-only t)
  (aln2 nil :type (or null sam-alignment) :read-only t))

(setf (documentation 'make-pair 'function)
      "Constructor for struct pair that takes a score and two sam-alignment instances as parameters."
      (documentation 'pair-p 'function)
      "Default predicate for struct pair."
      (documentation 'copy-pair 'function)
      "Default copier function for struct pair."
      (documentation 'pair-score 'function)
      "Read the pair score of type fixnum."
      (documentation 'pair-aln1 'function)
      "Read the first fragment of type sam-alignment from a pair."
      (documentation 'pair-aln2 'function)
      "Read the second fragment of type sam-alignment from a pair.")

(declaim (inline pair-rg pair-refid1 pair-pos1 pair-reversed1-p pair-refid2 pair-pos2 pair-reversed2-p))

(defun pair-rg (pair)
  "Access the read group of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-rg (the sam-alignment (pair-aln1 pair))))

(defun pair-refid1 (pair)
  "Access the first REFID of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-refid (the sam-alignment (pair-aln1 pair))))

(defun pair-pos1 (pair)
  "Access the first unclipped position of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-adapted-pos (the sam-alignment (pair-aln1 pair))))

(defun pair-reversed1-p (pair)
  "Access the first direction of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-reversed-p (the sam-alignment (pair-aln1 pair))))

(defun pair-refid2 (pair)
  "Access the second REFID of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-refid (the sam-alignment (pair-aln2 pair))))

(defun pair-pos2 (pair)
  "Access the second unclipped position of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-adapted-pos (the sam-alignment (pair-aln2 pair))))

(defun pair-reversed2-p (pair)
  "Access the second direction of a pair."
  (declare (pair pair) #.*optimization*)
  (sam-alignment-reversed-p (the sam-alignment (pair-aln2 pair))))

(defun handle-pair= (p1 p2)
  "Do the two handles refer to pair instances at the same unclipped positions and with the same directions?"
  (declare (handle p1 p2) #.*optimization*)
  (let ((p1 (handle-object p1))
        (p2 (handle-object p2)))
    (declare (pair p1 p2))
    (and (eq (pair-rg p1) (pair-rg p2))
         (=  (the int32 (pair-refid1 p1)) (the int32 (pair-refid1 p2)))
         (=  (the int32 (pair-pos1 p1)) (the int32 (pair-pos1 p2)))
         (eq (pair-reversed1-p p1) (pair-reversed1-p p2))
         (=  (the int32 (pair-refid2 p1)) (the int32 (pair-refid2 p2)))
         (=  (the int32 (pair-pos2 p1)) (the int32 (pair-pos2 p2)))
         (eq (pair-reversed2-p p1) (pair-reversed2-p p2)))))

(defun pair-hash (p)
  "Hash function that corresponds to handle-pair=, but operates on a pair directly."
  (declare (pair p) #.*optimization*)
  (logxor (sxhash (the base-string (pair-rg p)))
          (sxhash (the int32 (pair-refid1 p)))
          (sxhash (the int32 (pair-refid2 p)))
          (sxhash (+ (ash (the int32 (pair-pos1 p)) 32) (the int32 (pair-pos2 p))))
          (sxhash (the boolean (pair-reversed1-p p)))
          (sxhash (the boolean (pair-reversed2-p p)))))

(defun classify-pair (aln fragments pairs deterministic)
  "For each list of pairs with the same unclipped positions and directions, all except the one with the highest score are marked as duplicates."
  (when (true-pair-p aln)
    (let ((aln1 aln)
          (aln2 (let ((split (hash-table-split aln fragments)))
                  (with-hash-table-locked split
                    (let ((value (gethash aln split)))
                      (cond (value (remhash aln split) value)
                            (t (setf (gethash aln split) aln)
                               (return-from classify-pair))))))))
      (when (> (the int32 (sam-alignment-adapted-pos aln1))
               (the int32 (sam-alignment-adapted-pos aln2)))
        (rotatef aln1 aln2))
      (let* ((score (+ (sam-alignment-adapted-score aln1)
                       (sam-alignment-adapted-score aln2)))
             (keypair (make-pair score aln1 aln2))
             (hash (pair-hash keypair))
             (pairkey (make-handle keypair hash))
             (entry nil)
             (best (let ((split (hash-table-split pairkey pairs)))
                     (or (gethash pairkey split)
                         (progn
                           (setq entry (make-pair score aln1 aln2))
                           (let ((handle (make-handle entry hash)))
                             (with-hash-table-locked split
                               (let ((value (gethash pairkey split)))
                                 (cond (value value)
                                       (t (setf (gethash handle split) handle)
                                          (return-from classify-pair)))))))))))
        (declare (dynamic-extent keypair pairkey))
        (loop for best-pair = (handle-object best)
              for best-pair-score = (pair-score best-pair)
              until (cond ((> best-pair-score score)
                           (mark-as-duplicate aln1)
                           (mark-as-duplicate aln2))
                          ((= best-pair-score score)
                           (cond (deterministic ; code for correctness checks
                                  (cond ((string> (sam-alignment-qname aln1) (sam-alignment-qname (pair-aln1 best-pair)))
                                         (mark-as-duplicate aln1)
                                         (mark-as-duplicate aln2))
                                        ((compare-and-swap (handle-object best) best-pair
                                                           (or entry (setq entry (make-pair score aln1 aln2))))
                                         (mark-as-duplicate (pair-aln1 best-pair))
                                         (mark-as-duplicate (pair-aln2 best-pair)))))
                                 (t (mark-as-duplicate aln1)
                                    (mark-as-duplicate aln2))))
                          ((compare-and-swap (handle-object best) best-pair
                                             (or entry (setq entry (make-pair score aln1 aln2))))
                           (mark-as-duplicate (pair-aln1 best-pair))
                           (mark-as-duplicate (pair-aln2 best-pair)))))))))

(defun mark-duplicates (deterministic)
  "A filter for marking duplicate sam-alignment instances. Depends on the add-refid filter being called before to fill in the refid."
  (lambda (header)
    (declare (ignore header))
    (let ((splits (* 16 *number-of-threads*)))
      ; set up tables once header is parsed, tables will serve for all alignments
      (let ((fragments (make-split-hash-table splits :test #'handle-fragment= :hash-function #'handle-hash))
            (pairs-fragments (make-split-hash-table splits :test #'sam-alignment-pair= :hash-function #'sam-alignment-pair-hash))
            (pairs (make-split-hash-table splits :test #'handle-pair= :hash-function #'handle-hash))
            (rg-table (make-synchronized-hash-table :test #'string= :hash-function #'sxhash)))
        (lambda ()
          (lambda (alignment)
            (when (sam-alignment-flag-notany alignment #.(+ +unmapped+ +secondary+ +duplicate+ +supplementary+))
              ; mark duplicate checks
              (adapt-alignment rg-table alignment)
              (classify-fragment alignment fragments deterministic)
              (classify-pair alignment pairs-fragments pairs deterministic))
            t))))))
