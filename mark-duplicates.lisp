(in-package :elprep)

(declaim (inline sam-alignment-rg (setf sam-alignment-rg)))

(defun sam-alignment-rg (aln)
  "Access the read group optional field of a sam-alignment."
  (sam-alignment-tag aln :rg))

(defun (setf sam-alignment-rg) (new-value aln)
  "Access the read group optional field of a sam-alignment."
  (setf (sam-alignment-tag aln :rg) new-value))

(defun make-phred-score-table ()
  "Map Phred qualities to a reasonable range and an error flag indicating if it is outside a valid range."
  (let ((score-table (make-array 512 :element-type 'base-char :allocation :long-lived :single-thread t)))
    (loop for char below 256
          for pos = (ash char 1) do
          (if (or (< char 33) (> char 126))
            (setf (sbchar score-table pos)      #.(code-char 0)
                  (sbchar score-table (1+ pos)) #.(code-char 1))
            (let ((qual (- char 33)))
              (setf (sbchar score-table pos)      (if (>= qual 15) (code-char qual) #.(code-char 0))
                    (sbchar score-table (1+ pos)) #.(code-char 0)))))
    score-table))

(declaim (notinline compute-phred-score))

(defun compute-phred-score (aln)
  "Sum the adapted Phred qualities of a sam-alignment."
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (let ((string (sam-alignment-qual aln))
        (score-table (load-time-value (make-phred-score-table) t))
        (score 0) (error 0))
    (declare (simple-base-string string) (fixnum score error))
    (loop for i of-type fixnum below (length string)
          for pos of-type fixnum = (ash (char-code (sbchar string i)) 1) do
          (setq score (+ score (char-code (sbchar score-table pos)))
                error (logior error (char-code (sbchar score-table (1+ pos))))))
    (assert (= error 0))
    score))

(define-symbol-macro upcase-cigar-operations "MIDNSHPX=")

(defconstant +min-upcase-cigar-operation+ (reduce 'min upcase-cigar-operations :key 'char-code)
  "The smallest CIGAR operation, upper case only.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(defconstant +max-upcase-cigar-operation+ (reduce 'max upcase-cigar-operations :key 'char-code)
  "The largest CIGAR operation, upper case only.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(declaim (inline cigar-aux-pos))

(defun cigar-aux-pos (char)
  "Position in a vector of CIGAR operations, starting with the smallest CIGAR operation, upper case only."
  (declare (base-char char) #.*fixnum-optimization*)
  (- (char-code char) +min-upcase-cigar-operation+))

(defun make-unclipped-aux-tables ()
  "Map CIGAR operations to flags indicating whether they are clipped and/or reference operations, upper case only."
  (let* ((tablesize (1+ (- +max-upcase-cigar-operation+
                           +min-upcase-cigar-operation+)))
         (clipped   (make-array tablesize
                                :initial-element #.(code-char 0)
                                :element-type 'base-char
                                :allocation :long-lived
                                :single-thread t))
         (reference (make-array tablesize
                                :initial-element #.(code-char 0)
                                :element-type 'base-char
                                :allocation :long-lived
                                :single-thread t)))
    (flet ((set-clipped-flag (char)
             (setf (sbchar clipped (cigar-aux-pos char)) #.(code-char 1)))
           (set-reference-flag (char)
             (setf (sbchar reference (cigar-aux-pos char)) #.(code-char 1))))
      (set-clipped-flag #\S)
      (set-clipped-flag #\H)
      (set-reference-flag #\M)
      (set-reference-flag #\D)
      (set-reference-flag #\N)
      (set-reference-flag #\=)
      (set-reference-flag #\X))
    (cons clipped reference)))

(declaim (notinline compute-unclipped-position))

(defun compute-unclipped-position (aln)
  "Compute unclipped position of a sam-alignment, based on its POS and CIGAR string."
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (let ((cigar (scan-cigar-string 'vector (sam-alignment-cigar aln))))
    (declare (simple-vector cigar))
    (if (= 0 (length cigar))
      (sam-alignment-pos aln)
      (let ((tables (load-time-value (make-unclipped-aux-tables) t)))
        (declare (cons tables))
        (if (sam-alignment-reversed-p aln)
          (1- (+ (sam-alignment-pos aln)
                 (loop for i of-type fixnum = (1- (length cigar)) then (1- i)
                       for (key . value) of-type (symbol . fixnum) = (svref cigar i)
                       for p of-type fixnum = (cigar-aux-pos (sbchar (symbol-name key) 0))
                       for c of-type fixnum = (char-code (sbchar (car tables) p))
                       for r of-type fixnum = (char-code (sbchar (cdr tables) p)) 
                       for clipped of-type fixnum = c then (* c clipped)
                       sum (* (logior r clipped) value) fixnum
                       until (= i 0))))
          (- (sam-alignment-pos aln)
             (loop for i of-type fixnum below (length cigar)
                   for (key . value) of-type (symbol . fixnum) = (svref cigar i)
                   for p of-type fixnum = (cigar-aux-pos (sbchar (symbol-name key) 0))
                   until (= 0 (char-code (sbchar (car tables) p)))
                   sum value fixnum)))))))

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
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (setf (sam-alignment-flag aln) (logior (sam-alignment-flag aln) +duplicate+)))

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
    (and (eq (sam-alignment-rg          f1) (sam-alignment-rg          f2))
         (=  (sam-alignment-refid       f1) (sam-alignment-refid       f2))
         (=  (sam-alignment-adapted-pos f1) (sam-alignment-adapted-pos f2))
         (eq (sam-alignment-reversed-p  f1) (sam-alignment-reversed-p  f2)))))

(defun fragment-hash (f)
  "Hash function that corresponds to handle-fragment=, but operates an a sam-alignment directly."
  (declare (sam-alignment f) #.*optimization*)
  (logxor (sxhash (sam-alignment-rg f))
          (sxhash (sam-alignment-refid f))
          (sxhash (sam-alignment-adapted-pos f))
          (sxhash (sam-alignment-reversed-p f))))

(declaim (inline true-fragment-p true-pair-p))

(defun true-fragment-p (aln)
  "Is this sam-alignment definitely not part of a pair?"
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (/= (logand (sam-alignment-flag aln) #.(logior +multiple+ +next-unmapped+)) +multiple+))

(defun true-pair-p (aln)
  "Is this sam-alignment definitely part of a pair?"
  (declare (sam-alignment aln) #.*fixnum-optimization*)
  (= (logand (sam-alignment-flag aln) #.(logior +multiple+ +next-unmapped+)) +multiple+))

(defun classify-fragment (aln fragments) ; fragments is the hash table to store intermediate "best" alns, passed as a parameter
  "For each list of sam-alignment instances with the same unclipped position and direction, all except the one with the highest score are marked as duplicates.
   If there are fragments in such a list that are actually part of pairs, all the true fragments are marked as duplicates and the pairs are left untouched."
  (let* ((hash (fragment-hash aln))
         (key (make-handle aln hash))
         (split (hash-table-split key fragments))
         (best (or (gethash key split)
                   (let ((entry (make-handle aln hash)))
                     (hcl:with-hash-table-locked split
                       (let ((value (gethash key split)))
                         (cond (value value)
                               (t (setf (gethash entry split) entry)
                                  (return-from classify-fragment)))))))))
    (declare (dynamic-extent key))
    (cond ((true-fragment-p aln)
           (loop with aln-score = (sam-alignment-adapted-score aln)
                 for best-aln = (handle-object best)
                 until (cond ((or (true-pair-p best-aln)
                                  (>= (sam-alignment-adapted-score best-aln) aln-score))
                              (mark-as-duplicate aln) t)
                             ((sys:compare-and-swap (handle-object best) best-aln aln)
                              (mark-as-duplicate best-aln) t))))
          (t ; the aln is a true pair object, there may be a true fragment stored in the hash table which we then need to mark and swap out
           (loop for best-aln = (handle-object best)
                 until (cond ((true-pair-p best-aln) t) ; stop, the best in the hash tab is a pair, this is marked via mark-duplicates
                             ((sys:compare-and-swap (handle-object best) best-aln aln)
                              (mark-as-duplicate best-aln) t)))))))

(defun sam-alignment-pair= (aln1 aln2)
  "Are the two sam-alignment fragments part of the same pair?"
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  (and (eq      (sam-alignment-rg    aln1) (sam-alignment-rg    aln2))
       (string= (sam-alignment-qname aln1) (sam-alignment-qname aln2))))

(defun sam-alignment-pair-hash (aln)
  "Hash function that corresponds to sam-alignment-pair=."
  (declare (sam-alignment aln) #.*optimization*)
  (or (sam-alignment-temp aln :pair-hash)
      (let ((hash (logxor (sxhash (sam-alignment-rg aln))
                          (sxhash (sam-alignment-qname aln)))))
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
    (and (eq (pair-rg          p1) (pair-rg          p2))
         (=  (pair-refid1      p1) (pair-refid1      p2))
         (=  (pair-pos1        p1) (pair-pos1        p2))
         (eq (pair-reversed1-p p1) (pair-reversed1-p p2))
         (=  (pair-refid2      p1) (pair-refid2      p2)) 
         (=  (pair-pos2        p1) (pair-pos2        p2))
         (eq (pair-reversed2-p p1) (pair-reversed2-p p2)))))

(defun pair-hash (p)
  "Hash function that corresponds to handle-pair=, but operates on a pair directly."
  (declare (pair p) #.*optimization*)
  (logxor (sxhash (pair-rg p))
          (sxhash (pair-refid1 p))
          (sxhash (pair-refid2 p))
          (sxhash (+ (ash (pair-pos1 p) 32) (pair-pos2 p)))
          (sxhash (pair-reversed1-p p))
          (sxhash (pair-reversed2-p p))))

(defun classify-pair (aln fragments pairs)
  "For each list of pairs with the same unclipped positions and directions, all except the one with the highest score are marked as duplicates."
  (when (true-pair-p aln)
    (let ((aln1 aln)
          (aln2 (let ((split (hash-table-split aln fragments)))
                  (hcl:with-hash-table-locked split
                    (let ((value (gethash aln split)))
                      (cond (value (remhash aln split) value)
                            (t (setf (gethash aln split) aln)
                               (return-from classify-pair))))))))
      (when (> (sam-alignment-adapted-pos aln1)
               (sam-alignment-adapted-pos aln2))
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
                             (hcl:with-hash-table-locked split
                               (let ((value (gethash pairkey split)))
                                 (cond (value value)
                                       (t (setf (gethash handle split) handle)
                                          (return-from classify-pair)))))))))))
        (declare (dynamic-extent keypair pairkey))
        (loop for best-pair = (handle-object best)
              until (cond ((>= (pair-score best-pair) score)
                           (mark-as-duplicate aln1)
                           (mark-as-duplicate aln2) t)
                          ((sys:compare-and-swap (handle-object best) best-pair
                                                 (or entry (setq entry (make-pair score aln1 aln2))))
                           (mark-as-duplicate (pair-aln1 best-pair))
                           (mark-as-duplicate (pair-aln2 best-pair)) t)))))))

(defun mark-duplicates (header)
  "A filter for marking duplicate sam-alignment instances. Depends on the add-refid filter being called before to fill in the refid."
  (declare (ignore header))
  (let ((splits (* 8 *number-of-threads*)))
    ; set up tables once header is parsed, tables will serve for all alignments
    (let ((fragments (make-split-hash-table splits :test 'handle-fragment= :hash-function 'handle-hash))
          (pairs-fragments (make-split-hash-table splits :test 'sam-alignment-pair= :hash-function 'sam-alignment-pair-hash))
          (pairs (make-split-hash-table splits :test 'handle-pair= :hash-function 'handle-hash))
          (rg-table (make-hash-table :test 'string= :hash-function 'sxhash)))
      (lambda ()
        (lambda (alignment)
          (when (sam-alignment-flag-notany alignment #.(+ +unmapped+ +secondary+ +duplicate+ +supplementary+))
            ; mark duplicate checks
            (adapt-alignment rg-table alignment)
            (classify-fragment alignment fragments)
            (classify-pair alignment pairs-fragments pairs))
          t)))))
