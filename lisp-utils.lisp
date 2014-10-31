(in-package :elprep)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defparameter *optimization*
    '(optimize (speed 3) (space 0) (debug 1) (safety 0)
               (compilation-speed 0))
    "Standard optimization settings without fixnum optimizations.")

  (defparameter *fixnum-optimization*
    '(optimize (speed 3) (space 0) (debug 1) (safety 0) 
               (compilation-speed 0) (hcl:fixnum-safety 0))
    "Standard optimizations settings with fixnum optimizations."))

(hcl:defglobal-variable *keyword* (find-package :keyword)
  "The :keyword package.")

(declaim (inline intern-key))

(defun intern-key (string)
  "Intern a string in the :keyword package."
  (declare (simple-base-string string) #.*optimization*)
  (intern string *keyword*))

(hcl:defglobal-variable *unique-value* (gensym)
  "A unique value for use in the functions presentp and unique.")

(declaim (inline presentp))

(defun presentp (indicator plist)
  "Is the indicator present in the property list?"
  (declare (list plist) #.*optimization*)
  (let ((unique-value *unique-value*))
    (not (eq (getf plist indicator unique-value) unique-value))))

(declaim (inline unique))

(defun unique (indicator plist)
  "Assert that the indicator is not present in the property list."
  (assert (not (presentp indicator plist)))
  indicator)

(defmacro with-modify-hash ((key value found) (hash-table form) &body body)
  "Macro version of LispWorks's modify-hash function."
  `(hcl:modify-hash ,hash-table ,form (lambda (,key ,value ,found)
                                        (declare (ignorable ,key ,value ,found))
                                        (block nil ,@body))))

(defmacro unwind-protectn (&body forms)
  "Like unwind-protect, except that all forms but the last are protected, and only the last form is used for cleanup"
  `(unwind-protect (progn ,@(butlast forms)) ,@(last forms)))

(defun compose-thunks (thunks)
  "Return a single thunk that executes the given thunks in sequence."
  (declare (list thunks))
  (cond ((null thunks) (constantly nil))
        ((null (cdr thunks)) (car thunks))
        (t (lambda ()
             (declare #.*optimization*)
             (loop for fun in thunks do (funcall fun))))))

(declaim (inline mapfiltermap))

(defun mapfiltermap (inmap filter outmap list &optional tail)
  "Apply the following steps to each element in the list, optionally bounded by tail:
   - Apply inmap.
   - Filter out each element for which filter returns nil.
   - Apply outmap."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (locally (declare (cons list))
      (loop for car = (funcall inmap (car list))
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            when filtered collect (funcall outmap car)
            until (eq cdr tail)
            do (setq list cdr)))))

(declaim (inline nmapfiltermap))

(defun nmapfiltermap (inmap filter outmap list &optional tail)
  "Destructively apply the following steps to each element in the list, optionally bounded by tail:
   - Apply inmap.
   - Filter out each element for which filter returns nil.
   - Apply outmap."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (let ((head (cons nil list)))
      (declare (cons list head) (dynamic-extent head))
      (loop with prev of-type cons = head
            for car = (funcall inmap (car list))
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            until (eq cdr tail) do
            (locally (declare (cons cdr))
              (if filtered
                (setf (car list) (funcall outmap car) prev list list cdr)
                (setf (car list) (car cdr) (cdr list) (cdr cdr))))
            finally (if filtered
                      (setf (car list) (funcall outmap car))
                      (setf (cdr prev) tail)))
      (cdr head))))

(declaim (inline mapfilter))

(defun mapfilter (map filter list &optional tail)
  "Apply the following steps to each element in the list, optionally bounded by tail:
   - Apply map.
   - Filter out each element for which filter returns nil."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (locally (declare (cons list))
      (loop for car = (funcall map (car list))
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            when filtered collect car
            until (eq cdr tail)
            do (setq list cdr)))))

(declaim (inline nmapfilter))

(defun nmapfilter (map filter list &optional tail)
  "Destructively apply the following steps to each element in the list, optionally bounded by tail:
   - Apply map.
   - Filter out each element for which filter returns nil."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (let ((head (cons nil list)))
      (declare (cons list head) (dynamic-extent head))
      (loop with prev of-type cons = head
            for car = (funcall map (car list))
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            until (eq cdr tail) do
            (locally (declare (cons cdr))
              (if filtered
                (setf (car list) car prev list list cdr)
                (setf (car list) (car cdr) (cdr list) (cdr cdr))))
            finally (if filtered
                      (setf (car list) car)
                      (setf (cdr prev) tail)))
      (cdr head))))

(declaim (inline filtermap))

(defun filtermap (filter map list &optional tail)
  "Apply the following steps to each element in the list, optionally bounded by tail:
   - Filter out each element for which filter returns nil.
   - Apply map."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (locally (declare (cons list))
      (loop for car = (car list)
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            when filtered collect (funcall map car)
            until (eq cdr tail)
            do (setq list cdr)))))

(declaim (inline nfiltermap))

(defun nfiltermap (filter map list &optional tail)
  "Destructively apply the following steps to each element in the list, optionally bounded by tail:
   - Filter out each element for which filter returns nil.
   - Apply map."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (let ((head (cons nil list)))
      (declare (cons list head) (dynamic-extent head))
      (loop with prev of-type cons = head
            for car = (car list)
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            until (eq cdr tail) do
            (locally (declare (cons cdr))
              (if filtered
                (setf (car list) (funcall map car) prev list list cdr)
                (setf (car list) (car cdr) (cdr list) (cdr cdr))))
            finally (if filtered
                      (setf (car list) (funcall map car))
                      (setf (cdr prev) tail)))
      (cdr head))))

(declaim (inline filter))

(defun filter (filter list &optional tail)
  "Filter out each element from the list, optionally bounded by tail, for which filter returns nil."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (locally (declare (cons list))
      (loop for car = (car list)
            for cdr = (cdr list)
            when (funcall filter car) collect car
            until (eq cdr tail)
            do (setq list cdr)))))

(declaim (inline nfilter))

(defun nfilter (filter list &optional tail)
  "Destructively filter out each element from the list, optionally bounded by tail, for which filter returns nil."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (let ((head (cons nil list)))
      (declare (cons list head) (dynamic-extent head))
      (loop with prev of-type cons = head
            for car = (car list)
            for cdr = (cdr list)
            for filtered = (funcall filter car)
            until (eq cdr tail) do
            (locally (declare (cons cdr))
              (if filtered
                (setf prev list list cdr)
                (setf (car list) (car cdr) (cdr list) (cdr cdr))))
            finally (if filtered
                      ()
                      (setf (cdr prev) tail)))
      (cdr head))))

(declaim (inline mapcar*))

(defun mapcar* (map list &optional tail)
  "Like mapcar, except operates on only one list, optionally bounded by tail."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (locally (declare (cons list))
      (loop for car = (car list)
            for cdr = (cdr list)
            collect (funcall map car)
            until (eq cdr tail)
            do (setq list cdr)))))

(declaim (inline nmapcar*))

(defun nmapcar* (map list &optional tail)
  "Like mapcar, except destructively operates on only one list, optionally bounded by tail."
  (declare #.*optimization*)
  (if (eq list tail) tail
    (loop with cur of-type cons = list
          for car = (car cur)
          for cdr = (cdr cur)
          do (setf (car cur) (funcall map car))
          until (eq cdr tail)
          do (setq cur cdr)
          finally (return list))))

(defun nthdiff (n list)
  "Return a copy of the first n elements of list, and the nth cdr of the list."
  (declare (fixnum n) (list list) #.*optimization*)
  (loop for tail on list repeat n
        collect (car tail) into head
        finally (return (values head tail))))

;;; split hash tables

(defstruct (split-hash-table (:constructor %make-split-hash-table (hash-function vector)))
  "A collection of hash tables distributing the entries based on a permutation of the same hash function used for the splits.
   The struct split-hash-table has a constructor %make-split-hash-table that takes a hash function and a vector of splits as parameters.
   Read-only accessor split-hash-table-hash-function refers to the hash function.
   Read-only accessor split-hash-table-vector of type simple-vector refers to the splits.
   Primary use of this struct is to allow for locking splits separately to avoid lock contention when operating on a hash table in parallel."
  (hash-function nil :read-only t)
  (vector #() :type simple-vector :read-only t))

(setf (documentation '%make-split-hash-table 'function)
      "Constructor for struct split-hash-table that takes a hash function and a vector of splits as parameters."
      (documentation 'split-hash-table-p 'function)
      "Default predicate for struct split-hash-table."
      (documentation 'copy-split-hash-table 'function)
      "Default copier function for struct split-hash-table."
      (documentation 'split-hash-table-hash-function 'function)
      "Read the split-hash-table hash function."
      (documentation 'split-hash-table-vector 'function)
      "Read the split-hash-table vector of splits of type simple-vector.")

(defun make-split-hash-table (splits &rest args
                                     &key
                                     test size rehash-size rehash-threshold
                                     (hash-function (error "No hash function passed to make-split-hash-table."))
                                     weak-kind single-thread free-function)
  "Constructor for split-hash-table that takes the number of splits and initialization arguments as for make-hash-table as parameters.
   The :hash-function initialization argument must be explicitly provided."
  (declare (dynamic-extent args) (ignore test size rehash-size rehash-threshold weak-kind single-thread free-function))
  (loop with vector = (make-array splits :single-thread t)
        for i below splits do
        (setf (svref vector i) (apply 'make-hash-table args))
        finally (return (%make-split-hash-table hash-function vector))))

(defconstant +lowest-15-bits+  (1- (ash 1 15))
  "Bit mask for the lowest 15 bits of a fixnum, used for split hash tables.")
(defconstant +highest-45-bits+ (ash (1- (ash 1 45)) 15)
  "Bit mask for the highest 45 bits of a fixnum, used for split hash tables.")

(declaim (inline rotate-45))

(defun rotate-45 (n)
  "Rotate a fixnum by 45 bits, used for split hash tables."
  (declare (fixnum n) #.*fixnum-optimization*)
  (logior (ash (logand +lowest-15-bits+ n) 45)
          (ash (logand +highest-45-bits+ n) -15)))

(declaim (inline hash-table-split))

(defun hash-table-split (key table)
  "Return a split of a split-hash-table for a given key."
  (declare (split-hash-table table) #.*fixnum-optimization*)
  (let ((hash-function (split-hash-table-hash-function table))
        (vector (split-hash-table-vector table)))
    (declare (vector vector))
    (svref vector (rem (rotate-45 (funcall hash-function key)) (length vector)))))

(defgeneric copy-stream (input output)
  (:documentation "Efficient copying of the contents of an input stream to an output stream.")
  (:method ((input buffered-stream) (output stream))
   "Specialization for LispWorks's buffered-stream."
   (loop do (with-stream-input-buffer (buffer index limit) input
              (when (< index limit)
                (stream-write-sequence output buffer index limit)))
         while (stream-fill-buffer input))))

(defun elprep-debugger-hook (condition hook)
  (declare (ignore hook))
  (dbg:log-bug-form (format nil "An error occurred in elPrep: ~A" condition) :message-stream t)
  (lw:quit :status -1 :confirm nil :ignore-errors-p t))

(defun process-run (name function &rest arguments)
  (declare (dynamic-extent arguments))
  (apply #'mp:process-run-function name '() function arguments))
