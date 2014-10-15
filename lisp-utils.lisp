(in-package :elprep)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defparameter *optimization*
    '(optimize (speed 3) (space 0) (debug 1) (safety 0)
               (compilation-speed 0))
    "Standard optimization settings without fixnum optimizations.")

  (defparameter *fixnum-optimization*
    '(optimize (speed 3) (space 0) (debug 1) (safety 0)
               (compilation-speed 0) #+lispworks (hcl:fixnum-safety 0))
    "Standard optimizations settings with fixnum optimizations."))

;;; portability

(defmacro sbs (string)
  "Coerce (literal) string to 'simple-base-string at macro-expansion time."
  #+lispworks string
  #+sbcl (coerce string 'simple-base-string))

#+lispworks (define-symbol-macro empty-sbs "")
#+sbcl (sb-ext:defglobal empty-sbs (sbs "") "An empty simple-base-string.")

(defmacro defglobal (var value &optional (doc nil docp))
  "Define a global variable."
  #+lispworks
  (if docp
    `(hcl:defglobal-variable ,var ,value ,doc)
    `(hcl:defglobal-variable ,var ,value))
  #+sbcl
  (if docp
    `(sb-ext:defglobal ,var ,value ,doc)
    `(sb-ext:defglobal ,var ,value)))

(declaim (inline make-single-thread-hash-table make-synchronized-hash-table))

(defun make-single-thread-hash-table (&rest args &key test size rehash-size rehash-treshold hash-function)
  "Like make-hash-table, but ensure it is single-thread, not synchronized."
  (declare (dynamic-extent args) (ignore test size rehash-size rehash-treshold hash-function))
  #+lispworks (apply #'cl:make-hash-table :single-thread t args)
  #+sbcl (apply #'cl:make-hash-table :synchronized nil args))

(defun make-synchronized-hash-table (&rest args &key test size rehash-size rehash-treshold hash-function)
  "Like make-hash-table, but ensure it is synchronized, not single-thread."
  (declare (dynamic-extent args) (ignore test size rehash-size rehash-treshold hash-function))
  #+lispworks (apply #'cl:make-hash-table :single-thread nil args)
  #+sbcl (apply #'cl:make-hash-table :synchronized t args))

#+sbcl
(defun modify-hash (hash-table key function)
  "Similar to LispWorks's hcl:modify-hash."
  (if (sb-ext:hash-table-synchronized-p hash-table)
    (sb-ext:with-locked-hash-table (hash-table)
      (multiple-value-bind (value foundp)
          (gethash key hash-table)
        (values (setf (gethash key hash-table)
                      (locally (declare #.*optimization*)
                        (funcall (the function function) key value foundp)))
                key)))
    (multiple-value-bind (value foundp)
        (gethash key hash-table)
      (values (setf (gethash key hash-table)
                    (locally (declare #.*optimization*)
                      (funcall (the function function) key value foundp)))
              key))))

#+sbcl
(defmacro with-hash-table-locked (hash-table &body body)
  "Renamed sb-ext:with-locked-hash-table."
  `(sb-ext:with-locked-hash-table (,hash-table) ,@body))

(declaim (inline open-program close-program program-input program-output))

#+lispworks
(progn
  (defun open-program (command &rest args &key (direction :input))
    "Similar to LispWorks's sys:open-pipe, to enable portability between LispWorks and SBCL."
    (declare (dynamic-extent args) (ignore direction))
    (apply #'sys:open-pipe command args))

  (defun close-program (program)
    "Close a program opened with open-program."
    (close program))

  (defun program-input (program)
    "Get the stream of a prgoram that expects input."
    program)

  (defun program-output (program)
    "Get the stream of a program that produces output."
    program))

#+sbcl
(progn
  (defun open-program (command &key (direction :input) (external-format :default))
    "Similar to LispWorks's sys:open-pipe, to enable portability between LispWorks and SBCL."
    (etypecase command
      (string (sb-ext:run-program "/bin/sh" (list "-c" command)
                                  :wait nil (ecase direction (:input :output) (:output :input)) :stream
                                  :external-format external-format))
      (list   (sb-ext:run-program (first command) (rest command)
                                  :wait nil (ecase direction (:input :output) (:output :input)) :stream
                                  :external-format external-format))))

  (defun close-program (program)
    "Close a program opened with open-program."
    (sb-ext:process-close program))

  (defun program-input (program)
    "Get the stream of a program that expects input."
    (sb-ext:process-input program))

  (defun program-output (program)
    "Get the stream of a program that produces output."
    (sb-ext:process-output program)))

(defmacro with-open-program ((program command &rest args) &body body)
  "Open a program for a block of code, and ensure that it is closed again in the end."
  `(let ((,program (open-program ,command ,@args)))
     (unwind-protect (progn ,@body)
       (close-program ,program))))

(declaim (inline make-buffered-stream))

#+lispworks
(progn
  (defun make-buffered-stream (stream)
    "Create a wrapper around a stream to enable buffers that can be explicitly accessed.
     This is a no-op in LispWorks, but creates a separater buffer in SBCL."
    (etypecase stream (buffered-stream stream)))

  (declaim (inline buffered-listen buffered-peek buffered-read-line))

  (defun buffered-listen (stream)
    "Like listen, but on buffered-stream."
    (declare (buffered-stream stream) #.*optimization*)
    (loop (stream:with-stream-input-buffer (buffer index limit) stream
            (declare (ignore buffer) (fixnum index limit))
            (when (< index limit)
              (return-from buffered-listen t)))
          (unless (stream:stream-fill-buffer stream)
            (return-from buffered-listen nil))))

  (defun buffered-peek (stream)
    "Like (peek-char nil ... nil #\EOT), but on buffered-stream."
    (declare (buffered-stream stream) #.*optimization*)
    (loop (stream:with-stream-input-buffer (buffer index limit) stream
            (declare (simple-base-string buffer) (fixnum index limit))
            (when (< index limit)
              (return-from buffered-peek (lw:sbchar buffer index))))
          (unless (stream:stream-fill-buffer stream)
            (return-from buffered-peek #\EOT))))

  (defun buffered-read-line (stream)
    "Like read-line, but on buffered-stream."
    (stream:stream-read-line stream)))

#+sbcl
(progn
  (defconstant +buffered-stream-limit+ 8192)

  (deftype stream-buffer ()
    '(simple-array character (8192)))

  (defstruct (buffered-stream (:constructor make-buffered-stream (stream))
                              (:copier nil)
                              (:predicate nil))
    (buffer (make-array +buffered-stream-limit+ :element-type 'character) :type stream-buffer :read-only t)
    (index +buffered-stream-limit+ :type fixnum)
    (limit +buffered-stream-limit+ :type fixnum)
    (stream nil :type stream))

  (setf (documentation 'make-buffered-stream 'function)
        "Create a wrapper around a stream to enable buffers that can be explicitly accessed.
         This is a no-op in LispWorks, but creates a separater buffer in SBCL.")

  (defun buffered-listen-slow-path (stream)
    "Like listen, but on buffered-stream, slow path."
    (declare (buffered-stream stream) #.*optimization*)
    (let ((position (read-sequence
                     (buffered-stream-buffer stream)
                     (buffered-stream-stream stream))))
      (declare (fixnum position))
      (setf (buffered-stream-index stream) 0
            (buffered-stream-limit stream) position)
      (> position 0)))

  (declaim (inline buffered-listen))

  (defun buffered-listen (stream)
    "Like listen, but on buffered-stream."
    (declare (buffered-stream stream) #.*optimization*)
    (let ((index (buffered-stream-index stream))
          (limit (buffered-stream-limit stream)))
      (declare (fixnum index limit))
      (or (< index limit)
          (buffered-listen-slow-path stream))))

  (defun buffered-peek-slow-path (stream buffer)
    "Like (peek-char nil ... nil #\EOT), but on buffered-stream, slow path."
    (declare (buffered-stream stream) (simple-base-string buffer) #.*optimization*)
    (let ((position (read-sequence
                     (buffered-stream-buffer stream)
                     (buffered-stream-stream stream))))
      (declare (fixnum position))
      (setf (buffered-stream-index stream) 0
            (buffered-stream-limit stream) position)
      (if (= position 0) #\EOT (schar buffer 0))))

  (declaim (inline buffered-peek))

  (defun buffered-peek (stream)
    "Like (peek-char nil ... nil #\EOT), but on buffered-stream."
    (declare (buffered-stream stream) #.*optimization*)
    (let ((buffer (buffered-stream-buffer stream))
          (index (buffered-stream-index stream))
          (limit (buffered-stream-limit stream)))
      (declare (stream-buffer buffer) (fixnum index limit))
      (if (< index limit)
        (schar buffer index)
        (buffered-peek-slow-path stream buffer))))

  (defun buffered-read-line (stream)
    "Like read-line, but on buffered-stream."
    (declare (buffered-stream stream) #.*optimization*)
    (loop with strings of-type list = '()
          with buffer of-type stream-buffer = (buffered-stream-buffer stream) do
          (loop with index of-type fixnum = (buffered-stream-index stream)
                with limit of-type fixnum = (buffered-stream-limit stream)
                for end of-type fixnum from index below limit
                for flag of-type boolean = (char= (schar buffer end) #\Newline)
                until flag finally
                (let ((string (make-array (- end index) :element-type 'base-char)))
                  (declare (simple-base-string string))
                  (loop for j of-type fixnum from index below end
                        for i of-type fixnum from 0
                        do (setf (schar string i) (schar buffer j)))
                  (cond (flag (setf (buffered-stream-index stream) (1+ end))
                              (cond (strings (push string strings)
                                             (return-from buffered-read-line (values (apply #'concatenate 'simple-base-string (nreverse strings)) nil)))
                                    (t       (return-from buffered-read-line (values string nil)))))
                        (t    (setf (buffered-stream-index stream) end)
                              (push string strings)))))
          (let ((position (read-sequence buffer (buffered-stream-stream stream))))
            (declare (fixnum position))
            (setf (buffered-stream-index stream) 0
                  (buffered-stream-limit stream) position)
            (when (= position 0)
              (if strings
                (if (cdr strings)
                  (return-from buffered-read-line (values (apply #'concatenate 'simple-base-string (nreverse strings)) t))
                  (return-from buffered-read-line (values (car strings) t)))
                (return-from buffered-read-line (values empty-sbs t))))))))

#+sbcl
(progn
  (defglobal *working-directory-name* "" "Cached working directory name.")
  (defglobal *working-directory* nil "Cached working directory.")

  (defun get-working-directory ()
    "Similar to LispWorks's hcl:get-working-directory."
    (let ((cwd (sb-posix:getcwd)))
      (if (string= cwd *working-directory-name*)
        *working-directory*
        (setq *working-directory-name* cwd
              *working-directory*
              (parse-namestring (concatenate 'string (sb-posix:getcwd) "/")))))))

(defvar *number-of-threads* 1
  "The number of threads used by run-pipeline and run-pipeline-in-situ to parallelize the filtering process.
   Default is 1, which results in sequential execution. Usually, parallelization only occurs when this value is greater than 3.")

(defun process-run (name function &rest arguments)
  "Wrapper around mp:process-run-function / sb-thread:make-thread."
  (declare (dynamic-extent arguments))
  #+lispworks (apply #'mp:process-run-function name '() function arguments)
  #+sbcl (sb-thread:make-thread function :name name :arguments (copy-list arguments)))

#+sbcl
(progn
  (declaim (inline make-mailbox mailbox-send mailbox-read process-join))

  (defun make-mailbox ()
    "Similar to LispWorks's mp:make-mailbox."
    (sb-concurrency:make-mailbox))

  (defun mailbox-send (mailbox object)
    "Similar to LispWorks's mp:mailbox-send."
    (sb-concurrency:send-message mailbox object))

  (defun mailbox-read (mailbox)
    "Similar to LispWorks's mp:mailbox-read."
    (sb-concurrency:receive-message mailbox))

  (defun process-join (process)
    "Similar to LispWorks's mp:process-join."
    (sb-thread:join-thread process)))

;;; utilities

(defglobal *keyword* (find-package :keyword)
  "The :keyword package.")

(declaim (inline intern-key))

(defun intern-key (string)
  "Intern a string in the :keyword package."
  (declare (simple-base-string string) #.*optimization*)
  (intern string *keyword*))

(defun intern-key/copy (string)
  "Intern a string in the :keyword package; copy string if not found."
  (declare (simple-base-string string) #.*optimization*)
  (let ((keyword *keyword*))
    (or (find-symbol string keyword)
        (intern (copy-seq string) keyword))))

(defglobal *unique-value* (gensym)
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
  `(modify-hash ,hash-table ,form (lambda (,key ,value ,found)
                                    (declare (ignorable ,key ,value ,found))
                                    (block nil ,@body))))

(defmacro unwind-protectn (&body forms)
  "Like unwind-protect, except that all forms but the last are protected, and only the last form is used for cleanup"
  `(unwind-protect (progn ,@(butlast forms)) ,@(last forms)))

(defun get-function (object)
  "Get a function object from a function (returned as is) or a function name."
  (etypecase object
    (function object)
    (symbol (symbol-function object))
    (cons (fdefinition object))))

(defun compose-thunks (thunks)
  "Return a single thunk that executes the given thunks in sequence."
  (declare (list thunks))
  (cond ((null thunks) (constantly nil))
        ((null (cdr thunks)) (car thunks))
        (t (lambda ()
             (declare #.*optimization*)
             (loop for fun in thunks do (funcall (the function fun)))))))

(declaim (inline mapfiltermap))

(defun mapfiltermap (inmap filter outmap list &optional tail)
  "Apply the following steps to each element in the list, optionally bounded by tail:
   - Apply inmap.
   - Filter out each element for which filter returns nil.
   - Apply outmap."
  (declare (function inmap filter outmap) #.*optimization*)
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
  (declare (function inmap filter outmap) #.*optimization*)
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
  (declare (function map filter) #.*optimization*)
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
  (declare (function map filter) #.*optimization*)
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
  (declare (function filter map) #.*optimization*)
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
  (declare (function filter map) #.*optimization*)
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
  (declare (function filter) #.*optimization*)
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
  (declare (function filter) #.*optimization*)
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
  (declare (function map) #.*optimization*)
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
  (declare (function map) #.*optimization*)
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
  (hash-function nil :type function :read-only t)
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

(defun make-split-hash-table (splits &rest args &key
                                     test size rehash-size rehash-threshold
                                     (hash-function (error "No hash function passed to make-split-hash-table.")))
  "Constructor for split-hash-table that takes the number of splits and initialization arguments as for make-hash-table as parameters.
   The :hash-function initialization argument must be explicitly provided."
  (declare (dynamic-extent args) (ignore size rehash-size rehash-threshold))
  (setq test (get-function test)
        hash-function (get-function hash-function))
  (loop with vector = (make-array splits #+lispworks :single-thread #+lispworks t)
        for i below splits do
        (setf (svref vector i) (apply #'make-synchronized-hash-table :test test :hash-function hash-function args))
        finally (return (%make-split-hash-table hash-function vector))))

(defconstant +total-bits+ #.(integer-length most-positive-fixnum))
(defconstant +low-bits+ 15)
(defconstant +high-bits+ (- +total-bits+ +low-bits+))

(defconstant +lowest-bits+  (1- (ash 1 +low-bits+))
  "Bit mask for the lowest bits of a fixnum, used for split hash tables.")
(defconstant +highest-bits+ (ash (1- (ash 1 +high-bits+)) +low-bits+)
  "Bit mask for the highest bits of a fixnum, used for split hash tables.")

(declaim (inline rotate-15))

(defun rotate-15 (n)
  "Rotate a fixnum by 15 bits, used for split hash tables."
  (declare (fixnum n) #.*optimization*)
  (logior (the fixnum (ash (logand +lowest-bits+ n) +high-bits+))
          (the fixnum (ash (logand +highest-bits+ n) (- +low-bits+)))))

(declaim (inline hash-table-split))

(defun hash-table-split (key table)
  "Return a split of a split-hash-table for a given key."
  (declare (split-hash-table table) #.*optimization*)
  (let ((hash-function (split-hash-table-hash-function table))
        (vector (split-hash-table-vector table)))
    (declare (function hash-function) (vector vector))
    (svref vector (rem (rotate-15 (funcall hash-function key)) (length vector)))))

(defgeneric copy-stream (input output)
  (:documentation "Efficient copying of the contents of an input stream to an output stream.")
  #+lispworks
  (:method ((input buffered-stream) (output stream))
   (declare #.*optimization*)
   (loop do (stream:with-stream-input-buffer (buffer index limit) input
              (declare (simple-base-string buffer) (fixnum index limit))
              (when (< index limit)
                (stream:stream-write-sequence output buffer index limit)))
         while (stream:stream-fill-buffer input)))
  #+sbcl
  (:method ((input buffered-stream) (output stream))
   (declare #.*optimization*)
   (loop with buffer of-type stream-buffer = (buffered-stream-buffer input)
         with stream of-type stream = (buffered-stream-stream input)
         do (let ((index (buffered-stream-index input))
                  (limit (buffered-stream-limit input)))
              (declare (fixnum index limit))
              (when (< index limit)
                (write-sequence buffer output index limit)))
         until (let ((position (read-sequence buffer stream)))
                 (declare (fixnum position))
                 (setf (buffered-stream-index input) 0
                       (buffered-stream-limit input) position)
                 (= position 0))))
  #+sbcl
  (:method ((input stream) (output stream))
   (let ((buffered-input (make-buffered-stream input)))
     (declare (dynamic-extent buffered-input))
     (copy-stream buffered-input output))))
