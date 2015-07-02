(in-package :elprep)
(in-simple-base-string-syntax)

;;; input streams

#+lispworks
(defmacro with-buffer-dispatch (buffer &body body)
  "Specialize a body of code for different types of buffers in SBCL."
  (error "with-buffer-dispatch is not available in LispWorks: ~S."
         `(with-buffer-dispatch ,buffer ,@body)))

#+lispworks
(defmacro with-ascii-stream-input-buffer (buffer stream &body body)
  "Specialize a body of code for different types of a buffered-ascii-input-stream's input buffer in SBCL."
  (error "with-ascii-stream-input-buffer is not available in LispWorks: ~S."
         `(with-ascii-stream-input-buffer ,buffer ,stream ,@body)))

#+sbcl
(progn
  (defconstant +ascii-stream-buffer-size+ (expt 2 13)
    "Size of buffers in ascii-input-stream.")

  (declaim (inline make-buffered-ascii-input-stream))

  (defstruct (buffered-ascii-input-stream
              (:constructor make-buffered-ascii-input-stream
               (stream &optional (element-type 'base-char) &aux
                       (buffer (ecase element-type
                                 (base-char (make-array +ascii-stream-buffer-size+ :element-type 'octet))
                                 (character (make-array +ascii-stream-buffer-size+ :element-type 'character))))
                       (limit  (read-sequence buffer stream))))
              (:copier nil))
    "A wrapper for input streams that efficiently buffers ASCII input and allows for direct access of the buffer.
     This is somewhat similar to LispWorks's buffered-stream concept.
     The struct buffered-ascii-input-stream has a constructor make-buffered-ascii-input-stream that takes a stream and an optional element-type as parameters.
     Accessor buffered-ascii-input-stream-index refers to the current start position in the input buffer.
     Accessor buffered-ascii-input-stream-limit refers to the current end position in the input buffer.
     Accessor buffered-ascii-input-stream-buffer refers to the actual input buffer.
     Accessor buffered-ascii-input-stream-secondary-buffer refers to a secondary buffer that can be used to convert the primary buffer to a different type.
     Accessor buffered-ascii-input-stream-element-type refers to the element-type of the underlying stream.
     Accessor buffered-ascii-input-stream-stream refers to the underlying stream.
     The relevant Gray stream methods are specialized for this struct."
    (index 0 :type fixnum)
    (limit 0 :type fixnum)
    buffer
    (secondary-buffer nil)
    (element-type 'base-char :type symbol :read-only t)
    (stream nil :type (or null stream)))

  (setf (documentation 'make-buffered-ascii-input-stream 'function)
        "Constructor for struct buffered-ascii-input-stream that takes a stream and an optional element-type as parameters."
        (documentation 'buffered-ascii-input-stream-p 'function)
        "Default predicate for struct buffered-ascii-input-stream-p."
        (documentation 'buffered-ascii-input-stream-index 'function)
        "Current start position in the input buffer of a buffered-ascii-input-stream."
        (documentation 'buffered-ascii-input-stream-limit 'function)
        "Current end position in the input buffer of a buffered-ascii-input-stream."
        (documentation 'buffered-ascii-input-stream-buffer 'function)
        "The input buffer of a buffered-ascii-input-stream."
        (documentation 'buffered-ascii-input-stream-secondary-buffer 'function)
        "A secondary buffer for a buffered-ascii-input-stream used in copy-stream.
         May become obsolete when SBCL supports :element-type '(unsigned-byte 8) in sb-ext:run-program."
        (documentation 'buffered-ascii-input-stream-element-type 'function)
        "The element-type of stream underlying a buffered-ascii-input-stream.
         May become obsolete when SBCL supports :element-type '(unsigned-byte 8) in sb-ext:run-program."
        (documentation 'buffered-ascii-input-stream-stream 'function)
        "The stream underlying a buffered-ascii-input-stream.")

  (defmethod input-stream-p ((stream buffered-ascii-input-stream)) t)

  (declaim (inline stream-fill-buffer))

  (defun stream-fill-buffer (stream)
    "Fill the buffered-ascii-input-stream with data from the underlying stream.
     If data is available, update index and limit and return true."
    (declare (buffered-ascii-input-stream stream) #.*optimization*)
    (let ((position (read-sequence (buffered-ascii-input-stream-buffer stream)
                                   (buffered-ascii-input-stream-stream stream))))
      (declare (fixnum position))
      (when (> position 0)
        (setf (buffered-ascii-input-stream-index stream) 0
              (buffered-ascii-input-stream-limit stream) position)
        t)))

  (defmacro with-buffer-dispatch (buffer &body body)
    "Specialize a body of code for different types of buffers. In the code body, bchar can be used to access elements in the buffer."
    `(etypecase ,buffer
       ((simple-array octet (*))
        (let ((,buffer ,buffer))
          (declare ((simple-array octet (*)) ,buffer))
          (flet ((bchar (buffer index)
                   (declare ((simple-array octet (*)) buffer) (fixnum index))
                   (code-char (the octet (aref buffer index))))
                 ((setf bchar) (new-value buffer index)
                   (declare (base-char new-value) ((simple-array octet (*)) buffer) (fixnum index))
                   (setf (aref buffer index) (char-code new-value))
                   new-value))
            (declare (inline bchar (setf bchar)) (ignorable #'bchar #'(setf bchar)))
            ,@body)))
       ((simple-array character (*))
        (let ((,buffer ,buffer))
          (declare ((simple-array character (*)) ,buffer))
          (flet ((bchar (buffer index)
                   (declare ((simple-array character (*)) buffer) (fixnum index))
                   (aref buffer index))
                 ((setf bchar) (new-value buffer index)
                   (declare (character new-value) ((simple-array character (*)) buffer) (fixnum index))
                   (setf (aref buffer index) new-value)))
            (declare (inline bchar (setf bchar)) (ignorable #'bchar #'(setf bchar)))
            ,@body)))))

  (defmacro with-ascii-stream-input-buffer (buffer stream &body body)
    "Specialize a body of code for different types of buffers in a buffered-ascii-input-stream.
     In the code body, bchar can be used to access elements in the buffer."
    `(let ((,buffer (buffered-ascii-input-stream-buffer ,stream)))
       (with-buffer-dispatch ,buffer ,@body)))

  (defmethod stream-element-type ((stream buffered-ascii-input-stream))
    (buffered-ascii-input-stream-element-type stream))

  (defmethod close ((stream buffered-ascii-input-stream) &key abort)
    (declare (ignore abort))
    (setf (buffered-ascii-input-stream-stream stream) nil
          (buffered-ascii-input-stream-buffer stream) nil
          (buffered-ascii-input-stream-secondary-buffer stream) nil))

  (defmethod stream-file-position ((stream buffered-ascii-input-stream) &optional position)
    (declare #.*optimization*)
    (if position
      (when (file-position (buffered-ascii-input-stream-stream stream) position)
        (stream-fill-buffer stream) t)
      (+ (file-position (buffered-ascii-input-stream-stream stream))
         (the fixnum (- (the fixnum (buffered-ascii-input-stream-index stream))
                        (the fixnum (buffered-ascii-input-stream-limit stream)))))))

  (defmethod stream-clear-input ((stream buffered-ascii-input-stream)))

  (defmethod stream-read-sequence ((stream buffered-ascii-input-stream) sequence &optional (start 0) end)
    (declare #.*optimization*)
    (assert (typep start 'fixnum))
    (if end (assert (typep end 'fixnum))
      (setq end (length sequence)))
    (let ((start start) (end end))
      (declare (fixnum start end))
      (when (> end start)
        (with-ascii-stream-input-buffer buffer stream
          (loop (let* ((index (buffered-ascii-input-stream-index stream))
                       (limit (buffered-ascii-input-stream-limit stream))
                       (size  (min (the fixnum (- end start))
                                   (the fixnum (- limit index)))))
                  (declare (fixnum index limit size))
                  (setf (subseq sequence start end)
                        (subseq buffer index limit))
                  (setf (buffered-ascii-input-stream-index stream) (the fixnum (+ index size)))
                  (setq start (the fixnum (+ start size))))
                (unless (stream-fill-buffer stream)
                  (return-from stream-read-sequence start)))))))

  (defmethod stream-peek-char ((stream buffered-ascii-input-stream))
    (declare #.*optimization*)
    (with-ascii-stream-input-buffer buffer stream
      (loop (let ((index (buffered-ascii-input-stream-index stream))
                  (limit (buffered-ascii-input-stream-limit stream)))
              (declare (fixnum index limit))
              (when (< index limit)
                (return-from stream-peek-char (bchar buffer index))))
            (unless (stream-fill-buffer stream)
              (return-from stream-peek-char :eof)))))

  (defmethod stream-read-char-no-hang ((stream buffered-ascii-input-stream))
    (declare #.*optimization*)
    (with-ascii-stream-input-buffer buffer stream
      (loop (let ((index (buffered-ascii-input-stream-index stream))
                  (limit (buffered-ascii-input-stream-limit stream)))
              (declare (fixnum index limit))
              (when (< index limit)
                (let ((char (bchar buffer index)))
                  (setf (buffered-ascii-input-stream-index stream) (the fixnum (1+ index)))
                  (return-from stream-read-char-no-hang char))))
            (unless (stream-fill-buffer stream)
              (return-from stream-read-char-no-hang :eof)))))

  (defmethod stream-read-char ((stream buffered-ascii-input-stream))
    (declare #.*optimization*)
    (with-ascii-stream-input-buffer buffer stream
      (loop (let ((index (buffered-ascii-input-stream-index stream))
                  (limit (buffered-ascii-input-stream-limit stream)))
              (declare (fixnum index limit))
              (when (< index limit)
                (let ((char (bchar buffer index)))
                  (setf (buffered-ascii-input-stream-index stream) (the fixnum (1+ index)))
                  (return-from stream-read-char char))))
            (unless (stream-fill-buffer stream)
              (return-from stream-read-char :eof)))))

  (defmethod stream-read-line ((stream buffered-ascii-input-stream))
    (declare #.*optimization*)
    (let ((strings '()))
      (declare (list strings))
      (with-ascii-stream-input-buffer buffer stream
        (loop (let ((index (buffered-ascii-input-stream-index stream))
                    (limit (buffered-ascii-input-stream-limit stream)))
                (declare (fixnum index limit))
                (loop for end of-type fixnum from index below limit
                      for flag of-type boolean = (char= (bchar buffer end) #\Newline)
                      until flag finally
                      (let ((string (make-array (the fixnum (- end index)) :element-type 'base-char)))
                        (declare (simple-base-string string))
                        (loop for j of-type fixnum from index below end
                              for i of-type fixnum from 0
                              do (setf (schar string i) (bchar buffer j)))
                        (cond (flag (setf (buffered-ascii-input-stream-index stream) (the fixnum (1+ end)))
                                    (cond (strings (push string strings)
                                                   (return-from stream-read-line
                                                     (values (apply #'concatenate 'simple-base-string (nreverse strings)) nil)))
                                          (t       (return-from stream-read-line
                                                     (values string nil)))))
                              (t   (setf (buffered-ascii-input-stream-index stream) limit)
                                   (push string strings))))))
              (unless (stream-fill-buffer stream)
                (if strings
                  (if (cdr strings)
                    (return-from stream-read-line (values (apply #'concatenate 'simple-base-string (nreverse strings)) t))
                    (return-from stream-read-line (values (car strings) t)))
                  (return-from stream-read-line (values :eof t))))))))

  (defun skip-line (stream)
    "Skip characters from stream until a newline is reached.
     Can be used in place of read-line when its return value is discarded."
    (declare (buffered-ascii-input-stream stream) #.*optimization*)
    (with-ascii-stream-input-buffer buffer stream
      (loop (let ((index (buffered-ascii-input-stream-index stream))
                  (limit (buffered-ascii-input-stream-limit stream)))
              (declare (fixnum index limit))
              (loop for pos of-type fixnum from index below limit do
                    (when (char= (bchar buffer pos) #\Newline)
                      (setf (buffered-ascii-input-stream-index stream) (the fixnum (1+ pos)))
                      (return-from skip-line))
                    finally (setf (buffered-ascii-input-stream-index stream) pos)))
            (unless (stream-fill-buffer stream)
              (return-from skip-line)))))

  (defmethod stream-listen ((stream buffered-ascii-input-stream))
    (declare #.*optimization*)
    (or (let ((index (buffered-ascii-input-stream-index stream))
              (limit (buffered-ascii-input-stream-limit stream)))
          (declare (fixnum index limit))
          (< index limit))
        (listen (buffered-ascii-input-stream-stream stream))))

  (defmethod stream-unread-char ((stream buffered-ascii-input-stream) character)
    (declare (ignore character) #.*optimization*)
    (setf (buffered-ascii-input-stream-index stream)
          (the fixnum (1- (the fixnum (buffered-ascii-input-stream-index stream)))))
    nil)

  (defun copy-stream (input output)
    "Efficient copying of the contents of a buffered-ascii-input-stream to an output stream."
    (declare (buffered-ascii-input-stream input) (stream output) #.*optimization*)
    (ecase (buffered-ascii-input-stream-element-type input)
      (base-char (let ((buffer (buffered-ascii-input-stream-buffer input))
                       (second (or (buffered-ascii-input-stream-secondary-buffer input)
                                   (setf (buffered-ascii-input-stream-secondary-buffer input)
                                         (make-array +ascii-stream-buffer-size+ :element-type 'base-char)))))
                   (declare ((simple-array octet (*)) buffer) (simple-base-string second))
                   (loop (let ((index (buffered-ascii-input-stream-index input))
                               (limit (buffered-ascii-input-stream-limit input)))
                           (declare (fixnum index limit))
                           (when (< index limit)
                             (loop for i of-type fixnum from index below limit do
                                   (setf (schar second i) (code-char (the octet (aref buffer i)))))
                             (write-string second output :start index :end limit)
                             (setf (buffered-ascii-input-stream-index input) limit)))
                         (unless (stream-fill-buffer input)
                           (return-from copy-stream (values))))))
      (character (let ((buffer (buffered-ascii-input-stream-buffer input)))
                   (declare ((simple-array character (*)) buffer))
                   (loop (let ((index (buffered-ascii-input-stream-index input))
                               (limit (buffered-ascii-input-stream-limit input)))
                           (declare (fixnum index limit))
                           (when (< index limit)
                             (write-string buffer output :start index :end limit)
                             (setf (buffered-ascii-input-stream-index input) limit)))
                         (unless (stream-fill-buffer input)
                           (return-from copy-stream (values)))))))))

#+lispworks
(progn
  (defun skip-line (stream)
    "Skip characters from stream until a newline is reached.
     Can be used in place of read-line when its return value is discarded."
    (declare (buffered-stream stream) #.*fixnum-optimization*)
    (loop (with-stream-input-buffer (source index limit) stream
            (declare (simple-base-string source) (fixnum index limit))
            (loop for pos of-type fixnum from index below limit do
                  (when (char= (lw:sbchar source pos) #\Newline)
                    (setq index (1+ pos))
                    (return-from skip-line))
                  finally (setq index pos)))
          (unless (stream-fill-buffer stream)
            (return-from skip-line))))

  (defun copy-stream (input output)
    "Efficient copying of the contents of a buffered-stream to an output stream."
    (declare (buffered-stream input) (stream output) #.*optimization*)
    (loop (with-stream-input-buffer (buffer index limit) input
            (declare (simple-base-string buffer) (fixnum index limit))
            (when (< index limit)
              (stream-write-string output buffer index limit)
              (setq index limit)))
          (unless (stream-fill-buffer input)
            (return-from copy-stream (values))))))


;;; output streams

(declaim (inline writec write-newline write-tab))

#+sbcl
(defun writec (out c)
  "Write a character to an output stream."
  (write-char c out))

#+lispworks
(defun writec (out c)
  "Write a character to an output stream."
  (declare (buffered-stream out) (base-char c) #.*optimization*)
  (loop (with-stream-output-buffer (buffer index limit) out
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (setf (lw:sbchar buffer index) c)
            (setq index (the fixnum (1+ index)))
            (return-from writec c)))
        (stream-flush-buffer out)))

(defun write-newline (out)
  "Write a newline to an output stream."
  (writec out #\Newline))

(defun write-tab (out)
  "Write a tabulator to an output stream."
  (writec out #\Tab))

#+sbcl
(progn
  (declaim (inline writestr writeln))
  
  (defun writestr (out string)
    "Write a string to an output stream."
    (write-string string out))

  (defun writeln (out string &aux (length (length string)))
    "Write a simple-base-string to an output stream, but only up to a #\Newline."
    (declare (simple-base-string string) (fixnum length) #.*optimization*)
    (write-string string out :end (the fixnum (1+ (the fixnum (loop for i of-type fixnum below length
                                                                    when (char= (schar string i) #\Newline) return i
                                                                    finally (return length))))))))

#+lispworks
(progn 
  (defun writestr (out string &aux (end (length string)))
    "Write a base-string to an output stream."
    (declare (buffered-stream out) (base-string string) (fixnum end) #.*optimization*)
    (multiple-value-bind (string* start) (unwrap-displaced-array string)
      (declare (simple-base-string string*) (fixnum start))
      (setq end (the fixnum (+ start end)))
      (loop (with-stream-output-buffer (buffer index limit) out
              (declare (simple-base-string buffer) (fixnum index limit))
              (loop for i of-type fixnum from index below limit do
                    (setf (lw:sbchar buffer i) (lw:sbchar string* start))
                    (when (= (setq start (the fixnum (1+ start))) end)
                      (setq index (the fixnum (1+ i)))
                      (return-from writestr string))
                    finally (setq index limit)))
            (stream-flush-buffer out))))

  (defun writeln (out string)
    "Write a simple-base-string to an output stream, but only up to a #\Newline."
    (declare (buffered-stream out) (simple-base-string string) #.*optimization*)
    (let ((start 0))
      (declare (fixnum start))
      (loop (with-stream-output-buffer (buffer index limit) out
              (declare (simple-base-string buffer) (fixnum index limit))
              (loop for i of-type fixnum from index below limit do
                    (cond ((char= (setf (lw:sbchar buffer i) (lw:sbchar string start)) #\Newline)
                           (setq index (the fixnum (1+ i)))
                           (return-from writeln string))
                          (t (setq start (the fixnum (1+ start)))))
                    finally (setq index limit)))
            (stream-flush-buffer out)))))

(declaim (inline make-sim-stream))

(defstruct (sim-stream (:constructor make-sim-stream
                        (&optional (initial-size 128) &aux
                                   (string (make-array initial-size :element-type 'base-char
                                                       #+lispworks :single-thread #+lispworks t)))))
  "A pseudo string stream to allow for more low-level efficient implementations of formatting functions speficic to elPrep.
   The struct sim-stream has a constructor make-sim-stream that takes an optional initial size for the target string as a parameter.
   Accessor sim-stream-string refers to the target string for output.
   Accessor sim-stream-index refers to the current index into sim-stream-string.
   Accessor sim-stream-%floats refers to an optional actual string-stream for floating point output."
  (string "" :type simple-base-string)
  (index   0 :type fixnum)
  (%floats nil))

(setf (documentation 'make-sim-stream 'function)
      "Constructor for struct sim-stream that takes an optional initial size for the target string as a parameter."
      (documentation 'sim-stream-p 'function)
      "Default predicate for struct sim-stream."
      (documentation 'copy-sim-stream 'function)
      "Default copier function for struct sim-stream."
      (documentation 'sim-stream-string 'function)
      "Access the sim-stream target string."
      (documentation 'sim-stream-index 'function)
      "Access the sim-stream index into sim-stream-string."
      (documentation 'sim-stream-%floats 'function)
      "Access the optional sim-stream string-stream for formatting floating point numbers.")

(declaim (inline sim-stream-floats))

(defun sim-stream-floats (s)
  "Read the sim-stream string-stream for formatting floating point numbers, creating a new one if not already available."
  (declare (sim-stream s) #.*optimization*)
  (or (sim-stream-%floats s)
      (setf (sim-stream-%floats s)
            (make-string-output-stream :element-type 'base-char))))

(defun slow-ensure-sim-space (out string index required-length)
  "Ensure there is enough space in the target string of sim-stream out at index with a required-length of additional characters, slow path."
  (declare (sim-stream out) (simple-base-string string) (fixnum index required-length) #.*optimization*)
  (let ((new-string (make-array (the fixnum (+ (the fixnum (* 128 (the fixnum (floor required-length 128)))) 128))
                                :element-type 'base-char #+lispworks :single-thread #+lispworks t)))
    (declare (simple-base-string new-string))
    (loop for i of-type fixnum below index
          do (setf (schar new-string i) (schar string i)))
    (setf (sim-stream-string out) new-string)))

(declaim (inline ensure-sim-space))

(defun ensure-sim-space (out length)
  "Ensure there is enough space in the target string of sim-stream out, with length additional characters required, fast path."
  (declare (sim-stream out) (fixnum length) #.*optimization*)
  (let ((string (sim-stream-string out))
        (index  (sim-stream-index out)))
    (declare (simple-base-string string) (fixnum index))
    (let ((required-length (+ index length)))
      (declare (fixnum required-length))
      (if (< required-length (length string)) string
        (slow-ensure-sim-space out string index required-length)))))

(declaim (inline sim-writec sim-write-newline sim-write-tab))

(defun sim-writec (out c)
  "Write a character c to sim-stream out."
  (declare (sim-stream out) (base-char c) #.*optimization*)
  (let ((string (ensure-sim-space out 1))
        (index  (sim-stream-index out)))
    (declare (simple-base-string string) (fixnum index))
    (setf (schar string index) c)
    (setf (sim-stream-index out) (the fixnum (1+ index)))))

(defun sim-write-newline (out)
  "Write a newline to sim-stream out."
  (sim-writec out #\Newline))

(defun sim-write-tab (out)
  "Write a tab to sim-stream out."
  (sim-writec out #\Tab))

(defun sim-writestr (out string &aux (length (length string)))
  "Write a string to sim-stream out."
  (declare (sim-stream out) (string string) (fixnum length) #.*optimization*)
  (multiple-value-bind (source start) (unwrap-displaced-array string)
    (declare (simple-string source) (fixnum start))
    (let ((target (ensure-sim-space out length))
          (index  (sim-stream-index out))
          (end (the fixnum (+ start length))))
      (declare (simple-base-string target) (fixnum index end))
      (loop for i of-type fixnum from start below end
            for j of-type fixnum from index
            do (setf (schar target j) (schar source i)))
      (setf (sim-stream-index out) (the fixnum (+ index length))))))

(defun sim-write-integer (out integer)
  "Write a decimal integer to sim-stream out."
  (declare (sim-stream out) (integer integer) #.*optimization*)
  (when (< integer 0)
    (sim-writec out #\-)
    (setq integer (- integer)))
  (let ((start (sim-stream-index out)))
    (declare (fixnum start))
    (flet ((truncate-fixnum-loop (fixnum)
             (declare (fixnum fixnum))
             (multiple-value-bind (div rem) (truncate fixnum 10)
               (declare (fixnum div rem))
               (loop do (sim-writec out (code-char (the fixnum (+ #.(char-code #\0) rem))))
                     until (= div 0)
                     do (setf (values div rem) (truncate div 10))))))
      (declare (inline truncate-fixnum-loop))
      (if (typep integer 'fixnum)
        (truncate-fixnum-loop integer)
        (multiple-value-bind (div rem) (truncate integer 10)
          (declare (integer div) (fixnum rem))
          (loop do (sim-writec out (code-char (the fixnum (+ #.(char-code #\0) rem))))
                until (typep div 'fixnum)
                do (setf (values div rem) (truncate div 10))
                finally (truncate-fixnum-loop div)))))
    (let* ((string (sim-stream-string out))
           (end (sim-stream-index out))
           (half (floor (the fixnum (- end start)) 2)))
      (declare (fixnum end half))
      (setq end (the fixnum (1- end)))
      (loop for i of-type fixnum below half
            do (rotatef (schar string (the fixnum (+ start i)))
                        (schar string (the fixnum (- end i))))))))

(declaim (inline sim-write-byte))

(defun sim-write-byte (out byte)
  "Write a byte in hex notation to sim-stream out."
  (declare (sim-stream out) (fixnum byte) #.*optimization*)
  (flet ((sim-write-nibble (nibble)
           (declare (fixnum nibble))
           (sim-writec out (code-char (the fixnum (+ (if (< nibble 10)
                                                       #.(char-code #\0)
                                                       #.(- (char-code #\A) 10))
                                                     nibble))))))
    (declare (inline sim-write-nibble))
    (multiple-value-bind (hi lo) (floor byte 16)
      (declare (fixnum hi lo))
      (sim-write-nibble hi)
      (sim-write-nibble lo))))

(defun sim-write-fixed-size-fixnum (out fixnum size)
  "Write a positive fixnum of known maximum width to sim-stream out."
  (declare (sim-stream out) (fixnum fixnum size) #.*optimization*)
  (let* ((string (ensure-sim-space out size))
         (index  (sim-stream-index out))
         (new-index (+ index size)))
    (declare (simple-base-string string) (fixnum index new-index))
    (multiple-value-bind (div rem) (truncate fixnum 10)
      (declare (fixnum div rem))
      (loop for i of-type fixnum from (the fixnum (1- new-index)) downto index do
            (setf (schar string i) (code-char (the fixnum (+ #.(char-code #\0) rem))))
            (setf (values div rem) (truncate div 10))))
    (setf (sim-stream-index out) new-index)))

(defun sim-write-float (out float)
  "Write a floating point number to sim-stream out."
  (declare (sim-stream out) (single-float float) #.*optimization*)
  (let ((stream (sim-stream-floats out)))
    (format stream "~E" float)
    (sim-writestr out (get-output-stream-string stream))))

#+lispworks
(progn
  (defstruct process
    "This struct represents the return values of sys:run-shell-command in LispWorks.
     It has default constructor, predicate, and copier.
     Read-only accessor process-stream refers to the stream of a process.
     Read-only accessor process-error refers to the error-stream of a process.
     Read-only accessor process-pid refers to the PID of a process."
    (stream nil :read-only t)
    (error nil :read-only t)
    #+lispworks6 (pid nil :read-only t))

  (setf (documentation 'make-process 'function)
        "Default constructor for struct process."
        (documentation 'copy-process 'function)
        "Default copier for struct process."
        (documentation 'process-p 'function)
        "Default predicate for struct process."
        (documentation 'process-stream 'function)
        "The stream of an external process."
        (documentation 'process-error 'function)
        "The error stream of an external process."
        #+lispworks6 (documentation 'process-pid 'function)
        #+lispworks6 "The PID of an external process.")

  (declaim (inline process-input process-output))

  (defun process-input (process)
    "Returns the result of process-stream. Exists for compatibility with SBCL's sb-ext:process-input."
    (declare (process process) #.*optimization*)
    (process-stream process))

  (defun process-output (process)
    "Returns the result of process-stream. Exists for compatibility with SBCL's sb-ext:process-output."
    (declare (process process) #.*optimization*)
    (process-stream process))

  (defun run-program (command command-args &key
                              (input nil)
                              (output nil)
                              (error nil)
                              (wait t)
                              (if-input-does-not-exist nil)
                              (if-output-exists :error)
                              (if-error-exists :error)
                              (external-format :latin-1))
    "A wrapper for sys:run-shell-command in LispWorks to make it behave more like SBCL's sb-ext:run-program."
    (declare (ignore external-format))
    (multiple-value-bind
        (stream error-stream #+lispworks6 pid)
        (sys:run-shell-command (cons command command-args)
                               :input input :output output :error-output error :wait wait
                               :if-input-does-not-exist if-input-does-not-exist
                               :if-output-exists if-output-exists
                               :if-error-output-exists if-error-exists
                               :save-exit-status t)
      (if wait stream
        (make-process :stream stream
                      :error error-stream
                      #+lispworks6 :pid #+lispworks6 pid))))

  (defun process-close (program &key (wait t))
    "Close all streams of an external process, wait for the process to finish if requested, and return its exit code if available."
    (let (#-lispworks6 (exit-code (sys:pipe-exit-status (process-stream program) :wait wait)))
      (let ((stream (process-stream program)))
        (when stream (close stream)))
      (let ((error (process-error program)))
        (when error (close error)))
      #+lispworks6 (sys:pid-exit-status (process-pid program) :wait wait)
      #-lispworks6 exit-code)))

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

(declaim (inline command-line-arguments))

(defun command-line-arguments ()
  "Wrapper for sys:*line-arguments-list* in LispWorks, and sb-ext:*posix-argv* in SBCL."
  #+lispworks sys:*line-arguments-list*
  #+sbcl sb-ext:*posix-argv*)
