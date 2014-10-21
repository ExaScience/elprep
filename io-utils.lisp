(in-package :elprep)

(defmacro with-string-dispatch (string &body body)
  `(etypecase ,string
     (simple-base-string
      (let ((,string ,string))
        (declare (simple-base-string ,string))
        (macrolet ((dchar (string index) `(the base-char (schar ,string ,index))))
          ,@body)))
     #+lispworks
     (lw:simple-text-string
      (let ((,string ,string))
        (declare (lw:simple-text-string ,string))
        (macrolet ((dchar (string index) `(the base-char (schar ,string ,index))))
          ,@body)))
     ((simple-array character (*))
      (let ((,string ,string))
        (declare ((simple-array character (*)) ,string))
        (macrolet ((dchar (string index) `(the base-char (schar ,string ,index))))
          ,@body)))
     (string
      (let ((,string ,string))
        (declare (string ,string))
        (macrolet ((dchar (string index) `(the base-char (char ,string ,index))))
          ,@body)))))

#+sbcl
(progn
  (defconstant +ascii-stream-buffer-size+ (expt 2 13))
  (deftype ascii-stream-octet-buffer () `(simple-array octet (,+ascii-stream-buffer-size+)))
  (deftype ascii-stream-character-buffer () `(simple-array character (,+ascii-stream-buffer-size+)))
  (deftype ascii-stream-buffer () '(or ascii-stream-octet-buffer ascii-stream-character-buffer))

  (declaim (inline make-ascii-stream-buffer))
  
  (defun make-ascii-stream-buffer (element-type)
    (ecase element-type
      (octet (make-array +ascii-stream-buffer-size+ :element-type 'octet))
      (character (make-array +ascii-stream-buffer-size+ :element-type 'character))))
                   
  (defmacro with-buffer-dispatch (buffer &body body)
    `(etypecase ,buffer
       (ascii-stream-octet-buffer
        (flet ((bchar (buffer index)
                 (declare (ascii-stream-octet-buffer buffer) (fixnum index))
                 (the base-char (code-char (the octet (aref buffer index)))))
               ((setf bchar) (char buffer index)
                 (declare (base-char char) (ascii-stream-octet-buffer buffer) (fixnum index))
                 (setf (aref buffer index) (the octet (char-code char)))))
          (declare (inline bchar (setf bchar)) (ignorable #'bchar #'(setf bchar)))
          ,@body))
       (ascii-stream-character-buffer
        (flet ((bchar (buffer index)
                 (declare (ascii-stream-character-buffer buffer) (fixnum index))
                 (the base-char (schar buffer index)))
               ((setf bchar) (char buffer index)
                 (declare (base-char char) (ascii-stream-character-buffer buffer) (fixnum index))
                 (setf (schar buffer index) char)))
          (declare (inline bchar (setf bchar)) (ignorable #'bchar #'(setf bchar)))
          ,@body))))

  ;;; input streams

  (declaim (inline %make-ascii-stream))
  
  (defstruct (ascii-stream (:constructor %make-ascii-stream
                            (base element-type &aux (buffer (make-ascii-stream-buffer element-type))))
                           (:copier nil))
    (base nil :type stream)
    (index 0 :type fixnum)
    (limit +ascii-stream-buffer-size+ :type fixnum)
    (buffer #() :type ascii-stream-buffer :read-only t))
  
  (declaim (inline ascii-stream-fill-buffer))
  
  (defun ascii-stream-fill-buffer (stream)
    (declare (ascii-stream stream) #.*optimization*)
    (let ((position (read-sequence (ascii-stream-buffer stream)
                                   (ascii-stream-base stream))))
      (declare (ascii-stream stream) (fixnum position))
      (setf (ascii-stream-index stream) 0
            (ascii-stream-limit stream) position)
      (> position 0)))

  (declaim (inline make-ascii-stream))

  (defun make-ascii-stream (base &optional (element-type 'octet))
    (assert (input-stream-p base))
    (assert (subtypep (stream-element-type base) element-type))
    (let ((stream (%make-ascii-stream base element-type)))
      (ascii-stream-fill-buffer stream)
      stream))

  (declaim (inline close-ascii-stream))

  (defun close-ascii-stream (stream)
    (close (ascii-stream-base stream)))

  (defun ascii-stream-peek-char (stream eof-value)
    (declare (ascii-stream stream) #.*optimization*)
    (let ((buffer (ascii-stream-buffer stream)))
      (declare (ascii-stream-buffer buffer))
      (with-buffer-dispatch buffer
        (loop (let ((index (ascii-stream-index stream))
                    (limit (ascii-stream-limit stream)))
                (declare (fixnum index limit))
                (when (< index limit)
                  (return-from ascii-stream-peek-char (bchar buffer index))))
              (unless (ascii-stream-fill-buffer stream)
                (return-from ascii-stream-peek-char eof-value))))))

  (defun ascii-stream-read-line (stream)
    (declare (ascii-stream stream) #.*optimization*)
    (let ((buffer (ascii-stream-buffer stream))
          (strings '()))
      (declare (ascii-stream-buffer buffer) (list strings))
      (with-buffer-dispatch buffer
        (loop (let ((index (ascii-stream-index stream))
                    (limit (ascii-stream-limit stream)))
                (declare (fixnum index limit))
                (loop for end of-type fixnum from index below limit
                      for flag of-type boolean = (char= (bchar buffer end) #\Newline)
                      until flag finally
                      (let ((string (make-array (the fixnum (- end index)) :element-type 'base-char)))
                        (declare (simple-base-string string))
                        (loop for j of-type fixnum from index below end
                              for i of-type fixnum from 0
                              do (setf (schar string i) (bchar buffer j)))
                        (cond (flag (setf (ascii-stream-index stream) (the fixnum (1+ end)))
                                    (cond (strings (push string strings)
                                                   (return-from ascii-stream-read-line
                                                     (values (apply #'concatenate 'simple-base-string (nreverse strings)) nil)))
                                          (t       (return-from ascii-stream-read-line
                                                     (values string nil)))))
                              (t    (setf (ascii-stream-index stream) end)
                                    (push string strings))))))
              (unless (ascii-stream-fill-buffer stream)
                (if strings
                  (if (cdr strings)
                    (return-from ascii-stream-read-line (values (apply #'concatenate 'simple-base-string (nreverse strings)) t))
                    (return-from ascii-stream-read-line (values (car strings) t)))
                  (return-from ascii-stream-read-line (values empty-sbs t))))))))

  (defun ascii-stream-listen (stream)
    (declare (ascii-stream stream) #.*optimization*)
    (loop (let ((index (ascii-stream-index stream))
                (limit (ascii-stream-limit stream)))
            (declare (fixnum index limit))
            (when (< index limit)
              (return-from ascii-stream-listen t)))
          (unless (ascii-stream-fill-buffer stream)
            (return-from ascii-stream-listen nil))))

  ;;; output streams

  (declaim (inline ascii-stream-write-string))

  (defun ascii-stream-write-string (stream string)
    (write-string string stream))

  (declaim (inline ascii-stream-write-char))

  (defun ascii-stream-write-char (stream char)
    (write-char char stream)))

#+lispworks
(progn
  (declaim (inline ascii-stream-peek-char))
  
  (defun ascii-stream-peek-char (stream eof-value)
    (declare (buffered-stream stream) #.*optimization*)
    (loop (with-stream-input-buffer (buffer index limit) stream
            (declare (simple-base-string buffer) (fixnum index limit))
            (when (< index limit)
              (return-from ascii-stream-peek-char (lw:sbchar buffer index))))
          (unless (stream-fill-buffer stream)
            (return-from ascii-stream-peek-char eof-value))))

  (declaim (inline ascii-stream-read-line))

  (defun ascii-stream-read-line (stream)
    (stream-read-line stream))

  (declaim (inline ascii-stream-listen))
  
  (defun ascii-stream-listen (stream)
    (declare (buffered-stream stream) #.*optimization*)
    (loop (with-stream-input-buffer (buffer index limit) stream
            (declare (ignore buffer) (fixnum index limit))
            (when (< index limit)
              (return-from ascii-stream-listen t)))
          (unless (stream-fill-buffer stream)
            (return-from ascii-stream-listen nil))))

  (defun ascii-stream-write-string (stream string)
    (declare (buffered-stream stream) (string string) (fixnum start) #.*optimization*)
    (let ((start 0))
      (declare (fixnum start))
      (with-string-dispatch string
        (let ((end (length string)))
          (declare (fixnum end))
          (loop (with-stream-output-buffer (buffer index limit) stream
                  (declare (simple-base-string buffer) (fixnum index limit))
                  (loop for i of-type fixnum from index below limit do
                        (setf (lw:sbchar buffer i) (dchar string start))
                        (when (= (setq start (the fixnum (1+ start))) end)
                          (setq index (the fixnum (1+ i)))
                          (return-from ascii-stream-write-string string))
                        finally (setq index limit)))
                (stream-flush-buffer stream))))))

  (declaim (inline ascii-stream-write-char))

  (defun ascii-stream-write-char (stream char)
    (declare (buffered-stream stream) (base-char char) #.*optimization*)
    (loop (with-stream-output-buffer (buffer index limit) stream
            (declare (simple-base-string buffer) (fixnum index limit))
            (when (< index limit)
              (setf (lw:sbchar buffer index) char)
              (setq index (the fixnum (1+ index)))
              (return-from ascii-stream-write-char char)))
          (stream-flush-buffer stream))))

(declaim (inline open-program close-program program-stream))

#+lispworks
(progn
  (defun open-program (command &rest args &key (direction :input) (buffered t))
    "Similar to LispWorks's sys:open-pipe, to enable portability between LispWorks and SBCL."
    (declare (dynamic-extent args) (ignore direction buffered))
    (apply #'sys:open-pipe command args))

  (defun close-program (program)
    "Close a program opened with open-program."
    (close program))

  (defun program-stream (program)
    "Get the stream of a program."
    program))

#+sbcl
(progn
  (defstruct (program
              (:constructor open-program
               (command &key (direction :input) (external-format :default) (buffered t) &aux
                        (process 
                         (etypecase command
                           (string (sb-ext:run-program
                                    "/bin/sh" (list "-c" command)
                                    :wait nil (ecase direction (:input :output) (:output :input)) :stream
                                    :external-format external-format))
                           (list   (sb-ext:run-program
                                    (first command) (rest command)
                                    :wait nil (ecase direction (:input :output) (:output :input)) :stream
                                    :external-format external-format))))
                        (stream (ecase direction
                                  (:input (let ((stream (sb-ext:process-output process)))
                                            (if buffered (make-ascii-stream stream 'character) stream)))
                                  (:output (sb-ext:process-input process))))))
              (:copier nil))
    (process nil :read-only t)
    (stream nil :read-only t))

  (defun close-program (program)
    "Close a program opened with open-program."
    (let ((program-stream (program-stream program)))
      (when (and (streamp program-stream) (output-stream-p program-stream))
        (force-output program-stream)))
    (sb-ext:process-close (program-process program))))

(defmacro with-open-program ((program command &rest args) &body body)
  "Open a program for a block of code, and ensure that it is closed again in the end."
  `(let ((,program (open-program ,command ,@args)))
     (unwind-protect (progn ,@body)
       (close-program ,program))))

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

#+lispworks
(defun copy-stream (input output)
  "Efficient copying of the contents of an input stream to an output stream."
  (declare (buffered-stream input) (stream output) #.*optimization*)
  (loop do (with-stream-input-buffer (buffer index limit) input
             (declare (simple-base-string buffer) (fixnum index limit))
             (when (< index limit)
               (stream-write-sequence output buffer index limit)))
        while (stream-fill-buffer input))
  (values))

#+sbcl
(defun copy-stream (input output)
  "Efficient copying of the contents of an input stream to an output stream."
  (declare (ascii-stream input) (stream output) #.*optimization*)
  (let ((buffer (ascii-stream-buffer input)))
    (declare (ascii-stream-buffer buffer))
    (loop do (let ((index (ascii-stream-index input))
                   (limit (ascii-stream-limit input)))
               (declare (fixnum index limit))
               (when (< index limit)
                 (write-sequence buffer output :start index :end limit)))
          while (ascii-stream-fill-buffer input))
    (values)))
