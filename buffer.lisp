(in-package :elprep)

(defconstant +buffer-chunk-size+ 1024)

(declaim (inline make-buffer buffer-p))

(defstruct (buffer (:constructor make-buffer ()))
  (pos 0 :type fixnum)
  (str #() :type simple-vector)
  (hash-value -1 :type fixnum))

(declaim (inline reinitialize-buffer))

(defun reinitialize-buffer (buf)
  (declare (buffer buf) #.*fixnum-optimization*)
  (setf (buffer-pos buf) 0)
  (setf (buffer-hash-value buf) -1)
  buf)

(declaim (inline buffer-emptyp))

(defun buffer-emptyp (buf)
  (declare (buffer buf) #.*fixnum-optimization*)
  (= (buffer-pos buf) 0))

(defun ensure-str (buf old-str n)
  "internal"
  (declare (buffer buf) (simple-vector old-str) (fixnum n) #.*fixnum-optimization*)
  (let ((new-str (make-array n :single-thread t)))
    (declare (simple-vector new-str))
    (loop for i of-type fixnum below (length old-str)
          do (setf (svref new-str i) (svref old-str i)))
    (loop for i of-type fixnum from (length old-str) below (length new-str)
          do (setf (svref new-str i)
                   (make-array +buffer-chunk-size+ :element-type 'base-char :single-thread t)))
    (setf (buffer-str buf) new-str)))

(declaim (inline ensure-chunk))

(defun ensure-chunk (buf hi)
  "internal"
  (declare (buffer buf) (fixnum hi) #.*fixnum-optimization*)
  (let ((str (buffer-str buf)))
    (declare (simple-vector str))
    (svref (if (< hi (length str)) str
             (ensure-str buf str (1+ hi))) hi)))

(defmethod print-object ((buf buffer) stream)
  "internal"
  (print-unreadable-object (buf stream :type t :identity t)
    (format stream ":POS ~S :STR ~S"
            (buffer-pos buf)
            (unless (zerop (length (buffer-str buf))) "..."))))

(declaim (inline buffer-push))

(defun buffer-push (buf char)
  "add a character to a buffer"
  (declare (buffer buf) (base-char char) #.*fixnum-optimization*)
  (let ((pos (buffer-pos buf)))
    (declare (fixnum width pos))
    (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
      (declare (fixnum hi lo))
      (let ((chunk (ensure-chunk buf hi)))
        (declare (simple-base-string chunk))
        (setf (schar chunk lo) char)
        (setf (buffer-pos buf) (1+ pos))))))

(declaim (notinline buffer-push))

(declaim (notinline slow-buffer-extend))

(defun slow-buffer-extend (buf pos hi lo chunk string start end length)
  "internal"
  (declare (buffer buf) (fixnum pos hi lo) (simple-base-string chunk)
           (simple-base-string string) (fixnum start end length)
           #.*fixnum-optimization*)
  (loop with source of-type fixnum = start do
        (loop for target of-type fixnum from lo below +buffer-chunk-size+ do
              (setf (schar chunk target)
                    (schar string source))
              (when (= (incf source) end)
                (setf (buffer-pos buf) (+ pos length))
                (return-from slow-buffer-extend)))
        (incf hi) (setq lo 0)
        (setq chunk (ensure-chunk buf hi))))

(declaim (inline buffer-extend))

(defun buffer-extend (buf string &optional (start 0) end)
  "add a string to a buffer"
  (declare (buffer buf) (simple-base-string string) (fixnum start)
           #.*fixnum-optimization*)
  (let* ((end (or end (length string)))
         (length (- end start))
         (pos (buffer-pos buf)))
    (declare (fixnum end length pos))
    (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
      (declare (fixnum hi lo))
      (let ((chunk (ensure-chunk buf hi)))
        (declare (simple-base-string chunk))
        (if (<= (+ lo length) +buffer-chunk-size+)
          (loop for i of-type fixnum from start below end
                for j of-type fixnum from lo do
                (setf (schar chunk j)
                      (schar string i))
                finally
                (setf (buffer-pos buf) (+ pos length))
                (return (values)))
          (slow-buffer-extend buf pos hi lo chunk string start end length))))))

(declaim (notinline buffer-extend))

(defun write-buffer (buf stream)
  "write a buffer to a stream"
  (declare (buffer buf) (stream:buffered-stream stream)
           #.*fixnum-optimization*)
  (let ((pos (buffer-pos buf))
        (str (buffer-str buf)))
    (declare (fixnum pos) (simple-vector str))
    (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
      (declare (fixnum hi lo))
      (loop for i of-type fixnum below hi
            do (stream:stream-write-string stream (svref str i) 0 +buffer-chunk-size+)) ; buffer
      (when (> lo 0) (stream:stream-write-string stream (svref str hi) 0 lo)))) ; buffer
  (values))

(defun read-line-into-buffer (stream buf)
  "read a line into a buffer"
  (declare (stream:buffered-stream stream) (buffer buf)
           #.*fixnum-optimization*)
  (loop (stream:with-stream-input-buffer (buffer index limit) stream
          (declare (simple-base-string buffer) (fixnum index limit))
          (loop for i of-type fixnum from index below limit do
                (when (char= (schar buffer i) #\Newline)
                  (buffer-extend buf buffer index (1+ i)) ;make sure to include the newline
                  (setq index (1+ i))
                  (return-from read-line-into-buffer buf))
                finally
                (buffer-extend buf buffer index limit)
                (setq index limit)))
        (unless (stream:stream-fill-buffer stream)
          (return-from read-line-into-buffer buf))))

(defun buffer-partition (buf separator &rest targets)
  "get substrings from a buffer and feed them to target buffers;
   separator is a character, like #\Tab;
   targets is a property list with numbers as keys and buffers as values;
   the targets need to be sorted by key;
   for example (buffer-partition buf #\Tab 3 buf1 6 buf2)"
  (declare (buffer buf) (base-char separator) (dynamic-extent targets) ;
           #.*fixnum-optimization*)
  (let ((current-target 0))
    (declare (fixnum current-target))
    (flet ((get-target-buf ()
             (if targets
               (when (= current-target (the fixnum (car targets)))
                 (pop targets)
                 (pop targets))
               (return-from buffer-partition))))
      (declare (inline get-target-buf))
      (let ((target-buf (get-target-buf)))
        (declare ((or buffer null) target-buf))
        (flet ((next-target ()
                 (incf current-target)
                 (setq target-buf (get-target-buf))))
          (declare (inline next-target))
          (let ((pos (buffer-pos buf))
                (str (buffer-str buf)))
            (declare (fixnum pos) (simple-vector str))
            (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
              (declare (fixnum hi lo))
              (loop for i of-type fixnum below hi
                    for chunk of-type simple-base-string = (svref str i)
                    for start of-type fixnum = 0 do
                    (loop for end of-type fixnum below +buffer-chunk-size+ do
                          (when (char= (schar chunk end) separator)
                            (when target-buf
                              (buffer-extend target-buf chunk start end))
                            (next-target)
                            (setq start (+ end 1)))
                          finally
                          (when target-buf
                            (buffer-extend target-buf chunk start +buffer-chunk-size+))))
              (when (> lo 0)
                (loop with chunk of-type simple-base-string = (svref str hi)
                      with start of-type fixnum = 0
                      for end of-type fixnum below +buffer-chunk-size+ do
                      (when (char= (schar chunk end) separator)
                        (when target-buf
                          (buffer-extend target-buf chunk start end))
                        (next-target)
                        (setq start (+ end 1)))
                      finally
                      (when target-buf
                        (buffer-extend target-buf chunk start lo))))))))))
  (values))

(defun buffer-string (buf)
  "turn a buffer into a string representation; only for debugging"
  (declare (buffer buf) #.*fixnum-optimization*)
  (let ((pos (buffer-pos buf))
        (str (buffer-str buf)))
    (declare (fixnum pos) (simple-vector str))
    (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
      (declare (fixnum hi lo))
      (let ((result (make-array pos :element-type 'base-char :single-thread t))
            (target -1))
        (declare (simple-base-string result) (fixnum target))
        (loop for i of-type fixnum below hi
              for chunk of-type simple-base-string = (svref str i) do
              (loop for j of-type fixnum below +buffer-chunk-size+ do
                    (setf (schar result (incf target))
                          (schar chunk j))))
        (when (> lo 0)
          (loop with chunk of-type simple-base-string = (svref str hi)
                for j of-type fixnum below lo do
                (setf (schar result (incf target))
                      (schar chunk j))))
        result))))

(defun buffer= (buf1 buf2)
  "compare the contents of two buffers"
  (declare (buffer buf1 buf2) #.*fixnum-optimization*)
  (or (eq buf1 buf2)
      (let ((pos1 (buffer-pos buf1))
            (str1 (buffer-str buf1))
            (pos2 (buffer-pos buf2))
            (str2 (buffer-str buf2)))
        (declare (fixnum pos1 pos2) (simple-vector str1 str2))
        (when (= pos1 pos2)
          (multiple-value-bind (hi lo) (floor pos1 +buffer-chunk-size+)
            (declare (fixnum hi lo))
            (loop for i of-type fixnum below hi
                  for chunk1 of-type simple-base-string = (svref str1 i)
                  for chunk2 of-type simple-base-string = (svref str2 i) do
                  (loop for j of-type fixnum below +buffer-chunk-size+ do
                        (when (char/= (schar chunk1 j)
                                      (schar chunk2 j))
                          (return-from buffer= nil))))
            (when (> lo 0)
              (loop with chunk1 of-type simple-base-string = (svref str1 hi)
                    with chunk2 of-type simple-base-string = (svref str2 hi)
                    for j of-type fixnum below lo do
                    (when (char/= (schar chunk1 j)
                                  (schar chunk2 j))
                      (return-from buffer= nil)))))
          t))))

(defun buffer-parse-integer (buf)
  "convert a buffer to an integer"
  (declare (buffer buf) #.*optimization*)
  (let ((pos (buffer-pos buf))
        (str (buffer-str buf))
        (sign +1)
        (result 0))
    (declare (fixnum pos) (simple-vector str) (fixnum sign) (integer result))
    (flet ((update-result (char)
             (declare (base-char char))
             (assert (and (char<= #\0) (char<= #\9)))
             (let ((digit (- (char-code char) #.(char-code #\0))))
               (declare (fixnum digit))
               (if (and (typep result 'fixnum) (< (the fixnum result) #.(floor most-positive-fixnum 10)))
                 (setq result (the fixnum (+ (the fixnum (* (the fixnum result) 10)) digit)))
                 (setq result (+ (* result 10) digit))))))
      (declare (inline update-result))
      (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
        (declare (fixnum hi lo))
        (cond ((= hi 0)
               (assert (> lo 0))
               (let* ((chunk (svref str 0)) (char (schar chunk 0)) (start 0))
                 (declare (simple-base-string chunk) (base-char char) (fixnum start))
                 (cond ((char= char #\+)                (setq start 1) (assert (> lo 1)))
                       ((char= char #\-) (setq sign -1) (setq start 1) (assert (> lo 1))))
                 (loop for j of-type fixnum from start below lo
                       do (update-result (schar chunk j)))))
              (t (let* ((chunk (svref str 0)) (char (schar chunk 0)) (start 0))
                   (declare (simple-base-string chunk) (base-char char) (fixnum start))
                   (cond ((char= char #\+)                (setq start 1))
                         ((char= char #\-) (setq sign -1) (setq start 1)))
                   (loop for j of-type fixnum from start below +buffer-chunk-size+
                         do (update-result (schar chunk j))))
                 (loop for i of-type fixnum from 1 below hi
                       for chunk of-type simple-base-string = (svref str i) do
                       (loop for j of-type fixnum below +buffer-chunk-size+
                             do (update-result (schar chunk j))))
                 (when (> lo 0)
                   (loop with chunk of-type simple-base-string = (svref str hi)
                         for j of-type fixnum below lo
                         do (update-result (schar chunk j))))))))
    (* sign result)))

(declaim (inline rotate-1))

(defun rotate-1 (n)
  "internal"
  (declare (fixnum n) #.*fixnum-optimization*)
  (logior (ash n -1) (ash (logand n 1) #.(1- (integer-length most-positive-fixnum)))))

(defun buffer-hash (buf)
  "get the hash code for a buffer; once a hash code is computed, the buffer shouldn't change anymore!
   can be used for hash tables, like in (make-hash-table :test #'buffer= :hash-function #'buffer-hash)"
  (declare (buffer buf) #.*fixnum-optimization*)
  (let ((pos (buffer-pos buf))
        (str (buffer-str buf))
        (hash (buffer-hash-value buf)))
    (declare (fixnum pos) (simple-vector str) (fixnum hash))
    (if (> hash -1)
      (return-from buffer-hash hash)
      (setq hash 0))
    (multiple-value-bind (hi lo) (floor pos +buffer-chunk-size+)
      (declare (fixnum hi lo))
      (loop for i of-type fixnum below hi
            for chunk of-type simple-base-string = (svref str i) do
            (loop for j of-type fixnum below +buffer-chunk-size+ do
                  (setq hash (logxor (rotate-1 hash) (char-code (schar chunk j))))))
      (when (> lo 0)
        (loop with chunk of-type simple-base-string = (svref str hi)
              for j of-type fixnum below lo do
              (setq hash (logxor (rotate-1 hash) (char-code (schar chunk j)))))))
    (setf (buffer-hash-value buf) hash)))
