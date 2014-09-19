(in-package :elprep)

;;; input

(defconstant +eof+ (code-char 255)
  "Constant character value for indicating an end of file.")

(declaim (inline advance
                 peekc readc scanc
                 end-of-entry-p
                 end-of-line-p
                 end-of-entry-or-line-p
                 end-of-file-p))

(defun advance (stream)
  "Advance stream by one character."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop (with-stream-input-buffer (buffer index limit) stream
          (declare (ignore buffer) (fixnum index limit))
          (when (< index limit)
            (incf index)
            (return-from advance)))
        (unless (stream-fill-buffer stream)
          (return-from advance))))

(defun peekc (stream)
  "Peek current character in stream."
  (declare (buffered-stream stream) #.*optimization*)
  (loop (with-stream-input-buffer (buffer index limit) stream
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (return-from peekc (sbchar buffer index))))
        (unless (stream-fill-buffer stream)
          (return-from peekc +eof+))))

(defun readc (stream)
  "Read current character in stream."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop (with-stream-input-buffer (buffer index limit) stream
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (return-from readc (prog1 (sbchar buffer index)
                                 (incf index)))))
        (unless (stream-fill-buffer stream)
          (return-from readc +eof+))))

(defun scanc (stream char)
  "Read current character in stream, ensuring it's the given one."
  (declare (buffered-stream stream) (base-char char) #.*fixnum-optimization*)
  (assert (char= (readc stream) char))
  char)

(defun end-of-entry-p (char)
  "Is the character a tab, delimiting an entry in a tab-delimited text file?"
  (declare (base-char char) #.*fixnum-optimization*)
  (char= char #\Tab))

(defun end-of-line-p (char)
  "Is the character an end of line in a text file?"
  (declare (base-char char) #.*fixnum-optimization*)
  (or (char= char #\Newline)
      (char= char +eof+)))

(defun end-of-entry-or-line-p (char)
  "Is the character a tab or end of line, delimiting an entry in a tab-delimited text file?"
  (declare (base-char char) #.*fixnum-optimization*)
  (or (char= char #\Tab)
      (char= char #\Newline)
      (char= char +eof+)))

(defun end-of-file-p (char)
  "Is the character and end of file in a text file?"
  (declare (base-char char) #.*fixnum-optimization*)
  (char= char +eof+))

(declaim (inline %scan-string))

(defun %scan-string (stream predicate)
  "Read characters from stream into a string until predicate on a character returns true."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (let ((strings '()))
    (declare (list strings))
    (loop (with-stream-input-buffer (source index limit) stream
            (declare (simple-base-string source) (fixnum index limit))
            (loop for end of-type fixnum from index below limit
                  for flag of-type boolean = (funcall predicate (sbchar source end))
                  until flag finally
                  (let ((string (make-array (- end index) :element-type 'base-char)))
                    (declare (simple-base-string string))
                    (loop for j of-type fixnum from index below end
                          for i of-type fixnum from 0
                          do (setf (sbchar string i) (sbchar source j)))
                    (setq index end)
                    (if flag
                      (cond (strings (push string strings)
                                     (return-from %scan-string (apply 'string-append (nreverse strings))))
                            (t       (return-from %scan-string string)))
                      (push string strings)))))
          (unless (stream-fill-buffer stream)
            (if strings
              (if (cdr strings)
                (return-from %scan-string (apply 'string-append (nreverse strings)))
                (return-from %scan-string (car strings)))
              (return-from %scan-string ""))))))

(defun scan-string (stream)
  "Read a tab-delimited entry as string from stream."
  (declare (buffered-stream stream) #.*optimization*)
  (%scan-string stream 'end-of-entry-p))

(defun scan-stringn (stream)
  "Read a tab- or end-of-line-delimited entry as string from stream."
  (declare (buffered-stream stream) #.*optimization*)
  (%scan-string stream 'end-of-entry-or-line-p))

(declaim (inline not-digit-p))

(defun not-digit-p (char)
  "Is the character not a digit?"
  (declare (base-char char) #.*fixnum-optimization*)
  (not (and (char<= #\0 char) (char<= char #\9))))

(defun scan-digits (stream)
  "Read a string of digits from stream."
  (declare (buffered-stream stream) #.*optimization*)
  (%scan-string stream 'not-digit-p))

(defun slow-scan-integer (stream)
  "Scan an integer from stream, slow path."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (let* ((negative (let ((char (peekc stream)))
                     (declare (base-char char))
                     (cond ((char= char #\+) (advance stream) nil)
                           ((char= char #\-) (advance stream) t))))
         (string   (scan-digits stream))
         (value    (parse-integer string)))
    (if negative (- value) value)))

(defun scan-integer (stream)
  "Scan an integer from stream."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop (with-stream-input-buffer (buffer index limit) stream
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (let* ((pos index) (char (sbchar buffer pos)))
              (declare (fixnum pos) (base-char char))
              (flet ((nextc () (if (= (incf pos) limit)
                                 (return-from scan-integer (slow-scan-integer stream))
                                 (setq char (sbchar buffer pos)))))
                (declare (inline nextc))
                (when (or (char= char #\-)
                          (char= char #\+))
                  (nextc))
                (assert (and (char<= #\0 char) (char<= char #\9)))
                (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9)))
                (return-from scan-integer (prog1 (parse-integer buffer :start index :end pos)
                                            (setq index pos)))))))
        (unless (stream-fill-buffer stream)
          (error "Stream overflow in scan-integer."))))

(defun slow-scan-float (stream)
  "Scan a floating point number from stream, slow path. Returns either an integer or a single-float."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (let ((char (peekc stream)))
    (declare (base-char char))
    (flet ((peekc () (setq char (peekc stream)))
           (advance () (advance stream)))
      (declare (inline peekc advance))
      (let* ((sign        (cond ((char= char #\+) (advance) "+")
                                ((char= char #\-) (advance) "-")
                                (t "")))
             (digits      (prog1 (scan-digits stream) (peekc)))
             (dot         (cond ((char= char #\.) (advance) ".")
                                (t "")))
             (dot-digits  (if (= 0 (length dot)) ""
                            (prog1 (scan-digits stream) (peekc))))
             (exp         (cond ((char= char #\e) (advance) (peekc) "e")
                                ((char= char #\E) (advance) (peekc) "E")
                                (t "")))
             (exp-sign    (if (= 0 (length exp)) ""
                            (cond ((char= char #\+) (advance) "+")
                                  ((char= char #\-) (advance) "-")
                                  (t ""))))
             (exp-digits  (if (= 0 (length exp)) ""
                            (scan-digits stream)))
             (string      (string-append sign digits dot dot-digits exp exp-sign exp-digits))
             (result      (read-from-string string)))
        (declare (string sign digits dot dot-digits exp exp-sign exp-digits string))
        (etypecase result
          (integer result)
          (single-float result)
          (number (coerce result 'single-float)))))))

(defun scan-float (stream)
  "Scan a floating point number from stream. Returns either an integer or a single-float."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop (with-stream-input-buffer (buffer index limit) stream
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (let* ((pos index) (char (sbchar buffer pos)))
              (declare (fixnum pos) (base-char char))
              (flet ((nextc () (if (= (incf pos) limit)
                                 (return-from scan-float (slow-scan-float stream))
                                 (setq char (sbchar buffer pos)))))
                (declare (inline nextc))
                (when (or (char= char #\+)
                          (char= char #\-))
                  (nextc))
                (let ((digits (and (char<= #\0 char) (char<= char #\9))))
                  (declare (boolean digits))
                  (when digits
                    (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9))))
                  (cond ((char= char #\.)
                         (nextc)
                         (assert (and (char<= #\0 char) (char<= char #\9)))
                         (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9))))
                        (t (assert digits))))
                (when (or (char= char #\e)
                          (char= char #\E))
                  (nextc)
                  (when (or (char= char #\+)
                            (char= char #\-))
                    (nextc))
                  (assert (and (char<= #\0 char) (char<= char #\9)))
                  (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9))))
                (let ((result (read-from-string buffer t nil :start index :end pos)))
                  (setq index pos)
                  (return-from scan-float
                    (etypecase result
                      (integer result)
                      (single-float result)
                      (number (coerce result 'single-float)))))))))
        (unless (stream-fill-buffer stream)
          (error "Stream overflow in scan-float."))))

(defun parse-sam-tag (stream)
  "Parse the TAG: portion of a SAM file tag."
  (declare (buffered-stream stream) #.*optimization*)
  (let ((tag-string (make-array 2 :element-type 'base-char)))
    (declare (simple-base-string tag-string))
    (setf (sbchar tag-string 0) (readc stream))
    (setf (sbchar tag-string 1) (readc stream))
    (scanc stream #\:)
    tag-string))

(defun parse-sam-byte-array (stream)
  "Parse a byte array in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type H."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (flet ((hex-value (char)
           (declare (base-char char))
           (cond ((and (char<= #\0 char) (char<= char #\9))
                  (- (char-int char) #.(char-int #\0)))
                 ((and (char<= #\a char) (char<= char #\f))
                  (- (char-int char) #.(- (char-int #\a) 10)))
                 ((and (char<= #\A char) (char<= char #\F))
                  (- (char-int char) #.(- (char-int #\A) 10)))
                 (t (error "Not a hex digit ~S." char)))))
    (declare (inline hex-value))
    (loop for count of-type fixnum from 1
          collect (+ (ash (hex-value (readc stream)) 4)
                     (hex-value (readc stream)))
          into list until (end-of-entry-or-line-p (peekc stream))
          finally (return (make-array count :element-type '(unsigned-byte 8) :initial-contents list)))))

(defun parse-sam-numeric-array (stream)
  "Parse a numeric array in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type B."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (advance stream) ;; ignore type tag
  (loop do (scanc stream #\,)
        collect (scan-float stream)
        until (end-of-entry-or-line-p (peekc stream))))

(defun parse-sam-header-line (stream)
  "Parse an @HD line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop until (end-of-line-p (peekc stream))
        nconc (let ((tag (progn
                             (scanc stream #\Tab)
                             (parse-sam-tag stream))))
                  (declare (simple-base-string tag))
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key tag) record) (scan-stringn stream))
                                               (error "Unknown tag ~A in SAM header line." tag)))
                    ("VN" (list (unique :VN record) (scan-stringn stream)))
                    ("SO" (list (unique :SO record) (scan-stringn stream)))))
        into record finally
        (advance stream)
        (unless (presentp :VN record)
          (cerror "Ignore absence of VN tag and continue."
                  "VN tag missing in @HD line when reading ~A." stream))
        (return record)))

(defun parse-sam-reference-sequence-dictionary-entry (stream)
  "Parse an @SQ line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop until (end-of-line-p (peekc stream))
        nconc (let ((tag (progn
                             (scanc stream #\Tab)
                             (parse-sam-tag stream))))
                  (declare (simple-base-string tag))
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key tag) record) (scan-stringn stream))
                                               (error "Unknown tag ~A in SAM reference sequence dictionary line." tag)))
                    ("SN" (list (unique :SN record) (scan-stringn stream)))
                    ("LN" (list (unique :LN record) (scan-integer stream)))
                    ("AS" (list (unique :AS record) (scan-stringn stream)))
                    ("M5" (list (unique :M5 record) (parse-sam-byte-array stream)))
                    ("SP" (list (unique :SP record) (scan-stringn stream)))
                    ("UR" (list (unique :UR record) (scan-stringn stream)))))
        into record finally
        (advance stream)
        (unless (presentp :SN record)
          (cerror "Ignore absence of SN tag and continue."
                  "SN tag missing in @SQ line when reading ~A." stream))
        (unless (presentp :LN record)
          (cerror "Ignore absence of LN tag and continue."
                  "LN tag missing in @SQ line when reading ~A." stream))
        (return record)))

(defun parse-sam-read-group (stream)
  "Parse an @RG line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop until (end-of-line-p (peekc stream))
        nconc (let ((tag (progn
                             (scanc stream #\Tab)
                             (parse-sam-tag stream))))
                  (declare (simple-base-string tag))
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key tag) record) (scan-stringn stream))
                                               (error "Unknown tag ~A in SAM read group line." tag)))
                    ("ID" (list (unique :ID record) (scan-stringn stream)))
                    ("CN" (list (unique :CN record) (scan-stringn stream)))
                    ("DS" (list (unique :DS record) (scan-stringn stream)))
                    ("DT" (list (unique :DT record) (parse-date-time (scan-stringn stream))))
                    ("FO" (list (unique :FO record) (scan-stringn stream)))
                    ("KS" (list (unique :KS record) (scan-stringn stream)))
                    ("LB" (list (unique :LB record) (scan-stringn stream)))
                    ("PG" (list (unique :PG record) (scan-stringn stream)))
                    ("PI" (list (unique :PI record) (scan-integer stream)))
                    ("PL" (list (unique :PL record) (scan-stringn stream)))
                    ("PU" (list (unique :PU record) (scan-stringn stream)))
                    ("SM" (list (unique :SM record) (scan-stringn stream)))))
        into record finally
        (advance stream)
        (unless (presentp :ID record)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @RG line when reading ~A." stream))
        (return record)))

(defun parse-sam-program (stream)
  "Parse an @PG line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop until (end-of-line-p (peekc stream))
        nconc (let ((tag (progn
                             (scanc stream #\Tab)
                             (parse-sam-tag stream))))
                  (declare (simple-base-string tag))
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key tag) record) (scan-stringn stream))
                                               (error "Unknown tag ~A in SAM program line." tag)))
                    ("ID" (list (unique :ID record) (scan-stringn stream)))
                    ("PN" (list (unique :PN record) (scan-stringn stream)))
                    ("CL" (list (unique :CL record) (scan-stringn stream)))
                    ("PP" (list (unique :PP record) (scan-stringn stream)))
                    ("DS" (list (unique :DS record) (scan-stringn stream)))
                    ("VN" (list (unique :VN record) (scan-stringn stream)))))
        into record finally
        (advance stream)
        (unless (presentp :ID record)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @PG line when reading ~A." stream))
        (return record)))

(defun parse-sam-comment (stream)
  "Parse an @CO line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (when (char= (peekc stream) #\Tab)
    (advance stream))
  (prog1 (%scan-string stream 'end-of-line-p)
    (advance stream)))

(defun parse-sam-header (stream)
  "Parse a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (let (hd sq rg pg co user-tags)
    (loop until (let ((char (peekc stream)))
                  (declare (base-char char))
                  (or (end-of-file-p char) (char/= char #\@)))
          for code of-type simple-base-string = (scan-string stream) do
          (string-case (code :default (if (sam-header-user-tag-p code)
                                        (push (parse-sam-comment stream) (getf user-tags (intern-key code)))
                                        (error "Unknown SAM record type code ~A." code)))
            ("@HD" (progn 
                     (assert (null hd))
                     (setq hd (parse-sam-header-line stream))))
            ("@SQ" (push (parse-sam-reference-sequence-dictionary-entry stream) sq))
            ("@RG" (push (parse-sam-read-group stream) rg))
            ("@PG" (push (parse-sam-program stream) pg))
            ("@CO" (push (parse-sam-comment stream) co)))
          finally (return (make-sam-header
                           :hd hd
                           :sq (nreverse sq)
                           :rg (nreverse rg)
                           :pg (nreverse pg)
                           :co (nreverse co)
                           :user-tags (loop for cons on user-tags by 'cddr
                                            do (setf (cdr cons) (nreverse (cdr cons)))
                                            finally (return user-tags)))))))

(define-symbol-macro optional-field-type-tags "AifZHB")
(defconstant +min-optional-field-type-tag+ (reduce 'min optional-field-type-tags :key 'char-code)
  "The smallest optional field type tag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.")
(defconstant +max-optional-field-type-tag+ (reduce 'max optional-field-type-tags :key 'char-code)
  "The largest optional field type tag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.")

(defun scan-error (stream)
  "Signal invalid type tag for optional field in SAM alignment."
  (declare (ignore stream))
  (error "Invalid type tag for optional field in SAM alignment."))

(defun make-optional-field-scan-table ()
  "Create a dispatch table for scanning optional fields in a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5."
  (let ((table (make-array (1+ (- +max-optional-field-type-tag+
                                  +min-optional-field-type-tag+))
                           :initial-element #'scan-error
                           :allocation :long-lived
                           :single-thread t)))
    (flet ((s (index value) (setf (svref table (- (char-code index) +min-optional-field-type-tag+)) value)))
      (s #\A #'readc)
      (s #\i #'scan-integer)
      (s #\f #'scan-float)
      (s #\Z #'scan-stringn)
      (s #\H #'parse-sam-byte-array)
      (s #\B #'parse-sam-numeric-array))
    table))

(declaim (notinline parse-sam-alignment))

(defun parse-sam-alignment (stream)
  "Parse a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (unless (end-of-file-p (peekc stream))
    (flet ((do-stringn ()
             (scan-stringn stream))
           (do-string ()
             (prog1 (scan-string stream)
               (scanc stream #\Tab)))
           (do-int32 ()
             (prog1 (scan-integer stream)
               (scanc stream #\Tab))))
      (declare (inline do-stringn do-string do-int32))
      (make-sam-alignment
       :qname (do-string)
       :flag  (do-int32)
       :rname (do-string)
       :pos   (do-int32)
       :mapq  (do-int32)
       :cigar (do-string)
       :rnext (do-string)
       :pnext (do-int32)
       :tlen  (do-int32)
       :seq   (do-string)
       :qual  (do-stringn)
       :tags  (loop until (end-of-line-p (peekc stream)) do
                    (scanc stream #\Tab) nconc
                    (let ((tag (intern-key (parse-sam-tag stream)))
                          (type (readc stream)))
                      (declare (symbol tag) (base-char type))
                      (scanc stream #\:)
                      (list tag (funcall (the function (svref (load-time-value (make-optional-field-scan-table) t)
                                                              (- (char-code type) +min-optional-field-type-tag+)))
                                         stream)))
                    finally (advance stream))))))

(defun parse-sam (stream)
  "Parse a complete SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1."
  (make-sam :header     (parse-sam-header stream)
            :alignments (loop for aln = (parse-sam-alignment stream)
                              while aln collect aln)))

;;; output

(declaim (inline %writec))

(defun %writec (out c)
  "Write a character to output stream."
  (declare (buffered-stream out) (base-char c) #.*fixnum-optimization*)
  (loop (with-stream-output-buffer (buffer index limit) out
          (declare (simple-base-string buffer) (fixnum index limit))
          (when (< index limit)
            (setf (sbchar buffer index) c)
            (incf index)
            (return-from %writec)))
        (stream-flush-buffer out)))

(declaim (inline writec write-tab write-newline))

(defun writec (out c)
  "Write a character to output stream."
  (declare (buffered-stream out) (base-char c) #.*optimization*)
  (%writec out c))

(defun write-tab (out)
  "Write a tabulator to output stream."
  (declare (buffered-stream out) #.*optimization*)
  (%writec out #\Tab))

(defun write-newline (out)
  "Write an end-of-line to output stream."
  (declare (buffered-stream out) #.*optimization*)
  (%writec out #\Newline))

(declaim (inline writestr))

(defun writestr (out s)
  "Write a string to output stream."
  (declare (buffered-stream out) (simple-base-string s) #.*fixnum-optimization*)
  (let ((sindex 0) (length (length s)))
    (declare (fixnum sindex length))
    (when (> length 0)
      (loop (with-stream-output-buffer (buffer tindex limit) out
              (declare (simple-base-string buffer) (fixnum tindex limit))
              (loop with rlen of-type fixnum = (min (- length sindex) (- limit tindex))
                    for i of-type fixnum below rlen
                    do (setf (sbchar buffer (+ tindex i))
                             (sbchar s (+ sindex i)))
                    finally (incf sindex rlen) (incf tindex rlen)))
            (if (< sindex length)
              (stream-flush-buffer out)
              (return-from writestr))))))

(defun format-sam-string (out tag string)
  "Write a SAM file TAG of type string."
  (declare (buffered-stream out) (simple-base-string tag string) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (writec out #\:)
  (writestr out string))

(defun format-sam-integer (out tag value)
  "Write a SAM file TAG of type integer."
  (declare (buffered-stream out) (simple-base-string tag) (integer value) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (format out ":~D" value))

(defun format-sam-byte-array (out tag byte-array)
  "Write a SAM file TAG of type byte array.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type H."
  (declare (buffered-stream out) (simple-base-string tag) (sequence byte-array) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (writec out #\:)
  (etypecase byte-array
    (list   (loop for byte in (the list byte-array)
                  do (format out "~2,'0X" byte)))
    (vector (loop for byte across (the vector byte-array)
                  do (format out "~2,'0X" byte)))))

(defun format-sam-datetime (out tag datetime)
  "Write a SAM file TAG of type date/time.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag @RG, DT."
  (declare (buffered-stream out) (simple-base-string tag) (integer datetime) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (multiple-value-bind
      (sec min hour day month year)
      (decode-universal-time datetime)
    (format out ":~4,'0D-~2,'0D-~2,'0DT~2,'0D:~2,'0D:~2,'0DZ"
            year month day hour min sec)))

(defun format-sam-header-user-tag (out tag value)
  "Write a user-defined SAM file TAG of type string."
  (declare (buffered-stream out) (symbol tag) #.*optimization*)
  (let ((tag-string (symbol-name tag)))
    (declare (simple-base-string tag-string))
    (if (sam-header-user-tag-p tag-string)
      (format-sam-string out tag-string value)
      (error "Unknown tag ~A in a SAM header line." tag))))

(defun format-sam-header-line (out list)
  "Write an @HD line in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list list) #.*optimization*)
  (when list
    (unless (presentp :VN list)
      (cerror "Ignore absence of VN tag and continue."
              "VN tag missing in @HD line when writing ~A." out))
    (writestr out "@HD")
    (loop for (tag value) of-type (symbol t) on list by 'cddr do
          (case tag 
            (:VN (format-sam-string out "VN" value))
            (:SO (format-sam-string out "SO" value))
            (t   (format-sam-header-user-tag out tag value))))
    (write-newline out)))

(defun format-sam-reference-sequence-dictionary (out list)
  "Write the @SQ lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :SN plist)
          (cerror "Ignore absence of SN tag and continue."
                  "SN tag missing in @SQ line when writing ~A." out))
        (unless (presentp :LN plist)
          (cerror "Ignore absence of LN tag and continue."
                  "LN tag missing in @SQ line when writing ~A." out))
        (writestr out "@SQ")
        (loop for (tag value) on plist by 'cddr do
              (case tag
                (:SN (format-sam-string out "SN" value))
                (:LN (format-sam-integer out "LN" value))
                (:AS (format-sam-string out "AS" value))
                (:M5 (format-sam-byte-array out "M5" value))
                (:SP (format-sam-string out "SP" value))
                (:UR (format-sam-string out "UR" value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-read-groups (out list)
  "Write the @RG lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :ID plist)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @RG line when writing ~A." out))
        (writestr out "@RG")
        (loop for (tag value) on plist by 'cddr do
              (case tag
                (:ID (format-sam-string out "ID" value))
                (:CN (format-sam-string out "CN" value))
                (:DS (format-sam-string out "DS" value))
                (:DT (format-sam-datetime out "DT" value))
                (:FO (format-sam-string out "FO" value))
                (:KS (format-sam-string out "KS" value))
                (:LB (format-sam-string out "LB" value))
                (:PG (format-sam-string out "PG" value))
                (:PI (format-sam-integer out "PI" value))
                (:PL (format-sam-string out "PL" value))
                (:PU (format-sam-string out "PU" value))
                (:SM (format-sam-string out "SM" value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-programs (out list)
  "Write the @PG lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :ID plist)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @PG line when writing ~A." out))
        (writestr out "@PG")
        (loop for (tag value) on plist by 'cddr do
              (case tag
                (:ID (format-sam-string out "ID" value))
                (:PN (format-sam-string out "PN" value))
                (:CL (format-sam-string out "CL" value))
                (:PP (format-sam-string out "PP" value))
                (:DS (format-sam-string out "DS" value))
                (:VN (format-sam-string out "VN" value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-comments (out list)
  "Write the @CO lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list list) #.*optimization*)
  (loop for string in list do
        (writestr out "@CO")
        (write-tab out)
        (writestr out string)
        (write-newline out)))

(defun format-sam-user-tags (out tags)
  "Write the user-defined header lines.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (list tags) #.*optimization*)
  (loop for (code list) of-type (symbol list) on tags by 'cddr
        for code-string = (symbol-name code) do
        (loop for string in list do
              (writestr out code-string)
              (write-tab out)
              (writestr out string)
              (write-newline out))))

(defun format-sam-header (out header)
  "Write a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (buffered-stream out) (sam-header header) #.*optimization*)
  (format-sam-header-line                   out (sam-header-hd header))
  (format-sam-reference-sequence-dictionary out (sam-header-sq header))
  (format-sam-read-groups                   out (sam-header-rg header))
  (format-sam-programs                      out (sam-header-pg header))
  (format-sam-comments                      out (sam-header-co header))
  (format-sam-user-tags                     out (sam-header-user-tags header)))

(eval-when (:compile-toplevel :execute)
  (define-symbol-macro   int8 #\c)
  (define-symbol-macro  uint8 #\C)
  (define-symbol-macro  int16 #\s)
  (define-symbol-macro uint16 #\S)
  (define-symbol-macro  int32 #\i)
  (define-symbol-macro uint32 #\I)
  (define-symbol-macro    num #\f)

  (defvar *integer-types* 
    (list int8 uint8 int16 uint16 int32 uint32)
    "Integer types supported in SAM file numeric arrays.
     See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type B.")

  (defun make-integer-type-descriptors (&rest types)
    "Create descriptors for the integer types in SAM file numeric arrays.
     Types with smaller indexes are preferred over types with larger indexes.
     See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type B."
    (declare (dynamic-extent types))
    (loop for type in *integer-types*
          for index from 0
          when (member type types :test 'char=) 
          collect (cons type index))))

(declaim (inline common-number-type))

(defun common-number-type (numbers)
  "Determine the smallest numeric type supported by SAM file numeric arrays that can represent all numbers in the given list.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type B."
  (declare (list numbers) #.*optimization*)
  (flet ((integer-types (number)
           (declare (number number))
           (etypecase number
             (float '())
             (integer (cond ((< number  #.(- (expt 2 31))) '())
                            ((< number  #.(- (expt 2 15))) '#.(make-integer-type-descriptors int32))
                            ((< number  #.(- (expt 2  7))) '#.(make-integer-type-descriptors int16 int32))
                            ((< number 0)                  '#.(make-integer-type-descriptors int8 int16 int32))
                            ((> number #.(1- (expt 2 32))) '())
                            ((> number #.(1- (expt 2 31))) '#.(make-integer-type-descriptors uint32))
                            ((> number #.(1- (expt 2 16))) '#.(make-integer-type-descriptors int32 uint32))
                            ((> number #.(1- (expt 2 15))) '#.(make-integer-type-descriptors uint16 int32 uint32))
                            ((> number #.(1- (expt 2  8))) '#.(make-integer-type-descriptors int16 uint16 int32 uint32))
                            ((> number #.(1- (expt 2  7))) '#.(make-integer-type-descriptors uint8 int16 uint16 int32 uint32))
                            (t                             '#.(make-integer-type-descriptors int8 uint8 int16 uint16 int32 uint32)))))))
    (loop for number in numbers
          for types = (copy-list (integer-types number))
          then (nintersection types (integer-types number) :key 'car :test 'char=)
          unless types return #\f
          finally (return (caar (stable-sort types '< :key 'cdr))))))

(defun format-sam-tag (out tag value)
  "Write a SAM file TAG, dispatching on actual type of the given value.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5."
  (declare (buffered-stream out) (symbol tag) #.*optimization*)
  (write-tab out)
  (writestr out (symbol-name tag))
  (etypecase value
    (character    (writestr out ":A:") (writec out value))           ; printable character
    (int32        (format out ":i:~D" value))                        ; signed 32-bit integer
    (integer      (format out ":f:~D" value))                        ; single-precision floating point number
    (single-float (format out ":f:~E" value))                        ; single-precision floating point number
    (double-float (format out ":f:~E" (coerce value 'single-float))) ; single-precision floating point number
    (string       (writestr out ":Z:") (writestr out value))         ; printable string
    (vector       (writestr out ":H:")                               ; byte array in hex format
                  (loop for byte across (the vector value)
                        do (format out "~2,'0X" byte)))
    (list         (writestr out ":B:")                               ; integer or numeric array
                  (let ((type (common-number-type value)))
                    (writec out type)
                    (if (char= type num)
                      (loop for number in value do
                            (if (integerp number) (format out ",~D" number)
                              (format out ",~E" (coerce number 'single-float))))
                      (loop for number in value do (format out ",~D" number)))))))

(defun format-sam-alignment (out aln)
  "Write a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and 1.5."
  (declare (buffered-stream out) (sam-alignment aln) #.*optimization*)
  (writestr out    (sam-alignment-qname aln)) (write-tab out)
  (format out "~D" (sam-alignment-flag aln))  (write-tab out)
  (writestr out    (sam-alignment-rname aln)) (write-tab out)
  (format out "~D" (sam-alignment-pos aln))   (write-tab out)
  (format out "~D" (sam-alignment-mapq aln))  (write-tab out)
  (writestr out    (sam-alignment-cigar aln)) (write-tab out)
  (writestr out    (sam-alignment-rnext aln)) (write-tab out)
  (format out "~D" (sam-alignment-pnext aln)) (write-tab out)
  (format out "~D" (sam-alignment-tlen aln))  (write-tab out)
  (writestr out    (sam-alignment-seq aln))   (write-tab out)
  (writestr out    (sam-alignment-qual aln))
  (loop for (key value) on (sam-alignment-tags aln) by 'cddr
        do (format-sam-tag out key value))
  (write-newline out))

(defun format-sam (out sam)
  "Write a complete SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1."
  (declare (buffered-stream out) (sam sam) #.*optimization*)
  (format-sam-header out (sam-header sam))
  (loop for aln in (sam-alignments sam)
        do (format-sam-alignment out aln)))


(hcl:defglobal-variable *samtools* nil
  "Location of the samtools binary.")

(defun get-samtools ()
  "Determine location of the samtools binary."
  (or *samtools* (with-open-stream (stream (sys:open-pipe "command -v samtools"))
                   (multiple-value-bind (line eof) (stream-read-line stream)
                     (if eof
                       (error "samtools not found. Please download it from http://samtools.sourceforge.net and make sure that its binary is present in your PATH.")
                       (setq *samtools* line))))))

(defun check-stdout (pathname)
  (ignore-errors
    (string= (namestring (truename pathname))
             (load-time-value (namestring (truename #p"/dev/stdout"))))))

(defun check-stdin (pathname)
  (ignore-errors
    (string= (namestring (truename pathname))
             (load-time-value (namestring (truename #p"/dev/stdin"))))))

(defvar *reference-fasta* nil)
(defvar *reference-fai* nil)

(defun %open-sam (pathname direction header-only kind)
  "Open a SAM file for :input or :output. If kind is :bam or :cram, use samtools view for input or output.
   Tell samtools view to only return the header for :input when :header-only is true."
  (ecase direction
    (:input (cond ((eq kind :sam)
                   (if (check-stdin pathname)
                       (make-synonym-stream '*terminal-io*)
                     (open pathname :direction :input :element-type 'base-char :if-does-not-exist :error)))
                  (t (open pathname :direction :probe :element-type 'base-char :if-does-not-exist :error)
                     (sys:open-pipe (list (get-samtools) "view" (if header-only "-H" "-h") "-@" (write-to-string *number-of-threads*)
                                          (namestring (translate-logical-pathname pathname)))
                                    :direction :input))))
    (:output (cond ((eq kind :sam)
                    (if (check-stdout pathname)
                        (make-synonym-stream '*terminal-io*)
                      (open pathname :direction :output :element-type 'base-char :if-exists :supersede)))
                   ((eq kind :cram)
                    (let ((reference-fasta *reference-fasta*)
                          (reference-fai *reference-fai*))
                      (cond (reference-fasta
                             (sys:open-pipe (list (get-samtools) "view" "-C" "-@" (write-to-string *number-of-threads*) "-T" reference-fasta
                                                  "-o" (namestring (translate-logical-pathname pathname)) "-") :direction :output))
                            (reference-fai
                             (sys:open-pipe (list (get-samtools) "view" "-C" "-@" (write-to-string *number-of-threads*) "-t" reference-fai
                                                  "-o" (namestring (translate-logical-pathname pathname)) "-") :direction :output)))))
                   (t (sys:open-pipe (list (get-samtools) "view" (ecase kind (:bam "-Sb") (:cram "-C")) "-@" (write-to-string *number-of-threads*)
                                           "-o" (namestring (translate-logical-pathname pathname)) "-")
                                     :direction :output))))))

(defun sam-file-kind (filename)
  "Determine whether the file is :bam for .bam, :cram for .cram, or :sam in all other cases."
  (let ((type (pathname-type (pathname filename))))
    (etypecase type
      (string (cond ((string-equal "bam" type) :bam)
                    ((string-equal "cram" type) :cram)
                    (t :sam)))
      (symbol :sam))))

(declaim (inline open-sam))

(defun open-sam (filename &key (direction :input) header-only)
  "Open a SAM file for :input or :output. If the file is .bam or .cram, use samtools view for input or output.
   Tell samtools view to only return the header for :input when :header-only is true."
  (%open-sam filename direction header-only (sam-file-kind filename)))

(defun open-temporary-sam (sibling)
  "Open a temporary SAM file for :output in the same folder as the sibling file. If the file is .bam or .cram, use samtools view for output."
  (let ((location (lw:pathname-location (pathname sibling)))
        (kind (sam-file-kind sibling)))
    (cond ((eq kind :sam)
           (let ((stream (hcl:open-temp-file :directory location)))
             (values stream (pathname stream))))
          (t (let ((pathname (hcl:create-temp-file :directory location)))
               (assert (delete-file pathname))
               (values (%open-sam pathname :output nil kind) pathname))))))
