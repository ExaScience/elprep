(in-package :elprep)

;;; input

(declaim (inline make-scanner))

(defstruct (scanner (:constructor make-scanner (&optional (string empty-sbs)))
                    (:copier nil)
                    (:predicate nil))
  (index 0 :type fixnum)
  (string empty-sbs :type simple-base-string))

(setf (documentation 'make-scanner 'function)
      "Create a scanner to scan/parse strings of type simple-base-string."
      (documentation 'scanner-index 'function)
      "The current index into the string while scanning/parsing it."
      (documentation 'scanner-string 'function)
      "The string being scanned/parsed.")

(declaim (inline reinitialize-scanner))
 
(defun reinitialize-scanner (scanner string)
  "Reinitialize the scanner so it can be reused."
  (declare (scanner scanner) (simple-base-string string))
  (setf (scanner-index scanner) 0
        (scanner-string scanner) string)
  scanner)

(declaim (inline advance peekc readc scanc end-of-entry-p))

(defun advance (scanner)
  "Advance scanner by one character if possible."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (when (< index (length string))
      (setf (scanner-index scanner) (the fixnum (1+ index))))))

(defun peekc (scanner)
  "Peek current character in scanner if possible, or return #\EOT."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (if (< index (length string))
      (schar string index)
      #\EOT)))

(defun readc (scanner)
  "Read current character in scanner if possible, or return #\EOT."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (if (< index (length string))
      (prog1 (schar string index)
        (setf (scanner-index scanner) (the fixnum (1+ index))))
      #\EOT)))

(defun scanc (scanner char)
  "Read current character in scanner, ensuring it is the given one."
  (declare (scanner scanner) (base-char char) #.*optimization*)
  (assert (char= (readc scanner) char))
  char)

(defun end-of-entry-p (char)
  "Is the character a tab or end of line, delimiting an entry in a tab-delimited text file?"
  (declare (base-char) #.*optimization*)
  (or (char= char #\Tab)
      (char= char #\EOT)))

(defun scan-string (scanner)
  "Read a tab-delimited entry as string from scanner."
  (declare (scanner scanner) #.*optimization*)
  (loop with index of-type fixnum = (scanner-index scanner)
        with string of-type simple-base-string = (scanner-string scanner)
        for end of-type fixnum from index below (length string)
        until (char= (schar string end) #\Tab) finally
        (let ((result (make-array (the fixnum (- end index)) :element-type 'base-char)))
          (declare (simple-base-string result))
          (loop for j of-type fixnum from index below end
                for i of-type fixnum from 0
                do (setf (schar result i) (schar string j)))
          (setf (scanner-index scanner) end)
          (return result))))

(defun scan-integer (scanner)
  "Scan an integer from scanner."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (let* ((pos index) (char (schar string pos)))
      (declare (fixnum pos) (base-char char))
      (flet ((nextc () (if (= (incf pos) (length string))
                         (setq char #\EOT)
                         (setq char (schar string pos)))))
        (declare (inline nextc))
        (when (or (char= char #\-)
                  (char= char #\+))
          (nextc))
        (assert (and (char<= #\0 char) (char<= char #\9)))
        (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9))))
      (setf (scanner-index scanner) pos)
      (parse-integer string :start index :end pos))))

(defun scan-float (scanner)
  "Scan a floating point number from scanner. Returns either an integer or a single-float."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (let* ((pos index) (char (schar string pos)))
      (declare (fixnum pos) (base-char char))
      (flet ((nextc () (if (= (incf pos) (length string))
                         (setq char #\EOT)
                         (setq char (schar string pos)))))
        (declare (inline nextc))
        (when (or (char= char #\-)
                  (char= char #\+))
          (nextc))
        (cond ((and (char<= #\0 char) (char<= char #\9))
               (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9)))
               (when (char= char #\.)
                 (nextc)
                 (assert (and (char<= #\0 char) (char<= char #\9)))
                 (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9)))))
              (t (when (char= char #\.)
                   (nextc))
                 (assert (and (char<= #\0 char) (char<= char #\9)))
                 (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9)))))
        (when (or (char= char #\e)
                  (char= char #\E))
          (nextc)
          (when (or (char= char #\-)
                    (char= char #\+))
            (nextc))
          (assert (and (char<= #\0 char) (char<= char #\9)))
          (loop do (nextc) while (and (char<= #\0 char) (char<= char #\9)))))
      (setf (scanner-index scanner) pos)
      (with-input-from-string (stream string :start index :end pos)
        (let ((*read-base* 10) (*read-default-float-format* 'single-float) (*read-eval* nil))
          (let ((value (read stream)))
            (etypecase value
              (integer value)
              (single-float value)
              (number (coerce value 'single-float)))))))))

(defun parse-sam-tag (scanner &optional (tag-string (make-array 2 :element-type 'base-char)))
  "Parse the TAG: portion of a SAM file tag."
  (declare (scanner scanner) (simple-base-string tag-string) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (assert (<= index (- (length string) 3)))
    (setf (schar tag-string 0) (schar string index))
    (setf (schar tag-string 1) (schar string (the fixnum (+ index 1))))
    (assert (char= (schar string (the fixnum (+ index 2))) #\:))
    (setf (scanner-index scanner) (the fixnum (+ index 3)))
    tag-string))

(defun parse-sam-byte-array (scanner)
  "Parse a byte array in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type H."
  (declare (scanner scanner) #.*fixnum-optimization*)
  (flet ((hex-value (char)
           (declare (base-char char))
           (cond ((and (char<= #\0 char) (char<= char #\9))
                  (- (char-int char) #.(char-int #\0)))
                 ((and (char<= #\a char) (char<= char #\f))
                  (- (char-int char) #.(- (char-int #\a) 10)))
                 ((and (char<= #\A char) (char<= char #\F))
                  (- (char-int char) #.(- (char-int #\A) 10)))
                 (t (error "Not a hex digit ~S in ~S." char scanner)))))
    (declare (inline hex-value))
    (loop for count of-type fixnum from 1
          collect (+ (ash (hex-value (readc scanner)) 4)
                     (hex-value (readc scanner)))
          into list until (end-of-entry-p (peekc scanner))
          finally (return (make-array count :element-type 'octet :initial-contents list)))))

(defun parse-sam-numeric-array (scanner)
  "Parse a numeric array in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type B."
  (declare (scanner scanner) #.*optimization*)
  (advance scanner) ;; ignore type tag
  (loop do (scanc scanner #\,)
        collect (scan-float scanner)
        until (end-of-entry-p (peekc scanner))))

(defun parse-sam-header-line (scanner)
  "Parse an @HD line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (scanner scanner) #.*optimization*)
  (let ((tag (make-array 2 :element-type 'base-char)))
    (declare (simple-base-string tag) (dynamic-extent tag))
    (loop until (char= (peekc scanner) #\EOT)
          nconc (progn
                  (scanc scanner #\Tab)
                  (parse-sam-tag scanner tag)
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key/copy tag) record) (scan-string scanner))
                                               (error "Unknown tag ~A in SAM header line ~S." tag scanner)))
                    (#.(sbs "VN") (list (unique :VN record) (scan-string scanner)))
                    (#.(sbs "SO") (list (unique :SO record) (scan-string scanner)))))
          into record finally
          (advance scanner)
          (unless (presentp :VN record)
            (cerror "Ignore absence of VN tag and continue."
                    "VN tag missing in @HD line when reading ~S." scanner))
          (return record))))

(defun parse-sam-reference-sequence-dictionary-entry (scanner)
  "Parse an @SQ line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (scanner scanner) #.*optimization*)
  (let ((tag (make-array 2 :element-type 'base-char)))
    (declare (simple-base-string tag) (dynamic-extent tag))
    (loop until (char= (peekc scanner) #\EOT)
          nconc (progn
                  (scanc scanner #\Tab)
                  (parse-sam-tag scanner tag)
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key/copy tag) record) (scan-string scanner))
                                               (error "Unknown tag ~A in SAM reference sequence dictionary line ~S." tag scanner)))
                    (#.(sbs "SN") (list (unique :SN record) (scan-string scanner)))
                    (#.(sbs "LN") (list (unique :LN record) (scan-integer scanner)))
                    (#.(sbs "AS") (list (unique :AS record) (scan-string scanner)))
                    (#.(sbs "M5") (list (unique :M5 record) (parse-sam-byte-array scanner)))
                    (#.(sbs "SP") (list (unique :SP record) (scan-string scanner)))
                    (#.(sbs "UR") (list (unique :UR record) (scan-string scanner)))))
          into record finally
          (advance scanner)
          (unless (presentp :SN record)
            (cerror "Ignore absence of SN tag and continue."
                    "SN tag missing in @SQ line when reading ~S." scanner))
          (unless (presentp :LN record)
            (cerror "Ignore absence of LN tag and continue."
                    "LN tag missing in @SQ line when reading ~S." scanner))
          (return record))))

(defun parse-sam-read-group (scanner)
  "Parse an @RG line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (scanner scanner) #.*optimization*)
  (let ((tag (make-array 2 :element-type 'base-char)))
    (declare (simple-base-string tag) (dynamic-extent tag))
    (loop until (char= (peekc scanner) #\EOT)
          nconc (progn
                  (scanc scanner #\Tab)
                  (parse-sam-tag scanner tag)
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key/copy tag) record) (scan-string scanner))
                                               (error "Unknown tag ~A in SAM read group line ~S." tag scanner)))
                    (#.(sbs "ID") (list (unique :ID record) (scan-string scanner)))
                    (#.(sbs "CN") (list (unique :CN record) (scan-string scanner)))
                    (#.(sbs "DS") (list (unique :DS record) (scan-string scanner)))
                    (#.(sbs "DT") (list (unique :DT record) (parse-date-time (scan-string scanner))))
                    (#.(sbs "FO") (list (unique :FO record) (scan-string scanner)))
                    (#.(sbs "KS") (list (unique :KS record) (scan-string scanner)))
                    (#.(sbs "LB") (list (unique :LB record) (scan-string scanner)))
                    (#.(sbs "PG") (list (unique :PG record) (scan-string scanner)))
                    (#.(sbs "PI") (list (unique :PI record) (scan-integer scanner)))
                    (#.(sbs "PL") (list (unique :PL record) (scan-string scanner)))
                    (#.(sbs "PU") (list (unique :PU record) (scan-string scanner)))
                    (#.(sbs "SM") (list (unique :SM record) (scan-string scanner)))))
          into record finally
          (advance scanner)
          (unless (presentp :ID record)
            (cerror "Ignore absence of ID tag and continue."
                    "ID tag missing in @RG line when reading ~S." scanner))
          (return record))))

(defun parse-sam-program (scanner)
  "Parse an @PG line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (scanner scanner) #.*optimization*)
  (let ((tag (make-array 2 :element-type 'base-char)))
    (declare (simple-base-string tag) (dynamic-extent tag))
    (loop until (char= (peekc scanner) #\EOT)
          nconc (progn
                  (scanc scanner #\Tab)
                  (parse-sam-tag scanner tag)
                  (string-case (tag :default (if (sam-header-user-tag-p tag)
                                               (list (unique (intern-key/copy tag) record) (scan-string scanner))
                                               (error "Unknown tag ~A in SAM program line ~S." tag scanner)))
                    (#.(sbs "ID") (list (unique :ID record) (scan-string scanner)))
                    (#.(sbs "PN") (list (unique :PN record) (scan-string scanner)))
                    (#.(sbs "CL") (list (unique :CL record) (scan-string scanner)))
                    (#.(sbs "PP") (list (unique :PP record) (scan-string scanner)))
                    (#.(sbs "DS") (list (unique :DS record) (scan-string scanner)))
                    (#.(sbs "VN") (list (unique :VN record) (scan-string scanner)))))
          into record finally
          (advance scanner)
          (unless (presentp :ID record)
            (cerror "Ignore absence of ID tag and continue."
                    "ID tag missing in @PG line when reading ~S." scanner))
          (return record))))

(defun parse-sam-comment (scanner)
  "Parse an @CO line in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (scanner scanner) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (when (char= (schar string index) #\Tab)
      (incf index))
    (let ((result (make-array (- (length string) index) :element-type 'base-char)))
      (declare (simple-base-string result))
      (loop for j of-type fixnum from index below (length string)
            for i of-type fixnum from 0
            do (setf (schar result i) (schar string j)))
      (setf (scanner-index scanner) (length string))
      result)))

(defun parse-sam-header-code (scanner &optional (code-string (make-array 3 :element-type 'base-char)))
  "Parse the record type code of a SAM header line."
  (declare (scanner scanner) (simple-base-string code-string) #.*optimization*)
  (let ((index (scanner-index scanner))
        (string (scanner-string scanner)))
    (declare (fixnum index) (simple-base-string string))
    (assert (<= index (- (length string) 3)))
    (setf (schar code-string 0) (schar string index))
    (setf (schar code-string 1) (schar string (the fixnum (+ index 1))))
    (setf (schar code-string 2) (schar string (the fixnum (+ index 2))))
    (setf (scanner-index scanner) (the fixnum (+ index 3)))
    code-string))

(defun parse-sam-header (stream)
  "Parse a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream stream) #.*optimization*)
  (let ((scanner (make-scanner))
        (code (make-array 3 :element-type 'base-char))
        hd sq rg pg co user-tags)
    (declare (scanner scanner) (simple-base-string code) (dynamic-extent scanner code))
    (loop until (char/= (the base-char (ascii-stream-peek-char stream #\EOT)) #\@) do
          (reinitialize-scanner scanner (ascii-stream-read-line stream))
          (parse-sam-header-code scanner code)
          (string-case (code :default (if (sam-header-user-tag-p code)
                                        (push (parse-sam-comment scanner) (getf user-tags (intern-key/copy code)))
                                        (error "Unknown SAM record type code ~A in header line ~S." code scanner)))
            (#.(sbs "@HD") (progn
                             (assert (null hd))
                             (setq hd (parse-sam-header-line scanner))))
            (#.(sbs "@SQ") (push (parse-sam-reference-sequence-dictionary-entry scanner) sq))
            (#.(sbs "@RG") (push (parse-sam-read-group scanner) rg))
            (#.(sbs "@PG") (push (parse-sam-program scanner) pg))
            (#.(sbs "@CO") (push (parse-sam-comment scanner) co)))
          finally (return (make-sam-header
                           :hd hd
                           :sq (nreverse sq)
                           :rg (nreverse rg)
                           :pg (nreverse pg)
                           :co (nreverse co)
                           :user-tags (loop for cons on user-tags by #'cddr
                                            do (setf (cdr cons) (nreverse (cdr cons)))
                                            finally (return user-tags)))))))

(declaim (inline %skip-line))

#+lispworks
(defun %skip-line (stream)
  "Skip characters from stream until a newline is reached."
  (declare (buffered-stream stream) #.*fixnum-optimization*)
  (loop (with-stream-input-buffer (source index limit) stream
          (declare (simple-base-string source) (fixnum index limit))
          (loop for pos of-type fixnum from index below limit do
                (when (char= (lw:sbchar source pos) #\Newline)
                  (setq index (1+ pos))
                  (return-from %skip-line))
                finally (setq index pos)))
        (unless (stream-fill-buffer stream)
          (return-from %skip-line))))

#+sbcl
(defun %skip-line (stream)
  "Skip characters from stream until a newline is reached."
  (declare (ascii-stream stream) #.*optimization*)
  (let ((buffer (ascii-stream-buffer stream)))
    (declare (ascii-stream-buffer buffer))
    (with-buffer-dispatch buffer
      (loop (let ((index (ascii-stream-index stream))
                  (limit (ascii-stream-limit stream)))
              (declare (fixnum index limit))
              (loop for pos of-type fixnum from index below limit do
                    (when (char= (bchar buffer pos) #\Newline)
                      (setf (ascii-stream-index stream) (the fixnum (1+ pos)))
                      (return-from %skip-line))
                    finally (setf (ascii-stream-index stream) pos)))
            (unless (ascii-stream-fill-buffer stream)
              (return-from %skip-line))))))

(defun skip-sam-header (stream)
  "Skip the SAM file header section."
  (declare #.*optimization*)
  (loop until (char/= (the base-char (ascii-stream-peek-char stream #\EOT)) #\@)
        do (%skip-line stream))
  (values))

(define-symbol-macro optional-field-type-tags (sbs "AifZHB"))

(defconstant +min-optional-field-type-tag+ (reduce #'min optional-field-type-tags :key #'char-code)
  "The smallest optional field type tag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.")

(defconstant +max-optional-field-type-tag+ (reduce #'max optional-field-type-tags :key #'char-code)
  "The largest optional field type tag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.")

(defun scan-error (scanner)
  "Signal invalid type tag for optional field in SAM alignment."
  (declare (scanner scanner))
  (error "Invalid type tag for optional field in SAM alignment ~S." scanner))

(defun make-optional-field-scan-table ()
  "Create a dispatch table for scanning optional fields in a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5."
  (let ((table (make-array (1+ (- +max-optional-field-type-tag+
                                  +min-optional-field-type-tag+))
                           :initial-element #'scan-error
                           #+lispworks :allocation #+lispworks :long-lived
                           #+lispworks :single-thread #+lispworks t)))
    (flet ((s (index value) (setf (svref table (- (char-code index) +min-optional-field-type-tag+)) value)))
      (s #\A #'readc)
      (s #\i #'scan-integer)
      (s #\f #'scan-float)
      (s #\Z #'scan-string)
      (s #\H #'parse-sam-byte-array)
      (s #\B #'parse-sam-numeric-array))
    table))

(declaim (notinline parse-sam-alignment))

(defun parse-sam-alignment (string)
  "Parse a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4."
  (declare (simple-base-string string) #.*optimization*)
  (let ((scanner (make-scanner string)))
    (declare (scanner scanner) (dynamic-extent scanner))
    (flet ((do-stringn ()
             (scan-string scanner))
           (do-string ()
             (prog1 (scan-string scanner)
               (scanc scanner #\Tab)))
           (do-int32 ()
             (prog1 (scan-integer scanner)
               (scanc scanner #\Tab))))
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
       :tags  (let ((tag-string (make-array 2 :element-type 'base-char)))
                (declare (simple-base-string tag-string) (dynamic-extent tag-string))
                (loop until (char= (peekc scanner) #\EOT)
                      nconc (progn
                              (scanc scanner #\Tab)
                              (let ((tag (intern-key/copy (parse-sam-tag scanner tag-string)))
                                    (type (readc scanner)))
                                (declare (symbol tag) (base-char type))
                                (scanc scanner #\:)
                                (list tag (funcall (the function (svref (load-time-value (make-optional-field-scan-table) t)
                                                                        (- (char-code type) +min-optional-field-type-tag+)))
                                                   scanner))))))))))

(defun parse-sam (stream)
  "Parse a complete SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1."
  (make-sam :header     (parse-sam-header stream)
            :alignments (loop while (ascii-stream-listen stream)
                              collect (parse-sam-alignment (ascii-stream-read-line stream)))))


;;; output

(declaim (inline writec write-newline write-tab))

(defun writec (out c)
  "Write a character to output stream."
  (ascii-stream-write-char out c))

(defun write-newline (out)
  "Write a newline to output stream."
  (writec out #\Newline))

(defun write-tab (out)
  "Write a tabulator to output stream."
  (writec out #\Tab))

(declaim (inline writestr))

(defun writestr (out string)
  (ascii-stream-write-string out string))

(defun format-sam-string (out tag string)
  "Write a SAM file TAG of type string."
  (declare (stream out) (string tag string) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (writec out #\:)
  (writestr out string))

(defun format-sam-integer (out tag value)
  "Write a SAM file TAG of type integer."
  (declare (stream out) (simple-base-string tag) (integer value) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (format out (sbs ":~D") value))

(defun format-sam-byte-array (out tag byte-array)
  "Write a SAM file TAG of type byte array.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5, Type H."
  (declare (stream out) (simple-base-string tag) (sequence byte-array) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (writec out #\:)
  (etypecase byte-array
    (list   (loop for byte in (the list byte-array)
                  do (format out (sbs "~2,'0X") byte)))
    (vector (loop for byte across (the vector byte-array)
                  do (format out (sbs "~2,'0X") byte)))))

(defun format-sam-datetime (out tag datetime)
  "Write a SAM file TAG of type date/time.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag @RG, DT."
  (declare (stream out) (simple-base-string tag) (integer datetime) #.*optimization*)
  (write-tab out)
  (writestr out tag)
  (multiple-value-bind
      (sec min hour day month year)
      (decode-universal-time datetime)
    (format out (sbs ":~4,'0D-~2,'0D-~2,'0DT~2,'0D:~2,'0D:~2,'0DZ")
            year month day hour min sec)))

(defun format-sam-header-user-tag (out tag value)
  "Write a user-defined SAM file TAG of type string."
  (declare (stream out) (symbol tag) #.*optimization*)
  (let ((tag-string (symbol-name tag)))
    (if (sam-header-user-tag-p tag-string)
      (format-sam-string out tag-string value)
      (error "Unknown tag ~A in a SAM header line." tag))))

(defun format-sam-header-line (out list)
  "Write an @HD line in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list list) #.*optimization*)
  (when list
    (unless (presentp :VN list)
      (cerror "Ignore absence of VN tag and continue."
              "VN tag missing in @HD line when writing ~A." out))
    (writestr out (sbs "@HD"))
    (loop for (tag value) of-type (symbol t) on list by #'cddr do
          (case tag 
            (:VN (format-sam-string out (sbs "VN") value))
            (:SO (format-sam-string out (sbs "SO") value))
            (t   (format-sam-header-user-tag out tag value))))
    (write-newline out)))

(defun format-sam-reference-sequence-dictionary (out list)
  "Write the @SQ lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :SN plist)
          (cerror "Ignore absence of SN tag and continue."
                  "SN tag missing in @SQ line when writing ~A." out))
        (unless (presentp :LN plist)
          (cerror "Ignore absence of LN tag and continue."
                  "LN tag missing in @SQ line when writing ~A." out))
        (writestr out (sbs "@SQ"))
        (loop for (tag value) on plist by #'cddr do
              (case tag
                (:SN (format-sam-string out (sbs "SN") value))
                (:LN (format-sam-integer out (sbs "LN") value))
                (:AS (format-sam-string out (sbs "AS") value))
                (:M5 (format-sam-byte-array out (sbs "M5") value))
                (:SP (format-sam-string out (sbs "SP") value))
                (:UR (format-sam-string out (sbs "UR") value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-read-groups (out list)
  "Write the @RG lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :ID plist)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @RG line when writing ~A." out))
        (writestr out (sbs "@RG"))
        (loop for (tag value) on plist by #'cddr do
              (case tag
                (:ID (format-sam-string out (sbs "ID") value))
                (:CN (format-sam-string out (sbs "CN") value))
                (:DS (format-sam-string out (sbs "DS") value))
                (:DT (format-sam-datetime out (sbs "DT") value))
                (:FO (format-sam-string out (sbs "FO") value))
                (:KS (format-sam-string out (sbs "KS") value))
                (:LB (format-sam-string out (sbs "LB") value))
                (:PG (format-sam-string out (sbs "PG") value))
                (:PI (format-sam-integer out (sbs "PI") value))
                (:PL (format-sam-string out (sbs "PL") value))
                (:PU (format-sam-string out (sbs "PU") value))
                (:SM (format-sam-string out (sbs "SM") value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-programs (out list)
  "Write the @PG lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list list) #.*optimization*)
  (loop for plist in list do
        (unless (presentp :ID plist)
          (cerror "Ignore absence of ID tag and continue."
                  "ID tag missing in @PG line when writing ~A." out))
        (writestr out (sbs "@PG"))
        (loop for (tag value) on plist by #'cddr do
              (case tag
                (:ID (format-sam-string out (sbs "ID") value))
                (:PN (format-sam-string out (sbs "PN") value))
                (:CL (format-sam-string out (sbs "CL") value))
                (:PP (format-sam-string out (sbs "PP") value))
                (:DS (format-sam-string out (sbs "DS") value))
                (:VN (format-sam-string out (sbs "VN") value))
                (t   (format-sam-header-user-tag out tag value))))
        (write-newline out)))

(defun format-sam-comments (out list)
  "Write the @CO lines in a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list list) #.*optimization*)
  (loop for string in list do
        (writestr out (sbs "@CO"))
        (write-tab out)
        (writestr out string)
        (write-newline out)))

(defun format-sam-user-tags (out tags)
  "Write the user-defined header lines.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (list tags) #.*optimization*)
  (loop for (code list) of-type (symbol list) on tags by #'cddr
        for code-string = (symbol-name code) do
        (loop for string in list do
              (writestr out code-string)
              (write-tab out)
              (writestr out string)
              (write-newline out))))

(defun format-sam-header (out header)
  "Write a SAM file header section.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3."
  (declare (stream out) (sam-header header) #.*optimization*)
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
          when (member type types :test #'char=) 
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
             (fixnum  (locally (declare (fixnum number))
                        (cond ((< number  #.(- (expt 2 31))) '())
                              ((> number #.(1- (expt 2 32))) '())
                              ((> number #.(1- (expt 2 31))) '#.(make-integer-type-descriptors uint32))

                              ((< number  #.(- (expt 2 15))) '#.(make-integer-type-descriptors int32))
                              ((> number #.(1- (expt 2 16))) '#.(make-integer-type-descriptors int32 uint32))
                              ((> number #.(1- (expt 2 15))) '#.(make-integer-type-descriptors uint16 int32 uint32))

                              ((< number  #.(- (expt 2  7))) '#.(make-integer-type-descriptors int16 int32))
                              ((> number #.(1- (expt 2  8))) '#.(make-integer-type-descriptors int16 uint16 int32 uint32))
                              ((> number #.(1- (expt 2  7))) '#.(make-integer-type-descriptors uint8 int16 uint16 int32 uint32))

                              ((< number 0)                  '#.(make-integer-type-descriptors int8 int16 int32))
                              (t                             '#.(make-integer-type-descriptors int8 uint8 int16 uint16 int32 uint32)))))
             (integer (locally (declare (integer number))
                        (cond ((< number  #.(- (expt 2 31))) '())
                              ((> number #.(1- (expt 2 32))) '())
                              ((> number #.(1- (expt 2 31))) '#.(make-integer-type-descriptors uint32))

                              ((< number  #.(- (expt 2 15))) '#.(make-integer-type-descriptors int32))
                              ((> number #.(1- (expt 2 16))) '#.(make-integer-type-descriptors int32 uint32))
                              ((> number #.(1- (expt 2 15))) '#.(make-integer-type-descriptors uint16 int32 uint32))

                              ((< number  #.(- (expt 2  7))) '#.(make-integer-type-descriptors int16 int32))
                              ((> number #.(1- (expt 2  8))) '#.(make-integer-type-descriptors int16 uint16 int32 uint32))
                              ((> number #.(1- (expt 2  7))) '#.(make-integer-type-descriptors uint8 int16 uint16 int32 uint32))

                              ((< number 0)                  '#.(make-integer-type-descriptors int8 int16 int32))
                              (t                             '#.(make-integer-type-descriptors int8 uint8 int16 uint16 int32 uint32))))))))
    (loop for number in numbers
          for types = (copy-list (integer-types number))
          then (nintersection types (integer-types number) :key #'car :test #'char=)
          unless types return #\f
          finally (return (caar (stable-sort types #'< :key #'cdr))))))

(defun format-sam-tag (out tag value)
  "Write a SAM file TAG, dispatching on actual type of the given value.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5."
  (declare (stream out) (symbol tag) #.*optimization*)
  (write-tab out)
  (writestr out (symbol-name tag))
  (etypecase value
    (character    (writestr out (sbs ":A:")) (writec out value))           ; printable character
    (int32        (format out (sbs ":i:~D") value))                        ; signed 32-bit integer
    (integer      (format out (sbs ":f:~D") value))                        ; single-precision floating point number
    (single-float (format out (sbs ":f:~E") value))                        ; single-precision floating point number
    (double-float (format out (sbs ":f:~E") (coerce value 'single-float))) ; single-precision floating point number
    (string       (writestr out (sbs ":Z:"))                               ; printable string
                  (writestr out value))
    (vector       (writestr out (sbs ":H:"))                               ; byte array in hex format
                  (loop for byte across (the vector value)
                        do (format out (sbs "~2,'0X") byte)))
    (list         (writestr out (sbs ":B:"))                               ; integer or numeric array
                  (let ((type (common-number-type value)))
                    (writec out type)
                    (if (char= type num)
                      (loop for number in value do
                            (if (integerp number) (format out (sbs ",~D") number)
                              (format out (sbs ",~E") (coerce number 'single-float))))
                      (loop for number in value do (format out (sbs ",~D") number)))))))

(defun format-sam-alignment (out aln)
  "Write a SAM file read alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and 1.5."
  (declare (stream out) (sam-alignment aln) #.*optimization*)
  (writestr out (sam-alignment-qname aln)) (write-tab out)
  (format out (sbs "~D") (sam-alignment-flag aln))  (write-tab out)
  (writestr out          (sam-alignment-rname aln)) (write-tab out)
  (format out (sbs "~D") (sam-alignment-pos aln))   (write-tab out)
  (format out (sbs "~D") (sam-alignment-mapq aln))  (write-tab out)
  (writestr out          (sam-alignment-cigar aln)) (write-tab out)
  (writestr out          (sam-alignment-rnext aln)) (write-tab out)
  (format out (sbs "~D") (sam-alignment-pnext aln)) (write-tab out)
  (format out (sbs "~D") (sam-alignment-tlen aln))  (write-tab out)
  (writestr out          (sam-alignment-seq aln))   (write-tab out)
  (writestr out          (sam-alignment-qual aln))
  (loop for (key value) on (sam-alignment-tags aln) by #'cddr
        do (format-sam-tag out key value))
  (write-newline out))

(defun format-sam (out sam)
  "Write a complete SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1."
  (declare (stream out) (sam sam) #.*optimization*)
  (format-sam-header out (sam-header sam))
  (loop for aln in (sam-alignments sam)
        do (format-sam-alignment out aln)))

(defglobal *stderr* #+lispworks (sys:make-stderr-stream) #+sbcl sb-sys:*stderr*)

(defun setup-standard-streams ()
  (setq *standard-output* *stderr*
        *error-output* *stderr*
        *trace-output* *stderr*))

(defglobal *samtools* nil
  "Location of the samtools binary.")

(defun get-samtools ()
  "Determine location of the samtools binary."
  (or *samtools* 
      (with-open-program (program "command -v samtools" :buffered nil)
        (let ((line (read-line (program-stream program) nil)))
          (if line
            (setq *samtools* line)
            (error "samtools not found. Please download it from http://samtools.sourceforge.net and make sure that its binary is present in your PATH."))))))

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
  ;
  ; In LispWorks, :element-type 'base-char and :external-format :latin-1 is the default.
  ; This is both fast and fully compatible with all ASCII I/O to and from files and pipes.
  ;
  ; In SBCL, we use :element-type 'octet for input files and wrap them with our own ASCII conversion.
  ; For fast ASCII output, we use :element-type 'base-char and :external-format :utf-8.
  ; See open-program for details on pipes.
  ;
  (ecase direction
    (:input (cond ((eq kind :sam)
                   #+lispworks
                   (if (check-stdin pathname) *terminal-io*
                     (open pathname :direction :input :element-type 'base-char :if-does-not-exist :error))
                   #+sbcl
                   (make-ascii-stream (open pathname :direction :input :element-type 'octet :if-does-not-exist :error)))
                  (t (open pathname :direction :probe :if-does-not-exist :error)
                     (open-program (list (get-samtools) "view" (if header-only "-H" "-h") "-@" (write-to-string *number-of-threads*)
                                         (namestring (translate-logical-pathname pathname))) :direction :input))))
    (:output (cond ((eq kind :sam)
                    #+lispworks
                    (if (check-stdout pathname) *terminal-io*
                      (open pathname :direction :output :element-type 'base-char :if-exists :supersede))
                    #+sbcl
                    (open pathname :direction :output :element-type 'base-char :external-format :utf-8 :if-exists :supersede))

                   ((eq kind :cram)
                    (let ((reference-fasta *reference-fasta*)
                          (reference-fai *reference-fai*))
                      (cond ((and (null reference-fasta) (null reference-fai))
                             (error "When creating CRAM output, either a reference-fasta or a reference-fai must be provided."))
                            ((and reference-fasta reference-fai)
                             (error "When creating CRAM output, only either a reference-fasta or a reference-fai must be provided, but not both.")))
                      (open-program `(,(get-samtools) "view" "-C"
                                      "-@" ,(write-to-string *number-of-threads*)
                                      ,@(cond (reference-fasta
                                               `("-T" ,reference-fasta))
                                              (reference-fai
                                               `("-t" ,reference-fai)))
                                      "-o" ,(namestring (translate-logical-pathname pathname)) "-")
                                    :direction :output)))
                   (t (open-program (list (get-samtools) "view" "-Sb" "-@" (write-to-string *number-of-threads*)
                                          "-o" (namestring (translate-logical-pathname pathname)) "-") :direction :output))))))

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
  #+lispworks
  (let ((location (lw:pathname-location (pathname sibling)))
        (kind (sam-file-kind sibling)))
    (cond ((eq kind :sam)
           (let ((stream (hcl:open-temp-file :directory location)))
             (values stream (pathname stream))))
          (t (let ((pathname (hcl:create-temp-file :directory location)))
               (assert (delete-file pathname))
               (values (%open-sam pathname :output nil kind) pathname)))))
  #+sbcl
  (let* ((stream (cl-fad:open-temporary
                  :template (merge-pathnames "%" (merge-pathnames sibling (get-working-directory)))
                  :element-type 'octet))
         (pathname (pathname stream))
         (kind (sam-file-kind sibling)))
    (cond ((eq kind :sam)
           (values stream pathname))
          (t (close stream)
             (assert (delete-file pathname))
             (values (%open-sam pathname :output nil kind) pathname)))))

#+lispworks
(progn
  (declaim (inline close-sam sam-input sam-output))

  (defun close-sam (sam)
    "Close a SAM file, no matter whether it is an actual file or a pipe."
    (when (output-stream-p sam) (stream-flush-buffer sam))
    (unless (eq sam *terminal-io*) (close sam)))

  (defun sam-stream (sam)
    "Get the stream for a SAM file."
    sam))

#+sbcl
(progn
  (defun close-sam (sam)
    "Close a SAM file, no matter whether it is an actual file or a pipe."
    (cond ((streamp sam) (close sam))
          ((ascii-stream-p sam) (close-ascii-stream sam))
          ((program-p sam) (close-program sam))
          (t (error "Not a proper sam file connection: ~S." sam))))

  (defun sam-stream (sam)
    "Get the stream for a SAM file."
    (cond ((streamp sam) sam)
          ((ascii-stream-p sam) sam)
          ((program-p sam) (program-stream sam))
          (t (error "Not a proper sam file connection: ~S." sam)))))

(defun invoke-with-open-sam (function filename &rest args &key (direction :input) header-only)
  "Call a function and pass it a stream for reading or writing a SAM file."
  (declare (dynamic-extent args) (ignorable direction header-only))
  (let ((sam (apply #'open-sam filename args)))
    (unwind-protect
        (funcall function (sam-stream sam))
      (close-sam sam))))

(defmacro with-open-sam ((stream filename &rest args &key (direction :input) header-only) &body body)
  "Macro version of invoke-with-open-sam."
  (declare (ignore direction header-only))
  `(invoke-with-open-sam (lambda (,stream) ,@body) ,filename ,@args))
