(in-package :elprep)

(deftype  octet () '(unsigned-byte 8))
(deftype  int16 () '(signed-byte 16))
(deftype uint16 () '(unsigned-byte 16))
(deftype  int32 () '(signed-byte 32))

(defvar *sam-file-format-version* (sbs "1.4")
  "The SAM file format version string supported by this library.
   This is entered by default in a @HD line in the header section of a SAM file, unless user code explicitly asks for a different version number.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
   Default is \"1.4\".")

(defstruct sam-header
  "The information stored in the header section of a SAM file. 
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
   The struct sam-header has a default constructor make-sam-header.
   Accessor sam-header-hd of type property list refers to the @HD line.
   Accessor sam-header-sq of type list of property lists refers to the @SQ lines.
   Accessor sam-header-rg of type list of property lists refers to the @RG lines.
   Accessor sam-header-pg of type list of property lists refers to the @PG lines.
   Accessor sam-header-co of type list of strings refers to the @CO lines.
   Accessor sam-header-user-tags of type property list refers to lines defined by end users.
   All information in property lists are stored in the same order as they occur in a SAM file."
  (hd '() :type list)
  (sq '() :type list)
  (rg '() :type list)
  (pg '() :type list)
  (co '() :type list)
  (user-tags '() :type list))

(setf (documentation 'make-sam-header 'function)
      "Default constructor for struct sam-header."
      (documentation 'sam-header-p 'function)
      "Default predicate for struct sam-header."
      (documentation 'copy-sam-header 'function)
      "Default copier function for struct sam-header."
      (documentation 'sam-header-hd 'function)
      "Access the sam-header @HD line of type property list."
      (documentation 'sam-header-sq 'function)
      "Access the sam-header @SQ lines of type list of property lists."
      (documentation 'sam-header-rg 'function)
      "Access the sam-header @RG lines of type list of property lists."
      (documentation 'sam-header-pg 'function)
      "Access the sam-header @PG lines of type list of property lists."
      (documentation 'sam-header-co 'function)
      "Access the sam-header @CO lines of type list of strings."
      (documentation 'sam-header-user-tags 'function)
      "Access the sam-header user-defined header lines of type property list.")

(defun sam-header-ensure-hd (hdr &key (vn *sam-file-format-version*) (so :unknown so-p))
  "Ensure an @HD line is present in the given sam-header.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag @HD."
  (declare (sam-header hdr) #.*optimization*)
  (setf (getf (sam-header-hd hdr) :VN) vn)
  (when so-p
    (setf (getf (sam-header-hd hdr) :SO)
          (string-downcase (string so))))
  hdr)

(declaim (inline sam-header-user-tag (setf sam-header-user-tag)))

(defun sam-header-user-tag (hdr tag)
  "Access a sam-header user tag."
  (declare (sam-header hdr) (symbol tag) #.*optimization*)
  (getf (sam-header-user-tags hdr) tag))

(defun (setf sam-header-user-tag) (value hdr tag)
  "Access a sam-header user tag."
  (declare (sam-header hdr) (symbol tag) #.*optimization*)
  (setf (getf (sam-header-user-tags hdr) tag) value))

(declaim (inline sam-header-user-tag-p))

(defun sam-header-user-tag-p (code)
  "Does this tag string represent a user-defined tag?"
  (declare (simple-string code) #.*optimization*)
  (loop for c of-type character across code
        thereis (and (char<= #\a c) (char<= c #\z))))

(defstruct sam-alignment
  "A single read alignment with mandatory and optional fields that can be contained in a SAM file alignment line.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and 1.5.
   The struct sam-alignment has a default constructor make-sam-alignment.
   Accessor sam-alignment-qname of type simple-base-string refers to the Query template NAME.
   Accessor sam-alignment-flag of type fixnum refers to the bitwise FLAG.
   Accessor sam-alignment-rname of type simple-base-string refers to the Reference sequence NAME.
   Accessor sam-alignment-pos of type int32 refers to the 1-based leftmost mapping POSition.
   Accessor sam-alignment-mapq of type fixnum refers to the MAPping quality.
   Accessor sam-alignment-cigar of type simple-base-string refers to the CIGAR string.
   Accessor sam-alignment-rnext of type simple-base-string refers to the Reference sequence name of the mate/NEXT read.
   Accessor sam-alignment-pnext of type int32 refers to the 1-based leftmost mapping Position of the mate/NEXT read.
   Accessor sam-alignment-tlen of type int32 refers to the observed Template LENgth.
   Accessor sam-alignment-seq of type simple-base-string refers to the segment SEQuence.
   Accessor sam-alignment-qual of type simple-base-string refers to the ASCII of Phred-scaled base QUALity+33.
   Accessor sam-alignment-tags of type property list refers to the optional fields in a read alignment.
   Accessor sam-alignment-xtags of type property list refers to additional optional fields not stored in SAM files, but reserved for other storage formats.
   Accessor sam-alignment-temps of type property list refers to additional optional fields not stored in any storage format, but reserved for temporary values in filters."
  (qname   empty-sbs :type simple-base-string)
  (flag    0         :type uint16)
  (rname   empty-sbs :type simple-base-string)
  (pos     0         :type int32)
  (mapq    0         :type octet)
  (cigar   empty-sbs :type simple-base-string)
  (rnext   empty-sbs :type simple-base-string)
  (pnext   0         :type int32)
  (tlen    0         :type int32)
  (seq     empty-sbs :type simple-base-string)
  (qual    empty-sbs :type simple-base-string)
  (tags   '()        :type list)
  (xtags  '()        :type list)
  (temps  '()        :type list))

(setf (documentation 'make-sam-alignment 'function)
      "Default constructor for struct sam-alignment."
      (documentation 'sam-alignment-p 'function)
      "Default predicate for struct sam-alignment."
      (documentation 'copy-sam-alignment 'function)
      "Default copier function for struct sam-header."
      (documentation 'sam-alignment-qname 'function)
      "Access the sam-alignment Query template NAME of type simple-base-string."
      (documentation 'sam-alignment-flag 'function)
      "Access the sam-alignment bitwise FLAG of type fixnum."
      (documentation 'sam-alignment-rname 'function)
      "Access the sam-alignment Reference sequence NAME of type simple-base-string."
      (documentation 'sam-alignment-pos 'function)
      "Access the sam-alignment 1-based leftmost mapping POSition of type int32."
      (documentation 'sam-alignment-mapq 'function)
      "Access the sam-alignment MAPping quality of type fixnum."
      (documentation 'sam-alignment-cigar 'function)
      "Access the sam-alignment CIGAR string of type simple-base-string."
      (documentation 'sam-alignment-rnext 'function)
      "Access the sam-alignment Reference sequence name of the mate/NEXT read of type simple-base-string."
      (documentation 'sam-alignment-pnext 'function)
      "Access the sam-alignment 1-based leftmost mapping Position of the mate/NEXT read of type int32."
      (documentation 'sam-alignment-tlen 'function)
      "Access the sam-alignment observed Template LENgth of type int32."
      (documentation 'sam-alignment-seq 'function)
      "Access the sam-alignment segment SEQuence of type simple-base-string."
      (documentation 'sam-alignment-qual 'function)
      "Access the sam-alignment ASCII of Phred-scaled base QUALity+33 of type simple-base-string."
      (documentation 'sam-alignment-tags 'function)
      "Access the sam-alignment optional fields of type property list."
      (documentation 'sam-alignment-xtags 'function)
      "Access the sam-alignment extended optional fields of type property list."
      (documentation 'sam-alignment-temps 'function)
      "Access the sam-alignment temporary values of type property list.")

(declaim (inline sam-alignment-tag (setf sam-alignment-tag)))

(defun sam-alignment-tag (aln tag)
  "Access a sam-alignment optional field in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (getf (sam-alignment-tags aln) tag))

(defun (setf sam-alignment-tag) (value aln tag)
  "Access a sam-alignment optional field in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (setf (getf (sam-alignment-tags aln) tag) value))

(declaim (inline sam-alignment-xtag (setf sam-alignment-xtag)))

(defun sam-alignment-xtag (aln tag)
  "Access a sam-alignment extended optional field in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (getf (sam-alignment-xtags aln) tag))

(defun (setf sam-alignment-xtag) (value aln tag)
  "Access a sam-alignment extended optional field in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (setf (getf (sam-alignment-xtags aln) tag) value))

(declaim (inline sam-alignment-temp (setf sam-alignment-temp)))

(defun sam-alignment-temp (aln tag)
  "Access a sam-alignment temporary value in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (getf (sam-alignment-temps aln) tag))

(defun (setf sam-alignment-temp) (value aln tag)
  "Access a sam-alignment temporary value in the sam-alignment."
  (declare (sam-alignment aln) (symbol tag) #.*optimization*)
  (setf (getf (sam-alignment-temps aln) tag) value))

(declaim (inline sam-alignment-refid (setf sam-alignment-refid)))

(defun sam-alignment-refid (aln)
  "Access the commonly used REFID temporary field in the sam-alignment."
  (sam-alignment-temp aln :refid))

(defun (setf sam-alignment-refid) (new-value aln)
  "Access the commonly used REFID temporary field in the sam-alignment."
  (setf (sam-alignment-temp aln :refid) new-value))

(declaim (inline check-refid-type))

(defun check-refid-type (value)
  "Ensure that value is probably an int32."
  (check-type value
              #+x86 integer
              #+x86-64 fixnum
              "a 32-bit integer")
  value)

(defun coordinate< (aln1 aln2)
  "Compare two alignments according to their coordinate.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag @HD, SO."
  (declare (sam-alignment aln1 aln2) #.*optimization*)
  (let ((refid1 (check-refid-type (sam-alignment-refid aln1)))
        (refid2 (check-refid-type (sam-alignment-refid aln2))))
    (declare (int32 refid1 refid2))
    (or (< refid1 refid2)
        (and (= refid1 refid2)
             (< (sam-alignment-pos aln1)
                (sam-alignment-pos aln2))))))

(defconstant +multiple+        #x1
  "Bit value for sam-alignment-flag: template having multiple segments in sequencing.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +proper+          #x2
  "Bit value for sam-alignment-flag: each segment properly aligned according to the aligner.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +unmapped+        #x4
  "Bit value for sam-alignment-flag: segment unmapped.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +next-unmapped+   #x8
  "Bit value for sam-alignment-flag: next segment in the template unmapped.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +reversed+       #x10
  "Bit value for sam-alignment-flag: SEQ being reversed complemented.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +next-reversed+  #x20
  "Bit value for sam-alignment-flag: SEQ of the next segment in the template being reversed.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +first+          #x40
  "Bit value for sam-alignment-flag: the first segment in the template.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +last+           #x80
  "Bit value for sam-alignment-flag: the last segment in the template.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +secondary+     #x100
  "Bit value for sam-alignment-flag: secondary alignment.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +qc-failed+     #x200
  "Bit value for sam-alignment-flag: not passing quality controls.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +duplicate+     #x400
  "Bit value for sam-alignment-flag: PCR or optical duplicate.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")
(defconstant +supplementary+ #x800
  "Bit value for sam-alignment-flag: supplementary alignment.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.")

(declaim (inline sam-alignment-multiple-p
                 sam-alignment-proper-p
                 sam-alignment-unmapped-p
                 sam-alignment-next-unmapped-p
                 sam-alignment-reversed-p
                 sam-alignment-next-reversed-p
                 sam-alignment-first-p
                 sam-alignment-last-p
                 sam-alignment-secondary-p
                 sam-alignment-qc-failed-p
                 sam-alignment-duplicate-p
                 sam-alignment-supplementary-p))

(defun sam-alignment-multiple-p (aln)
  "Check for template having multiple segments in sequencing in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +multiple+) 0))

(defun sam-alignment-proper-p (aln)
  "Check for each segment being properly aligned according to the aligner in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +proper+) 0))

(defun sam-alignment-unmapped-p (aln)
  "Check for segment unmapped in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +unmapped+) 0))

(defun sam-alignment-next-unmapped-p (aln)
  "Check for next segment in the template unmapped in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +next-unmapped+) 0))

(defun sam-alignment-reversed-p (aln)
  "Check for SEQ being reversed complemented in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +reversed+) 0))

(defun sam-alignment-next-reversed-p (aln)
  "Check for SEQ of the next segment in the template being reversed in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +next-reversed+) 0))

(defun sam-alignment-first-p (aln)
  "Check for being the first segment in the template in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +first+) 0))

(defun sam-alignment-last-p (aln)
  "Check for being the last segment in the template in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +last+) 0))

(defun sam-alignment-secondary-p (aln)
  "Check for secondary alignment in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +secondary+) 0))

(defun sam-alignment-qc-failed-p (aln)
  "Check for not passing quality controls in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +qc-failed+) 0))

(defun sam-alignment-duplicate-p (aln)
  "Check for PCR or optical duplicate in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +duplicate+) 0))

(defun sam-alignment-supplementary-p (aln)
  "Check for supplementary alignment in sam-alignment-flag.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2."
  (declare (sam-alignment aln) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) +supplementary+) 0))

(declaim (inline sam-alignment-flag-every
                 sam-alignment-flag-some
                 sam-alignment-flag-notevery
                 sam-alignment-flag-notany))

(defun sam-alignment-flag-every (aln flag)
  "Check for every bit in flag being set in sam-alignment-flag."
  (declare (sam-alignment aln) (fixnum flag) #.*optimization*)
  (= (logand (sam-alignment-flag aln) flag) flag))

(defun sam-alignment-flag-some (aln flag)
  "Check for some bits in flag being set in sam-alignment-flag."
  (declare (sam-alignment aln) (fixnum flag) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) flag) 0))

(defun sam-alignment-flag-notevery (aln flag)
  "Check for not every bit in flag being set in sam-alignment-flag."
  (declare (sam-alignment aln) (fixnum flag) #.*optimization*)
  (/= (logand (sam-alignment-flag aln) flag) flag))

(defun sam-alignment-flag-notany (aln flag)
  "Check for not any bit in flag being set in sam-alignment-flag."
  (declare (sam-alignment aln) (fixnum flag) #.*optimization*)
  (= (logand (sam-alignment-flag aln) flag) 0))

(declaim (inline make-sam))

(defstruct sam
  "A complete SAM data set that can be contained in a SAM file.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.
   The struct sam has a default constructor make-sam.
   Accessor sam-header of type sam-header refers to the header.
   Accessor sam-alignments of type list of sam-alignment refers to the read alignments."
  (header (make-sam-header) :type sam-header)
  (alignments '() :type list))

(setf (documentation 'make-sam 'function)
      "Default constructor for struct sam."
      (documentation 'sam-p 'function)
      "Default predicate for struct sam."
      (documentation 'copy-sam 'function)
      "Default copier function for struct sam."
      (documentation 'sam-header 'function)
      "Access the sam header of type sam-header."
      (documentation 'sam-alignments 'function)
      "Access the sam list of sam-alignment instances.")

;;; mapping cigar strings to alists or avectors

(define-symbol-macro cigar-operations (sbs "MmIiDdNnSsHhPpXx="))

(defconstant +min-cigar-operation+ (reduce #'min cigar-operations :key #'char-code)
  "The smallest CIGAR operation.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(defconstant +max-cigar-operation+ (reduce #'max cigar-operations :key #'char-code)
  "The largest CIGAR operation.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(eval-when (#+sbcl :compile-toplevel :load-toplevel :execute)
  (defun make-cigar-operations-table ()
    "Map CIGAR operations represented as characters to upcase characters.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
    (let ((table (make-array (1+ (- +max-cigar-operation+
                                    +min-cigar-operation+))
                             :initial-element nil
                             #+lispworks :allocation #+lispworks :long-lived
                             #+lispworks :single-thread #+lispworks t)))
      (loop for char across cigar-operations do
            (setf (svref table (- (char-code char) +min-cigar-operation+))
                  (char-upcase char)))
      table)))

(defglobal *cigar-operations*
  (make-cigar-operations-table)
  "Map CIGAR operations represented as characters to upcase characters.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(declaim (inline make-cigar-operation))

(defun make-cigar-operation (table cigar i)
  "Parse a CIGAR length + operation from position i in the cigar string.
   Return a cons cell with upcase character + length, and a next position in the cigar string.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (simple-base-string cigar) (fixnum i) #.*optimization*)
  (loop for j of-type fixnum = i then (1+ j)
        for char of-type fixnum = (char-code (schar cigar j))
        when (or (< char #.(char-code #\0)) (> char #.(char-code #\9)))
        return (values (let ((length (parse-integer cigar :start i :end j))
                             (operation (svref table (the fixnum (- char +min-cigar-operation+)))))
                         (if operation (cons operation length)
                           (error "Invalid CIGAR operation ~S." (code-char char))))
                       (the fixnum (1+ j)))))

(defglobal *cigar-list-cache*
  (let ((table (make-synchronized-hash-table :test #'equal)))
    (setf (gethash (sbs "*") table) '())
    table)
  "Cache CIGAR Strings to association lists.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(declaim (notinline slow-scan-cigar-string-to-list))

(defun slow-scan-cigar-string-to-list (cigar)
  "Convert a cigar string to an association list, slow path.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (simple-base-string cigar) #.*optimization*)
  (loop with table = *cigar-operations*
        with i of-type fixnum = 0
        with cigar-operation do
        (setf (values cigar-operation i) (make-cigar-operation table cigar i))
        collect cigar-operation into list
        until (= i (length cigar))
        finally (return (with-modify-hash (key value found) 
                            ((the hash-table *cigar-list-cache*) cigar)
                          (if found value list)))))

(declaim (inline fast-scan-cigar-string-to-list))

(defun fast-scan-cigar-string-to-list (cigar)
  "Convert a cigar string to an association list, fast path.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (simple-base-string cigar) #.*optimization*)
  (multiple-value-bind (value found)
      (gethash cigar (the hash-table *cigar-list-cache*))
    (if found value (slow-scan-cigar-string-to-list cigar))))

(defglobal *cigar-vector-cache*
  (let ((table (make-synchronized-hash-table :test #'equal)))
    (setf (gethash (sbs "*") table) #())
    table)
  "Cache CIGAR strings to assocation vectors.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.")

(declaim (notinline slow-scan-cigar-string-to-vector))

(defun slow-scan-cigar-string-to-vector (cigar)
  "Convert a cigar string to an assocation vector, slow path.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (simple-base-string cigar) #.*optimization*)
  (loop with table = *cigar-operations*
        with vector = (make-array (loop for i of-type fixnum = 0 then (1+ i)
                                        until (= i (length cigar))
                                        count (let ((char (char-code (schar cigar i))))
                                                (declare (fixnum char))
                                                (or (< char #.(char-code #\0)) (> char #.(char-code #\9))))))
        with i of-type fixnum = 0
        with cigar-operation
        for index of-type fixnum from 0 do
        (setf (values cigar-operation i) (make-cigar-operation table cigar i))
        (setf (svref vector index) cigar-operation)
        until (= i (length cigar))
        finally (return (with-modify-hash (key value found)
                            ((the hash-table *cigar-vector-cache*) cigar)
                          (if found value vector)))))

(declaim (inline fast-scan-cigar-string-to-vector))

(defun fast-scan-cigar-string-to-vector (cigar)
  "Convert a cigar string to an association vector, fast path.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (simple-base-string cigar) #.*optimization*)
  (or (gethash cigar (the hash-table *cigar-vector-cache*))
      (slow-scan-cigar-string-to-vector cigar)))

(declaim (inline scan-cigar-string))

(defun scan-cigar-string (type cigar)
  "Convert a cigar string to an association 'list or 'vector.
   See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6."
  (declare (symbol type))
  (ecase type
    (list   (fast-scan-cigar-string-to-list cigar))
    (vector (fast-scan-cigar-string-to-vector cigar))))

(define-compiler-macro scan-cigar-string (&whole form type cigar)
  (case type
    (list   `(fast-scan-cigar-string-to-list ,cigar))
    (vector `(fast-scan-cigar-string-to-vector ,cigar))
    (t form)))
