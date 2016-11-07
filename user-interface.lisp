(in-package :elprep)
(in-simple-base-string-syntax)

(defun run-best-practices-pipeline-intermediate-list (file-in file-out &key (sorting-order :keep) (filters '()) (filters2 '()) (gc-on 0) (timed nil) (split-file nil))
  "Run the best practices pipeline. Version that uses an intermediate list so that sorting and mark-duplicates are supported."
  #+lispworks-32bit
  (declare (ignore gc-on))
  (let* ((filtered-reads (make-sam)) ; where the result of run-pipeline will be stored
         (working-directory (get-working-directory))
         (filename (merge-pathnames file-in working-directory)))
    #+lispworks-64bit
    (unless (= gc-on 2) (system:set-blocking-gen-num 0 :do-gc nil))
    #+sbcl
    (unless (= gc-on 2) (setf (sb-ext:bytes-consed-between-gcs)
                              (sb-ext:dynamic-space-size)))
    (flet ((first-phase ()
             (run-pipeline
              filename filtered-reads
              :filters filters
              :sorting-order sorting-order
              :split-file split-file)))
      (cond (timed
             (format t "Reading SAM into memory and applying filters.~%")
             (time (first-phase)))
            (t (first-phase))))
    #+lispworks-64bit
    (when (= gc-on 1) (system:marking-gc 0 :max-size 0))
    #+sbcl
    (when (= gc-on 1) (sb-ext:gc :full t))
    ; write to file
    (let ((file-out-name (merge-pathnames file-out working-directory)))
      (unless (check-stdout file-out)
        (ensure-directories-exist file-out-name))
      (flet ((second-phase ()
               (run-pipeline
                filtered-reads file-out-name
                :sorting-order (if (eq sorting-order :unsorted) :unsorted :keep)
                :filters filters2)))
        (cond (timed
               (format t "Write to file.~%")
               (time (second-phase)))
              (t (second-phase)))))))

(defun run-best-practices-pipeline (file-in file-out &key (sorting-order :keep) (filters '()) (gc-on 0) (timed nil) (split-file nil))
  "Run the best practices pipeline. Version that doesn't use an intermediate list when neither sorting nor mark-duplicates are needed."
  #+lispworks-32bit
  (declare (ignore gc-on))
  #+lispworks-64bit
  (unless (= gc-on 2) (system:set-blocking-gen-num 0 :do-gc nil))
  #+sbcl
  (unless (= gc-on 2) (setf (sb-ext:bytes-consed-between-gcs)
                            (sb-ext:dynamic-space-size)))
  (flet ((all-phases ()
           (let* ((working-directory (get-working-directory))
                  (file-in-name (merge-pathnames file-in working-directory))
                  (file-out-name (merge-pathnames file-out working-directory)))
             (run-pipeline file-in-name file-out-name
                           :filters filters
                           :sorting-order sorting-order
                           :split-file split-file))))
    (cond (timed
           (format t "Running pipeline.~%")
           (time (all-phases)))
          (t (all-phases)))))

(defvar *program-name* "elPrep"
  "Name of the elprep binary.")

(defvar *program-version* "2.5"
  "Version of the elprep binary.")

(defvar *program-url* "http://github.com/exascience/elprep"
  "URL for more information about elprep.")

(defvar *program-help* "sam-file sam-output-file ~% [--replace-reference-sequences sam-file] ~% [--filter-unmapped-reads [strict]] ~% [--replace-read-group read-group-string]~% [--mark-duplicates [remove] [deterministic]] ~% [--sorting-order [keep | unknown | unsorted | queryname | coordinate]] ~% [--clean-sam] ~% [--nr-of-threads nr] ~% [--gc-on [0 | 1 | 2]] ~% [--timed] ~% [--reference-t fai-file] ~% [--reference-T fasta-file] ~% [--split-file] ~%"
  "Help string for the elprep-script binary.")

;;; error handling

(defun create-log-filename ()
  "Create a log filename for writing error messages from within the elPrep binary."
  (multiple-value-bind
      (second minute hour date month year day daylight-p timezone)
      (get-decoded-time)
    (declare (ignore day daylight-p))
    (format nil "logs/elprep/elprep-~D-~2,'0D-~2,'0D-~2,'0D-~2,'0D-~2,'0D-GMT~D.log"
            year month date hour minute second timezone)))

#+lispworks
(defun elprep-debugger-hook (condition hook)
  "Write an error report to a log file and exist the elPrep binary."
  (declare (ignore hook))
  (dbg:log-bug-form (format nil "An error occurred in ~A ~A: ~A" *program-name* *program-version* condition)
                    :log-file (merge-pathnames (create-log-filename) (user-homedir-pathname))
                    :message-stream t)
  (lw:quit :status 1 :ignore-errors-p t))

#+sbcl
(defun elprep-debugger-hook (condition hook)
  "Write an error report to a log file and exist the elPrep binary."
  (declare (ignore hook))
  (let ((log-path (merge-pathnames (create-log-filename) (user-homedir-pathname))))
    (ensure-directories-exist log-path)
    (with-open-file (out log-path :direction :output :if-exists :rename)
      (format out "An error occurred in ~A ~A: ~A~%~%" *program-name* *program-version* condition)
      (format out "Lisp implementation: ~A ~A~%" (lisp-implementation-type) (lisp-implementation-version))
      (format out "Site: ~A~%" (long-site-name))
      (format out "Machine: ~A ~A ~A~%" (machine-instance) (machine-type) (machine-version))
      (format out "Software: ~A ~A~%~%" (software-type) (software-version))
      (let ((*standard-output* out)) (room t))
      (terpri out)
      (sb-debug:print-backtrace :stream out :print-thread t :print-frame-source t :method-frame-style :normal))
    (format t "Wrote error log to ~A~%" log-path)
    (sb-ext:exit :abort t)))

(defun exit-script (help-string &optional error-string &rest error-args)
  (when error-string (apply #'format t error-string error-args))
  (format t help-string)
  #+lispworks
  (lw:quit :status (if error-string 1 0) :ignore-errors-p t)
  #+sbcl
  (sb-ext:exit :abort error-string))

(defun elprep-filter-script ()
  "Command line script for elprep filter script."
  (let ((cmd-line (rest (command-line-arguments)))
        (sorting-order :keep)
        (nr-of-threads 1)
        (mark-duplicates-p nil)
        (gc-on nil)
        (timed nil)
        ; filters
        (replace-ref-seq-dct-filter nil)
        (ref-seq-dct nil)
        (remove-unmapped-reads-filter nil)
        (filter-unmapped-arg nil)
        (replace-read-group-filter nil)
        (read-group-string nil)
        (mark-duplicates-filter nil)
        (remove-duplicates-filter nil)
        (mark-duplicates-deterministic nil)
        (clean-sam-filter nil)
        (rename-chromosomes-filter nil)
        (reference-fai nil)
        (reference-fasta nil)
        (split-file nil))
    (loop with entry while cmd-line do (setq entry (pop cmd-line))
          if (string= entry "-h") do
          (exit-script *program-help*)
          else if (string= entry "--replace-reference-sequences")
          do (setf replace-ref-seq-dct-filter (list (replace-reference-sequence-dictionary-from-sam-file (setf ref-seq-dct (pop cmd-line)))))
          else if (string= entry "--filter-unmapped-reads")
          do (setf remove-unmapped-reads-filter (if (string= "strict" (first cmd-line))
                                                  (progn (setf filter-unmapped-arg (pop cmd-line)) (list #'filter-unmapped-reads-strict))
                                                  (list #'filter-unmapped-reads)))
          else if (string= entry "--replace-read-group")
          do (setf replace-read-group-filter (list (add-or-replace-read-group (parse-read-group-from-string (setf read-group-string (pop cmd-line))))))
          else if (string= entry "--mark-duplicates")
          do (setf mark-duplicates-filter (progn
                                            (setf mark-duplicates-p t)
                                            (loop (cond ((string= "remove" (first cmd-line))
                                                         (pop cmd-line)
                                                         (setf remove-duplicates-filter (list #'filter-duplicate-reads)))
                                                        ((string= "deterministic" (first cmd-line))
                                                         (pop cmd-line)
                                                         (setf mark-duplicates-deterministic t))
                                                        (t (return (list (mark-duplicates mark-duplicates-deterministic))))))))
          else if (string= entry "--sorting-order")
          do (let ((so (first cmd-line))) ; peek it
               (if (or (not so) (search "--" so)) ; if no sorting order or next option, use default sorting order = keep
                 (setq sorting-order :keep)
                 (progn
                   (pop cmd-line)
                   (setq sorting-order (intern (string-upcase so) :keyword))
                   (unless (member sorting-order '(:keep :unknown :unsorted :queryname :coordinate))
                     (exit-script *program-help* "Invalid sorting-order ~A.~%" so)))))
          else if (string= entry "--clean-sam")
          do (setf clean-sam-filter (list #'clean-sam))
          else if (string= entry "--nr-of-threads")
          do (setf nr-of-threads (parse-integer (pop cmd-line)))
          else if (string= entry "--gc-on")
          do (let ((lvl (first cmd-line)))
               (if (or (not lvl) (search "--" lvl)) ; no sorting level given
                 (setf gc-on 0)
                 (progn (setf lvl (parse-integer (pop cmd-line) :junk-allowed t))
                   (if (and lvl (or (= lvl 0) (= lvl 1) (= lvl 2)))
                     (setf gc-on lvl)
                     (exit-script *program-help* "Invalid gc-on option ~A.~%" lvl)))))
          else if (string= entry "--timed")
          do (setf timed t)
          else if (string= entry "--reference-t")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *program-help* "Please provide reference file with --reference-t.~%")
                 (setf *reference-fai* (setf reference-fai (pop cmd-line)))))
          else if (string= entry "--reference-T")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *program-help* "Please provide reference file with --reference-T.~%")
                 (setf *reference-fasta* (setf reference-fasta (pop cmd-line)))))
          else if (string= entry "--rename-chromosomes")
          do (setf rename-chromosomes-filter (list #'rename-chromosomes))
          else if (string= entry "--split-file")
          do (setf split-file t)
          else collect entry into conversion-parameters ; main required parameters, input and output sam files
          finally
          (when (/= (length conversion-parameters) 2)
            (exit-script *program-help* "Incorrect number of parameters: Expected 2, got ~S~@[ ~A~].~%"
                         (length conversion-parameters) conversion-parameters))
          ; print a warning when replacing ref seq dictionary and trying to keep the order
          (when (and replace-ref-seq-dct-filter (eq sorting-order :keep))
            (format t "WARNING: Requesting to keep the order of the input file while replacing the reference sequence dictionary may force an additional sorting phase to ensure the original sorting order is respected."))
          ; check that when cram is used, the reference dictionary is passed
          (when (and (eq (sam-file-kind (second conversion-parameters)) :cram)
                     (not reference-fai) (not reference-fasta))
            (exit-script *program-help* "ERROR: Attempting to output to cram without specifying a reference file. Please add --reference-t or --reference-T to your call.~%"))
          ; set the default gc-settings, unless the user set them
          (unless gc-on
            (if (or mark-duplicates-filter 
                    (and replace-ref-seq-dct-filter (eq sorting-order :keep)) 
                    (eq sorting-order :coordinate) (eq sorting-order :queryname))
              (setf gc-on 0)
              (setf gc-on 2)))
          (let* ((cmd-string
                  (with-output-to-string (s nil :element-type 'base-char)
                    (format s "~a ~a ~a"
                            #+lispworks (first (command-line-arguments))
                            #+sbcl (first sb-ext:*posix-argv*)
                            (first conversion-parameters)
                            (second conversion-parameters))
                    (when remove-unmapped-reads-filter (format s " --filter-unmapped-reads~@[ ~a~]" filter-unmapped-arg))
                    (when clean-sam-filter (format s " --clean-sam"))
                    (when replace-ref-seq-dct-filter (format s " --replace-reference-sequences ~a" ref-seq-dct))
                    (when replace-read-group-filter (format s " --replace-read-group ~s" read-group-string))
                    (when mark-duplicates-filter
                      (format s " --mark-duplicates")
                      (when remove-duplicates-filter
                        (format s " remove"))
                      (when mark-duplicates-deterministic
                        (format s " deterministic")))
                    (format s " --sorting-order ~(~a~) --gc-on ~a --nr-of-threads ~a" sorting-order gc-on nr-of-threads)
                    (when timed (format s " --timed"))
                    (when reference-fai (format s " --reference-t ~a" reference-fai))
                    (when reference-fasta (format s " --reference-T ~a" reference-fasta))
                    (when split-file (format s " --split-file"))))
                 ; optimal order for filters
                 (filters (nconc (list (add-pg-line (format nil "~A ~A" *program-name* *program-version*)
                                                    :pn *program-name*
                                                    :vn *program-version*
                                                    :ds *program-url*
                                                    :cl cmd-string))
                                 remove-unmapped-reads-filter
                                 rename-chromosomes-filter
                                 clean-sam-filter
                                 replace-ref-seq-dct-filter       
                                 replace-read-group-filter
                                 (when (or replace-ref-seq-dct-filter mark-duplicates-filter (member sorting-order '(:coordinate :queryname))) (list #'add-refid))
                                 mark-duplicates-filter
                                 (list #'filter-optional-reads)))
                 (filters2 remove-duplicates-filter))
            (format t "Executing command:~%  ~a~%" cmd-string)
            (setq *number-of-threads* nr-of-threads
                  *reference-fasta* reference-fasta
                  *reference-fai* reference-fai)
            (if (or mark-duplicates-p
                    (member sorting-order '(:coordinate :queryname))
                    (and replace-ref-seq-dct-filter (eq sorting-order :keep)))
              (let ((conversion-parameters
                     (nconc conversion-parameters
                            `(:sorting-order ,sorting-order
                              :filters ,filters :filters2 ,filters2 :gc-on ,gc-on :timed ,timed :split-file ,split-file))))
                (apply #'run-best-practices-pipeline-intermediate-list conversion-parameters))
              (let ((conversion-parameters
                     (nconc conversion-parameters
                            `(:sorting-order ,sorting-order
                              :filters ,filters :gc-on ,gc-on :timed ,timed :split-file ,split-file))))
                (apply #'run-best-practices-pipeline conversion-parameters)))))))

(defvar *split-program-help* "split [sam-file | /path/to/input/] /path/to/output/ ~% [--output-prefix name] ~% [--output-type [sam | bam | cram]] ~% [--nr-of-threads nr] ~% [--reference-t fai-file] ~% [--reference-T fasta-file] ~%"
  "Help string for the elprep-split-script binary.")

(defun elprep-split-script ()
  "Command line script for elprep split script."
  (let ((cmd-line (rest (rest (command-line-arguments)))) ; skip elprep split part of the command
        (input nil) (output-path nil) (output-prefix nil) (output-type :sam) (output-extension nil) (nr-of-threads 1) (reference-fai nil) (reference-fasta nil))
    (loop with entry while cmd-line do (setq entry (pop cmd-line))
          if (string= entry "-h") do
          (exit-script *split-program-help*)
          else if (string= entry "--output-type")
          do (let ((output-kind (intern (string-upcase (first cmd-line)) :keyword)))
               (cond ((or (not output-kind) (search "--" cmd-line) (not (member output-kind '(:sam :bam :cram)))) ; no correct output type given
                      (exit-script *split-program-help* "Please provide a valid output type.~%"))
                     (t (pop cmd-line)
                        (setf output-type output-kind)
                        (setf output-extension (ecase output-type (:bam "bam") (:sam "sam") (:cram "cram"))))))
          else if (string= entry "--output-prefix")
          do (let ((prefix (first cmd-line)))
               (cond ((or (not prefix) (search "--" cmd-line)) ; no prefix given
                      (exit-script *split-program-help* "Please provide a valid output prefix.~%"))
                     (t 
                      (pop cmd-line) 
                      (setf output-prefix prefix))))
          else if (string= entry "--nr-of-threads")
          do (setf nr-of-threads (parse-integer (pop cmd-line)))
          else if (string= entry "--reference-t")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *split-program-help* "Please provide reference file with --reference-t.~%")
                 (setf *reference-fai* (setf reference-fai (pop cmd-line)))))
          else if (string= entry "--reference-T")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *split-program-help* "Please provide reference file with --reference-T.~%")
                 (setf *reference-fasta* (setf reference-fasta (pop cmd-line)))))
          else collect entry into io-parameters
          finally
          (when (/= (length io-parameters) 2) ; checks on input parameters
            (exit-script *split-program-help* "Incorrect number of parameters: Expected 2, got ~S~@[ ~A~].~%"
                         (length io-parameters) io-parameters))
          ; fill in defaults
          (setf input (first io-parameters))
          (setf output-path (second io-parameters))
          ; check that output path is really a path
          (when (or (pathname-name (pathname output-path)) (not (pathname-directory (pathname output-path))))
            (exit-script *split-program-help* "Given output path is not a path: ~a ~%" output-path))
          (unless output-prefix (setf output-prefix (pathname-name input)))
          (unless output-extension (setf output-extension (ecase (sam-file-kind input) (:bam "bam") (:sam "sam") (:cram "cram"))))
          ; print feedback
          (let ((cmd-string
                 (with-output-to-string (s nil :element-type 'base-char)
                   (format s "~a split ~a ~a " (first (command-line-arguments)) (first io-parameters) (second io-parameters))
                   (format s "--output-prefix ~a " output-prefix)
                   (format s "--output-type ~(~a~) " output-type)
                   (when reference-fai
                     (format s "--reference-t ~a" reference-fai))
                   (when reference-fasta
                     (format s "--reference-T ~a" reference-fasta)))))
            (format t "Executing command:~%  ~a~%" cmd-string))
          (setq *number-of-threads* nr-of-threads)
          (ensure-directories-exist output-path)
          (let ((working-directory (get-working-directory)))
            (split-file-per-chromosome (merge-pathnames input working-directory) (merge-pathnames output-path working-directory) output-prefix output-extension)))))

(defvar *merge-program-help* "merge /path/to/input/ sam-output-file ~% [--nr-of-threads nr] ~% [--reference-t fai-file] ~% [--reference-T fasta-file] ~%"
  "Help string for the elprep-merge-script binary.")

(defun elprep-merge-script ()
  "Command line script for elprep merge script."
  (let ((cmd-line (rest (rest (command-line-arguments)))) ; skip elprep merge part of the command
        (input-path nil) (output nil) (nr-of-threads 1) (reference-fai nil) (reference-fasta nil))
    (loop with entry while cmd-line do (setq entry (pop cmd-line))
          if (string= entry "-h") do
          (exit-script *merge-program-help*)
          else if (string= entry "--nr-of-threads")
          do (setf nr-of-threads (parse-integer (pop cmd-line)))
          else if (string= entry "--reference-t")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *merge-program-help* "Please provide reference file with --reference-t.~%")
                 (setf *reference-fai* (setf reference-fai (pop cmd-line)))))
          else if (string= entry "--reference-T")
          do (let ((ref (first cmd-line)))
               (if (or (not ref) (search "--" cmd-line)) ; no file given
                 (exit-script *merge-program-help* "Please provide reference file with --reference-T.~%")
                 (setf *reference-fasta* (setf reference-fasta (pop cmd-line)))))
          else collect entry into io-parameters
          finally         
          (when (/= (length io-parameters) 2)
            (exit-script *merge-program-help* "Incorrect number of parameters: Expected 2, got ~S~@[ ~A~].~%"
                         (length io-parameters) io-parameters))
          (setf input-path (first io-parameters))
          (setf output (second io-parameters))
          (let ((files-to-merge (cl-fad:list-directory input-path)))
            (cond ((not files-to-merge)
                   (exit-script *merge-program-help* "Given directory ~a does not exist. ~%" input-path))
                  ((= (length files-to-merge) 1)
                   (exit-script *merge-program-help* "Given directory ~a does not contain any files to merge. ~%" input-path)))
            ; extract the input prefix
            (let* ((header (with-open-sam (in (first files-to-merge) :direction :input) (parse-sam-header in)))
                   (input-prefix
                    (loop for file in files-to-merge
                          do (let* ((ffile (file-namestring file))
                                    (idx (search "-unmapped" ffile :from-end t))) ; there is at least a file with unmapped tag
                               (when idx (return (subseq ffile 0 idx))))))
                   (first-file-name (file-namestring (first files-to-merge)))                   
                   (input-extension (ecase (sam-file-kind first-file-name) (:bam "bam") (:sam "sam") (:cram "cram"))))
              ; print feedback
              (let ((cmd-string
                     (with-output-to-string (s nil :element-type 'base-char)
                       (format s "~a merge ~a ~a" (first (command-line-arguments)) input-path output)
                       (when reference-fai
                         (format s "--reference-t ~a" reference-fai))
                       (when reference-fasta
                         (format s "--reference-T ~a" reference-fasta)))))
                (format t "Executing command:~%  ~a~%" cmd-string))
              (setq *number-of-threads* nr-of-threads)
              (let ((working-directory (get-working-directory)))
                (let ((sorting-order (getf (sam-header-hd header) :so)))
                  (if (string= sorting-order "coordinate")
                      (merge-sorted-files-split-per-chromosome (merge-pathnames input-path working-directory) 
                                                               (merge-pathnames output working-directory) input-prefix input-extension header)
                    (merge-unsorted-files-split-per-chromosome (merge-pathnames input-path working-directory) 
                                                               (merge-pathnames output working-directory) input-prefix input-extension header)))))))))

(defun elprep-script ()
  "Command line script for elPrep."
  (setf *debugger-hook* #'elprep-debugger-hook)
  (setup-standard-streams)
  (format t "~A version ~A. See ~A for more information.~%"
          *program-name* *program-version* *program-url*)
  (let ((cmd-line (rest (command-line-arguments))))
    (when (null cmd-line)
      (format t "Incorrect number of parameters. ~%")
      (format t "Filter parameters: ~%")
      (format t *program-help*)
      (format t "Split parameters: ~%")
      (format t *split-program-help*)
      (format t "Merge parameters: ~%")
      (format t *merge-program-help*)
      (exit-script ""))
    (cond ((string= (first cmd-line) "split") (elprep-split-script))
          ((string= (first cmd-line) "merge") (elprep-merge-script))
          (t 
           (elprep-filter-script)))))
