(in-package :elprep)
(in-simple-base-string-syntax)

(defgeneric run-pipeline (input output &key &allow-other-keys)
  (:documentation
   "Reads a SAM data set from input, applies filters to the header and the alignments,
    and writes the result to output. Eventually returns the output.

    Uses *number-of-threads* to execute the pipeline in parallel. (Default: 1.)

    Accepted keywords:
    - :filters: A list of filters. Default is the empty list.
    - :destructive: Operate destructively on the input or not. One of :default (the default), false, or true.
    - :sorting-order: One of :keep (the default), :unsorted, :unkown, :coordinate, or :queryname.
    - :chunk-size: Number of alignments to read from input at a time. Default is +default-chunk-size+.
    - :split-file: Operate on intermediate split files or not. One of false (the default), or true.

    A filter is a function that accepts zero or more arguments and returns zero or more values:
    - The first return value is nil or a next-level filter.
    - The second value is nil or a thunk to initialize filtering at the next level.
    - The third value is nil or a thunk to finalize filtering at the next level.

    The levels are: global, then thread-local, then local.

    For operating on SAM data sets:
    - A global filter receives the SAM header that it can modify.
      If an alignment filter is needed, the global filter needs to return a thread-local filter.
    - The thread-local filter receives no arguments.
      If an alignment filters is needed, the thread filter needs to return one.
    - The local filter finally receives an alignment that it can modify, and returns a boolean value,
      indicating whether the alignment is to be included in the final result or not.

    run-pipeline can be destructive to the input.
    This is indicated by the user with the :destructive keyword parameter that can have on of three values:
    - t: Operation is destructive.
    - nil: Operation is non-destructive. This may require copying of the input alignments
           but this only happens if there are any alignment-filters that may actually modify the alignments.
    - :default: Operation is destructive or non-destructive depending on what makes most sense for the given input:
                - For in-memory, :default means destructive.
                - For files, :default means non-destructive.
                - For databases:, :default means non-destructive.

    Specialize this generic function for new kinds of input sources.
    It is recommended to define two methods: One where only input is specialized, and another one where output
    is specialized as well for the same kind of input source, to take advantage of knowing the best way to run
    pipelines over the same input and output sources, and to identify cases where the pipeline has to effectively
    run in situ because input and output are the same."))

(defgeneric run-pipeline-in-situ (sam &key &allow-other-keys)
  (:documentation
   "Applies filters destructively and returns the first argument.
    This is like run-pipeline, except that input and output are the same.

    Specialize this generic function for new kinds of input/output sources."))

;; filter composition

(defun %compose-filters (filters &rest args)
  "Low-level filter composition, used by compose-global-filters and compose-thread-filters."
  (declare (list filters) (dynamic-extent args))
  (loop with local-filter with local-init with local-exit
        with local-exits of-type list
        for filter in filters do
        (setf (values local-filter local-init local-exit)
              (apply filter args))
        when local-filter collect (get-function local-filter) into local-filters
        when local-init collect (get-function local-init) into local-inits
        when local-exit do (push (get-function local-exit) local-exits)
        finally (return (values local-filters
                                (compose-thunks local-inits)
                                (compose-thunks local-exits)))))

(declaim (inline compose-global-filters))

(defun compose-global-filters (filters header)
  "Compose global (header) filters, and return thread-local filters, a global initialization thunk, and a global finalization thunk."
  (%compose-filters filters header))

(defmacro with-thread-filters ((thread-filters global-init global-exit) (filters header) &body body)
  "Macro version of compose-global-filters."
  `(multiple-value-bind
       (,thread-filters ,global-init ,global-exit)
       (compose-global-filters ,filters ,header)
     (flet ((,global-init () (funcall ,global-init))
            (,global-exit () (funcall ,global-exit)))
       (declare (inline ,global-init ,global-exit))
       ,@body)))

(declaim (inline compose-thread-filters))

(defun compose-thread-filters (thread-filters)
  "Compose thread-local filters, and return alignment filters, a thread-local initialization thunk, and a thread-local finalization thunk."
  (%compose-filters thread-filters))

(defmacro with-alignment-filters ((aln-filters local-init local-exit) (thread-filters) &body body)
  "Macro version of compose-thread-filters."
  `(multiple-value-bind
       (,aln-filters ,local-init ,local-exit)
       (compose-thread-filters ,thread-filters)
     (flet ((,local-init () (funcall ,local-init))
            (,local-exit () (funcall ,local-exit)))
       (declare (inline ,local-init ,local-exit))
       ,@body)))

(defgeneric create-chunk-filter (input-filter aln-filters output-filter destructive)
  (:documentation
   "Combines an optional alignment input mapper, a list of alignment filters, and an optional output mapper into a single function.
    Creates best combination based on which mappers and filters are present, and whether they are allowed to be destructive or not.
    May return nil if alignments don't need to be processed at all."))

(defmethod create-chunk-filter :around (input-filter (aln-filters cons) output-filter destructive)
  "Combines all alignment-filters into a single one."
  (let ((aln-filter (if (cdr aln-filters)
                      (lambda (aln)
                        (declare #.*optimization*)
                        (loop for filter in aln-filters
                              always (funcall (the function filter) aln)))
                      (car aln-filters))))
    (call-next-method input-filter aln-filter output-filter destructive)))

(defmethod create-chunk-filter ((input-filter null) (aln-filter null) (output-filter null) destructive)
  (declare (ignore input-filter aln-filter output-filter destructive))
  nil)

(defmethod create-chunk-filter ((input-filter t) (aln-filter null) (output-filter null) destructive)
  (declare (ignore aln-filter output-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nmapcar* input-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapcar* input-filter head tail))))

(defmethod create-chunk-filter ((input-filter null) (aln-filter t) (output-filter null) destructive)
  (declare (ignore input-filter output-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nfilter aln-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapfilter #'copy-structure aln-filter head tail))))

(defmethod create-chunk-filter ((input-filter null) (aln-filter null) (output-filter t) destructive)
  (declare (ignore input-filter aln-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nmapcar* output-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapcar* output-filter head tail))))

(defmethod create-chunk-filter ((input-filter t) (aln-filter t) (output-filter null) destructive)
  (declare (ignore output-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nmapfilter input-filter aln-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapfilter input-filter aln-filter head tail))))

(defmethod create-chunk-filter ((input-filter t) (aln-filter null) (output-filter t) destructive)
  (declare (ignore aln-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nmapcar* (lambda (aln)
                  (funcall (the function output-filter)
                           (funcall (the function input-filter) aln)))
                head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapcar* (lambda (aln)
                 (funcall (the function output-filter)
                          (funcall (the function input-filter) aln)))
               head tail))))

(defmethod create-chunk-filter ((input-filter null) (aln-filter t) (output-filter t) destructive)
  (declare (ignore input-filter))
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nfiltermap aln-filter output-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapfiltermap #'copy-structure aln-filter output-filter head tail))))

(defmethod create-chunk-filter ((input-filter t) (aln-filter t) (output-filter t) destructive)
  (if destructive
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (nmapfiltermap input-filter aln-filter output-filter head tail))
    (lambda (head &optional tail)
      (declare #.*optimization*)
      (mapfiltermap input-filter aln-filter output-filter head tail))))

;;; sorting order

(defun effective-sorting-order (sorting-order header original-sorting-order)
  "Determine effective sorting order: Some filters may destroy the sorting order recorded in the input.
   If this happens, and requested sorting order is :keep, then we need to effectively sort the result
   according to the original sorting order."
  (cond ((eq sorting-order :keep)
         (cond ((string= (getf (sam-header-hd header) :so "unknown") original-sorting-order) :keep)
               (t (sam-header-ensure-hd header :so original-sorting-order)
                  (intern-key (string-upcase original-sorting-order)))))
        (t (sam-header-ensure-hd header :so sorting-order)
           sorting-order)))

(declaim (inline sorting-criterion))

(defun sorting-criterion (sorting-order)
  "Determine sorting criterion for given sorting order, which is a list of parameters to be passed to sort/stable-sort."
  (case sorting-order
    (:coordinate (load-time-value (list #'coordinate<)))
    (:queryname (load-time-value (list #'string< :key #'sam-alignment-qname)))))

;;; generic filter thread setup

(defun call-with-threads (threads name filter-thread source-thread)
  "Set up the filter threads for a filter pipeline, and execute the input thread."
  (declare (fixnum threads) #.*optimization*)
  (loop repeat threads
        for mailbox = (make-mailbox 2)
        for thread = (thread-run name filter-thread mailbox)
        collect mailbox into mailboxes
        collect thread into threads
        finally (unwind-protect
                    (funcall (the function source-thread) mailboxes)
                  (loop for mailbox in mailboxes do (mailbox-send mailbox :eop))
                  (mapc #'thread-join threads))))

(defstruct (temporary-file (:constructor make-temporary-file (sibling)))
  "An internal representation for temporary files.
   The struct temporary-file has a constructor make-temporary-file that takes a sibling file as a parameter.
   Accessor temporary-file-sibling refers to the sibling.
   Primary use of this struct is to create temporary intermediate result files in run-pipeline-in-situ."
  sibling)

(setf (documentation 'make-temporary-file 'function)
      "Constructor for struct temporary-file that takes a sibling file as a parameter."
      (documentation 'temporary-file-p 'function)
      "Default predicate for struct temporary-file."
      (documentation 'copy-temporary-file 'function)
      "Default copier function for struct temporary-file."
      (documentation 'temporary-file-sibling 'function)
      "Access the temporary-file sibling.")

(declaim (inline check-file-sorting-order))

(defun check-file-sorting-order (sorting-order)
  "Assert that no sorting is asked for files."
  (ecase sorting-order ((:keep :unknown :unsorted))))

(defgeneric get-output-functions (output header original-sq &key sorting-order &allow-other-keys)
  (:documentation
   "Sets up the output thread, feeds the given header to the output data set, and returns four functions:
    One for potential additional output mapping of alignments (output-filter, can be NIL, which defaults to #'IDENTITY),
    one for invoking the function that filters chunks (wrap-thread, can be NIL, which defaults to #'FUNCALL),
    one for receiving chunks of alignments in filter threads (receive-chunk),
    and one for finalizing alignments in the main thread (finalize).
    The function wrap-thread is called once per thread, receives the main function for processing chunks as a thunk,
    and has to eventually call that thunk. The purpose of wrap-thread is to enable setting up additional resources
    to be used in output-filter and/or receive-chunk.
    When the (already effective) :sorting-order is :coordinate or :queryname, sort the received chunks before
    writing to output. Sorting is not supported for files, but by default only for in-memory representations.
    The parameter original-sq should be used to query the reference sequence dictionary,
    because the one in header may have already been filtered and modified.

    Specialize this generic function for new kinds of output sources."))

(defun chunk-output-loop (mailbox ordered output-chunk)
  "Abstract loop for writing chunks of alignments. Each chunk has a sequence number prepended to it.
   If ordered is true, ensure that chunks are written in the order of the sequence numbers.
   Otherwise ignore the sequence numbers."
  (declare (mailbox mailbox) (function output-chunk) #.*optimization*)
  (if ordered
    (loop with stash = '() with run of-type integer = 0
          for chunk = (mailbox-read mailbox) until (eq chunk :eop)
          if (= (car chunk) run) do
          (loop do (funcall output-chunk (cdr chunk))
                (incf run)
                (setq chunk (car stash))
                while (and chunk (= (car chunk) run))
                do (setq stash (cdr stash)))
          else do (setq stash (merge 'list (list chunk) stash #'< :key #'car))
          finally (assert (null stash)))
    (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop)
          do (funcall output-chunk (cdr chunk)))))

(defmacro with-chunk-output ((chunk) (mailbox ordered) &body body)
  "Macro version of chunk-output-loop."
  `(chunk-output-loop ,mailbox ,ordered (lambda (,chunk) ,@body)))

(defmethod get-output-functions ((output sam) header original-sq &key (sorting-order :keep) (split-file nil))
  (ecase sorting-order
    ((:keep :unknown :unsorted)
     (let* ((head (cons nil nil))
            (tail head)
            (mailbox (make-mailbox *number-of-threads*))
            (thread (thread-run
                     "in-memory output thread"
                     (lambda ()
                       (setf (sam-header output) header)
                       (with-chunk-output (chunk) (mailbox (not (eq sorting-order :unsorted)))
                         (when chunk (setf (cdr tail) chunk tail (last chunk))))))))
       (values nil nil
               (lambda (chunk)
                 (mailbox-send mailbox chunk))
               (lambda ()
                 (mailbox-send mailbox :eop)
                 (thread-join thread)
                 (setf (sam-alignments output) (cdr head))
                 output))))
    ((:coordinate)
     (let* ((threads *number-of-threads*)
            (table   (make-array (if split-file threads                  ; distribute evenly across positions
                                   (1+ (length (sam-header-sq header)))) ; or distribute over reference sequences
                                 :initial-element '() #+lispworks :single-thread #+lispworks t))
            (mailbox (make-mailbox threads))
            (thread  (thread-run
                      "in-memory output thread with sorting by coordinate"
                      (lambda ()
                        (declare #.*optimization*)
                        (setf (sam-header output) header)
                        (if (= (length table) 1) ; only one bucket, so no bucket distribution necessary (for example, for unmapped split files)
                          (with-chunk-output (chunk) (mailbox nil)
                            (when chunk (setf (cdr (last chunk)) (svref table 0) (svref table 0) chunk)))
                          (macrolet ((receive-chunks (index-form) ; distribute over buckets according to index-form
                                       `(with-chunk-output (chunk) (mailbox nil)
                                          (loop while chunk do
                                                (let* ((cons chunk)
                                                       (aln (car cons))
                                                       (index ,index-form))
                                                  (declare (cons chunk cons) (sam-alignment aln) (fixnum index))
                                                  (setf chunk (cdr chunk)
                                                        (cdr cons) (svref table index)
                                                        (svref table index) cons))))))
                            (if split-file
                              (let ((min-pos -1) (bucket-size 1))
                                (declare (int32 min-pos) (rational bucket-size))
                                (flet ((get-table-index (aln)
                                         (declare (sam-alignment aln) #.*optimization*)
                                         (when (< min-pos 0)
                                           (let ((refid (sam-alignment-refid aln)))
                                             (declare (fixnum refid))
                                             (if (< refid 0)
                                               (return-from get-table-index 0)
                                               (let* ((sq (loop with sn = (getf (nth refid (sam-header-sq header)) :sn)
                                                                for entry in original-sq
                                                                when (string= sn (getf entry :sn)) return entry))
                                                      (len (or (getf sq :ln) #.(1- (expt 2 31)))))
                                                 (declare (int32 len))
                                                 (setq min-pos 1 bucket-size (/ len threads))))))
                                         (floor (the int32 (- (sam-alignment-pos aln) min-pos)) bucket-size)))
                                  (declare (inline get-table-index))
                                  (receive-chunks (the fixnum (get-table-index aln)))))
                              (receive-chunks (the fixnum (1+ (the int32 (sam-alignment-refid aln))))))))))))
       (assert (> (length table) 0))
       (values nil nil
               (lambda (chunk)
                 (mailbox-send mailbox chunk))
               (lambda ()
                 (mailbox-send mailbox :eop)
                 (thread-join thread)
                 (cond ((= (length table) 1) ; only one bucket
                        (let ((first-aln (car (svref table 0))))
                          (when first-aln
                            (unless (< (sam-alignment-refid first-aln) 0) ; only unmapped reads => no sorting necessary
                              (setf (svref table 0)
                                    (stable-sort (svref table 0) #'< :key #'sam-alignment-pos)))))
                        (setf (sam-alignments output) (svref table 0)))
                       (t (claws:reset-workers threads)
                          (labels ((recur (min max)
                                     (declare (fixnum min max) #.*optimization*)
                                     (let ((length (- max min)))
                                       (declare (fixnum length))
                                       (cond ((= length 0))
                                             ((= length 1)
                                              (setf (svref table min)
                                                    (stable-sort (svref table min) #'< :key #'sam-alignment-pos)))
                                             (t (let* ((half (ash length -1))
                                                       (mid (+ min half)))
                                                  (declare (fixnum half mid))
                                                  (claws:spawn () (recur mid max))
                                                  (recur min mid)
                                                  (claws:sync)))))))
                            (recur 0 (length table)))
                          (claws:reset-workers 1)
                          (setf (sam-alignments output)
                                (if split-file
                                  (loop for list across table nconc list)
                                  (reduce #'nconc table :from-end t :start 1 :initial-value (svref table 0))))))
                 output))))
    ((:queryname)
     (let* ((tree    (make-simple-tree 16))
            (mailbox (make-mailbox *number-of-threads*))
            (thread  (thread-run
                      "in-memory output thread with sorting by queryname"
                      (lambda ()
                        (setf (sam-header output) header)
                        (with-chunk-output (chunk) (mailbox nil)
                          (when chunk (setq tree (insert-node tree chunk))))))))
       (values nil nil
               (lambda (chunk)
                 (mailbox-send mailbox chunk))
               (lambda ()
                 (mailbox-send mailbox :eop)
                 (thread-join thread)
                 (let ((sorting-criterion (sorting-criterion sorting-order)))
                   (setf (sam-alignments output)
                         (tree-reduce tree *number-of-threads*
                                      (lambda (chunk)
                                        (apply #'stable-sort chunk sorting-criterion))
                                      (lambda (left right)
                                        (apply #'merge 'list left right sorting-criterion)))))
                 output))))))

(defmethod get-output-functions ((output pathname) header original-sq &key (sorting-order :keep))
  (declare (ignore original-sq))
  (check-file-sorting-order sorting-order)
  (let* ((mailbox (make-mailbox *number-of-threads*))
         (thread  (thread-run
                   "file output thread"
                   (lambda ()
                     (with-open-sam (out output :direction :output)
                       (format-sam-header out header)
                       (with-chunk-output (chunk) (mailbox (not (eq sorting-order :unsorted)))
                         (loop for aln in chunk do (writeln out aln))))))))
    (values (lambda (aln)
              (let ((sim (make-sim-stream (estimate-sam-alignment-output-length aln))))
                (declare (dynamic-extent sim))
                (sim-format-sam-alignment sim aln)
                (sim-stream-string sim)))
            nil
            (lambda (chunk)
              (mailbox-send mailbox chunk))
            (lambda ()
              (mailbox-send mailbox :eop)
              (thread-join thread)
              output))))

(defmethod get-output-functions ((output temporary-file) header original-sq &key (sorting-order :keep))
  (declare (ignore original-sq))
  (check-file-sorting-order sorting-order)
  (let ((mailbox (make-mailbox *number-of-threads*))
        (outbox  (make-mailbox)))
    (thread-run
     "temporary file output thread"
     (lambda ()
       (multiple-value-bind
           (out program pathname)
           (open-temporary-sam (temporary-file-sibling output))
         (unwind-protectn
           (format-sam-header out header)
           (with-chunk-output (chunk) (mailbox (not (eq sorting-order :unsorted)))
             (loop for aln in chunk do (writeln out aln)))
           (close-sam out program))
         (mailbox-send outbox pathname))))
    (values (lambda (aln)
              (let ((sim (make-sim-stream (estimate-sam-alignment-output-length aln))))
                (declare (dynamic-extent sim))
                (sim-format-sam-alignment sim aln)
                (sim-stream-string sim)))
            nil
            (lambda (chunk)
              (mailbox-send mailbox chunk))
            (lambda ()
              (mailbox-send mailbox :eop)
              (mailbox-read outbox)))))

(defmethod get-output-functions ((output null) header original-sq &key sorting-order)
  (declare (ignore header original-sq sorting-order))
  (values nil (constantly nil) (constantly nil)))

(defmacro with-output-functions ((output-filter wrap-thread receive-chunk) form &body body)
  "Macro for receiving results from get-output-functions."
  (let ((finalize (copy-symbol 'finalize)))
    `(multiple-value-bind
         (,output-filter ,wrap-thread ,receive-chunk ,finalize) ,form
       (when ,output-filter (setq ,output-filter (get-function ,output-filter)))
       (unless ,wrap-thread (setq ,wrap-thread #'funcall))
       (setq ,receive-chunk (get-function ,receive-chunk))
       (flet ((,output-filter (aln)
                (declare #.*optimization*)
                (funcall (the function ,output-filter) aln))
              (,wrap-thread (thunk) (funcall ,wrap-thread thunk))
              (,receive-chunk (chunk)
                (declare #.*optimization*)
                (funcall (the function ,receive-chunk) chunk)))
         (declare (inline ,output-filter ,wrap-thread ,receive-chunk))
         (locally ,@body)
         (funcall ,finalize)))))

;;; filters and pipelines

(defconstant +default-chunk-size+ #.(expt 2 13)
  "The number of alignments to read from an input source at a time.
   Default is 8192.")

(defmethod run-pipeline-in-situ :before ((sam t) &key (destructive t))
  (assert (eq destructive t)))

(defmethod run-pipeline-in-situ ((sam sam) &rest args)
  (declare (dynamic-extent args))
  (let ((input (make-sam :header (sam-header sam)
                         :alignments (sam-alignments sam))))
    (declare (dynamic-extent input))
    (apply #'run-pipeline input sam :destructive t args)))

(defmethod run-pipeline-in-situ ((sam pathname) &rest args &key header filters (sorting-order :keep))
  (declare (dynamic-extent args))
  (if (and (null header) (null filters))
    (check-file-sorting-order sorting-order)
    (rename-file (apply #'run-pipeline sam (make-temporary-file sam) :destructive t args)
                 sam))
  sam)

(defmacro with-prepared-header ((header original-sorting-order original-sq alns) (input destructive) &body body)
  "Prepare the header for further processing, when input is in memory."
  (let ((inputv (copy-symbol 'input)))
    `(let* ((,inputv ,input)
            (,header (sam-header ,inputv))
            (,alns (sam-alignments ,inputv))
            (,original-sorting-order (getf (sam-header-hd ,header) :so "unknown"))
            (,original-sq (sam-header-sq ,header)))
       (declare (ignorable ,original-sq))
       (if ,destructive
         (setf (sam-alignments ,inputv) '())
         (setq ,header (copy-structure ,header)))
       ,@body)))

(defmethod run-pipeline ((input sam) (output sam) &key filters (sorting-order :keep) (destructive :default))
  "Optimize when both input and output are in memory."
  (if (> *number-of-threads* 3)
    ;; enable multithreading
    (call-next-method)
    ;; simplified case: everything in a single thread when copying completely in-memory
    (with-prepared-header (header original-sorting-order original-sq alns) (input destructive)
      (with-thread-filters (thread-filters global-init global-exit) (filters header)
        (setf (sam-header output) header)
        (if (null alns)
          (setf (sam-alignments output) '())
          (let* ((sorting-order (effective-sorting-order sorting-order header original-sorting-order))
                 (sorting-criterion (sorting-criterion sorting-order)))
            (flet ((no-filter-case () (setf (sam-alignments output)
                                            (if sorting-criterion
                                              (apply #'stable-sort
                                                     (if destructive alns (copy-list alns))
                                                     sorting-criterion)
                                              alns))))
              (cond ((null thread-filters) (no-filter-case))
                    (t (global-init)
                       (unwind-protect
                           (with-alignment-filters (aln-filters local-init local-exit) (thread-filters)
                             (let ((chunk-filter (create-chunk-filter nil aln-filters nil destructive)))
                               (cond ((null chunk-filter) (no-filter-case))
                                     (t (local-init)
                                        (unwind-protectn
                                          (setq alns (funcall chunk-filter alns))
                                          (setf (sam-alignments output)
                                                (if sorting-criterion
                                                  (apply #'stable-sort alns sorting-criterion)
                                                  alns))
                                          (local-exit))))))
                         (global-exit))))))))
      output)))

(defmacro with-sam-chunk (sam-chunk &body body)
  "Access the components of a chunk of sam-alignment instances."
  (let ((cons (copy-symbol 'cons)))
    `(let* ((,cons (cdr ,sam-chunk))
            (head (car ,cons))
            (tail (cdr ,cons)))
       ,@body)))

(defmethod run-pipeline ((input sam) output &rest args &key filters (sorting-order :keep) (destructive :default) (chunk-size +default-chunk-size+))
  (declare (dynamic-extent args))
  (with-prepared-header (header original-sorting-order original-sq alns) (input destructive)
    (with-thread-filters (thread-filters global-init global-exit) (filters header)
      (setq sorting-order (effective-sorting-order sorting-order header original-sorting-order))
      (with-output-functions
          (output-filter wrap-thread receive-chunk)
          (apply #'get-output-functions output header original-sq :sorting-order sorting-order args)
        (when alns
          (cond ((and (null thread-filters)
                      (null output-filter))
                 (wrap-thread
                  (lambda ()
                    (ecase sorting-order
                      ((:keep :unknown :unsorted)
                       (receive-chunk (cons 0 alns)))
                      ((:coordinate :queryname)
                       (if destructive
                         (loop with chunk-size-1 = (1- chunk-size)
                               do (let ((head alns) (tail (nthcdr chunk-size-1 alns)))
                                    (if tail
                                      (setf alns (cdr tail) (cdr tail) nil)
                                      (setq alns nil))
                                    (receive-chunk (cons 0 head)))
                               while alns)
                         (loop do (multiple-value-bind
                                      (head tail) (nthdiff chunk-size alns)
                                    (setq alns tail)
                                    (receive-chunk (cons 0 head)))
                               while alns)))))))
                (t
                 ; source thread:
                 ;   get chunks out of input list; send chunks with head and tail - no copies are created in source thread
                 ; filter threads:
                 ;   if destructive, mark last cdr as nil, modify in-situ, send to target
                 ;   if non-destructive, copy then modify, send to target
                 (global-init)
                 (unwind-protect
                     (call-with-threads
                      (max 1 (- *number-of-threads* 2))
                      "in-memory filter thread"
                      (lambda (mailbox)
                        (wrap-thread
                         (lambda ()
                           (with-alignment-filters (aln-filters local-init local-exit) (thread-filters)
                             (let ((chunk-filter (create-chunk-filter nil aln-filters output-filter destructive)))
                               (cond (chunk-filter
                                      (macrolet ((chunk-filter (&rest args)
                                                   `(locally (declare #.*optimization*)
                                                      (funcall (the function chunk-filter) ,@args))))
                                        (local-init)
                                        (unwind-protect
                                            (if destructive
                                              (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do
                                                    (with-sam-chunk chunk
                                                      (when tail (setf (cdr tail) nil))
                                                      (setf (cdr chunk) (chunk-filter head))
                                                      (receive-chunk chunk)))
                                              (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do
                                                    (with-sam-chunk chunk
                                                      (setf (cdr chunk) (chunk-filter head (cdr tail)))
                                                      (receive-chunk chunk))))
                                          (local-exit))))
                                     (t (if destructive
                                          (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do
                                                (with-sam-chunk chunk
                                                  (when tail (setf (cdr tail) nil))
                                                  (setf (cdr chunk) head)
                                                  (receive-chunk chunk)))
                                          (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do
                                                (with-sam-chunk chunk
                                                  (setf (cdr chunk) (ldiff head (cdr tail)))
                                                  (receive-chunk chunk)))))))))))
                      (lambda (mailboxes)
                        (loop with serial = -1
                              with chunk-size-1 = (1- chunk-size)
                              for mailbox-ring = mailboxes then (or (cdr mailbox-ring) mailboxes)
                              do (let ((head alns) (tail (nthcdr chunk-size-1 alns)))
                                   (setq alns (cdr tail))
                                   (mailbox-send (car mailbox-ring) (cons (incf serial) (cons head tail))))
                              while alns)))
                   (global-exit)))))))))

(defmethod run-pipeline ((input pathname) (output pathname) &rest args &key filters (sorting-order :keep) (destructive :default))
  "Optimize when both input and output are files."
  (declare (dynamic-extent args))
  (let ((input-name (truename input)))
    (let ((output-name (probe-file output)))
      (when output-name
        (when (equal input-name output-name)
          (return-from run-pipeline (apply #'run-pipeline-in-situ input :destructive t args))))))
  (when (string-equal (pathname-type input) (pathname-type output))
    (unless filters
      (check-file-sorting-order sorting-order)
      (if (or (not destructive) (eq destructive :default))
        #+lispworks (lw:copy-file input output)
        #+sbcl (cl-fad:copy-file input output :overwrite t)
        (rename-file input output))
      (return-from run-pipeline output)))
  (call-next-method))

(defmethod run-pipeline ((input pathname) output &rest args &key filters (sorting-order :keep) (destructive :default) (chunk-size +default-chunk-size+))
  (declare (dynamic-extent args))
  (unwind-protect
      (with-open-sam (in input :direction :input)
        (let* ((header (parse-sam-header in))
               (original-sorting-order (getf (sam-header-hd header) :so "unknown"))
               (original-sq (sam-header-sq header)))
          (with-thread-filters (thread-filters global-init global-exit) (filters header)
            (setq sorting-order (effective-sorting-order sorting-order header original-sorting-order))
            (when (null thread-filters)
              (typecase output
                (pathname       (check-file-sorting-order sorting-order)
                                (with-open-sam (out output :direction :output)
                                  (format-sam-header out header)
                                  (copy-stream in out))
                                (return-from run-pipeline output))
                (temporary-file (check-file-sorting-order sorting-order)
                                (multiple-value-bind
                                    (out program pathname)
                                    (open-temporary-sam (temporary-file-sibling output))
                                  (unwind-protectn
                                    (format-sam-header out header)
                                    (copy-stream in out)
                                    (close-sam out program))
                                  (return-from run-pipeline pathname)))))
            (with-output-functions
                (output-filter wrap-thread receive-chunk)
                (apply #'get-output-functions output header original-sq :sorting-order sorting-order args)
              (global-init)
              (unwind-protect
                  (call-with-threads
                   (max 1 (- *number-of-threads* 2))
                   "file filter thread"
                   (lambda (mailbox)
                     (wrap-thread
                      (lambda ()
                        (with-alignment-filters (aln-filters local-init local-exit) (thread-filters)
                          (if (and (null aln-filters) (or (typep output 'pathname)
                                                          (typep output 'temporary-file)))
                            (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do (receive-chunk chunk))
                            (let ((chunk-filter (create-chunk-filter #'parse-sam-alignment aln-filters output-filter t)))
                              (macrolet ((chunk-filter (&rest args)
                                           `(locally (declare #.*optimization*)
                                              (funcall (the function chunk-filter) ,@args))))
                                (local-init)
                                (unwind-protect
                                    (loop for chunk = (mailbox-read mailbox) until (eq chunk :eop) do
                                          (setf (cdr chunk) (chunk-filter (cdr chunk)))
                                          (receive-chunk chunk))
                                  (local-exit)))))))))
                   (lambda (mailboxes)
                     (loop with serial = -1
                           for mailbox-ring = mailboxes then (or (cdr mailbox-ring) mailboxes)
                           for alns = (loop repeat chunk-size
                                            for aln = (read-line in nil)
                                            while aln collect aln)
                           while alns do (mailbox-send (car mailbox-ring) (cons (incf serial) alns)))))
                (global-exit))))))
    (unless (or (not destructive) (eq destructive :default))
      (delete-file input))))
