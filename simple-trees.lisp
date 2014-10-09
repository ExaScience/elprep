(in-package :elprep)

(declaim (inline make-simple-tree))

(defstruct (simple-tree 
            (:constructor make-simple-tree
             (rank &optional depth &aux (nodes (make-array rank :initial-element nil
                                                           #+lispworks :single-thread #+lispworks t)))))
  "A simple balanced tree.
   The struct simple-tree has a constructor that takes the rank and optionally the depth as parameters.
   Accessor simple-tree-index of type fixnum points to the current node in simple-tree-nodes.
   Read-only accessor simple-tree-nodes of type simple-vector refers to the nodes or leaves.
   Read-only accessor simple-tree-depth of type fixnum refers to the depth of the tree.
   Read-only accessor simple-tree-rank of type fixnum refers to the rank of the tree."
  (index 0 :type fixnum)
  (nodes #() :type simple-vector :read-only t)
  (depth 0 :type fixnum :read-only t)
  (rank 0 :type fixnum :read-only t))

(setf (documentation 'make-simple-tree 'function)
      "Constructor for struct simple-tree that takes the rank and optionally the depth as parameters."
      (documentation 'simple-tree-p 'function)
      "Default predicate for struct simple-tree."
      (documentation 'copy-simple-tree 'function)
      "Default copier function for struct simple-tree."
      (documentation 'simple-tree-index 'function)
      "Access the simple-tree index of type fixnum."
      (documentation 'simple-tree-nodes 'function)
      "Read the simple-tree nodes of type simple-vector."
      (documentation 'simple-tree-depth 'function)
      "Read the simple-tree depth of type fixnum."
      (documentation 'simple-tree-rank 'function)
      "Read the simple-tree rank of type fixnum.")

(declaim (inline make-subtree))

(defun make-subtree (tree)
  "Create a subtree for the given simple-tree."
  (declare (simple-tree tree) #.*optimization*)
  (make-simple-tree
   (simple-tree-rank tree)
   (the fixnum (1- (simple-tree-depth tree)))))

(declaim (inline make-super-tree))

(defun make-super-tree (tree)
  "Create a super tree for the given simple-tree, and add the given tree to the new super tree as its first node."
  (declare (simple-tree tree) #.*optimization*)
  (let ((super-tree
         (make-simple-tree
          (simple-tree-rank tree)
          (the fixnum (1+ (simple-tree-depth tree))))))
    (declare (simple-tree super-tree))
    (setf (simple-tree-index super-tree) 1)
    (setf (svref (simple-tree-nodes super-tree) 0) tree)
    super-tree))

(defun insert-node (top-tree node)
  "Insert a node into the given simple-tree. Returns the tree unless it is fully occupied.
   If it is fully occupied, create a fresh super tree, insert both the given tree and the node there,
   and return that super tree."
  (declare (simple-tree top-tree) #.*optimization*)
  (labels ((recur (tree)
             (declare (simple-tree tree))
             (let ((index (simple-tree-index tree))
                   (nodes (simple-tree-nodes tree)))
               (declare (fixnum index) (simple-vector nodes))
               (when (< index (length nodes))
                 (when (= (simple-tree-depth tree) 0)
                   (setf (svref nodes index) node)
                   (setf (simple-tree-index tree) (the fixnum (1+ index)))
                   ;; done: jump out of the recursion!
                   (return-from insert-node top-tree))
                 ;; depth /= 0
                 (recur (or (svref nodes index) (setf (svref nodes index) (make-subtree tree))))
                 ;; tree is full, try to create a new subtree
                 (setf (simple-tree-index tree) (incf index))
                 (when (< index (length nodes))
                   (recur (setf (svref nodes index) (make-subtree tree))))))))
    (recur top-tree)
    ;; tree is full, try from a new top tree
    (recur (setq top-tree (make-super-tree top-tree)))
    ;; new top tree is also full, which should never happen
    (error "This code should not be reached.")))

(defun tree-reduce (tree threads map reduce)
  "Perform a parallel map/reduce traversal over the given simple-tree."
  (declare (simple-tree tree) (fixnum threads) (function map reduce) #.*fixnum-optimization*)
  (claws:reset-workers threads)
  (unwind-protect
      (labels ((reduce-vector (vector start end map)
                 (declare (simple-vector vector) (fixnum start end) (function map))
                 (let ((length (- end start)))
                   (declare (fixnum length))
                   (cond ((= length 0))
                         ((= length 1)
                          (funcall map (svref vector start)))
                         (t (let* ((half (ash length -1))
                                   (middle (+ start half))
                                   left right)
                              (declare (fixnum half middle))
                              (claws:spawn (left) (reduce-vector vector start middle map))
                              (setq right (reduce-vector vector middle end map))
                              (claws:sync)
                              (funcall reduce left right))))))
               (recur-tree (tree)
                 (declare (simple-tree tree))
                 (let ((index (simple-tree-index tree))
                       (nodes (simple-tree-nodes tree)))
                   (declare (fixnum index) (simple-vector nodes))
                   (if (= (simple-tree-depth tree) 0)
                     (reduce-vector nodes 0 index map)
                     (reduce-vector nodes 0 (if (and (< index (length nodes))
                                                     (svref nodes index))
                                              (1+ index) index)
                                    #'recur-tree)))))
        (recur-tree tree))
    (claws:reset-workers 1)))
