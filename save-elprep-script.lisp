(in-package :cl-user)

(let ((quicklisp-init (print (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname)))))
  (when (probe-file quicklisp-init)
	(load quicklisp-init)))

(asdf:load-system :elprep)

#+lispworks
(deliver #'elprep:elprep-script "elprep" 4 :multiprocessing t)

#+sbcl
(sb-ext:save-lisp-and-die "elprep" :toplevel #'elprep:elprep-script :executable t :save-runtime-options t)
