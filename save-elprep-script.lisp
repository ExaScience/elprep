(in-package :cl-user)

(let ((quicklisp-init (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname))))
  (when (probe-file quicklisp-init)
	(load quicklisp-init)))

(asdf:load-system :elprep)

(deliver 'elprep:elprep-script "elprep" 4 :multiprocessing t)
