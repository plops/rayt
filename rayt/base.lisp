(provide :base)
(unless (find-package :base)
    (make-package :base :use '(:cl))) 

(in-package :base)
(export 
 '(num vec with-arrays v copy-vec
   .+ .- .* ./ dot norm normalize
   cross .s vx vy vz v-spherical
   check-unit-vector check-range
   req +pif+ +2pif+))

(declaim (inline v copy-vec .+ .- .* ./ dot norm normalize .s))

(deftype num ()
  `single-float)

(defconstant +pif+ #.(coerce pi 'num))
(defconstant +2pif+ #.(coerce (* 2 pi) 'num))

(deftype vec ()
  `(simple-array num 1))

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(eval-when (:compile-toplevel)
 (defun num (x)
   (declare (number x))
   (coerce x 'single-float)))

(defmacro vec (&rest rest)
  (let ((a (gensym)))
   `(let ((,a (make-array ,(list-length rest)
			 :element-type 'num)))
      ,@(let ((i 0))
	     (declare (type (integer 0 #.(1- (expt 2 16))) i))
	     (loop for e in rest collect
		  (prog1
		       `(setf (aref ,a ,i)
			     ,(typecase e ;; FIXME once-only e?
					(fixnum (* 1s0 e))
					(num e)
					(t `(* 1s0 ,e))))
		     (incf i))))
      ,a)))

#+nil 
(vec 2 2 1)
#+nil
(vec .2 .3 1)
(defun v (&optional
	  (z 0s0)
	  (y 0s0)
	  (x 0s0))
  (the vec (vec z y x)))


(declaim (ftype (function (vec) (values vec &optional)) copy-vec))
(defun copy-vec (a)
  (let* ((n (length a))
	 (b (make-array n
			:element-type (array-element-type a))))
    (dotimes (i n)
      (setf (aref b i) (aref a i)))
    b))
#+nil
(copy-vec (v))

(progn
  (declaim (ftype (function (vec &rest t) (values vec &optional))
		  .+ .- .* ./))
  (defun .+ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (vec e))
	(dotimes (i (length r))
	  (incf (aref r i) (aref e i))))
      r))
  (defun .- (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (vec e))
	(dotimes (i (length r))
	  (decf (aref r i) (aref e i))))
     r))
  (defun .* (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (vec e))
	(dotimes (i (length r))
	 (setf (aref r i) (* (aref r i) (aref e i)))))
     r))
  (defun ./ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(declare (vec e))
	(dotimes (i (length r))
	  (setf (aref r i) (/ (aref r i) (aref e i)))))
      r)))
#+nil
(./ (vec 1 2 3) (vec 3 32 2) (vec 32 4 2))

(defun dot (a b)
  (declare (vec a b))
  (let ((r 0s0))
    (declare (num r))
    (dotimes (i (length a))
      (incf r (* (aref a i) (aref b i))))
    (the num r)))
#+nil
(dot (v 1 3 0) (v 3))

(defun norm (v)
  (declare (vec v))
  (let ((l2 (dot v v)))
    (declare (type (single-float 0s0) l2)) ;; FIXME: write num here
    (the num (sqrt l2))))
#+nil
(norm (v 1 1 0))

(declaim (ftype (function (num vec) (values vec &optional)) .s))
(defun .s (s a)
  "scalar multiplication"
  (let ((r (v)))
    (dotimes (i (length a))
      (setf (aref r i) (* s (aref a i))))
    r))
#+nil
(.s .3 (v 1 1))


(defun normalize (v)
  (declare (vec v))
  (let ((b (copy-vec v)))
    (the vec (.s (/ (norm v)) b))))

(declaim (ftype (function (vec vec) (values vec &optional)) cross))
(defun cross (a b)
  (let ((r (v)))
    (with-arrays (r a b)
      (setf (r 2) (- (* (a 1) (b 0))
		     (* (a 0) (b 1)))
	    (r 1) (- (* (a 0) (b 2))
		     (* (a 2) (b 0)))
	    (r 0) (- (* (a 2) (b 1))
		     (* (a 1) (b 2)))))
    r))
#+nil
(cross (v 0 0 1) (v 0 1 0))


(declaim (inline vx vy vz))
(defun vx (v)
  (declare (vec v))
  (the num (aref v 2)))

(defun vy (v)
  (declare (vec v))
  (the num (aref v 1)))

(defun vz (v)
  (declare (vec v))
  (the num (aref v 0)))

(defun v-spherical (theta phi)
  "Convert spherical coordinates into cartesian."
  (declare (type (single-float 0s0 #.(/ (coerce pi 'num) 4)) theta)
	   (type (single-float #.(coerce (- pi) 'num) #.(coerce pi 'num)) phi))
  (let* ((st (sin theta)))
    (the vec
      (v (cos theta)
	 (* st (sin phi))
	 (* st (cos phi))))))


(defun check-unit-vector (&rest rest)
  (declare (ignore rest)))

(defun check-range (min max &rest rest)
  (declare (ignore min max rest)))


;; (defun check-unit-vector (&rest rest)
;;   ;; turn off for faster speed
;;   (dolist (e rest)
;;     (unless (< (abs (- (norm e) 1)) 1s-6)
;;       (error "vector isn't normalized"))))

;; (defun check-range (min max &rest rest)
;; #+nil  (declare (num min max))
;;   (dolist (e rest)
;;     (declare (num e))
;;     (unless (< min e max)
;;       (error "range check failed"))))

(defun req (&optional name)
  (error "Required argument ~@[~S~] missing" name))
