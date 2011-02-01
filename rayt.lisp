(deftype num ()
  `single-float)

(deftype vec ()
  `(simple-array num 1))

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(defun num (x)
  (coerce x 'num))

(defmacro vec (&rest rest)
  (let ((a (gensym)))
   `(let ((,a (make-array ,(length rest)
			 :element-type 'num)))
      ,@(let ((i -1))
	     (loop for e in rest collect
		  `(setf (aref ,a ,(incf i))
			 ,(typecase e ;; FIXME once-only e?
				    (number (coerce e 'num))
				    (t `(* (num 1) ,e))))))
      ,a)))
#+nil
(vec 2 2 1)

(defun v (&optional
	  (z (num 0))
	  (y (num 0))
	  (x (num 0)))
  (vec z y x))

(declaim (ftype (function (vec) (values vec &optional)) copy-vec))
(defun copy-vec (a)
  (let* ((n (length a))
	 (b (make-array n
			:element-type (array-element-type a))))
    (dotimes (i n)
      (setf (aref b i) (aref a i)))
    b))
#+nil
(copy-vec (vec3))

(progn
  (declaim (ftype (function (vec &rest t) (values vec &optional))
		  .+ .- .* ./))
  (defun .+ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(dotimes (i (length r))
	  (incf (aref r i) (aref e i))))
      r))
  (defun .- (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(dotimes (i (length r))
	  (decf (aref r i) (aref e i))))
     r))
  (defun .* (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(dotimes (i (length r))
	 (setf (aref r i) (* (aref r i) (aref e i)))))
     r))
  (defun ./ (a &rest rest)
    (let ((r (copy-vec a)))
      (dolist (e rest)
	(dotimes (i (length r))
	  (setf (aref r i) (/ (aref r i) (aref e i)))))
      r)))
#+nil
(./ (vec 1 2 3) (vec 3 32 2) (vec 32 4 2))

(progn
  (declaim (ftype (function (vec vec) (values num &optional)) dot))
  (defun dot (a b)
    (let ((r (num 0)))
      (dotimes (i (length a))
	(incf r (* (aref a i) (aref b i))))
      r)))
#+nil
(dot (v 1 3 0) (v 3))

(defun norm (v)
  (declare (vec v))
  (the num (sqrt (dot v v))))
#+nil
(norm (v 1 1 0))

(defun normalize (v)
  (declare (vec v))
  (let ((b (copy-vec v)))
    (the vec (.s (/ (norm v)) b))))

(progn
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
      r)))
#+nil
(cross (v 0 0 1) (v 0 1 0))

(progn
  (declaim (ftype (function (num vec) (values vec &optional)) .s))
  (defun .s (s a)
    "scalar multiplication"
    (let ((r (v)))
      (dotimes (i (length a))
	(setf (aref r i) (* s (aref a i))))
      r)))
#+nil
(.s .3 (v 1 1))

(progn
  (declaim (ftype (function (num num num) (values num &optional))
		  quadratic-root-dist))
  (declaim (optimize (speed 3) (safety 0) (debug 0)))
  (defun quadratic-root-dist (a b c)
    (let ((eps (coerce 1s-12 'num))
	  (zero (coerce 0 'num))
	  (det2 (- (* b b) (* 4 a c))))
     (if (<= det2 zero)
       zero
       (let* ((q (* .5 (+ b (* (signum b) (sqrt det2)))))
	      (aa (abs a))
	      (aq (abs q)))
	 (cond ((or (< aq eps) (< aa eps)) zero)
	       (t (- (/ q a) (/ c q)))))))))
#+nil
(quadratic-root-dist 1s0 2s0 -3s0)

(defun find-inverse-ray-angle (height focal-length)
  "HEIGHT meridional distance of a point in sample space towards optic
axis. Coordinates in mm."
  (declare (num rho focal-length))
  (the num (asin (/ height focal-length))))
#+nil
(find-inverse-ray-angle 2.2 2.61)

(defun find-bfp-radius (na f)
  (declare (num na f)) (the num (* f na)))

(defun find-focal-length (mag)
  (declare (num mag))  (the num (/ 164.5 mag)))

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
  (declare ((double-float 0d0 #.(/ pi 4)) theta)
	   ((double-float 0d0 #.(* 2 pi)) phi)
	   (values vec &optional))
  (let* ((st (sin theta)))
    (the vec
      (v (* st (cos phi))
	 (* st (sin phi))
	 (cos theta)))))

(defun check-unit-vector (&rest rest)
  ;; turn off for faster speed
  (dolist (e rest)
    (unless (< (abs (- (norm e) 1)) 1d-12)
      (error "vector isn't normalized"))))

(defun check-range (min max &rest rest)
  (dolist (e rest)
    (unless (< min e max)
      (error "range check failed"))))

(defun intersect (start dir c n)
  "intersection of ray and plane"
  (declare (vec start dir c n))
  (let* ((hess-dist (dot c n))
	 (div (dot n dir)))
    (when (< (abs div) 1d-12)
      (error 'ray-and-plane-parallel))
    (let ((eta (/ (- hess-dist (dot n start))
		div)))
      (the vec (.+ start (.s eta dir))))))
#+nil
(intersect (v -1 1) (v 1) (v) (v 1))


(defun refract (start dir c n f na ri)
  (declare (vec start dir c n)
	   (num f na ri))
  (check-unit-vector start dir n)
  (check-range 1 4 ri)
  (check-range 0 100 f)
  (check-range 0 1.6 na)
  (let* ((i (intersect start dir c n)) ;; intersection with lens plane
 	 (rho (.- i c))
	 (cosphi (dot n dir))
	 (r (.- (.s (/ f cosphi) dir)
		rho)) ;; unscaled direction of non-immersed lens
	 (a (.s (* f (- ri 1s0)) n))
	 (ru (.+ r a))
	 (rho2 (dot rho rho))
	 (nf (* ri f))
	 (nf2 (* nf nf))
	 (rat (let ((rat (- nf2 rho2)))
		(when (<= rat 0)
		  (error 'ray-doesnt-hit-principal-sphere))
		rat))
	 (s (.s (- nf (sqrt rat)) dir))
	 (ro (.- ru s))
	 (nro (normalize ro))
	 (cosu (dot nro n)))
    (let ((sinu2 (- 1 (* cosu cosu))))
      (let ((sinu-max (/ na ri)))
	(when (<= (* sinu-max sinu-max)
		  sinu2)
	  (error 'angle-too-steep))
	sinu2))
    (values 
     (the vec ro) ;; direction, not normalized, points to object
     (the vec (.+ s i)) ;; intersection of gaussian sphere and ray
     )))

(defun ray-behind-objective (obj bfp/r c n f na ri)
  "OBJ point in object space (in mm). BFP/R 3D point on BFP, 1 is on
the border of the BFP (values z y x, z is ignored).  Focal length
F, numerical aperture NA and center of objective C."
  (declare (vec obj bfp/r c n)
	   (num f na ri))
  (let* ((theta (find-inverse-ray-angle (norm obj) f))
	 (phi (atan (vy obj) (vx obj)))
	 (dir (v-spherical theta phi))
	 (r (find-bfp-radius na f))
	 (start (.+ c (v (- f) 
			 (* r (vy bfp/r))
			 (* r (vx bfp/r))))))
    (refract start dir c n f na ri)))


(defparameter *centers* ; in mm
  (let* ((l '((6 79 177)
	      (12 125 111)
	      (13 101 135)
	      (13 155 89)
	      (11 219 120)
	      (3 212 168)))
	 (a (make-array 
	     (length l)
	     :element-type 'vec
	     :initial-contents
	     (mapcar #'(lambda (q) 
			 (destructuring-bind (z y x) q
			   (v (* .01 z) (* .001 y) (* .001 x))))
		     l))))
    a))



(defun req (&optional name)
  (error "Required argument ~@[~S~] missing" name))

(defclass model ()
  ((centers :accessor centers :type (simple-array vec 1)
	    :initarg :centers
	    :initform (req 'centers))
   ))

#+nil
(make-instance 'model :centers *centers*)