(provide :rayt)
#.(unless (find-package :rayt)
  (make-package :rayt :use '(:cl)))
(in-package :rayt)
(export '(refract))
#.(unless (find-package :poisson)
  (load "poisson.lisp"))
#.(use-package :base)
#.(use-package :poisson)

(progn
  (declaim (ftype (function (num num num) (values num &optional))
		  quadratic-root-dist))
  (defun quadratic-root-dist (a b c)
    (let ((eps (coerce 1s-7 'num))
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
  (declare (num height focal-length))
  (when (< focal-length height)
    (error 'height-smaller-than-focal-length))
  (let ((rat (/ height focal-length)))
    (declare ((single-float 0s0 1s0) rat))
   (the num (asin rat))))
#+nil
(find-inverse-ray-angle 2.2 2.61)

(defun find-bfp-radius (na f)
  (declare (num na f)) (the num (* f na)))

(defun find-focal-length (mag)
  (declare (num mag))  (the num (/ 164.5 mag)))


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
  (check-unit-vector dir n)
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
     (the vec ro)	;; direction, not normalized, points to object
     (the vec (.+ s i)) ;; intersection of gaussian sphere and ray
     )))

(defun ray-behind-objective (obj bfp/r c n f na ri)
  "OBJ point in object space (in mm). BFP/R 3D point on BFP, 1 is on
the border of the BFP (values z y x, z is ignored).  Focal length F,
numerical aperture NA and center of objective C (position where
Gaussian sphere cuts optic axis). Normal N of objective (points in
direction of excitation light)."
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

#+nil
(with-open-file (s "asy/obj.asy" :direction :output
		   :if-exists :supersede
		   :if-does-not-exist :create)
  (macrolet ((asy (str &rest rest)
	       `(progn
		  (format s ,str ,@rest)
		  (terpri s))))
    (flet ((coord (v)
	     (format nil "(~f,~f,~f)"
		     (vx v) (vy v) (vz v))))
      (let* ((f (find-focal-length 63s0))
	    (na 1.4s0)
	    (ri 1.515s0)
	    (rif (* ri f))
	    (r (find-bfp-radius na f)))
       (asy "import three;")
       (asy "size(1000,1000);")
       (dotimes (i (length *centers*))
	(asy "draw(shift(~a)*scale3(~f)*unitsphere,~a);"
	     (coord (.s ri (aref *centers* i)))
	     (* ri (* .001 8))
	     (if (= i 0)
		 "red"
		 "blue")))
       (asy "draw((0,0,0)--(0,0,1));")
       (asy "draw(shift(~f,~f,~f)*scale3(~f)*unitsquare3);"
	    (- r) (- r)
	    (* -1 (+ rif f))
	    (* 2 r))
       (let ((field .1))
	 (asy "draw(shift(~f,~f,0)*scale3(~f)*unitsquare3);"
	      (- field) (- field)
	      (* 2 field)))
       (asy "draw(shift(0,0,~f)*scale3(~f)*unitcircle3);"
	    (- rif)
	    r)
       (let ((nk 3)
	     (nj 5))
	 (dotimes (j nj)
	  (dotimes (k nk)
	    (let ((bfp/r (v 0 (/ j nj)))
		  (obj (v 0 (* .1 (/ k nk)))))
	      (multiple-value-bind (dir hit)
		  (ray-behind-objective obj
					bfp/r
					(v (- rif))
					(v 1 0 0) 
					f na ri)
		(let ((bfp (.s r bfp/r)))
		  (asy "draw(~a--~a--~a);"
		       (coord (.+ bfp
				  (v (- (+ rif f)))))
		       (coord hit) (coord obj))))))))
       (asy "clip(unitsphere);")))))

(defparameter *centers* ; in mm
  (let* ((l '((3 212 168)
	      (6 79 177)
	      (12 125 111)
	      (13 101 135)
	      (13 155 89)
	      (11 219 120)))
	 (central 0)
	 (offset (- (first (elt l central))))
	 (a (make-array 
	     (length l)
	     :element-type 'vec
	     :initial-contents
	     (mapcar #'(lambda (q) 
			 (destructuring-bind (z y x) q
			   (declare (fixnum z y x))
			   (v (* .01 (+ z offset))
			      (* .001 (- y 150))
			      (* .001 (- x 150)))))
		     l))))
    a))


(defclass model ()
  ((centers :accessor centers :type (simple-array vec 1)
	    :initarg :centers
	    :initform (req 'centers))
   ))

#+nil
(make-instance 'model :centers *centers*)
;; for process take random element
;; replace with the last one and shrink array by one

