(declaim (optimize (speed 0) (safety 3) (debug 3)))

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

(eval-when (:compile-toplevel)
 (defun num (x)
   (declare (number x))
   (coerce x 'single-float)))

(defmacro vec (&rest rest)
  (let ((a (gensym)))
   `(let ((,a (make-array ,(list-length rest)
			 :element-type 'num)))
      ,@(let ((i 0))
	     (declare ((integer 0 #.(1- (expt 2 16))) i))
	     (loop for e in rest collect
		  (prog1
		       `(setf (aref ,a ,i)
			     ,(typecase e ;; FIXME once-only e?
					(fixnum (* 1s0 e))
					(num e)
					(t `(* #.(num 1) ,e))))
		     (incf i))))
      ,a)))
#+nil
(vec 2 2 1)
#+nil
(vec .2 .3 1)
(defun v (&optional
	  (z #.(num 0))
	  (y #.(num 0))
	  (x #.(num 0)))
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
  (let ((r (num 0)))
    (declare (num r))
    (dotimes (i (length a))
      (incf r (* (aref a i) (aref b i))))
    (the num r)))
#+nil
(dot (v 1 3 0) (v 3))

(defun norm (v)
  (declare (vec v))
  (let ((l2 (dot v v)))
    (declare ((single-float 0s0) l2)) ;; FIXME: write num here
    (the num (sqrt l2))))
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
  (declare ((single-float 0s0 #.(/ (coerce pi 'num) 4)) theta)
	   ((single-float 0s0 #.(* 2 (coerce pi 'num))) phi))
  (let* ((st (sin theta)))
    (the vec
      (v (cos theta)
	 (* st (sin phi))
	 (* st (cos phi))))))

(defun check-unit-vector (&rest rest)
  ;; turn off for faster speed
  (dolist (e rest)
    (unless (< (abs (- (norm e) 1)) 1s-7)
      (error "vector isn't normalized"))))

(defun check-range (min max &rest rest)
#+nil  (declare (num min max))
  (dolist (e rest)
    (declare (num e))
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



(defun req (&optional name)
  (error "Required argument ~@[~S~] missing" name))

(defclass model ()
  ((centers :accessor centers :type (simple-array vec 1)
	    :initarg :centers
	    :initform (req 'centers))
   ))

#+nil
(make-instance 'model :centers *centers*)

(defun aref-circ (a j i)
  (declare ((simple-array fixnum 2) a)
	   (fixnum j i))
  (destructuring-bind (y x) (array-dimensions a)
    (let ((ii (cond ((< i 0) (+ i x))
		    ((<= x i) (- i x))
		    (t i)))
	  (jj (cond ((< j 0) (+ j y))
		    ((<= y j) (- j y))
		    (t j))))
      (aref a jj ii))))


(defun in-rectangle-p (p pmin pmax)
  (declare (vec p pmin pmax))
  (the boolean
   (progn
     (loop for i from 1 below 3 do
	  (unless (< (aref pmin i) (aref p i) (aref pmax i))
	    (return-from in-rectangle-p nil)))
     t)))

#+nil
(in-rectangle-p (v 0 7.4 9.2) (v) (v 0 12 12))

(defun generate-poisson (h w min-dist &optional (new-point-count 30))
  (let* ((cell-size (/ min-dist (sqrt 2)))
	(grid (make-array (list (ceiling h cell-size)
				(ceiling w cell-size))
			  :element-type 'fixnum
			  :initial-element -1))
	(points (make-array 10000
			    :element-type 'vec
			    :initial-element (v)
			    :adjustable t
			    :fill-pointer 0))
	(process ())
	(output (make-array 1000
			    :element-type 'fixnum
			    :initial-element -1
			    :adjustable t
			    :fill-pointer 0)))
    (labels ((image->grid (p)
	       (declare (vec p))
	       (the vec 
		 (v 0
		    (floor (vy p) cell-size)
		    (floor (vx p) cell-size))))
	     (make-random-point (p &optional (mindist nil))
	       (declare (vec p)
			((or null num) mindist))
	       (let ((radius (+ (* 2 min-dist) 
				(if mindist
				    (random mindist)
				    0s0)))
		     (angle (random (coerce (* 2 pi) 'num))))
		 (the vec
		   (v 0
		      (+ (vy p) (* radius (sin angle)))
		      (+ (vx p) (* radius (cos angle)))))))
	     (point-collides-p (p)
	       (declare (vec p))
	       (the boolean
		 (let ((gp (image->grid p)))
		   (loop for j from -2 upto 2 do
			(loop for i from -2 upto 2 do
			     (let* ((x (floor (vx gp)))
				    (y (floor (vy gp)))
				    (ind (aref-circ grid (+ y j) (+ x i))))
			       (unless (= -1 ind)
				 (format t "found circ ~a~%"
					 (list ind x y process))
				 (when (< min-dist (norm (.- p (aref points ind))))
				   (return-from point-collides-p t))))))))))
      (let ((first-point (v 0 
			    (random (* 1s0 h)) 
			    (random (* 1s0 w))))
	    (ind 0)) 
	;; store the first point in all containers
	(vector-push-extend first-point points)
	(vector-push-extend ind output)
	(push ind process)
	(let ((q (image->grid first-point)))
	  (format t "first point ~a~%" (list first-point q ind))
	  (setf (aref grid (floor (vy q)) (floor (vx q))) 
	       ind)))
      (loop while process do
	   (let* ((np (length process))
		  (point (elt process (if (= 0 np) 
					 0
					 (random np)))))
	     (setf process (remove point process))
	     (dotimes (i new-point-count)
	       (let ((new-point (make-random-point (aref points point))))
		 (when (and (in-rectangle-p new-point (v) (v 0 h w))
			    (not (point-collides-p new-point)))
		   (vector-push-extend new-point points)
		   (let ((cur (fill-pointer points)))
		     (push cur process)
		     (vector-push-extend cur output)
		     (let ((q (image->grid new-point)))
		       (setf (aref grid (floor (vy q))
				   (floor (vx q)))
			     cur))))))))
      (let ((m (make-array (length output)
			   :element-type 'vec
			   :initial-element (v))))
	(dotimes (i (length m))
	  (setf (aref m i) (aref points (aref output i))))
	m))))


(let ((m (make-array 10 :element-type 'single-float
		     :fill-pointer 0 :adjustable t)))
  (vector-push-extend 3s0 m)
  (vector-push-extend 33s0 m)
  (list (fill-pointer m) m))

(progn
 (defun run ()
   (with-open-file (s "/dev/shm/o.asy" :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
     (macrolet ((asy (str &rest rest)
		  `(progn
		     (format s ,str ,@rest)
		     (terpri s))))
       (asy "import graph;size(400,400);")
       (let* ((r 10s0)
	      (m (generate-poisson 400 400 r 10)))
	 (dotimes (i (length m))
	   (asy "draw(Circle((~f,~f),~f));"
		(vx (aref m i)) (vy (aref m i))
		r))))))
 (run))
