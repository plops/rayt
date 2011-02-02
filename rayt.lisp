(declaim (optimize (speed 3) (safety 1) (debug 1)))

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
;; for process take random element
;; replace with the last one and shrink array by one

(defun make-random-queue ()
  (the (array fixnum 1)
   (make-array 1000 :element-type 'fixnum
	       :adjustable t
	       :fill-pointer 0)))

(defun push-random (e queue)
  (declare (fixnum e)
	   ((array fixnum 1) queue))
  (vector-push-extend e queue))

(defun pop-random (queue)
  (declare ((array fixnum 1) queue))
  (let* ((n (length queue))
	 (i (if (= 0 n) 
		(return-from pop-random nil)
		(random n)))
	 (e (aref queue i)))
    (setf (aref queue i) (aref queue (1- n)))
    (decf (fill-pointer queue))
    (the fixnum e)))

#+nil
(let ((m (make-random-queue)))
  (dotimes (i 10)
   (push-random i m))
  (dotimes (i 8)
    (format t "~a~%" (pop-random m)))
  m)

(defclass ranges ()
  ((mi :accessor mi :initarg :mi
       :type (array num 1))
   (ma :accessor ma :initarg :ma 
       :type (array num 1))))

(defun make-ranges (&optional (min 0s0) (max (coerce (* 2 pi) 'num)))
  (let ((mi (make-array 8 :element-type 'num ;; NOTE: this conses a bit 
			:initial-element 0s0 
			:fill-pointer 0
			:adjustable t))
	(ma (make-array 8 :element-type 'num 
			:initial-element 0s0 
			:fill-pointer 0
			:adjustable t)))
    (vector-push-extend min mi)
    (vector-push-extend max ma)
    (make-instance 'ranges :ma ma :mi mi)))

(defmethod print-object ((r ranges) (s stream))
  (format s "#<ranges ~{~{~4,2f ~}~}>" (loop for i below (length (mi r)) collect
				(list (aref (mi r) i) (aref (ma r) i)))))

(defgeneric delete-range (ranges pos))
(defmethod delete-range ((r ranges) pos)
  (declare (fixnum pos))
  (let ((n (length (mi r))))
    (when (< pos (- n 1))
      (loop for i from pos below n do
	   (setf (aref (mi r) i) (aref (mi r) (1+ i)))
	   (setf (aref (ma r) i) (aref (ma r) (1+ i)))))
    (vector-pop (mi r))
    (vector-pop (ma r))
    nil))
#+nil
(let ((r (make-ranges)))
  (subtract r .2 .3)
  (delete-range r 0)
  (delete-range r 0)
  (length (mi r)))
#+NIL
(let ((r (make-ranges)))
  (subtract r 0s0 (coerce (* 2 pi) 'num)))

(defgeneric insert (ranges pos min max))
(defmethod insert ((r ranges) pos min max)
  (declare (num min max)
	   (fixnum pos))
  (let ((n (length (mi r))))
    (cond
      ((< pos n) ;; make space for one and move towards right
       (vector-push-extend (aref (mi r) (1- n)) (mi r))
       (vector-push-extend (aref (ma r) (1- n)) (ma r))
       (loop for i from (- n 2) downto pos do
	    (setf (aref (mi r) (1+ i)) (aref (mi r) i)
		  (aref (ma r) (1+ i)) (aref (ma r) i)))
       (setf (aref (mi r) pos) min
	     (aref (ma r) pos) max))
      ((= pos n)
       (vector-push-extend min (mi r))
       (vector-push-extend max (ma r)))
      (t (error "~a" (list 'trying-to-write-behind-ranges 'pos pos 'n n)))))
  pos)

(defgeneric subtract (ranges min max))
(defmethod subtract ((r ranges) min max)
  (declare (num min max))
  (assert (<= min max))
  (let ((2pi (coerce (* 2 pi) 'num))
	(eps .000001s0))
   (cond ((< 2pi min) 
	  (subtract r (- min 2pi) (- max 2pi)))
	 ((< max 0s0)
	  (subtract r (+ min 2pi) (+ max 2pi)))
	 ((< min 0s0)
	  (subtract r 0s0 max)
	  (subtract r (+ min 2pi) 2pi))
	 ((< 2pi max)
	  (subtract r min 2pi)
	  (subtract r 0s0 (- max 2pi)))
	 ((= 0 (fill-pointer (mi r)))
	  nil)
	 (t 
	  (let ((pos 0))
	    (if (< min (aref (mi r) 0))
		(setf pos -1) ;; left has to go to front
		(let ((lo 0)
		      (mid 0)
		      (hi (fill-pointer (mi r))))
		  ;; binary search for segment whose mi is just smaller than min
		  (loop while (< lo (- hi 1)) do
		       (setf mid (floor (+ lo hi) 2))
		       (if (< (aref (mi r) mid) min)
			   (setf lo mid)
			   (setf hi mid)))
		  (setf pos lo)))
	    (cond 
	      ((= -1 pos) ;; left of all segments
	       (setf pos 0))
	      ((< min (aref (ma r) pos)) ;; new interval starts inside old segment 
	       (let ((mip (aref (mi r) pos))
		     (map (aref (ma r) pos)))
		 (if (< (- min mip) eps) ;; new interval starts at beginning of old segment
		     (if (< max map) ;; is shorter than old segment
			 (setf (aref (mi r) pos) max) ; so let it start at new end
			 (delete-range r pos))
		     ;; there is a gap between new and old segment
		     (progn
		       (setf (aref (ma r) pos) min) ;; max of old segment is beginning of new
		       (when (< max map) ;; new segment within old one
			 (insert r (1+ pos) max map))
		       (incf pos)))))
	      ((and (< pos (1- (fill-pointer (mi r)))) ;; there are more old segments on the right
		    (< (aref (mi r) (1+ pos)) max)) ;; new interval is extending into (or over) next old segment
	       (incf pos)))
	    ;; iterate towards right over old segments
	    (loop while (and (< pos (length (mi r))) ;; the loop decreases the fill-pointer!
			     (< (aref (mi r) pos) max))
	       do
		 (if (< (- (aref (ma r) pos) max) eps)
		     (delete-range r pos) ;; new segment overlaps this old one
		     (setf (aref (mi r) pos) max) ;; new segment ends within this old one
		     ))))))
  nil)

#+nil
(let ((r (make-ranges)))
    (subtract r 1.1s0 1.2s0)
    (subtract r 0s0 2s0)
    r)

(defmacro make-adjustable (type init)
  `(make-array 1000 :element-type ',type :initial-element ,init
	       :adjustable t :fill-pointer 0))

(defun make-random-point ()
  (v 0s0 (- (random 2s0) 1) (- (random 2s0) 1)))

(defparameter *candidates* (make-random-queue))
(defparameter *points* (make-adjustable vec (v)))
(defparameter *grid* (make-array (list 1 1 1) :element-type 'fixnum))
(defparameter *points-per-cell* 9)
(declaim ((array fixnum 1) *candidates)
	 ((array vec 1) *points*)
	 ((simple-array fixnum 3) *grid*)
	 (fixnum *points-per-cell*))

(defun print-grid (&key count)
  (destructuring-bind (h w n) (array-dimensions *grid*)
    (dotimes (j h)
     (dotimes (i w)
       (format t "~1a" (if (= -1 (aref *grid* j i 0))
			   "." (if count
				   (loop for k below n count
					(/= -1 (aref *grid* j i k)))
				   "x"))))
     (terpri))))

#+nil
(print-grid :count t)

(defun get-grid-size (radius)
  (declare (num radius))
  (the fixnum
    (max 2 (ceiling (/ 2s0 (* 4s0 radius))))))

(defun make-grid (radius)
  (let ((grid-size (get-grid-size radius)))
   (setf *grid* (make-array (list grid-size grid-size *points-per-cell*)
			    :element-type 'fixnum
			    :initial-element -1))))

(declaim (inline get-point get-grid-point))
#+sbcl
(declaim (ftype (function (vec) (values fixnum fixnum &optional))
		get-grid-point))
(defun get-grid-point (v)
  (declare (vec v))
  (destructuring-bind (h w n) (array-dimensions *grid*)
    (declare (ignore n)
	     (fixnum h w))
    (let* ((x (vx v))
	   (y (vy v)))
      (declare (num x y))
      (values (floor (* (+ 1s0 y) h) 2)
	      (floor (* (+ 1s0 x) w) 2)))))

(defun get-point (index)
  (aref *points* index))

(defun add-to-grid (v ind)
  (declare (vec v))
  (multiple-value-bind (j i) (get-grid-point v)
    (dotimes (k *points-per-cell*)
      (when (= -1 (aref *grid* j i k))
	(setf (aref *grid* j i k) ind)
	(return))
      (when (= k (1- *points-per-cell*))
	(error "grid overflowed max points per cell"))))
  t)

#+nil
(progn
  (run 0.1)
  (print-grid :count t))

(defun add-point (p)
  (declare (vec p))
  (let ((ind (vector-push-extend p *points*)))
    (add-to-grid p ind)
    (push-random ind *candidates*)))

(defun get-tiled (v)
  (declare (vec v))
  (let* ((x (vx v))
	 (y (vy v))
	 (xx (cond ((< x -1s0) (+ 2s0 x))
		   ((< 1s0 x) (- x 2s0))
		   (t x)))
	 (yy (cond ((< y -1s0) (+ 2s0 y))
		   ((< 1s0 y) (- y 2s0))
		   (t y))))
    (the vec (v 0s0 yy xx))))

(defgeneric find-neighbour-ranges (ranges index radius))
(defmethod find-neighbour-ranges ((r ranges) index radius)
  (declare (fixnum index)
	   (num radius))
  (let* ((candidate (get-point index))
	 (range2 (* 4s0 4s0 radius radius))
	 (grid-size (array-dimension *grid* 1)) 
	 (grid-cell-size (/ 2s0 grid-size))
	 (n (max (floor grid-size 2s0) 
		 (ceiling (* 4s0 radius) grid-cell-size))))
    (multiple-value-bind (gy gx) (get-grid-point candidate)
      (let ((xside (if (< (* .5 grid-cell-size)
		       (- (vx candidate) (- (* gx grid-cell-size) 1)))
		       1 0))
	    (yside (if (< (* .5 grid-cell-size)
			  (- (vy candidate) (- (* gy grid-cell-size) 1)))
		       1 0)))
	(declare (fixnum n gx gy)
		 (num grid-cell-size))
	(loop for j from (- n) upto n do
	    (let ((iy (cond ((= j 0) yside)
			    ((= j 1) 0)
			    (t 1))))
	      (declare (fixnum iy))
	      (loop for i from (- n) upto n do
		   (let* ((ix (cond ((= i 0) xside)
				    ((= i 1) 0)
				    (t 1)))
			  ;; offset to closes cell
			  (dx (- (vx candidate)
				 (* grid-cell-size (+ 0s0 gx ix i))
				 -1s0))
			  (dy (- (vy candidate)
				 (* grid-cell-size (+ 0s0 gy iy j))
				 -1s0)))
		     (declare (fixnum ix iy i j))
		     (when (< (+ (* dx dx) (* dy dy))
			      range2)
		       (let ((cx (floor (mod (+ gx i grid-size) grid-size)))
			     ;; make sure the range of cx is 0 .. grid-size-1
			     (cy (floor (mod (+ gy j grid-size) grid-size))))
			 (dotimes (k *points-per-cell*)
			   (declare ((simple-array fixnum 3) *grid*))
			   (let ((ind (aref *grid* cy cx k)))
			     (declare (fixnum ind))
			     (if (= -1 ind)
				 (return)
				 (if (/= ind index)
				     (let* ((pt (get-point ind))
					    (v (get-tiled (.- pt candidate)))
					    (dist2 (dot v v)))
				       (declare ((single-float 0s0) dist2))
				       (when (< dist2 range2)
					 (let* ((d (sqrt dist2))
						(angle (atan (vy v) (vx v)))
						(arg (* .25s0 (/ d radius)))
						(theta (acos arg)))
					   (declare ((single-float 0s0) d)
						    ((single-float 0s0 1s0) arg)
						    (num angle theta))
					   (subtract r
						     (- angle theta) 
						     (+ angle theta)))))))))))))))))))

(defgeneric make-random-angle (ranges))
(defmethod make-random-angle ((r ranges))
  (let ((n (length (mi r))))
    (declare (fixnum n))
    (when (= n 0)
      (error "range is empty"))
    (let* ((i (random n))
	   (min (aref (mi r) i))
	   (max (aref (ma r) i)))
      (declare (fixnum i)
	       (num min max))
      (the num (+ min (random (- max min)))))))


(defun make-periphery-point (v angle radius)
  (declare (vec v)
	   (num angle radius))
  (let* ((x (* 2 radius (cos angle)))
	 (y (* 2 radius (sin angle))))
    (the vec (get-tiled (.+ v (v 0 y x))))))

(defun generate-poisson (radius)
  (setf
   *candidates* (make-random-queue)
   *points* (make-adjustable vec (v)))
  (make-grid radius)
  (add-point (make-random-point))
  (loop while (< 0 (length *candidates*)) do
       (let* ((index (pop-random *candidates*))
	      (candidate (get-point index))
	      (rl (make-ranges))
	      (pi/3 #.(coerce (/ pi 3) 'num)))
	 (find-neighbour-ranges rl index radius)
	 (loop while (< 0 (length (mi rl)))
	    do
	      (let* ((angle (make-random-angle rl)) 
		     (pt (make-periphery-point
			  candidate angle radius)))
		(declare (num angle))
		(add-point pt)
		(subtract rl (- angle pi/3) (+ angle pi/3)))))))
(require :sb-sprof)
#+nil
(sb-sprof:with-profiling (:max-samples 1000
				       :report :flat
				       :loop nil)
  (generate-poisson .01))

(defun run (&optional (radius .1s0))
  (with-open-file (s "/dev/shm/o.asy" :direction :output
		     :if-exists :supersede
		     :if-does-not-exist :create)
    (macrolet ((asy (str &rest rest)
		 `(progn
		    (format s ,str ,@rest)
		    (terpri s))))
      (asy "import graph;size(400,400);")
      (let* ((r radius)
;	     (gs (get-grid-size r))
;	     (cs (/ 2 gs))
	     )
	(generate-poisson r)
	;; (dotimes (i (length m))
	;;   (let ((e (aref m i)))
	;;     (asy "fill(Circle((~f,~f),~f),gray(.8));"
	;; 	 (vx e) (vy e)
	;; 	 (* 2 r))))
	;; (dotimes (i (length m))
	;;   (let ((e (aref m i)))
	;;     (asy "draw(Circle((~f,~f),~f),gray(.6));"
	;; 	 (vx e) (vy e)
	;; 	 (* 2 r))))
 	;; (loop for i from -1 upto 1 by cs do
	;;      (asy "draw((-1,~f)--(1,~f));" i i))
	;; (loop for i from -1 upto 1 by cs do
	;;      (asy "draw((~f,-1)--(~f,1));" i i))
	
	(dotimes (i (length *points*))
	  (let ((e (aref *points* i)))
	   (asy "dot((~f,~f));" (vx e) (vy e))))
	
	#+nil (let ((n (length m)))
	 (dotimes (i n)
	   (let ((e (aref m i)))
	     (asy "fill(Circle((~f,~f),~f),gray(.4));"
		  (vx e) (vy e) r)
	     (asy "label(\"~d\",(~f,~f),NW);"
		  i (vx e) (vy e)))))))))

#+nil
(time (generate-poisson .01s0))
#+nil
(time
 (progn
   (run .1s0)
   #+nil (print-grid :count t)
   (terpri)))
