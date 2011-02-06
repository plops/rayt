(in-package :rayt)

(declaim (optimize (speed 1) (safety 2) (debug 2)))


(progn
  (declaim (ftype (function (num num num) (values num &optional))
		  quadratic-root-dist))
  (defun quadratic-root-dist (a b c)
    (let ((eps (coerce 1s-7 'num))
	  (zero 0s0)
	  (det2 (- (* b b) (* 4 a c))))
     (if (<= det2 zero)
       zero
       (let* ((q (* .5s0 (+ b (* (signum b) (sqrt det2)))))
	      (aa (abs a))
	      (aq (abs q)))
	 (cond ((or (< aq eps) (< aa eps)) zero)
	       (t (- (/ q a) (/ c q)))))))))
#+nil
(quadratic-root-dist 1s0 2s0 -3s0)

(declaim (ftype (function (vec vec vec num)
			  #+ecl num #-ecl (values num &optional))
		ray-sphere-intersection-length)
	 (inline ray-sphere-intersection-length))
(defun ray-sphere-intersection-length (start dir center radius)
  (declare (type vec start dir center)
	   (type num radius))
 ; (check-unit-vector dir)
  (let* ((l (.- center start))
	 (c (- (dot l l) (* radius radius)))
	 (b (* -2s0 (dot l dir))))
    (the num (abs (quadratic-root-dist 1s0 b c)))))

(defun find-inverse-ray-angle (height focal-length)
  "HEIGHT meridional distance of a point in sample space towards optic
axis. Coordinates in mm."
  (declare (type num height focal-length))
  (when (< focal-length height)
    (error 'height-smaller-than-focal-length))
  (let ((rat (/ height focal-length)))
    (declare (type (single-float 0s0 1s0) rat))
   (the num (asin rat))))
#+nil
(find-inverse-ray-angle 2.2 2.61)

(defun find-bfp-radius (na f)
  (declare (type num na f)) (the num (* f na)))

(defun find-focal-length (mag)
  (declare (type num mag))  (the num (/ 164.5 mag)))


(defun intersect (start dir c n)
  "intersection of ray and plane"
  (declare (type vec start dir c n))
  (let* ((hess-dist (dot c n))
	 (div (dot n dir)))
    (when (< (abs div) 1d-12)
      (error 'ray-and-plane-parallel))
    (let ((eta (/ (- hess-dist (dot n start))
		div)))
      (the vec (.+ start (.s eta dir))))))
#+nil
(intersect (v -1 1) (v 1) (v) (v 1))

(define-condition ray-does-not-hit-principal-sphere () ())
(define-condition ray-too-steep () ())

(declaim (ftype (function (vec vec vec vec num num num)
			  (values vec vec &optional))
		refract))
(defun refract (start dir c n f na ri)
  (declare (type vec start dir c n)
	   (type num f na ri))
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
		  (error 'ray-does-not-hit-principal-sphere))
		rat))
	 (s (.s (- nf (sqrt rat)) dir))
	 (ro (.- ru s))
	 (nro (normalize ro))
	 (cosu (dot nro n)))
    (let ((sinu2 (- 1 (* cosu cosu))))
      (let ((sinu-max (/ na ri)))
	(when (<= (* sinu-max sinu-max)
		  sinu2)
	  (error 'ray-too-steep))
	sinu2))
    (values 
     (the vec nro)	;; direction, points to object, either nro or ro
     (the vec (.+ s i)) ;; intersection of gaussian sphere and ray
     )))

(defun ray-behind-objective (obj bfp/r c n f na ri)
  "OBJ point in object space (in mm). BFP/R 3D point on BFP, 1 is on
the border of the BFP (values z y x, z is ignored).  Focal length F,
numerical aperture NA and center of objective C (position where
Gaussian sphere cuts optic axis). Normal N of objective (points in
direction of excitation light)."
  (declare (type vec obj bfp/r c n)
	   (type num f na ri))
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
       (let ((pd (generate-poisson .04s0)))
	 (dotimes (j (length pd))
	   (let ((bfp/r (aref pd j))
		 (obj (v 0 .01 .02)))
	     (handler-case 
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
			  (coord hit) (coord obj))))
	       (ray-too-steep ())
	       (ray-does-not-hit-principal-sphere ())))))
       (asy "clip(unitsphere);")))))

(defparameter *centers* nil) ;; in mm
(defparameter *centers-fix* nil)
(defparameter *data-dz-mm* 1s-3)
(defparameter *data-dx-mm* .1s-3)
(defparameter *data-width-px* 300)
			
(let* ((l '((3 212 168)
	    (6 79 177)
	    (12 125 111)
	    (13 101 135)
	    (13 155 89)
	    (11 219 120)))
      ; (central 0)
       (offset 0 #+nil (- (first (elt l central))))
       (a (make-array 
	   (length l)
	   :element-type 'vec
	   :initial-contents
	   (mapcar #'(lambda (q) 
		       (destructuring-bind (z y x) q
			 (declare (type fixnum z y x))
			 (v (* .01 (+ z offset))
			    (* .001 (- y 150))
			    (* .001 (- x 150)))))
		   l)))
       (afix (make-array (length l)
			 :element-type 'vec
			 :initial-contents
			 (mapcar #'(lambda (q) 
				     (destructuring-bind (z y x) q
				       (declare (type fixnum z y x))
				       (v z y x)))
				 l))))
  (setf *centers* a
	*centers-fix* afix))


(defclass model ()
  ((centers :accessor centers :type (simple-array vec 1)
	    :initarg :centers
	    :initform (req 'centers))
   ))

#+nil
(make-instance 'model :centers *centers*)

(declaim (ftype (function ((simple-array num 2) fixnum vec num)
			  (values num &optional))
		incf-vec)
	 (inline incf-vec))
(defun incf-vec (img n v val)
  "Increase the value in an image at a position given by V."
  (declare (type (simple-array num 2) img)
	   (type fixnum n) (type vec v)
	   (type num val))
  (let* ((x (vx v))
	 (y (vy v))
	 (i (floor (* (+ 1s0 x) n) 2))
	 (j (floor (* (+ 1s0 y) n) 2)))
    (declare (type num x y)
	     (type fixnum i j))
    (incf (aref img j i) val)))

(defun write-pgm (filename img)
  (declare (type simple-string filename)
	   (type (array (unsigned-byte 8) 2) img))
  (destructuring-bind (h w) (array-dimensions img)
    (declare (type (integer 0 65535) w h))
    (with-open-file (s filename
		       :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
      (declare (type stream s))
      (format s "P5~%~D ~D~%255~%" w h))
    (with-open-file (s filename 
		       :element-type '(unsigned-byte 8)
		       :direction :output
		       :if-exists :append)
      (let ((data-1d (make-array (* h w)
		      :element-type '(unsigned-byte 8)
		      :displaced-to img)))
	(write-sequence data-1d s)))
    nil))

(defun linear-array (img)
  (let ((img1 (make-array (reduce #'* (array-dimensions img))
			  :element-type (array-element-type img)
			  :displaced-to img)))
    img1))

(defvar *pd* nil)

(defun normalize-im (img)
  (declare (type (simple-array num 2) img))
  (destructuring-bind (h w) (array-dimensions img)
   (let* ((img1 (linear-array img))
	  (buf (make-array (list h w) 
			   :element-type '(unsigned-byte 8)
			   :adjustable nil))
	  (buf1 (linear-array buf))
	  (ma (let ((e (reduce #'max img1)))
		(unless e
		  (return-from normalize-im buf))
		e))
	  (s (/ 255s0 ma)))
     (declare (type (simple-array (unsigned-byte 8) 2) buf)
	      (type (array (unsigned-byte 8) 1) buf1))
     (dotimes (i (length buf1))
       (setf (aref buf1 i) (floor (* s (aref img1 i)))))
     (the (simple-array (unsigned-byte 8) 2) buf))))

(defun trace-from-bfp (obj nucleus &key (radius 2s-3) (w 100))
  (declare (type vec obj)
	   (type fixnum w)
	   (type num radius))
  (unless *pd*
    (setf *pd* (generate-poisson .01s0)))
  (let* ((f (find-focal-length 63s0))
	(na 1.4s0)
	(ri 1.515s0)
	(rif (* ri f))
;	(r (find-bfp-radius na f))
	(pd *pd*)
	(img (make-array (list w w) :element-type 'num
			 :adjustable nil
			 :initial-element 0s0)))
   (declare (type (simple-array vec 1) pd))
   (dotimes (j (length pd))
     (let ((bfp/r (aref pd j))
	   (offset-mm (v 0 0 (vz (aref *centers* nucleus)))))
       (handler-case 
	   (multiple-value-bind (dir hit)
	       (ray-behind-objective obj bfp/r
				     (v (- rif))
				     (v 1 0 0) 
				     f na ri)
	     (declare (ignore hit))
	     (loop for i below (length *centers*) do
		  (unless (= i nucleus)
		   (incf-vec img w bfp/r 
			     (ray-sphere-intersection-length 
			      obj dir
			      (.- (aref *centers* i) offset-mm)
			      radius)))))
	 (ray-too-steep ())
	 (ray-does-not-hit-principal-sphere ()))))
   (the (simple-array num 2) img)))


(defmacro tmp (str &rest rest)
  `(concatenate 'string "/home/martin/tmp/t0204/"
		(format nil ,str ,@rest)))


(defun ..+ (dst src)
  (let ((a1 (linear-array dst))
	(b1 (linear-array src)))
   (dotimes (i (length a1))
     (incf (aref a1 i) (aref b1 i))))
  dst)

#+nil
(progn
  (setf *pd* (generate-poisson .003s0))
  (length *pd*))

#+nil
(time
 (let* ((n 8)
	(field .1s0)			; diameter
	(s (/ field n))
	(w 200)
	(sum (make-array (list w w) :element-type 'num
			 :initial-element 0s0)))
   (dotimes (j n)
     (dotimes (i n)
        
       (let ((bfp (trace-from-bfp 
		   (.- (.s s (v 0 j i))
		       (.s (* .5s0 field) (v 0 1 1)))
		   0
		   :radius 32s-3
		   :w w)))
	 (..+ sum bfp)
	 (write-pgm (tmp "bfp~2,'0d-~2,'0d.pgm" j i)
		    (normalize-im
		     bfp)))))
   (write-pgm (tmp "sum.pgm")
	      (normalize-im sum))))

(defun get-nuclei-in-slice (slice-px r-px)
  "Return a list of nuclei that should be visible in the given slice."
  (declare (type fixnum slice-px r-px))
  (let ((res ()))
   (dotimes (nucleus (length *centers-fix*))
     (let ((c (aref *centers-fix* nucleus)))
       (when (<= (abs (- slice-px (vz c))) r-px)
	 (push nucleus res))))
   res))

(defun make-image (h &key (w h) (type '(unsigned-byte 8)))
  (make-array (list h w)
	      :element-type type))


(defun draw-nucleus (img nucleus &key 
		     (radius-mm 1.5s-3)
		     (dx-mm 1s0))
  (declare (type (simple-array (unsigned-byte 8) 2) img)
	   (type fixnum nucleus)
	   (type num radius-mm dx-mm))
  (destructuring-bind (h w) (array-dimensions img)
    (declare (ignore h))
    (let ((s (/ (* 1s0 w) *data-width-px*)))
      (setf dx-mm (* s *data-dx-mm*))
      (let* ((c (aref *centers-fix* nucleus))
	     ;(s (/ dx-mm *data-dx-mm*))
	     (q (.s s c)) ;; note: vz wrong
	     (x (floor (vx q)))
	     (y (floor (vy q)))
	     (rad (floor (* s radius-mm) *data-dx-mm*)))
       (format t "~a~%" (list dx-mm 'dx-mm c q rad 's s ))
       (raster-disk img y x rad))))
  img)

(defun draw-slice-through-point (slice-px &key 
				 (radius-mm .003s0)
				 (h 100) (w h)
				 (dx-mm (* *data-dx-mm* (/ *data-width-px* w))))
  (declare (type fixnum slice-px w h)
	   (type num radius-mm dx-mm))
  (let ((img (make-image h :w w)))
    (dolist (e (get-nuclei-in-slice slice-px
				    (ceiling radius-mm *data-dz-mm*)))
      (draw-nucleus img e :radius-mm radius-mm :dx-mm dx-mm))
    img))

#+nil
(dotimes (k 20)
 (write-pgm (tmp "lcos~3,'0d.pgm" k)
	    (draw-slice-through-point k :radius-mm 2.5e-3 :h 300)))

(defun sum-bfp (nucleus &key (w-bfp 100) (w-ffp 30) (radius-ffp-mm 1.5s-3)
		(radius-bfp-mm radius-ffp-mm))
  "Average BFP exposure pattern for an area with radius RADIUS-FFP-MM
on LCoS that illuminates the give NUCLEUS. The size of the projected
spheres is defined by RADIUS-BFP-MM."
  (let* ((ffp (make-image w-ffp))
	 (bfp (make-image w-bfp :type 'num))
	 (field (* *data-dx-mm* *data-width-px*))
	 (s (/ field w-bfp)))
    (draw-nucleus ffp nucleus :radius-mm radius-ffp-mm)
    (dotimes (j w-ffp)
      (dotimes (i w-ffp)
	(when (< 0 (aref ffp j i))
	  (..+ bfp
	       (trace-from-bfp 
		(.- (.s s (v 0 j i))
		    (.s (* .5s0 field) (v 0 1 1)))
		nucleus
		:radius radius-bfp-mm
		:w w-bfp)))))
    (write-pgm (tmp "sumffp~3,'0d.pgm" nucleus) ffp)
    (write-pgm (tmp "sumbfp~3,'0d.pgm" nucleus) (normalize-im bfp))))

#+nil
(dotimes (k (length *centers*))
  (time (sum-bfp k 
	    :radius-ffp-mm 2.4s-3
	    :radius-bfp-mm 10s-3
	    :w-bfp 300)))



;; direct (and therefore faster) projection of the sphere through
;; objective onto bfp

(defun project-nucleus-into-bfp (bfp nucleus ffp-pos 
				 &key (mag 63s0) (f (find-focal-length mag))
				 (ri 1.515s0) (radius-mm 1.5s-3) (triangles 13) (na 1.45))
  (declare (type fixnum nucleus triangles)
	   (type (simple-array (unsigned-byte 8) 2) bfp)
	   (type vec ffp-pos)
	   (type num radius-mm ri mag f))
  (let* ((radius (* radius-mm ri))
	 (nuc-center (.s ri (.- (aref *centers* nucleus) ffp-pos)))
	 (dist (norm nuc-center))
	 (c0 (.s (/ dist) nuc-center))
	 ;; the billboard bounded by the tangents from ffp-pos to
	 ;; nucleus.  2D construction of diameter: find difference
	 ;; between the two intersections of a circle with radius
	 ;; DIST=R and the nucleus' circle RADIUS-MM=r: (\vec x-(R
	 ;; 0)^T)^2=R^2 and (\vec x)^2 = r^2, \vec x=(x,y)^T,
	 ;; x..distance from nucleus center to bill-board. y..radius of
	 ;; billboard
	 (bb-x (/ (* radius radius)
		  (* 2 dist)))
	 (bb-y (* radius
		  (let ((arg (- 1s0 (/ (* 2s0 dist)))))
		    (if (<= 0 arg)
			(sqrt arg)
			(error "distance from nucleus center to billboard is negative~% ~a" 
			       (list 'dist dist 'arg arg 'bb-x bb-x))))))
	 (meridian (.s bb-y (cross c0 (v 0 0 1)))) ;; prependicular to ray
	 (rif (* ri f))
	 (bfp-rad (find-bfp-radius na f))
	 (s (/ (array-dimension bfp 0) bfp-rad))
	 (vertices (loop for i below triangles collect
			(let* ((m (rotation-matrix 
				   (* i (/ +2*pi+ triangles)) c0))
			       (cc (.+ (.s (- 1 (/ bb-x dist)) nuc-center) 
				       (m* m meridian)))
			       (gauss-hit (.s (* s rif) (normalize cc))))
			  gauss-hit))))
    (dotimes (i triangles)
      (let ((e (elt vertices i))
	    (f (elt vertices (if (< i (1- triangles))
				 (1+ i)
				 0))))
	(macrolet ((tri (a b c)
		     `(raster-triangle bfp ,@(let ((r ()))
						  (dolist (e (list a b c))
						       (push `(vy ,e) r)
						       (push `(vx ,e) r))
						  (reverse r)) 1)))
	  (tri c0 e f))))))


(defun run ()
 (let ((bfp (make-array (list 200 200) :element-type '(unsigned-byte 8))))
   (project-nucleus-into-bfp bfp 0 (v))))

#+nil
(run)