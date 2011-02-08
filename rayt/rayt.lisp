(in-package :rayt)

(declaim (optimize (speed 1) (safety 3) (debug 3)))


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
	   (offset-mm (v (vz (aref *centers* nucleus)) 0 0)))
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
  `(concatenate 'string "/home/martin/tmp/t0207/"
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


(defun make-image (h &key (w h) (type '(unsigned-byte 8)))
  (make-array (list h w)
	      :element-type type))




(defun sum-bfp (nucleus &key (w-bfp 100) (w-ffp 30) (radius-ffp-mm 1.5s-3)
		(radius-bfp-mm radius-ffp-mm))
  "Average BFP exposure pattern for an area with radius RADIUS-FFP-MM
on LCoS that illuminates the given NUCLEUS. The size of the projected
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


 #+nil
 (format t "~a~%"
	 (run))



 (defun run ()
   (let ((bfp (make-array (list 200 200) :element-type '(unsigned-byte 8)))
	 (res ())
	 (nuc 0))
     (dotimes (i (length *centers*))
       (unless (= i nuc)
	 (sum-bfp-raster bfp nuc i :radius-ffp-mm 16s-3) 
	 #+nil(push (project-nucleus-into-bfp bfp nuc i (v) :radius-mm 16s-3 :triangles 23)
		    res)))
     (write-pgm "/dev/shm/bfp.pgm" bfp)
     ;; compare with raytracer
     (write-pgm "/dev/shm/bfp-rt.pgm" (normalize-im (trace-from-bfp
						     (v) nuc :w 200 :radius 16s-3)))
     (reverse res)))


 #+nil
 (require :rayt)

 #+nil
 (let* ((ri 1.515s0)
	(f (find-focal-length 63s0))
	(dir (ray-behind-objective (v 0 0 0) 
				   (v 0 0 .999) 
				   (v (* -1 ri f))
				   (v 1) f 1.45 ri)))
   (project-ray-into-bfp dir ri f))

 #+nil
 (project-ray-into-bfp (normalize (v (sqrt 2) 1 0)) 1.515 (find-focal-length 63s0)))