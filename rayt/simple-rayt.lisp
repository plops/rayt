(in-package :simple-rayt)

(defun project-ray-into-bfp (dir ri f)
  "Given a direction in sample space, return the point where ray will
  intersect the BFP."
  (declare (type vec dir)
	   (num ri f))
  (check-unit-vector dir)
  (let* ((x (vx dir))
	 (y (vy dir))
	 (z (vz dir))
	 ;(q (/ y x))
	 (atanq (atan y x))
	 (cosatanq (cos atanq) #+nil (/ (sqrt (+ 1s0 (* q q)))))
	 (sinatanq (sin atanq) #+nil (* x cosatanq))
	 (sinacosz (sin (acos z)) #+nil (sqrt (- 1s0 (* z z))))
	 (nf (* ri f))
	 (rho (* nf sinacosz)))
    (format t "project into bfp ~%~a~%" (list 
					    'cosatanq cosatanq
					    'sinatanq sinatanq
					    'sinacosz sinacosz
					    'nf nf
					    'rho rho
					    'y (* rho sinatanq)
					    'x (* rho cosatanq)
					    ))
    (the vec
     (v (- (+ nf f))
	(* rho sinatanq)
	(* rho cosatanq)))))

#+nil
(* 1.515 (find-focal-length 63s0)
   (sin (atan 1 1)))

#+nil
(project-ray-into-bfp (normalize (v 1 1 0)) 1.515 (find-focal-length 63s0))

;; direct (and therefore faster) projection of the sphere through
;; objective onto bfp


(defun project-nucleus-into-bfp (bfp illum-nucleus protected-nucleus ffp-pos 
				 &key (mag 63s0) (f (find-focal-length mag))
				 (ri 1.515s0) (radius-mm 1.5s-3) (triangles 13) (na 1.45))
  (declare (type fixnum illum-nucleus protected-nucleus triangles)
	   (type (simple-array (unsigned-byte 8) 2) bfp)
	   (type vec ffp-pos)
	   (type num radius-mm ri mag f))
  (with-open-file (stream (tmp "obj.asy") :direction :output
		   :if-exists :supersede
		   :if-does-not-exist :create)
  (macrolet ((asy (str &rest rest)
	       `(progn
		  (format stream ,str ,@rest)
		  (terpri stream))))
    (flet ((coord (v)
	     (format nil "(~f,~f,~f)"
		     (vx v) (vy v) (vz v))))
      (let* ((f (find-focal-length 63s0))
	    (na 1.4s0)
	    (ri 1.515s0)
	    (rif (* ri f))
	    (r (find-bfp-radius na f)))
       (asy "import three;")
       (asy "size(200,200);")
       (let* ((radius (* radius-mm ri))
	      (offset-mm (v (vz (aref *centers* illum-nucleus)) 0 0)) ;; focus on nucleus
	      (nuc-center (let ((c (.- (aref *centers* protected-nucleus)
				       (.+ ffp-pos offset-mm))))
			    (.s (* (signum (- (vz c))) ri) c)))
	 (dist (norm nuc-center))
	 (c0 (.s (/ dist) nuc-center))
	 ;; the billboard bounded by the tangents from ffp-pos to
	 ;; nucleus.  2D construction of diameter: find difference
	 ;; between the two intersections of a circle with radius
	 ;; DIST=R and the nucleus' circle RADIUS-MM=r: (\vec x-(R
	 ;; 0)^T)^2=R^2 and (\vec x)^2 = r^2, \vec x=(x,y)^T,
	 ;; x..distance from nucleus center to bill-board. y..radius of
	 ;; billboard
	 ;; solve([(x-R)^2+y^2=R^2,x^2+y^2=r^2],[x,y]);
	      (bb-x (/ (* radius radius)
		       (* 2 dist)))
	      (bb-y (let ((y2 (- (* radius radius) (* bb-x bb-x)))) 
		      (if (<= 0 y2) (sqrt y2)
			  (error "arg neg ~a" (list 'bb-x bb-x)))))
	      (meridian (.s bb-y (cross c0 (v 0 0 1)))) ;; prependicular to ray
	      (rif (* ri f))
	      (bfp-rad (find-bfp-radius na f))
	      (s (/ (array-dimension bfp 0) (* 2 bfp-rad))))
	 (flet ((scale (q)
		  (.s s (.+ q (v 0 bfp-rad bfp-rad)))))
	   (let* ((ccs ())
		  (vc0 (scale (project-ray-into-bfp c0 ri f)))
		 (vertices (loop for i below triangles collect
				(let* ((m (rotation-matrix 
					   (* i (/ +2*pi+ triangles)) c0))
				       (cc (.+ (.s (- 1 (/ bb-x dist)) nuc-center) 
					       (m* m meridian)))
				       (b (project-ray-into-bfp (normalize cc) ri f)))
				  (push cc ccs) 
				  (scale b)))))
	    (progn	 
	      #+nil (asy "currentprojection=perspective(
camera=(4.98788094302443,4.00955657639435,2.01112380474988),
up=(0.000392603828345185,0.000309874871399901,-0.00154884955268845),
target=(0.0409059444797579,0.105001572683979,-0.0240156297638063),
zoom=0.696946091482217,
angle=3.99957957637564,
autoadjust=false);")
	      (asy "draw(Label(\"$z$\",1),(0,0,0)--(0,0,.1),Arrow3);")
	      (asy "draw(Label(\"$y$\",1,NW),(0,0,0)--(0,.1,0),Arrow3);")
	      (asy "draw(Label(\"$x$\",1),(0,0,0)--(.1,0,0),Arrow3);")
	      (asy "draw((0,0,0)..(0,.1,0));")
	      (asy "draw((0,0,0)..(.1,0,0));")
	      (dotimes (i (length *centers*))
			   (let ((c (let ((c (.- (aref *centers* i)
						 (.+ ffp-pos offset-mm))))
				      (.s (* (signum (- (vz c))) ri) c))))
			     (asy "draw(shift(~a)*scale3(~f)*unitsphere,lightgreen);"
				  (coord c)
				  (* .9 radius))
			     (asy "label(~s,~a);" (format nil "~d" i)
				  (coord (.- c (v 0 (- radius) (- radius)))))))
			 (asy "draw((0,0,0)..~a);" (coord nuc-center))
			 
			 (asy "draw((0,0,0)..~a);" (coord (.+ (.s (- 1 (/ bb-x dist))
								  nuc-center) 
							      meridian)))
			 (dolist (e ccs)
			   (asy "draw((0,0,0)--~a);" (coord e)))
		       
#+nil
			 (progn (asy "draw(~a" (coord vc0))
			   (dolist (e vertices)
			     (asy "--~a" (coord (v 0 (vy e) (vx e)))))
			   (asy ");")
			   )			 
			 (format t "~a~%" (list ccs (map-into vertices 
						     #'(lambda (v) (v (floor (vz v))
								 (floor (vy v)) 
								 (floor (vx v)))) vertices))))
	    (dotimes (i triangles)
	      (let ((e (elt vertices i))
		    (f (elt vertices (if (< i (1- triangles))
					 (1+ i)
					 0))))
		(macrolet ((tri (a b c)
			     `(raster-triangle bfp ,@(let ((r ()))
							  (dolist (e (list a b c))
							    (push `(floor (vy ,e)) r)
							    (push `(floor (vx ,e)) r))
							  (reverse r)) 1)))
		  (tri vc0 e f))))
	    vc0))))))))


#+nil ;asy -tex pdflatex obj.asy
(let ((bfp (make-image 256)))
 (project-nucleus-into-bfp bfp 0 2 (let ((e (aref *centers* 0)))
				     (v 0 (vy e) (vx e)))
			   :radius-mm 16s-3)


 (defun sum-bfp-raster (bfp illum-nucleus protected-nucleus 
			&key (w-ffp 30) (radius-ffp-mm 1.5s-3)
			(radius-project-mm 1.5s-3))
   "Average BFP exposure pattern of PROTECTED-NUCLEUS for an area with
radius RADIUS-FFP-MM on LCoS that illuminates the given
ILLUM-NUCLEUS. The size of the projected spheres is defined by
RADIUS-BFP-MM."
   (declare (type fixnum illum-nucleus protected-nucleus)
	    (type (simple-array (unsigned-byte 8) 2) bfp)
	    (type num radius-ffp-mm))
   (destructuring-bind (h w) (array-dimensions bfp)
     (declare (ignore h))
     (let* ((ffp (make-image w-ffp))
	    (field (* *data-dx-mm* *data-width-px*))
	    (s (/ field w)))
       (draw-nucleus ffp illum-nucleus :radius-mm radius-ffp-mm)
       (dotimes (j w-ffp)
	 (dotimes (i w-ffp)
	   (when (< 0 (aref ffp j i))
	     (project-nucleus-into-bfp bfp illum-nucleus protected-nucleus
				       (.- (.s s (v 0 j i))
					   (.s (* .5s0 field) (v 0 1 1)))
				       :radius-mm radius-project-mm :triangles 23))))
       (format t "~a~%" ffp)
       (write-pgm (tmp "sumffp~3,'0d.pgm" illum-nucleus) ffp)
       (write-pgm (tmp "sumbfp~3,'0d.pgm" illum-nucleus) bfp))))
