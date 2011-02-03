(provide :poisson)
#+ecl (make-package :poisson :use '(:cl))
#.(unless (find-package :poisson)
  (make-package :poisson :use '(:cl)))
(in-package :poisson)
(export '(generate-poisson))
#.(unless (find-package :base)
    (load "base.lisp"))
#+ecl (load "base")
#.(use-package :base)

(declaim (optimize (speed 3) (safety 2) (debug 2)))
#+nil (declaim (optimize (speed 1) (safety 3) (debug 3)))


(defclass ranges ()
  ((mi :accessor mi :initarg :mi
       :type (simple-array num 1))
   (ma :accessor ma :initarg :ma 
       :type (simple-array num 1))
   (n :accessor n :initarg :n :type fixnum :initform 0)))

(defun make-ranges (&optional (min 0s0) (max +2pif+))
  (let ((mi (make-array 8 :element-type 'num ;; NOTE: this conses a bit 
			:initial-element 0s0 
			;:adjustable t
			;:fill-pointer 0
			))
	(ma (make-array 8 :element-type 'num 
			:initial-element 0s0 
			;:adjustable t
			;:fill-pointer 0
			)))
   #+nil (vector-push-extend min mi)
   #+nil (vector-push-extend max ma)
   (setf (aref mi 0) min
	 (aref ma 0) max)
   (make-instance 'ranges :ma ma :mi mi :n 1)))

#+nil
(make-ranges)

(defgeneric reset (ranges &optional min max))
(defmethod reset ((r ranges) &optional (min 0s0) (max +2pif+))
  (with-slots (mi ma n) r
    (declare (fixnum n))
    #+nil (setf (fill-pointer mi) 0
	  (fill-pointer ma) 0)
    (setf n 1)
    (setf (aref mi 0) min
	  (aref ma 0) max)
    #+nil (vector-push-extend min mi)
    #+nil (vector-push-extend max ma)))

#+nil
(let ((r (make-ranges)))
  (subtract r .1 .2)
  (subtract r .4 .5)
  r)

(defmethod print-object ((r ranges) s)
  (format s "#<ranges ~{~{~4,2f ~}~}>"
	  (loop for i below (n r) collect
	       (list (aref (mi r) i) (aref (ma r) i)))))

(defgeneric delete-range (ranges pos))
(defmethod delete-range ((r ranges) pos)
  (declare (fixnum pos))
  (with-slots (mi ma n) r
    (declare (fixnum n))
    (when (< pos (- n 1))
      (loop for i from pos below n do
	   (setf (aref mi i) (aref mi (1+ i)))
	   (setf (aref ma i) (aref ma (1+ i)))))
    (decf n)
    nil))
#+nil
(let ((r (make-ranges)))
  (subtract r .2 .3)
  (delete-range r 0)
 ; (delete-range r 0)
  (subseq (mi r) 0 (n r)))
#+NIL
(let ((r (make-ranges)))
  (subtract r 0s0 (coerce (* 2 pi) 'num)))

(defgeneric insert (ranges pos min max))
(defmethod insert ((r ranges) pos min max)
  (declare (num min max)
	   (fixnum pos))
  (with-slots (mi ma n) r
    (declare (fixnum n))
    (cond
      ((< pos n) ;; make space for one and move towards right
       (loop for i from (- n 1) downto pos do
	    (setf (aref mi (1+ i)) (aref mi i)
		  (aref ma (1+ i)) (aref ma i)))
       (setf (aref mi pos) min
	     (aref ma pos) max)
       (incf n))
      ((= pos n)
       (setf (aref (mi r) n) min
	     (aref (ma r) n) max)
       (incf n))
      (t (error "~a" (list 'trying-to-write-behind-ranges 'pos pos 'n n)))))
  pos)

(defgeneric subtract (ranges min max))
(defmethod subtract ((r ranges) min max)
  (declare (num min max))
  (assert (<= min max))
  (with-slots (mi ma n) r
    (declare (fixnum n)
	     ((simple-array num 1) mi ma))
   (let ((eps .000001s0))
     (cond ((< +2pif+ min) 
	    (subtract r (- min +2pif+) (- max +2pif+)))
	   ((< max 0s0)
	    (subtract r (+ min +2pif+) (+ max +2pif+)))
	   ((< min 0s0)
	    (subtract r 0s0 max)
	    (subtract r (+ min +2pif+) +2pif+))
	   ((< +2pif+ max)
	    (subtract r min +2pif+)
	    (subtract r 0s0 (- max +2pif+)))
	   ((= 0 n)
	    nil)
	   (t 
	    (let ((pos 0))
	      (if (< min (aref mi 0))
		  (setf pos -1) ;; left has to go to front
		  (let ((lo 0)
			(mid 0)
			(hi n))
		    (declare (fixnum hi))
		    ;; binary search for segment whose mi is just smaller than min
		    (loop while (< lo (- hi 1)) do
			 (setf mid (floor (+ lo hi) 2))
			 (if (< (aref mi mid) min)
			     (setf lo mid)
			     (setf hi mid)))
		    (setf pos lo)))
	      (cond 
		((= -1 pos) ;; left of all segments
		 (setf pos 0))
		((< min (aref ma pos)) ;; new interval starts inside old segment 
		 (let ((mip (aref mi pos))
		       (map (aref ma pos)))
		   (declare (num mip map))
		   (if (< (- min mip) eps) ;; new interval starts at beginning of old segment
		       (if (< max map) ;; is shorter than old segment
			   (setf (aref mi pos) max) ; so let it start at new end
			   (delete-range r pos))
		       ;; there is a gap between new and old segment
		       (progn
			 (setf (aref ma pos) min) ;; max of old segment is beginning of new
			 (when (< max map) ;; new segment within old one
			   (insert r (1+ pos) max map))
			 (incf pos)))))
		((and (< pos (1- n)) ;; there are more old segments on the right
		      (< (aref mi (1+ pos)) max)) ;; new interval is extending into (or over) next old segment
		 (incf pos)))
	      ;; iterate towards right over old segments
	      (loop while (and (< pos n) ;; the loop decreases the fill-pointer!
			       (< (aref mi pos) max))
		 do
		 (if (< (- (aref ma pos) max) eps)
		     (delete-range r pos) ;; new segment overlaps this old one
		     (setf (aref mi pos) max) ;; new segment ends within this old one
		     )))))))
  nil)

#+nil
(let ((r (make-ranges)))
    (subtract r 1.1s0 1.2s0)
    (subtract r 0s0 2s0)
    r)



(defparameter *candidates* (make-array 30000 :element-type 'fixnum))
(defparameter *candidates-n* 0)

(defun make-candidates ()
  (setf *candidates* (make-array 30000 :element-type 'fixnum)
	*candidates-n* 0))

(defparameter *points* (make-array 100000 
				   :element-type '(simple-array single-float 1) 
				   :initial-element (v)))
(defparameter *points-n* 0)
(defparameter *grid* (make-array (list 1 1 1) :element-type 'fixnum))
(defparameter *points-per-cell* 7)
(declaim ((simple-array fixnum 1) *candidates*)
	 ((simple-array vec 1) *points*)
	 ((simple-array fixnum 3) *grid*)
	 (fixnum *points-per-cell* *points-n* *candidates-n*))

(defun pop-candidate ()
  (let* ((n *candidates-n*)
	 (i (if (= 0 n) 
		(return-from pop-candidate nil)
		(random n)))
	 (e (aref *candidates* i)))
    (declare ((simple-array fixnum 1) *candidates*))
    (setf (aref *candidates* i)
	  (aref *candidates* (1- n)))
    (decf *candidates-n*)
    #+nil (decf (fill-pointer queue))
    (the fixnum e)))


(defun push-point (p)
  (declare (vec p))
  #+nil (vector-push-extend p *points*)
  (setf (aref *points* *points-n*) p)
  (prog1
      *points-n*
    (incf *points-n*)))

(defun push-candidate (i)
  (declare (fixnum i))
  (setf (aref *candidates* *candidates-n*) i)
  (prog1
      *candidates-n*
    (incf *candidates-n*)))

(defmacro make-adjustable (type init)
  `(make-array 100000 
	       :element-type ',type 
	       :initial-element ,init
	       ;:adjustable t 
	       ;:fill-pointer 0
	       ))

(defun make-random-point ()
  (v 0s0 (- (random 2s0) 1) (- (random 2s0) 1)))


#+nil
(defun make-random-queue ()
  (the (array fixnum 1)
   (make-array 30000 :element-type 'fixnum
	       ;:adjustable t
	       ;:fill-pointer 0
	       )))



#+nil
(defun push-random (e queue)
  (declare (fixnum e)
	   ((array fixnum 1) queue))
  (vector-push-extend e queue))

#+nil
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

#+nil
(progn 
  (make-candidates)
  (dotimes (i 10)
    (push-candidate i))
  (dotimes (i 8)
    (format t "~a~%" (pop-candidate)))
  (subseq *candidates* 0 *candidates-n*))



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
	(error "grid overflowed max points per cell k=~d" k))))
  t)

#+nil
(progn
  (run 0.1)
  (print-grid :count t))


(defun add-point (p)
  (declare (vec p))
  (let ((ind (push-point p)))
    (add-to-grid p ind)
    (push-candidate ind)))

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
	 (distance-pixel-f (* 2 radius grid-size))
	 (distance-pixel (ceiling distance-pixel-f)
	   #+nil (ceiling (* 4 radius) grid-cell-size))
	 (n (min (floor grid-size 2) 
		 distance-pixel)))
    (declare (fixnum n)
	     (num distance-pixel-f)
	     ((signed-byte 64) distance-pixel))
    (multiple-value-bind (gy gx) (get-grid-point candidate)
      (let ((xside (if (< (* .5 grid-cell-size)
		       (- (vx candidate) (- (* gx grid-cell-size) 1)))
		       1 0))
	    (yside (if (< (* .5 grid-cell-size)
			  (- (vy candidate) (- (* gy grid-cell-size) 1)))
		       1 0)))
	(loop for j from (- n) upto n do
	    (let ((iy (cond ((= j 0) yside)
			    ((= j 1) 0)
			    (t 1))))
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
		     (declare (fixnum i j))
		     (when (< (+ (* dx dx) (* dy dy))
			      range2)
		       (let* ((ax (+ gx i grid-size))
			      (cx (mod ax grid-size))
			      ;; make sure the range of cx is 0 .. grid-size-1
			      (ay (+ gy j grid-size))
			      (cy (mod ay grid-size)))
			 (dotimes (k *points-per-cell*)
			   (let ((ind (aref *grid* cy cx k)))
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
					   (declare ((single-float 0s0 1s0) arg))
					   (subtract r
						     (- angle theta) 
						     (+ angle theta)))))))))))))))))))

(defgeneric make-random-angle (ranges))
(defmethod make-random-angle ((r ranges))
  (with-slots (mi ma n) r
    (declare (fixnum n))
    (when (= n 0)
      (error "range is empty"))
    (let* ((i (random n))
	   (min (aref mi i))
	   (max (aref ma i)))
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
  (make-candidates)
  (setf *points* (make-adjustable vec (v))
	*points-n* 0)
  (make-grid radius)
  (add-point (make-random-point))
  (let ((r (make-ranges)))
    (loop while (< 0 *candidates-n*) do
	(let* ((index (pop-candidate))
	       (candidate (get-point index))
	       (pi/3 (/ +pif+ 3s0)))
	  (reset r)
	  (find-neighbour-ranges r index radius)
	  (loop while (let ((n (n r)))
			(declare (fixnum n))
			(< 0 n))
	     do
	     (let* ((angle (make-random-angle r)) 
		    (pt (make-periphery-point
			 candidate angle radius)))
	       (declare (num angle))
	       (add-point pt)
	       (subtract r (- angle pi/3) (+ angle pi/3)))))))
  (subseq *points* 0 *points-n*))
#+nil
(progn
  (generate-poisson .01s0)
  (print-grid :count t)
  nil)

#+nil
(progn
  (time (progn (generate-poisson .008s0) nil))
#+nil (print-grid :count t))

#+nil
(require :sb-sprof)

#+nil
(progn
 (sb-sprof:with-profiling (:max-samples 1000
					:report :graph ;:flat
					:loop nil)
   (progn (generate-poisson .008)
	  nil)))
#+nil
(progn (generate-poisson .008)
       nil)
#+nil
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
(time
 (progn
   (run .1s0)
   #+nil (print-grid :count t)
   (terpri)))
