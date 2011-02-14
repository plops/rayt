(in-package :raster)

(deftype bfp-type ()
  `(simple-array (unsigned-byte 64) 2))

(declaim (inline draw-point))
(defun draw-point (img j i val)
  (declare (type bfp-type img)
	   (type fixcoord j i)
	   (type (unsigned-byte 64) val))
  (destructuring-bind (h w) (array-dimensions img)
    (declare (type fixnum w h))
    (when (and (<= 0 i (1- w))
	       (<= 0 j (1- h)))
      (incf (aref img j i) val))))

(defun raster-line (img y x y1 x1 &optional (val 255))
  ;; wikipedia Bresenham's_line_algorithm
  (declare (type fixnum y x y1 x1)
	   (type bfp-type img)
	   (type (unsigned-byte 8) val))
  (let* ((dx (abs (- x1 x)))
	 (dy (abs (- y1 y)))
	 (sx (if (< x x1) 1 -1))
	 (sy (if (< y y1) 1 -1))
	 (e (- dx dy)))
    (loop
       (draw-point img y x val) 
       (when (and (= x x1)
		  (= y y1))
	 (return-from raster-line
	   (the bfp-type img)))
       (let ((e2 (* 2 e)))
	 (when (< (- dy) e2)
	   (decf e dy)
	   (incf x sx))
	 (when (< e2 dx)
	   (incf e dx)
	   (incf y sy))))))


(defun raster-circle (img y0 x0 r &optional (val 255))
  ;; wikipedia Midpoint_circle_algorithm
  (declare (type bfp-type img) 
	   (type fixnum x0 y0 r)
	   (type (unsigned-byte 64) val))
  (let ((f (- 1 r))
	(dx 1)
	(dy (* -2 r))
	(x 0)
	(y r))
    (declare (type fixnum f dx dy x y))
    (macrolet ((q (j i)
		 `(progn (draw-point img (+ y0 ,j) (+ x0 ,i) val)
			 (draw-point img (+ y0 ,i) (+ x0 ,j) val))))
      (q r 0)
      (q (- r) 0)
      (loop while (< x y) do
	   ;; dx = 2x+1
	   ;; dy = -2y
	   ;; f=x^2+y^2-r^2+2x-y+1
	   (when (<= 0 f)
	     (decf y)
	     (incf dy 2)
	     (incf f dy))
	   (incf x)
	   (incf dx 2)
	   (incf f dx)
	   (q y x)
	   (q y (- x))
	   (q (- y) x)
	   (q (- y) (- x)))))
  (the bfp-type img))

(defun raster-disk (img y0 x0 r &optional (val 255))
  (declare (type bfp-type img) (type (unsigned-byte 64) val) 
	   (type fixnum x0 y0 r))
  (let ((f (- 1 r))
	(dx 1)
	(dy (* -2 r))
	(x 0)
	(y r))
    (declare (type fixnum f dx dy x y))
    (macrolet ((q (y i i1)
		 `(raster-line img (+ y0 ,y) (+ x0 ,i1) (+ y0 ,y) (+ x0 ,i) val)))
      (q 0 r (- r))
      (loop while (< x y) do
	   ;; dx = 2x+1
	   ;; dy = -2y
	   ;; f=x^2+y^2-r^2+2x-y+1
	   (when (<= 0 f)
	     (decf y)
	     (incf dy 2)
	     (incf f dy))
	   (incf x)
	   (incf dx 2)
	   (incf f dx)
	   (q y (- x) x)
	   (q (- y) (- x) x)
	   (q x y (- y))
	   (q (- x) y (- y)))))
  (the bfp-type img))

(declaim (optimize (speed 3) (safety 0) (debug 1)))

;; graphics gems III pp. 599
(deftype fixcoord ()
  `(integer -10000 10000))
(defstruct edge
  (ymin 0 :type fixcoord)
  (ymax 0 :type fixcoord)
  (xi 0 :type fixcoord)
  (si 0 :type fixcoord)
  (r 0 :type fixcoord)
  (inc 0 :type fixcoord)
  (dec 0 :type fixcoord))

(defun edge-setup (y0 x0 y1 x1)
  (declare (type fixcoord y0 x0 y1 x1))
  (the edge
   (let ((dx (- x1 x0))
	 (dy (- y1 y0)))
     (declare (type fixcoord dx dy))
     (if (= 0 dy)
	 (make-edge :ymin y0 :ymax y1)
	 (let* ((si (floor dx dy))
		(sf (- dx (* si dy))))
	   (declare (type fixcoord si) 
		    (type fixcoord sf))
	   (make-edge :ymin y0 :ymax y1
		      :si si :xi (+ x0 si)
		      :r (- (* 2 sf) dy)
		      :inc sf
		      :dec (- sf dy)))))))

(defun edge-scan (edge)
  (declare (type edge edge))
  (the fixcoord 
    (with-slots (xi r si dec inc) edge
      (let ((x xi))
	(if (<= 0 r)
	    (progn (incf xi (1+ si))
		   (incf r dec))
	    (progn (incf xi si)
		   (incf r inc)))
	x))))


(declaim (inline draw-span))
(defun draw-span (img y x1 x2 val)
  (declare (type fixcoord y x1 x2)
	   (type bfp-type img) (type (unsigned-byte 64) val) )
  (loop for x from (1+ x1) upto x2 do 
       (draw-point img y x val)))

;; triangles must be sorted y0 <= y1 <= y2
(defun sorted-triangle (img y0 x0 y1 x1 y2 x2 val)
  (declare (type fixcoord y0 x0 y1 x1 y2 x2)
	   (type bfp-type img) (type (unsigned-byte 64) val) )
  (let* ((handedness (- (* (- y1 y0)
			   (- x2 x0))
			(* (- x1 x0)
			   (- y2 y0))))
	 (left (if (< handedness 0)
		   (edge-setup y0 x0 y2 x2)
		   (edge-setup y0 x0 y1 x1)))
	 (right (if (< handedness 0)
		    (edge-setup y0 x0 y1 x1)
		    (edge-setup y0 x0 y2 x2))))
    (declare (type fixcoord handedness))
    (loop for y from (1+ (edge-ymin left)) upto (min (edge-ymax left)
						     (edge-ymax right)) do
	 (draw-span img y (edge-scan left) (edge-scan right) val))
    (if (<= 0 handedness)
	(setf left (edge-setup y1 x1 y2 x2))
	(setf right (edge-setup y1 x1 y2 x2)))
    (loop for y from (1+ (max (edge-ymin left)
			      (edge-ymin right)))
	 upto (edge-ymax left) do
	 (draw-span img y (edge-scan left) (edge-scan right) val))
    img))

(defun raster-triangle (img y0 x0 y1 x1 y2 x2 &optional (val 255))
  (declare (type fixcoord y0 x0 y1 x1 y2 x2)
	   (type bfp-type img) (type (unsigned-byte 64) val))
  (when (< y1 y0)
    (rotatef y0 y1)
    (rotatef x0 x1))
  (when (< y2 y0)
    (rotatef y0 y2)
    (rotatef x0 x2))
  (when (< y2 y1)
    (rotatef y1 y2)
    (rotatef x1 x2))
  (sorted-triangle img y0 x0 y1 x1 y2 x2 val))


#+nil
(let ((m (make-array (list 300 300)
		     :element-type '(unsigned-byte 8))))
  (dotimes (i 30)
   (raster-circle m 150 150 (* 3 i)))
  (raster-line m 12 13 150 190)
  (raster-disk m 200 200 40)
  (raster-triangle m 10 10 30 30 300 400)
  (raster-triangle m 12 20 30 40 60 12)
  (rayt::write-pgm "/dev/shm/o.pgm" m))

