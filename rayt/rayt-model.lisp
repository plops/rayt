(in-package :rayt-model)

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


(defun find-bfp-radius (na f)
  (declare (type num na f)) (the num (* f na)))

(defun find-focal-length (mag)
  (declare (type num mag))  (the num (/ 164.5 mag)))

(defun get-nuclei-in-slice (slice-px r-px)
  "Return a list of nuclei that should be visible in the given slice."
  (declare (type fixnum slice-px r-px))
  (let ((res ()))
   (dotimes (nucleus (length *centers-fix*))
     (let ((c (aref *centers-fix* nucleus)))
       (when (<= (abs (- slice-px (vz c))) r-px)
	 (push nucleus res))))
   res))

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
       (raster-disk img y x rad 56))))
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
