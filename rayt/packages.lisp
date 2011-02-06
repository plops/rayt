(defpackage :base
  (:use :cl)
  (:export 
   #:num #:vec #:with-arrays #:v #:copy-vec
   #:.+ #:.- #:.* #:./ #:dot #:norm #:normalize
   #:cross #:.s #:vx #:vy #:vz #:v-spherical
   #:check-unit-vector #:check-range
   #:req #:+pi+ #:+2*pi+ #:+pi/2+
   #:mat
   #:m
   #:rotation-matrix
   #:m*))

(defpackage :poisson
  (:use :cl :base)
  (:export #:generate-poisson))

(defpackage :raster
  (:use :cl :base)
  (:export
   #:raster-line
   #:raster-circle
   #:raster-disk
   #:raster-triangle))

(defpackage :rayt
  (:use :cl :base :poisson :raster)
  (:export
   #:*centers-fix*
   #:*centers*
   #:*data-dz-mm*
   #:*data-dx-mm*
   #:*data-width-px*
   #:raster-circle
   #:raster-disk
   #:raster-line
   #:draw-nucleus
   #:draw-slice-through-point
   #:sum-bfp))