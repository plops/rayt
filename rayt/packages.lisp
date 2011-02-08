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

(defpackage :rayt-model
  (:use :cl :base :raster)
  (:export
   #:*centers-fix*
   #:*centers*
   #:*data-dz-mm*
   #:*data-dx-mm*
   #:*data-width-px*
   #:draw-nucleus
   #:draw-slice-through-point
   #:find-bfp-radius
   #:find-focal-length))

(defpackage :simple-rayt
  (:use :cl :base :raster :rayt-model)
  (:export   
   #:project-ray-into-bfp
   #:project-nucleus-into-bfp
   #:sum-bfp-raster))

(defpackage :rayt
  (:use :cl :base :poisson :raster :rayt-model)
  (:export
   #:sum-bfp))