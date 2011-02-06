(defpackage :base
  (:use :cl)
  (:export 
   #:num #:vec #:with-arrays #:v #:copy-vec
   #:.+ #:.- #:.* #:./ #:dot #:norm #:normalize
   #:cross #:.s #:vx #:vy #:vz #:v-spherical
   #:check-unit-vector #:check-range
   #:req #:+pif+ #:+2pif+))

(defpackage :poisson
  (:use :cl :base)
  (:export #:generate-poisson))

(defpackage :raster
  (:use :cl :base)
  (:export ))

(defpackage :rayt
  (:use :cl :base :poisson)
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