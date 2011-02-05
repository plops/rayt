(require :asdf)
(require :rayt)

(defpackage :run
  (:use :cl :rayt))
(in-package :run)

(sum-bfp 0)