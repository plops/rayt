(asdf:defsystem rayt
  :components ((:module "rayt"
                        :components
                        ((:file "packages")
                         (:file "base" :depends-on ("packages"))
                         (:file "raster" :depends-on ("packages" "base"))
			 (:file "poisson" :depends-on ("packages" "base"))
                         (:file "rayt" :depends-on ("packages" "base" "poisson"))))))
