(asdf:defsystem rayt
  :components ((:module "rayt"
                        :components
                        ((:file "packages")
                         (:file "base" :depends-on ("packages"))
                         (:file "raster" :depends-on ("packages" "base"))
			 (:file "poisson" :depends-on ("packages" "base"))
			 (:file "rayt-model" :depends-on ("packages" "base" "raster"))
			 (:file "simple-rayt" :depends-on ("packages" "base" "raster"
								      "rayt-model"))
                         (:file "rayt" :depends-on ("packages" "base" "poisson"
							       "raster" "rayt-model"))))))
