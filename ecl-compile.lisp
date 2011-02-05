(require :asdf)
(asdf:make-build :rayt :type :program :epilogue-code '(ext:quit 0)
		 :move-here t)
(quit)