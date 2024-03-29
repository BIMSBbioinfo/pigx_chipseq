;;; PiGx_chipseq - reports pipeline for reads from ChIPseq experiments.
;;; Copyright © 2017, 2021 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx_chipseq.
;;;
;;; PiGx_chipseq is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment for
;;; PiGx_chipseq:
;;;
;;;  $ guix environment -l guix.scm
;;;
;;; To install the package from a release tarball do this:
;;;
;;;  $ guix package --with-source=pigx_chipseq-0.0.1.tar.gz -f guix.scm
;;;
;;; This environment file was developed for Guix version
;;; v0.14.0-3177-gbcddf30af


(use-modules (guix packages) (gnu packages))

(define %pigx-chipseq-version
  (symbol->string (with-input-from-file "VERSION" read)))

(define-public pigx-chipseq-development
  (let ((parent (specification->package "pigx-chipseq")))
    (package (inherit parent)
      (version %pigx-chipseq-version)
      (source (string-append (getcwd) "/pigx_chipseq-" version ".tar.gz"))
      (native-inputs
       `(("autoconf" ,(specification->package "autoconf"))
         ("automake" ,(specification->package "automake"))
         ,@(package-native-inputs parent))))))

pigx-chipseq-development
