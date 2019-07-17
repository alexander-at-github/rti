#! /usr/bin/env racket
#lang racket/base
(require
  (only-in racket/cmdline
           command-line)
  (only-in racket/file
           file->lines)
  (only-in racket/list
           eighth
           empty?)
  (only-in racket/string
           string-prefix?
           string-split
           string-trim))

;; This script reads an input file, skips lines starting with a hash symbol, and
;; returns a list of distinct elements within a certain column of the input file.
;; The particular column can be set by defining the option "position".

;; Options
(define position eighth)

;; parameters
(define infile-name (make-parameter null))

;; read command line arguments
(command-line
  #:args (one-arg) (infile-name one-arg))

(define (handle-line line acc)
  (let* ([tokens (string-split line)]
         [tok (position tokens)])
    (if (not (member tok acc))
      (cons tok acc)
      acc)))

(define (handle-lines lines acc)
  (if (empty? lines)
    acc
    (handle-lines (cdr lines) (handle-line (car lines) acc))))

(define (filter-comments lines)
  (filter
    (lambda (aa)
      ;(cond [(string-prefix? (string-trim aa) "#") (printf "~a\n" (string-trim aa))]) ;; debug
      (not (string-prefix? (string-trim aa) "#")))
    lines))

(define (handle file-name)
  (handle-lines (filter-comments (file->lines file-name)) (list)))

(define result (handle (infile-name)))

(printf "~a\n" result)
