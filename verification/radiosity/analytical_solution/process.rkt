#! /usr/bin/env racket
#lang racket/base

(require (only-in racket/function
                  identity)
         (only-in racket/list
                  argmax)
         (only-in racket/list
                  first
                  second))
(require csv-reading)

(define make-reader
  (make-csv-reader-maker
    '((separator-chars #\ ))))

(define next-row-x
  (make-reader (open-input-file "original/position.sp.0.1.ar.45.txt")))
(define xx
  (let ([nrx (map string->number (next-row-x))])
    (map (lambda (oo) (/ oo (argmax identity nrx))) nrx)))

(define next-row-y
  (make-reader (open-input-file "original/flux.sp.0.1.ar.45.txt")))
(define yy
  (let ([nry (map string->number (next-row-y))])
    (map (lambda (oo) (/ oo (argmax identity nry))) nry)))

(define zip (lambda (l1 l2) (map list l1 l2)))

(for ([pair (zip xx yy)])
  (printf "~a ~a\n" (first pair) (second pair)))
