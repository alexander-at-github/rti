#! /usr/bin/env racket
#lang racket/base
(require (only-in racket/cmdline
                  command-line)
         (only-in racket/file
                  file->lines)
         (only-in racket/function
                  identity)
         (only-in racket/list
                  append*
                  append-map
                  drop
                  empty?
                  first
                  second)
         (only-in math/statistics
                  mean
                  variance))

;; Options
;; Number of measurments droped from each sample
(define drop-num 3)

;; define parameters
(define file-paths-a (make-parameter null))
(define machine-readable-mode (make-parameter #f))

;; read command line argument and write it to the parameter
(command-line
  #:multi [("-m" "--machine-readable") "Produce machine reable output" (machine-readable-mode #t)]
  #:args args (file-paths-a args))

;; Filter file paths
(define file-paths
  (let ([exists (lambda (file-path)
                  (if (file-exists? file-path)
                      #t
                      ((lambda () (printf "Cannot find file ~a\n" file-path) #f))))])
    (filter exists (file-paths-a))))

(define (m-match input mstr)
  (regexp-match* (pregexp mstr) input)) ;; regexp-match* matches multiple occurences in input

(define (get-matching lines str)
  (let ([match (lambda (line) (m-match line str))])
    (append-map match lines)))

(define (m-mean lst)
  (exact->inexact (mean lst)))

(define (m-variance lst)
  (exact->inexact (variance lst #:bias #t))) ;; the #:bias #t argument turns bias correction on (n vs n-1).

(define (format-human-summary num-threads file-path filter-str sample-size sample)
  (let ([mu-hat (m-mean sample)]
        [var-hat (m-variance sample)])
    (format "number of threads: ~a file path: ~a filter: ~a sample mean: ~a sample variance: ~a confidence interval of mean: ~a"
            num-threads
            file-path
            filter-str
            mu-hat
            var-hat
            (mu-confidence-95 sample-size mu-hat var-hat)
            )))

(define (format-machine-summary num-threads file-path filter-str sample-size sample)
  (let* ([mu-hat (m-mean sample)]
         [var-hat (m-variance sample)]
         [mu-hat-confidence (mu-confidence-95 sample-size mu-hat var-hat)]
         [mu-hat-low (first mu-hat-confidence)]
         [mu-hat-high (second mu-hat-confidence)])
    (string-append
     (format "# numbr of threads: ~a\n" num-threads)
     (format "# file-path: ~a\n" file-path)
     (format "# fitler string: ~a\n" filter-str)
     (format "#\n")
     (format "# number of threads, sample mean, 95 confidence interval low, 95 confidence interval high\n")
     (format "~a ~a ~a ~a\n" num-threads mu-hat mu-hat-low mu-hat-high))
  ))

(define (m-drop lst num)
  (if (> (length lst) num)
      (drop lst num)
      '())) ;; otherwise return empty list

;; Defining t-0.025 for several degrees of freedom (95% confidence)
(define (t-0.025 freedom)
  ;;(printf "t-0.025 freedom: ~a\n" freedom)
  (hash-ref
   ;; t-distribution values from a table
   #hash((29 . 2.045) (30 . 2.042) (31 . 2.040) (32 . 2.037) (33 . 2.035) (34 . 2.032) (35 . 2.030))
   freedom))

(define (mu-confidence-95 sample-size mu-hat var-hat)
  (let ([range (* (t-0.025 (sub1 sample-size)) (sqrt (/ var-hat sample-size)))])
    (list (- mu-hat range) (+ mu-hat range))
    ))

(for* ([filter-str '(
                     "rti::triangle_geometry_from_gmsh"
                     "rti::sphere_geometry_from_gmsh"
                     "rti::oriented_disc_geometry_from_gmsh"
                     "rti::disc_geometry_from_gmsh"
                     )]
       [file file-paths]
       )
  (let* (;; parse file into lines
         [lines (file->lines file)]
         ;; extract the datum specifying the number of threads used
         [num-threads1 (get-matching (get-matching lines "Using\\s\\d+\\sthreads") "\\d+")]
         ;; If there is no information about the number of threads, then return zero.
         [num-threads2 (if (null? num-threads1) "0" (car num-threads1))]
         ;; extract the running times numbers
         [m-lines (get-matching lines (format "\\(.*~a.*\\)" filter-str))]
         [time-list1 (get-matching m-lines "\\s\\d+ns")]
         [time-list2 (get-matching time-list1 "\\d+")]
         [time-list3 (m-drop time-list2 drop-num)]
         [samples (map string->number time-list3)]
         ;; size of sample
         [sample-size (length samples)]
         ;; Script output format
         [format-function (if (machine-readable-mode) format-machine-summary format-human-summary)]
         )
    ;; (println "m-lines")
    ;; (println m-lines)
    ;; (println "time-list1")
    ;; (println time-list1)
    ;; (println "time-list2")
    ;; (println time-list2)
    ;; (println num-threads2)
    ;; (printf "running times for file ~a: ~a\n" file samples)
    (printf "# droped ~a measurments from each sample\n" drop-num)
    (display (format-function num-threads2 file filter-str sample-size samples))
    (printf "#\n")
    ))

