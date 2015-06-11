(ns mro.util
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]))

;; This file defines some utilities we need,
;; mostly for reading and writing CSV files.

(defn to-identifier
  "Take a string and return a properly formatted identifier string."
  [s]
  (-> s
      string/trim
      string/lower-case
      (string/replace "_" "-")
      (string/replace "/" "-")
      (string/replace #"\s+" "-")))

(defn to-keyword
  "Take a string and return a properly formatted keyword."
  [s]
  (keyword (to-identifier s)))

(defn to-int
  "Try to convert a value to an integer. Nil is simply returned."
  [x]
  (cond
    (nil? x) nil
    (instance? Integer x) x
    (number? x) (int x)
    (and (string? x) (string/blank? x)) nil
    (and (string? x) (re-matches #"\d+" x)) (Integer/parseInt x)
    :else
    (throw (Exception. (str "Input is not an integer: " x)))))

(defn read-rows
  "Given sequence of sequences where the first sequence contains headers,
   return a sequence of maps with the headers as keys."
  [rows]
  (map (partial zipmap (map to-keyword (first rows)))
       (rest rows)))

(defn clean-map
  "Given a map, return a map for which none of the values
   is nil or a blank string."
  [coll]
  (->> coll
       (remove (comp string/blank? second))
       (into {})))

(defn read-csv
  "Load data from the comma-separated values file at the given path."
  [path]
  (->> path
       slurp
       csv/read-csv
       read-rows
       (map clean-map)))

(defn write-csv
  "Given a path, a list of headers, and a vector of maps,
   translate the headers to keywords, then write a CSV file
   to the path using the headers as the list of columns.
   Return the data."
  [path headers rows]
  (with-open [w (io/writer path)]
    (csv/write-csv w [headers])
    (doseq [row rows]
      (csv/write-csv w [(map #(get row (to-keyword %) "") headers)])))
  rows)

;; These functions are help fetch sequences and merge them into chain data.

(defn read-fasta
  "Read a FASTA file into a map from descriptions to sequences."
  [path]
  (->> path
       slurp
       string/split-lines
       (partition-by #(.startsWith % ">"))
       (map (partial string/join ""))
       (apply hash-map)))

(defn add-sequences
  "Add sequences from a FASTA file to a table of chains."
  [fasta-path chain-path]
  (let [seqs (read-fasta fasta-path)]
    (->> chain-path
         slurp
         csv/read-csv
         read-rows
         (map (fn [{:keys [accession-protein] :as row}]
                (assoc
                 row
                 :sequence
                 (->> seqs
                      (filter #(.contains (first %) accession-protein))
                      first
                      second))))
         (write-csv "mro-chains-seqs.csv"
                    ["Name" "Resource Name"
                     "Source" "Accession DNA/RNA" "Accession Protein"
                     "Sequence"]))))

(defn add-ids
  "Add properly-formated IDs to CSV file for chains."
  [path]
  (->> path
       slurp
       csv/read-csv
       read-rows
       (map (fn [{:keys [name] :as row}]
              (assoc row :id (str "MRO:" (to-identifier name)))))
       (write-csv "mro-chains.csv"
                  ["ID" "Name" "Resource Name"
                   "Source" "Accession DNA/RNA" "Accession Protein"
                   "Sequence"])))
