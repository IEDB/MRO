(ns mro.core
  (:require [clojure.string :as string]
            [clojure.set :as set]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [ontodev.excel :as xls]
            [mro.util :refer [read-csv write-csv to-identifier to-int
                              to-keyword clean-map split-list]]))

;; IEDB stores MHC information in a very dense SQL table.
;; This code takes that table and creates several new tables
;; that express the structure of MRO much more explicitly.


;; First we need to define some names and prefixes for taxa.

(def mhc-taxon-table
  [[1     nil    "organism"           "Any organism"]
   [9031  "BF"   "chicken"            "Gallus gallus"]
   [9483  "Caja" "marmoset"           "Callithrix jacchus"]
   [9490  "Saoe" "cotton-top tamarin" "Saguinus oedipus"]
   [9544  "Mamu" "rhesus macaque"     "Macaca mulatta"]
   [9593  "Gogo" "gorilla"            "Gorilla gorilla"]
   [9597  "Papa" "bonobo"             "Pan paniscus"]
   [9598  "Patr" "chimpanzee"         "Pan troglodytes"]
   [9606  "HLA"  "human"              "Homo sapiens"]
   [9615  "DLA"  "dog"                "Canis lupus familiaris"]
   [9796  "Eqca" "horse"              "Equus caballus"]
   [9796  "ELA"  "horse"              "Equus caballus"]
   [9823  "SLA"  "pig"                "Sus scrofa"]
   [9913  "BoLA" "cattle"             "Bos taurus"]
   [9940  "Ovar" "sheep"              "Ovis aries"]
   [10090 "H-2"  "mouse"              "Mus musculus"]
   [10116 "RT1"  "rat"                "Rattus novegicus"]])

(def taxon-mhc-prefixes
  (->> mhc-taxon-table
       (map (juxt first second))
       (into {})))

(def taxon-mhc-names
  (->> mhc-taxon-table
       (map (juxt first #(nth % 2)))
       (into {})))

(def taxon-mhc-ids
  (->> mhc-taxon-table
       (map (juxt #(nth % 2) first))
       (into {})))

(defn get-long-class
  "Given a class identifier string, return a nicer string for use in labels."
  [class]
  (when class
    (get {"i"  "MHC class I"
          "ii" "MHC class II"}
         (string/lower-case (name class))
         "non-classical MHC")))

(defn get-short-class
  "Given a class identifier string, return a short string for use in labels."
  [class]
  (when class
    (get {"i"  "class I"
          "ii" "class II"}
         (string/lower-case (name class))
         "non-classical")))

;; We use these definitions to adjust rows of the alleles.csv table.

(defn format-row
  "Make adjustments to a row map."
  [{:keys [synonyms class taxon-id] :as row}]
  (assoc row
         :synonyms    (split-list synonyms)
         :taxon-id    (to-int taxon-id)
         :class       (get-long-class class)
         :micro-class class
         :short-class (get-short-class class)
         :taxon-label (get taxon-mhc-names (to-int taxon-id))
         :prefix      (get taxon-mhc-prefixes (to-int taxon-id))))

;; The core of the system is a set of templates
;; for reformatting the data.
;; Some of the templates require special test functions
;; that we define here.
;; Each test function takes a row and determines
;; whether the template applies (a boolean test).

(defn chain-i-sublocus?
  "True only for rat An loci, where n is a number."
  [{:keys [taxon-id class chain-i-locus]}]
  (boolean
   (and (= "MHC class I" class)
        chain-i-locus
        (= taxon-id 10116)
        (re-find #"^A\d$" chain-i-locus))))

(defn chain-ii-sublocus?
  "True only for cow or human DRBn loci, where n is a number."
  [{:keys [taxon-id class chain-ii-locus]}]
  (boolean
   (and (= "MHC class II" class)
        chain-ii-locus
        (or (= taxon-id 9606) (= taxon-id 9913))
        (re-find #"^DRB\d$" chain-ii-locus))))

(defn mutant?
  [row]
  (or (:chain-i-mutation row) (:chain-ii-mutation row)))

;; This map holds all the special test functions.

(def special-functions
  {"chain-i sublocus?"
   chain-i-sublocus?
   "chain-ii sublocus?"
   chain-ii-sublocus?
   "not chain-i sublocus?"
   (comp not chain-i-sublocus?)
   "not chain-ii sublocus?"
   (comp not chain-ii-sublocus?)
   "chain-i sublocus parent"
   (fn [{:keys [prefix chain-i-locus]}]
     (format "%s-%s locus"
             prefix
             (string/replace chain-i-locus #"\d$" "")))
   "chain-ii sublocus parent"
   (fn [{:keys [prefix chain-ii-locus]}]
     (format "%s-%s locus"
             prefix
             (string/replace chain-ii-locus #"\d$" "")))
   "not prefix chain-i-locus?"
   (fn [{:keys [prefix chain-i-locus chain-i]}]
     (not (= (str prefix "-" chain-i-locus) chain-i)))
   "not prefix chain-ii-locus?"
   (fn [{:keys [prefix chain-ii-locus chain-ii]}]
     (not (= (str prefix "-" chain-ii-locus) chain-ii)))
   "human?"
   (fn [{:keys [prefix]}] (= "HLA" prefix))
   "not human?"
   (fn [{:keys [prefix]}] (not= "HLA" prefix))
   "alpha chain?"
   (fn [{:keys [prefix class chain-i-locus]}]
     (or (= "MHC class I" class)
         (and (= "HLA" prefix) (= "E" chain-i-locus))
         (and (= "HLA" prefix) (= "G" chain-i-locus))))
   "class II not sublocus?"
   (fn [{:keys [class] :as row}]
     (and (= "MHC class II" class)
          (not (chain-i-sublocus? row))
          (not (chain-ii-sublocus? row))))})

;; The hard work is applying templates to rows from the alleles table.
;; We define the templates in an Excel spreadsheet,
;; which lets us refer to other cells and keep our configuration DRY.

(def templates
  (->> (-> "src/mro/templates.xlsx"
           xls/load-workbook
           (xls/read-sheet "Sheet1"))
       (map clean-map)
       (map (fn [{:keys [depth test parent] :as row}]
              (merge
               row
               {:depth (to-int depth)}
               (when test
                 (if (find special-functions test)
                   {:test-name test
                    :test      (get special-functions test)}
                   (throw (Exception. (str "Could not find special function " test)))))
               (when (find special-functions parent)
                 {:parent-fn (get special-functions parent)}))))))

;; To locate specific templates
;; we use their branch, depth, and test-name (if given).
;; We create a map from these indices to templates.

(defn get-template-index
  "Given a template, return its index vector."
  [{:keys [branch depth test-name] :as template}]
  (if test-name
    [branch depth test-name]
    [branch depth]))

;; We use the index to create a map, e.g.
;; {["molecules" 1] {:branch "molecules" :depth 1 ...}}

(def template-map
  (zipmap (map get-template-index templates) templates))

(defn template-applies?
  "Given a template and a row,
   if the template has no test then return true,
   otherwise return the value of the test."
  [{:keys [test]} row]
  (if (nil? test) true (test row)))

;; The most important operation is substituting values into template strings.
;; This should work as you expect:
;;
;; - variables start with "$" and name a key
;; - every variable is replaces with the corresponding key from the row map

(defn fill-template
  "Given a template string and a map,
   replace `$keys` in the string with values from the map,
   and return the new string."
  [template-string row]
  (try
    (-> template-string
        (string/replace "-$" " @DASH@ $")
        (string/replace #"\$(\S+)" #(get row (keyword (second %))))
        (string/replace " @DASH@ " "-")
        (string/replace "Beta-2-microglobulin chain" "Beta-2-microglobulin"))
    (catch Exception e nil)))

;; We also need to be able to find the most-specific template
;; for a given row of data.

(defn assign-template
  "Given a row from the alleles.csv table,
   return an IRI for this node."
  [{:keys [level taxon-name] :as row}]
  (cond
    (= taxon-name "organism (all species)")
    (cond
      (= level "class")     ["molecules" 1]
      (= level "haplotype") ["haplotype-molecules" 1]
      (= level "serotype")  ["serotype-molecules" 1])
    (mutant? row)
    (cond
      (= level "partial molecule")  ["mutant-molecules" 3]
      (= level "complete molecule") ["mutant-molecules" 3])
    :else
    (cond
      (= level "class")             ["molecules" 2]
      (= level "locus")             ["molecules" 3 "alpha chain?"]
      (= level "haplotype")         ["haplotype-molecules" 2]
      (= level "serotype")          ["serotype-molecules" 2]
      (= level "partial molecule")  ["molecules" 4]
      (= level "complete molecule") ["molecules" 4])))

(defn assign-iri
  "Given a row, rename :id to :iedb-id,
   and use the most specific available template
   to assign a new :id IRI."
  [row]
  (merge
   (dissoc row :id)
   {:iedb-id (:id row)}
   (when-let [index (assign-template row)]
     (when-let [value (fill-template
                       (get-in template-map [index :label])
                       row)]
       {:id (str "MRO:" (to-identifier value))}))))

;; Now we define functions for applying a single template to a single row,
;; and for applying sets of templates to sequences of rows.

(defn apply-template
  "Given a template and a row map,
   apply the template and return the updated row map."
  [{:keys [branch level label test class-type parent parent-fn] :as template}
   row]
  (when-not label
    (throw (Exception. (str ":label key is required: " template))))
  (when-let [value (fill-template label row)]
    (apply
     merge
     (dissoc row :synonyms)
     (select-keys template [:branch :depth])
     {:label      value
      :id         (str "MRO:" (to-identifier value))
      :class-type (if class-type class-type "subclass")}
     (when (and level
                (:level row)
                (.contains (:level row) level)
                (if (mutant? row) (= branch "mutant-molecules") true)
                (seq (:synonyms row)))
       {:synonyms (string/join "|" (:synonyms row))})
     (when (and parent (not parent-fn))
       {:parent (fill-template parent row)})
     (when parent-fn
       {:parent (parent-fn row)})
     (for [key [:iedb-label :in-taxon :gene :alpha-chain :beta-chain
                :with-haplotype :with-serotype :mutant-of]]
       (when (find template key)
         {key (fill-template (get template key) row)})))))

(defn select-cols
  "Given a list of columns pairs and a row,
   return the row with just those column keys selected."
  [cols row]
  (select-keys
   row
   (->> cols
        (map (comp keyword to-identifier first))
        (concat [:branch :depth]))))

(defn merge-synonyms
  "Given a sequence of rows,
   if any two rows are identical except for their :synonyms and :iedb-id,
   merge them into one row."
  [rows]
  (->> rows
       (reduce
        (fn [coll row]
          (update-in coll [(:id row)] (fnil conj []) row))
        {})
       (reduce
        (fn [coll [id rows]]
          (cond
            (= 1 (count rows))
            (conj coll (first rows))
            (and (= 2 (count rows))
                 (= (dissoc (first rows)  :synonyms :iedb-id)
                    (dissoc (second rows) :synonyms :iedb-id)))
            (conj coll (apply merge rows))
            :else
            (println "BAD DUPLICATE FOR" id)))
        [])))

(defn warn-on-duplicate-ids
  [rows]
  (let [ids (->> rows (map :id) set)]
    (->> rows
         (map :id)
         frequencies
         (filter #(> (second %) 1))
         (map first)
         (map (partial println "DUPLICATE ID:"))
         doall))
  rows)

(defn apply-templates
  "Given a template group name, a vector of column-name/ROOBT-template pairs,
   and a sequence of rows,
   apply all the templates in the group to all the rows,
   add the ROBOT-template to the top,
   and write to a CSV file."
  [name cols rows]
  (->> (for [row      rows
             template (filter #(= name (:branch %)) templates)]
         (when (template-applies? template row)
           (apply-template template row)))
       (remove nil?)
       (map (partial select-cols cols))
       ; filter out the B2M chain
       (remove #(= (:id %) "MRO:beta-2-microglobulin"))
       set
       merge-synonyms
       warn-on-duplicate-ids
       (sort-by :label)
       (concat [(->> cols
                     (map (juxt (comp keyword to-identifier first) second))
                     (into {}))])
       (write-csv (format "mro-%s.csv" name) (map first cols))))

;; Finally, we specify which groups of templates run on which rows,
;; and how the resulting CSV file for each template group should be formatted.

(def default-cols
  [["ID"         "ID"]
   ["Label"      "A rdfs:label"]
   ["Synonyms"   "A IAO:0000118 SPLIT=|"]
   ["Class Type" "CLASS_TYPE"]
   ["Parent"     "C %"]
   ["In Taxon"   "C 'in taxon' some %"]])

(def default-iedb-cols
  (let [[before after] (split-at 2 default-cols)]
    (vec (concat before [["IEDB Label" "A OBI:9991118"]] after))))

(defn process-loci
  [rows]
  (apply-templates
   "loci"
   default-cols
   (->> rows
        (remove #(= 1 (:taxon-id %))))))

(defn process-haplotypes
  [rows]
  (apply-templates
   "haplotypes"
   default-cols
   (->> rows
        (filter :haplotype)
        (remove #(= 1 (:taxon-id %))))))

(defn process-serotypes
  [rows]
  (apply-templates
   "serotypes"
   default-cols
   (->> rows
        (filter :serotype)
        (remove #(= 1 (:taxon-id %))))))

(defn process-chains
  [rows]
  (apply-templates
   "chains"
   (conj (vec (drop-last default-cols))
         ["Gene" "C 'gene product of' some %"])
   (->> rows
        (remove #(= 1 (:taxon-id %))))))

(defn process-molecules
  [rows]
  (apply-templates
   "molecules"
   (conj default-iedb-cols
         ["Alpha Chain"    "C 'has part' some %"]
         ["Beta Chain"     "C 'has part' some %"]
         ["With Haplotype" "C 'haplotype member of' some %"]
         ["With Serotype"  "C 'serotype member of' some %"])
   (->> rows
        (filter #(contains? #{"complete molecule" "partial molecule"}
                            (:level %)))
        (remove #(= 1 (:taxon-id %)))
        ; TODO: rat data is broken
        (remove #(= 10116 (:taxon-id %))))))

(defn process-haplotype-molecules
  [rows]
  (apply-templates
   "haplotype-molecules"
   (conj default-iedb-cols
         ["With Haplotype" "C 'haplotype member of' some %"])
   (->> rows
        (filter :haplotype)
        (remove #(= 1 (:taxon-id %))))))

(defn process-serotype-molecules
  [rows]
  (apply-templates
   "serotype-molecules"
   (conj default-iedb-cols
         ["With Serotype" "C 'serotype member of' some %"])
   (->> rows
        (filter :serotype)
        (remove #(= 1 (:taxon-id %))))))

(defn process-mutant-molecules
  [rows]
  (apply-templates
   "mutant-molecules"
   (conj default-iedb-cols
         ["Mutant Of" "C 'is synthetic protein mutant of' some %"])
   (->> rows
        (filter mutant?))))

;; We generate an mro-iedb.csv file that updates alleles.csv
;; by assigning IRIs to rows,
;; and removing any rows that have not been used.

(def iedb-columns
  [["ID"         "ID"]
   ["IEDB ID"    "A MRO:has-iedb-mhc-id"]
   ["Label"      "A MRO:has-mhc-full-label"]
   ["Level"      "A MRO:has-mhc-restriction"]
   ["Taxon ID"   "A MRO:has-taxon-id"]
   ["Taxon Name" "A MRO:has-taxon-label"]
   ["Class"      "A MRO:has-mhc-class"]
   ["Locus"      "A MRO:has-mhc-locus"]
   ["Haplotype"  "A MRO:has-mhc-haplotype"]
   ["Serotype"   "A MRO:has-mhc-serotype"]])

(def iedb-keys (->> iedb-columns (map first) (map to-keyword)))

(def iedb-template (zipmap iedb-keys (map second iedb-columns)))

;; The next part is not elegant.
;; We use the information in the processed rows
;; to generated a bunch of extra IEDB annotations.
;; It would be cleaner if we were more consistent with naming.

(defn process-iedb-row
  "Given a map from IRIs to IEDB ID numbers and a row,
   return an updated row with extra IEDB annotations."
  [id-map
   {:keys [id branch depth parent in-taxon with-haplotype with-serotype]
    :as row}]
  (merge
   row
   {:level
    (condp = branch
      "molecules" (condp = depth
                    1 "class"
                    2 "class"
                    3 "locus"
                    "complete molecule")
      "haplotype-molecules" "haplotype"
      "serotype-molecules" "serotype"
      "mutant-molecules" (condp = depth
                           1 "class"
                           2 "class"
                           "complete molecule")
      nil)
    :class
    (condp = parent
      "MHC class I protein complex"       "MHC class I"
      "MHC class II protein complex"      "MHC class II"
      "non-classical MHC protein complex" "non-classical MHC"
      nil)}
   (when (find id-map id)
     (->> [:iedb-id :level :class :locus]
          (map (juxt identity #(get-in id-map [id %])))
          (into {})))
   (when in-taxon
     {:taxon-name in-taxon
      :taxon-id   (get taxon-mhc-ids in-taxon)})
   (when with-haplotype
     {:haplotype with-haplotype})
   (when with-serotype
     {:serotype with-serotype})))

(defn process-iedb
  "Given the original alleles.csv rows,
   and the newly generated rows,
   write a table of extra IEDB annotations."
  [old-rows new-rows]
  (->> new-rows
       (remove #(= "ID" (:id %)))
       (map (partial process-iedb-row
                     (->> old-rows
                          (map assign-iri)
                          (map (juxt :id identity))
                          (into {}))))
       (sort-by :label)
       (concat [iedb-template])
       (write-csv "mro-iedb.csv" (map first iedb-columns))))

;; Now we run all of these specific processing functions.

(defn process-table
  "Given a path to the alleles.csv file,
   generate a bunch of CSV files for various branches."
  [path]
  (let [rows (map format-row (read-csv path))]
    (->> (concat
          (process-loci rows)
          (process-haplotypes rows)
          (process-serotypes rows)
          (process-chains rows)
          (process-molecules rows)
          (process-haplotype-molecules rows)
          (process-serotype-molecules rows)
          (process-mutant-molecules rows))
         (process-iedb rows))))

(defn -main
  "Run process-table on alleles.csv."
  [& args]
  (process-table "alleles.csv"))
