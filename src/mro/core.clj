(ns mro.core
  (:require [clojure.string :as string]
            [clojure.set :as set]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [ontodev.excel :as xls]
            [mro.util :refer [read-csv write-csv to-identifier to-int
                              clean-map split-list]]))

;; IEDB stores MHC information in a very dense SQL table.
;; This code takes that table and creates several new tables
;; that express the structure of MRO much more explicitly.)


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


;; We use these definitions to adjust rows of the alleles table.

(defn format-row
  "Make adjustments to a row map."
  [{:keys [synonyms includes class taxon-id] :as row}]
  (assoc row
         :synonyms    (set/union
                       (split-list synonyms)
                       (split-list includes))
         :taxon-id    (to-int taxon-id)
         :class       (get-long-class class)
         :micro-class class
         :short-class (get-short-class class)
         :taxon-label (get taxon-mhc-names (to-int taxon-id))
         :prefix      (get taxon-mhc-prefixes (to-int taxon-id))))

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

; Define some fnctions to use in the template sheet
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
       (map (fn [{:keys [test parent] :as row}]
              (merge
               row
               (when test
                 (if (find special-functions test)
                   {:test (get special-functions test)}
                   (throw (Exception. (str "Could not find special function " test)))))
               (when (find special-functions parent)
                 {:parent-fn (get special-functions parent)}))))))

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
        (concat []))))

(defn merge-synonyms
  "Given a sequence of rows,
   if any two rows are identical except for their :synonyms,
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
                 (= (dissoc (first rows)  :synonyms)
                    (dissoc (second rows) :synonyms)))
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
         ["Alpha Chain" "C 'has part' some %"]
         ["Beta Chain"  "C 'has part' some %"]
         ["With Haplotype" "C 'haplotype member of' some %"]
         ["With Serotype" "C 'serotype member of' some %"])
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
         ["Mutant Of" "C 'mutant of' some %"])
   (->> rows
        (filter mutant?))))

(defn process-table
  "Given a path to the alleles.csv file,
   generate a bunch of CSV files for various branches."
  [path]
  (let [rows (map format-row (read-csv path))]
    (concat
     (process-loci rows)
     (process-haplotypes rows)
     (process-serotypes rows)
     (process-chains rows)
     (process-molecules rows)
     (process-haplotype-molecules rows)
     (process-serotype-molecules rows)
     (process-mutant-molecules rows))))

(defn -main
  "Run process-table on alleles.csv."
  [& args]
  (process-table "alleles.csv"))
