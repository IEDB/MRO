(ns mro.core
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [ontodev.excel :as xls]
            [mro.util :refer [read-csv write-csv to-identifier to-int
                              clean-map]]))

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

(defn clean-class
  "Given a class identifier string, return a nicer string for us in labels."
  [class]
  (when class
    (get {"i"  "MHC class I"
          "ii" "MHC class II"}
         (string/lower-case (name class))
         "non-classical MHC")))


;; We use these definitions to adjust rows of the alleles table.

(defn format-row
  "Make adjustments to a row map."
  [{:keys [class taxon-id] :as row}]
  (assoc row
         :taxon-id    (to-int taxon-id)
         :class       (clean-class class)
         :taxon-label (get taxon-mhc-names (to-int taxon-id))
         :prefix      (get taxon-mhc-prefixes (to-int taxon-id))))

(defn special-locus?
  "True only if the chain-ii locus is DRBn where n is a number"
  [{:keys [class chain-ii-locus]}]
  (and (= "MHC class II" class)
       chain-ii-locus
       (re-find #"DRB\d$" chain-ii-locus)))

; Define some fnctions to use in the template sheet
(def special-functions
  {"special locus?"
   special-locus?
   "not special locus?"
   (comp not special-locus?)
   "special locus parent"
   (fn [{:keys [prefix chain-ii-locus]}]
     (format "%s-%s locus"
             prefix
             (string/replace chain-ii-locus #"\d$" "")))
   "human?"
   (fn [{:keys [prefix]}] (= "HLA" prefix))
   "not human?"
   (fn [{:keys [prefix]}] (not= "HLA" prefix))
   "alpha chain?"
   (fn [{:keys [prefix class chain-i-locus]}]
     (or (= "MHC class I" class)
         (and (= "HLA" prefix) (= "E" chain-i-locus))
         (and (= "HLA" prefix) (= "G" chain-i-locus))))
   "class II not special?"
   (fn [{:keys [class] :as row}]
     (and (= "MHC class II" class)
          (not (special-locus? row))))})


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
               (when test {:test (get special-functions test)})
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
  [{:keys [label test class-type parent parent-fn] :as template}
   row]
  (when-not label
    (throw (Exception. (str ":label key is required: " template))))
  (when-let [value (fill-template label row)]
    (apply
     merge
     row
     {:label      value
      :id         (str "MRO:" (to-identifier value))
      :class-type (if class-type class-type "subclass")}
     (when (and parent (not parent-fn))
       {:parent (fill-template parent row)})
     (when parent-fn
       {:parent (parent-fn row)})
     (for [key [:in-taxon :gene :alpha-chain :beta-chain
                :with-haplotype :with-serotype :mutant-of]]
       (when (find template key)
         {key (fill-template (get template key) row)})))))

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
       (map #(select-keys % (map (comp keyword to-identifier first) cols)))
       set
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
   ["Class Type" "CLASS_TYPE"]
   ["Parent"     "C %"]
   ["In Taxon"   "C 'in taxon' some %"]])

(defn process-table
  "Given a path to an "
  [path]
  (let [rows (map format-row (read-csv path))]
    (apply-templates
     "loci"
     default-cols
     (->> rows
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "haplotypes"
     default-cols
     (->> rows
          (filter :haplotype)
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "serotypes"
     default-cols
     (->> rows
          (filter :serotype)
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "chains"
     (conj (vec (drop-last default-cols))
           ["Gene" "C 'gene product of' some %"])
     (->> rows
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "molecules"
     (conj default-cols
           ["Alpha Chain" "C 'has part' some %"]
           ["Beta Chain"  "C 'has part' some %"]
           ["With Haplotype" "C 'haplotype member of' some %"]
           ["With Serotype" "C 'serotype member of' some %"])
     (->> rows
          (filter #(contains? #{"complete molecule" "partial molecule"}
                              (:level %)))
          (remove #(= 1 (:taxon-id %))) 
          ; TODO: rat data is broken
          (remove #(= 10116 (:taxon-id %)))))
    (apply-templates
     "haplotype-molecules"
     (conj default-cols
           ["With Haplotype" "C 'haplotype member of' some %"])
     (->> rows
          (filter :haplotype)
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "serotype-molecules"
     (conj default-cols
           ["With Serotype" "C 'serotype member of' some %"])
     (->> rows
          (filter :serotype)
          (remove #(= 1 (:taxon-id %)))))
    (apply-templates
     "mutant-molecules"
     (conj default-cols
           ["Mutant Of" "C 'mutant of' some %"])
     (->> rows
          (filter #(or (:chain-i-mutation %) (:chain-ii-mutation %)))))))

(defn -main
  "Run process-table on alleles.csv."
  [& args]
  (process-table "alleles.csv"))
