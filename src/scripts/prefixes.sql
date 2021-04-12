CREATE TABLE IF NOT EXISTS prefix (
  prefix TEXT PRIMARY KEY,
  base TEXT NOT NULL
);

INSERT OR IGNORE INTO prefix VALUES
("BFO",       "http://purl.obolibrary.org/obo/BFO_"),
("ECO",       "http://purl.obolibrary.org/obo/ECO_"),
("GO",        "http://purl.obolibrary.org/obo/GO_"),
("IAO",       "http://purl.obolibrary.org/obo/IAO_"),
("OBI",       "http://purl.obolibrary.org/obo/OBI_"),
("oio",       "http://www.geneontology.org/formats/oboInOwl#"),
("owl",       "http://www.w3.org/2002/07/owl#"),
("PR",        "http://purl.obolibrary.org/obo/PR_"),
("REO",       "http://purl.obolibrary.org/obo/REO_"),
("rdf",       "http://www.w3.org/1999/02/22-rdf-syntax-ns#"),
("rdfs",      "http://www.w3.org/2000/01/rdf-schema#"),
("RO",        "http://purl.obolibrary.org/obo/RO_"),
("HANCESTRO", "http://purl.obolibrary.org/obo/HANCESTRO_"),
("NCIT",      "http://purl.obolibrary.org/obo/NCIT_"),
("IDO",       "http://purl.obolibrary.org/obo/IDO_");
