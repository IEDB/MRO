CREATE TABLE IF NOT EXISTS prefix (
  prefix TEXT PRIMARY KEY,
  base TEXT NOT NULL
);

INSERT OR IGNORE INTO prefix VALUES
('rdf', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'),
('rdfs', 'http://www.w3.org/2000/01/rdf-schema#'),
('xsd', 'http://www.w3.org/2001/XMLSchema#'),
('owl', 'http://www.w3.org/2002/07/owl#'),
('oio', 'http://www.geneontology.org/formats/oboInOwl#'),
('dce', 'http://purl.org/dc/elements/1.1/'),
('dct', 'http://purl.org/dc/terms/'),
('foaf', 'http://xmlns.com/foaf/0.1/'),
('obo',  'http://purl.obolibrary.org/obo/'),

('BFO', 'http://purl.obolibrary.org/obo/BFO_'),
('ECO', 'http://purl.obolibrary.org/obo/ECO_'),
('GO', 'http://purl.obolibrary.org/obo/GO_'),
('IAO', 'http://purl.obolibrary.org/obo/IAO_'),
('OBI', 'http://purl.obolibrary.org/obo/OBI_'),
('MRO', 'http://purl.obolibrary.org/obo/MRO_'),
('PR', 'http://purl.obolibrary.org/obo/PR_'),
('RO', 'http://purl.obolibrary.org/obo/RO_');
