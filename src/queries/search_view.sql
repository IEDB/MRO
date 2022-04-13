-- Create a view with search labels as LABEL - SYNONYM [TERM ID]
DROP VIEW IF EXISTS mro_search_view;
CREATE VIEW mro_search_view AS
WITH term_ids AS (
    SELECT * FROM (
        SELECT DISTINCT subject AS subject FROM mro
        UNION
        SELECT DISTINCT predicate FROM mro
    )
),
labels AS (
    SELECT DISTINCT subject, object
    FROM mro WHERE predicate = 'rdfs:label'
),
synonyms AS (
    SELECT * FROM (
        SELECT DISTINCT subject, object FROM mro
        WHERE predicate = 'IAO:0000118'
        UNION
        SELECT DISTINCT subject, object FROM mro
        WHERE predicate = 'OBI:9991118'
    )
)
SELECT
    t.subject AS subject,
    COALESCE(l.object, "") || COALESCE(" - " || s.object, "") || " [" || t.subject || "]" AS label
FROM term_ids t
LEFT JOIN labels l ON t.subject = l.subject
LEFT JOIN synonyms s ON t.subject = s.subject;
