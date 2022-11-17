use FSM_Analytics;

SELECT DISTINCT bc.*
				, m.medication_name
				, m.generic_name
				, med.order_start_date_key
FROM [edw_research_dm].[luo_khan_nlp_breast_cancer_export] bc
LEFT JOIN [NM_BI].[dim].[vw_patient_current] p
	ON bc.patient_ir_id = p.ir_id
LEFT JOIN [NM_BI].[fact].[vw_medication_order] med
	ON p.patient_key = med.patient_key
	AND (med.order_start_datetime < bc.date_of_diagnosis)
	AND YEAR(med.order_start_datetime) < 2017
LEFT JOIN [nm_bi].[dim].[vw_medication] m                          --adding medication dimension table
	ON m.medication_key = med.medication_key
	AND (m.medication_name LIKE '%metformin%' or m.generic_name  LIKE '%metformin%' ) --search on medication and generic names for better results
order by bc.patient_ir_id;

SELECT *
FROM [edw_research_dm].[luo_khan_nlp_breast_cancer_export]
