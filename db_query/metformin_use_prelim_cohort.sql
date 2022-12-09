use FSM_Analytics;

SELECT DISTINCT bc.*
				, m.medication_name
				, m.generic_name
				, mo.order_start_date_key
				, mo.order_end_date_key
FROM [edw_research_dm].[luo_khan_nlp_breast_cancer_export] bc
LEFT JOIN [NM_BI].[dim].[vw_patient_current] p
	ON bc.patient_ir_id = p.ir_id
LEFT JOIN [NM_BI].[fact].[vw_medication_order] mo
	ON p.patient_key = mo.patient_key
	AND (mo.order_start_datetime < bc.date_of_diagnosis)
	AND YEAR(mo.order_start_datetime) < 2017
LEFT JOIN [NM_BI].[dim].[vw_medication_order_profile] mop
	ON mo.medication_order_profile_key = mop.medication_order_profile_key
	AND (
		mop.order_status NOT IN ('Suspend', 'Order Canceled', 'Voided', 'Canceled', 'On Hold, Med Student', 'Voided With Results', 'Suspended') 
		OR mop.order_status IS NULL
	)
LEFT JOIN [nm_bi].[dim].[vw_medication] m                          --adding medication dimension table
	ON m.medication_key = mo.medication_key
	AND (m.medication_name LIKE '%metformin%' or m.generic_name  LIKE '%metformin%' ) --search on medication and generic names for better results
order by bc.patient_ir_id;


SELECT DISTINCT bc.*
				, m.medication_name
				, m.generic_name
				, mo.order_start_date_key
				, mo.order_end_date_key
FROM [edw_research_dm].[luo_khan_nlp_breast_cancer_export] bc
LEFT JOIN [NM_BI].[dim].[vw_patient_current] p
	ON bc.patient_ir_id = p.ir_id
LEFT JOIN [NM_BI].[fact].[vw_medication_order] mo
	ON p.patient_key = mo.patient_key
	AND (mo.order_start_datetime < bc.date_of_diagnosis)
	AND YEAR(mo.order_start_datetime) < 2017
LEFT JOIN [NM_BI].[dim].[vw_medication_order_profile] mop
	ON mo.medication_order_profile_key = mop.medication_order_profile_key
	AND (
		mop.order_status NOT IN ('Suspend', 'Order Canceled', 'Voided', 'Canceled', 'On Hold, Med Student', 'Voided With Results', 'Suspended') 
		OR mop.order_status IS NULL
	)
LEFT JOIN [nm_bi].[dim].[vw_medication] m                          --adding medication dimension table
	ON m.medication_key = mo.medication_key
	AND (m.medication_name LIKE '%statin%' or m.generic_name  LIKE '%statin%' ) --search on medication and generic names for better results
order by bc.patient_ir_id;

