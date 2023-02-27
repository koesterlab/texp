use std::path::{Path, PathBuf};

use csv;
use anyhow::Result;


use crate::common::Outdir;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::{ProbDistribution1d, self};


pub(crate) fn write_fold_changes(
    preprocessing: &Path,
    differential_expression_path: &Path,
    output: &Path,
)-> Result<()> {

    let in_dir = Outdir::open(&differential_expression_path)?;
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190432).collect();

    // Open a csv file for writing
    let mut wtr = csv::Writer::from_path(output)?;
    wtr.serialize(("feature_id", "foldchange", "probability")).unwrap();

    feature_ids.iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            // println!("feature_id {:?}", feature_id);
            //test if file for feature_id exists at differential_expression_path before deserializing
            //build string for file path and feature_id with ending ".mpk"
            let mut file_path = differential_expression_path.to_path_buf();
            file_path.push(feature_id);
            file_path.set_extension("mpk");
            // println!("file_path {:?}", file_path);
            if !file_path.exists() {
                // println!("file does not exist");
                return Ok(());
            }  
            
            let diff_exp_distribution: ProbDistribution1d = in_dir.deserialize_value(feature_id)?;
            for (fold_change, (prob, _, _)) in diff_exp_distribution.points {
                wtr.serialize((feature_id, fold_change, prob.exp())).unwrap();
            }
            Ok(())
    })?;
    Ok(())
}


pub(crate) fn write_kallisto_fold_changes(
        preprocessing: &Path,
        sample_ids: Vec<String>,
        output: &Path,
    )-> Result<()> {
        // println!("differential_expression_path {:?}", differential_expression_path);
        println!("sample_ids {:?}", sample_ids);
        // let in_dir = Outdir::open(&differential_expression_path)?;
        let preprocessing = Preprocessing::from_path(preprocessing)?;
        let scale_factors = preprocessing.scale_factors();
        let means_disp = preprocessing.mean_disp_estimates();
        let means1 = means_disp.get(&sample_ids[0]).unwrap().means(); 
        let s1 = scale_factors.get(&sample_ids[0]).unwrap();
        let means2 = means_disp.get(&sample_ids[1]).unwrap().means();
        let s2 = scale_factors.get(&sample_ids[1]).unwrap();
        let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190432).collect();
    
        // Open a csv file for writing
        let mut wtr = csv::Writer::from_path(output)?;
        wtr.serialize(("feature_id", "foldchange", "probability")).unwrap();
    
        feature_ids.iter()
            .try_for_each(|(i, feature_id)| -> Result<()> {

                let mean1 = means1.get(*i).unwrap() * s1;
                let mean2 = means2.get(*i).unwrap() * s2;
                let fold_change = mean1 / mean2;
                let disp1 = means_disp.get(&sample_ids[0]).unwrap().dispersions().get(*i).unwrap();
                let disp2 = means_disp.get(&sample_ids[1]).unwrap().dispersions().get(*i).unwrap();
                println!("feature_id {:?}, mean1 {:?} disp1 {:?} mean2 {:?} disp2 {:?}fold_change {:?}", feature_id, mean1, disp1, mean2, disp2, fold_change);
                
                // println!("feature_id {:?}", feature_id);
                //test if file for feature_id exists at differential_expression_path before deserializing
                //build string for file path and feature_id with ending ".mpk"
                // let mut file_path = differential_expression_path.to_path_buf();
                // file_path.push(feature_id);
                // file_path.set_extension("mpk");
                // // println!("file_path {:?}", file_path);
                // if !file_path.exists() {
                //     // println!("file does not exist");
                //     return Ok(());
                // }  
                
                // let diff_exp_distribution: ProbDistribution1d = in_dir.deserialize_value(feature_id)?;
                // for (fold_change, (prob, _, _)) in diff_exp_distribution.points {
                wtr.serialize((feature_id, fold_change, 1.)).unwrap();
                // }
                Ok(())
        })?;
       
    Ok(())
}



pub(crate) fn write_kallisto_counts(
    preprocessing: &Path,
    sample_ids: Vec<String>,
    output: &Path,
)-> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
        let scale_factors = preprocessing.scale_factors();
        let means_disp = preprocessing.mean_disp_estimates();

        let mut wtr = csv::Writer::from_path(output)?;
        wtr.serialize(("sample_id", "feature_id", "normalized_count", "dispersion")).unwrap();
        let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190432).collect();

        for sample_id in sample_ids {
            let means = means_disp.get(&sample_id).unwrap().means();
            let dispersions = means_disp.get(&sample_id).unwrap().dispersions();
            let s = scale_factors.get(&sample_id).unwrap();
            
            feature_ids.iter()
                .try_for_each(|(i, feature_id)| -> Result<()> {
                    let mean = means.get(*i).unwrap() * s;
                    let dispersion = dispersions.get(*i).unwrap();
                    wtr.serialize((sample_id.clone(), feature_id, mean, dispersion)).unwrap();
                    Ok(())
            })?;
        }

    Ok(())
}
    