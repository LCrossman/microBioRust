// src/bin/gen_stubs.rs

use _microbiorust::get_stub_info; 
use pyo3_stub_gen::Result;

//from pyo3-stub-gen github page
fn main() -> Result<()> {
    println!("Gathering Python stub information...");
    
    //fetch the gathered info from your main library
    let info = get_stub_info()?;

    //generate the module structures
    info.generate()?;
    let out_dir = std::path::PathBuf::from("microbiorust");

    // Create the alias stub for the gbk submodule
    let gbk_stub_content = "from . import gbk_to_faa, parse_gbk, gbk_to_fna, gbk_to_ffn, gbk_to_faa_count, gbk_to_gff\n";
    std::fs::write(out_dir.join("gbk.pyi"), gbk_stub_content)?;

    // Create the alias stub for embl (if you have one)
    let embl_stub_content = "from . import embl_to_faa, parse_embl, embl_to_fna, embl_to_ffn, embl_to_gff\n";
    std::fs::write(out_dir.join("embl.pyi"), embl_stub_content)?;
    let align_stub_content = "from . import purge_gaps, subset_msa_alignment, get_consensus\n";
    std::fs::write(out_dir.join("align.pyi"), align_stub_content)?;
    let blast_stub_content = "from . import parse_tabular, parse_xml\n";
    std::fs::write(out_dir.join("blast.pyi"), blast_stub_content)?;
    let seqmetrics_stub_content = "from . import hydrophobicity, amino_counts, amino_percentage\n";
    std::fs::write(out_dir.join("seqmetrics.pyi"), seqmetrics_stub_content)?;
    println!("Submodule alias stubs generated successfully!");
    println!("All stubs generated successfully!");
    Ok(())
}
