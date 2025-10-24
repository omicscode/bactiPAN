use glob::glob;
use rayon::prelude::*;
use std::error::Error;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process::Command;

/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() -> io::Result<()> {
    println!("Starting analyzing the pangenome");
    println!("Alignment tools options are macse, prank, muscle, or all");
    let dirpath = read_input("Enter the directory path to the sequences: ")?;
    let species: usize =
        read_input("Enter the number of the species involved in the orthology search: ")?
            .trim()
            .parse()
            .expect("Invalid number of species");
    let threads: usize = read_input("Enter the number of the threads: ")?
        .trim()
        .parse()
        .expect("Invalid number of threads");
    let tool = read_input("Enter the alignment tools: ")?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Failed to initialize rayon thread pool");
    std::env::set_current_dir(&dirpath)?;
    println!("Checking whether all the files are correct");
    println!("Implementing a regular expression for the file search and the number enumeration");

    let mut file_count_correct = true;
    for entry in glob("*.fa")? {
        let path = entry?;
        let content = fs::read_to_string(&path)?;
        let species_count = content.lines().filter(|line| line.starts_with('>')).count();
        if species_count != species {
            file_count_correct = false;
            break;
        }
    }

    if file_count_correct {
        println!(
            "Files are all correct and the number of the species in all the files are present"
        );
    }

    println!("Starting the analysis");

    match tool.trim().to_lowercase().as_str() {
        "macse" => run_macse(&dirpath, threads)?,
        "prank" => run_prank(&dirpath, threads)?,
        "muscle" => run_muscle(&dirpath, threads)?,
        "all" => {
            run_macse(&dirpath, threads)?;
            run_prank(&dirpath, threads)?;
            run_muscle(&dirpath, threads)?;
        }
        _ => println!("Invalid tool specified. Options are: macse, prank, muscle, all"),
    }

    Ok(())
}

fn read_input(prompt: &str) -> io::Result<String> {
    print!("{}", prompt);
    io::stdout().flush()?;
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    Ok(input.trim().to_string())
}

fn run_macse(dirpath: &str, threads: usize) -> io::Result<()> {
    println!("Aligning the pangenome using macse");
    let paths: Vec<_> = glob(&format!("{}/*.fa", dirpath))?.collect::<Result<Vec<_>, _>>()?;
    paths.par_iter().try_for_each(|path| {
        let filename = path.file_name().unwrap().to_str().unwrap();
        let output_aa = format!("{}.AA", filename);
        let output_nt = format!("{}.NT", filename);
        let log_file = format!("{}.macse.run.log.txt", filename);

        Command::new("java")
            .args(&[
                "-jar",
                "-Xmx100g",
                "macse",
                "-prog",
                "alignSequences",
                "-gc_def",
                "12",
                "-seq",
                path.to_str().unwrap(),
                "-out_AA",
                &output_aa,
                "-out_NT",
                &output_nt,
            ])
            .output()
            .and_then(|output| {
                fs::write(&log_file, output.stdout)?;
                Ok(())
            })
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))
    })?;

    let nt_paths: Vec<_> = glob("*.NT")?.collect::<Result<Vec<_>, _>>()?;
    nt_paths.par_iter().try_for_each(|nt_path| {
        let new_name = nt_path.with_extension("ntaligned.fasta");
        fs::rename(nt_path, &new_name)?;
        Ok(())
    })?;

    println!("Renamed the aligned files as ntaligned");

    let aligned_paths: Vec<_> = glob("*.ntaligned.fasta")?.collect::<Result<Vec<_>, _>>()?;
    aligned_paths.par_iter().try_for_each(|aligned_path| {
        let trimmed_name = aligned_path.with_extension("trimmed.fasta");
        Command::new("trimal")
            .args(&[
                "-in",
                aligned_path.to_str().unwrap(),
                "-out",
                trimmed_name.to_str().unwrap(),
                "-nogaps",
            ])
            .output()?;
        Ok(())
    })?;
    Command::new("wget")
        .arg("https://github.com/marekborowiec/AMAS/AMAS.py")
        .output()?;
    Command::new("chmod").args(&["755", "AMAS.py"]).output()?;
    Command::new("python3")
        .args(&["AMAS.py", "-in", "*.fasta", "-f", "fasta", "-d", "dna"])
        .output()?;
    fs::rename("concatenated.out", "macsealignmentconcatenated.fasta")?;
    fs::rename("partitions.txt", "macsealignmentpartitions.txt")?;
    Command::new("iqtree")
        .args(&[
            "--seqtype",
            "DNA",
            "-s",
            "macsealignmentconcatenated.fasta",
            "--alrt",
            "1000",
            "-b",
            "1000",
            "-T",
            &threads.to_string(),
        ])
        .output()?;
    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "macsealignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "macsephylogeny_GAMMA",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
        ])
        .output()?;

    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "macsealignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "macsephylogeny_GTRCAT",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
            "-b",
            "1000",
        ])
        .output()?;

    println!("Finishing up the analysis");
    Ok(())
}

fn run_prank(dirpath: &str, threads: usize) -> io::Result<()> {
    println!("Aligning the pangenome using prank probabilistic alignment");
    let paths: Vec<_> = glob(&format!("{}/*.fa", dirpath))?.collect::<Result<Vec<_>, _>>()?;
    paths.par_iter().try_for_each(|path| {
        let output_name = path.with_extension("prankaligned.fasta");
        Command::new("prank")
            .args(&[
                "-d",
                path.to_str().unwrap(),
                "-o",
                output_name.to_str().unwrap(),
            ])
            .output()?;
        Ok(())
    })?;
    let aligned_paths: Vec<_> = glob("*.prankaligned.fasta")?.collect::<Result<Vec<_>, _>>()?;
    aligned_paths.par_iter().try_for_each(|aligned_path| {
        let trimmed_name = aligned_path.with_extension("pranktrimmed.fasta");
        Command::new("trimal")
            .args(&[
                "-in",
                aligned_path.to_str().unwrap(),
                "-out",
                trimmed_name.to_str().unwrap(),
                "-nogaps",
            ])
            .output()?;
        Ok(())
    })?;
    Command::new("wget")
        .arg("https://github.com/marekborowiec/AMAS/AMAS.py")
        .output()?;
    Command::new("chmod").args(&["755", "AMAS.py"]).output()?;
    Command::new("python3")
        .args(&["AMAS.py", "-in", "*.fasta", "-f", "fasta", "-d", "dna"])
        .output()?;

    fs::rename("concatenated.out", "prankalignmentconcatenated.fasta")?;
    fs::rename("partitions.txt", "prankalignmentpartitions.txt")?;
    Command::new("iqtree")
        .args(&[
            "--seqtype",
            "DNA",
            "-s",
            "prankalignmentconcatenated.fasta",
            "--alrt",
            "1000",
            "-b",
            "1000",
            "-T",
            &threads.to_string(),
        ])
        .output()?;
    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "prankalignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "prankphylogeny_GAMMA",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
        ])
        .output()?;
    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "prankalignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "prankphylogeny_GTRCAT",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
            "-b",
            "1000",
        ])
        .output()?;

    Ok(())
}

fn run_muscle(dirpath: &str, threads: usize) -> io::Result<()> {
    println!("Aligning the pangenome using muscle alignment");
    let paths: Vec<_> = glob(&format!("{}/*.fa", dirpath))?.collect::<Result<Vec<_>, _>>()?;
    paths.par_iter().try_for_each(|path| {
        let output_name = path.with_extension("musclealigned.fasta");
        Command::new("muscle")
            .args(&[
                "-in",
                path.to_str().unwrap(),
                "-out",
                output_name.to_str().unwrap(),
            ])
            .output()?;
        Ok(())
    })?;
    let aligned_paths: Vec<_> = glob("*.musclealigned.fasta")?.collect::<Result<Vec<_>, _>>()?;
    aligned_paths.par_iter().try_for_each(|aligned_path| {
        let trimmed_name = aligned_path.with_extension("muscletrimmed.fasta");
        Command::new("trimal")
            .args(&[
                "-in",
                aligned_path.to_str().unwrap(),
                "-out",
                trimmed_name.to_str().unwrap(),
                "-nogaps",
            ])
            .output()?;
        Ok(())
    })?;
    Command::new("wget")
        .arg("https://github.com/marekborowiec/AMAS/AMAS.py")
        .output()?;
    Command::new("chmod").args(&["755", "AMAS.py"]).output()?;
    Command::new("python3")
        .args(&["AMAS.py", "-in", "*.fasta", "-f", "fasta", "-d", "dna"])
        .output()?;
    fs::rename("concatenated.out", "musclealignmentconcatenated.fasta")?;
    fs::rename("partitions.txt", "musclealignmentpartitions.txt")?;
    Command::new("iqtree")
        .args(&[
            "--seqtype",
            "DNA",
            "-s",
            "musclealignmentconcatenated.fasta",
            "--alrt",
            "1000",
            "-b",
            "1000",
            "-T",
            &threads.to_string(),
        ])
        .output()?;
    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "musclealignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "musclephylogeny_GAMMA",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
        ])
        .output()?;
    Command::new("raxmlHPC-PTHREADS")
        .args(&[
            "-s",
            "musclealignmentconcatenated.fasta",
            "--no-seq-check",
            "-O",
            "-m",
            "GTRGAMMA",
            "-p",
            "12345",
            "-n",
            "musclephylogeny_GTRCAT",
            "-T",
            &threads.to_string(),
            "-N",
            "50",
            "-b",
            "1000",
        ])
        .output()?;
    Ok(())
}
