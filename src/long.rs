use std::env;
use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::Path;
use std::process::Command;

/*
Gaurav Sablok
codeprog@icloud.com
*/

pub fn longread() -> io::Result<()> {
    println!("Generating or updating the genome assembly from long reads such as PacBio");
    println!(
        "If you have an already assembled genome, use the update option or else use the genome assembly option"
    );
    println!(
        "Requires: canu, bamtools, bowtie2, blasr, lastz, flye, mecatref, pilon, jasper, and busco"
    );
    let checkavail = prompt("Do you have all these tools installed? (yes/no): ")?;
    let tool_paths = if checkavail.to_lowercase() == "yes" {
        println!("Please provide paths to the tools:");
        let paths = vec![
            prompt("Path to canu: ")?,
            prompt("Path to blasr: ")?,
            prompt("Path to lastz: ")?,
            prompt("Path to flye: ")?,
            prompt("Path to mecatref: ")?,
            prompt("Path to bamtools: ")?,
            prompt("Path to bowtie2: ")?,
            prompt("Path to pilon: ")?,
        ];
        for (i, path) in paths.iter().enumerate() {
            println!("Provided path {}: {}", i + 1, path);
        }
        paths
    } else {
        println!("Conda environment creation is not directly supported in Rust.");
        println!("Please ensure the tools are installed manually or via conda before proceeding.");
        vec![]
    };

    for path in &tool_paths {
        if !path.is_empty() {
            let current_path = env::var("PATH").unwrap_or_default();
            env::set_var("PATH", format!("{}:{}", path, current_path));
        }
    }
    let species = prompt("Please provide the species name: ")?;
    let first = prompt("Are you assembling the genome for the first time? (yes/no): ")?;
    let (reads, read_removal) = if first.to_lowercase() == "yes" {
        let reads = prompt("Please provide the path for the genome assembly reads: ")?;
        let removal = prompt(
            "Do you have a reference genome or a contaminant file for read removal? (yes/no): ",
        )?;
        let read_removal = if removal.to_lowercase() == "yes" {
            prompt("Please provide the reference genome or contaminant file: ")?
        } else {
            String::new()
        };
        (reads, read_removal)
    } else {
        (String::new(), String::new())
    };

    let mapping = prompt("How many reads do you want to keep while blasr mapping: ")?;
    let coverage = prompt("Do you want to calculate the coverage also? (yes/no): ")?;
    let polish = prompt("Do you want to polish the genome also? (yes/no): ")?;
    let polish_option = if polish.to_lowercase() == "yes" {
        prompt("Which polishing method do you want to use? (pilon/jasper): ")?
    } else {
        String::new()
    };

    let (previous_assembly, new_reads_directory) = if first.to_lowercase() != "yes" {
        let update = prompt("Are you updating an existing assembly? (yes/no): ")?;
        if update.to_lowercase() == "yes" {
            (
                prompt("Please provide the path for the previous assembly: ")?,
                prompt("Please provide the path for the directory containing the new reads: ")?,
            )
        } else {
            (String::new(), String::new())
        }
    } else {
        (String::new(), String::new())
    };

    let illumina_reads = if coverage.to_lowercase() == "yes" {
        prompt("Please provide the directory containing the Illumina reads: ")?
    } else {
        String::new()
    };

    let assembler = prompt("Please select the choice of the assembler (canu/flye/mecatref): ")?;
    let (threads, calibrate, sensitivity, read_selection, meryl_memory, flye_overlap, ref_genome) =
        match assembler.as_str() {
            "canu" => (
                prompt("Please provide the threads you want to use: ")?,
                prompt("Please provide the genome size for calibration: ")?,
                prompt("Please provide the sensitivity for the genome assembly: ")?,
                prompt("Please provide the minimum length for read selection: ")?,
                prompt("Please provide the memory bindings for meryl: ")?,
                String::new(),
                String::new(),
            ),
            "flye" => (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                prompt("Please provide the minimum overlap for the long reads: ")?,
                String::new(),
            ),
            "mecatref" => (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                prompt("Please provide the name of the reference genome FASTA file: ")?,
            ),
            _ => (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            ),
        };

    let alignment_tracks =
        prompt("Do you want to make genome alignments for comparison? (yes/no): ")?;
    let (comparison, comparison_gff, threshold, align_coverage, file_format) =
        if alignment_tracks.to_lowercase() == "yes" {
            (
                prompt("Please provide the reference genome for comparison: ")?,
                prompt("Please provide the genome GFF files for comparison: ")?,
                prompt("Please provide the threshold sensitivity (default 70% if empty): ")
                    .unwrap_or("70".to_string()),
                prompt("Please provide the coverage for alignments (default 60% if empty): ")
                    .unwrap_or("60".to_string()),
                prompt("Which format of alignments do you want? (sam/general): ")
                    .unwrap_or("sam".to_string()),
            )
        } else {
            (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            )
        };

    println!("Thank you for choosing the workflow options.");
    println!("Using working directory: {}", env::current_dir()?.display());
    if !species.is_empty()
        && first.to_lowercase() == "yes"
        && ["canu", "flye", "mecatref"].contains(&assembler.as_str())
        && alignment_tracks.to_lowercase() == "yes"
        && !mapping.is_empty()
        && coverage.to_lowercase() == "yes"
        && polish.to_lowercase() == "yes"
        && !polish_option.is_empty()
    {
        process_new_assembly(
            &species,
            &reads,
            &read_removal,
            &mapping,
            &assembler,
            &threads,
            &calibrate,
            &sensitivity,
            &read_selection,
            &meryl_memory,
            &flye_overlap,
            &ref_genome,
            &illumina_reads,
            &comparison,
            &comparison_gff,
            &threshold,
            &align_coverage,
            &file_format,
            &polish_option,
        )?;
    } else if !species.is_empty()
        && first.to_lowercase() == "no"
        && !previous_assembly.is_empty()
        && !new_reads_directory.is_empty()
        && ["flye", "mecatref"].contains(&assembler.as_str())
        && alignment_tracks.to_lowercase() == "yes"
        && !mapping.is_empty()
        && coverage.to_lowercase() == "yes"
        && polish.to_lowercase() == "yes"
        && !polish_option.is_empty()
    {
        process_updated_assembly(
            &species,
            &previous_assembly,
            &new_reads_directory,
            &read_removal,
            &mapping,
            &assembler,
            &threads,
            &calibrate,
            &sensitivity,
            &read_selection,
            &meryl_memory,
            &flye_overlap,
            &ref_genome,
            &illumina_reads,
            &comparison,
            &comparison_gff,
            &threshold,
            &align_coverage,
            &file_format,
            &polish_option,
        )?;
    } else {
        println!("You haven't selected enough options to run this workflow.");
    }

    println!("Processing finished.");
    Ok(())
}

fn prompt(question: &str) -> io::Result<String> {
    print!("{}", question);
    io::stdout().flush()?;
    let mut input = String::new();
    io::stdin().read_line(&mut input)?;
    Ok(input.trim().to_string())
}

fn process_new_assembly(
    species: &str,
    reads: &str,
    read_removal: &str,
    mapping: &str,
    assembler: &str,
    threads: &str,
    calibrate: &str,
    sensitivity: &str,
    read_selection: &str,
    meryl_memory: &str,
    flye_overlap: &str,
    ref_genome: &str,
    illumina_reads: &str,
    comparison: &str,
    comparison_gff: &str,
    threshold: &str,
    align_coverage: &str,
    file_format: &str,
    polish_option: &str,
) -> io::Result<()> {
    // Create directories
    let dirs = vec![
        "reads_assembly",
        "genome_assembly",
        "illumina_reads",
        "polished_genome",
        "sam_file_illumina",
    ];
    for dir in &dirs {
        fs::create_dir_all(dir)?;
    }

    let longreads = format!("{}/reads_assembly", env::current_dir()?.display());
    let genome_assembly = format!("{}/genome_assembly", env::current_dir()?.display());
    let illumina_dir = format!("{}/illumina_reads", env::current_dir()?.display());
    let polished_genome = format!("{}/polished_genome", env::current_dir()?.display());
    let coverage_dir = format!("{}/sam_file_illumina", env::current_dir()?.display());
    process_reads(
        reads,
        &longreads,
        read_removal,
        mapping,
        species,
        &genome_assembly,
    )?;
    process_illumina_reads(illumina_reads, &illumina_dir)?;
    let genome_assembly_fasta = run_assembler(
        assembler,
        species,
        &genome_assembly,
        &longreads,
        threads,
        calibrate,
        sensitivity,
        read_selection,
        meryl_memory,
        flye_overlap,
        ref_genome,
    )?;

    run_alignments(
        &genome_assembly_fasta,
        comparison,
        comparison_gff,
        threshold,
        align_coverage,
        file_format,
        species,
        &genome_assembly,
    )?;

    process_coverage(
        &genome_assembly_fasta,
        &illumina_dir,
        species,
        threads,
        &polished_genome,
        &coverage_dir,
    )?;

    polish_genome(
        &genome_assembly_fasta,
        &polished_genome,
        species,
        polish_option,
        &illumina_dir,
        threads,
    )?;

    Ok(())
}

fn process_updated_assembly(
    species: &str,
    previous_assembly: &str,
    new_reads_directory: &str,
    read_removal: &str,
    mapping: &str,
    assembler: &str,
    threads: &str,
    calibrate: &str,
    sensitivity: &str,
    read_selection: &str,
    meryl_memory: &str,
    flye_overlap: &str,
    ref_genome: &str,
    illumina_reads: &str,
    comparison: &str,
    comparison_gff: &str,
    threshold: &str,
    align_coverage: &str,
    file_format: &str,
    polish_option: &str,
) -> io::Result<()> {
    // Create directories
    let dirs = vec![
        "new_genomic_reads",
        "combined_assembly",
        "illumina_reads",
        "polished_genome",
        "sam_file_illumina",
    ];
    for dir in &dirs {
        fs::create_dir_all(dir)?;
    }

    let new_longreads = format!("{}/new_genomic_reads", env::current_dir()?.display());
    let combined_assembly = format!("{}/combined_assembly", env::current_dir()?.display());
    let illumina_dir = format!("{}/illumina_reads", env::current_dir()?.display());
    let polished_genome = format!("{}/polished_genome", env::current_dir()?.display());
    let coverage_dir = format!("{}/sam_file_illumina", env::current_dir()?.display());

    process_updated_reads(
        previous_assembly,
        new_reads_directory,
        &new_longreads,
        &combined_assembly,
    )?;
    process_illumina_reads(illumina_reads, &illumina_dir)?;

    let unaligned_reads = format!("{}/{}.unaligned.fasta", combined_assembly, species);
    Command::new("blasr")
        .args(&[
            "all_reads.fasta",
            read_removal,
            "--best",
            mapping,
            "--bam",
            &format!("{}/{}.bam", combined_assembly, species),
            "--unaligned",
            &unaligned_reads,
        ])
        .current_dir(&combined_assembly)
        .status()?;

    println!("Finished cleaning of the reads");

    let genome_assembly_fasta = run_assembler(
        assembler,
        species,
        &combined_assembly,
        &unaligned_reads,
        threads,
        calibrate,
        sensitivity,
        read_selection,
        meryl_memory,
        flye_overlap,
        ref_genome,
    )?;

    run_alignments(
        &genome_assembly_fasta,
        comparison,
        comparison_gff,
        threshold,
        align_coverage,
        file_format,
        species,
        &combined_assembly,
    )?;

    // Process coverage
    process_coverage(
        &genome_assembly_fasta,
        &illumina_dir,
        species,
        threads,
        &polished_genome,
        &coverage_dir,
    )?;

    polish_genome(
        &genome_assembly_fasta,
        &polished_genome,
        species,
        polish_option,
        &illumina_dir,
        threads,
    )?;

    Ok(())
}

fn process_reads(
    reads: &str,
    longreads: &str,
    read_removal: &str,
    mapping: &str,
    species: &str,
    genome_assembly: &str,
) -> io::Result<()> {
    for entry in fs::read_dir(reads)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "gz") {
            Command::new("gzip").arg(entry.path()).status()?;
        }
    }
    Command::new("sh")
        .arg("-c")
        .arg("cat *.fasta > all_reads.fasta")
        .current_dir(reads)
        .status()?;
    fs::rename(
        format!("{}/all_reads.fasta", reads),
        format!("{}/all_reads.fasta", longreads),
    )?;

    if !read_removal.is_empty() {
        Command::new("blasr")
            .args(&[
                "all_reads.fasta",
                read_removal,
                "--best",
                mapping,
                "--bam",
                &format!("{}/{}.bam", longreads, species),
                "--unaligned",
                &format!("{}/{}.unaligned.fasta", longreads, species),
            ])
            .current_dir(longreads)
            .status()?;
        println!("Finished cleaning of the reads");
        fs::copy(
            format!("{}/{}.unaligned.fasta", longreads, species),
            format!("{}/{}.unaligned.fasta", genome_assembly, species),
        )?;
    }
    Ok(())
}

fn process_illumina_reads(illumina_reads: &str, illumina_dir: &str) -> io::Result<()> {
    for entry in fs::read_dir(illumina_reads)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "gz") {
            Command::new("gunzip").arg(entry.path()).status()?;
        }
    }

    Command::new("sh")
        .arg("-c")
        .arg("cat *.R1.fastq > final.R1.fastq")
        .current_dir(illumina_reads)
        .status()?;
    Command::new("sh")
        .arg("-c")
        .arg("cat *.R2.fastq > final.R2.fastq")
        .current_dir(illumina_reads)
        .status()?;

    let mut r1_file = File::create(format!("{}/illumina.R1.txt", illumina_dir))?;
    let mut r2_file = File::create(format!("{}/illumina.R2.txt", illumina_dir))?;
    for entry in fs::read_dir(illumina_reads)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(false, |ext| ext == "fastq") {
            if path.to_string_lossy().contains(".R1.") {
                writeln!(r1_file, "{}", path.display())?;
            } else if path.to_string_lossy().contains(".R2.") {
                writeln!(r2_file, "{}", path.display())?;
            }
        }
    }
    let mut r1_final = File::create(format!("{}/illumina.R1.final.txt", illumina_dir))?;
    writeln!(r1_final, "final.R1.fastq")?;
    let mut r2_final = File::create(format!("{}/illumina.R1.final.txt", illumina_dir))?;
    writeln!(r2_final, "final.R2.fastq")?;
    Ok(())
}

fn process_updated_reads(
    previous_assembly: &str,
    new_reads_directory: &str,
    new_longreads: &str,
    combined_assembly: &str,
) -> io::Result<()> {
    for entry in fs::read_dir(previous_assembly)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "fasta") {
            fs::copy(
                entry.path(),
                format!(
                    "{}/{}",
                    combined_assembly,
                    entry.path().file_name().unwrap().to_string_lossy()
                ),
            )?;
        }
    }
    for entry in fs::read_dir(new_reads_directory)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "gz") {
            fs::copy(
                entry.path(),
                format!(
                    "{}/{}",
                    new_longreads,
                    entry.path().file_name().unwrap().to_string_lossy()
                ),
            )?;
        }
    }
    Command::new("tar")
        .args(&["zxvf", "*.gz"])
        .current_dir(new_longreads)
        .status()?;
    for entry in fs::read_dir(new_longreads)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "fasta") {
            fs::copy(
                entry.path(),
                format!(
                    "{}/{}",
                    combined_assembly,
                    entry.path().file_name().unwrap().to_string_lossy()
                ),
            )?;
        }
    }
    Command::new("sh")
        .arg("-c")
        .arg("cat *.fasta > update_genome_assembly.fasta")
        .current_dir(combined_assembly)
        .status()?;
    Ok(())
}

fn run_assembler(
    assembler: &str,
    species: &str,
    genome_assembly: &str,
    longreads: &str,
    threads: &str,
    calibrate: &str,
    sensitivity: &str,
    read_selection: &str,
    meryl_memory: &str,
    flye_overlap: &str,
    ref_genome: &str,
) -> io::Result<String> {
    let unaligned_reads = format!("{}/{}.unaligned.fasta", longreads, species);
    let genome_assembly_fasta = format!("{}/{}.fasta", genome_assembly, species);

    match assembler {
        "canu" => {
            Command::new("canu")
                .args(&[
                    "gridOptions=--time=24:00:00",
                    &format!("-corMhapSensitivity={}", sensitivity),
                    &format!("-p={}", species),
                    &format!("-d={}", genome_assembly),
                    &format!("-genomeSize={}", calibrate),
                    &format!("-minReadLength={}", read_selection),
                    &format!("-merylMemory={}", meryl_memory),
                    "-gnuplotImageFormat=png",
                    &format!("-ovsThreads={}", threads),
                    &format!("-ovbThreads={}", threads),
                    &format!("-pacbio-raw={}", unaligned_reads),
                ])
                .current_dir(genome_assembly)
                .status()?;
        }
        "flye" => {
            Command::new("flye")
                .args(&[
                    "--pacbio-raw",
                    &unaligned_reads,
                    "--genome-size",
                    calibrate,
                    "--threads",
                    threads,
                    "--out-dir",
                    genome_assembly,
                    "--min-overlap",
                    flye_overlap,
                ])
                .current_dir(genome_assembly)
                .status()?;
        }
        "mecatref" => {
            Command::new("mecat2ref")
                .args(&[
                    "-d",
                    &unaligned_reads,
                    "-r",
                    ref_genome,
                    "-o",
                    &format!("{}/{}.reference.sam", genome_assembly, species),
                    "-w",
                    &format!("{}/{}_intermediate", genome_assembly, species),
                    "-t",
                    "24",
                    "-n",
                    "20",
                    "-n",
                    "50",
                    "-m",
                    "2",
                    "-x",
                    "0",
                ])
                .current_dir(genome_assembly)
                .status()?;
        }
        _ => println!("Unsupported assembler: {}", assembler),
    }

    for entry in fs::read_dir(genome_assembly)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "fasta") {
            return Ok(entry.path().to_string_lossy().to_string());
        }
    }
    Ok(genome_assembly_fasta)
}

fn run_alignments(
    genome_assembly_fasta: &str,
    comparison: &str,
    comparison_gff: &str,
    threshold: &str,
    align_coverage: &str,
    file_format: &str,
    species: &str,
    genome_assembly: &str,
) -> io::Result<()> {
    let output_sam = format!("{}/{}.refalignments.sam", genome_assembly, species);
    if threshold.is_empty() && align_coverage.is_empty() && file_format.is_empty() {
        Command::new("lastz")
            .args(&[
                &format!("{}[multiple]", comparison),
                genome_assembly_fasta,
                "--identity=75",
                "--coverage=60",
                "--format=sam",
                "--ambiguous=iupac",
                "--output",
                &output_sam,
            ])
            .current_dir(genome_assembly)
            .status()?;
    } else {
        Command::new("lastz")
            .args(&[
                &format!("{}[multiple]", comparison),
                genome_assembly_fasta,
                &format!("--identity={}", threshold),
                &format!("--coverage={}", align_coverage),
                &format!("--format={}", file_format),
                "--ambiguous=iupac",
                "--output",
                &output_sam,
            ])
            .current_dir(genome_assembly)
            .status()?;
    }
    println!("Genome alignments have been written");

    Command::new("quast")
        .args(&[
            genome_assembly_fasta,
            "-r",
            comparison,
            "-g",
            comparison_gff,
            "--fragmented",
            "-o",
            &format!(
                "{}/{}.quast.report",
                genome_assembly,
                Path::new(genome_assembly_fasta)
                    .file_stem()
                    .unwrap()
                    .to_string_lossy()
            ),
            "-m",
            "500",
            "-t",
            "40",
            "-e",
        ])
        .current_dir(genome_assembly)
        .status()?;
    println!("Performing genome completeness analysis");
    Ok(())
}

fn process_coverage(
    genome_assembly_fasta: &str,
    illumina_dir: &str,
    species: &str,
    threads: &str,
    polished_genome: &str,
    coverage_dir: &str,
) -> io::Result<()> {
    fs::copy(
        genome_assembly_fasta,
        format!(
            "{}/{}",
            illumina_dir,
            Path::new(genome_assembly_fasta)
                .file_name()
                .unwrap()
                .to_string_lossy()
        ),
    )?;

    Command::new("bowtie2-build")
        .args(&[genome_assembly_fasta, species])
        .current_dir(illumina_dir)
        .status()?;

    let r1_file = format!("{}/illumina.R1.txt", illumina_dir);
    let r2_file = format!("{}/illumina.R2.txt", illumina_dir);
    let r1_lines: Vec<String> = fs::read_to_string(&r1_file)?
        .lines()
        .map(String::from)
        .collect();
    let r2_lines: Vec<String> = fs::read_to_string(&r2_file)?
        .lines()
        .map(String::from)
        .collect();

    for (col1, col2) in r1_lines.iter().zip(r2_lines.iter()) {
        Command::new("bowtie2")
            .args(&[
                "-t",
                "-x",
                species,
                "-p",
                threads,
                "--very-sensitive-local",
                "-1",
                col1,
                "-2",
                col2,
                "-S",
                &format!("{}/{}.sam", illumina_dir, species),
                "--no-unal",
                "--al-conc",
                &format!("{}/{}.aligned.fastq", illumina_dir, species),
            ])
            .current_dir(illumina_dir)
            .status()?;
    }

    fs::copy(
        format!("{}/{}.sam", illumina_dir, species),
        format!("{}/{}.sam", polished_genome, species),
    )?;

    for entry in fs::read_dir(illumina_dir)? {
        let entry = entry?;
        if entry.path().extension().map_or(false, |ext| ext == "sam") {
            let f = entry.path();
            let f_bam = format!("{}.bam", f.to_string_lossy());
            Command::new("samtools")
                .args(&["view", "-b", &f.to_string_lossy(), "-o", &f_bam])
                .current_dir(coverage_dir)
                .status()?;
            Command::new("bamtools")
                .args(&[
                    "coverage",
                    "-in",
                    &f_bam,
                    "-out",
                    &format!("{}/{}.bam.coverage.txt", coverage_dir, species),
                ])
                .current_dir(coverage_dir)
                .status()?;
            Command::new("bamtools")
                .args(&[
                    "stats",
                    "-in",
                    &f_bam,
                    "-insert",
                    ">",
                    &format!(
                        "{}/{}.bam.stats.coverage.illumina.txt",
                        coverage_dir, species
                    ),
                ])
                .current_dir(coverage_dir)
                .status()?;
            Command::new("bamtools")
                .args(&[
                    "sort",
                    "-in",
                    &f_bam,
                    "-out",
                    &format!("{}.sorted.bam", f.to_string_lossy()),
                    "-n",
                    "500000",
                    "-mem",
                    "1024",
                ])
                .current_dir(coverage_dir)
                .status()?;
            fs::copy(
                format!("{}.sorted.bam", f.to_string_lossy()),
                format!(
                    "{}/{}.sorted.bam",
                    polished_genome,
                    Path::new(&f_bam).file_stem().unwrap().to_string_lossy()
                ),
            )?;
        }
    }
    Ok(())
}

fn polish_genome(
    genome_assembly_fasta: &str,
    polished_genome: &str,
    species: &str,
    polish_option: &str,
    illumina_dir: &str,
    threads: &str,
) -> io::Result<()> {
    fs::copy(
        genome_assembly_fasta,
        format!(
            "{}/{}",
            polished_genome,
            Path::new(genome_assembly_fasta)
                .file_name()
                .unwrap()
                .to_string_lossy()
        ),
    )?;

    if polish_option == "jasper" {
        let jasper_url =
            "https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz";
        let jasper_tar = format!("{}/jasper-1.0.3.tar.gz", polished_genome);
        if !Path::new(&jasper_tar).exists() {
            println!("Downloading JASPER...");
            Command::new("wget")
                .args(&[jasper_url, "-O", &jasper_tar])
                .current_dir(polished_genome)
                .status()?;
            if !Path::new(&jasper_tar).exists() {
                println!("JASPER download failed, retrying...");
                Command::new("wget")
                    .args(&[jasper_url, "-O", &jasper_tar])
                    .current_dir(polished_genome)
                    .status()?;
            }
        }
        Command::new("tar")
            .args(&["zxvf", "jasper-1.0.3.tar.gz"])
            .current_dir(polished_genome)
            .status()?;
        Command::new("./configure")
            .arg(&format!("--prefix={}/jasper-1.0.3", polished_genome))
            .current_dir(format!("{}/jasper-1.0.3", polished_genome))
            .status()?;
        Command::new("make")
            .arg("install")
            .current_dir(format!("{}/jasper-1.0.3", polished_genome))
            .status()?;
        env::set_var("PYTHONPATH", format!("{}/jasper-1.0.3", polished_genome));
        env::set_var(
            "PATH",
            format!(
                "{}/jasper-1.0.3:{}",
                polished_genome,
                env::var("PATH").unwrap_or_default()
            ),
        );

        Command::new("jasper.sh")
            .args(&[
                "-t",
                threads,
                genome_assembly_fasta,
                "-r",
                &format!("{}/final.R1.fastq", illumina_dir),
                &format!("{}/final.R1.fastq", illumina_dir), // Note: R1 used twice as per original script
            ])
            .current_dir(polished_genome)
            .status()?;
    } else {
        Command::new("java")
            .args(&[
                "-Xmx100g",
                "-jar",
                "pilon-1.28.jar",
                "--genome",
                genome_assembly_fasta,
                "--frags",
                &format!("{}/{}.sorted.bam", polished_genome, species),
                "--output",
                &format!("{}/{}.polished_genome", polished_genome, species),
                "--threads",
                "60",
                "--changes",
                "--fix",
                "all",
                "--mindepth",
                "10",
            ])
            .current_dir(polished_genome)
            .status()?;
    }
    Ok(())
}
