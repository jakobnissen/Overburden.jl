# Overburden.jl
Overburden is a command-line tool and Julia package which does depth-aware subsampling of reads.

Specifically, given a sorted BAM file and a minimum depth M, it will produce a new BAM file where all positions P, which had depth D have at least coverage coverage min(D, M). In order words, it will remove reads that are redundant w.r.t. reaching coverage M at any position.
Overburden can be used to reduce the size of read files without removing valuable information, by discarding redundant reads.

## Installation
Make sure you have a recent version of Julia installed (Julia 1.13+).
Install it as a Julia app:

```shell
julia --startup=no -e 'using Pkg; Pkg.Apps.add(url="https://github.com/jakobnissen/Overburden.jl")'
```

You should then be able to run with `overburden`. If not, check that the installation directory
(default: `~/.julia/bin`) is on your `PATH`.

### Experimental: Static compilation:
Try running `sh build.sh` to statically compile `overburden` with JuliaC on Julia 1.13+.
If it works, the result will be an executable binary in `build/bin/overburden`.
The executable must not be moved away from its dynamic libraries that are also created
in the build directory on compilation, as it links to them with relative paths, but you can move the whole build directory.

## Usage
For now, control the number of threads with the env var, e.g. `JULIA_NUM_THREADS=4`. I have omitted that below.

Usage: `overburden [-d INT] [-l INT] INPUT OUTPUT`
Example: `overburden -d 15 ecoli.sorted.bam ecoli.subsampled.bam`

The input must be a BAM file sorted by coordinate.

## Walkthrough
Typically people will want to subsample FASTQ reads based on the reads' alignment to some reference.

1. Map your reads to the reference, e.g. using `strobealign`, then sort by coordinate with `samtools`
    - `strobealign -t 8 ref.fa fw.fq.gz rv.fq.gz | samtools sort -@ 3 -m 3G -o sorted.bam`
2. Subsample with `overburden` (see usage above, command simplified here)
    - `overburden sorted.bam subsampled.bam`
3. Extract the FASTQ from the BAM files
    - `samtools fastq -@ 3 -1 fw.ss.fq.gz -2 rv.ss.fq.gz subsampled.bam`

## Algorithm and implementation details
Overburden makes two passes over the BAM file. The first pass sweeps the alignment blocks of each reference from left to right.
Wherever the reads selected so far cover a position less than min(M, depth) times, it selects more templates (a read and its mate together), preferring those whose alignment reaches furthest to the right.
The second pass writes all records of the selected templates.

This greedy approach is fast and always meets the depth target, but it may keep more reads than necessary.
Finding the smallest such subset is NP-hard.
See [algorithm.md](algorithm.md) for details.
