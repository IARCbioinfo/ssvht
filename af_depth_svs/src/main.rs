use rust_htslib::{bam, bam::Read};
use bio::io::bed;
use std::env;


fn main() {

    let args: Vec<String> = env::args().collect();

    let bam = &args[1];
    let bed_file = &args[2];
    
    println!("bam file {}", bam);
    println!("bed file {}", bed_file);

//go to specific site
let mut bam_index = bam::IndexedReader::from_path(&bam).unwrap();
let mut reader = bed::Reader::from_file(&bed_file).unwrap();

for record in reader.records() {
//we fetch any record 
	let rec = record.ok().expect("Error reading record.");
	let _fbam = bam_index.fetch((rec.chrom(), rec.start(), rec.end())); // coordinates 10000..20000 on reference named "chrX"
	
	let mut qv = 0;
	let mut bases=0;
	for p in bam_index.pileup() {
		let pileup = p.unwrap();	
		if u64::from(pileup.pos()) >= rec.start() && u64::from(pileup.pos()) <= rec.end() {
			for alignment in pileup.alignments() {
        			if !alignment.is_del() && !alignment.is_refskip() && alignment.record().mapq() > 20 {
            			//println!("Base {} ALNQ {}", alignment.record().seq()[alignment.qpos().unwrap()], alignment.record().mapq());
					qv+=1;
        			}
			}	
    		//println!("{}:{} depth {} MQ {}", pileup.tid(), pileup.pos(), pileup.depth(), qv);
   		//qv=0;
       		 bases+=1;	
		}
	}
    	println!("{}:{} {} {:?} {} {} {}", rec.chrom(),rec.start(),rec.end(),rec.name(), qv/bases, qv, bases);

}

}


/*
let bam = bam::Reader::from_path(&"test/test.bam").unwrap();
let header = bam::Header::from_template(bam.header());

// print header records to the terminal, akin to samtool
for (key, records) in header.to_hashmap() {
    for record in records {
	if key== "SQ" {
         println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
	}
    }
}
*/
/*let mut qv = 0;
for p in bam.pileup() {
    let pileup = p.unwrap();
	for alignment in pileup.alignments() {
        	if !alignment.is_del() && !alignment.is_refskip() &&  alignment.record().mapq() > 40 {
            		//println!("Base {} ALNQ {}", alignment.record().seq()[alignment.qpos().unwrap()], alignment.record().mapq());
			qv+=1;
        	}
	}	
    println!("{}:{} depth {} MQ {}", pileup.tid(), pileup.pos(), pileup.depth(), qv);
   qv=0;
}
*/

