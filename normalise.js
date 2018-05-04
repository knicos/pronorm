#!/usr/bin/node
/**
 * Progenesis version (normalise-to-all).
 *
 * - Load data from CSV, filter out rows with missing data
 * - For each combination of runs in the two sets:
 *   - Calculate the log ratios between them
 *   - Iteratively filter outliers using median and MAD, 3 x deviation
 *   - Calc scale factor in abundance space as 2^mean (mean in log space)
 *   - Scale corresponding abundance set by factor 
 * - Group peptides by protein
 * - For each combination of runs in the two protein sets:
 *   - Calculate log 2 fold change for proteins
 * - Calculate average fold across runs/normalisations for each protein
 * - Print results as CSV
 */

const fs = require('fs');
let yargs = require('yargs');
const ss = require('simple-statistics');

let argv = yargs
.describe("in", "Input data file")
.describe("out", "Output file")
.describe("set1", "Data columns to use, eg. '2,3'")
.describe("set2", "Data columns to use, eg. '2,3'")
.describe("info-cols", "Info columns for export")
.describe("gene-col", "Gene label column")
//.describe("dump","Output all normalisation errors")
//.describe("sort", "Sort by 'error' or 'gradient'")
.describe("group", "Group by protein")
//.describe("means", "Print final fold change means")
.describe("median", "Use median instead of mean")
.describe("noskip", "Include missing data")
//.describe("log2", "Use log base 2 of normalisation")
//.describe("log10", "Use log base 10 of normalisation")
//.describe("dump-input", "Process CSV only")
.demandOption("in", "Must specify an input file")
.demandOption("set1", "Missing first set of runs")
.demandOption("set2", "Missing second set of runs")
.demandOption("gene-col", "Missing protein or gene column number")
.boolean(["group","median", "noskip"])
.argv;

const infile = argv["in"];
const outfile = argv["out"];

function parseCSV(data) {
	let out = {
		size: 0,
		labels: [],
		A: [],
		B: []
	};
	let indata = data;
	let inlines = indata.split("\n");

	let set1 = (typeof argv.set1 == "string") ? argv.set1.split(",") : [argv.set1];
	for (var i=0; i<set1.length; i++) {
		set1[i] = parseInt(set1[i])-1;
		out.A.push([]);
	}
	let set2 = (typeof argv.set2 == "string") ? argv.set2.split(",") : [argv.set2];
	for (var i=0; i<set2.length; i++) {
		set2[i] = parseInt(set2[i])-1;
		out.B.push([]);
	}

	//let set1 = parseInt(colsplit[0])-1;
	//let col2 = parseInt(colsplit[1])-1;
	const genecol = parseInt(argv["gene-col"])-1;
	const infocol = parseInt(argv["info-cols"])-1;

	// Extract the data
	for (var i=1; i<inlines.length; i++) {
		let skip = false;

		let line = inlines[i].split(",");
		let gene = line[genecol];
		if (!isNaN(infocol)) gene += ","+line[infocol];
		let A = []; //parseFloat(line[col1]);
		let B = []; //parseFloat(line[col2]);
		for (var j=0; j<set1.length; j++) {
			let a = parseFloat(line[set1[j]]);
			if (a == 0 || isNaN(a)) {
				skip = true;
				a = 1;
			}
			A.push(a);
		}
		for (var j=0; j<set2.length; j++) {
			let a = parseFloat(line[set2[j]]);
			if (a == 0 || isNaN(a)) {
				skip = true;
				a = 1;
			}
			B.push(a);
		}

		if (!argv.noskip && skip == true) continue;

		for (var j=0; j<set1.length; j++) out.A[j].push(A[j]);
		for (var j=0; j<set2.length; j++) out.B[j].push(B[j]);

		// Skip missing data.
		//if(isNaN(A)) console.log("NaN: ", line[col1]); 
		//if(isNaN(B)) console.log("NaN: ", line[col2]); 
		//if (A == 0 || B == 0 || isNaN(A) || isNaN(B)) continue;

		out.labels.push(gene);
		//out.A.push(A);
		//out.B.push(B);
	}

	return out;
}

/**
 * Normalise array A (outer index x) to reference value valA in set A and reference
 * value valB in another set. ValA and valB must be corresponding variables.
 */
function normaliseSet(A, x, valA, valB) {
	let res = [];
	//const log = (argv.log2) ? Math.log2 : (argv.log10) ? Math.log10 : dummylog;

	for (var i=0; i<A.length; i++) {
		let n = (A[i][x] / valA * valB);
		res.push(n);
	}
	return res;
}

function logRatio(A,B) {
	var res = new Array(A.length);
	for (var i=0; i<A.length; i++) {
		res[i] = Math.log2(A[i] / B[i]);
	}
	return res;
}

function mad(A) {
	let m = median(A);
	let t = [];
	for (var i=0; i<A.length; i++) t.push(Math.abs(A[i] - m));
	return median(t);
}

function outlierFilter(A) {
	let t = A;
	let o;
	let s = 0;

	do {
		o = t;
		let m = median(o);
		let M = mad(o);
		let lim1 = m + 3 * 1.4826 * M;
		let lim2 = m - 3 * 1.4826 * M;
		t = [];
		for (var i=0; i<o.length; i++) {
			if (o[i] >= lim2 && o[i] <= lim1) t.push(o[i]);
			else s += Math.abs(o[i]);
		}
	} while (t.length < o.length);
	//console.log("Filtered",A.length - t.length, s / (A.length - t.length));
	return t;
}

function averageExp(A) {
	let a = 0.0;
	for (var i=0; i<A.length; i++) a+=A[i];
	return Math.pow(2,a / A.length);
	//return a/A.length;
}

function median(A) {
	let s = A.sort();
	let a = s[Math.floor(s.length/2)];
	return a;
	//return a/A.length;
}

function medianExp(A) {
	let s = A.sort();
	let a = s[Math.floor(s.length/2)];
	return Math.pow(2,a);
	//return a/A.length;
}

function average(A) {
	let a = 0.0;
	for (var i=0; i<A.length; i++) a+=A[i];
	return a / A.length;
	//return a/A.length;
}

function scaleSet(A,x) {
	return A.map((a) => x*a);
}

fs.readFile(infile, 'utf8', (err, data) => {
	if (err) throw err;

	let distDiff = [];
	let dSets = [];
	let input = parseCSV(data);

	let fExp = (argv.median) ? medianExp : averageExp;

	let r = [];
	for (var x=0; x<input.A.length; x++) {
	for (var y=0; y<input.B.length; y++) {
	r.push(input.A[x]);
	let ratio = fExp(outlierFilter(logRatio(input.A[x],input.B[y])));
	//console.log(ratio);
	r.push(scaleSet(input.B[y],ratio));
	}
	}

	if (argv.means) {
		for (var i=0; i<r.length; i+=2) {
			let avg = median(logRatio(r[i],r[i+1]));
			console.log(avg);
		}
		return;
	}

	//return;

	//let i = parseInt(argv.nrs);
	//let normB = normaliseBy(input, i);

	if (argv.group) {
		// Now group by protein and calculate fold change.
		let grouping = {};
		for (var j=0; j<input.labels.length; j++) {
			let l = input.labels[j];
			if (!grouping.hasOwnProperty(l)) {
				grouping[l] = {sums: [], name: l};
				for (var k=0; k<r.length; k++) grouping[l].sums.push(r[k][j]);
			} else {
				for (var k=0; k<r.length; k++) grouping[l].sums[k] += r[k][j];
			}
		}

		let grouparray = [];
		for (var x in grouping) {
			let g = grouping[x];

			// Now calculate fold changes
			g.folds = [];
			for (var ii=0; ii<g.sums.length; ii+=2) {
				var f = Math.log2(g.sums[ii] / g.sums[ii+1]);
				g.folds.push(f);
			}

			grouparray.push(g);
		}

		//console.log(grouping['slr1452']);

		process.stdout.write("Protein,");
		for (var k=0; k<input.A.length*input.B.length; k++) {
			process.stdout.write("Fold"+k+",");			
		}
		process.stdout.write("AvgFold\n");

		for (var j=0; j<grouparray.length; j++) {
			process.stdout.write(grouparray[j].name+",");
			let foldsum = 0;
			for (var k=0; k<grouparray[j].folds.length; k++) {
				foldsum += grouparray[j].folds[k];
				process.stdout.write(grouparray[j].folds[k]+",");
			}
			process.stdout.write((foldsum / grouparray[j].folds.length)+"\n");
		}

	} else {
		process.stdout.write("Protein,");
		for (var k=0; k<r.length; k+=2) {
			process.stdout.write("Fold"+(k/2)+",");			
		}
		process.stdout.write("AvgFold\n");

		for (var k=0; k<input.labels.length; k++) {
			process.stdout.write(input.labels[k]+",");
			let foldsum = 0;
			for (var j=0; j<r.length; j+=2) {
				let fold = Math.log2(r[j+1][k] / r[j][k]);
				foldsum += fold;
				process.stdout.write(fold+",");
			}
			process.stdout.write((foldsum / r.length)+"\n");
		}
	}
});

