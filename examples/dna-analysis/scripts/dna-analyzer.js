/**
 * DNA Deep Analysis using RuVector
 *
 * This script performs comprehensive DNA sequence analysis including:
 * - K-mer frequency analysis
 * - GC content calculation
 * - Sequence similarity using vector embeddings
 * - Codon usage analysis
 * - Motif detection
 * - Structural feature prediction
 */

const fs = require('fs');
const path = require('path');

// DNA Analysis Constants
const NUCLEOTIDES = ['A', 'T', 'G', 'C'];
const CODON_TABLE = {
  'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
  'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
  'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
  'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
  'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
  'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
  'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
  'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
  'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
  'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
};

// Known regulatory motifs
const REGULATORY_MOTIFS = {
  'TATAAA': 'TATA Box (Promoter)',
  'CAAT': 'CAAT Box (Promoter)',
  'GCCGCC': 'GC Box (Promoter)',
  'AATAAA': 'Polyadenylation Signal',
  'CTCF': 'CTCF Binding Site',
  'CACGTG': 'E-box (Enhancer)',
  'GATA': 'GATA Binding Site',
  'CCAAT': 'CCAAT Box',
  'GCGCGC': 'CpG Island Marker'
};

class DNAAnalyzer {
  constructor(sequence, name = 'Unknown') {
    this.sequence = sequence.toUpperCase().replace(/[^ATGC]/g, '');
    this.name = name;
    this.length = this.sequence.length;
  }

  // Basic Statistics
  getBasicStats() {
    const counts = { A: 0, T: 0, G: 0, C: 0 };
    for (const base of this.sequence) {
      if (counts[base] !== undefined) counts[base]++;
    }

    const gcContent = ((counts.G + counts.C) / this.length) * 100;
    const atContent = ((counts.A + counts.T) / this.length) * 100;

    return {
      length: this.length,
      baseCounts: counts,
      gcContent: gcContent.toFixed(2),
      atContent: atContent.toFixed(2),
      gcSkew: ((counts.G - counts.C) / (counts.G + counts.C)).toFixed(4),
      atSkew: ((counts.A - counts.T) / (counts.A + counts.T)).toFixed(4)
    };
  }

  // K-mer Frequency Analysis
  getKmerFrequencies(k = 3) {
    const kmers = {};
    for (let i = 0; i <= this.sequence.length - k; i++) {
      const kmer = this.sequence.substring(i, i + k);
      kmers[kmer] = (kmers[kmer] || 0) + 1;
    }

    // Sort by frequency
    const sorted = Object.entries(kmers)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 20);

    return {
      k,
      totalKmers: Object.keys(kmers).length,
      topKmers: sorted,
      uniqueRatio: (Object.keys(kmers).length / (this.length - k + 1)).toFixed(4)
    };
  }

  // Generate k-mer embedding vector for similarity search
  generateKmerEmbedding(k = 4) {
    const possibleKmers = [];
    const generateKmers = (prefix, len) => {
      if (len === 0) {
        possibleKmers.push(prefix);
        return;
      }
      for (const n of NUCLEOTIDES) {
        generateKmers(prefix + n, len - 1);
      }
    };
    generateKmers('', k);

    // Count k-mers in sequence
    const kmerCounts = {};
    for (let i = 0; i <= this.sequence.length - k; i++) {
      const kmer = this.sequence.substring(i, i + k);
      kmerCounts[kmer] = (kmerCounts[kmer] || 0) + 1;
    }

    // Create normalized frequency vector
    const total = this.length - k + 1;
    const vector = possibleKmers.map(kmer => (kmerCounts[kmer] || 0) / total);

    return {
      dimensions: possibleKmers.length,
      vector,
      kmerMap: possibleKmers
    };
  }

  // Codon Usage Analysis
  getCodonUsage() {
    const codonCounts = {};
    const aminoAcidCounts = {};

    for (let i = 0; i <= this.sequence.length - 3; i += 3) {
      const codon = this.sequence.substring(i, i + 3);
      if (CODON_TABLE[codon]) {
        codonCounts[codon] = (codonCounts[codon] || 0) + 1;
        const aa = CODON_TABLE[codon];
        aminoAcidCounts[aa] = (aminoAcidCounts[aa] || 0) + 1;
      }
    }

    // Calculate codon bias (RSCU - Relative Synonymous Codon Usage)
    const rscu = {};
    const synonymousCodons = {};

    for (const [codon, aa] of Object.entries(CODON_TABLE)) {
      if (!synonymousCodons[aa]) synonymousCodons[aa] = [];
      synonymousCodons[aa].push(codon);
    }

    for (const [aa, codons] of Object.entries(synonymousCodons)) {
      const total = codons.reduce((sum, c) => sum + (codonCounts[c] || 0), 0);
      const expected = total / codons.length;
      for (const codon of codons) {
        if (total > 0) {
          rscu[codon] = ((codonCounts[codon] || 0) / expected).toFixed(3);
        }
      }
    }

    return {
      codonCounts,
      aminoAcidCounts,
      rscu,
      totalCodons: Math.floor(this.length / 3)
    };
  }

  // Find Regulatory Motifs
  findMotifs() {
    const found = [];

    for (const [motif, description] of Object.entries(REGULATORY_MOTIFS)) {
      let index = this.sequence.indexOf(motif);
      while (index !== -1) {
        found.push({
          motif,
          description,
          position: index,
          context: this.sequence.substring(Math.max(0, index - 10), Math.min(this.length, index + motif.length + 10))
        });
        index = this.sequence.indexOf(motif, index + 1);
      }
    }

    return found.sort((a, b) => a.position - b.position);
  }

  // Open Reading Frame Detection
  findORFs(minLength = 100) {
    const startCodon = 'ATG';
    const stopCodons = ['TAA', 'TAG', 'TGA'];
    const orfs = [];

    for (let frame = 0; frame < 3; frame++) {
      let i = frame;
      let orfStart = -1;

      while (i <= this.sequence.length - 3) {
        const codon = this.sequence.substring(i, i + 3);

        if (codon === startCodon && orfStart === -1) {
          orfStart = i;
        } else if (stopCodons.includes(codon) && orfStart !== -1) {
          const orfLength = i - orfStart + 3;
          if (orfLength >= minLength) {
            orfs.push({
              frame: frame + 1,
              start: orfStart,
              end: i + 3,
              length: orfLength,
              proteinLength: Math.floor(orfLength / 3),
              sequence: this.sequence.substring(orfStart, i + 3)
            });
          }
          orfStart = -1;
        }
        i += 3;
      }
    }

    return orfs.sort((a, b) => b.length - a.length);
  }

  // CpG Island Detection
  findCpGIslands(windowSize = 200, minGC = 50, minObsExp = 0.6) {
    const islands = [];

    for (let i = 0; i <= this.sequence.length - windowSize; i += 50) {
      const window = this.sequence.substring(i, i + windowSize);

      let gc = 0, cpg = 0, c = 0, g = 0;
      for (let j = 0; j < window.length; j++) {
        if (window[j] === 'G') { g++; gc++; }
        if (window[j] === 'C') { c++; gc++; }
        if (j < window.length - 1 && window[j] === 'C' && window[j + 1] === 'G') {
          cpg++;
        }
      }

      const gcPercent = (gc / windowSize) * 100;
      const expectedCpG = (c * g) / windowSize;
      const obsExp = expectedCpG > 0 ? cpg / expectedCpG : 0;

      if (gcPercent >= minGC && obsExp >= minObsExp) {
        islands.push({
          start: i,
          end: i + windowSize,
          gcContent: gcPercent.toFixed(2),
          obsExpRatio: obsExp.toFixed(3),
          cpgCount: cpg
        });
      }
    }

    // Merge overlapping islands
    const merged = [];
    for (const island of islands) {
      if (merged.length === 0 || island.start > merged[merged.length - 1].end) {
        merged.push(island);
      } else {
        merged[merged.length - 1].end = Math.max(merged[merged.length - 1].end, island.end);
      }
    }

    return merged;
  }

  // Repeat Sequence Detection
  findRepeats(minLength = 6, minOccurrences = 3) {
    const repeats = {};

    for (let len = minLength; len <= 20; len++) {
      for (let i = 0; i <= this.sequence.length - len; i++) {
        const pattern = this.sequence.substring(i, i + len);
        if (!repeats[pattern]) {
          repeats[pattern] = [];
        }
        repeats[pattern].push(i);
      }
    }

    // Filter by minimum occurrences
    const significantRepeats = Object.entries(repeats)
      .filter(([_, positions]) => positions.length >= minOccurrences)
      .map(([pattern, positions]) => ({
        pattern,
        length: pattern.length,
        occurrences: positions.length,
        positions: positions.slice(0, 10)
      }))
      .sort((a, b) => (b.length * b.occurrences) - (a.length * a.occurrences))
      .slice(0, 20);

    return significantRepeats;
  }

  // Nucleotide Complexity (Shannon Entropy)
  calculateEntropy(windowSize = 100) {
    const entropies = [];

    for (let i = 0; i <= this.sequence.length - windowSize; i += windowSize / 2) {
      const window = this.sequence.substring(i, i + windowSize);
      const counts = { A: 0, T: 0, G: 0, C: 0 };

      for (const base of window) {
        if (counts[base] !== undefined) counts[base]++;
      }

      let entropy = 0;
      for (const count of Object.values(counts)) {
        if (count > 0) {
          const p = count / windowSize;
          entropy -= p * Math.log2(p);
        }
      }

      entropies.push({
        position: i,
        entropy: entropy.toFixed(4),
        complexity: entropy > 1.8 ? 'high' : entropy > 1.5 ? 'medium' : 'low'
      });
    }

    const avgEntropy = entropies.reduce((sum, e) => sum + parseFloat(e.entropy), 0) / entropies.length;

    return {
      windowSize,
      averageEntropy: avgEntropy.toFixed(4),
      maxEntropy: 2.0,
      windows: entropies
    };
  }

  // Translate to Protein
  translate(frame = 0) {
    let protein = '';
    for (let i = frame; i <= this.sequence.length - 3; i += 3) {
      const codon = this.sequence.substring(i, i + 3);
      const aa = CODON_TABLE[codon] || 'X';
      protein += aa;
    }
    return protein;
  }

  // Run Complete Analysis
  runCompleteAnalysis() {
    console.log(`\n${'='.repeat(80)}`);
    console.log(`DEEP DNA ANALYSIS: ${this.name}`);
    console.log(`${'='.repeat(80)}\n`);

    const results = {
      geneName: this.name,
      timestamp: new Date().toISOString(),
      basicStats: this.getBasicStats(),
      kmerAnalysis: {
        k3: this.getKmerFrequencies(3),
        k4: this.getKmerFrequencies(4),
        k6: this.getKmerFrequencies(6)
      },
      codonUsage: this.getCodonUsage(),
      regulatoryMotifs: this.findMotifs(),
      orfs: this.findORFs(100),
      cpgIslands: this.findCpGIslands(),
      repeats: this.findRepeats(),
      entropy: this.calculateEntropy(),
      embedding: this.generateKmerEmbedding(4)
    };

    // Print Summary
    console.log('1. BASIC STATISTICS');
    console.log('-'.repeat(40));
    console.log(`   Sequence Length: ${results.basicStats.length} bp`);
    console.log(`   GC Content: ${results.basicStats.gcContent}%`);
    console.log(`   AT Content: ${results.basicStats.atContent}%`);
    console.log(`   GC Skew: ${results.basicStats.gcSkew}`);
    console.log(`   AT Skew: ${results.basicStats.atSkew}`);
    console.log(`   Base Composition: A=${results.basicStats.baseCounts.A}, T=${results.basicStats.baseCounts.T}, G=${results.basicStats.baseCounts.G}, C=${results.basicStats.baseCounts.C}`);

    console.log('\n2. K-MER ANALYSIS');
    console.log('-'.repeat(40));
    console.log(`   3-mers: ${results.kmerAnalysis.k3.totalKmers} unique (ratio: ${results.kmerAnalysis.k3.uniqueRatio})`);
    console.log(`   4-mers: ${results.kmerAnalysis.k4.totalKmers} unique (ratio: ${results.kmerAnalysis.k4.uniqueRatio})`);
    console.log(`   Top 3-mers: ${results.kmerAnalysis.k3.topKmers.slice(0, 5).map(k => `${k[0]}(${k[1]})`).join(', ')}`);

    console.log('\n3. CODON USAGE');
    console.log('-'.repeat(40));
    console.log(`   Total Codons: ${results.codonUsage.totalCodons}`);
    const topAA = Object.entries(results.codonUsage.aminoAcidCounts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 5);
    console.log(`   Most Common Amino Acids: ${topAA.map(([aa, count]) => `${aa}(${count})`).join(', ')}`);

    console.log('\n4. REGULATORY MOTIFS');
    console.log('-'.repeat(40));
    if (results.regulatoryMotifs.length > 0) {
      for (const motif of results.regulatoryMotifs.slice(0, 10)) {
        console.log(`   ${motif.motif} at pos ${motif.position}: ${motif.description}`);
      }
    } else {
      console.log('   No known regulatory motifs found');
    }

    console.log('\n5. OPEN READING FRAMES');
    console.log('-'.repeat(40));
    console.log(`   Total ORFs (>100bp): ${results.orfs.length}`);
    for (const orf of results.orfs.slice(0, 5)) {
      console.log(`   Frame ${orf.frame}: ${orf.start}-${orf.end} (${orf.length}bp, ${orf.proteinLength}aa)`);
    }

    console.log('\n6. CpG ISLANDS');
    console.log('-'.repeat(40));
    console.log(`   Total CpG Islands: ${results.cpgIslands.length}`);
    for (const island of results.cpgIslands.slice(0, 5)) {
      console.log(`   ${island.start}-${island.end}: GC=${island.gcContent}%, Obs/Exp=${island.obsExpRatio}`);
    }

    console.log('\n7. REPEAT SEQUENCES');
    console.log('-'.repeat(40));
    console.log(`   Significant Repeats Found: ${results.repeats.length}`);
    for (const repeat of results.repeats.slice(0, 5)) {
      console.log(`   "${repeat.pattern}" (${repeat.length}bp): ${repeat.occurrences} occurrences`);
    }

    console.log('\n8. SEQUENCE COMPLEXITY (Shannon Entropy)');
    console.log('-'.repeat(40));
    console.log(`   Average Entropy: ${results.entropy.averageEntropy} / ${results.entropy.maxEntropy}`);
    const complexityDist = results.entropy.windows.reduce((acc, w) => {
      acc[w.complexity] = (acc[w.complexity] || 0) + 1;
      return acc;
    }, {});
    console.log(`   Complexity Distribution: High=${complexityDist.high || 0}, Medium=${complexityDist.medium || 0}, Low=${complexityDist.low || 0}`);

    console.log('\n9. VECTOR EMBEDDING');
    console.log('-'.repeat(40));
    console.log(`   Embedding Dimensions: ${results.embedding.dimensions}`);
    console.log(`   Non-zero Features: ${results.embedding.vector.filter(v => v > 0).length}`);

    return results;
  }
}

// Parse FASTA file
function parseFasta(content) {
  const sequences = [];
  const lines = content.split('\n');
  let currentName = '';
  let currentSeq = '';

  for (const line of lines) {
    if (line.startsWith('>')) {
      if (currentName && currentSeq) {
        sequences.push({ name: currentName, sequence: currentSeq });
      }
      currentName = line.substring(1).trim();
      currentSeq = '';
    } else {
      currentSeq += line.trim();
    }
  }

  if (currentName && currentSeq) {
    sequences.push({ name: currentName, sequence: currentSeq });
  }

  return sequences;
}

// Calculate sequence similarity using cosine similarity
function cosineSimilarity(vec1, vec2) {
  let dotProduct = 0;
  let norm1 = 0;
  let norm2 = 0;

  for (let i = 0; i < vec1.length; i++) {
    dotProduct += vec1[i] * vec2[i];
    norm1 += vec1[i] * vec1[i];
    norm2 += vec2[i] * vec2[i];
  }

  return dotProduct / (Math.sqrt(norm1) * Math.sqrt(norm2));
}

// Main Analysis Function
async function main() {
  console.log('\n' + '═'.repeat(80));
  console.log('                    RUVECTOR DNA DEEP ANALYSIS SYSTEM');
  console.log('                    Public Gene Sequence Analysis');
  console.log('═'.repeat(80));

  // Load FASTA file
  const fastaPath = path.join(__dirname, '..', 'data', 'genes.fasta');
  const fastaContent = fs.readFileSync(fastaPath, 'utf-8');
  const sequences = parseFasta(fastaContent);

  console.log(`\nLoaded ${sequences.length} sequences from FASTA file`);
  console.log('Genes: ' + sequences.map(s => s.name.split(' ')[1] || s.name.split(' ')[0]).join(', '));

  const allResults = [];
  const embeddings = [];

  // Analyze each sequence
  for (const seq of sequences) {
    const analyzer = new DNAAnalyzer(seq.sequence, seq.name);
    const results = analyzer.runCompleteAnalysis();
    allResults.push(results);
    embeddings.push({
      name: results.geneName,
      vector: results.embedding.vector
    });
  }

  // Similarity Analysis
  console.log('\n' + '═'.repeat(80));
  console.log('                    SEQUENCE SIMILARITY MATRIX');
  console.log('═'.repeat(80));
  console.log('\nCosine Similarity based on 4-mer frequency vectors:\n');

  const geneNames = embeddings.map(e => {
    const parts = e.name.split(' ');
    return parts[1] ? parts[1].replace('(', '').replace(')', '').substring(0, 8) : parts[0].substring(0, 8);
  });

  // Print header
  console.log('          ' + geneNames.map(n => n.padStart(10)).join(''));

  for (let i = 0; i < embeddings.length; i++) {
    let row = geneNames[i].padEnd(10);
    for (let j = 0; j < embeddings.length; j++) {
      const sim = cosineSimilarity(embeddings[i].vector, embeddings[j].vector);
      row += (sim * 100).toFixed(1).padStart(10) + '%';
    }
    console.log(row);
  }

  // Comparative Analysis
  console.log('\n' + '═'.repeat(80));
  console.log('                    COMPARATIVE GENE ANALYSIS');
  console.log('═'.repeat(80));

  console.log('\n1. GC CONTENT COMPARISON');
  console.log('-'.repeat(60));
  for (const result of allResults) {
    const shortName = result.geneName.split(' ')[1] || result.geneName;
    const gcBar = '█'.repeat(Math.round(parseFloat(result.basicStats.gcContent) / 2));
    console.log(`   ${shortName.padEnd(15)} ${result.basicStats.gcContent}% ${gcBar}`);
  }

  console.log('\n2. SEQUENCE LENGTH COMPARISON');
  console.log('-'.repeat(60));
  const maxLen = Math.max(...allResults.map(r => r.basicStats.length));
  for (const result of allResults) {
    const shortName = result.geneName.split(' ')[1] || result.geneName;
    const lenBar = '█'.repeat(Math.round((result.basicStats.length / maxLen) * 40));
    console.log(`   ${shortName.padEnd(15)} ${result.basicStats.length.toString().padStart(6)} bp ${lenBar}`);
  }

  console.log('\n3. ORF DENSITY');
  console.log('-'.repeat(60));
  for (const result of allResults) {
    const shortName = result.geneName.split(' ')[1] || result.geneName;
    const orfDensity = (result.orfs.length / (result.basicStats.length / 1000)).toFixed(2);
    console.log(`   ${shortName.padEnd(15)} ${result.orfs.length} ORFs (${orfDensity} per kb)`);
  }

  console.log('\n4. REGULATORY MOTIF DENSITY');
  console.log('-'.repeat(60));
  for (const result of allResults) {
    const shortName = result.geneName.split(' ')[1] || result.geneName;
    const motifDensity = (result.regulatoryMotifs.length / (result.basicStats.length / 1000)).toFixed(2);
    console.log(`   ${shortName.padEnd(15)} ${result.regulatoryMotifs.length} motifs (${motifDensity} per kb)`);
  }

  console.log('\n5. SEQUENCE COMPLEXITY');
  console.log('-'.repeat(60));
  for (const result of allResults) {
    const shortName = result.geneName.split(' ')[1] || result.geneName;
    const entropy = parseFloat(result.entropy.averageEntropy);
    const complexityBar = '█'.repeat(Math.round(entropy * 20));
    console.log(`   ${shortName.padEnd(15)} Entropy: ${result.entropy.averageEntropy} ${complexityBar}`);
  }

  // Save results for ruvector
  const ruvectorData = allResults.map(result => ({
    id: result.geneName.split(' ')[0],
    vector: result.embedding.vector,
    metadata: {
      geneName: result.geneName,
      length: result.basicStats.length,
      gcContent: parseFloat(result.basicStats.gcContent),
      orfCount: result.orfs.length,
      cpgIslands: result.cpgIslands.length,
      entropy: parseFloat(result.entropy.averageEntropy),
      motifCount: result.regulatoryMotifs.length,
      repeatCount: result.repeats.length
    }
  }));

  const outputPath = path.join(__dirname, '..', 'results', 'ruvector-data.json');
  fs.writeFileSync(outputPath, JSON.stringify(ruvectorData, null, 2));
  console.log(`\nRuVector-compatible data saved to: ${outputPath}`);

  // Save full analysis report
  const reportPath = path.join(__dirname, '..', 'results', 'analysis-report.json');
  fs.writeFileSync(reportPath, JSON.stringify(allResults, null, 2));
  console.log(`Full analysis report saved to: ${reportPath}`);

  console.log('\n' + '═'.repeat(80));
  console.log('                         ANALYSIS COMPLETE');
  console.log('═'.repeat(80));
  console.log('\nKey Findings:');

  // Find most similar pair
  let maxSim = 0;
  let simPair = ['', ''];
  for (let i = 0; i < embeddings.length; i++) {
    for (let j = i + 1; j < embeddings.length; j++) {
      const sim = cosineSimilarity(embeddings[i].vector, embeddings[j].vector);
      if (sim > maxSim) {
        maxSim = sim;
        simPair = [geneNames[i], geneNames[j]];
      }
    }
  }
  console.log(`  • Most similar sequences: ${simPair[0]} and ${simPair[1]} (${(maxSim * 100).toFixed(1)}% similarity)`);

  // Find gene with highest GC content
  const highestGC = allResults.reduce((max, r) =>
    parseFloat(r.basicStats.gcContent) > parseFloat(max.basicStats.gcContent) ? r : max
  );
  const gcName = highestGC.geneName.split(' ')[1] || highestGC.geneName;
  console.log(`  • Highest GC content: ${gcName} (${highestGC.basicStats.gcContent}%)`);

  // Find gene with most ORFs
  const mostORFs = allResults.reduce((max, r) => r.orfs.length > max.orfs.length ? r : max);
  const orfName = mostORFs.geneName.split(' ')[1] || mostORFs.geneName;
  console.log(`  • Most ORFs: ${orfName} (${mostORFs.orfs.length} ORFs)`);

  // Clinical relevance
  console.log('\nClinical Relevance:');
  console.log('  • BRCA1: DNA repair gene; mutations linked to breast/ovarian cancer');
  console.log('  • TP53: Tumor suppressor; most frequently mutated gene in human cancers');
  console.log('  • APOE: Lipid metabolism; ε4 variant is major risk factor for Alzheimer\'s');

  return { allResults, ruvectorData, embeddings };
}

// Run analysis
main().catch(console.error);
