/**
 * RuVector DNA Database Integration
 *
 * This script demonstrates how to use RuVector for DNA sequence
 * similarity search and genomic data management.
 */

const fs = require('fs');
const path = require('path');

// Import analysis results
const resultsPath = path.join(__dirname, '..', 'results', 'ruvector-data.json');

async function main() {
  console.log('\n' + '═'.repeat(80));
  console.log('              RUVECTOR DNA SEQUENCE DATABASE DEMO');
  console.log('═'.repeat(80));

  // Check if analysis results exist
  if (!fs.existsSync(resultsPath)) {
    console.log('\nError: Run dna-analyzer.js first to generate analysis data');
    console.log('Usage: node scripts/dna-analyzer.js');
    process.exit(1);
  }

  const data = JSON.parse(fs.readFileSync(resultsPath, 'utf-8'));
  console.log(`\nLoaded ${data.length} gene embeddings from analysis results`);

  // Simulate RuVector operations
  console.log('\n' + '-'.repeat(60));
  console.log('SIMULATING RUVECTOR DATABASE OPERATIONS');
  console.log('-'.repeat(60));

  // 1. Database Creation
  console.log('\n1. Creating RuVector Database');
  console.log(`   Dimensions: ${data[0].vector.length} (4-mer frequency vectors)`);
  console.log('   Index Type: HNSW (Hierarchical Navigable Small World)');
  console.log('   Distance Metric: Cosine Similarity');

  const dbConfig = {
    dimensions: data[0].vector.length,
    indexType: 'hnsw',
    distanceMetric: 'cosine',
    efConstruction: 200,
    M: 16
  };
  console.log(`   Config: ${JSON.stringify(dbConfig)}`);

  // 2. Insert Vectors
  console.log('\n2. Inserting Gene Vectors');
  for (const gene of data) {
    console.log(`   ✓ Inserted: ${gene.id}`);
    console.log(`     - Length: ${gene.metadata.length} bp`);
    console.log(`     - GC Content: ${gene.metadata.gcContent}%`);
    console.log(`     - ORFs: ${gene.metadata.orfCount}`);
  }

  // 3. Similarity Search Demo
  console.log('\n3. Similarity Search Results');
  console.log('   Query: Finding genes similar to TP53 variants');

  // Find TP53 variants
  const tp53Vectors = data.filter(d => d.id.includes('NM_000546') || d.id.includes('NM_001126112'));

  if (tp53Vectors.length >= 2) {
    const sim = cosineSimilarity(tp53Vectors[0].vector, tp53Vectors[1].vector);
    console.log(`\n   TP53 Variant Comparison:`);
    console.log(`   - ${tp53Vectors[0].id} vs ${tp53Vectors[1].id}`);
    console.log(`   - Similarity: ${(sim * 100).toFixed(2)}%`);
    console.log(`   - Interpretation: ${sim > 0.99 ? 'Nearly identical (expected for transcript variants)' :
                                       sim > 0.95 ? 'Highly similar' :
                                       sim > 0.8 ? 'Moderately similar' : 'Different sequences'}`);
  }

  // 4. Cross-gene comparison
  console.log('\n4. Cross-Gene Similarity Analysis');
  console.log('   Comparing different gene families:\n');

  const brca1 = data.find(d => d.id.includes('NM_007294'));
  const tp53 = data.find(d => d.id.includes('NM_000546'));
  const apoe = data.find(d => d.id.includes('NM_000041'));

  if (brca1 && tp53 && apoe) {
    const comparisons = [
      { name: 'BRCA1 vs TP53', sim: cosineSimilarity(brca1.vector, tp53.vector) },
      { name: 'BRCA1 vs APOE', sim: cosineSimilarity(brca1.vector, apoe.vector) },
      { name: 'TP53 vs APOE', sim: cosineSimilarity(tp53.vector, apoe.vector) }
    ];

    for (const comp of comparisons) {
      const bar = '█'.repeat(Math.round(comp.sim * 40));
      console.log(`   ${comp.name.padEnd(15)} ${(comp.sim * 100).toFixed(1).padStart(5)}% ${bar}`);
    }
  }

  // 5. Metadata Filtering
  console.log('\n5. Metadata-Based Filtering');
  console.log('   Query: Genes with GC content > 50% AND > 5 ORFs\n');

  const filtered = data.filter(d =>
    d.metadata.gcContent > 50 && d.metadata.orfCount > 5
  );

  if (filtered.length > 0) {
    for (const gene of filtered) {
      console.log(`   ✓ ${gene.id}`);
      console.log(`     GC: ${gene.metadata.gcContent}%, ORFs: ${gene.metadata.orfCount}`);
    }
  } else {
    console.log('   No genes match criteria');
  }

  // 6. Export for RuVector CLI
  console.log('\n6. RuVector CLI Commands');
  console.log('-'.repeat(60));

  const cliCommands = `
# Create the genomic database
ruvector create --dimensions ${data[0].vector.length} --path ./genomics.db

# Insert gene vectors
ruvector insert --db ./genomics.db --input results/ruvector-data.json --format json

# Search for similar sequences
ruvector search --db ./genomics.db --query "[${data[0].vector.slice(0, 5).join(', ')}...]" --top-k 3

# Database info
ruvector info --db ./genomics.db

# Export all data
ruvector export --db ./genomics.db --output backup.json --format json
`;
  console.log(cliCommands);

  // 7. Generate Summary Statistics
  console.log('7. Database Statistics');
  console.log('-'.repeat(60));

  const stats = {
    totalSequences: data.length,
    avgLength: Math.round(data.reduce((s, d) => s + d.metadata.length, 0) / data.length),
    avgGC: (data.reduce((s, d) => s + d.metadata.gcContent, 0) / data.length).toFixed(2),
    totalORFs: data.reduce((s, d) => s + d.metadata.orfCount, 0),
    totalMotifs: data.reduce((s, d) => s + d.metadata.motifCount, 0),
    vectorDimensions: data[0].vector.length
  };

  console.log(`   Total Sequences:    ${stats.totalSequences}`);
  console.log(`   Average Length:     ${stats.avgLength} bp`);
  console.log(`   Average GC Content: ${stats.avgGC}%`);
  console.log(`   Total ORFs:         ${stats.totalORFs}`);
  console.log(`   Total Motifs:       ${stats.totalMotifs}`);
  console.log(`   Vector Dimensions:  ${stats.vectorDimensions}`);

  console.log('\n' + '═'.repeat(80));
  console.log('                    RUVECTOR INTEGRATION COMPLETE');
  console.log('═'.repeat(80));
  console.log('\nData is ready for import into RuVector database.');
  console.log('Use the CLI commands above to create and populate the database.\n');
}

// Cosine similarity helper
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

main().catch(console.error);
