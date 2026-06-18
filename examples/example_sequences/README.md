# Example Sequences

This directory contains four single-record FASTA files used as example protein sequences throughout the protpy documentation, tutorials, and test suite. All four sequences are coronavirus proteins sourced from [UniProtKB/Swiss-Prot](https://www.uniprot.org/).

---

## Files

### `test_fasta1.fasta`

**SARS-CoV-1 Spike glycoprotein (S)**

- **UniProt accession:** [P59594](https://www.uniprot.org/uniprot/P59594)
- **Organism:** *Severe acute respiratory syndrome coronavirus* (SARS-CoV-1, OX=694009)
- **Length:** 1,255 amino acids
- **Function:** Surface glycoprotein responsible for host cell attachment and membrane fusion. Binds the ACE2 receptor to mediate viral entry. The longest sequence in this set and a useful benchmark for descriptor computation time.

---

### `test_fasta2.fasta`

**SARS-CoV-1 Nucleoprotein (N)**

- **UniProt accession:** [P59595](https://www.uniprot.org/uniprot/P59595)
- **Organism:** *Severe acute respiratory syndrome coronavirus* (SARS-CoV-1, OX=694009)
- **Length:** 422 amino acids
- **Function:** RNA-binding protein that packages the viral genome into a ribonucleoprotein complex. Highly expressed during infection and widely used as a diagnostic antigen.

---

### `test_fasta3.fasta`

**SARS-CoV-1 Envelope small membrane protein (E)**

- **UniProt accession:** [P59637](https://www.uniprot.org/uniprot/P59637)
- **Organism:** *Severe acute respiratory syndrome coronavirus* (SARS-CoV-1, OX=694009)
- **Length:** 76 amino acids
- **Function:** Small integral membrane protein involved in viral assembly, budding, and release. The shortest sequence in this set; useful for fast tests and progress-bar demonstrations.

---

### `test_fasta4.fasta`

**SARS-CoV-2 Nucleoprotein (N)**

- **UniProt accession:** [P0DTC9](https://www.uniprot.org/uniprot/P0DTC9)
- **Organism:** *Severe acute respiratory syndrome coronavirus 2* (SARS-CoV-2, OX=2697049)
- **Length:** 419 amino acids
- **Function:** RNA-binding nucleoprotein of SARS-CoV-2. Closely related to the SARS-CoV-1 nucleoprotein (`test_fasta2.fasta`); the pair is useful for comparing descriptors between orthologous proteins.

---

## Usage

```python
import protpy

# Load a single sequence from a FASTA file
seqs = protpy.read_sequences_from_file("examples/example_sequences/test_fasta1.fasta")
protein_seq = seqs[0]

# Compute all descriptors for that sequence
df = protpy.all_descriptors(protein_seq)
# Shape: 1 × 10,757

# Compute descriptors for all four sequences at once
import os
fasta_files = [
    f"examples/example_sequences/test_fasta{i}.fasta" for i in range(1, 5)
]
sequences = [protpy.read_sequences_from_file(f)[0] for f in fasta_files]
result = protpy.encode_dataset(sequences, output_file=None)
# Shape: 4 × 10,758  (first column: sequence_id = seq1…seq4)
```

---

## Summary

| File | Protein | Organism | UniProt | Length |
|------|---------|----------|---------|--------|
| `test_fasta1.fasta` | Spike glycoprotein (S) | SARS-CoV-1 | [P59594](https://www.uniprot.org/uniprot/P59594) | 1,255 aa |
| `test_fasta2.fasta` | Nucleoprotein (N) | SARS-CoV-1 | [P59595](https://www.uniprot.org/uniprot/P59595) | 422 aa |
| `test_fasta3.fasta` | Envelope protein (E) | SARS-CoV-1 | [P59637](https://www.uniprot.org/uniprot/P59637) | 76 aa |
| `test_fasta4.fasta` | Nucleoprotein (N) | SARS-CoV-2 | [P0DTC9](https://www.uniprot.org/uniprot/P0DTC9) | 419 aa |
