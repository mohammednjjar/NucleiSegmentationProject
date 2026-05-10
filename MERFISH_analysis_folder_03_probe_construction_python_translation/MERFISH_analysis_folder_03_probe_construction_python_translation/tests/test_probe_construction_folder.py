from pathlib import Path
import tempfile
import numpy as np

from probe_construction import OTMap, OTMap2, OTTable, Transcriptome, TRDesigner, TargetRegions, PrimerDesigner


def test_otmap_accumulates_and_reads_values():
    m = OTMap([[1, 2, 1], [10, 20, 5]])
    assert m.length() == 2
    assert np.allclose(m.GetValues([1, 2, 3]), [15, 20, 0])
    m2 = OTMap2([[1, 2, 1], [10, 20, 5]])
    assert np.allclose(m2.GetValues([1, 2, 3]), [15, 20, 0])


def test_transcriptome_lookup_and_slice():
    t = Transcriptome([['id1', 'id2'], ['geneA', 'geneB'], ['ACGTACGT', 'TTTTCCCC'], [1.0, 2.0]])
    assert t.numTranscripts == 2
    assert t.GetSequenceByID('id1') == 'ACGTACGT'
    assert t.GetIDsByName('geneA') == ['id1']
    s = t.Slice(geneName=['geneB'])
    assert s.numTranscripts == 1
    assert s.GetSequenceByID('id2') == 'TTTTCCCC'


def test_ottable_penalty_save_load():
    t = Transcriptome([['id1'], ['geneA'], ['ACGTACGT'], [1.0]])
    ot = OTTable(t, 3)
    penalty, hashes = ot.CalculatePenalty('ACGT')
    assert penalty.shape[0] == 2
    assert np.all(penalty > 0)
    with tempfile.TemporaryDirectory() as tmp:
        ot.Save(tmp)
        loaded = OTTable.Load(tmp)
        p2, _ = loaded.CalculatePenalty('ACGT')
        assert np.allclose(penalty, p2)


def test_trdesigner_region_properties_and_design():
    t = Transcriptome([['id1'], ['geneA'], ['ACGTACGTACGT'], [1.0]])
    ot = OTTable(t, 3)
    trd = TRDesigner(transcriptome=t, OTTables=[ot], OTTableNames=['self'], specificityTable=ot, verbose=False)
    gc, ids, names = trd.GetRegionGC(4)
    assert len(gc[0]) == 9
    regs = trd.DesignTargetRegions(regionLength=4, GC=[0.0, 1.0], Tm=[-100, 100], specificity=[0, 10])
    assert len(regs) == 1
    assert isinstance(regs[0], TargetRegions)
    assert regs[0].numRegions > 0


def test_targetregions_fasta_and_save_load():
    tr = TargetRegions(geneName='geneA', id='id1', geneSequence='ACGTACGT', startPos=[1, 5], regionLength=[4, 4], GC=[0.5, 0.5], Tm=[10, 11], specificity=[1, 1], isoSpecificity=[1, 1])
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / 'regions.fasta'
        tr.fastawrite(fasta, overwrite=True)
        assert fasta.read_text().count('>') == 2
        tr.Save(Path(tmp) / 'tr')
        loaded = TargetRegions.Load(Path(tmp) / 'tr')
        assert loaded.numRegions == 2


def test_primerdesigner_filters_and_fasta():
    pd = PrimerDesigner(seqs=['ACGTACGT', 'AAAACCCC', 'TGCATGCA'], primerLength=8, numPrimersToGenerate=0, verbose=False)
    assert pd.numPrimers == 3
    keep = pd.RemoveForbiddenSeqs(['AAAA'])
    assert keep.tolist() == [True, False, True]
    pd.RemoveSelfCompPrimers(homologyMax=4)
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / 'primers.fasta'
        pd.WriteFasta(fasta, namePrefix='p')
        assert fasta.exists()
