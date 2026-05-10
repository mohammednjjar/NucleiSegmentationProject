from pathlib import Path
import sys, tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from deprecated.misc.PageBreak import PageBreak
from deprecated.misc.StringFind import StringFind
from deprecated.misc.PlotCorr2 import PlotCorr2
from deprecated.probe_construction.WriteFasta import WriteFasta
from deprecated.probe_construction.OligoArrayCmd import OligoArrayCmd
from deprecated.probe_construction.BatchLaunchOligoArray import BatchLaunchOligoArray
from deprecated.probe_construction.CompileOligoArrayOutput import CompileOligoArrayOutput
from deprecated.probe_construction.AddCntrlSeqs import AddCntrlSeqs
from deprecated.probe_construction.GenProbe import GenProbe
from deprecated.reports.GenerateFPKMReport import GenerateFPKMReport
from deprecated.reports.GenerateHammingSphereReport import GenerateHammingSphereReport
from deprecated.reports.GenerateBitFlipReport import GenerateBitFlipReport
from deprecated.reports.GenerateMoleculeStatsReport import GenerateMoleculeStatsReport
from deprecated.reports.GenerateOnBitHistograms import GenerateOnBitHistograms


def test_misc_functions():
    assert PageBreak(display=False).startswith('---')
    idx, not_idx = StringFind(['abc','def','ab'], ['a','z'])
    assert idx == [[0, 2], []]
    assert not_idx == [1]
    corr = PlotCorr2([1,2,3], [1,2,4])
    assert corr['rho'] > 0.9


def test_fasta_and_commands(tmp_path):
    fasta = tmp_path / 'genes.fasta'
    WriteFasta(str(fasta), ['gene1','gene2'], ['ACGT'*400, 'TTTT'])
    assert fasta.read_text().count('>') == 2
    cmd, save = OligoArrayCmd(savePath=str(tmp_path) + '/', fastaName='abc.fasta')
    assert 'abc_oligos.txt' in cmd
    jobs = BatchLaunchOligoArray(cmd.replace('abc','genename'), str(fasta), savePath=str(tmp_path) + '/', execute=False, runExternal=False)
    assert len(jobs) == 2
    assert (tmp_path / 'gene1.fasta').exists()


def test_compile_and_controls(tmp_path):
    f = tmp_path / 'GeneA_oligos.txt'
    f.write_text('x\t1\t2\t3\t4\t5\t6\ttarget\tACGTAC\textra\n')
    data = CompileOligoArrayOutput(str(tmp_path) + '/', zeroOffTarget=False)
    assert data['GeneName'] == ['GeneA']
    picked = {'CommonName':['G1'], 'Sequence':[['AAAA','CCCC','GGGG','TTTT']], 'FivePrimeEnd':[[1,2,3,4]], 'IsoformName':['iso'], 'GeneName':['G1'], 'Nprobes':[4], 'FPKM':[1]}
    controls = [{'Header':'c1','Sequence':'ACAC'}, {'Header':'c2','Sequence':'TGTG'}, {'Header':'c3','Sequence':'CCCC'}, {'Header':'c4','Sequence':'GGGG'}]
    out = AddCntrlSeqs(picked, controls, 1, 1, 1, seed=1)
    assert len(out['CommonName']) == 3


def test_genprobe_and_reports(tmp_path):
    codebook = np.array([[1,1,0,0],[1,0,1,0]], dtype=np.uint8)
    genes = {'CommonName':['GeneA','blank001'], 'IsoformName':['iso',''], 'Sequence':[['ATGC','GGCC'], []], 'FivePrimeEnd':[[10,20], [0]]}
    secs = [{'Header':f'B{i+1}', 'Sequence':'ACGTACGT'} for i in range(4)]
    primers = [{'Header':f'P{i+1}', 'Sequence':'AAAACCCC'} for i in range(5)]
    seqs, names, pnum, params = GenProbe(codebook, genes, secs, primers, 2, 1, str(tmp_path) + '/', seed=2)
    assert seqs and names
    words = [{'geneName':'GeneA','isExactMatch':True,'isCorrectedMatch':False,'cellID':1,
              'intCodeword':12,'numHyb':4,'a':[1,2,3,4],'bg':[1,1,1,1],'h':[2,2,2,2],
              'numOnBits':2,'xc':[1,2,3,4],'yc':[2,3,4,5],'wordCentroidX':2.5,'wordCentroidY':3.5}]
    fpkm = [{'geneName':'GeneA','FPKM':5}]
    report, _ = GenerateFPKMReport(words, fpkm)
    assert report['countsWOUnknown'] == [1]
    hs, _ = GenerateHammingSphereReport(words, {'1100':'GeneA'}, reportsToGenerate=None)
    assert hs['hammingSphereCounts'][0,0] == 1
    bf, _ = GenerateBitFlipReport(words, {'1100':'GeneA'}, numHybs=4)
    assert 'scaledHybProb' in bf
    ms, _ = GenerateMoleculeStatsReport(words)
    assert ms['aN'][0] == 1
    figs, _ = GenerateOnBitHistograms(words, reportsToGenerate=[])
    assert figs == []
