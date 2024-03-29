package org.exbio.tfprio.configs.modules;

import org.exbio.pipejar.configs.ConfigModule;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.ExternalConfig;
import org.exbio.pipejar.configs.ConfigTypes.InputTypes.InternalConfig;

import java.io.File;

public class TEPIC extends ConfigModule {
    public final ExternalConfig<Boolean> randomize = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<File> referenceGenome = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> PWMs = new ExternalConfig<>(File.class);
    public final ExternalConfig<File> bedChromatinSignal = new ExternalConfig<>(File.class);
    public final ExternalConfig<Integer> columnBedfile = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Integer> windowSize = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<File> onlyDNasePeaks = new ExternalConfig<>(File.class);
    public final ExternalConfig<Boolean> exponentialDecay = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Boolean> doNotNormalizePeakLength = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Boolean> originalDecay = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<File> psemsLengthFile = new ExternalConfig<>(File.class);
    public final ExternalConfig<Boolean> entireGeneBody = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Boolean> doZip = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<File> twoBitFile = new ExternalConfig<>(File.class);
    public final ExternalConfig<Double> pValue = new ExternalConfig<>(Double.class);
    public final ExternalConfig<Integer> maxMinutesPerChromosome = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Boolean> chromosomePrefix = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Boolean> transcriptBased = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<File> loopListFile = new ExternalConfig<>(File.class);
    public final ExternalConfig<Integer> loopWindows = new ExternalConfig<>(Integer.class);
    public final ExternalConfig<Boolean> onlyPeakFeatures = new ExternalConfig<>(Boolean.class);
    public final ExternalConfig<Double> affinityCutoff = new ExternalConfig<>(Double.class);
    public final InternalConfig<String> sequenceFileName = new InternalConfig<>("trap_sequences.tsv");
    public final InternalConfig<String> geneRegionFileName = new InternalConfig<>("regions_to_target_genes.tsv");
}
