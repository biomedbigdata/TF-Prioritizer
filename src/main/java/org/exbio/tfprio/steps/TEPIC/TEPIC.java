package org.exbio.tfprio.steps.TEPIC;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.OptionalConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.UsageConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.*;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class TEPIC extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final Map<String, Map<String, Collection<InputFile>>> chipSeqFiles = new HashMap<>();
    private final RequiredConfig<String> sequenceFileName = new RequiredConfig<>(configs.tepic.sequenceFileName);

    private final RequiredConfig<File> f_pwms = new RequiredConfig<>(configs.tepic.PWMs);
    private final OptionalConfig<File> bedChromatinSignal =
            new OptionalConfig<>(configs.tepic.bedChromatinSignal, false);
    private final OptionalConfig<Integer> columnBedfile = new OptionalConfig<>(configs.tepic.columnBedfile, false);
    private final OptionalConfig<File> geneAnnotationFile =
            new OptionalConfig<>(configs.inputConfigs.geneAnnotationFile, false);
    private final OptionalConfig<Integer> windowSize = new OptionalConfig<>(configs.tepic.windowSize, false);
    private final OptionalConfig<File> onlyDNasePeaks = new OptionalConfig<>(configs.tepic.onlyDNasePeaks, false);
    private final OptionalConfig<Boolean> exponentialDecay =
            new OptionalConfig<>(configs.tepic.exponentialDecay, false);
    private final OptionalConfig<Boolean> doNotNormalizePeakLength =
            new OptionalConfig<>(configs.tepic.doNotNormalizePeakLength, false);
    private final OptionalConfig<Boolean> originalDecay = new OptionalConfig<>(configs.tepic.originalDecay, false);
    private final OptionalConfig<File> psemsLengthFile = new OptionalConfig<>(configs.tepic.psemsLengthFile, false);
    private final OptionalConfig<Boolean> entireGeneBody = new OptionalConfig<>(configs.tepic.entireGeneBody, false);
    private final OptionalConfig<Boolean> doZip = new OptionalConfig<>(configs.tepic.doZip, false);
    private final OptionalConfig<File> twoBitFile = new OptionalConfig<>(configs.tepic.twoBitFile, false);
    private final OptionalConfig<Double> pValue = new OptionalConfig<>(configs.tepic.pValue, false);
    private final OptionalConfig<Integer> maxMinutesPerChromosome =
            new OptionalConfig<>(configs.tepic.maxMinutesPerChromosome, false);
    private final OptionalConfig<Boolean> chromosomePrefix =
            new OptionalConfig<>(configs.tepic.chromosomePrefix, false);
    private final OptionalConfig<Boolean> transcriptBased = new OptionalConfig<>(configs.tepic.transcriptBased, false);
    private final OptionalConfig<File> loopListFile = new OptionalConfig<>(configs.tepic.loopListFile, false);
    private final OptionalConfig<Integer> loopWindows = new OptionalConfig<>(configs.tepic.loopWindows, false);
    private final OptionalConfig<Boolean> onlyPeakFeatures =
            new OptionalConfig<>(configs.tepic.onlyPeakFeatures, false);
    private final InputFile referenceGenome;

    private final File executable;

    public TEPIC(Configs configs, Map<String, Map<String, Collection<OutputFile>>> latestChipSeq,
                 OutputFile referenceGenome) {
        super(configs, false,
                latestChipSeq.values().stream().flatMap(x -> x.values().stream()).flatMap(Collection::stream).collect(
                        Collectors.toSet()), referenceGenome);

        InputFile tepicDirectory = new InputFile(inputDirectory, "TEPIC");
        this.referenceGenome = addInput(referenceGenome);

        String path = "/org/exbio/tfprio/steps/TEPIC/Code";
        URL resource = getClass().getResource(path);
        if (resource == null || !new File(resource.getFile()).exists()) {
            // Usually entered when running via Jar
            logger.info("Jar mode");
            try {
                copyResources(path, tepicDirectory.toPath());
            } catch (URISyntaxException | IOException e) {
                throw new RuntimeException(e);
            }
        } else {
            // Usually entered when running via IDE
            logger.info("IDE mode");
            try {
                copyDirectory(new File(resource.getFile()), tepicDirectory,
                        pathname -> !pathname.getName().endsWith(".class"));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        makeAllChildrenExecutable(tepicDirectory);

        executable = extend(tepicDirectory, "TEPIC.sh");

        latestChipSeq.forEach((group, hmMap) -> {
            OutputFile dGroupOut = new OutputFile(outputDirectory, group);
            chipSeqFiles.put(group, new HashMap<>());
            outputFiles.put(group, new HashMap<>());
            OutputFile dGroupIn = new OutputFile(inputDirectory, group);
            hmMap.forEach((hm, samples) -> {
                OutputFile dHmOut = new OutputFile(dGroupOut, hm);
                chipSeqFiles.get(group).put(hm, new HashSet<>());
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile dHmIn = new OutputFile(dGroupIn, hm);
                samples.forEach(sample -> {
                    OutputFile sampleDir = addOutput(dHmOut,
                            sample.getName().substring(0, sample.getName().lastIndexOf('.')).replace(".", "_"));
                    outputFiles.get(group).get(hm).add(sampleDir);
                    InputFile inputFile = addInput(dHmIn, sample);
                    bridge.put(inputFile, sampleDir);
                    chipSeqFiles.get(group).get(hm).add(inputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        String order = "gbopScdnawfmervkiqhs";
        return new HashSet<>() {{
            File linkedReferenceGenome = new File(outputDirectory, referenceGenome.getName());
            try {
                softLink(linkedReferenceGenome, referenceGenome);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            chipSeqFiles.forEach((group, hmMap) -> hmMap.forEach((hm, samples) -> {
                samples.forEach(sample -> add(() -> {
                    Map<Character, String> stringConfigs = new HashMap<>() {{
                        put('b', sample.getAbsolutePath());
                        put('o', bridge.get(sample).getAbsolutePath() + "/");
                        put('S', new File(bridge.get(sample), sequenceFileName.get()).getAbsolutePath());
                        put('g', linkedReferenceGenome.getAbsolutePath());
                    }};

                    Map<Character, UsageConfig<?>> otherConfigs = new HashMap<>() {{
                        put('p', f_pwms);
                        put('d', bedChromatinSignal);
                        put('n', columnBedfile);
                        put('a', geneAnnotationFile);
                        put('w', windowSize);
                        put('f', onlyDNasePeaks);
                        put('e', exponentialDecay);
                        put('l', doNotNormalizePeakLength);
                        put('x', originalDecay);
                        put('m', psemsLengthFile);
                        put('y', entireGeneBody);
                        put('z', doZip);
                        put('r', twoBitFile);
                        put('v', pValue);
                        put('i', maxMinutesPerChromosome);
                        put('j', chromosomePrefix);
                        put('t', transcriptBased);
                        put('h', loopListFile);
                        put('s', loopWindows);
                        put('q', onlyPeakFeatures);
                    }};

                    otherConfigs.forEach((key, config) -> {
                        if (config.isSet()) {
                            if (config.get().getClass().equals(Boolean.class)) {
                                if (List.of('e', 'q').contains(key)) {
                                    // Boolean parameters
                                    stringConfigs.put(key, config.toString().toUpperCase());
                                } else {
                                    // Flags
                                    if ((boolean) config.get()) {
                                        stringConfigs.put(key, "");
                                    }
                                }
                            } else {
                                stringConfigs.put(key, config.get().toString());
                            }
                        }
                    });
                    String command = executable.getAbsolutePath() + " " +
                            order.chars().mapToObj(c -> (char) c).filter(stringConfigs::containsKey).map(c -> {
                                if (stringConfigs.containsKey(c)) {
                                    return "-" + c + " " + stringConfigs.get(c);
                                } else {
                                    return "";
                                }
                            }).collect(Collectors.joining(" "));

                    executeAndWait(command, true);
                    return true;
                }));
            }));
        }};
    }
}
