package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;
import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class DeSeq2 extends ExecutableStep<Configs> {
    public final Map<String, OutputFile> outputFiles = new HashMap<>();
    public final Map<String, OutputFile> normalizedCounts = new HashMap<>();
    private final InputFile script;
    private final InputFile batchFile;
    private final Set<String> pairings;
    private final Map<InputFile, OutputFile> filteredBatchFiles = new HashMap<>();
    private final RequiredConfig<Integer> countFilter = new RequiredConfig<>(configs.deSeq2.countFilter);
    private final Map<String, InputFile> inputFiles = new HashMap<>();

    public DeSeq2(Configs configs, Map<String, OutputFile> pairingCounts, OutputFile batchFile) {
        super(configs, false, pairingCounts.values(), batchFile);


        pairings = pairingCounts.keySet();

        pairings.forEach(pairing -> {
            OutputFile input = pairingCounts.get(pairing);
            InputFile inputFile = addInput(input);
            inputFiles.put(pairing, inputFile);

            OutputFile outputFile = addOutput(inputFile.getName());
            outputFiles.put(pairing, outputFile);

            OutputFile normalizedCount = addOutput(pairing + "_normalized.tsv");
            normalizedCounts.put(pairing, normalizedCount);

            filteredBatchFiles.put(inputFile,
                    addOutput(inputFile.getName().substring(0, inputFile.getName().lastIndexOf(".")) + "_meta.tsv"));
        });

        this.batchFile = addInput(batchFile);

        this.script = addInput(getClass().getResourceAsStream("deseq2.R"), "deseq2.R");
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            try {
                List<String> batchLines = readLines(batchFile);

                Map<String, String> sampleGroup = batchLines.stream().skip(1).map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line[1]));
                Map<String, String> sampleBatch = batchLines.stream().skip(1).map(line -> line.split("\t")).collect(
                        Collectors.toMap(line -> line[0], line -> line[2]));

                pairings.forEach(pairing -> add(() -> {
                    InputFile input = inputFiles.get(pairing);
                    OutputFile output = outputFiles.get(pairing);
                    OutputFile normalizedCount = normalizedCounts.get(pairing);

                    OutputFile filteredBatchFile = filteredBatchFiles.get(input);

                    final List<String> order;
                    try (BufferedReader reader = new BufferedReader(new FileReader(input))) {
                        order = Arrays.stream(reader.readLine().split("\t")).skip(1).toList();
                    }

                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(filteredBatchFile))) {
                        writer.write(batchLines.get(0));
                        writer.newLine();

                        order.forEach(sample -> {
                            try {
                                writer.write(
                                        String.join("\t", sample, sampleGroup.get(sample), sampleBatch.get(sample)));
                                writer.newLine();
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });
                    }

                    executeAndWait("Rscript " + script.getAbsolutePath() + " --metadata " +
                            filteredBatchFile.getAbsolutePath() + " --input " + input.getAbsolutePath() + " --output " +
                            output.getAbsolutePath() + " --out_normalized " + normalizedCount.getAbsolutePath() +
                            " --minCount " + countFilter.get(), true);
                    return true;
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }
}
