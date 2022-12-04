package org.exbio.tfprio.steps.rnaSeq;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class FilterTPM extends ExecutableStep {
    private final HashMap<InputFile, OutputFile> bridge;

    private final RequiredConfig<Double> tpmFilter = new RequiredConfig<>(Configs.deSeq2.tpmFilter);

    public FilterTPM(Collection<OutputFile> dependencies) {
        super(false, dependencies);

        bridge = dependencies.stream().collect(HashMap::new,
                (map, dependency) -> map.put(addInput(dependency), addOutput(dependency.getName())), HashMap::putAll);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            bridge.forEach((inputFile, outputFile) -> add(() -> {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                    List<String> inputLines = readLines(inputFile);
                    writer.write(inputLines.get(0));
                    writer.newLine();
                    inputLines.stream().skip(1).map(line -> line.split("\t")).forEachOrdered(split -> {
                        try {
                            writer.write(split[0] + "\t" + String.join("\t",
                                    Arrays.stream(split).skip(1).mapToDouble(Double::parseDouble).map(
                                            value -> value >= tpmFilter.get() ? value : 0).mapToObj(
                                            String::valueOf).toArray(String[]::new)));
                            writer.newLine();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    });
                }

                return true;
            }));
        }};
    }
}
