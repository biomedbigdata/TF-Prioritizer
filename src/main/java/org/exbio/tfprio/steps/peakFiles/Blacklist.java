package org.exbio.tfprio.steps.peakFiles;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.configs.ConfigTypes.UsageTypes.RequiredConfig;
import org.exbio.pipejar.pipeline.ExecutableStep;
import org.exbio.tfprio.configs.Configs;
import org.exbio.tfprio.lib.Region;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.readLines;

public class Blacklist extends ExecutableStep<Configs> {
    public final Map<String, Map<String, Collection<OutputFile>>> outputFiles = new HashMap<>();
    private final Map<InputFile, OutputFile> bridge = new HashMap<>();
    private final RequiredConfig<File> blacklistConfig = new RequiredConfig<>(configs.mixOptions.blackListPath);
    private final InputFile blacklist;

    public Blacklist(Configs configs, Map<String, Map<String, Collection<OutputFile>>> dependencies) {
        super(configs, false, dependencies.values().stream().flatMap(
                stringCollectionMap -> stringCollectionMap.values().stream()).flatMap(Collection::stream).collect(
                Collectors.toSet()));
        blacklist = addInput(blacklistConfig);

        dependencies.forEach((group, hmMap) -> {
            outputFiles.put(group, new HashMap<>());
            OutputFile d_groupOut = addOutput(group);
            hmMap.forEach((hm, sampleFiles) -> {
                outputFiles.get(group).put(hm, new HashSet<>());
                OutputFile d_hmOut = addOutput(d_groupOut, hm);

                sampleFiles.forEach(sampleFile -> {
                    InputFile inputFile = addInput(sampleFile);
                    OutputFile outputFile = addOutput(d_hmOut, inputFile.getName());
                    bridge.put(inputFile, outputFile);
                    outputFiles.get(group).get(hm).add(outputFile);
                });
            });
        });
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            try {
                logger.trace("Reading blacklist file");
                Collection<Region> blacklistRegions = readLines(blacklist).stream().map(Region::new).toList();
                logger.trace("Found " + blacklistRegions.size() + " blacklisted regions");

                bridge.forEach((inputFile, outputFile) -> add(() -> {
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                        // Read all regions from the input file
                        long count = readLines(inputFile).stream().map(Region::new)
                                                         // Filter out all blacklisted regions
                                                         .filter(region -> blacklistRegions.stream().noneMatch(
                                                                 region::overlaps))
                                                         // Get string representations of remaining regions
                                                         .map(Region::toString)
                                                         // Write them to the output file
                                                         .peek(line -> {
                                                             try {
                                                                 writer.write(line);
                                                                 writer.newLine();
                                                             } catch (IOException e) {
                                                                 throw new RuntimeException(e);
                                                             }
                                                         }).count();
                        logger.trace("Wrote " + count + " regions to " + outputFile);
                    }
                    return true;
                }));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }};
    }
}
