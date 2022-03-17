package lib;

import tfprio.TFPRIO;
import util.FileFilters.DirectoryFilter;
import util.FileFilters.FileFilter;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

import static util.FileManagement.*;

public class CheckChromosomes extends ExecutableStep
{
    /**
     * Make sure that chromosomes are annotated without a "chr" prefix
     */
    public void execute()
    {
        File file_root_input = TFPRIO.configs.tepic.inputDirectory.get();
        File chr_annotation_output = TFPRIO.configs.mixOptions.fileStructure.d_preprocessingCheckChr.get();

        // I do not really know why this assignments have to take place
        //options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        //options_intern.tepic_input_directory = chr_annotation_output.getAbsolutePath();

        for (File fileDir : Objects.requireNonNull(file_root_input.listFiles(new DirectoryFilter())))
        {
            String group = fileDir.getName();
            TFPRIO.groupsToHms.put(group, new HashSet<>());

            File d_outputGroup = extend(chr_annotation_output, group);

            for (File d_hm : Objects.requireNonNull(fileDir.listFiles(new DirectoryFilter())))
            {
                TFPRIO.groupsToHms.get(group).add(d_hm.getName());
                File d_outputGroupHm = extend(d_outputGroup, d_hm.getName());

                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(new FileFilter())))
                {
                    executorService.execute(() ->
                    {
                        List<String> inputLines = null;
                        try
                        {
                            inputLines = readLines(f_sample);
                        } catch (IOException e)
                        {
                            logger.error("Could not read content from " + f_sample.getAbsolutePath());
                            System.exit(1);
                        }
                        StringBuilder sb_output = new StringBuilder();

                        for (String inputLine : inputLines)
                        {
                            String processedLine;
                            if (inputLine.startsWith("chr"))
                            {
                                processedLine = inputLine.substring("chr".length());
                            } else
                            {
                                processedLine = inputLine;
                            }
                            sb_output.append(processedLine);
                            sb_output.append("\n");
                        }
                        File f_outputGroupHmSample = extend(d_outputGroupHm, f_sample.getName());

                        try
                        {
                            writeFile(f_outputGroupHmSample, sb_output.toString());
                        } catch (IOException e)
                        {
                            logger.error("Could not write content to " + f_outputGroupHmSample);
                            System.exit(1);
                        }
                    });

                }
            }
        }
        try
        {
            TFPRIO.configs.general.latestInputDirectory.setValue(chr_annotation_output);
        } catch (IllegalAccessException e)
        {
            logger.error(e.getMessage());
        }
    }
}