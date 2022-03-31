package lib.Dynamite;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class RunDynamite extends ExecutableStep
{
    private final Config<File> f_dynamite_script = TFPRIO.configs.tepic.f_dynamite;
    private final Config<File> d_input = TFPRIO.configs.dynamite.fileStructure.d_preprocessing_prepareClassification;

    private final Config<File> d_output = TFPRIO.configs.dynamite.fileStructure.d_output;

    private final Config<String> outVar = TFPRIO.configs.dynamite.outVar;
    private final Config<Integer> oFolds = TFPRIO.configs.dynamite.oFolds;
    private final Config<Integer> iFolds = TFPRIO.configs.dynamite.iFolds;
    private final Config<Boolean> performance = TFPRIO.configs.dynamite.performance;
    private final Config<Double> alpha = TFPRIO.configs.dynamite.alpha;
    private final Config<Integer> cores = TFPRIO.configs.dynamite.cores;
    private final Config<Boolean> randomize = TFPRIO.configs.dynamite.randomize;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_dynamite_script, d_input));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(outVar, oFolds, iFolds, performance, alpha, cores, randomize));
    }

    @Override protected void execute()
    {
        for (Map.Entry<String, Set<String>> pairingEntry : TFPRIO.groupCombinationsToHms.entrySet())
        {
            String pairing = pairingEntry.getKey();

            for (String hm : pairingEntry.getValue())
            {
                File d_in = extend(d_input.get(), pairing, hm);
                int count_line_br = 0;

                for (File f_x : Objects.requireNonNull(d_in.listFiles(Filters.fileFilter)))
                {
                    try
                    {
                        count_line_br += readLines(f_x).size();
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }

                if (count_line_br < 2)
                {
                    continue;
                }

                File d_out = extend(d_output.get(), pairing, hm);
                makeSureDirectoryExists(d_out, logger);

                String command_edited = "Rscript " + f_dynamite_script.get().getAbsolutePath() + " --dataDir=" +
                        d_in.getAbsolutePath() + " --outDir=" + d_out.getAbsolutePath() + " --out_var=" + outVar.get() +
                        " --Ofolds=" + oFolds.get() + " --Ifolds=" + iFolds.get() + " --performance=" +
                        performance.get() + " --alpha=" + alpha.get() + " --cores=" + cores.get() + " --randomise=" +
                        randomize.get();

                executorService.submit(() ->
                {
                    try
                    {
                        executeAndWait(command_edited, logger);
                    } catch (IOException | InterruptedException e)
                    {
                        logger.error(e.getMessage());
                    }
                });
            }
        }
    }
}
