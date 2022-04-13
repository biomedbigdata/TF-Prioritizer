package tfprio;

import util.Configs.Config;
import util.Configs.Configs;
import util.ExecutionTimeMeasurement;
import util.Logger;
import util.MapSymbolAndEnsg;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static util.FileManagement.extend;


public class TFPRIO
{
    /**
     * The {@link Configs} object which is referenced by all the {@link lib.ExecutableStep} instances
     */
    public static Configs configs;

    /**
     * Mapping the existing groups to their associated histone modifications.
     */
    public static Map<String, Set<String>> groupsToHms = new HashMap<>();

    /**
     * Mapping the existing group combinations to their associated histone modifications.
     */
    public static Map<String, Set<String>> groupCombinationsToHms = new HashMap<>();

    /**
     * Contains all the file configs whose allocated files or directories will be created during pipeline execution.
     * <p>
     * Is being filled during simulation.
     * In development mode, also existing file structures are considered as created by the pipeline. Otherwise,
     * commenting out executableSteps inside the {@link Workflow} would lead to a simulator error.
     */
    public static Set<Config<File>> createdFileStructure = new HashSet<>();

    /**
     * Allows disk read efficient mapping of geneSymbols, geneIDs and their descriptions.
     */
    public static MapSymbolAndEnsg mapSymbolAndEnsg;


    /**
     * Stores the file config pointing to the latest modification of the input files.
     * <p>
     * Will be replaced by the {@link Workflow} getLatestInputDirectory method.
     */
    @Deprecated public static Config<File> latestInputDirectory;

    /**
     * Contains all the histone modifications available in the input data.
     */
    public static Set<String> existingHms = new HashSet<>();

    static File workingDirectory;
    static File sourceDirectory;
    static File configFile;

    public static void main(String[] args) throws Exception
    {
        ExecutionTimeMeasurement timer = new ExecutionTimeMeasurement();
        ArgParser.parseArguments(args);

        configs = new Configs(workingDirectory, sourceDirectory);

        configs.merge(extend(sourceDirectory, "config_templates", "configsTemplate.json"));
        configs.merge(configFile);
        configs.validate();

        Logger logger = new Logger("TFPRIO");

        Workflow workflow = new Workflow();

        if (workflow.simulationSuccessful())
        {
            workflow.run();
        }

        logger.info("Finished. Execution took " + timer.stopAndGetDeltaFormatted() + ".");
    }
}
