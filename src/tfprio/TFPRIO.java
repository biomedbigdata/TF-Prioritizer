package tfprio;

import lib.Blacklist.Blacklist;
import lib.CheckChromosomes;
import lib.CreateDeseq2Scripts;
import lib.CreateTpmMappings;
import lib.FootprintIntervals.CreateFootprintsBetweenPeaks;
import lib.MixOptions.MixMutuallyExclusive;
import lib.MixOptions.MixOptions;
import org.apache.commons.cli.*;

import util.Configs.Configs;
import util.MapSymbolAndEnsg;
import util.Options_intern;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;


public class TFPRIO
{
    public static Configs configs;
    public static Map<String, Set<String>> groupsToHms = new HashMap<>();
    public static MapSymbolAndEnsg mapSymbolAndEnsg;

    public static void main(String[] args) throws Exception
    {
        long startTime = System.currentTimeMillis();
        Options_intern options_intern = new Options_intern();
        parseArguments(args, options_intern);

        configs = new Configs(new File(options_intern.com2pose_working_directory),
                new File(options_intern.path_to_COM2POSE));

        configs.merge(new File(options_intern.path_to_COM2POSE + File.separator + "config_templates" + File.separator +
                "configsTemplate.json"));

        configs.validate();

        new CheckChromosomes();

        //mix histone modifications
        if (configs.mixOptions.level.isSet())
        {
            new MixOptions();
        }

        //create footprints if needed
        if (configs.tepic.tfBindingSiteSearch.isSet() && (configs.tepic.tfBindingSiteSearch.get().equals("BETWEEN") ||
                configs.tepic.tfBindingSiteSearch.get().equals("EXCL_BETWEEN")))
        {
            new CreateFootprintsBetweenPeaks();
        }


        //black listed regions filter
        if (configs.blacklist.bedFilePath.isSet())
        {
            new Blacklist();
        }

        if (configs.mixOptions.mutuallyExclusive.get())
        {
            new MixMutuallyExclusive();
        }


        //DESeq2
        if ((!configs.tepic.ensgSymbolFile.isSet() || configs.deSeq2.biomartDatasetSpecies.isSet()) &&
                configs.general.ensgMappingEnabled.get())
        {
            mapSymbolAndEnsg = new MapSymbolAndEnsg();
        }

        /*
        Kann BatchVariablen akzeptieren
        Wenn entsprechende Config gesetzt ist, soll File übergeben werden, wo dann samples Batches zugeordnet sind.
        Diese Zuordnungen sollen dann DESeq2 übergeben werden
        */


        new CreateDeseq2Scripts();

        new CreateTpmMappings();
        /*
        com2pose_lib.create_gene_positions();

        if (options_intern.deseq2_tpm_filter > 0)
        {
            com2pose_lib.preprocess_deseq2_input_tpm();
        }
        com2pose_lib.run_and_postprocess_DESeq2();


        if (!options_intern.path_tgen.equals(""))
        {
            com2pose_lib.preprocess_tgen();
            com2pose_lib.run_tgen();
            com2pose_lib.merge_tgen();
        }


        //TEPIC
        com2pose_lib.run_tepic();
        if (options_intern.tepic_randomize_tf_gene_matrix)
        {
            com2pose_lib.randomize_tepic();
        }

        com2pose_lib.postprocess_tepic_output();
        com2pose_lib.create_open_regions_violin_plots();

        //TGen
        if (!options_intern.path_tgen.equals(""))
        {
            if (!options_intern.mix_mutually_exclusive)
            {
                com2pose_lib.create_tgen_groups();
            }

            com2pose_lib.filter_target_genes_tgen();

            if (options_intern.tgen_self_regulatory)
            {
                com2pose_lib.integrate_self_regulatory_tgen();
            }
        }

        //DYNAMITE
        com2pose_lib.preprocess_dynamite();
        //install Rscripts needed for DYNAMITE
        com2pose_lib.install_required_packages();
        com2pose_lib.run_DYNAMITE();

        //PLOTS
        com2pose_lib.create_tp_plots();
        com2pose_lib.analyze_plots_data();
        com2pose_lib.get_top_k_target_genes_plots();


        //PREPARE DISTR ANALYSIS
        com2pose_lib.prepare_distribution_analysis();

        //DISTRIBUTION ANALYSIS
        com2pose_lib.perform_distribution_analysis();
        com2pose_lib.create_distribution_plots();
        //com2pose_lib.create_overview_html_report_distribution();
        com2pose_lib.calculate_discounted_cumulative_gain_rank_distribution_analysis();

        //get target genes of TFs based on dcg analysis
        com2pose_lib.get_top_k_target_genes_dcg();

        //generate DCG target genes heatmap
        com2pose_lib.generateHeatmaps();

        //create energy motifs of dcg genes
        com2pose_lib.create_tf_binding_logos_biophysical_tfSequence();
        com2pose_lib.create_tf_binding_logos_predicted_binding_sites();

        //check for ChIP-seq data for TFs on ChIP-Atlas
        if (options_intern.chip_atlas_activated_chip_atlas)
        {
            com2pose_lib.get_chip_atlas_data_list();
            com2pose_lib.get_chip_atlas_data();
            //com2pose_lib.run_igv_chip_atlas_data();
            //com2pose_lib.run_igv_chr_wide_data();

        }

        //shoot IGV for own TF ChIP-seq data
        if (!options_intern.igv_path_to_igv.equals("") && !options_intern.igv_path_to_tf_chip_seq.equals(""))
        {
            //com2pose_lib.run_igv_own_data();
        }

        //shot IGV for dcg target gene heatmaps
        com2pose_lib.run_igv_on_dcg_targetgenes_heatmaps();

        if (options_intern.igv_important_locus_all_prio_tf.size() > 0)
        {
            com2pose_lib.run_igv_on_important_loci();
        }

        if (options_intern.igv_top_log2fc > 0)
        {
            com2pose_lib.run_igv_top_log2gc();
        }

        //run co-occurence analysis of ChIP-ATLAS .bed files
        if (options_intern.chip_atlas_activated_chip_atlas)
        {
            com2pose_lib.run_cooccurence_analysis();
        }

        new Report(options_intern).generate();
        */
        double deltaSeconds = (double) (System.currentTimeMillis() - startTime) / 1e3;
        System.out.println("COM2POSE finished. Execution took " + deltaSeconds + " seconds.");
    }


    private static void parseArguments(String[] args, Options_intern options_intern)
    {
        Options options = new Options();

        Option opt_cfg_file = new Option("c", "com2pose-config", true,
                "[REQ]: COM2POSE config file. Example in /COM2POSE/config_templates/com2pose_template.cfg");
        opt_cfg_file.setRequired(true);
        options.addOption(opt_cfg_file);

        Option opt_working_dir = new Option("w", "working-directory", true,
                "[REQ]: working directory where COM2POSE can create, remove and edit files");
        opt_working_dir.setRequired(true);
        options.addOption(opt_working_dir);

        Option opt_path_to_COM2POSE = new Option("p", "path-com2pose", true, "[REQ]: filepath to COM2POSE folder");
        opt_path_to_COM2POSE.setRequired(true);
        options.addOption(opt_path_to_COM2POSE);

        Option opt_tgen = new Option("t", "tgen-dir", true,
                "[OPT]: Use a consensus approach of TGen and TEPIC, provide directory to MEME suite installation folder here");
        options.addOption(opt_tgen);

        Option opt_write_logfile = new Option("l", "write-log-file", false,
                "[OPT]: If flag is set no logfile will be written, default: logfile will be written");
        options.addOption(opt_write_logfile);

        Option opt_do_ensg_mapping = new Option("m", "do-ensg-mapping", false,
                "[OPT]: If flag is set no ensg mapping will be done (only use when biomart error)");
        options.addOption(opt_do_ensg_mapping);

        Option opt_do_tpm_length_calculation = new Option("a", "do-tpm-length-calculation", false,
                "[OPT]: If flag is set no tpm length calculation will be done (only use when biomart error)");
        options.addOption(opt_do_tpm_length_calculation);

        Option opt_do_position_calculation = new Option("b", "do-position-calculations", false,
                "[OPT]: If flag is set no gene position retrieval via biomart will be done (only use when biomart error)");
        options.addOption(opt_do_position_calculation);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;


        try
        {
            cmd = parser.parse(options, args);


            options_intern.config_data_path = cmd.getOptionValue("com2pose-config");
            options_intern.com2pose_working_directory = cmd.getOptionValue("working-directory");
            options_intern.path_to_COM2POSE = cmd.getOptionValue("path-com2pose");

            if (cmd.hasOption("tgen-dir"))
            {
                options_intern.path_tgen = cmd.getOptionValue("tgen-dir");
            }

            if (cmd.hasOption("write-log-file"))
            {
                options_intern.write_to_logfile = false;
            }
            if (cmd.hasOption("do-ensg-mapping"))
            {
                options_intern.do_ensg_mapping = false;
            }
            if (cmd.hasOption("do-tpm-length-calculation"))
            {
                options_intern.calculate_tpm_lengths = false;
            }
            if (cmd.hasOption("do-position-calculations"))
            {
                options_intern.calculcate_gene_positions = false;
            }


        } catch (ParseException e)
        {
            System.out.println(e.getMessage());
            formatter.printHelp(
                    "-c <com2pose-config> -w <working-directory> -p <path-com2pose> [-t <tgen-dir>] [-l] [-m] [-a] [-t] [-b]",
                    options);
            System.exit(1);
        }
    }
}