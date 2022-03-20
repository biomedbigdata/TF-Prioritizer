package util;

import util.Configs.Config;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class FileManagement
{
    public static String readFile(File file) throws IOException
    {
        return Files.readString(file.toPath());
    }

    public static List<String> readLines(File file) throws IOException
    {
        return Files.readAllLines(file.toPath());
    }

    public static void copyFile(File source, File target) throws IOException
    {
        copyFile(source, target, false);
    }

    public static void copyFile(File source, File target, boolean compression) throws IOException
    {
        if (target.exists())
        {
            return;
        }

        if (source.exists())
        {
            if (!target.getParentFile().exists())
            {
                target.getParentFile().mkdirs();
            }

            if (target.getName().endsWith(".png") && compression)
            {
                String command =
                        "pngtopnm " + source.getAbsolutePath().replace(" ", "\\ ") + " | pnmquant 16 | pnmtopng >" +
                                " " + target.getAbsolutePath().replace(" ", "\\ ");
                String[] cmd = {"/bin/sh", "-c", command};
                Process child = Runtime.getRuntime().exec(cmd);
                try
                {
                    child.waitFor(5, TimeUnit.SECONDS);
                } catch (InterruptedException ignored)
                {
                }
            } else
            {
                Files.copy(source.toPath(), target.toPath(), REPLACE_EXISTING);
            }
        }
    }

    public static void copyDirectory(File source, File target, boolean compression)
    {
        copyDirectory(source, target, compression, ".*", new ArrayList<>());
    }

    public static void copyDirectory(File source, File target, boolean compression, String fileNameRegex,
                                     List<String> removables)
    {
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        for (File sourceFile : Objects.requireNonNull(source.listFiles()))
        {
            if (sourceFile.isFile())
            {
                if (sourceFile.getName().matches(fileNameRegex))
                {
                    String cleanName = sourceFile.getName();
                    String fileExtension = cleanName.substring(cleanName.lastIndexOf("."));

                    cleanName = cleanName.replace(fileExtension, "");

                    for (String removable : removables)
                    {
                        cleanName = cleanName.replace(removable, "");
                    }

                    while (cleanName.contains("__"))
                    {
                        cleanName = cleanName.replace("__", "_");
                    }
                    while (cleanName.startsWith("_"))
                    {
                        cleanName = cleanName.substring(1);
                    }
                    while (cleanName.endsWith("_"))
                    {
                        cleanName = cleanName.substring(0, cleanName.length() - 1);
                    }
                    cleanName = cleanName + fileExtension;
                    File targetFile = new File(target.getAbsolutePath() + File.separator + cleanName);
                    executorService.execute(() ->
                    {
                        try
                        {
                            copyFile(sourceFile, targetFile, compression);
                        } catch (IOException ignored)
                        {
                        }
                    });
                }
            } else
            {
                File subDir = new File(target.getAbsolutePath() + File.separator + sourceFile.getName());
                copyDirectory(sourceFile, subDir, compression, fileNameRegex, removables);
            }
        }
        try
        {
            executorService.shutdown();
            executorService.awaitTermination(30, TimeUnit.SECONDS);
        } catch (InterruptedException e)
        {
            System.out.println("[REPORT] Error during copy process: " + source.getAbsolutePath() + " to " +
                    target.getAbsolutePath());
        }
    }

    public static void writeHTML(File file, String content, int relativationDepth) throws IOException
    {
        content = content.replace("{RELATIVATION}", (".." + File.separator).repeat(relativationDepth));

        writeFile(file, content);
    }

    public static void writeFile(File file, String content) throws IOException
    {
        makeSureFileExists(file);
        if (file.isFile())
        {
            Files.writeString(file.toPath(), content);
        } else
        {
            System.out.println(file.getAbsolutePath());
            throw new IllegalArgumentException("Can only process files, not directories!");
        }
    }

    public static String findValueInTable(String term, int searchIndex, int resultIndex, File file, String sep,
                                          boolean ignoreCase) throws FileNotFoundException, NoSuchFieldException
    {
        if (file.isFile())
        {
            try (Scanner scanner = new Scanner(file))
            {
                while (scanner.hasNextLine())
                {
                    String line = scanner.nextLine();
                    if (line.split(sep).length > searchIndex && line.split(sep).length > resultIndex)
                    {

                        if (term.equals(line.split(sep)[searchIndex]) ||
                                ignoreCase && term.equalsIgnoreCase(line.split(sep)[searchIndex]))
                        {
                            return line.split(sep)[resultIndex];
                        }
                    }
                }
            }
        }
        throw new NoSuchFieldException();
    }

    public static File getFileIfInDirectory(File directory, String fileNameRegex, boolean lookingForFiles)
    {
        if (directory == null)
        {
            return null;
        }

        if (directory.listFiles() == null)
        {
            return null;
        }

        for (File entry : directory.listFiles())
        {
            if (lookingForFiles != entry.isFile())
            {
                continue;
            }
            if (entry.getName().matches(fileNameRegex))
            {
                return entry;
            }
        }
        return null;
    }

    public static synchronized void makeSureFileExists(File file) throws IOException
    {
        if (!file.exists())
        {
            try
            {
                if (!file.getParentFile().exists())
                {
                    if (!file.getParentFile().mkdirs())
                    {
                        throw new IOException("parent directory");
                    }
                }
                if (!file.createNewFile())
                {
                    throw new IOException("file");
                }
            } catch (IOException e)
            {
                throw new IOException("Exception during " + e.getMessage() + " creation: " + file.getAbsolutePath());
            }
        }
    }

    public static synchronized void makeSureDirectoryExists(File directory) throws IOException
    {
        if (directory.isFile())
        {
            throw new IllegalArgumentException("Can not handle files. Received: " + directory.getAbsolutePath());
        }
        if (!directory.exists())
        {
            if (!directory.mkdirs())
            {
                throw new IOException("Could not create directory: " + directory.getAbsolutePath());
            }
        }
    }

    public static void appendToFile(File file, String content) throws IOException
    {
        makeSureFileExists(file);

        try
        {
            Files.write(file.toPath(), content.getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e)
        {
            throw new IOException("Could not append content to file: " + file.getAbsolutePath());
        }
    }

    public static Config<File> extend(Config<File> fileConfig, String extension)
    {
        return new Config<>(extend(fileConfig.get(), extension));
    }

    public static File extend(File file, String... extensions)
    {
        File extended = new File(file.getAbsolutePath());

        for (String extension : extensions)
        {
            extended = new File(extended.getAbsolutePath() + File.separator + extension);
        }
        return extended;
    }

    public static String hashFileContent(File file) throws IOException
    {
        String content = readFile(file);
        return String.valueOf(content.hashCode());
    }

    public static void createSourceFile(File file, String hash) throws IOException
    {
        File sourceFile = getSourceFile(file);
        writeFile(sourceFile, hash);
    }

    public static boolean hasValidSourceFile(File file, String hash)
    {
        File sourceFile = getSourceFile(file);

        if (file.exists() && sourceFile.exists() && sourceFile.isFile() && sourceFile.canRead())
        {
            try
            {
                if (Files.readAttributes(file.toPath(), BasicFileAttributes.class).lastModifiedTime().compareTo(
                        Files.readAttributes(sourceFile.toPath(), BasicFileAttributes.class).lastModifiedTime()) > 0)
                {
                    // File has been modified since last generation
                    return false;
                }
                if (readFile(sourceFile).equals(hash))
                {
                    return true;
                }
            } catch (IOException e)
            {
                return false;
            }
        }
        return false;
    }

    public static boolean isEmpty(File file) throws IOException
    {
        String content = readFile(file);
        return content.isBlank();
    }

    private static File getSourceFile(File file)
    {
        return extend(file.getParentFile(), "." + file.getName());
    }
}
