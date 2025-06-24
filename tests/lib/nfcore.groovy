/**
 * Creates the modules directory and .nf-core.yml configuration file
 * @param libDir The directory path to initialise an nf-core library at
 */
def nfcore_setup(String libDir) {
    new File("${libDir}/modules").mkdirs()
    def nfcore_yml = new File("${libDir}/.nf-core.yml")
    nfcore_yml.write(
    """
    repository_type: "pipeline"
    template:
        name: test
    """.stripIndent()
    )
}

/**
 * Installs nf-core modules from a list
 * @param libDir An nf-core library initialised by nfcore_setup()
 * @param modules List of module names to install (e.g., ["minimap2/index", "samtools/view"])
 */
def nfcore_install(String libDir, List<String> modules) {
    modules.each { module ->
        def command = ['bash', '-c', "cd ${libDir} && nf-core modules install ${module}"]
        def command_result = command.execute()
        command_result.waitFor()
        if (command_result.exitValue() != 0) {
            println "Error installing module ${module}: ${command_result.err.text}"
        } else {
            println "Successfully installed module: ${module}"
        }
    }
}

/**
 * Creates a symbolic link from the installed nf-core modules to the base directory
 * @param libDir An nf-core library initialised by nfcore_setup()
 * @param modulesDir The root directory of an nf-core style modules repository
 */
def nfcore_link(String libDir, String modulesDir) {
    def sourceDir = new File("${libDir}/modules/nf-core")
    def targetDir = new File("${modulesDir}/modules/nf-core")

    try {
        java.nio.file.Files.createSymbolicLink(
            targetDir.toPath(),
            sourceDir.toPath()
        )
        println "Successfully created symlink: ${targetDir} -> ${sourceDir}"
    } catch (java.nio.file.FileAlreadyExistsException e) {
        println "Symlink already exists: ${targetDir}"
    } catch (Exception e) {
        println "Error creating symlink: ${e.message}"
        throw e
    }
}

/**
 * Cleanup function to remove created directories and files
 * @param libDir The launch directory path
 * @param modulesDir The base directory path
 */
def nfcore_cleanup(String libDir, String modulesDir) {
    new File("${libDir}/modules").deleteDir()
    new File("${libDir}/.nf-core.yml").delete()
    new File("${libDir}/modules.json").delete()
    new File("${modulesDir}/modules/nf-core").delete()
}
