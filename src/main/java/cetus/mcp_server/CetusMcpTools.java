package cetus.mcp_server;

import cetus.exec.Driver;
// NEW (Use these)
import org.springframework.ai.tool.annotation.Tool;
import org.springframework.ai.tool.annotation.ToolParam;
import org.springframework.stereotype.Component;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@Component
public class CetusMcpTools {

    @Tool(name = "optimize_c_code", 
          description = "Uses the Cetus compiler to analyze and optimize C source code for parallelization.")
    public String runCetusOptimization(
        @ToolParam(description = "The specific Cetus flags to use") List<String> flags,
        @ToolParam(description = "The path to the C file to be analyzed") String filePath
    ) {
        try {
            List<String> argsList = new ArrayList<>(flags);
            argsList.add(filePath);
            
            System.out.println("Running MCP server Cetus with the following arguments: " + String.join(" ", argsList));
            String[] args = argsList.toArray(new String[0]);

            // Call the same Driver you used in your Controller
            Driver.main(args);

            return "Cetus successfully processed " + filePath + " with flags: " + String.join(" ", flags);
        } catch (Exception e) {
            return "Cetus Error: " + e.getMessage();
        }
    }

    @Tool(name = "profile_application", 
          description = "profile the application using tools.")
          public String profileApplication(@ToolParam(description = "The path to the C binary file to be analyzed") String filePath, @ToolParam(description = "This is where the profile infomration will be stored at") String locationresult) throws IOException {
           
            File binaryfile = new File(filePath);
            File workingDirectory = new File(binaryfile.getParent());

            File destinationDirec = new File(locationresult);
            if(!destinationDirec.exists())
            {
                boolean created=destinationDirec.mkdir();
                if(!created){
                    throw new IOException("We could no create the directory");
                }
            }

            File outputFile = new File(destinationDirec, "profile_output.txt");
            
            ProcessBuilder processBuilder = new ProcessBuilder("perf", "stat", filePath);
            processBuilder.directory(workingDirectory);
            processBuilder.redirectOutput(outputFile);

            Process process = processBuilder.start();

            int exitCode;
            try {
                exitCode = process.waitFor();
                 if (exitCode != 0) {
                    return "Profiling failed. Check output file: " + outputFile.getAbsolutePath();
                }

                return "Profiling completed successfully. Results stored at: "
                        + outputFile.getAbsolutePath();

            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
                 return "There was an error trying to profile the application";

            }

          }

}