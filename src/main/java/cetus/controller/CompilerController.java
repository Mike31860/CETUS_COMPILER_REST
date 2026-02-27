package cetus.controller;


import org.springframework.web.bind.annotation.*;
import cetus.exec.Driver;
import java.util.List;

@RestController
@RequestMapping("/api")
public class CompilerController {

    @GetMapping("/hello")
    public String sayHello() {
        return "The API is responding at /api/hello!";
    }

    // This endpoint will receive a list of command line arguments for Cetus
    @PostMapping("/run-cetus")
    public String runCetus(@RequestBody List<String> commandArgs) {
        try {
            // 1. Convert the List from JSON into a String array
            // This is exactly what the command line used to provide to main()
            String[] args = commandArgs.toArray(new String[0]);

            // 2. Call the Cetus Driver directly
            // Driver.main will handle the parsing of -alias, -omp2gpu, and the file list
            Driver.main(args);

            return "Cetus execution finished with arguments: " + String.join(" ", args);
        } catch (Exception e) {
            return "Execution Error: " + e.getMessage();
        }
    }
}