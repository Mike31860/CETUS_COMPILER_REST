package cetus;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.ComponentScan;

@SpringBootApplication
// This tells Spring to look for your Controller inside the cetus.controller package
@ComponentScan(basePackages = {"cetus.controller", "cetus.exec"}) 
public class RestServer {
    public static void main(String[] args) {
        SpringApplication.run(RestServer.class, args);
    }
}