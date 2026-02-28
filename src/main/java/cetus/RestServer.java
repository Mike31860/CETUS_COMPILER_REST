package cetus;
import cetus.mcp_server.CetusMcpTools;
import cetus.mcp_server.CetusMcpTools;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.ai.tool.ToolCallbackProvider;
import org.springframework.ai.tool.method.MethodToolCallbackProvider; // Note the .method subpackage
import org.springframework.context.annotation.Bean;

@SpringBootApplication
// This tells Spring to look for your Controller inside the cetus.controller package
@ComponentScan(basePackages = {"cetus.controller", "cetus.exec", "cetus.mcp_server"}) 
public class RestServer {
    public static void main(String[] args) {
        SpringApplication.run(RestServer.class, args);
    }

    @Bean
    public ToolCallbackProvider cetusToolProvider(CetusMcpTools cetusMcpTools) {
        // In 2.0.0-M2, use builder() then .toolObjects()
        return MethodToolCallbackProvider.builder()
                .toolObjects(cetusMcpTools)
                .build();
    }
}

