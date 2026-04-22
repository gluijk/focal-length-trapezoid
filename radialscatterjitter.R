# Radial scatter: statistical exercise confirming the law of large numbers
# Originally to add jitter noise in focal length/aspect ratio calculations, but it's too concentrated
# www.overfitting.net
# https://www.overfitting.net/2026/03/tres-micro-ejercicios-de-vibe-coding.html


library(Cairo)

radial_scatter <- function(R = 1, n = 200000, bins = 60, alpha = 0.03, seed = 1){
    
    # --- Generate samples ---
    set.seed(seed)
    r <- runif(n, 0, R)
    theta <- runif(n, 0, 2*pi)
    x <- r * cos(theta)
    y <- r * sin(theta)
    
    name=sprintf("scatter_%06d.png", n)
    CairoPNG(name, width=512, height=800)
        op <- par(mfrow = c(2,1))
    
        # --- 1. 2D density map ---
        # Scatterplot
        plot(
            x, y, asp = 1,
            pch = 16, cex = 0.3,
            col = rgb(1, 0, 0, alpha),
            xlim = c(-R, R), ylim = c(-R, R),
            xlab = "x", ylab = "y",
            main = paste0("Scatter distribution (", as.integer(n), " samples)")
        )
        
        # Draw boundary circle
        t <- seq(0, 2*pi, length.out = 1000)
        lines(R*cos(t), R*sin(t), lwd = 1, col='gray')
        
        
        # --- 2. Radial density estimation ---
        edges <- seq(0, R, length.out = bins + 1)
        counts <- hist(r, breaks = edges, plot = FALSE)$counts
        
        r_mid <- (edges[-1] + edges[-length(edges)]) / 2
        ring_area <- pi * (edges[-1]^2 - edges[-length(edges)]^2)
        density_area <- counts / ring_area
        density_area <- density_area / max(density_area) # normalize for visual comparison
        
        # Theoretical 1/r curve (scaled)
        theory <- 1 / r_mid
        theory <- theory / max(theory)  # normalize for visual comparison
        
        plot(r_mid, density_area, col='red', type="l", xlim=c(0, 1), ylim=c(0.01, 1),
             xlab="Radius", ylab="Normalized spatial density (log)",
             main="Radial density (area corrected)",
             log="y")  # y-axis on log scale
        lines(r_mid, theory, lty=2)
        legend("topright", legend=c("Measured","1/r theory"), lty=c(1,2))
        
        par(op)
    dev.off()
    
    invisible(data.frame(x=x, y=y, r=r))
}

SEED=7
radial_scatter(R = 1, n = 200, alpha = 0.5, seed = SEED)
radial_scatter(R = 1, n = 1000, alpha = 0.5, seed = SEED)
radial_scatter(R = 1, n = 5000, alpha = 0.25, seed = SEED)
radial_scatter(R = 1, n = 25000, alpha = 0.15, seed = SEED)
radial_scatter(R = 1, n = 100000, alpha = 0.05, seed = SEED)
radial_scatter(R = 1, n = 400000, alpha = 0.01, seed = SEED)


