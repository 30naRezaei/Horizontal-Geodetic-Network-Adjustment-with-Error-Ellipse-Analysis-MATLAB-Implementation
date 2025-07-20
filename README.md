Horizontal Geodetic Network Adjustment with Error Ellipse Analysis – MATLAB Implementation
This repository presents a comprehensive MATLAB code for 2D (horizontal) network adjustment of a geodetic control network using the least squares method, enriched with precision analysis and error ellipse computation for each estimated station.

 Key Features
✅ Observation Processing and Accuracy Estimation
•	Calculates the mean observational accuracy based on the nominal precision of the theodolite (total station) and the field measurements.
•	Constructs the variance-covariance matrix of the observations and derives the weight matrix accordingly.
✅ Network Adjustment and Statistical Analysis
•	Performs least squares adjustment to estimate the unknown station coordinates.
•	Determines the degrees of freedom of the adjustment.
•	Computes the residuals (corrections) of observations and estimates the accuracy of adjusted variables (e.g., coordinates, directions).
•	Evaluates the posterior variance factor and performs a statistical test on it to verify the reliability of the network adjustment.
✅ Error Ellipses and Precision Visualization
•	Computes and draws error ellipses for each network station based on the covariance matrix of adjusted coordinates.
•	Calculates the area of each error ellipse, providing a quantitative representation of positional precision.
•	Notes that larger elongation of an ellipse along a specific direction indicates greater uncertainty in that direction.
✅ Visualization
•	Graphically represents the survey network, including stations, connections, and associated error ellipses using MATLAB plotting functions.
________________________________________
🗂 Input & Assumptions
•	Input includes observed distances and azimuths between stations for a defined 2D network.
•	Each dataset ends with a given azimuth (bearing) for one edge to resolve orientation.
•	Station 1 is fixed with known coordinates (x₁ = 1000, y₁ = 1000) and serves as the control point.
________________________________________
💻 Technical Details
•	Programming Language: MATLAB
•	Adjustment Method: Constrained Least Squares
•	Output: Adjusted coordinates, residuals, variance factor, error ellipses, and ellipse areas
________________________________________
🔍 Applications
•	Teaching geodetic network adjustment and precision analysis
•	Professional surveying quality control
•	Positional uncertainty visualization in 2D geodetic networks

