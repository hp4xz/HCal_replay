#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include <TCanvas.h>
#include <TText.h>

std::vector<std::string> wrapLine(const std::string& text, unsigned lineLength) {
    std::vector<std::string> result;
    std::istringstream iss(text);
    std::string word;
    
    std::string currentLine;
    while (iss >> word) {
        if (currentLine.length() + word.length() > lineLength) {
            result.push_back(currentLine);
            currentLine = word;
        } else {
            if (!currentLine.empty()) {
                currentLine += " ";
            }
            currentLine += word;
        }
    }
    if (!currentLine.empty()) {
        result.push_back(currentLine);
    }
    return result;
}

void txtToPdf(const char* inputFileName, const char* outputFileName) {
    // Open the text file
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Couldn't open the input file!" << std::endl;
        return;
    }

    // Determine the number of lines and wrap lines
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(inputFile, line)) {
        std::vector<std::string> wrapped = wrapLine(line, 60);  // wrap lines exceeding 60 characters
        lines.insert(lines.end(), wrapped.begin(), wrapped.end());
    }

    int numLines = lines.size();

    // Determine canvas size based on the number of lines
    int canvasHeight = 20 * numLines;  // assuming roughly 20 pixels per line, adjust if needed
    int canvasWidth = 800;  // default width

    // Create a canvas
    TCanvas c("c", "Text to PDF", canvasWidth, canvasHeight);

    // Use TText to handle the text drawing
    TText t;
    t.SetTextFont(42);  // Choose a font
    t.SetTextSize(20.0 / canvasHeight);  // Adjust the text size

    float yPosition = 0.98;  // Start from the top
    const float yDecrement = 2.0 / numLines;  // Line spacing

    for (const std::string& wrappedLine : lines) {
        t.DrawText(0.05, yPosition, wrappedLine.c_str());  // Draw text on the canvas
        yPosition -= yDecrement;
    }

    // Save the canvas to a PDF file
    c.Print(outputFileName, "pdf");

    inputFile.close();
}

int dbToPdf() {
    const char* inputFileName = "/w/halla-scshelf2102/sbs/seeds/sbsoffline/SBS-replay/DB/db_sbs.hcal.dat";
    const char* outputFileName = "hcal_db_prior_to_pass2_final.pdf";
    txtToPdf(inputFileName, outputFileName);
    return 0;
}
