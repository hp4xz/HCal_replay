#include <TCanvas.h>
#include <TLatex.h>
#include <TText.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <sstream>

// Define the structs
struct CommonData {
    std::string headerText;
    double someValue;
    // Add other members as needed
};

struct MyData {
    std::string text;
    double value1;
    int value2;
};

void DrawStringsOnCanvas(const std::vector<MyData>& data, const CommonData& commonData) {
    // Create a new ROOT file to save the canvases
    TFile *file = new TFile("output.root", "RECREATE");

    for (size_t i = 0; i < data.size(); ++i) {
        // Create a new canvas for each MyData struct
        TCanvas *canvas = new TCanvas(Form("c%zu", i+1), Form("Canvas %zu", i+1), 600, 400);
        canvas->cd();

        double yPos = 0.9; // Starting Y position for text, top of the canvas
        double yStep = 0.05; // Step size for each line of text

        // Display common data
        TLatex *headerLatex = new TLatex(0.1, yPos, commonData.headerText.c_str());
        headerLatex->SetTextColor(kRed); // Set the text color to red
        headerLatex->SetTextSize(0.03);
        headerLatex->Draw();
        yPos -= yStep;

        std::stringstream ss;
        ss << "Value: " << commonData.someValue;
        TLatex *commonValueLatex = new TLatex(0.1, yPos, ss.str().c_str());
        commonValueLatex->SetTextColor(kRed);
        commonValueLatex->SetTextSize(0.03);
        commonValueLatex->Draw();
        yPos -= yStep;

        // Now display the data from the MyData struct
        ss.str("");
        ss << "String: " << data[i].text;
        TLatex *myDataTextLatex = new TLatex(0.1, yPos, ss.str().c_str());
        myDataTextLatex->SetTextSize(0.03);
        myDataTextLatex->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Double: " << data[i].value1;
        TLatex *myDataDoubleLatex = new TLatex(0.1, yPos, ss.str().c_str());
        myDataDoubleLatex->SetTextSize(0.03);
        myDataDoubleLatex->Draw();
        yPos -= yStep;

        ss.str("");
        ss << "Int: " << data[i].value2;
        TLatex *myDataIntLatex = new TLatex(0.1, yPos, ss.str().c_str());
        myDataIntLatex->SetTextSize(0.03);
        myDataIntLatex->Draw();

        // Write the canvas to the file
        canvas->Write();
        canvas->Update();
    }

    // Close the file
    file->Close();
}

void scanv() {
    CommonData commonData = {"Common Header", 123.45};

    std::vector<MyData> data = {
        {"String 1", 1.1, 2},
        {"String 2", 3.3, 4},
        {"String 3", 5.5, 6}
        // Add more MyData structs as needed
    };

    DrawStringsOnCanvas(data, commonData);
}
