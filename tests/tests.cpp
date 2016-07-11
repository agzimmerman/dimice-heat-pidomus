#include <gtest/gtest.h>
#include "../src/step-14.h"
#include "test_tools.h"
std::string TEST_DIR;
TEST(Step14, ResultsMatch) {
	std::cout << std::endl << "Running step14." << std::endl;
    step14_main(); // This creates the output file which will be verified.
    std::string file_name = "log.txt";
	std::cout << "Comparing " << file_name << " to " << TEST_DIR+"/"+file_name << std::endl;
    EXPECT_TRUE(are_files_equal(file_name, TEST_DIR+"/"+file_name));
}
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    assert(argc == 2);
    TEST_DIR = argv[1];
	std::cout << "TEST_DIR = " << TEST_DIR << std::endl;
    return RUN_ALL_TESTS();
}