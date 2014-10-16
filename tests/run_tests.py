import unittest
from matrix_builder_test import *
from state_freqs_test import *

if __name__ == '__main__':

    run_tests = unittest.TextTestRunner()
    print "Running unit tests for pyvolve"
    
    ##################### MATRIX_BUILDER TESTS #######################
    print "Testing assemble_matrix function of matrixBuilder for codon model"
    test_suite_assemble_matrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_assemble_matrix_tests)
    run_tests.run(test_suite_assemble_matrix)

    print "Testing the functions in the base class matrixBuilder"
    test_suite_baseMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_baseClass_tests)
    run_tests.run(test_suite_baseMatrix)
    
    print "Testing mechCodon_Matrix, a subclass of the parent matrixBuilder"
    test_suite_mechCodonMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_mechCodon_Matrix_tests)
    run_tests.run(test_suite_mechCodonMatrix)

    print "Testing aminoAcids_Matrix, a subclass of the parent matrixBuilder"
    test_suite_aminoAcidMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_aminoAcid_Matrix_tests)
    run_tests.run(test_suite_aminoAcidMatrix)

    print "Testing ECM_Matrix, a subclass of the parent matrixBuilder"
    test_suite_empCodonMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_ECM_Matrix_tests)
    run_tests.run(test_suite_empCodonMatrix)
 
    print "Testing nucleotide_Matrix, a subclass of the parent matrixBuilder"
    test_suite_nucleotideMatrix = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_nucleotide_Matrix_tests)
    run_tests.run(test_suite_nucleotideMatrix)

    print "Testing mutSel_Matrix, a subclass of the parent matrixBuilder, when unit of evolution is codons."
    test_suite_mutSelMatrix_codon = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_mutSel_codon_Matrix_tests)
    run_tests.run(test_suite_mutSelMatrix_codon)

    print "Testing mutSel_Matrix, a subclass of the parent matrixBuilder, when unit of evolution is nucleotides."
    test_suite_mutSelMatrix_nuc = unittest.TestLoader().loadTestsFromTestCase(matrixBuilder_mutSel_nuc_Matrix_tests)
    run_tests.run(test_suite_mutSelMatrix_nuc) 
    
    ############################## STATE_FREQS TESTS ###########################
    
    print "Testing the RandomFrequencies subclass of StateFrequencies"
    test_suite_Rand = unittest.TestLoader().loadTestsFromTestCase(state_freqs_RandomFrequencies_Tests)
    run_tests.run(test_suite_Rand)

    print "Testing the EqualFrequencies subclass of StateFrequencies"
    test_suite_Equal = unittest.TestLoader().loadTestsFromTestCase(state_freqs_EqualFrequencies_Tests)
    run_tests.run(test_suite_Equal)

    print "Testing the CustomFrequencies subclass of StateFrequencies"
    test_suite_User = unittest.TestLoader().loadTestsFromTestCase(state_freqs_CustomFrequencies_Tests)
    run_tests.run(test_suite_User)

    print "Testing the ReadFrequencies subclass of StateFrequencies"
    test_suite_Read = unittest.TestLoader().loadTestsFromTestCase(state_freqs_ReadFrequencies_Tests)
    run_tests.run(test_suite_Read)
