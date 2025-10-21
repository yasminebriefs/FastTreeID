test_that ("the program throws errors on invalid inputs", {
	directed <- c(1, 2, 3, 4)

	bidirected <- matrix(c(1,2 , 1,3, 1,4, 1,5), byrow = TRUE, ncol = 2)

	expect_no_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181"))

	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "288230376151711743")) # prime too small
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "576460752303423489")) # prime too large
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "47646075230342348a")) # prime invalid
	expect_error(fasttreeid_identify(bidirected, directed, "0x00900146909", "533993101691002181")) # seed invalid

	directed <- c(1, 2, 3, 5)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # parent too large

	directed <- c(0, 2, 3, 4)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # parent too small

	directed <- c(1, 2, 3, 4)

	bidirected <- matrix(c(2,2 , 1,3, 1,4, 1,5), byrow = TRUE, ncol = 2)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # u = v

	bidirected <- matrix(c(1,2 , 1,3, 1,4, 0,5), byrow = TRUE, ncol = 2)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # bidirected edge index u too small

	bidirected <- matrix(c(1,2 , 1,3, 1,4, 6,5), byrow = TRUE, ncol = 2)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # bidirected edge index u too large

	bidirected <- matrix(c(1,2 , 1,0, 1,4, 1,5), byrow = TRUE, ncol = 2)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # bidirected edge index v too small

	bidirected <- matrix(c(1,2 , 1,6, 1,4, 1,5), byrow = TRUE, ncol = 2)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # bidirected edge index v too large

	bidirected <- matrix(c(1,2,5), byrow = TRUE, ncol = 3)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # too many columns

	bidirected <- matrix(c(1, 2, 5), byrow = TRUE, ncol = 1)
	expect_error(fasttreeid_identify(bidirected, directed, "16630300901090146909", "533993101691002181")) # too few columns
})
