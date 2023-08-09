
class demonstratorPartitionChecker:
    """
        Evaluate the quality of a partition
    """

    def getCoveredBits(self, bit_range, lst_sec, showFalse = False, id_range="bit range"):
        """
            Counts:
                - number of secret bits lieing in bit range
                - number of bit ranges affected
        """
        parts_used = {}
        num_found_bits = 0
        dict_bitRange_used = {}
        
        for j in range(len(lst_sec)):
            if lst_sec[j] == 0:
                continue   

            foundPart = False
            for r in bit_range:
                if j >= r[0] and j < r[1]:
                    num_found_bits += 1
                    foundPart = True
                    dict_bitRange_used[r[0]] = True
                    break

            if showFalse and not foundPart:
                print ("Estimates -- Secret position in " + str(id_range) + " not found: " + str(j))

        return int(num_found_bits), len(dict_bitRange_used)

    def getNumSampleSpace(self, bit_range):
        """
            Return size of sample space
        """
        size_ss = 0
        for r in bit_range:
            size_ss += (r[1] - r[0])
        return size_ss