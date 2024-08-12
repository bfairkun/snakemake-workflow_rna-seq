#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : AggregateSortedCattedJuncBeds
# @created     : Monday Mar 18, 2024 20:58:25 CDT
#
# @description : 
######################################################################

import sys

def process_input():
    try:
        current_group = None
        group_sum = 0

        for i, line in enumerate(sys.stdin):
            columns = line.strip().split('\t')

            # Extract relevant columns
            column1 = columns[0]
            column6 = columns[5]
            column13 = int(columns[12])
            column14 = int(columns[13])
            column5 = int(columns[4])

            # Check if current line belongs to the same group or a new one
            if current_group is None or (column1, column6, column13, column14) != current_group:
                # If it's a new group, output the sum of column5 for the previous group
                if current_group is not None:
                    print('\t'.join(previous_group + [str(group_sum)]))
                # Start processing the new group
                current_group = (column1, column6, column13, column14)
                group_sum = 0
                previous_group = columns
            group_sum += column5

        # Output the sum for the last group
        if current_group is not None:
            print('\t'.join(previous_group + [str(group_sum)]))
    except BrokenPipeError:
        # Ignore broken pipe error
        return

if __name__ == "__main__":
    process_input()

