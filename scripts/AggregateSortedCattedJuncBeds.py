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
        previous_line_name = None
        group_sum = 0

        for i, line in enumerate(sys.stdin):
            columns = line.strip().split('\t')

            # Extract relevant columns
            name = '_'.join([columns[0], columns[12], columns[13], columns[5]])
            score = int(columns[4])

            # Check if current line junc (name) is same as previous
            if previous_line_name and previous_line_name != name:
                previous_line_columns[3:5] = (previous_line_name, str(group_sum))
                print('\t'.join(previous_line_columns))
                group_sum = 0
                # Start processing the new group
            previous_line_name = name
            previous_line_columns = columns[0:12]
            group_sum += score

        # Output the sum for the last group
        previous_line_columns[3:5] = (previous_line_name, str(group_sum))
        print('\t'.join(previous_line_columns))
    except BrokenPipeError:
        # Ignore broken pipe error
        return

if __name__ == "__main__":
    process_input()

