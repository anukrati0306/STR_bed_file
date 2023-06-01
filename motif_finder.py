def find_repeating_unit(string, unit_length):
    frequency = {}
    unit_start = 0
    unit_end = unit_length

    while unit_end <= len(string):
        unit = string[unit_start:unit_end]

        if unit in frequency:
            frequency[unit] += 1
        else:
            frequency[unit] = 1

        unit_start += 1
        unit_end += 1

    repeating_units = [unit for unit, count in frequency.items() if count > 4]
    return repeating_units