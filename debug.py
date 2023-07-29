def format_sample_fields(format_field, sample_field):
    format_field = format_field.split(':')
    output = {}
    for sample_id, sample_info in sample_field.items():
        sample_info_list = sample_info.split(':')
        sample_info_dict = {}
        for i, info in enumerate(sample_info_list):
            if i >= len(format_field):
                break
            field_name = format_field[i]
            sample_info_dict[field_name] = info
        output[sample_id] = sample_info_dict
    print(output)
    return output

