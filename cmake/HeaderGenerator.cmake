# HeaderGenerator.cmake Generates a formatted header 
# *****************************************
# *            <DESCRIPTION>              *
# *          Version <VERSION>            *
# * Copyright (C) <YEAR> <NAME> <<EMAIL>> *
# *****************************************
# and saves to the output variable. Formatting is done to fit the width of the lines.
# Usage is:
# generate_header(DESCRIPTION description
#   VERSION version
#   YEAR year
#   NAME name
#   EMAIL email
#   OUTPUT_VARIABLE output_variable)

function(generate_header)
    # Initialize variables to empty
    set(description "")
    set(version "")
    set(year "")
    set(name "")
    set(email "")
    set(output_variable "")

    # Parse key-value arguments
    while(ARGC GREATER 0)
        list(GET ARGV 0 key)
        list(GET ARGV 1 value)
        list(REMOVE_AT ARGV 0 1)
        math(EXPR ARGC "${ARGC} - 2")

        if(key STREQUAL "DESCRIPTION")
            set(description "${value}")
        elseif(key STREQUAL "VERSION")
            set(version "${value}")
        elseif(key STREQUAL "YEAR")
            set(year "${value}")
        elseif(key STREQUAL "NAME")
            set(name "${value}")
        elseif(key STREQUAL "EMAIL")
            set(email "${value}")
        elseif(key STREQUAL "OUTPUT_VARIABLE")
            set(output_variable "${value}")
        else()
            message(FATAL_ERROR "Unknown key: ${key}")
        endif()
    endwhile()

    # Ensure all required fields are set
    if(NOT description)
        message(FATAL_ERROR "DESCRIPTION is required")
    endif()

    if(NOT version)
        message(FATAL_ERROR "VERSION is required")
    endif()

    if(NOT year)
        message(FATAL_ERROR "YEAR is required")
    endif()

    if(NOT name)
        message(FATAL_ERROR "NAME is required")
    endif()

    if(NOT email)
        message(FATAL_ERROR "EMAIL is required")
    endif()

    if(NOT output_variable)
        message(FATAL_ERROR "OUTPUT_VARIABLE is required")
    endif()

    # Create the lines for the header
    set(lines
        "${description}"
        "Version ${version}"
        "Copyright (C) ${year} ${name} <${email}>"
    )

    # Find the max length of all lines
    set(max_length 0)
    foreach(line IN LISTS lines)
        string(LENGTH "${line}" line_length)
        if(line_length GREATER max_length)
            set(max_length ${line_length})
        endif()
    endforeach()

    # Add spaces to center the lines
    math(EXPR max_length "${max_length} + 2") # Add 2 for padding spaces
    math(EXPR max_length_plus2 "${max_length} + 2") # Add 2 for additonall star
    set(star_line "")
    string(REPEAT "*" ${max_length_plus2} star_line)

    set(formatted_lines)
    list(APPEND formatted_lines "${star_line}")

    foreach(line IN LISTS lines)
        string(LENGTH "${line}" line_length)
        math(EXPR left_spaces "(${max_length} - ${line_length}) / 2")
        math(EXPR right_spaces "${max_length} - ${line_length} - ${left_spaces}")

        # Add spaces and '*' to the left and right
        string(REPEAT " " ${left_spaces} left_padding)
        string(REPEAT " " ${right_spaces} right_padding)
        set(formatted_line "*${left_padding}${line}${right_padding}*")

        list(APPEND formatted_lines "${formatted_line}")
    endforeach()

    list(APPEND formatted_lines "${star_line}")

    # Now join the list with newlines
    string(REPLACE ";" ";" formatted_lines "${formatted_lines}")
    string(JOIN "\n" header ${formatted_lines})

    # Save the header into the output_variable variable
    set(${output_variable} "${header}" PARENT_SCOPE)

endfunction()