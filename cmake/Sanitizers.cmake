function (add_address_sanitizer TARGET_NAME)
    target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=address)
    target_link_options(${TARGET_NAME} PRIVATE -fsanitize=address)
endfunction ()

function (add_undefined_sanitizer TARGET_NAME)
    target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=undefined)
    target_link_options(${TARGET_NAME} PRIVATE -fsanitize=undefined)
endfunction ()
