.data
    input_file:        .asciiz "input.txt"
    output_file:       .asciiz "output.txt"
    desired_file:      .asciiz "desired.txt"

    str_filtered:      .asciiz "Filtered output: "
    str_mmse_lbl:      .asciiz "\nMMSE: "
    str_space:         .asciiz " "
    str_newline:       .asciiz "\n"
    str_size_error:    .asciiz "Error: size not match"

    float_0:           .float 0.0
    float_0p5:         .float 0.5
    float_1:           .float 1.0
    float_10:          .float 10.0

    input_count:       .word 0
    desired_count:     .word 0
    buf_ptr:           .word 0

    const_N:           .word 10
    const_M:           .word 10

    # --- Result arrays ---
    desired_signal:        .space 40
    input_signal:          .space 40
    output_signal:         .space 40      # y[n]
    optimize_coefficient:  .space 40      # h[k]

    gamma_xx:              .space 40      # autocorrelation
    gamma_dx:              .space 40      # cross-correlation

    R_matrix:              .space 400     # Toeplitz matrix 10x10
    mmse:                  .float 0.0

    # --- Buffers ---
    buf_write:         .space 1024
    input_buffer:      .space 1024
    desire_buffer:     .space 1024


.text
.globl main

main:
    jal     read_input_file
    jal     parse_input_buffer

    jal     read_desire_file
    jal     parse_desire_buffer

    lw      $t0, input_count
    lw      $t1, desired_count

    bne     $t0, $t1, size_error

    sw      $t0, const_N
    sw      $t0, const_M

    # Step 1: Compute correlations
    jal     compute_correlation

    # Step 2: Build Toeplitz matrix
    jal     build_toeplitz

    # Step 3: Solve R * h = gamma_dx
    jal     gaussian_elimination

    # Step 4: Compute output signal
    jal     compute_output

    # Step 5: Compute MMSE
    jal     compute_mmse

    # Step 6: Print and save result
    jal     print_results
    jal     save_output

    li      $v0, 10
    syscall


# =========================================================
# Compute Correlation
# =========================================================
compute_correlation:
    addiu   $sp, $sp, -4
    sw      $ra, 0($sp)

    lw      $t8, const_N
    lw      $t9, const_M

    l.s     $f20, float_10

    li      $s0, 0

cc_k_loop:
    bge     $s0, $t9, cc_done

    l.s     $f8, float_0
    l.s     $f9, float_0

    move    $s1, $s0

cc_n_loop:
    bge     $s1, $t8, cc_store

    # x[n]
    la      $t0, input_signal
    sll     $t1, $s1, 2
    add     $t1, $t0, $t1
    l.s     $f0, 0($t1)

    # x[n-k]
    sub     $t2, $s1, $s0
    sll     $t2, $t2, 2
    add     $t2, $t0, $t2
    l.s     $f1, 0($t2)

    mul.s   $f2, $f0, $f1
    add.s   $f8, $f8, $f2

    # d[n]
    la      $t3, desired_signal
    sll     $t4, $s1, 2
    add     $t4, $t3, $t4
    l.s     $f3, 0($t4)

    mul.s   $f4, $f3, $f1
    add.s   $f9, $f9, $f4

    addiu   $s1, $s1, 1
    j       cc_n_loop

cc_store:
    div.s   $f8, $f8, $f20
    div.s   $f9, $f9, $f20

    la      $t5, gamma_xx
    sll     $t6, $s0, 2
    add     $t5, $t5, $t6
    s.s     $f8, 0($t5)

    la      $t5, gamma_dx
    add     $t5, $t5, $t6
    s.s     $f9, 0($t5)

    addiu   $s0, $s0, 1
    j       cc_k_loop

cc_done:
    lw      $ra, 0($sp)
    addiu   $sp, $sp, 4
    jr      $ra


# =========================================================
# Build Toeplitz Matrix
# =========================================================
build_toeplitz:
    lw      $t8, const_M

    li      $s0, 0

bt_l_loop:
    bge     $s0, $t8, bt_done

    li      $s1, 0

bt_k_loop:
    bge     $s1, $t8, bt_next_row

    sub     $t0, $s0, $s1

    bgez    $t0, bt_abs_done
    neg     $t0, $t0

bt_abs_done:
    la      $t1, gamma_xx

    sll     $t2, $t0, 2
    add     $t1, $t1, $t2

    l.s     $f0, 0($t1)

    li      $t3, 10

    mul     $t4, $s0, $t3
    add     $t4, $t4, $s1

    sll     $t4, $t4, 2

    la      $t5, R_matrix
    add     $t5, $t5, $t4

    s.s     $f0, 0($t5)

    addiu   $s1, $s1, 1
    j       bt_k_loop

bt_next_row:
    addiu   $s0, $s0, 1
    j       bt_l_loop

bt_done:
    jr      $ra
