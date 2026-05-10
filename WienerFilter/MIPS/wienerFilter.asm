.data
    # --- File paths and strings ---
    input_file:     .asciiz "input.txt"
    output_file:    .asciiz "output.txt"
    desired_file:   .asciiz "desired.txt"
    str_filtered:   .asciiz "Filtered output: "
    str_mmse_lbl:   .asciiz "\nMMSE: "
    str_space:      .asciiz " "
    str_newline:    .asciiz "\n"
    str_size_error: .asciiz "Error: size not match"

    # --- Floating point constants ---
    float_0:        .float 0.0
    float_0p5:      .float 0.5
    float_1:        .float 1.0
    float_10:       .float 10.0

    # --- System variables ---
    input_count:    .word 0
    desired_count:  .word 0
    buf_ptr:        .word 0
    const_N:        .word 10           
    const_M:        .word 10           

    # --- Result arrays ---
    desired_signal: .space 40
    input_signal:   .space 40
    output_signal:  .space 40    # y[n]
    optimize_coefficient: .space 40    # h[k]
    gamma_xx:       .space 40    # autocorrelation
    gamma_dx:       .space 40    # cross-correlation
    R_matrix:       .space 400   # Toeplitz matrix 10x10
    mmse:           .float 0.0   

    # --- Buffers for file I/O ---
    buf_write:      .space 1024
    input_buffer:   .space 1024
    desire_buffer:  .space 1024

.text
.global main

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
    
    jal     compute_correlation     # Step 2: gamma_xx, gamma_dx
    jal     build_toeplitz          # Step 3: matrix R_M
    jal     gaussian_elimination    # Step 4: Solve R * h = gamma_dx
    jal     compute_output          # Step 5: Calculate y[n] = sum(h*x)
    jal     compute_mmse            # Step 6: Calculate MMSE

    jal     print_results           # Print to Console
    jal     save_output             # Write to file output.txt

    li      $v0, 10
    syscall

# =========================================================
# PROCEDURES
# =========================================================

# --- Correlation Calculation ---
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
    la      $t0, input_signal
    sll     $t1, $s1, 2
    add     $t1, $t0, $t1
    l.s     $f0, 0($t1)
    sub     $t2, $s1, $s0
    sll     $t2, $t2, 2
    add     $t2, $t0, $t2
    l.s     $f1, 0($t2)     
    mul.s   $f2, $f0, $f1
    add.s   $f8, $f8, $f2
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

# --- Build Toeplitz matrix ---
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

# --- Solve linear system (Gaussian Elimination) ---
gaussian_elimination:
    addiu   $sp, $sp, -4
    sw      $ra, 0($sp)
    lw      $t8, const_M
    li      $s0, 0
ge_fwd_i:
    bge     $s0, $t8, ge_backward
    li      $t3, 10
    mul     $t4, $s0, $t3
    add     $t4, $t4, $s0
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f10, 0($t5)
    addiu   $s1, $s0, 1
ge_fwd_k:
    bge     $s1, $t8, ge_next_i
    mul     $t4, $s1, $t3
    add     $t4, $t4, $s0
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f11, 0($t5)
    div.s   $f12, $f11, $f10
    la      $t6, gamma_dx
    sll     $t7, $s0, 2
    add     $t7, $t6, $t7
    l.s     $f13, 0($t7)
    sll     $t0, $s1, 2
    add     $t0, $t6, $t0
    l.s     $f14, 0($t0)
    mul.s   $f15, $f12, $f13
    sub.s   $f14, $f14, $f15
    s.s     $f14, 0($t0)
    move    $s2, $s0
ge_fwd_j:
    bge     $s2, $t8, ge_next_k
    mul     $t4, $s0, $t3
    add     $t4, $t4, $s2
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f16, 0($t5)
    mul     $t4, $s1, $t3
    add     $t4, $t4, $s2
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f17, 0($t5)
    mul.s   $f18, $f12, $f16
    sub.s   $f17, $f17, $f18
    s.s     $f17, 0($t5)
    addiu   $s2, $s2, 1
    j       ge_fwd_j
ge_next_k:
    addiu   $s1, $s1, 1
    j       ge_fwd_k
ge_next_i:
    addiu   $s0, $s0, 1
    j       ge_fwd_i
ge_backward:
    li      $s0, 9
ge_back_i:
    bltz    $s0, ge_done
    la      $t6, gamma_dx
    sll     $t7, $s0, 2
    add     $t7, $t6, $t7
    l.s     $f0, 0($t7)
    addiu   $s1, $s0, 1
ge_back_sum:
    bge     $s1, $t8, ge_back_div
    li      $t3, 10
    mul     $t4, $s0, $t3
    add     $t4, $t4, $s1
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f1, 0($t5)
    la      $t0, optimize_coefficient
    sll     $t1, $s1, 2
    add     $t1, $t0, $t1
    l.s     $f2, 0($t1)
    mul.s   $f3, $f1, $f2
    sub.s   $f0, $f0, $f3
    addiu   $s1, $s1, 1
    j       ge_back_sum
ge_back_div:
    li      $t3, 10
    mul     $t4, $s0, $t3
    add     $t4, $t4, $s0
    sll     $t4, $t4, 2
    la      $t5, R_matrix
    add     $t5, $t5, $t4
    l.s     $f4, 0($t5)
    div.s   $f5, $f0, $f4
    la      $t0, optimize_coefficient
    sll     $t1, $s0, 2
    add     $t1, $t0, $t1
    s.s     $f5, 0($t1)
    addiu   $s0, $s0, -1
    j       ge_back_i
ge_done:
    lw      $ra, 0($sp)
    addiu   $sp, $sp, 4
    jr      $ra

# --- Calculate output y[n] ---
compute_output:
    lw      $t8, const_N
    li      $s0, 0
co_n_loop:
    bge     $s0, $t8, co_done
    l.s     $f0, float_0
    li      $s1, 0
co_k_loop:
    bge     $s1, $t8, co_store
    sub     $t0, $s0, $s1
    bltz    $t0, co_k_next
    la      $t1, optimize_coefficient
    sll     $t2, $s1, 2
    add     $t1, $t1, $t2
    l.s     $f1, 0($t1)
    la      $t3, input_signal
    sll     $t4, $t0, 2
    add     $t3, $t3, $t4
    l.s     $f2, 0($t3)
    mul.s   $f3, $f1, $f2
    add.s   $f0, $f0, $f3
co_k_next:
    addiu   $s1, $s1, 1
    j       co_k_loop
co_store:
    la      $t5, output_signal
    sll     $t6, $s0, 2
    add     $t5, $t5, $t6
    s.s     $f0, 0($t5)
    addiu   $s0, $s0, 1
    j       co_n_loop
co_done:
    jr      $ra

# --- Calculate MMSE ---
compute_mmse:
    lw      $t8, const_N
    l.s     $f20, float_10
    l.s     $f0, float_0
    li      $s0, 0
cm_loop:
    bge     $s0, $t8, cm_store
    la      $t0, desired_signal
    sll     $t1, $s0, 2
    add     $t0, $t0, $t1
    l.s     $f1, 0($t0)
    la      $t2, output_signal
    add     $t2, $t2, $t1
    l.s     $f2, 0($t2)
    sub.s   $f3, $f1, $f2
    mul.s   $f4, $f3, $f3
    add.s   $f0, $f0, $f4
    addiu   $s0, $s0, 1
    j       cm_loop
cm_store:
    div.s   $f0, $f0, $f20
    s.s     $f0, mmse
    jr      $ra

# --- Print Result Procedure ---
print_results:
    la      $t0, buf_write
    sw      $t0, buf_ptr
    addiu   $sp, $sp, -4
    sw      $ra, 0($sp)

    li      $v0, 4
    la      $a0, str_filtered
    syscall
    
    la      $a0, str_filtered
    jal     append_string
    
    lw      $t8, const_N
    li      $s0, 0
pr_loop:
    bge     $s0, $t8, pr_mmse
    la      $t0, output_signal
    sll     $t1, $s0, 2
    add     $t0, $t0, $t1
    l.s     $f12, 0($t0)
    jal     print_float_1dp
    addiu   $t2, $s0, 1
    bge     $t2, $t8, pr_next
    
    # Print space
    li      $v0, 11
    li      $a0, ' '
    syscall

    lw      $t0, buf_ptr
    sb      $a0, 0($t0)
    addiu   $t0, $t0, 1
    sw      $t0, buf_ptr
pr_next:
    addiu   $s0, $s0, 1
    j       pr_loop
pr_mmse:
    li      $v0, 4
    la      $a0, str_mmse_lbl
    syscall
    
    la      $a0, str_mmse_lbl
    jal     append_string
    
    l.s     $f12, mmse
    jal     print_float_1dp
    
    li      $v0, 4
    la      $a0, str_newline
    syscall
    
    la      $a0, str_newline
    jal     append_string
    
    lw      $t0, buf_ptr
    sb      $zero, 0($t0)
    
    lw      $ra, 0($sp)
    addiu   $sp, $sp, 4
    jr      $ra

# --- Print Float with 1 decimal place ---
print_float_1dp:
    addiu   $sp, $sp, -4
    sw      $ra, 0($sp)
    l.s     $f0, float_0
    li      $t9, 0
    c.lt.s  $f12, $f0
    bc1f    pf1_pos
    li      $t9, 1
    neg.s   $f12, $f12
pf1_pos:
    l.s     $f1, float_10
    mul.s   $f2, $f12, $f1
    l.s     $f3, float_0p5
    add.s   $f2, $f2, $f3
    cvt.w.s $f4, $f2
    mfc1    $t0, $f4
    li      $t1, 10
    div     $t0, $t1
    mflo    $t2
    mfhi    $t3
    beqz    $t9, pf1_print
    
    li      $v0, 11
    li      $a0, '-'
    syscall
    lw      $t0, buf_ptr
    sb      $a0, 0($t0)
    addiu   $t0, $t0, 1
    sw      $t0, buf_ptr
pf1_print:
    li      $v0, 1
    move    $a0, $t2
    syscall
    
    addiu   $a0, $t2, 48
    lw      $t0, buf_ptr
    sb      $a0, 0($t0)
    addiu   $t0, $t0, 1
    sw      $t0, buf_ptr
    
    li      $v0, 11
    li      $a0, '.'
    syscall
    
    lw      $t0, buf_ptr
    sb      $a0, 0($t0)
    addiu   $t0, $t0, 1
    sw      $t0, buf_ptr
    
    li      $v0, 1
    move    $a0, $t3
    syscall
    
    addiu   $a0, $t3, 48
    lw      $t0, buf_ptr
    sb      $a0, 0($t0)
    addiu   $t0, $t0, 1
    sw      $t0, buf_ptr
    
    lw      $ra, 0($sp)
    addiu   $sp, $sp, 4
    jr      $ra

# --- Write Result to File ---
save_output:
    addiu   $sp, $sp, -4
    sw      $ra, 0($sp)

    li      $v0, 13
    la      $a0, output_file
    li      $a1, 1
    li      $a2, 0
    syscall

    move    $s7, $v0

    li      $v0, 15
    move    $a0, $s7
    la      $a1, buf_write

    la      $t0, buf_write
    li      $t1, 0
len_loop:
    lb      $t2, 0($t0)
    beqz    $t2, len_done
    addiu   $t1, $t1, 1
    addiu   $t0, $t0, 1
    j       len_loop
len_done:
    move    $a2, $t1
    
    li      $v0, 15
    syscall

    li      $v0, 16
    move    $a0, $s7
    syscall

    lw      $ra, 0($sp)
    addiu   $sp, $sp, 4
    jr      $ra

# --- Helper: Append string to buffer ---
append_string:
    lw      $t0, buf_ptr
as_loop:
    lb      $t1, 0($a0)
    beqz    $t1, as_done
    sb      $t1, 0($t0)
    addiu   $t0, $t0, 1
    addiu   $a0, $a0, 1
    j       as_loop
as_done:
    sw      $t0, buf_ptr
    jr      $ra

# --- File Reading Procedures ---
read_input_file:
    li      $v0, 13
    la      $a0, input_file
    li      $a1, 0
    li      $a2, 0
    syscall
    move    $s6, $v0
    li      $v0, 14
    move    $a0, $s6
    la      $a1, input_buffer
    li      $a2, 1024
    syscall
    la      $t0, input_buffer
    add     $t0, $t0, $v0
    sb      $zero, 0($t0)
    li      $v0, 16
    move    $a0, $s6
    syscall
    jr      $ra

read_desire_file:
    li      $v0, 13
    la      $a0, desired_file
    li      $a1, 0
    li      $a2, 0
    syscall
    move    $s6, $v0
    li      $v0, 14
    move    $a0, $s6    
    la      $a1, desire_buffer
    li      $a2, 1024
    syscall    
    la      $t0, desire_buffer
    add     $t0, $t0, $v0
    sb      $zero, 0($t0)
    li      $v0, 16    
    move    $a0, $s6    
    syscall
    jr      $ra

# --- Parsing Procedures ---
parse_float:
    li      $t7, 0              
pf_skip:
    lb      $t0, 0($a0)
    li      $t1, ' '
    beq     $t0, $t1, pf_skip_next
    li      $t1, '\r'
    beq     $t0, $t1, pf_skip_next
    li      $t1, '\n'
    beq     $t0, $t1, pf_skip_next
    j       pf_check_sign
pf_skip_next:
    addiu   $a0, $a0, 1
    j       pf_skip
pf_check_sign:
    li      $t1, '-'
    bne     $t0, $t1, pf_integer_start
    li      $t7, 1
    addiu   $a0, $a0, 1
pf_integer_start:
    li      $t2, 0
pf_integer_loop:
    lb      $t0, 0($a0)
    beqz    $t0, pf_done
    li      $t1, '.'
    beq     $t0, $t1, pf_decimal
    addiu   $t0, $t0, -48
    mul     $t2, $t2, 10
    add     $t2, $t2, $t0
    addiu   $a0, $a0, 1
    j       pf_integer_loop
pf_decimal:
    addiu   $a0, $a0, 1
    lb      $t0, 0($a0)
    addiu   $t0, $t0, -48
    move    $t3, $t0
    addiu   $a0, $a0, 1

    mtc1    $t2, $f1
    cvt.s.w $f1, $f1
    mtc1    $t3, $f2
    cvt.s.w $f2, $f2
    l.s     $f3, float_10
    div.s   $f2, $f2, $f3
    add.s   $f0, $f1, $f2
    
    beqz    $t7, pf_done
    neg.s   $f0, $f0
pf_done:
    move    $v0, $a0
    jr      $ra

parse_desire_buffer:
    addiu   $sp, $sp, -16
    sw      $ra, 0($sp)
    sw      $s0, 4($sp)
    sw      $s1, 8($sp)
    sw      $s2, 12($sp)

    la      $s0, desire_buffer
    la      $s1, desired_signal
    li      $s2, 0
pdb_loop:
    lb      $t0, 0($s0)
    beqz    $t0, pdb_done
    move    $a0, $s0
    jal     parse_float
    s.s     $f0, 0($s1)
    move    $s0, $v0
    addiu   $s1, $s1, 4
    addiu   $s2, $s2, 1
    j       pdb_loop
pdb_done:
    sw      $s2, desired_count
    lw      $ra, 0($sp)
    lw      $s0, 4($sp)
    lw      $s1, 8($sp)
    lw      $s2, 12($sp)
    addiu   $sp, $sp, 16
    jr      $ra

parse_input_buffer:
    addiu   $sp, $sp, -16
    sw      $ra, 0($sp)
    sw      $s0, 4($sp)
    sw      $s1, 8($sp)
    sw      $s2, 12($sp)

    la      $s0, input_buffer
    la      $s1, input_signal
    li      $s2, 0
pib_loop:
    lb      $t0, 0($s0)
    beqz    $t0, pib_done
    move    $a0, $s0
    jal     parse_float
    s.s     $f0, 0($s1)
    move    $s0, $v0
    addiu   $s1, $s1, 4
    addiu   $s2, $s2, 1
    j       pib_loop
pib_done:
    sw      $s2, input_count
    lw      $ra, 0($sp)
    lw      $s0, 4($sp)
    lw      $s1, 8($sp)
    lw      $s2, 12($sp)
    addiu   $sp, $sp, 16
    jr      $ra

# --- Error Handling ---
size_error:
    li      $v0, 4
    la      $a0, str_size_error
    syscall
    li      $v0, 10
    syscall
