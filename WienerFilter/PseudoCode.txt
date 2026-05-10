// =========================================================================
// MAIN ROUTINE
// =========================================================================

INIT  input signal array, desired signal array, and output signal array
INIT  optimize coefficient array, autocorrelation array, and cross correlation
INIT  R matrix for Toeplitz 

SET   max supported size to 10
SET   float normalizer to 10.0

CALL read input file into input buffer from "input.txt"
CALL parse input buffer into input signal array
SET input count to parsed element count

CALL read desired file into desired buffer from "desired.txt"
CALL parse desired buffer into desired signal array
SET desired count to parsed element count

IF input count does not equal desired count THEN
    PRINT "Error: size not match"
    CALL  exit program
END IF

SET   signal length to input count
SET   filter length to input count

CALL  compute correlation
CALL  build Toeplitz matrix
CALL  gaussian elimination
CALL  compute output
CALL  compute MMSE
CALL  output results


// =========================================================================
// SUBROUTINE: Parse Input Buffer
// =========================================================================

SET   buffer pointer to start of input buffer
SET   array pointer to input signal array
SET   parsed count to 0

WHILE current buffer character is not null terminator

    CALL    parse float using current buffer pointer
    
    STORE   parsed float into input signal array
    
    UPDATE  buffer pointer to returned position
    
    MOVE    array pointer to next float slot
    
    INCREMENT parsed count

END WHILE

SAVE parsed count into input count


// =========================================================================
// SUBROUTINE: Parse Desired Buffer
// =========================================================================

SET   buffer pointer to start of desired buffer
SET   array pointer to desired signal array
SET   parsed count to 0

WHILE current buffer character is not null terminator

    CALL    parse float using current buffer pointer
    
    STORE   parsed float into desired signal array
    
    UPDATE  buffer pointer to returned position
    
    MOVE    array pointer to next float slot
    
    INCREMENT parsed count

END WHILE

SAVE parsed count into desired count


// =========================================================================
// SUBROUTINE: Compute Correlation
// =========================================================================

FOR each shift from 0 to filter length minus 1
    SET   sum xx to 0.0
    SET   sum dx to 0.0
    
    FOR each sample from shift to signal length minus 1
        COMPUTE sum xx by adding input signal at sample multiplied by 
                input signal at (sample minus shift) // n - k
                
        COMPUTE sum dx by adding desired signal at sample multiplied by 
                input signal at (sample minus shift) // n - k
    ENDFOR
    
    COMPUTE autocorrelation at shift by dividing sum xx by float normalizer
    COMPUTE cross correlation at shift by dividing sum dx by float normalizer
ENDFOR


// =========================================================================
// SUBROUTINE: Build Toeplitz Matrix
// =========================================================================

FOR each row on the matrix from 0 to filter length minus 1
    FOR each column on the matrix from 0 to filter length minus 1
        
        COMPUTE lag as the absolute difference between row and column
        SET     R matrix position (row, column) to autocorrelation at lag
        
    ENDFOR
ENDFOR


// =========================================================================
// SUBROUTINE: Gaussian Elimination
// =========================================================================

FOR each pivot row from 0 to filter length minus 1
    FOR each target row from pivot row plus 1 to filter length minus 1
        
        COMPUTE elimination factor by dividing R matrix position (target row, pivot row) 
                by R matrix position (pivot row, pivot row) // div.s   $f12, $f11, $f10       # factor = R[k][i] / R[i][i]
                
        COMPUTE cross correlation at target row by subtracting 
                (elimination factor multiplied by cross correlation at pivot row) // sub.s   $f14, $f14, $f15       # gamma_dx[k] -= factor * gamma_dx[i]
        
        FOR each target column from pivot row to filter length minus 1
            COMPUTE R matrix position (target row, target column) by subtracting 
                    (elimination factor multiplied by R matrix position (pivot row, target column)) // sub.s   $f17, $f17, $f18       # R[k][j] -= factor * R[i][j]
        ENDFOR
        
    ENDFOR
ENDFOR

FOR each row backwards from filter length - 1 to 0
    SET     back substitution sum to cross correlation at row
    
    FOR each column from row plus 1 to filter length minus 1
        COMPUTE back substitution sum by subtracting 
                (R matrix position (row, column) multiplied by optimize coefficient at column) // sub.s   $f0, $f0, $f3          # sum -= R_matrix[i][j] * h[j]
    ENDFOR
    
    COMPUTE optimize coefficient at row by dividing back substitution sum 
            by R matrix position (row, row) // div.s   $f5, $f0, $f4          # $f5 = sum / R_matrix[i][i]
ENDFOR


// =========================================================================
// SUBROUTINE: Compute Output
// =========================================================================

FOR each sample from 0 to signal length minus 1
    SET   sum y to 0.0
    
    FOR each filter tap from 0 to filter length minus 1
        IF sample minus filter tap is greater than or equal to 0 THEN
            COMPUTE sum y by adding optimize coefficient at filter tap 
                    multiplied by input signal at (sample minus filter tap) // add.s   $f0, $f0, $f3          # sum_y += $f3
        END IF
    ENDFOR
    
    SET   output signal at sample to sum y
ENDFOR


// =========================================================================
// SUBROUTINE: Compute MMSE
// =========================================================================

SET   sum error squared to 0.0

FOR each sample from 0 to signal length minus 1
    COMPUTE error by subtracting output signal at sample 
            from desired signal at sample // sub.s   $f3, $f1, $f2          # diff = d[n] - y[n]
            
    COMPUTE sum error squared by adding error multiplied by error // add.s   $f0, $f0, $f4          # error_sum += diff_squared
ENDFOR

COMPUTE final MMSE by dividing sum error squared by float normalizer // div.s   $f0, $f0, $f20         # final_mmse = error_sum / 10.0


// =========================================================================
// SUBROUTINE: Output Results
// =========================================================================

SET   output text buffer to "Filtered output: "
PRINT "Filtered output: "

FOR each sample from 0 to signal length minus 1
    CALL    format float to 1 decimal place with output signal at sample
    PRINT   formatted float
    COMPUTE output text buffer by adding formatted float
    
    IF sample is less than signal length minus 1 THEN
        PRINT   empty space
        COMPUTE output text buffer by adding empty space
    END IF
ENDFOR

PRINT   newline and "MMSE: "
COMPUTE output text buffer by adding newline and "MMSE: "

ROUND float to 1 decimal place manually
PRINT   formatted MMSE
COMPUTE output text buffer by adding formatted MMSE and newline

CALL    save output to file with "output.txt" and output text buffer
