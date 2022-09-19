    program PipelineFeatureMatching2
    include 'link_fnl_static.h'
    !DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
    use linear_operators ! use IMSL linear operator functions
    implicit none

    ! Fixed Variables
    integer, parameter :: nX = 62 ! Number of features from ILI 1
    integer, parameter :: nY = 37 ! Number of features from ILI 2
    integer, parameter :: d = 3 ! Coordinate dimensions
    character (len=*), parameter :: fileName1 = 'EE_ILI1_031.csv'
    character (len=*), parameter :: fileName2 = 'EE_ILI2_037.csv'
    double precision, parameter :: annealingRate = 0.96d0 ! Annealing rate
    double precision, parameter :: zeta = 0.0001d0
    double precision, parameter :: multiplier_lambda2 = 200.d0
    double precision :: T_final ! Final temperature to stop the annealing process

    ! Variables
    integer :: i, j, counter
    double precision :: lambda2
    double precision :: T, T_old ! Temperature in annealing process
    double precision :: T_init ! Initial temperature to start the annealing process
    double precision :: T_outlier ! Temperature for outliers
    double precision, dimension(nX, d) :: x ! Coordinates of ILI 1 features
    double precision, dimension(nY, d) :: y ! Coordinates of ILI 2 features
    double precision, dimension(nX, nY) :: distFeatPairs ! Distance of all feature pairs 
    double precision, dimension(nX+1, nY+1) :: p, p_old ! Correspondence matrix
    double precision, dimension(d, d) :: affine ! Affine transformation matrix
    double precision, dimension(nX, d) :: x_trans ! Transformed x-coordinates
    double precision, dimension(2000, 2) :: E ! Energy values
    double precision :: start_time, end_time, elapsed_time ! Variables to determine the duration of the analysis
    logical :: stopAnalysis

    ! Set start time
    call cpu_time( start_time )
    
    ! Set initial values
    x = 1.0d0 ! Initialize coordinate matrix 
    y = 1.0d0 ! Initialize coordinate matrix
    affine = eye(d) ! Initilize the affine transformation matrix
    E = 0.d0 ! Initialize energy values
    
    ! Read data from csv files
    call ReadFeatureLocations(nX, d, fileName1, x) ! Read feature locations ILI 1
    call ReadFeatureLocations(nY, d, fileName2, y) ! Read feature locations ILI 2
    
    ! Set temperatures for analysis
    call DetermineDistanceFeaturePairs(nX, nX, d, x, x, distFeatPairs) ! Distance of all feature pairs
    T_init = ( maxval( distFeatPairs ) )**2 ! Set the initial temperature to start the annealing process
    T_final = 100.d0
    do i = 1, nX
        do j = 1, nY
            if (distFeatPairs(i,j)**2 > 0.d0) then
                if (T_final > (distFeatPairs(i,j)**2)) then
                    T_final = distFeatPairs(i,j)**2
                end if
            end if
        end do
    end do
    T_final = T_final
    
    T = T_init
    T_outlier = T_init
    counter = 0
    
    do while (T > T_final) ! Annealing process
        write(*,'(F10.6)') T ! Write current temperature on screen
        lambda2 = multiplier_lambda2 * T ! Set lambda2 for analysis 
        x_trans = matmul(x, affine) ! Transformed x-coordinates
        call DetermineCorrespondenceMatrix(nX, nY, d, x_trans, y, T, T_outlier, zeta, p, stopAnalysis) ! Determine correspondence matrix
        if (stopAnalysis == .false.) then
            counter = counter + 1 ! Count the number of steps in the annealing process
            E(counter, 1) = T
            call DetermineOptimalTransformation(nX, nY, d, x, y, p, lambda2, T, affine, E(counter, 2)) ! Determine the optimal affine transformation
            !call DetermineOptimalLinearTransformation(nX, nY, d, x, y, p, lambda2, T, affine, E(counter, 2)) ! Determine the optimal affine transformation
        else
            exit
        end if   
        T_old = T
        p_old = p        
        T = T * annealingRate ! Adjust temperature
    end do
    call cpu_time( end_time )
    elapsed_time = end_time - start_time
    x_trans = matmul(x, affine) ! Transformed x-coordinates
    call ReportResults("Results096.csv", nX, nY, d, counter, elapsed_time, T_init, T_old, annealingRate, zeta, lambda2, p_old, x_trans, affine, E)
    end program PipelineFeatureMatching2

    
    ! Subroutine to read the ILI coordinates of the features
    subroutine ReadFeatureLocations(n, d, fileName, location)
        integer, intent(in) :: n, d
        character (len=*), intent(in) :: fileName
        double precision, dimension(n, d), intent(out) :: location
        
        integer :: i
        open(10, FILE = fileName) ! Open csv file
        do i = 1, n
            read(10,*) location(i,1), location(i,2) ! Read coordinates from csv file
        end do
        close(10) ! Close csv file
    end subroutine ReadFeatureLocations
    
    
    ! Subroutine to determine the distance of all feature pairs
    subroutine DetermineDistanceFeaturePairs(nX, nY, d, x, y, distance)
        integer, intent(in) :: nX, nY, d
        double precision, dimension(nX, d), intent(in) :: x
        double precision, dimension(nY, d), intent(in) :: y
        double precision, dimension(nX, nY), intent(out) :: distance
        
        integer :: i, j
        
        do i = 1, nX
            do j = 1, nY
                distance(i, j) = sqrt( sum( (x(i,:) - y(j,:))**2 ) )
            end do
        end do
    end subroutine DetermineDistanceFeaturePairs
    
    
    ! Subroutine to determine the correspondence matrix
    subroutine DetermineCorrespondenceMatrix(nX, nY, d, x, y, T, T_outlier, zeta, p, stopAnalysis)
        use isnan_int
        integer, intent(in) :: nX, nY, d
        double precision, dimension(nX, d), intent(in) :: x
        double precision, dimension(nY, d), intent(in) :: y
        double precision, intent(in) :: T, T_outlier, zeta
        double precision, dimension(nX+1, nY+1), intent(out) :: p
        logical, intent(out) :: stopAnalysis
        
        integer :: i, j
        double precision, dimension(nX, nY) :: distFeatPairs
        double precision :: squreDist_ij
        double precision, dimension(d) :: outlierCenter
    
        ! Normal assignment
        call DetermineDistanceFeaturePairs(nX, nY, d, x, y, distFeatPairs) ! Distance of all feature pairs
        p(1:nX, 1:nY) = ( zeta - (distFeatPairs**2) )/ T - log( T )
        
        ! Outliers in x
        outlierCenter = sum(y, 1) / nY
        do i = 1, nX
            j = nY + 1
            squareDist_ij = sum( (x(i,1:d) - outlierCenter)**2 )
            p(i, j) = - squareDist_ij / T_outlier - log( T_outlier )
        end do
        
        ! Outliers in y
        outlierCenter = sum(x, 1) / nX
        do j = 1, nY
            i = nX + 1
            squareDist_ij = sum( (outlierCenter - y(j,1:d))**2 )
            p(i, j) = - squareDist_ij / T_outlier - log( T_outlier )
        end do
        
        p = exp( p )
        p(nX+1, nY+1) = 0.d0
        call NormalizeCorrespondenceMatrix(nX, nY, p)
        
        if ( isNan(p) == .true. ) then
            stopAnalysis = .true.
        else
            stopAnalysis = .false.
        end if
        
    end subroutine DetermineCorrespondenceMatrix

    
    ! Subroutine to normalize correspondence matrix
    subroutine NormalizeCorrespondenceMatrix(nX, nY, p)
        integer, intent(in) :: nX, nY
        double precision, dimension(nX+1, nY+1), intent(inout) :: p
        
        integer :: i, j, k
        double precision, parameter :: threshold = 1.0d-4
        double precision, dimension(nX) :: rowsum
        double precision, dimension(nY) :: colsum
        logical :: converged 

        do       
            ! Normalize columns
            colsum = sum( p(:,1:nY), 1 )
            do j = 1, nY
                if (colsum(j) > 0.d0) then
                    p(:,j) = p(:,j) / colsum(j)
                end if
            end do            
            ! Normalize rows
            rowsum = sum( p(1:nX,:), 2 )
            do i = 1, nX
                rowsum(i) = sum( p(i,:) )
                if (rowsum(i) > 0.d0) then
                    p(i,:) = p(i,:) / rowsum(i)
                end if
            end do            
            
            ! Check if normalization was successful
            converged = .true.
            colsum = sum( p(:,1:nY), 1 )
            do j = 1, nY
                if ( abs(colsum(j) - 1.d0) > threshold ) then 
                    converged = .false.
                end if
            end do    
         
            if (converged == .true.) exit             
        end do
    end subroutine NormalizeCorrespondenceMatrix
    
    
    ! Subroutine to determine the optimal transformation
    subroutine DetermineOptimalTransformation(nX, nY, d, x, y, p, lambda2, T, affine, E)
        use linear_operators
        integer, intent(in) :: nX, nY, d
        double precision, dimension(nX, d), intent(in) :: x
        double precision, dimension(nY, d), intent(in) :: y
        double precision, dimension(nX+1, nY+1), intent(in) :: p
        double precision, intent(in) :: lambda2, T
        double precision, dimension(d, d), intent(inout) :: affine ! Affine transformation matrix
        double precision, intent(out) :: E ! Energy value
        
        integer :: i, j
        double precision, dimension(nX, nY) :: x1, x2, y1, y2, p_red
        double precision, dimension(d, d) :: A
        double precision, dimension(d) :: B
        
        double precision, dimension(nX, d) :: x_trans
        double precision, dimension(nX, nY) :: distFeatPairs
        
        p_red = p(1:nX, 1:nY) ! Reduced p-matrix
        do i = 1, nX
            y1(i,:) = y(:,1)
            y2(i,:) = y(:,2)
        end do
        do j = 1, nY
            x1(:,j) = x(:,1)
            x2(:,j) = x(:,2)
        end do
        
        A(1,1) = sum( p_red * x1 * x1 ) + lambda2
        A(1,2) = sum( p_red * x1 * x2 )
        A(1,3) = sum( p_red * x1 )
        A(2,1) = A(1, 2)
        A(2,2) = sum( p_red * x2 * x2 ) + lambda2
        A(2,3) = sum( p_red * x2 )
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)
        A(3,3) = sum( p ) + lambda2
        B(1) = sum( p_red * y1 * x1 ) + lambda2
        B(2) = sum( p_red * y1 * x2 )
        B(3) = sum( p_red * y1 )
        affine(:,1) = A .ix. B
        
        B(1) = sum( p_red * y2 * x1 ) 
        B(2) = sum( p_red * y2 * x2 ) + lambda2
        B(3) = sum( p_red * y2 )
        affine(:,2) = A .ix. B
        
        ! Determine the energy function
        x_trans = matmul(x, affine)
        call DetermineDistanceFeaturePairs(nX, nY, d, x_trans, y, distFeatPairs) ! Distance of all feature pairs
        E = sum( p_red * distFeatPairs * distFeatPairs ) 
      !  E = E + lambda2 * ((affine(1,1)-1.d0)**2 + affine(2,1)**2 + + affine(1,2)**2 + (affine(2,2)-1.d0)**2) + affine(3,1)**2 + affine(3,2)**2
    end subroutine DetermineOptimalTransformation

    
    ! Subroutine to determine the optimal transformation: a21 and a12 are zero
    subroutine DetermineOptimalLinearTransformation(nX, nY, d, x, y, p, lambda2, T, affine, E)
        use linear_operators
        integer, intent(in) :: nX, nY, d
        double precision, dimension(nX, d), intent(in) :: x
        double precision, dimension(nY, d), intent(in) :: y
        double precision, dimension(nX+1, nY+1), intent(in) :: p
        double precision, intent(in) :: lambda2, T
        double precision, dimension(d, d), intent(inout) :: affine ! Affine transformation matrix
        double precision, intent(out) :: E ! Energy value

        integer :: i, j
        double precision, dimension(nX, nY) :: x1, x2, y1, y2, p_red
        double precision, dimension(d-1, d-1) :: A
        double precision, dimension(d-1) :: B
        double precision, dimension(d-1) :: help_matrix
        
        double precision, dimension(nX, d) :: x_trans
        double precision, dimension(nX, nY) :: distFeatPairs
        
        p_red = p(1:nX, 1:nY) ! Reduced p-matrix
        do i = 1, nX
            y1(i,:) = y(:,1)
            y2(i,:) = y(:,2)
        end do
        do j = 1, nY
            x1(:,j) = x(:,1)
            x2(:,j) = x(:,2)
        end do
        
        A(1,1) = sum( p_red * x1 * x1 ) + lambda2
        A(1,2) = sum( p_red * x1 )
        A(2,1) = A(1,2)
        A(2,2) = sum( p ) !+ lambda2
        B(1) = sum( p_red * y1 * x1 ) + lambda2
        B(2) = sum( p_red * y1 )
        help_matrix = A .ix. B
        affine(1,1) = help_matrix(1)
        affine(3,1) = help_matrix(2)
        
        A(1,1) = sum( p_red * x2 * x2 ) + lambda2
        A(1,2) = sum( p_red * x2 )
        A(2,1) = A(1,2)
        A(2,2) = sum( p ) !+ lambda2 
        B(1) = sum( p_red * y2 * x2 ) + lambda2
        B(2) = sum( p_red * y2 )
        help_matrix = A .ix. B
        affine(2,2) = help_matrix(1)
        affine(3,2) = help_matrix(2)
        
        ! Determine the energy function
        x_trans = matmul(x, affine)
        call DetermineDistanceFeaturePairs(nX, nY, d, x_trans, y, distFeatPairs) ! Distance of all feature pairs
        E = sum( p_red * distFeatPairs * distFeatPairs )
        E = E + lambda2 * ((affine(1,1)-1.d0)**2 + affine(2,1)**2 + + affine(1,2)**2 + (affine(2,2)-1.d0)**2) + affine(3,1)**2 + affine(3,2)**2
    end subroutine DetermineOptimalLinearTransformation
    
    ! Subroutine to deliver the results in a csv files
    subroutine ReportResults(fileNameResults, nX, nY, d, counter, time, T_init, T_final, annealingRate, zeta, lambda2, p, x_trans, affine, E)
        character (len=*), intent(in) :: fileNameResults
        integer, intent(in) :: nX, nY, d, counter
        double precision, intent(in) :: time, T_init, T_final, annealingRate, zeta, lambda2
        double precision, dimension(nX+1, nY+1), intent(in) :: p ! Correspondence Matrix
        double precision, dimension(nX, d), intent(in) :: x_trans ! Transformed x-coordinates
        double precision, dimension(d, d), intent(in) :: affine ! Affine transformation matrix
        double precision, dimension(2000, 2), intent(in) :: E ! Energy function
        
        integer :: i
    
        open(30, FILE = fileNameResults)
        write(30, *) 'Number of features ILI 1:'
        write(30, '(I5)') nX
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Number of features ILI 2:'
        write(30, '(I5)') nY
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Number of annealing steps:'
        write(30, '(I5)') counter
        write(30, '(1x, F, 100(",", F))')  
        write(30, *) 'Duration of the analysis in seconds:'
        write(30, '(F10.2)') time
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Initial temperature:'
        write(30, '(F10.4)') T_init
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Final temperature:'
        write(30, '(F10.6)') T_final
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Annealing rate:'
        write(30, '(F10.2)') annealingRate
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Zeta:'
        write(30, '(F10.4)') zeta
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Lambda 2:'
        write(30, '(F10.4)') lambda2
        write(30, '(1x, F, 100(",", F))')
        write(30, *) 'Correspondence matrix:'
        do i = 1, nX+1
            write(30, '(1x, F, 200(",", F))') p(i,:)
        end do
        write(30, '(1x, F, 200(",", F))')
        write(30, *) 'Transformed coordinates for ILI 1:'
        do i = 1, nX
            write(30, '(1x, F, 200(",", F))') x_trans(i,:)
        end do       
        write(30, '(1x, F, 200(",", F))')
        write(30, *) 'Affine transformation matrix:'
        do i = 1, d
            write(30, '(1x, F, 200(",", F))') affine(i,:)
        end do
        write(30, '(1x, F, 200(",", F))')
        write(30, *) 'Temperature, E-value:'        
        do i = 1, counter
            write(30, '(1x, F, 200(",", F))') E(i,:)
        end do
        close(30)
    end subroutine ReportResults

