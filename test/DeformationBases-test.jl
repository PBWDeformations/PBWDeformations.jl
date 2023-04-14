@testset ExtendedTestSet "All DeformationBases/*.jl tests" begin
    @testset "All DeformationBases/ArcDiagDeformBasis.jl tests" begin
        @testset "arcdiag_to_deformationmap__so(:exterior)" begin
            L = special_orthogonal_lie_algebra(QQ, 4)
            V = exterior_power(standard_module(L), 2)
            sp = smash_product(L, V)

            @testset "not all specialisations are zero" begin
                diag = ArcDiagram(4, 2, [5, 3, 2, 6, 1, 4])
                dm = PD.arcdiag_to_deformationmap__so(diag, sp)
                @test !iszero(dm)
            end

            @testset "deformation is equivariant, d = $deg" for deg in 1:3
                for diag in PD.pbw_arc_diagrams__so(V, deg)
                    dm = PD.arcdiag_to_deformationmap__so(diag, sp)
                    d = deform(sp, dm)
                    @test all(iszero, pbwdeform_eqs(d, disabled=[:b, :c, :d]))
                end
            end
        end

        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4)
                V = exterior_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = ArcDiagDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[0 x_2_3 x_2_4 -x_1_3 -x_1_4 0; -x_2_3 0 x_3_4 x_1_2 0 -x_1_4; -x_2_4 -x_3_4 0 0 x_1_2 x_1_3; x_1_3 -x_1_2 0 0 x_3_4 -x_2_4; x_1_4 0 -x_1_2 -x_3_4 0 x_2_3; 0 x_1_4 -x_1_3 x_2_4 -x_2_3 0]"
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = ArcDiagDeformBasis{QQFieldElem}(sp, 0:3)
                @test length(collect(b)) == 4
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 2, 2], [2, 3, 4], QQFieldElem[1, 3//2, -1//2], 4, 4)
            end

            @testset "SO_5, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5)
                V = exterior_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = ArcDiagDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[0 x_2_3 x_2_4 x_2_5 -x_1_3 -x_1_4 -x_1_5 0 0 0; -x_2_3 0 x_3_4 x_3_5 x_1_2 0 0 -x_1_4 -x_1_5 0; -x_2_4 -x_3_4 0 x_4_5 0 x_1_2 0 x_1_3 0 -x_1_5; -x_2_5 -x_3_5 -x_4_5 0 0 0 x_1_2 0 x_1_3 x_1_4; x_1_3 -x_1_2 0 0 0 x_3_4 x_3_5 -x_2_4 -x_2_5 0; x_1_4 0 -x_1_2 0 -x_3_4 0 x_4_5 x_2_3 0 -x_2_5; x_1_5 0 0 -x_1_2 -x_3_5 -x_4_5 0 0 x_2_3 x_2_4; 0 x_1_4 -x_1_3 0 x_2_4 -x_2_3 0 0 x_4_5 -x_3_5; 0 x_1_5 0 -x_1_3 x_2_5 0 -x_2_3 -x_4_5 0 x_3_4; 0 0 x_1_5 -x_1_4 0 x_2_5 -x_2_4 x_3_5 -x_3_4 0]"
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)
            end

            @testset "SO_4, S²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4)
                V = symmetric_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = ArcDiagDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 2
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[2], Int64[2], QQFieldElem[1], 2, 2)
            end

            @testset "SO_5, S²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5)
                V = symmetric_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = ArcDiagDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 2
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[2], Int64[2], QQFieldElem[1], 2, 2)
            end
        end

    end

    @testset "All DeformationBases/PseudographDeformBasis.jl tests" begin
        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4)
                V = exterior_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, PseudographDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, PseudographDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = PseudographDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[0 x_2_3 x_2_4 -x_1_3 -x_1_4 0; -x_2_3 0 x_3_4 x_1_2 0 -x_1_4; -x_2_4 -x_3_4 0 0 x_1_2 x_1_3; x_1_3 -x_1_2 0 0 x_3_4 -x_2_4; x_1_4 0 -x_1_2 -x_3_4 0 x_2_3; 0 x_1_4 -x_1_3 x_2_4 -x_2_3 0]"
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = PseudographDeformBasis{QQFieldElem}(sp, 0:3)
                @test length(collect(b)) == 4
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 2, 2], [2, 3, 4], QQFieldElem[1, 3//2, -1//2], 4, 4)
            end

            @testset "SO_5, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5)
                V = exterior_power(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, PseudographDeformBasis{QQFieldElem}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 0, 0)
                @test all_pbwdeformations(sp, PseudographDeformBasis{QQFieldElem}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = PseudographDeformBasis{QQFieldElem}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}[0 x_2_3 x_2_4 x_2_5 -x_1_3 -x_1_4 -x_1_5 0 0 0; -x_2_3 0 x_3_4 x_3_5 x_1_2 0 0 -x_1_4 -x_1_5 0; -x_2_4 -x_3_4 0 x_4_5 0 x_1_2 0 x_1_3 0 -x_1_5; -x_2_5 -x_3_5 -x_4_5 0 0 0 x_1_2 0 x_1_3 x_1_4; x_1_3 -x_1_2 0 0 0 x_3_4 x_3_5 -x_2_4 -x_2_5 0; x_1_4 0 -x_1_2 0 -x_3_4 0 x_4_5 x_2_3 0 -x_2_5; x_1_5 0 0 -x_1_2 -x_3_5 -x_4_5 0 0 x_2_3 x_2_4; 0 x_1_4 -x_1_3 0 x_2_4 -x_2_3 0 0 x_4_5 -x_3_5; 0 x_1_5 0 -x_1_3 x_2_5 0 -x_2_3 -x_4_5 0 x_3_4; 0 0 x_1_5 -x_1_4 0 x_2_5 -x_2_4 x_3_5 -x_3_4 0]"
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], QQFieldElem[], 1, 1)

                b = PseudographDeformBasis{QQFieldElem}(sp, 0:3)
                @test length(collect(b)) == 4
                @test all_pbwdeformations(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 3, 4], [2, 3, 4], QQFieldElem[1, 1, 1], 4, 4)
            end
        end

    end

end
