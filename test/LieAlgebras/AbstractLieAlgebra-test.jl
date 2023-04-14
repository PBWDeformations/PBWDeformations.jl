@testset ExtendedTestSet "All AbstractLieAlgebra.jl tests" begin

    sc = zeros(QQ, 3, 3, 3)
    sc[1, 2, 3] = 1
    sc[2, 1, 3] = -1
    sc[3, 1, 1] = 2
    sc[1, 3, 1] = -2
    sc[3, 2, 2] = -2
    sc[2, 3, 2] = 2
    sl2 = lie_algebra(QQ, sc, ["e", "f", "h"])

    @testset "conformance tests for $desc" for (desc, L, parentT, elemT) in [
        ("sl_2(QQ)", sl2, AbstractLieAlgebra{fmpq}, AbstractLieAlgebraElem{fmpq}),
        ("A_4(QQ)", lie_algebra(QQ, ('A', 4)), AbstractLieAlgebra{fmpq}, AbstractLieAlgebraElem{fmpq}),
        ("B_3(QQ)", lie_algebra(QQ, ('B', 3)), AbstractLieAlgebra{fmpq}, AbstractLieAlgebraElem{fmpq}),
        # ("A_4(CF(4))", lie_algebra(cyclotomic_field(4)[1], ('A', 4)), AbstractLieAlgebra{nf_elem}, AbstractLieAlgebraElem{nf_elem}),
        # ("B_3(CF(4))", lie_algebra(cyclotomic_field(4)[1], ('B', 3)), AbstractLieAlgebra{nf_elem}, AbstractLieAlgebraElem{nf_elem}),
    ]
        lie_algebra_conformance_test(L, parentT, elemT)
    end



end
