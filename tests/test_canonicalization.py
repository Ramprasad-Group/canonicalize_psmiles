from canonicalize_psmiles.canonicalize import canonicalize



def test_canonicalize_one():
    """Test canonicalization of one SMILES."""
    sm_init = [
        "N(=[*])c1ccc(NC=[*])cc1",
        "[*]=CNc1ccc(cc1)N=[*]",
        "[*]=Nc1ccc(NC=[*])cc1",
        "[*]=Nc1ccc(cc1)NC=[*]",
        "c1(ccc(NC=[*])cc1)N=[*]",
        "c1c(NC=[*])ccc(c1)N=[*]",
        "c1c(ccc(NC=[*])c1)N=[*]",
        "c1cc(N=[*])ccc1NC=[*]",
        "c1cc(N=[*])ccc1NC=[*]",
        "c1cc(ccc1N=[*])NC=[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    # assert one_smiles == ["[*]=CNc1ccc(N=[*])cc1"]
    assert one_smiles == ["[*]NC=Nc1ccc([*])cc1"]


def test_canonicalize_two():
    """Test canonicalization of one SMILES."""
    sm_init = [
        "[*]CCO[*]",
        "[*]OCC[*]",
        "C([*])CO[*]",
        "C(C[*])O[*]",
        "C(CO[*])[*]",
        "C(O[*])C[*]",
        "O([*])CC[*]",
        "O(CC[*])[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]COC[*]"]


def test_canonicalize_three():
    """Test canonicalization of one SMILES."""
    sm_init = [
        "*CCOCCO*",
        "*OCCOCC*",
        "C(*)COCCO*",
        "C(C*)OCCO*",
        "C(CO*)OCC*",
        "C(COCC*)O*",
        "C(COCCO*)*",
        "C(O*)COCC*",
        "C(OCCO*)C*",
        "O(*)CCOCC*",
        "O(CC*)CCO*",
        "O(CCO*)CC*",
        "O(CCOCC*)*",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))
    assert one_smiles == ["[*]COC[*]"]


def test_canonicalize_four():
    """Test canonicalization of one SMILES."""
    sm_init = [
        "C(=[*])(Nc1c(cc(cc1)-c1ccc(NC(=Nc2ccc(N=[*])cc2)C)c(C(=O)O)c1)C(=O)O)C",
        "C(Nc1c(cc(cc1)-c1cc(C(O)=O)c(NC(=[*])C)cc1)C(O)=O)(=Nc1ccc(cc1)N=[*])C",
        "CC(Nc1c(cc(cc1)-c1cc(C(O)=O)c(NC(C)=Nc2ccc(cc2)N=[*])cc1)C(=O)O)=[*]",
        "N(=[*])c1ccc(cc1)N=C(C)Nc1ccc(cc1C(O)=O)-c1ccc(c(C(=O)O)c1)NC(C)=[*]",
        "N(c1c(C(=O)O)cc(cc1)-c1ccc(c(C(=O)O)c1)NC(=[*])C)C(=Nc1ccc(N=[*])cc1)C",
        "c1(c(ccc(-c2cc(C(O)=O)c(NC(=[*])C)cc2)c1)NC(C)=Nc1ccc(N=[*])cc1)C(O)=O",
        "c1c(-c2cc(c(cc2)NC(=Nc2ccc(cc2)N=[*])C)C(O)=O)ccc(NC(=[*])C)c1C(=O)O",
        "c1cc(-c2ccc(NC(=[*])C)c(c2)C(O)=O)cc(C(=O)O)c1NC(C)=Nc1ccc(N=[*])cc1",
        "c1cc(-c2ccc(NC(C)=[*])c(c2)C(=O)O)cc(c1NC(=Nc1ccc(cc1)N=[*])C)C(=O)O",
        "c1cc(N=[*])ccc1N=C(Nc1c(cc(-c2ccc(c(C(=O)O)c2)NC(=[*])C)cc1)C(O)=O)C",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == [
        "[*]=Nc1ccc(N=C(C)Nc2ccc(-c3ccc(NC(=[*])C)c(C(=O)O)c3)cc2C(=O)O)cc1"
    ]


def test_canonicalize_five():
    """Test canonicalization of one SMILES."""

    # Last one is different.
    sm_init = [
        "[*]=COCCC=[*]",
        "[*]=CCCOC=[*]",
        "[*]=CCOCC=[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    two_smiles = list(set(canonicalize_smiles))
    two_smiles = sorted(two_smiles)

    assert two_smiles == ["[*]=CCCOC=[*]", "[*]=CCOCC=[*]"]


def test_canonicalize_six():
    """Test canonicalization of one SMILES."""
    sm_init = [
        "[*]C=CC(C)(C)C([*])(C)C",
        "[*]C=CC(C)(C)C([*])(C)C",
        "[*]C=CC(C)(C)C([*])(C)C",
    ]
    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]C=CC(C)(C)C([*])(C)C"]


def test_canonicalize_seven():
    """Test canonicalization of one SMILES."""

    # We cannot form a ring in this case because the
    # star neighbors are already connected. (polystyrene)
    # Just running RDKits canonicalization does the job though.

    sm_init = [
        "CC(C)(/C(=C/[*])[*])C",
        "C(=C\\[*])(\\[*])C(C)(C)C",
        "C(C)(C)(C)/C([*])=C/[*]",
        "C(=C(\\[*])C(C)(C)C)\\[*]",
        "[*]/C=C(/C(C)(C)C)[*]",
        "C(C)(C)(C)/C([*])=C/[*]",
        "C(=C\\[*])(\\[*])C(C)(C)C",
        "C(C)(C)(C)/C([*])=C/[*]",
        "C(C)(C)(/C([*])=C/[*])C",
        "CC(C)(C)/C([*])=C/[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == [r"[*]/C=C(\[*])C(C)(C)C"]


def test_canonicalize_polystyrene():
    """Test canonicalization of one SMILES."""

    sm_init = [
        "c1ccc(cc1)C(C[*])[*]",
        "c1ccc(cc1)C(C[*])[*]",
        "c1cc(C([*])C[*])ccc1",
        "c1ccccc1C(C[*])[*]",
        "[*]C(C[*])c1ccccc1",
        "c1ccc(C(C[*])[*])cc1",
        "C(C([*])c1ccccc1)[*]",
        "c1cccc(c1)C(C[*])[*]",
        "c1cc(C([*])C[*])ccc1",
        "C([*])C(c1ccccc1)[*]",
        "[*]CC(c1ccccc1)[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]CC([*])c1ccccc1"]


def test_canonicalize_cis_trans():
    """Test canonicalization of cis/trans SMILES."""

    sm_init = [
        "CC(/C([*])=C\\[*])(C)C",
        "[*]/C(C(C)(C)C)=C/[*]",
        "C(\\[*])=C(\\C(C)(C)C)[*]",
        "CC(/C([*])=C\\[*])(C)C",
        "C(C)(/C([*])=C\\[*])(C)C",
        "C(\\[*])(C(C)(C)C)=C/[*]",
        "CC(C)(C)/C(=C\\[*])[*]",
        "C(\\[*])=C(\\C(C)(C)C)[*]",
        "CC(C)(/C([*])=C\\[*])C",
        "CC(C)(C)/C(=C\\[*])[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]/C=C(/[*])C(C)(C)C"]


def test_canonicalize_polyethylene():
    """Test canonicalization of cis/trans SMILES."""

    sm_init = ["[*]CCCC[*]", "[*]CCC[*]", "[*]CC[*]", "[*]C[*]"]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]C[*]"]


def test_canonicalize_tricky():
    """Test canonicalization of tricky SMILES string."""

    sm_init = [
        "C(CCc1c([*])sc([*])c1)CCCO[Si](C)(C)C",
        "C[Si](OCCCCCCc1c(sc(c1)[*])[*])(C)C",
        "c1([*])sc([*])cc1CCCCCCO[Si](C)(C)C",
        "c1(c(cc([*])s1)CCCCCCO[Si](C)(C)C)[*]",
        "C[Si](OCCCCCCc1c(sc(c1)[*])[*])(C)C",
        "C(CCCCc1c([*])sc(c1)[*])CO[Si](C)(C)C",
        "C(CO[Si](C)(C)C)CCCCc1c([*])sc([*])c1",
        "C[Si](C)(C)OCCCCCCc1cc([*])sc1[*]",
        "c1(CCCCCCO[Si](C)(C)C)c([*])sc(c1)[*]",
        "C[Si](C)(OCCCCCCc1cc([*])sc1[*])C",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]c1cc(CCCCCCO[Si](C)(C)C)c([*])s1"]


def test_canonicalize_connected_neighbors():
    """Test canonicalization of SMILES strings that already have connected neighbors."""

    sm_init = [r"[*]/C=C(\[*])C(C)(C)C", r"[*]\C=C(/[*])C(C)(C)C"]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == [r"[*]/C=C(\[*])C(C)(C)C"]


def test_canonicalize_same_neighbors():
    """Test canonicalization of SMILES strings
    that star neighbors is are the same atom"""

    sm_init = [
        "[*]C(F)([*])F",
        "[*]C([*])(F)F",
        "C(F)([*])([*])F",
        "FC(F)([*])[*]",
        "C([*])([*])(F)F",
        "[*]C(F)([*])F",
        "C([*])(F)([*])F",
        "[*]C(F)([*])F",
        "[*]C([*])(F)F",
        "[*]C(F)([*])F",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]C([*])(F)F"]


def test_canonicalize_tricky_2():
    """That's a tricky one. All atoms are aromatic."""

    sm_init = [
        "[*]c1ccc(-c2ccc(N3C(=O)c4cc5c(cc4C3=O)C(C(F)(F)F)(C(F)(F)F)c3cc4c(cc3O5)C(=O)N([*])C4=O)cc2C(F)(F)F)c(C(F)(F)F)c1",  # noqa: E501
        "c1(C(F)(F)F)c(-c2c(C(F)(F)F)cc(cc2)N2C(=O)c3c(cc4Oc5c(C(C(F)(F)F)(c4c3)C(F)(F)F)cc3C(N(C(=O)c3c5)[*])=O)C2=O)ccc(c1)[*]",  # noqa: E501
        "C(F)(c1cc(ccc1-c1c(cc(N2C(c3cc4c(Oc5c(cc6C(=O)N(C(=O)c6c5)[*])C4(C(F)(F)F)C(F)(F)F)cc3C2=O)=O)cc1)C(F)(F)F)[*])(F)F",  # noqa: E501
        "O=C1c2c(C(=O)N1[*])cc1c(Oc3c(cc4C(=O)N(c5cc(c(-c6c(cc(cc6)[*])C(F)(F)F)cc5)C(F)(F)F)C(c4c3)=O)C1(C(F)(F)F)C(F)(F)F)c2",  # noqa: E501
        "C1(N(c2cc(c(cc2)-c2c(C(F)(F)F)cc([*])cc2)C(F)(F)F)C(=O)c2c1cc1c(c2)C(c2cc3C(N(C(c3cc2O1)=O)[*])=O)(C(F)(F)F)C(F)(F)F)=O",  # noqa: E501
        "FC(C1(c2cc3C(N(c4ccc(-c5c(C(F)(F)F)cc([*])cc5)c(c4)C(F)(F)F)C(=O)c3cc2Oc2cc3C(=O)N([*])C(c3cc21)=O)=O)C(F)(F)F)(F)F",  # noqa: E501
        "c1cc([*])cc(c1-c1ccc(N2C(=O)c3cc4c(C(c5c(cc6c(c5)C(=O)N([*])C6=O)O4)(C(F)(F)F)C(F)(F)F)cc3C2=O)cc1C(F)(F)F)C(F)(F)F",  # noqa: E501
        "FC(c1c(ccc(c1)N1C(=O)c2cc3Oc4c(cc5C(N(C(c5c4)=O)[*])=O)C(C(F)(F)F)(C(F)(F)F)c3cc2C1=O)-c1c(cc(cc1)[*])C(F)(F)F)(F)F",  # noqa: E501
        "FC(F)(F)c1c(-c2ccc(cc2C(F)(F)F)N2C(=O)c3c(cc4Oc5cc6C(N(C(=O)c6cc5C(C(F)(F)F)(c4c3)C(F)(F)F)[*])=O)C2=O)ccc(c1)[*]",  # noqa: E501
        "FC(c1cc(N2C(=O)c3cc4c(cc3C2=O)Oc2c(C4(C(F)(F)F)C(F)(F)F)cc3C(=O)N([*])C(=O)c3c2)ccc1-c1ccc(cc1C(F)(F)F)[*])(F)F",  # noqa: E501
        "C(C1(c2c(Oc3c1cc1C(=O)N(C(=O)c1c3)c1ccc(-c3ccc([*])cc3C(F)(F)F)c(C(F)(F)F)c1)cc1c(c2)C(=O)N([*])C1=O)C(F)(F)F)(F)(F)F",  # noqa: E501
    ]
    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == [
        "[*]c1ccc(-c2ccc(N3C(=O)c4cc5c(cc4C3=O)C(C(F)(F)F)(C(F)(F)F)c3cc4c(cc3O5)C(=O)N([*])C4=O)cc2C(F)(F)F)c(C(F)(F)F)c1"  # noqa: E501
    ]


def test_canonicalize_all_in_ring():
    """All atoms are part of a ring and there's only one ring"""
    sm_init = ["[*]C1C=CC([*])CC1"]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))

    assert one_smiles == ["[*]C1C=CC([*])CC1"]


def test_canonicalize_break_aromatic():
    sm_init = [
        "[*]C=Cc1sc([*])c(OCCCC)c1OCCCC",
        "[*]c1sc(c(OCCCC)c1OCCCC)C=C[*]",
        "c1(OCCCC)c(c(C=C[*])sc1[*])OCCCC",
        "C(=Cc1c(OCCCC)c(OCCCC)c([*])s1)[*]",
        "c1(OCCCC)c([*])sc(C=C[*])c1OCCCC",
        "s1c(c(c(OCCCC)c1C=C[*])OCCCC)[*]",
        "C([*])=Cc1sc(c(OCCCC)c1OCCCC)[*]",
        "C(Oc1c(OCCCC)c(sc1C=C[*])[*])CCC",
        "C(Oc1c(OCCCC)c(sc1C=C[*])[*])CCC",
        "c1(OCCCC)c(OCCCC)c(C=C[*])sc1[*]",
        "c1(sc(c(c1OCCCC)OCCCC)C=C[*])[*]",
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))
    assert one_smiles == ["[*]C=Cc1sc([*])c(OCCCC)c1OCCCC"]


def test_canonicalize_multiple():
    sm_init = [
        # Multiple paths between N1 and N2: remove length criteria for ring to break
        "[*]c1sc(-c2sc(-c3sc([*])c4nccnc34)c3c2OCCO3)c2c1OCCO2"
    ]

    canonicalize_smiles = []
    for s in sm_init:
        canonicalize_smiles.append(canonicalize(s))
    one_smiles = list(set(canonicalize_smiles))
    assert one_smiles == ["[*]c1sc(-c2sc(-c3sc([*])c4nccnc34)c3c2OCCO3)c2c1OCCO2"]