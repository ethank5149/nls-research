from figures import figures


def main():
    for fig in figures:
        fig.create()

    print("Placement Guide:")
    print("Figure 1 → Section 3 (Summary of Scenarios), near Table 1")
    print("Figure 2 → Section 3, opening (before detailed calculations)")
    print("Figure 3 → Section 3, after Example 5.4 (rotational DOF)")
    print("Figure 4 → Section 5 (Background), after eq. (5.3)")
    print("Figure 5 → Section 5, after Claim 5 / before Lemma 1")
    print("Figure 6 → Section 5, near definition of Q_R")
    print("Figure 7 → Section 1 (Introduction) or start of Section 4")


if __name__ == "__main__":
    main()