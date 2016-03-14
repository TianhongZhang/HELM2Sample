package org.helm2.sample;


import java.util.List;
import org.helm.notation2.MoleculeInfo;
import org.helm.notation2.MonomerFactory;
import org.helm.notation2.tools.HELM1Utils;
import org.helm.notation2.tools.HELM2NotationUtils;
import org.helm.notation2.calculation.MoleculeInformation;
import org.helm.notation2.tools.PeptideUtils;
import org.helm.notation2.tools.RNAUtils;
import org.helm.notation2.tools.SMILES;
import org.helm.notation2.tools.Validation;
import org.helm.notation2.parser.notation.HELM2Notation;
import org.helm.notation2.parser.notation.polymer.PeptideEntity;
import org.helm.notation2.parser.notation.polymer.PolymerNotation;
import org.helm.notation2.parser.notation.polymer.RNAEntity;

/**
 *
 * @author ZHANGTIANHONG
 */
public class HELM2Sample {

    public static void main(String[] args) {

        try {
            // initialize default monomer database
            MonomerFactory.getInstance();

            String notation;

            // cyclic peptide
            notation = "PEPTIDE1{A.A.C.G.K.[dK].C.H.A}$PEPTIDE1,PEPTIDE1,3:R3-7:R3$$$";
            validate(notation);
            validate(notation);
            getMonomerCount(notation);
            getCanonicalNotation(notation);
            getCanonicalSmiles(notation);
            getMoleculeInformation(notation);
            getAminoAcidSequence(notation);

            notation = "PEPTIDE1{A.A.C.G.[dK].E.C.H.A}$PEPTIDE1,PEPTIDE1,3:R3-7:R3$$$V2.0";
            validate(notation);
            getMonomerCount(notation);
            getCanonicalNotation(notation);
            getCanonicalSmiles(notation);
            getMoleculeInformation(notation);
            getAminoAcidSequence(notation);

            // peptide-chem conjugate
            notation = "PEPTIDE1{A.G.G.G.[seC].C.K.K.K.K}|CHEM1{MCC}$PEPTIDE1,CHEM1,9:R3-1:R1$$$";
            validate(notation);
            getMonomerCount(notation);
            getCanonicalNotation(notation);
            getCanonicalSmiles(notation);
            getMoleculeInformation(notation);
            getAminoAcidSequence(notation);

            // oligo-Chem Conjugate with dynamic chem modifier
            notation = "RNA1{R(A)P.[mR](A)P}|CHEM1{[*]OCCOCCOCCO[*] |$_R1;;;;;;;;;;;_R3$|}$RNA1,CHEM1,6:R2-1:R1$$$";
            validate(notation);
            getMonomerCount(notation);
            getCanonicalNotation(notation);
            getCanonicalSmiles(notation);
            getMoleculeInformation(notation);
            getNucleotideSequence(notation);

            // Invalid Notation with unknown Monomer xyz
            notation = "PEPTIDE1{A.A.C.G.[dK].[xyz].E.C.H.A}$$$$";
            validate(notation);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private static void validate(String notation) {
        System.out.println("/*********************");
        try {
            /*
             * All new methods in the HELM2NotationToolkit work with the
             * HELM2Notation object. It contains all necessary notation objects. Now
             * the class HELM2NotationUtils provides the method to parse HELM in
             * version 1 or 2 and generates the HELM2Notation object. In the past
             * HELM2NotationUtils was ContainerHELM2. Thank you very much for your
             * comment due to the InterConnections.
             */
            HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
            Validation.validateNotationObjects(helm2notation);

            System.out.println("simple polymer String: " + helm2notation.getListOfPolymers().toString());

            System.out.println("connection String: " + HELM2NotationUtils.getAllEdgeConnections(helm2notation.getListOfConnections()).toString());
            System.out.println("Base Pair String: " + HELM2NotationUtils.getAllBasePairConnections(helm2notation.getListOfConnections()).toString());

            System.out.println("Annotations section" + helm2notation.getListOfAnnotations().toString());

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

    private static void getMonomerCount(String notation) {
        System.out.println("/*********************");
        try {
            System.out.println("Testing getMonomerCount for: " + notation);
            System.out.println("Monomer Count: " + HELM2NotationUtils.getTotalMonomerCount(HELM2NotationUtils.readNotation(notation)));
            System.out.println("*********************");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

    private static void getCanonicalNotation(String notation) {
        System.out.println("/*********************");
        try {
            System.out.println("Testing getCanonicalNotation for: " + notation);
            System.out.println("Canonical Notation: " + HELM1Utils.getCanonical(HELM2NotationUtils.readNotation(notation)));
            System.out.println("*********************");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

    private static void getCanonicalSmiles(String notation) {

        System.out.println("/*********************");
        try {
            System.out.println("Testing getCanonicalSmiles for: " + notation);
            System.out.println("Canonical SMILES: " + SMILES.getCanonicalSMILESForAll(HELM2NotationUtils.readNotation(notation)));
            System.out.println("*********************");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

    private static void getMoleculeInformation(String notation) {
        System.out.println("/*********************");

        try {
            System.out.println("Testing testGetMoleculeInformation for: " + notation);
            HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);

            System.out.println("MW = " + MoleculeInformation.getMolecularWeight(helm2notation));

            System.out.println("MF = " + MoleculeInformation.getMolecularFormular(helm2notation));

            System.out.println("Mass = " + MoleculeInformation.getExactMass(helm2notation));

            MoleculeInfo mi = MoleculeInformation.getMoleculeProperties(helm2notation);

            System.out.println(mi.getMolecularWeight() + " " + mi.getMolecularFormula() + " " + mi.getExactMass() + " " + mi.getExtinctionCoefficient());

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

    private static void getNucleotideSequence(String notation) {
        /*
         * the class RNAUtils.java provides all methods to parse sequence from a rna
         * molecule and other rna dependent methods
         */
        System.out.println("/*********************");
        try {
            /* First generate HELM2Notation */
            HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
            System.out.println("Testing getNucleotideSequence for: " + notation);
            List<PolymerNotation> polymers = helm2notation.getListOfPolymers();
            for (PolymerNotation polymer : polymers) {
                if (polymer.getPolymerID() instanceof RNAEntity) {
                    String seq = RNAUtils.getNaturalAnalogSequence(polymer);
                    System.out.println("Nucleotide Sequence: " + seq);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");

    }

    private static void getAminoAcidSequence(String notation) {
        /*
         * the class PeptideUtils.java provides a method to parse sequence from a peptide
         * molecule
         */

        System.out.println("/*********************");
        try {
            HELM2Notation helm2notation = HELM2NotationUtils.readNotation(notation);
            System.out.println("Testing getAminoAcidSequence for: " + notation);
            List<PolymerNotation> polymers = helm2notation.getListOfPolymers();
            for (PolymerNotation polymer : polymers) {
                if (polymer.getPolymerID() instanceof PeptideEntity) {
                    String seq = PeptideUtils.getNaturalAnalogueSequence(polymer);
                    System.out.println("AA Sequence: " + seq);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        System.out.println("*********************/");
    }

}
