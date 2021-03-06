#ifndef _dgemm_param_tt_h_
#define _dgemm_param_tt_h_

#ifdef PastixDouble_PRECISION
//index, DIM_X, DIM_Y, BLK_M, BLK_N, BLK_K, dim_vec, DIM_XA, DIM_YA, DIM_XB, DIM_YB
#define TT_V_0 4, 8, 8, 24, 8, 1, 4, 8, 4, 8
#define TT_V_1 4, 8, 8, 32, 8, 1, 4, 8, 4, 8
#define TT_V_2 4, 8, 8, 40, 8, 1, 4, 8, 4, 8
#define TT_V_3 4, 8, 16, 16, 8, 1, 4, 8, 4, 8
#define TT_V_4 4, 8, 16, 24, 8, 1, 4, 8, 4, 8
#define TT_V_5 4, 8, 16, 32, 8, 1, 4, 8, 4, 8
#define TT_V_6 4, 8, 24, 16, 8, 1, 4, 8, 4, 8
#define TT_V_7 4, 8, 24, 24, 8, 1, 4, 8, 4, 8
#define TT_V_8 4, 8, 32, 16, 8, 1, 4, 8, 4, 8
#define TT_V_9 4, 16, 16, 32, 16, 1, 4, 16, 4, 16
#define TT_V_10 8, 4, 16, 16, 8, 1, 8, 4, 8, 4
#define TT_V_11 8, 4, 16, 24, 8, 1, 8, 4, 8, 4
#define TT_V_12 8, 4, 16, 32, 8, 1, 8, 4, 8, 4
#define TT_V_13 8, 4, 24, 8, 8, 1, 8, 4, 8, 4
#define TT_V_14 8, 4, 24, 16, 8, 1, 8, 4, 8, 4
#define TT_V_15 8, 4, 24, 24, 8, 1, 8, 4, 8, 4
#define TT_V_16 8, 4, 32, 8, 8, 1, 8, 4, 8, 4
#define TT_V_17 8, 4, 32, 16, 8, 1, 8, 4, 8, 4
#define TT_V_18 8, 4, 40, 8, 8, 1, 8, 4, 8, 4
#define TT_V_19 8, 8, 16, 24, 8, 1, 8, 8, 8, 8
#define TT_V_20 8, 8, 16, 24, 16, 1, 8, 8, 8, 8
#define TT_V_21 8, 8, 16, 32, 8, 1, 8, 8, 8, 8
#define TT_V_22 8, 8, 16, 32, 16, 1, 8, 8, 8, 8
#define TT_V_23 8, 8, 16, 40, 8, 1, 8, 8, 8, 8
#define TT_V_24 8, 8, 16, 48, 8, 1, 8, 8, 8, 8
#define TT_V_25 8, 8, 16, 56, 8, 1, 8, 8, 8, 8
#define TT_V_26 8, 8, 16, 64, 8, 1, 8, 8, 8, 8
#define TT_V_27 8, 8, 24, 16, 8, 1, 8, 8, 8, 8
#define TT_V_28 8, 8, 24, 16, 16, 1, 8, 8, 8, 8
#define TT_V_29 8, 8, 24, 24, 8, 1, 8, 8, 8, 8
#define TT_V_30 8, 8, 24, 24, 16, 1, 8, 8, 8, 8
#define TT_V_31 8, 8, 24, 32, 8, 1, 8, 8, 8, 8
#define TT_V_32 8, 8, 24, 40, 8, 1, 8, 8, 8, 8
#define TT_V_33 8, 8, 24, 48, 8, 1, 8, 8, 8, 8
#define TT_V_34 8, 8, 24, 56, 8, 1, 8, 8, 8, 8
#define TT_V_35 8, 8, 24, 64, 8, 1, 8, 8, 8, 8
#define TT_V_36 8, 8, 32, 16, 8, 1, 8, 8, 8, 8
#define TT_V_37 8, 8, 32, 16, 16, 1, 8, 8, 8, 8
#define TT_V_38 8, 8, 32, 24, 8, 1, 8, 8, 8, 8
#define TT_V_39 8, 8, 32, 32, 8, 1, 8, 8, 8, 8
#define TT_V_40 8, 8, 32, 40, 8, 1, 8, 8, 8, 8
#define TT_V_41 8, 8, 32, 48, 8, 1, 8, 8, 8, 8
#define TT_V_42 8, 8, 32, 56, 8, 1, 8, 8, 8, 8
#define TT_V_43 8, 8, 40, 16, 8, 1, 8, 8, 8, 8
#define TT_V_44 8, 8, 40, 24, 8, 1, 8, 8, 8, 8
#define TT_V_45 8, 8, 40, 32, 8, 1, 8, 8, 8, 8
#define TT_V_46 8, 8, 40, 40, 8, 1, 8, 8, 8, 8
#define TT_V_47 8, 8, 40, 48, 8, 1, 8, 8, 8, 8
#define TT_V_48 8, 8, 48, 16, 8, 1, 8, 8, 8, 8
#define TT_V_49 8, 8, 48, 24, 8, 1, 8, 8, 8, 8
#define TT_V_50 8, 8, 48, 32, 8, 1, 8, 8, 8, 8
#define TT_V_51 8, 8, 48, 40, 8, 1, 8, 8, 8, 8
#define TT_V_52 8, 8, 56, 16, 8, 1, 8, 8, 8, 8
#define TT_V_53 8, 8, 56, 24, 8, 1, 8, 8, 8, 8
#define TT_V_54 8, 8, 56, 32, 8, 1, 8, 8, 8, 8
#define TT_V_55 8, 8, 64, 16, 8, 1, 8, 8, 8, 8
#define TT_V_56 8, 8, 64, 24, 8, 1, 8, 8, 8, 8
#define TT_V_57 8, 16, 16, 48, 16, 1, 8, 16, 8, 16
#define TT_V_58 8, 16, 16, 64, 16, 1, 8, 16, 8, 16
#define TT_V_59 8, 16, 32, 32, 16, 1, 8, 16, 8, 16
#define TT_V_60 8, 16, 32, 48, 16, 1, 8, 16, 8, 16
#define TT_V_61 8, 16, 32, 64, 16, 1, 8, 16, 8, 16
#define TT_V_62 8, 16, 48, 32, 16, 1, 8, 16, 8, 16
#define TT_V_63 8, 16, 48, 48, 16, 1, 8, 16, 8, 16
#define TT_V_64 8, 16, 64, 32, 16, 1, 8, 16, 8, 16
#define TT_V_65 8, 24, 24, 48, 24, 1, 8, 24, 8, 24
#define TT_V_66 8, 32, 32, 64, 32, 1, 8, 32, 8, 32
#define TT_V_67 12, 24, 48, 48, 24, 1, 12, 24, 12, 24
#define TT_V_68 16, 4, 32, 16, 16, 1, 16, 4, 16, 4
#define TT_V_69 16, 8, 32, 32, 16, 1, 16, 8, 16, 8
#define TT_V_70 16, 8, 32, 48, 16, 1, 16, 8, 16, 8
#define TT_V_71 16, 8, 32, 64, 16, 1, 16, 8, 16, 8
#define TT_V_72 16, 8, 48, 16, 16, 1, 16, 8, 16, 8
#define TT_V_73 16, 8, 48, 32, 16, 1, 16, 8, 16, 8
#define TT_V_74 16, 8, 48, 48, 16, 1, 16, 8, 16, 8
#define TT_V_75 16, 8, 64, 16, 16, 1, 16, 8, 16, 8
#define TT_V_76 16, 8, 64, 32, 16, 1, 16, 8, 16, 8
#define TT_V_77 16, 16, 32, 48, 16, 1, 16, 16, 16, 16
#define TT_V_78 16, 16, 32, 48, 32, 1, 16, 16, 16, 16
#define TT_V_79 16, 16, 32, 64, 16, 1, 16, 16, 16, 16
#define TT_V_80 16, 16, 32, 64, 32, 1, 16, 16, 16, 16
#define TT_V_81 16, 16, 48, 32, 16, 1, 16, 16, 16, 16
#define TT_V_82 16, 16, 48, 32, 32, 1, 16, 16, 16, 16
#define TT_V_83 16, 16, 48, 48, 16, 1, 16, 16, 16, 16
#define TT_V_84 16, 16, 48, 48, 32, 1, 16, 16, 16, 16
#define TT_V_85 16, 16, 48, 64, 16, 1, 16, 16, 16, 16
#define TT_V_86 16, 16, 64, 32, 16, 1, 16, 16, 16, 16
#define TT_V_87 16, 16, 64, 32, 32, 1, 16, 16, 16, 16
#define TT_V_88 16, 16, 64, 48, 16, 1, 16, 16, 16, 16
#define TT_V_89 16, 16, 64, 64, 16, 1, 16, 16, 16, 16
#define TT_V_90 16, 32, 64, 64, 32, 1, 16, 32, 16, 32
#define TT_V_91 24, 8, 48, 24, 24, 1, 24, 8, 24, 8
#define TT_V_92 24, 12, 48, 48, 24, 1, 24, 12, 24, 12
#define TT_V_93 32, 8, 64, 32, 32, 1, 32, 8, 32, 8
#define TT_V_94 32, 16, 64, 64, 32, 1, 32, 16, 32, 16
#define TT_V_95 4, 8, 8, 24, 8, 2, 4, 8, 4, 8
#define TT_V_96 4, 8, 8, 32, 8, 2, 4, 8, 4, 8
#define TT_V_97 4, 8, 8, 40, 8, 2, 4, 8, 4, 8
#define TT_V_98 4, 8, 16, 16, 8, 2, 4, 8, 4, 8
#define TT_V_99 4, 8, 16, 24, 8, 2, 4, 8, 4, 8
#define TT_V_100 4, 8, 16, 32, 8, 2, 4, 8, 4, 8
#define TT_V_101 4, 8, 24, 16, 8, 2, 4, 8, 4, 8
#define TT_V_102 4, 8, 24, 24, 8, 2, 4, 8, 4, 8
#define TT_V_103 4, 8, 32, 16, 8, 2, 4, 8, 4, 8
#define TT_V_104 4, 16, 16, 32, 16, 2, 4, 16, 4, 16
#define TT_V_105 8, 8, 16, 32, 16, 2, 8, 8, 8, 8
#define TT_V_106 8, 8, 32, 16, 16, 2, 8, 8, 8, 8
#define TT_V_107 8, 16, 16, 48, 16, 2, 8, 16, 8, 16
#define TT_V_108 8, 16, 16, 64, 16, 2, 8, 16, 8, 16
#define TT_V_109 8, 16, 32, 32, 16, 2, 8, 16, 8, 16
#define TT_V_110 8, 16, 32, 48, 16, 2, 8, 16, 8, 16
#define TT_V_111 8, 16, 32, 64, 16, 2, 8, 16, 8, 16
#define TT_V_112 8, 16, 48, 32, 16, 2, 8, 16, 8, 16
#define TT_V_113 8, 16, 48, 48, 16, 2, 8, 16, 8, 16
#define TT_V_114 8, 16, 64, 32, 16, 2, 8, 16, 8, 16
#define TT_V_115 8, 32, 32, 64, 32, 2, 8, 32, 8, 32
#define TT_V_116 12, 24, 48, 48, 24, 2, 12, 24, 12, 24
#define TT_V_117 16, 16, 32, 64, 32, 2, 16, 16, 16, 16
#define TT_V_118 16, 16, 64, 32, 32, 2, 16, 16, 16, 16
#define TT_V_119 16, 32, 64, 64, 32, 2, 16, 32, 16, 32

#endif

#endif /* _dgemm_param_tt_h_ */
