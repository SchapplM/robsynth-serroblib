% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:45
% EndTime: 2019-01-03 10:25:49
% DurationCPUTime: 4.39s
% Computational Cost: add. (5091->293), mult. (15916->498), div. (0->0), fcn. (18069->18), ass. (0->185)
t1061 = cos(qJ(4));
t1058 = cos(pkin(6));
t1059 = sin(qJ(2));
t1008 = t1058 * t1059;
t950 = cos(qJ(1));
t1000 = t950 * t1008;
t1060 = sin(qJ(1));
t949 = cos(qJ(2));
t927 = t1060 * t949 + t1000;
t937 = sin(pkin(14));
t941 = cos(pkin(14));
t940 = sin(pkin(6));
t1039 = t940 * t950;
t1010 = t1060 * t1059;
t1017 = t949 * t1058;
t926 = -t1017 * t950 + t1010;
t939 = sin(pkin(7));
t943 = cos(pkin(7));
t998 = t1039 * t939 + t926 * t943;
t1064 = t927 * t937 + t941 * t998;
t918 = t1039 * t943 - t926 * t939;
t938 = sin(pkin(8));
t942 = cos(pkin(8));
t1075 = t1064 * t942 + t918 * t938;
t900 = -t927 * t941 + t937 * t998;
t946 = sin(qJ(4));
t862 = t900 * t1061 + t1075 * t946;
t882 = t1064 * t938 - t918 * t942;
t945 = sin(qJ(5));
t948 = cos(qJ(5));
t847 = t862 * t948 - t882 * t945;
t944 = sin(qJ(6));
t1088 = t847 * t944;
t947 = cos(qJ(6));
t1087 = t847 * t947;
t1036 = t950 * t949;
t996 = t1060 * t1008;
t913 = -qJD(1) * t996 - qJD(2) * t1010 + (qJD(2) * t1058 + qJD(1)) * t1036;
t1020 = qJD(1) * t1060;
t1009 = t940 * t1020;
t983 = t1060 * t1017 + t950 * t1059;
t912 = qJD(1) * t983 + qJD(2) * t927;
t987 = t1009 * t939 - t912 * t943;
t1063 = t913 * t937 - t941 * t987;
t1052 = t912 * t939;
t988 = t1009 * t943 + t1052;
t857 = t1063 * t938 + t942 * t988;
t1086 = qJD(5) * t847 + t857 * t948;
t1004 = t862 * t945 + t882 * t948;
t1085 = qJD(5) * t1004 + t857 * t945;
t879 = t913 * t941 + t937 * t987;
t1083 = -qJD(4) * t862 + t879 * t946;
t1077 = t900 * t946;
t1016 = t1058 * t939;
t1023 = t1059 * t937;
t1037 = t943 * t949;
t993 = t1037 * t941 - t1023;
t967 = t1016 * t941 + t940 * t993;
t1040 = t939 * t949;
t1029 = t940 * t1040;
t990 = t1058 * t943 - t1029;
t1070 = t990 * t938 + t967 * t942;
t982 = t996 - t1036;
t992 = t1060 * t939 * t940 - t943 * t983;
t974 = -t937 * t982 - t941 * t992;
t1025 = t943 * t1060;
t1013 = t940 * t1025;
t994 = t939 * t983 + t1013;
t1066 = t994 * t938 - t974 * t942;
t1068 = t1063 * t942 - t938 * t988;
t1074 = -t879 * t1061 + t1068 * t946;
t1062 = r_i_i_C(3) + pkin(13);
t911 = qJD(1) * t1000 + qJD(2) * t983 + t949 * t1020;
t976 = t982 * qJD(2);
t957 = qJD(1) * t998 + t943 * t976;
t954 = t911 * t937 + t941 * t957;
t958 = qJD(1) * t918 - t939 * t976;
t951 = -t954 * t938 + t958 * t942;
t1072 = t938 * t958 + t942 * t954;
t1071 = t1070 * t1061;
t1069 = t1066 * t1061;
t1067 = t1075 * t1061;
t997 = qJD(6) * (t944 * r_i_i_C(1) + t947 * r_i_i_C(2));
t1041 = t939 * t942;
t1038 = t941 * t943;
t959 = qJD(1) * t926 + t976;
t885 = t1038 * t911 - t937 * t959;
t867 = -t911 * t1041 - t885 * t938;
t901 = t937 * t992 - t941 * t982;
t864 = t901 * t1061 + t1066 * t946;
t884 = t938 * t974 + t942 * t994;
t1065 = -t864 * t945 + t884 * t948;
t1057 = qJ(3) * t939;
t887 = -t1038 * t913 + t912 * t937;
t1053 = t887 * t938;
t1035 = qJD(2) * t940;
t919 = t993 * t1035;
t1047 = t919 * t938;
t1022 = t1059 * t941;
t991 = t1022 * t943 + t937 * t949;
t923 = t991 * t940;
t1046 = t923 * t942;
t1043 = t937 * t943;
t1042 = t939 * t938;
t1034 = qJD(3) * t939;
t1033 = qJD(4) * t946;
t1032 = qJD(6) * t944;
t1031 = qJD(6) * t947;
t1028 = t938 * t1061;
t1027 = t939 * t1059;
t1026 = t942 * t1061;
t1019 = qJD(4) * t1061;
t1018 = pkin(11) * t942 + qJ(3);
t1015 = t939 * t1028;
t1014 = t940 * t1027;
t1012 = t942 * t1027;
t1011 = qJD(2) * t1029;
t1006 = t938 * t1014;
t1005 = t938 * t1011;
t849 = t864 * t948 + t884 * t945;
t905 = -t1038 * t927 + t926 * t937;
t906 = -t1043 * t927 - t926 * t941;
t870 = t906 * t1061 + (t1042 * t927 + t905 * t942) * t946;
t891 = t1041 * t927 - t905 * t938;
t850 = t870 * t948 + t891 * t945;
t907 = t1038 * t982 + t937 * t983;
t908 = t1043 * t982 - t941 * t983;
t872 = t908 * t1061 + (-t1042 * t982 + t907 * t942) * t946;
t892 = -t1041 * t982 - t907 * t938;
t851 = t872 * t948 + t892 * t945;
t916 = t940 * t1022 + (t1037 * t940 + t1016) * t937;
t877 = t916 * t1061 + t1070 * t946;
t897 = -t938 * t967 + t942 * t990;
t855 = t877 * t948 + t897 * t945;
t1003 = -t877 * t945 + t897 * t948;
t924 = (-t1023 * t943 + t941 * t949) * t940;
t890 = t924 * t1061 + (t1006 - t1046) * t946;
t909 = t1012 * t940 + t923 * t938;
t873 = t890 * t948 + t909 * t945;
t1002 = qJD(2) * t1014;
t1001 = r_i_i_C(1) * t947 - r_i_i_C(2) * t944 + pkin(5);
t995 = t1061 * t1006;
t978 = qJD(2) * t1046;
t975 = -t1001 * t948 - t1062 * t945 - pkin(4);
t973 = t1026 * t1064 + t1028 * t918 - t1077;
t971 = t1015 * t927 + t1026 * t905 - t906 * t946;
t970 = -t1015 * t982 + t1026 * t907 - t908 * t946;
t964 = -t1026 * t923 - t924 * t946 + t995;
t960 = t948 * t997 + (t1001 * t945 - t1062 * t948) * qJD(5);
t921 = qJD(2) * t924;
t920 = (-t1037 * t937 - t1022) * t1035;
t904 = (t938 * t991 + t1012) * t1035;
t903 = t1011 * t942 + t1047;
t888 = -t1043 * t913 - t912 * t941;
t886 = t1043 * t911 + t941 * t959;
t878 = -t911 * t941 + t937 * t957;
t876 = t916 * t946 - t1071;
t868 = t1041 * t913 - t1053;
t866 = t920 * t1061 + (-t919 * t942 + t1005) * t946 + t964 * qJD(4);
t865 = qJD(4) * t890 - t1005 * t1061 + t1026 * t919 + t920 * t946;
t863 = t901 * t946 - t1069;
t859 = t1067 - t1077;
t853 = -t916 * t1033 + t921 * t1061 + (t1002 * t938 - t978) * t946 + t1071 * qJD(4);
t852 = -qJD(2) * t995 + t916 * t1019 + t1033 * t1070 + t1061 * t978 + t921 * t946;
t844 = t866 * t948 + t903 * t945 + (-t890 * t945 + t909 * t948) * qJD(5);
t842 = t888 * t1061 + (t1042 * t913 + t887 * t942) * t946 + t971 * qJD(4);
t841 = qJD(4) * t870 - t1015 * t913 - t1026 * t887 + t888 * t946;
t840 = t886 * t1061 + (-t1042 * t911 + t885 * t942) * t946 + t970 * qJD(4);
t839 = qJD(4) * t872 + t1015 * t911 - t1026 * t885 + t886 * t946;
t838 = qJD(5) * t1003 + t853 * t948 + t904 * t945;
t836 = t973 * qJD(4) + t1074;
t835 = -t1026 * t1063 + t1028 * t988 - t1083;
t834 = -qJD(4) * t1067 + t1033 * t900 - t1074;
t833 = t1061 * t1068 + t1083;
t832 = qJD(4) * t1069 - t901 * t1033 + t878 * t1061 + t1072 * t946;
t831 = t901 * t1019 + t1033 * t1066 - t1061 * t1072 + t878 * t946;
t830 = t842 * t948 + t868 * t945 + (-t870 * t945 + t891 * t948) * qJD(5);
t828 = t840 * t948 + t867 * t945 + (-t872 * t945 + t892 * t948) * qJD(5);
t826 = t836 * t948 - t1085;
t824 = t834 * t948 + t1085;
t822 = qJD(5) * t1065 + t832 * t948 + t951 * t945;
t821 = qJD(5) * t849 + t832 * t945 - t948 * t951;
t820 = t822 * t947 + t831 * t944 + (-t849 * t944 + t863 * t947) * qJD(6);
t819 = -t822 * t944 + t831 * t947 + (-t849 * t947 - t863 * t944) * qJD(6);
t1 = [(t826 * t947 + t835 * t944) * r_i_i_C(1) + (-t826 * t944 + t835 * t947) * r_i_i_C(2) + t826 * pkin(5) + t836 * pkin(4) + t835 * pkin(12) - t879 * pkin(3) - t913 * pkin(2) - qJ(3) * t1052 + t1062 * (t836 * t945 + t1086) + ((-t947 * t973 - t1088) * r_i_i_C(1) + (t944 * t973 - t1087) * r_i_i_C(2)) * qJD(6) + t918 * qJD(3) - t857 * pkin(11) + (-t950 * pkin(1) + (-pkin(10) * t1060 - qJ(3) * t1025) * t940) * qJD(1) (t828 * t947 + t839 * t944 + (-t851 * t944 - t947 * t970) * qJD(6)) * r_i_i_C(1) + (-t828 * t944 + t839 * t947 + (-t851 * t947 + t944 * t970) * qJD(6)) * r_i_i_C(2) + t828 * pkin(5) + t840 * pkin(4) + t839 * pkin(12) + t886 * pkin(3) + t959 * pkin(2) - t911 * t1057 - t982 * t1034 + t1062 * (qJD(5) * t851 + t840 * t945 - t867 * t948) + t867 * pkin(11), t958 (t1031 * t864 + t832 * t944) * r_i_i_C(1) + (-t1032 * t864 + t832 * t947) * r_i_i_C(2) + t832 * pkin(12) + t975 * t831 + t960 * t863, -t1001 * t821 + t1062 * t822 - t1065 * t997, r_i_i_C(1) * t819 - r_i_i_C(2) * t820; -pkin(1) * t1020 - t911 * pkin(2) + t878 * pkin(3) + t832 * pkin(4) + t822 * pkin(5) + t831 * pkin(12) + t820 * r_i_i_C(1) + t819 * r_i_i_C(2) + qJD(3) * t1013 + t983 * t1034 - t959 * t1057 + t1062 * t821 + (qJ(3) * t943 + pkin(10)) * qJD(1) * t1039 + t951 * pkin(11) (t830 * t947 + t841 * t944) * r_i_i_C(1) + (-t830 * t944 + t841 * t947) * r_i_i_C(2) + t830 * pkin(5) + t842 * pkin(4) + t841 * pkin(12) + t888 * pkin(3) - pkin(11) * t1053 - t912 * pkin(2) + t1062 * (qJD(5) * t850 + t842 * t945 - t868 * t948) + ((-t850 * t944 - t947 * t971) * r_i_i_C(1) + (-t850 * t947 + t944 * t971) * r_i_i_C(2)) * qJD(6) + (t927 * qJD(3) + t1018 * t913) * t939, t988 (-t1031 * t862 + t834 * t944) * r_i_i_C(1) + (t1032 * t862 + t834 * t947) * r_i_i_C(2) + t834 * pkin(12) + t975 * t833 + t960 * t859, t1062 * t824 - t1004 * t997 + t1001 * (-t834 * t945 + t1086) (-t824 * t944 + t833 * t947) * r_i_i_C(1) + (-t824 * t947 - t833 * t944) * r_i_i_C(2) + ((-t859 * t944 + t1087) * r_i_i_C(1) + (-t859 * t947 - t1088) * r_i_i_C(2)) * qJD(6); 0 (t844 * t947 + t865 * t944) * r_i_i_C(1) + (-t844 * t944 + t865 * t947) * r_i_i_C(2) + t844 * pkin(5) + t866 * pkin(4) + t865 * pkin(12) + t920 * pkin(3) + pkin(11) * t1047 + t1062 * (qJD(5) * t873 + t866 * t945 - t903 * t948) + ((-t873 * t944 - t947 * t964) * r_i_i_C(1) + (-t873 * t947 + t944 * t964) * r_i_i_C(2)) * qJD(6) + (qJD(3) * t1027 + (-pkin(2) * t1059 + t1018 * t1040) * qJD(2)) * t940, t1002 (t1031 * t877 + t853 * t944) * r_i_i_C(1) + (-t1032 * t877 + t853 * t947) * r_i_i_C(2) + t853 * pkin(12) + t975 * t852 + t960 * t876, t1062 * t838 - t1003 * t997 + t1001 * (-qJD(5) * t855 - t853 * t945 + t904 * t948) (-t838 * t944 + t852 * t947) * r_i_i_C(1) + (-t838 * t947 - t852 * t944) * r_i_i_C(2) + ((-t855 * t947 - t876 * t944) * r_i_i_C(1) + (t855 * t944 - t876 * t947) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
