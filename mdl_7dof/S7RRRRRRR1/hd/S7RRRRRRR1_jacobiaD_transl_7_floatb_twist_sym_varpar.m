% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_7_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_7_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:21:08
% EndTime: 2018-11-26 21:21:13
% DurationCPUTime: 6.02s
% Computational Cost: add. (2918->345), mult. (8814->601), div. (0->0), fcn. (9592->14), ass. (0->214)
t978 = cos(qJ(2));
t1072 = qJD(2) * t978;
t979 = cos(qJ(1));
t1074 = qJD(1) * t979;
t971 = sin(qJ(2));
t972 = sin(qJ(1));
t1003 = t972 * t1072 + t971 * t1074;
t1066 = qJD(4) * t971;
t1075 = qJD(1) * t978;
t1033 = qJD(3) + t1075;
t1068 = qJD(3) * t978;
t970 = sin(qJ(3));
t1039 = t970 * t1068;
t1073 = qJD(2) * t971;
t1046 = t972 * t1073;
t1076 = qJD(1) * t972;
t977 = cos(qJ(3));
t1077 = t979 * t977;
t908 = t1033 * t1077 - t972 * t1039 - t977 * t1046 - t970 * t1076;
t1078 = t977 * t978;
t1086 = t970 * t979;
t941 = t972 * t1078 + t1086;
t969 = sin(qJ(4));
t976 = cos(qJ(4));
t868 = (t972 * t1066 + t908) * t976 + (-qJD(4) * t941 + t1003) * t969;
t1090 = t969 * t971;
t923 = t972 * t1090 + t941 * t976;
t1082 = t972 * t970;
t940 = t978 * t1082 - t1077;
t968 = sin(qJ(5));
t975 = cos(qJ(5));
t893 = t923 * t975 - t940 * t968;
t1006 = t978 * t1086 + t972 * t977;
t907 = t1006 * qJD(1) + t941 * qJD(3) - t970 * t1046;
t845 = t893 * qJD(5) + t868 * t968 + t907 * t975;
t966 = sin(qJ(7));
t1118 = t845 * t966;
t973 = cos(qJ(7));
t1117 = t845 * t973;
t1085 = t971 * t976;
t922 = -t972 * t1085 + t941 * t969;
t967 = sin(qJ(6));
t974 = cos(qJ(6));
t872 = t893 * t974 + t922 * t967;
t892 = t923 * t968 + t940 * t975;
t1116 = t872 * t966 + t892 * t973;
t1115 = t872 * t973 - t892 * t966;
t848 = t892 * qJD(5) - t868 * t975 + t907 * t968;
t867 = t923 * qJD(4) - t1003 * t976 + t908 * t969;
t1114 = t872 * qJD(6) - t848 * t967 - t867 * t974;
t1111 = t893 * t967 - t922 * t974;
t1110 = t867 * t967;
t1032 = qJD(4) * t977 - qJD(2);
t1014 = t1032 * t969;
t1069 = qJD(3) * t976;
t1108 = (t970 * t1069 + t1014) * t971;
t1031 = qJD(2) * t977 - qJD(4);
t1013 = t1031 * t978;
t1063 = qJD(5) * t970;
t1036 = t971 * t1063;
t1103 = -t976 * t1013 + t1036 + t1108;
t1079 = t976 * t978;
t904 = t1031 * t1079 - t1108;
t1101 = t904 - t1036;
t1071 = qJD(2) * t979;
t1002 = -t978 * t1071 + t971 * t1076;
t1070 = qJD(3) * t971;
t1041 = t970 * t1070;
t1067 = qJD(4) * t969;
t1100 = -t969 * t1041 - t978 * t1067 - t976 * t1073;
t1038 = t977 * t1070;
t1084 = t971 * t977;
t939 = t976 * t1084 - t969 * t978;
t935 = t939 * t979;
t1099 = -qJD(5) * t935 + t1002 * t970 - t979 * t1038;
t933 = t939 * t972;
t1098 = qJD(5) * t933 + t1003 * t970 + t972 * t1038;
t1097 = r_i_i_C(3) + pkin(4);
t1092 = t967 * t969;
t1091 = t968 * t976;
t1089 = t970 * t971;
t1088 = t970 * t975;
t1087 = t970 * t978;
t1083 = t971 * t979;
t1081 = t974 * t975;
t1080 = t975 * t976;
t1065 = qJD(4) * t976;
t1064 = qJD(5) * t968;
t1062 = qJD(5) * t975;
t1061 = qJD(5) * t976;
t1060 = qJD(6) * t967;
t1059 = qJD(6) * t969;
t1058 = qJD(6) * t975;
t1057 = qJD(7) * t966;
t1056 = qJD(7) * t968;
t1055 = qJD(7) * t973;
t1054 = t968 * t1089;
t1053 = t969 * t1089;
t1052 = t971 * t1088;
t1051 = t970 * t1083;
t1050 = t976 * t1078;
t1047 = t978 * t1074;
t1044 = t971 * t1071;
t1037 = t971 * t1065;
t1034 = qJD(1) + t1068;
t1030 = qJD(3) + t1061;
t1028 = r_i_i_C(1) * t973 - r_i_i_C(2) * t966;
t1025 = (-t1031 * t971 - t1039) * t976 + (-t1063 - t1014) * t978;
t988 = t1033 * t972 + t1044;
t906 = t1034 * t1086 + t988 * t977;
t1024 = -t1006 * t1061 - t906;
t1023 = t940 * t1061 - t908;
t905 = -t1034 * t1077 + t988 * t970;
t945 = t978 * t1077 - t1082;
t998 = -qJD(5) * t945 + t1006 * t1067 + t905 * t976;
t1022 = -t1006 * t1059 - t1024 * t968 + t998 * t975;
t997 = -qJD(5) * t941 + t940 * t1067 - t907 * t976;
t1021 = -t1023 * t968 + t940 * t1059 - t997 * t975;
t866 = (t979 * t1066 - t906) * t976 + (-qJD(4) * t945 - t1002) * t969;
t928 = -t976 * t1083 + t945 * t969;
t1020 = t928 * t1058 + t866;
t1019 = t922 * t1058 + t868;
t938 = t969 * t1084 + t1079;
t1018 = t938 * t1058 + t904;
t929 = t969 * t1083 + t945 * t976;
t865 = t929 * qJD(4) + t1002 * t976 - t906 * t969;
t994 = qJD(6) * t929 + t928 * t1064 - t865 * t975;
t1017 = t1020 * t967 + t928 * t1056 + t994 * t974;
t993 = qJD(6) * t923 + t922 * t1064 - t867 * t975;
t1016 = t1019 * t967 + t922 * t1056 + t993 * t974;
t903 = (t969 * t1072 + t1037) * t977 + t1100;
t992 = qJD(6) * t939 + t938 * t1064 - t903 * t975;
t1015 = t1018 * t967 + t938 * t1056 + t992 * t974;
t898 = -t1006 * t968 + t929 * t975;
t875 = t898 * t974 + t928 * t967;
t916 = t972 * t1054 - t933 * t975;
t932 = t938 * t972;
t886 = t916 * t974 - t932 * t967;
t918 = t968 * t1051 - t935 * t975;
t934 = t938 * t979;
t887 = t918 * t974 - t934 * t967;
t921 = t939 * t975 - t1054;
t891 = t921 * t974 + t938 * t967;
t943 = t1050 + t1090;
t927 = -t968 * t1087 + t943 * t975;
t942 = t969 * t1078 - t1085;
t896 = t927 * t974 + t942 * t967;
t897 = t1006 * t975 + t929 * t968;
t1012 = (-qJD(5) - t1069) * t977;
t1011 = t939 * t1076 + t1103 * t979;
t1010 = -qJD(1) * t935 + t1103 * t972;
t1008 = -t970 * t1080 - t968 * t977;
t1009 = qJD(6) * t1053 - t1008 * t1072 - (t975 * t1012 + (t1030 * t968 + t975 * t1067) * t970) * t971;
t1007 = -t970 * t1091 + t975 * t977;
t1005 = -t1006 * t1065 + t905 * t969;
t1004 = -t940 * t1065 - t907 * t969;
t1001 = -t970 * t1072 - t1038;
t913 = -t1006 * t1080 - t945 * t968;
t996 = -qJD(6) * t913 + t1005;
t911 = -t940 * t1080 - t941 * t968;
t995 = -qJD(6) * t911 + t1004;
t991 = -qJD(7) * (-t928 * t1081 + t929 * t967) + t928 * t1062 + t865 * t968;
t990 = -qJD(7) * (-t922 * t1081 + t923 * t967) + t922 * t1062 + t867 * t968;
t989 = -qJD(7) * (-t938 * t1081 + t939 * t967) + t938 * t1062 + t903 * t968;
t986 = qJD(5) * t939 - t1001;
t985 = -qJD(5) * t943 - t977 * t1068 + t970 * t1073;
t984 = -t1028 * t974 + t1097 * t967;
t983 = t1001 * t969 - t970 * t1037;
t931 = t1008 * t971;
t982 = -qJD(6) * t931 + t983;
t980 = (r_i_i_C(1) * t966 + r_i_i_C(2) * t973) * t974 * qJD(7) + (t1028 * t967 + t1097 * t974) * qJD(6);
t930 = t1007 * t971;
t926 = t975 * t1087 + t943 * t968;
t920 = t939 * t968 + t1052;
t917 = -t975 * t1051 - t935 * t968;
t915 = -t972 * t1052 - t933 * t968;
t914 = -t967 * t1053 + t931 * t974;
t912 = -t1006 * t1091 + t945 * t975;
t910 = -t940 * t1091 + t941 * t975;
t901 = t938 * qJD(2) - qJD(4) * t1050 + (t1039 - t1066) * t969;
t890 = -t921 * t967 + t938 * t974;
t885 = -t1006 * t1092 + t913 * t974;
t884 = -t940 * t1092 + t911 * t974;
t880 = t1003 * t969 * t977 + t976 * t1047 + (t1037 * t977 + t1100) * t972;
t878 = t938 * t1076 + (-t1032 * t1085 + (-t1013 + t1041) * t969) * t979;
t876 = t1007 * t1072 + (-t1030 * t1088 + (t970 * t1067 + t1012) * t968) * t971;
t874 = -t898 * t967 + t928 * t974;
t864 = t1101 * t975 - t986 * t968;
t863 = t1101 * t968 + t975 * t986;
t862 = t1025 * t975 + t985 * t968;
t861 = t1025 * t968 - t985 * t975;
t860 = t1010 * t975 + t1098 * t968;
t859 = t1010 * t968 - t1098 * t975;
t858 = t1011 * t975 - t1099 * t968;
t857 = t1011 * t968 + t1099 * t975;
t856 = -t1009 * t974 + t982 * t967;
t853 = -t1023 * t975 + t997 * t968;
t851 = t1024 * t975 + t998 * t968;
t844 = -t897 * qJD(5) + t866 * t975 + t905 * t968;
t843 = t898 * qJD(5) + t866 * t968 - t905 * t975;
t842 = t903 * t967 - t921 * t1060 + (qJD(6) * t938 + t864) * t974;
t841 = -t891 * qJD(6) - t864 * t967 + t903 * t974;
t840 = t862 * t974 - t901 * t967 + (-t927 * t967 + t942 * t974) * qJD(6);
t838 = -t1021 * t974 + t995 * t967;
t836 = t1022 * t974 + t996 * t967;
t830 = t860 * t974 - t880 * t967 + (-t916 * t967 - t932 * t974) * qJD(6);
t828 = t858 * t974 + t878 * t967 + (-t918 * t967 - t934 * t974) * qJD(6);
t826 = t1111 * qJD(6) + t848 * t974 - t1110;
t824 = t1110 - t893 * t1060 + (qJD(6) * t922 - t848) * t974;
t822 = t865 * t967 - t898 * t1060 + (qJD(6) * t928 + t844) * t974;
t821 = -t875 * qJD(6) - t844 * t967 + t865 * t974;
t820 = t822 * t973 - t843 * t966 + (-t875 * t966 - t897 * t973) * qJD(7);
t819 = -t822 * t966 - t843 * t973 + (-t875 * t973 + t897 * t966) * qJD(7);
t1 = [(t826 * t973 + t1118) * r_i_i_C(1) + (-t826 * t966 + t1117) * r_i_i_C(2) - t867 * pkin(3) + t1097 * t1114 + (t1116 * r_i_i_C(1) + t1115 * r_i_i_C(2)) * qJD(7) + t1003 * pkin(2) (t828 * t973 - t857 * t966) * r_i_i_C(1) + (-t828 * t966 - t857 * t973) * r_i_i_C(2) + t878 * pkin(3) + t1097 * (-t887 * qJD(6) - t858 * t967 + t878 * t974) + ((-t887 * t966 - t917 * t973) * r_i_i_C(1) + (-t887 * t973 + t917 * t966) * r_i_i_C(2)) * qJD(7) + (t972 * t1075 + t1044) * pkin(2) (t836 * t973 - t851 * t966) * r_i_i_C(1) + (-t836 * t966 - t851 * t973) * r_i_i_C(2) + t1097 * (-t1022 * t967 + t996 * t974) + ((-t885 * t966 - t912 * t973) * r_i_i_C(1) + (-t885 * t973 + t912 * t966) * r_i_i_C(2)) * qJD(7) + t1005 * pkin(3), t866 * pkin(3) + t1097 * (t1020 * t974 - t994 * t967) + (t1017 * r_i_i_C(1) + t991 * r_i_i_C(2)) * t973 + (t991 * r_i_i_C(1) - t1017 * r_i_i_C(2)) * t966 (-t898 * t1055 - t844 * t966) * r_i_i_C(1) + (t898 * t1057 - t844 * t973) * r_i_i_C(2) + t984 * t843 + t980 * t897, -t1097 * t822 + (-t874 * t1055 - t821 * t966) * r_i_i_C(2) + (-t874 * t1057 + t821 * t973) * r_i_i_C(1), r_i_i_C(1) * t819 - t820 * r_i_i_C(2); t1002 * pkin(2) + t865 * pkin(3) + t820 * r_i_i_C(1) + t819 * r_i_i_C(2) + t1097 * t821 (t830 * t973 - t859 * t966) * r_i_i_C(1) + (-t830 * t966 - t859 * t973) * r_i_i_C(2) - t880 * pkin(3) + t1097 * (-t886 * qJD(6) - t860 * t967 - t880 * t974) + ((-t886 * t966 - t915 * t973) * r_i_i_C(1) + (-t886 * t973 + t915 * t966) * r_i_i_C(2)) * qJD(7) + (t1046 - t1047) * pkin(2) (t838 * t973 - t853 * t966) * r_i_i_C(1) + (-t838 * t966 - t853 * t973) * r_i_i_C(2) + t1097 * (t1021 * t967 + t995 * t974) + ((-t884 * t966 - t910 * t973) * r_i_i_C(1) + (-t884 * t973 + t910 * t966) * r_i_i_C(2)) * qJD(7) + t1004 * pkin(3), t868 * pkin(3) + t1097 * (t1019 * t974 - t993 * t967) + (t1016 * r_i_i_C(1) + t990 * r_i_i_C(2)) * t973 + (t990 * r_i_i_C(1) - t1016 * r_i_i_C(2)) * t966 (-t893 * t1055 + t848 * t966) * r_i_i_C(1) + (t893 * t1057 + t848 * t973) * r_i_i_C(2) + t984 * t845 + t980 * t892, -t1097 * t824 + (t1055 * t1111 + t1114 * t966) * r_i_i_C(2) + (t1057 * t1111 - t1114 * t973) * r_i_i_C(1) (-t824 * t966 - t1117) * r_i_i_C(1) + (-t824 * t973 + t1118) * r_i_i_C(2) + (-t1115 * r_i_i_C(1) + t1116 * r_i_i_C(2)) * qJD(7); 0 (t840 * t973 - t861 * t966) * r_i_i_C(1) + (-t840 * t966 - t861 * t973) * r_i_i_C(2) - t901 * pkin(3) - pkin(2) * t1072 + t1097 * (-t896 * qJD(6) - t862 * t967 - t901 * t974) + ((-t896 * t966 - t926 * t973) * r_i_i_C(1) + (-t896 * t973 + t926 * t966) * r_i_i_C(2)) * qJD(7) (t856 * t973 - t876 * t966) * r_i_i_C(1) + (-t856 * t966 - t876 * t973) * r_i_i_C(2) + t1097 * (t1009 * t967 + t982 * t974) + ((-t914 * t966 - t930 * t973) * r_i_i_C(1) + (-t914 * t973 + t930 * t966) * r_i_i_C(2)) * qJD(7) + t983 * pkin(3), t904 * pkin(3) + t1097 * (t1018 * t974 - t992 * t967) + (t1015 * r_i_i_C(1) + t989 * r_i_i_C(2)) * t973 + (t989 * r_i_i_C(1) - t1015 * r_i_i_C(2)) * t966 (-t921 * t1055 - t864 * t966) * r_i_i_C(1) + (t921 * t1057 - t864 * t973) * r_i_i_C(2) + t984 * t863 + t980 * t920, -t1097 * t842 + (-t890 * t1055 - t841 * t966) * r_i_i_C(2) + (-t890 * t1057 + t841 * t973) * r_i_i_C(1) (-t842 * t966 - t863 * t973) * r_i_i_C(1) + (-t842 * t973 + t863 * t966) * r_i_i_C(2) + ((-t891 * t973 + t920 * t966) * r_i_i_C(1) + (t891 * t966 + t920 * t973) * r_i_i_C(2)) * qJD(7);];
JaD_transl  = t1;
