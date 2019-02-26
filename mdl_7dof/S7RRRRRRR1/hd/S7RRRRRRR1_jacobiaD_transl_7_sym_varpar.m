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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_7_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_7_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_7_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:44
% EndTime: 2019-02-26 22:54:50
% DurationCPUTime: 6.02s
% Computational Cost: add. (2918->345), mult. (8814->599), div. (0->0), fcn. (9592->14), ass. (0->212)
t978 = cos(qJ(2));
t1070 = qJD(2) * t978;
t979 = cos(qJ(1));
t1072 = qJD(1) * t979;
t971 = sin(qJ(2));
t972 = sin(qJ(1));
t1003 = t972 * t1070 + t971 * t1072;
t1064 = qJD(4) * t971;
t1073 = qJD(1) * t978;
t1033 = qJD(3) + t1073;
t1071 = qJD(2) * t971;
t1038 = t972 * t1071;
t1066 = qJD(3) * t978;
t970 = sin(qJ(3));
t1042 = t970 * t1066;
t1074 = qJD(1) * t972;
t977 = cos(qJ(3));
t1075 = t979 * t977;
t908 = t1033 * t1075 - t977 * t1038 - t972 * t1042 - t970 * t1074;
t1076 = t979 * t970;
t1081 = t972 * t978;
t941 = t977 * t1081 + t1076;
t969 = sin(qJ(4));
t976 = cos(qJ(4));
t868 = (t972 * t1064 + t908) * t976 + (-qJD(4) * t941 + t1003) * t969;
t1085 = t971 * t969;
t923 = t972 * t1085 + t941 * t976;
t940 = t970 * t1081 - t1075;
t968 = sin(qJ(5));
t975 = cos(qJ(5));
t893 = t923 * t975 - t940 * t968;
t1006 = t978 * t1076 + t972 * t977;
t907 = t1006 * qJD(1) + t941 * qJD(3) - t970 * t1038;
t845 = t893 * qJD(5) + t868 * t968 + t907 * t975;
t966 = sin(qJ(7));
t1116 = t845 * t966;
t973 = cos(qJ(7));
t1115 = t845 * t973;
t1084 = t971 * t976;
t922 = -t972 * t1084 + t941 * t969;
t967 = sin(qJ(6));
t974 = cos(qJ(6));
t872 = t893 * t974 + t922 * t967;
t892 = t923 * t968 + t940 * t975;
t1114 = t872 * t966 + t892 * t973;
t1113 = t872 * t973 - t892 * t966;
t848 = t892 * qJD(5) - t868 * t975 + t907 * t968;
t867 = t923 * qJD(4) - t1003 * t976 + t908 * t969;
t1112 = t872 * qJD(6) - t848 * t967 - t867 * t974;
t1109 = t893 * t967 - t922 * t974;
t1108 = t867 * t967;
t1032 = qJD(4) * t977 - qJD(2);
t1014 = t1032 * t969;
t1067 = qJD(3) * t976;
t1106 = (t970 * t1067 + t1014) * t971;
t1031 = qJD(2) * t977 - qJD(4);
t1013 = t1031 * t978;
t1061 = qJD(5) * t970;
t1040 = t971 * t1061;
t1101 = -t976 * t1013 + t1040 + t1106;
t1077 = t978 * t976;
t904 = t1031 * t1077 - t1106;
t1099 = t904 - t1040;
t1068 = qJD(3) * t971;
t1044 = t970 * t1068;
t1065 = qJD(4) * t969;
t1098 = -t969 * t1044 - t978 * t1065;
t1069 = qJD(2) * t979;
t1002 = -t978 * t1069 + t971 * t1074;
t1097 = t978 * t1072 - t1038;
t1041 = t977 * t1068;
t1078 = t978 * t969;
t1083 = t971 * t977;
t939 = t976 * t1083 - t1078;
t935 = t939 * t979;
t1096 = -qJD(5) * t935 + t1002 * t970 - t979 * t1041;
t933 = t939 * t972;
t1095 = qJD(5) * t933 + t1003 * t970 + t972 * t1041;
t1094 = r_i_i_C(3) + pkin(4);
t1089 = t967 * t969;
t1088 = t968 * t976;
t1087 = t970 * t975;
t1086 = t970 * t978;
t1082 = t971 * t979;
t1080 = t974 * t975;
t1079 = t975 * t976;
t1063 = qJD(4) * t976;
t1062 = qJD(5) * t968;
t1060 = qJD(5) * t975;
t1059 = qJD(5) * t976;
t1058 = qJD(6) * t967;
t1057 = qJD(6) * t969;
t1056 = qJD(6) * t975;
t1055 = qJD(7) * t966;
t1054 = qJD(7) * t968;
t1053 = qJD(7) * t973;
t1052 = t970 * t1085;
t1051 = t971 * t970 * t968;
t1050 = t971 * t1087;
t1049 = t977 * t1077;
t1039 = t971 * t1063;
t1035 = t971 * t1069;
t1034 = qJD(1) + t1066;
t1030 = qJD(3) + t1059;
t1028 = t973 * r_i_i_C(1) - t966 * r_i_i_C(2);
t1025 = (-t1031 * t971 - t1042) * t976 + (-t1061 - t1014) * t978;
t988 = t1033 * t972 + t1035;
t906 = t1034 * t1076 + t988 * t977;
t1024 = -t1006 * t1059 - t906;
t1023 = t940 * t1059 - t908;
t905 = -t1034 * t1075 + t988 * t970;
t945 = t978 * t1075 - t972 * t970;
t998 = -qJD(5) * t945 + t1006 * t1065 + t905 * t976;
t1022 = -t1006 * t1057 - t1024 * t968 + t998 * t975;
t997 = -qJD(5) * t941 + t940 * t1065 - t907 * t976;
t1021 = -t1023 * t968 + t940 * t1057 - t997 * t975;
t866 = (t979 * t1064 - t906) * t976 + (-qJD(4) * t945 - t1002) * t969;
t928 = -t976 * t1082 + t945 * t969;
t1020 = t928 * t1056 + t866;
t1019 = t922 * t1056 + t868;
t938 = t969 * t1083 + t1077;
t1018 = t938 * t1056 + t904;
t929 = t969 * t1082 + t945 * t976;
t865 = t929 * qJD(4) + t1002 * t976 - t906 * t969;
t994 = qJD(6) * t929 + t928 * t1062 - t865 * t975;
t1017 = t1020 * t967 + t928 * t1054 + t994 * t974;
t993 = qJD(6) * t923 + t922 * t1062 - t867 * t975;
t1016 = t1019 * t967 + t922 * t1054 + t993 * t974;
t903 = -t976 * t1071 + (t969 * t1070 + t1039) * t977 + t1098;
t992 = qJD(6) * t939 + t938 * t1062 - t903 * t975;
t1015 = t1018 * t967 + t938 * t1054 + t992 * t974;
t898 = -t1006 * t968 + t929 * t975;
t875 = t898 * t974 + t928 * t967;
t916 = t972 * t1051 - t933 * t975;
t932 = t938 * t972;
t886 = t916 * t974 - t932 * t967;
t918 = t979 * t1051 - t935 * t975;
t934 = t938 * t979;
t887 = t918 * t974 - t934 * t967;
t921 = t939 * t975 - t1051;
t891 = t921 * t974 + t938 * t967;
t943 = t1049 + t1085;
t927 = -t968 * t1086 + t943 * t975;
t942 = t977 * t1078 - t1084;
t896 = t927 * t974 + t942 * t967;
t897 = t1006 * t975 + t929 * t968;
t1012 = (-qJD(5) - t1067) * t977;
t1011 = t939 * t1074 + t1101 * t979;
t1010 = -qJD(1) * t935 + t1101 * t972;
t1008 = -t970 * t1079 - t968 * t977;
t1009 = qJD(6) * t1052 - t1008 * t1070 - (t975 * t1012 + (t1030 * t968 + t975 * t1065) * t970) * t971;
t1007 = -t970 * t1088 + t975 * t977;
t1005 = -t1006 * t1063 + t905 * t969;
t1004 = -t940 * t1063 - t907 * t969;
t1001 = -t970 * t1070 - t1041;
t913 = -t1006 * t1079 - t945 * t968;
t996 = -qJD(6) * t913 + t1005;
t911 = -t940 * t1079 - t941 * t968;
t995 = -qJD(6) * t911 + t1004;
t991 = -qJD(7) * (-t928 * t1080 + t929 * t967) + t928 * t1060 + t865 * t968;
t990 = -qJD(7) * (-t922 * t1080 + t923 * t967) + t922 * t1060 + t867 * t968;
t989 = -qJD(7) * (-t938 * t1080 + t939 * t967) + t938 * t1060 + t903 * t968;
t986 = qJD(5) * t939 - t1001;
t985 = -qJD(5) * t943 - t977 * t1066 + t970 * t1071;
t984 = -t1028 * t974 + t1094 * t967;
t983 = t1001 * t969 - t970 * t1039;
t931 = t1008 * t971;
t982 = -qJD(6) * t931 + t983;
t980 = (t966 * r_i_i_C(1) + t973 * r_i_i_C(2)) * t974 * qJD(7) + (t1028 * t967 + t1094 * t974) * qJD(6);
t930 = t1007 * t971;
t926 = t975 * t1086 + t943 * t968;
t920 = t939 * t968 + t1050;
t917 = -t979 * t1050 - t935 * t968;
t915 = -t972 * t1050 - t933 * t968;
t914 = -t967 * t1052 + t931 * t974;
t912 = -t1006 * t1088 + t945 * t975;
t910 = -t940 * t1088 + t941 * t975;
t901 = t938 * qJD(2) - qJD(4) * t1049 + (t1042 - t1064) * t969;
t890 = -t921 * t967 + t938 * t974;
t885 = -t1006 * t1089 + t913 * t974;
t884 = -t940 * t1089 + t911 * t974;
t880 = t1003 * t969 * t977 + t1097 * t976 + (t1039 * t977 + t1098) * t972;
t878 = t938 * t1074 + (-t1032 * t1084 + (-t1013 + t1044) * t969) * t979;
t876 = t1007 * t1070 + (-t1030 * t1087 + (t970 * t1065 + t1012) * t968) * t971;
t874 = -t898 * t967 + t928 * t974;
t864 = t1099 * t975 - t986 * t968;
t863 = t1099 * t968 + t986 * t975;
t862 = t1025 * t975 + t985 * t968;
t861 = t1025 * t968 - t985 * t975;
t860 = t1010 * t975 + t1095 * t968;
t859 = t1010 * t968 - t1095 * t975;
t858 = t1011 * t975 - t1096 * t968;
t857 = t1011 * t968 + t1096 * t975;
t856 = -t1009 * t974 + t982 * t967;
t853 = -t1023 * t975 + t997 * t968;
t851 = t1024 * t975 + t998 * t968;
t844 = -t897 * qJD(5) + t866 * t975 + t905 * t968;
t843 = t898 * qJD(5) + t866 * t968 - t905 * t975;
t842 = t903 * t967 - t921 * t1058 + (qJD(6) * t938 + t864) * t974;
t841 = -t891 * qJD(6) - t864 * t967 + t903 * t974;
t840 = t862 * t974 - t901 * t967 + (-t927 * t967 + t942 * t974) * qJD(6);
t838 = -t1021 * t974 + t995 * t967;
t836 = t1022 * t974 + t996 * t967;
t830 = t860 * t974 - t880 * t967 + (-t916 * t967 - t932 * t974) * qJD(6);
t828 = t858 * t974 + t878 * t967 + (-t918 * t967 - t934 * t974) * qJD(6);
t826 = qJD(6) * t1109 + t848 * t974 - t1108;
t824 = t1108 - t893 * t1058 + (qJD(6) * t922 - t848) * t974;
t822 = t865 * t967 - t898 * t1058 + (qJD(6) * t928 + t844) * t974;
t821 = -t875 * qJD(6) - t844 * t967 + t865 * t974;
t820 = t822 * t973 - t843 * t966 + (-t875 * t966 - t897 * t973) * qJD(7);
t819 = -t822 * t966 - t843 * t973 + (-t875 * t973 + t897 * t966) * qJD(7);
t1 = [(t826 * t973 + t1116) * r_i_i_C(1) + (-t826 * t966 + t1115) * r_i_i_C(2) - t867 * pkin(3) + t1094 * t1112 + (r_i_i_C(1) * t1114 + r_i_i_C(2) * t1113) * qJD(7) + t1003 * pkin(2) (t828 * t973 - t857 * t966) * r_i_i_C(1) + (-t828 * t966 - t857 * t973) * r_i_i_C(2) + t878 * pkin(3) + t1094 * (-t887 * qJD(6) - t858 * t967 + t878 * t974) + ((-t887 * t966 - t917 * t973) * r_i_i_C(1) + (-t887 * t973 + t917 * t966) * r_i_i_C(2)) * qJD(7) + (t972 * t1073 + t1035) * pkin(2) (t836 * t973 - t851 * t966) * r_i_i_C(1) + (-t836 * t966 - t851 * t973) * r_i_i_C(2) + t1094 * (-t1022 * t967 + t996 * t974) + ((-t885 * t966 - t912 * t973) * r_i_i_C(1) + (-t885 * t973 + t912 * t966) * r_i_i_C(2)) * qJD(7) + t1005 * pkin(3), t866 * pkin(3) + t1094 * (t1020 * t974 - t967 * t994) + (r_i_i_C(1) * t1017 + r_i_i_C(2) * t991) * t973 + (r_i_i_C(1) * t991 - r_i_i_C(2) * t1017) * t966 (-t1053 * t898 - t844 * t966) * r_i_i_C(1) + (t1055 * t898 - t844 * t973) * r_i_i_C(2) + t984 * t843 + t980 * t897, -t1094 * t822 + (-t1053 * t874 - t821 * t966) * r_i_i_C(2) + (-t1055 * t874 + t821 * t973) * r_i_i_C(1), r_i_i_C(1) * t819 - r_i_i_C(2) * t820; t1002 * pkin(2) + t865 * pkin(3) + t820 * r_i_i_C(1) + t819 * r_i_i_C(2) + t1094 * t821 (t830 * t973 - t859 * t966) * r_i_i_C(1) + (-t830 * t966 - t859 * t973) * r_i_i_C(2) - t880 * pkin(3) + t1094 * (-t886 * qJD(6) - t860 * t967 - t880 * t974) + ((-t886 * t966 - t915 * t973) * r_i_i_C(1) + (-t886 * t973 + t915 * t966) * r_i_i_C(2)) * qJD(7) - t1097 * pkin(2) (t838 * t973 - t853 * t966) * r_i_i_C(1) + (-t838 * t966 - t853 * t973) * r_i_i_C(2) + t1094 * (t1021 * t967 + t995 * t974) + ((-t884 * t966 - t910 * t973) * r_i_i_C(1) + (-t884 * t973 + t910 * t966) * r_i_i_C(2)) * qJD(7) + t1004 * pkin(3), t868 * pkin(3) + t1094 * (t1019 * t974 - t967 * t993) + (r_i_i_C(1) * t1016 + r_i_i_C(2) * t990) * t973 + (r_i_i_C(1) * t990 - r_i_i_C(2) * t1016) * t966 (-t1053 * t893 + t848 * t966) * r_i_i_C(1) + (t893 * t1055 + t848 * t973) * r_i_i_C(2) + t984 * t845 + t980 * t892, -t1094 * t824 + (t1053 * t1109 + t1112 * t966) * r_i_i_C(2) + (t1055 * t1109 - t1112 * t973) * r_i_i_C(1) (-t824 * t966 - t1115) * r_i_i_C(1) + (-t824 * t973 + t1116) * r_i_i_C(2) + (-r_i_i_C(1) * t1113 + r_i_i_C(2) * t1114) * qJD(7); 0 (t840 * t973 - t861 * t966) * r_i_i_C(1) + (-t840 * t966 - t861 * t973) * r_i_i_C(2) - t901 * pkin(3) - pkin(2) * t1070 + t1094 * (-t896 * qJD(6) - t862 * t967 - t901 * t974) + ((-t896 * t966 - t926 * t973) * r_i_i_C(1) + (-t896 * t973 + t926 * t966) * r_i_i_C(2)) * qJD(7) (t856 * t973 - t876 * t966) * r_i_i_C(1) + (-t856 * t966 - t876 * t973) * r_i_i_C(2) + t1094 * (t1009 * t967 + t982 * t974) + ((-t914 * t966 - t930 * t973) * r_i_i_C(1) + (-t914 * t973 + t930 * t966) * r_i_i_C(2)) * qJD(7) + t983 * pkin(3), t904 * pkin(3) + t1094 * (t1018 * t974 - t967 * t992) + (r_i_i_C(1) * t1015 + r_i_i_C(2) * t989) * t973 + (r_i_i_C(1) * t989 - r_i_i_C(2) * t1015) * t966 (-t1053 * t921 - t864 * t966) * r_i_i_C(1) + (t1055 * t921 - t864 * t973) * r_i_i_C(2) + t984 * t863 + t980 * t920, -t1094 * t842 + (-t1053 * t890 - t841 * t966) * r_i_i_C(2) + (-t1055 * t890 + t841 * t973) * r_i_i_C(1) (-t842 * t966 - t863 * t973) * r_i_i_C(1) + (-t842 * t973 + t863 * t966) * r_i_i_C(2) + ((-t891 * t973 + t920 * t966) * r_i_i_C(1) + (t891 * t966 + t920 * t973) * r_i_i_C(2)) * qJD(7);];
JaD_transl  = t1;
