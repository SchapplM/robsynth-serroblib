% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10_jacobiaD_transl_5_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobiaD_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:30
% EndTime: 2018-11-23 11:27:34
% DurationCPUTime: 3.84s
% Computational Cost: add. (10993->339), mult. (13147->516), div. (0->0), fcn. (10685->28), ass. (0->168)
t1035 = sin(qJ(5));
t1040 = cos(qJ(5));
t1081 = pkin(8) + qJ(4);
t1064 = cos(t1081) / 0.2e1;
t1082 = pkin(8) - qJ(4);
t1075 = cos(t1082);
t1001 = t1064 - t1075 / 0.2e1;
t1041 = cos(qJ(4));
t1083 = pkin(7) + qJ(3);
t1066 = cos(t1083) / 0.2e1;
t1084 = pkin(7) - qJ(3);
t1076 = cos(t1084);
t1003 = t1066 - t1076 / 0.2e1;
t1031 = sin(pkin(6));
t1042 = cos(qJ(3));
t1044 = cos(qJ(1));
t1037 = sin(qJ(3));
t1103 = qJD(3) * t1037;
t1039 = sin(qJ(1));
t1105 = qJD(1) * t1039;
t1085 = pkin(6) + qJ(2);
t1068 = cos(t1085) / 0.2e1;
t1086 = pkin(6) - qJ(2);
t1077 = cos(t1086);
t1052 = t1077 / 0.2e1 + t1068;
t1038 = sin(qJ(2));
t1088 = t1044 * t1038;
t1046 = t1039 * t1052 + t1088;
t1063 = sin(t1085) / 0.2e1;
t1074 = sin(t1086);
t1000 = t1063 - t1074 / 0.2e1;
t1051 = qJD(2) * t1000;
t1043 = cos(qJ(2));
t1089 = t1039 * t1043;
t966 = t1046 * qJD(1) + qJD(2) * t1089 + t1044 * t1051;
t1090 = t1039 * t1038;
t1017 = qJD(2) * t1090;
t1091 = t1039 * t1000;
t992 = t1052 * qJD(2);
t967 = -t1017 - qJD(1) * t1091 + (qJD(1) * t1043 + t992) * t1044;
t975 = -t1044 * t1052 + t1090;
t976 = t1044 * t1000 + t1089;
t1072 = sin(t1083);
t1061 = t1072 / 0.2e1;
t1073 = sin(t1084);
t1062 = t1073 / 0.2e1;
t997 = t1061 + t1062;
t987 = t997 * qJD(3);
t1004 = t1076 / 0.2e1 + t1066;
t989 = t1004 * qJD(3);
t998 = t1061 - t1073 / 0.2e1;
t1045 = (t1003 * t1105 + t1044 * t987) * t1031 - t967 * t1042 + t976 * t1103 + t966 * t998 + t975 * t989;
t1036 = sin(qJ(4));
t1101 = qJD(4) * t1036;
t1102 = qJD(3) * t1042;
t988 = (t1062 - t1072 / 0.2e1) * qJD(3);
t990 = t1003 * qJD(3);
t906 = -t966 * t1004 + (-t1044 * t990 + t997 * t1105) * t1031 - t967 * t1037 - t976 * t1102 - t975 * t988;
t1092 = t1031 * t1044;
t944 = t1003 * t1092 + t976 * t1042 - t975 * t998;
t946 = t975 * t1004 + t976 * t1037 + t997 * t1092;
t1033 = cos(pkin(7));
t1093 = t1031 * t1039;
t1079 = t1033 * t1093;
t1030 = sin(pkin(7));
t1106 = t966 * t1030;
t953 = qJD(1) * t1079 + t1106;
t1078 = t1033 * t1092;
t971 = -t975 * t1030 + t1078;
t1070 = sin(t1081);
t1059 = t1070 / 0.2e1;
t1071 = sin(t1082);
t1060 = t1071 / 0.2e1;
t995 = t1059 + t1060;
t983 = t995 * qJD(4);
t1002 = t1075 / 0.2e1 + t1064;
t985 = t1002 * qJD(4);
t996 = t1059 - t1071 / 0.2e1;
t881 = t953 * t1001 + t1041 * t1045 + t944 * t1101 - t906 * t996 + t946 * t985 + t971 * t983;
t1029 = sin(pkin(8));
t1032 = cos(pkin(8));
t899 = t906 * t1029 - t953 * t1032;
t1117 = t899 * t1035 + t881 * t1040;
t1116 = t881 * t1035 - t899 * t1040;
t912 = t971 * t1001 + t944 * t1041 - t946 * t996;
t929 = t946 * t1029 - t971 * t1032;
t1115 = t1035 * t929 + t1040 * t912;
t1114 = t1035 * t912 - t1040 * t929;
t1100 = qJD(4) * t1041;
t984 = (t1060 - t1070 / 0.2e1) * qJD(4);
t986 = t1001 * qJD(4);
t1113 = t906 * t1002 + t1036 * t1045 - t944 * t1100 - t946 * t984 + t953 * t995 - t971 * t986;
t1005 = t1068 - t1077 / 0.2e1;
t964 = t976 * qJD(1) + qJD(2) * t1088 + t1039 * t992;
t1104 = qJD(1) * t1044;
t1087 = t1044 * t1043;
t963 = qJD(1) * t975 - qJD(2) * t1087 + t1039 * t1051;
t980 = -t1087 + t1091;
t905 = (t1003 * t1104 - t1039 * t987) * t1031 + t964 * t1042 - t980 * t1103 - t963 * t998 + t1046 * t989;
t1109 = r_i_i_C(3) + pkin(13);
t1108 = pkin(11) * t1030;
t1107 = t963 * t1030;
t1099 = qJD(5) * t1035;
t1098 = qJD(5) * t1040;
t1097 = t1001 * t1030;
t1096 = t1005 * t1037;
t1095 = t1005 * t1042;
t1094 = t1030 * t1032;
t1080 = pkin(11) * t1033 + pkin(10);
t903 = t980 * t1102 + t963 * t1004 + t964 * t1037 - t1046 * t988 + (t1039 * t990 + t997 * t1104) * t1031;
t951 = qJD(1) * t1078 - t1107;
t897 = -t903 * t1029 + t951 * t1032;
t1058 = t1040 * r_i_i_C(1) - t1035 * r_i_i_C(2) + pkin(4);
t1057 = qJD(5) * (-t1035 * r_i_i_C(1) - t1040 * r_i_i_C(2));
t916 = t964 * t1004 - t963 * t1037 + t1046 * t1102 + t980 * t988;
t900 = -t916 * t1029 - t964 * t1094;
t968 = t980 * qJD(1) - t1044 * t992 + t1017;
t918 = t968 * t1004 + t966 * t1037 + t975 * t1102 - t976 * t988;
t901 = -t918 * t1029 - t968 * t1094;
t999 = t1063 + t1074 / 0.2e1;
t991 = t999 * qJD(2);
t993 = t1005 * qJD(2);
t937 = -t991 * t1004 + t1005 * t988 - t993 * t1037 - t999 * t1102;
t932 = -t937 * t1029 + t991 * t1094;
t1034 = cos(pkin(6));
t961 = t1034 * t1003 - t999 * t998 + t1095;
t950 = t1003 * t1093 + t1042 * t980 + t1046 * t998;
t934 = qJD(3) * t1096 + t1034 * t987 + t991 * t1042 + t999 * t989 + t993 * t998;
t948 = -t1004 * t1046 + t1037 * t980 + t997 * t1093;
t973 = t1030 * t1046 + t1079;
t875 = -t951 * t1001 - t1041 * t905 + t1101 * t950 + t903 * t996 + t948 * t985 + t973 * t983;
t933 = qJD(3) * t1095 + t993 * t1004 + t1034 * t990 - t991 * t1037 + t999 * t988;
t959 = t999 * t1004 + t1034 * t997 + t1096;
t974 = -t999 * t1030 + t1034 * t1033;
t891 = t934 * t1041 + t993 * t1097 + t1101 * t961 + t933 * t996 + t959 * t985 + t974 * t983;
t970 = t1005 * t998 + t999 * t1042;
t969 = t1005 * t1004 - t999 * t1037;
t958 = -t1042 * t1046 + t980 * t998;
t957 = t980 * t1004 + t1037 * t1046;
t956 = -t975 * t1042 - t976 * t998;
t955 = -t1004 * t976 + t975 * t1037;
t954 = -t1005 * t1094 - t969 * t1029;
t941 = -t959 * t1029 + t974 * t1032;
t940 = -t957 * t1029 - t980 * t1094;
t939 = -t955 * t1029 + t1094 * t976;
t938 = t1005 * t989 + t993 * t1042 - t999 * t1103 - t991 * t998;
t931 = -t948 * t1029 + t973 * t1032;
t928 = t1005 * t1097 + t970 * t1041 + t969 * t996;
t927 = t959 * t1041 + t961 * t996;
t926 = -t933 * t1029 - t993 * t1094;
t925 = -t974 * t1001 - t1041 * t961 + t959 * t996;
t923 = t948 * t1041 + t950 * t996;
t922 = -t1041 * t946 - t944 * t996;
t921 = t958 * t1041 + t980 * t1097 + t957 * t996;
t920 = t956 * t1041 - t1097 * t976 + t955 * t996;
t919 = -t966 * t1042 + t975 * t1103 + t968 * t998 - t976 * t989;
t917 = t963 * t1042 + t1046 * t1103 + t964 * t998 + t980 * t989;
t915 = -t973 * t1001 - t1041 * t950 + t948 * t996;
t896 = -t970 * t1101 + t938 * t1041 + t937 * t996 + t969 * t985 + (-t1001 * t991 - t1005 * t983) * t1030;
t894 = t933 * t1041 - t959 * t1101 - t934 * t996 + t961 * t985;
t889 = -t956 * t1101 + t919 * t1041 + t918 * t996 + t955 * t985 + (t1001 * t968 + t976 * t983) * t1030;
t887 = -t958 * t1101 + t917 * t1041 + t916 * t996 + t957 * t985 + (t1001 * t964 - t980 * t983) * t1030;
t885 = t906 * t1041 + t1045 * t996 + t1101 * t946 - t944 * t985;
t883 = t903 * t1041 - t948 * t1101 + t905 * t996 + t950 * t985;
t874 = -t903 * t1002 - t1036 * t905 - t1100 * t950 - t948 * t984 - t951 * t995 - t973 * t986;
t873 = t897 * t1035 + t875 * t1040 + (-t1035 * t915 + t1040 * t931) * qJD(5);
t872 = -t875 * t1035 + t897 * t1040 + (-t1035 * t931 - t1040 * t915) * qJD(5);
t1 = [t1117 * r_i_i_C(1) - t1116 * r_i_i_C(2) + t881 * pkin(4) + t1045 * pkin(3) - t967 * pkin(2) - pkin(11) * t1106 + t1109 * t1113 + (t1114 * r_i_i_C(1) + r_i_i_C(2) * t1115) * qJD(5) + t899 * pkin(12) + (-t1044 * pkin(1) - t1080 * t1093) * qJD(1) (t900 * t1035 + t887 * t1040) * r_i_i_C(1) + (-t887 * t1035 + t900 * t1040) * r_i_i_C(2) + t887 * pkin(4) + t917 * pkin(3) + t963 * pkin(2) - t964 * t1108 + t1109 * (t958 * t1100 - t916 * t1002 + t917 * t1036 - t957 * t984 + (t964 * t995 + t980 * t986) * t1030) + ((-t1035 * t921 + t1040 * t940) * r_i_i_C(1) + (-t1035 * t940 - t1040 * t921) * r_i_i_C(2)) * qJD(5) + t900 * pkin(12) (t883 * t1040 - t923 * t1099) * r_i_i_C(1) + (-t883 * t1035 - t923 * t1098) * r_i_i_C(2) + t883 * pkin(4) + t903 * pkin(3) + t1109 * (-t905 * t1002 + t903 * t1036 + t948 * t1100 - t950 * t984) + ((-t1035 * t905 - t950 * t1098) * r_i_i_C(1) + (-t1040 * t905 + t950 * t1099) * r_i_i_C(2) - t905 * pkin(12)) * t1029, t1109 * t875 + (t948 * t1002 + t1036 * t950 + t973 * t995) * t1057 - t1058 * t874, t872 * r_i_i_C(1) - t873 * r_i_i_C(2), 0; -pkin(11) * t1107 - t964 * pkin(2) - t905 * pkin(3) + t875 * pkin(4) + t873 * r_i_i_C(1) + t872 * r_i_i_C(2) + t1109 * t874 + t897 * pkin(12) + (-pkin(1) * t1039 + t1080 * t1092) * qJD(1) (t901 * t1035 + t889 * t1040) * r_i_i_C(1) + (-t889 * t1035 + t901 * t1040) * r_i_i_C(2) + t889 * pkin(4) + t919 * pkin(3) - t966 * pkin(2) - t968 * t1108 + t1109 * (t956 * t1100 - t918 * t1002 + t919 * t1036 - t955 * t984 + (t968 * t995 - t976 * t986) * t1030) + ((-t1035 * t920 + t1040 * t939) * r_i_i_C(1) + (-t1035 * t939 - t1040 * t920) * r_i_i_C(2)) * qJD(5) + t901 * pkin(12) (t885 * t1040 - t922 * t1099) * r_i_i_C(1) + (-t885 * t1035 - t922 * t1098) * r_i_i_C(2) + t885 * pkin(4) + t906 * pkin(3) + t1109 * (-t1002 * t1045 + t906 * t1036 - t1100 * t946 + t944 * t984) + ((-t1035 * t1045 + t1098 * t944) * r_i_i_C(1) + (-t1040 * t1045 - t1099 * t944) * r_i_i_C(2) - t1045 * pkin(12)) * t1029, -t1109 * t881 + (-t1002 * t946 - t944 * t1036 - t971 * t995) * t1057 + t1058 * t1113, t1116 * r_i_i_C(1) + t1117 * r_i_i_C(2) + (-r_i_i_C(1) * t1115 + t1114 * r_i_i_C(2)) * qJD(5), 0; 0 (t932 * t1035 + t896 * t1040) * r_i_i_C(1) + (-t896 * t1035 + t932 * t1040) * r_i_i_C(2) + t896 * pkin(4) + t938 * pkin(3) + t993 * pkin(2) + t991 * t1108 + t1109 * (t970 * t1100 - t937 * t1002 + t938 * t1036 - t969 * t984 + (t1005 * t986 - t991 * t995) * t1030) + ((-t1035 * t928 + t1040 * t954) * r_i_i_C(1) + (-t1035 * t954 - t1040 * t928) * r_i_i_C(2)) * qJD(5) + t932 * pkin(12) (t894 * t1040 - t927 * t1099) * r_i_i_C(1) + (-t894 * t1035 - t927 * t1098) * r_i_i_C(2) + t894 * pkin(4) + t933 * pkin(3) + t1109 * (t1002 * t934 + t933 * t1036 + t959 * t1100 - t961 * t984) + ((t1035 * t934 - t961 * t1098) * r_i_i_C(1) + (t1040 * t934 + t961 * t1099) * r_i_i_C(2) + t934 * pkin(12)) * t1029, t1109 * t891 + (t959 * t1002 + t1036 * t961 + t974 * t995) * t1057 + t1058 * (-t993 * t1030 * t995 + t933 * t1002 - t934 * t1036 + t1100 * t961 + t959 * t984 + t974 * t986) (-t891 * t1035 + t926 * t1040) * r_i_i_C(1) + (-t926 * t1035 - t891 * t1040) * r_i_i_C(2) + ((-t1035 * t941 - t1040 * t925) * r_i_i_C(1) + (t1035 * t925 - t1040 * t941) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
