% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:34
% EndTime: 2019-12-31 17:26:35
% DurationCPUTime: 1.46s
% Computational Cost: add. (5978->181), mult. (12460->246), div. (0->0), fcn. (8623->8), ass. (0->143)
t1057 = sin(qJ(3));
t1061 = cos(qJ(3));
t1062 = cos(qJ(2));
t1089 = qJD(1) * t1062;
t1058 = sin(qJ(2));
t1090 = qJD(1) * t1058;
t1023 = t1057 * t1090 - t1061 * t1089;
t1022 = qJD(4) + t1023;
t1097 = qJD(4) + t1022;
t1051 = t1062 * qJDD(1);
t1077 = qJD(2) * t1090;
t1032 = t1051 - t1077;
t1055 = t1062 ^ 2;
t1065 = qJD(1) ^ 2;
t1059 = sin(qJ(1));
t1063 = cos(qJ(1));
t1040 = t1059 * g(1) - t1063 * g(2);
t1070 = qJDD(1) * pkin(1) + t1040;
t1071 = qJD(2) * pkin(2) - pkin(6) * t1090;
t1004 = t1032 * pkin(2) + (t1055 * pkin(6) + pkin(5)) * t1065 - t1071 * t1090 + t1070;
t1025 = (t1057 * t1062 + t1058 * t1061) * qJD(1);
t1053 = qJD(2) + qJD(3);
t1056 = sin(qJ(4));
t1060 = cos(qJ(4));
t1010 = t1056 * t1025 - t1060 * t1053;
t1096 = t1010 ^ 2;
t1012 = t1060 * t1025 + t1056 * t1053;
t1095 = t1012 ^ 2;
t1094 = t1022 ^ 2;
t1093 = t1023 ^ 2;
t1092 = t1025 ^ 2;
t1091 = t1053 ^ 2;
t1041 = -t1063 * g(1) - t1059 * g(2);
t1027 = -t1065 * pkin(1) + qJDD(1) * pkin(5) + t1041;
t1019 = -t1058 * g(3) + t1062 * t1027;
t1086 = t1055 * t1065;
t1001 = -pkin(2) * t1086 + t1032 * pkin(6) - qJD(2) * t1071 + t1019;
t1076 = qJD(2) * t1089;
t1080 = t1058 * qJDD(1);
t1031 = t1076 + t1080;
t1084 = t1058 * t1065;
t1085 = t1058 * t1027;
t1066 = qJDD(2) * pkin(2) - t1031 * pkin(6) - t1085 + (qJD(2) * pkin(6) * qJD(1) + pkin(2) * t1084 - g(3)) * t1062;
t982 = t1061 * t1001 + t1057 * t1066;
t1088 = t1012 * t1010;
t1087 = t1025 * t1023;
t1083 = qJD(3) - t1053;
t1082 = qJD(4) - t1022;
t1054 = t1058 ^ 2;
t1081 = t1054 + t1055;
t1079 = qJDD(2) + qJDD(3);
t981 = -t1057 * t1001 + t1061 * t1066;
t1072 = -t1061 * t1031 - t1057 * t1032;
t1002 = -t1023 * qJD(3) - t1072;
t1075 = t1053 * t1023 - t1002;
t1074 = -t1056 * t1002 + t1060 * t1079;
t1073 = t1057 * t1031 - t1061 * t1032;
t1069 = -t1060 * t1002 - t1056 * t1079;
t1068 = -t1025 * qJD(3) - qJDD(4) - t1073;
t989 = (qJD(3) + t1053) * t1025 + t1073;
t1064 = qJD(2) ^ 2;
t1045 = t1062 * t1084;
t1043 = -t1064 - t1086;
t1042 = -t1054 * t1065 - t1064;
t1039 = -qJDD(2) + t1045;
t1038 = qJDD(2) + t1045;
t1037 = t1081 * t1065;
t1036 = -t1059 * qJDD(1) - t1063 * t1065;
t1035 = t1063 * qJDD(1) - t1059 * t1065;
t1034 = t1081 * qJDD(1);
t1033 = t1051 - 0.2e1 * t1077;
t1030 = 0.2e1 * t1076 + t1080;
t1026 = t1065 * pkin(5) + t1070;
t1018 = -t1062 * g(3) - t1085;
t1017 = -t1091 - t1092;
t1016 = t1062 * t1039 - t1058 * t1042;
t1015 = -t1058 * t1038 + t1062 * t1043;
t1014 = t1058 * t1039 + t1062 * t1042;
t1013 = t1062 * t1038 + t1058 * t1043;
t1009 = t1023 * pkin(3) - t1025 * pkin(7);
t1008 = -t1079 - t1087;
t1007 = t1079 - t1087;
t1006 = -t1091 - t1093;
t1003 = -t1092 - t1093;
t1000 = -t1058 * t1018 + t1062 * t1019;
t999 = t1062 * t1018 + t1058 * t1019;
t995 = -t1094 - t1095;
t994 = t1061 * t1008 - t1057 * t1017;
t993 = t1057 * t1008 + t1061 * t1017;
t992 = t1083 * t1023 + t1072;
t990 = -t1083 * t1025 - t1073;
t988 = -t1094 - t1096;
t987 = t1061 * t1006 - t1057 * t1007;
t986 = t1057 * t1006 + t1061 * t1007;
t985 = -t1095 - t1096;
t984 = t1068 - t1088;
t983 = -t1068 - t1088;
t980 = t1082 * t1010 + t1069;
t979 = -t1097 * t1010 - t1069;
t978 = -t1082 * t1012 + t1074;
t977 = t1097 * t1012 - t1074;
t976 = -t1058 * t993 + t1062 * t994;
t975 = t1058 * t994 + t1062 * t993;
t974 = -t1057 * t992 + t1061 * t990;
t973 = t1057 * t990 + t1061 * t992;
t972 = -t1058 * t986 + t1062 * t987;
t971 = t1058 * t987 + t1062 * t986;
t970 = -t1056 * t995 + t1060 * t984;
t969 = t1056 * t984 + t1060 * t995;
t968 = -t1091 * pkin(3) + t1079 * pkin(7) - t1023 * t1009 + t982;
t967 = -t1079 * pkin(3) - t1091 * pkin(7) + t1025 * t1009 - t981;
t966 = -t1056 * t983 + t1060 * t988;
t965 = t1056 * t988 + t1060 * t983;
t964 = t989 * pkin(3) + t1075 * pkin(7) - t1004;
t963 = -t1057 * t981 + t1061 * t982;
t962 = t1057 * t982 + t1061 * t981;
t961 = -t1056 * t980 + t1060 * t978;
t960 = t1056 * t978 + t1060 * t980;
t959 = -t1058 * t973 + t1062 * t974;
t958 = t1058 * t974 + t1062 * t973;
t957 = t1057 * t979 + t1061 * t970;
t956 = t1057 * t970 - t1061 * t979;
t955 = t1057 * t977 + t1061 * t966;
t954 = t1057 * t966 - t1061 * t977;
t953 = t1057 * t985 + t1061 * t961;
t952 = t1057 * t961 - t1061 * t985;
t951 = t1056 * t964 + t1060 * t968;
t950 = -t1056 * t968 + t1060 * t964;
t949 = -t1058 * t962 + t1062 * t963;
t948 = t1058 * t963 + t1062 * t962;
t947 = -t1058 * t956 + t1062 * t957;
t946 = t1058 * t957 + t1062 * t956;
t945 = -t1058 * t954 + t1062 * t955;
t944 = t1058 * t955 + t1062 * t954;
t943 = -t1058 * t952 + t1062 * t953;
t942 = t1058 * t953 + t1062 * t952;
t941 = -t1056 * t950 + t1060 * t951;
t940 = t1056 * t951 + t1060 * t950;
t939 = t1057 * t967 + t1061 * t941;
t938 = t1057 * t941 - t1061 * t967;
t937 = -t1058 * t938 + t1062 * t939;
t936 = t1058 * t939 + t1062 * t938;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1036, -t1035, 0, -t1059 * t1040 + t1063 * t1041, 0, 0, 0, 0, 0, 0, t1063 * t1015 - t1059 * t1033, t1063 * t1016 + t1059 * t1030, t1063 * t1034 - t1059 * t1037, t1063 * t1000 - t1059 * t1026, 0, 0, 0, 0, 0, 0, t1059 * t989 + t1063 * t972, -t1059 * t1075 + t1063 * t976, t1059 * t1003 + t1063 * t959, -t1059 * t1004 + t1063 * t949, 0, 0, 0, 0, 0, 0, t1059 * t965 + t1063 * t945, t1059 * t969 + t1063 * t947, t1059 * t960 + t1063 * t943, t1059 * t940 + t1063 * t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1035, t1036, 0, t1063 * t1040 + t1059 * t1041, 0, 0, 0, 0, 0, 0, t1059 * t1015 + t1063 * t1033, t1059 * t1016 - t1063 * t1030, t1059 * t1034 + t1063 * t1037, t1059 * t1000 + t1063 * t1026, 0, 0, 0, 0, 0, 0, t1059 * t972 - t1063 * t989, t1059 * t976 + t1063 * t1075, -t1063 * t1003 + t1059 * t959, t1063 * t1004 + t1059 * t949, 0, 0, 0, 0, 0, 0, t1059 * t945 - t1063 * t965, t1059 * t947 - t1063 * t969, t1059 * t943 - t1063 * t960, t1059 * t937 - t1063 * t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1013, t1014, 0, t999, 0, 0, 0, 0, 0, 0, t971, t975, t958, t948, 0, 0, 0, 0, 0, 0, t944, t946, t942, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1065, -qJDD(1), 0, t1041, 0, 0, 0, 0, 0, 0, t1015, t1016, t1034, t1000, 0, 0, 0, 0, 0, 0, t972, t976, t959, t949, 0, 0, 0, 0, 0, 0, t945, t947, t943, t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1065, 0, t1040, 0, 0, 0, 0, 0, 0, t1033, -t1030, t1037, t1026, 0, 0, 0, 0, 0, 0, -t989, t1075, -t1003, t1004, 0, 0, 0, 0, 0, 0, -t965, -t969, -t960, -t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1013, t1014, 0, t999, 0, 0, 0, 0, 0, 0, t971, t975, t958, t948, 0, 0, 0, 0, 0, 0, t944, t946, t942, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1043, t1039, t1051, t1019, 0, 0, 0, 0, 0, 0, t987, t994, t974, t963, 0, 0, 0, 0, 0, 0, t955, t957, t953, t939; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1038, t1042, -t1080, t1018, 0, 0, 0, 0, 0, 0, t986, t993, t973, t962, 0, 0, 0, 0, 0, 0, t954, t956, t952, t938; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1033, t1030, -t1037, -t1026, 0, 0, 0, 0, 0, 0, t989, -t1075, t1003, -t1004, 0, 0, 0, 0, 0, 0, t965, t969, t960, t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1006, t1008, t990, t982, 0, 0, 0, 0, 0, 0, t966, t970, t961, t941; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1007, t1017, t992, t981, 0, 0, 0, 0, 0, 0, -t977, -t979, -t985, -t967; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t989, -t1075, t1003, -t1004, 0, 0, 0, 0, 0, 0, t965, t969, t960, t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t988, t984, t978, t951; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t983, t995, t980, t950; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t977, t979, t985, t967;];
f_new_reg = t1;