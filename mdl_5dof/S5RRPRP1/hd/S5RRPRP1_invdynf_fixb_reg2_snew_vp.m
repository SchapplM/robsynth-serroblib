% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:28
% EndTime: 2020-01-03 11:59:30
% DurationCPUTime: 1.68s
% Computational Cost: add. (4741->139), mult. (6817->181), div. (0->0), fcn. (3888->8), ass. (0->106)
t1047 = sin(qJ(1));
t1050 = cos(qJ(1));
t1039 = qJD(1) + qJD(2);
t1037 = t1039 ^ 2;
t1038 = qJDD(1) + qJDD(2);
t1043 = sin(pkin(8));
t1044 = cos(pkin(8));
t1011 = t1043 * t1037 - t1044 * t1038;
t1046 = sin(qJ(2));
t1049 = cos(qJ(2));
t1058 = -t1044 * t1037 - t1043 * t1038;
t1073 = t1046 * t1011 + t1049 * t1058;
t994 = t1049 * t1011 - t1046 * t1058;
t1077 = t1047 * t994 + t1050 * t1073;
t1076 = t1047 * t1073 - t1050 * t994;
t1017 = t1046 * t1037 - t1049 * t1038;
t1057 = -t1049 * t1037 - t1046 * t1038;
t1072 = t1047 * t1017 + t1050 * t1057;
t1071 = -t1050 * t1017 + t1047 * t1057;
t1029 = -t1050 * g(2) - t1047 * g(3);
t1054 = qJDD(1) * pkin(1) + t1029;
t1028 = -t1047 * g(2) + t1050 * g(3);
t1052 = qJD(1) ^ 2;
t1055 = -t1052 * pkin(1) + t1028;
t999 = -t1046 * t1055 + t1049 * t1054;
t1053 = t1038 * pkin(2) + t999;
t1000 = t1046 * t1054 + t1049 * t1055;
t998 = -t1037 * pkin(2) + t1000;
t982 = t1043 * t1053 + t1044 * t998;
t1042 = -g(1) + qJDD(3);
t1045 = sin(qJ(4));
t1048 = cos(qJ(4));
t978 = -t1037 * pkin(3) + t1038 * pkin(7) + t982;
t975 = t1045 * t1042 + t1048 * t978;
t1066 = t1037 * t1048;
t1065 = t1039 * t1045;
t1041 = t1048 ^ 2;
t1064 = t1041 * t1037;
t1063 = t1045 * t1038;
t1062 = t1048 * t1038;
t1040 = t1045 ^ 2;
t1061 = t1040 + t1041;
t1060 = 0.2e1 * t1039 * t1048;
t1059 = qJD(4) * t1065;
t981 = -t1043 * t998 + t1044 * t1053;
t977 = -t1038 * pkin(3) - t1037 * pkin(7) - t981;
t1056 = -t1059 + t1062;
t1051 = qJD(4) ^ 2;
t1034 = t1048 * t1042;
t1027 = t1045 * t1066;
t1026 = -t1051 - t1064;
t1025 = -t1040 * t1037 - t1051;
t1024 = t1050 * qJDD(1) - t1047 * t1052;
t1023 = -t1047 * qJDD(1) - t1050 * t1052;
t1022 = -qJDD(4) + t1027;
t1021 = qJDD(4) + t1027;
t1020 = qJD(4) * pkin(4) - qJ(5) * t1065;
t1019 = t1061 * t1037;
t1014 = t1061 * t1038;
t1007 = -0.2e1 * t1059 + t1062;
t1006 = qJD(4) * t1060 + t1063;
t1004 = t1048 * t1022 - t1045 * t1025;
t1003 = -t1045 * t1021 + t1048 * t1026;
t1002 = t1045 * t1022 + t1048 * t1025;
t1001 = t1048 * t1021 + t1045 * t1026;
t997 = t1044 * t1014 - t1043 * t1019;
t996 = t1043 * t1014 + t1044 * t1019;
t988 = t1044 * t1004 + t1043 * t1006;
t987 = t1044 * t1003 - t1043 * t1007;
t986 = t1043 * t1004 - t1044 * t1006;
t985 = t1043 * t1003 + t1044 * t1007;
t984 = t1049 * t1000 - t1046 * t999;
t983 = t1046 * t1000 + t1049 * t999;
t980 = -t1046 * t996 + t1049 * t997;
t979 = t1046 * t997 + t1049 * t996;
t974 = -t1045 * t978 + t1034;
t973 = -t1046 * t986 + t1049 * t988;
t972 = -t1046 * t985 + t1049 * t987;
t971 = t1046 * t988 + t1049 * t986;
t970 = t1046 * t987 + t1049 * t985;
t969 = -t1056 * pkin(4) - qJ(5) * t1064 + t1020 * t1065 + qJDD(5) + t977;
t968 = -pkin(4) * t1064 + t1056 * qJ(5) - qJD(4) * t1020 + qJD(5) * t1060 + t975;
t967 = qJDD(4) * pkin(4) + t1034 + (pkin(4) * t1066 - t1038 * qJ(5) - 0.2e1 * qJD(5) * t1039 - t978) * t1045;
t966 = t1047 * t979 - t1050 * t980;
t965 = -t1043 * t981 + t1044 * t982;
t964 = t1047 * t980 + t1050 * t979;
t963 = t1043 * t982 + t1044 * t981;
t962 = -t1045 * t974 + t1048 * t975;
t961 = t1045 * t975 + t1048 * t974;
t960 = t1047 * t971 - t1050 * t973;
t959 = t1047 * t970 - t1050 * t972;
t958 = t1047 * t973 + t1050 * t971;
t957 = t1047 * t972 + t1050 * t970;
t956 = -t1045 * t967 + t1048 * t968;
t955 = t1045 * t968 + t1048 * t967;
t954 = t1043 * t977 + t1044 * t962;
t953 = t1043 * t962 - t1044 * t977;
t952 = -t1046 * t963 + t1049 * t965;
t951 = t1046 * t965 + t1049 * t963;
t950 = t1043 * t969 + t1044 * t956;
t949 = t1043 * t956 - t1044 * t969;
t948 = -t1046 * t953 + t1049 * t954;
t947 = t1046 * t954 + t1049 * t953;
t946 = -t1046 * t949 + t1049 * t950;
t945 = t1046 * t950 + t1049 * t949;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1042, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t961, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t955; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1024, t1023, 0, t1047 * t1028 + t1050 * t1029, 0, 0, 0, 0, 0, 0, t1071, t1072, 0, t1047 * t984 + t1050 * t983, 0, 0, 0, 0, 0, 0, t1076, t1077, 0, t1047 * t952 + t1050 * t951, 0, 0, 0, 0, 0, 0, t957, t958, t964, t1047 * t948 + t1050 * t947, 0, 0, 0, 0, 0, 0, t957, t958, t964, t1047 * t946 + t1050 * t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t1023, t1024, 0, -t1050 * t1028 + t1047 * t1029, 0, 0, 0, 0, 0, 0, -t1072, t1071, 0, t1047 * t983 - t1050 * t984, 0, 0, 0, 0, 0, 0, -t1077, t1076, 0, t1047 * t951 - t1050 * t952, 0, 0, 0, 0, 0, 0, t959, t960, t966, t1047 * t947 - t1050 * t948, 0, 0, 0, 0, 0, 0, t959, t960, t966, t1047 * t945 - t1050 * t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1052, -qJDD(1), 0, t1028, 0, 0, 0, 0, 0, 0, t1057, t1017, 0, t984, 0, 0, 0, 0, 0, 0, t1073, t994, 0, t952, 0, 0, 0, 0, 0, 0, t972, t973, t980, t948, 0, 0, 0, 0, 0, 0, t972, t973, t980, t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1052, 0, t1029, 0, 0, 0, 0, 0, 0, -t1017, t1057, 0, t983, 0, 0, 0, 0, 0, 0, -t994, t1073, 0, t951, 0, 0, 0, 0, 0, 0, t970, t971, t979, t947, 0, 0, 0, 0, 0, 0, t970, t971, t979, t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1042, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t961, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t955; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1037, -t1038, 0, t1000, 0, 0, 0, 0, 0, 0, t1058, t1011, 0, t965, 0, 0, 0, 0, 0, 0, t987, t988, t997, t954, 0, 0, 0, 0, 0, 0, t987, t988, t997, t950; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1038, -t1037, 0, t999, 0, 0, 0, 0, 0, 0, -t1011, t1058, 0, t963, 0, 0, 0, 0, 0, 0, t985, t986, t996, t953, 0, 0, 0, 0, 0, 0, t985, t986, t996, t949; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1042, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t961, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t955; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1037, -t1038, 0, t982, 0, 0, 0, 0, 0, 0, t1003, t1004, t1014, t962, 0, 0, 0, 0, 0, 0, t1003, t1004, t1014, t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1038, -t1037, 0, t981, 0, 0, 0, 0, 0, 0, t1007, -t1006, t1019, -t977, 0, 0, 0, 0, 0, 0, t1007, -t1006, t1019, -t969; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1042, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t961, 0, 0, 0, 0, 0, 0, t1001, t1002, 0, t955; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1026, t1022, t1062, t975, 0, 0, 0, 0, 0, 0, t1026, t1022, t1062, t968; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1021, t1025, -t1063, t974, 0, 0, 0, 0, 0, 0, t1021, t1025, -t1063, t967; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1007, t1006, -t1019, t977, 0, 0, 0, 0, 0, 0, -t1007, t1006, -t1019, t969; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1026, t1022, t1062, t968; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1021, t1025, -t1063, t967; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1007, t1006, -t1019, t969;];
f_new_reg = t1;