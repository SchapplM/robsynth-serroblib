% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:19
% EndTime: 2019-12-05 18:00:21
% DurationCPUTime: 1.58s
% Computational Cost: add. (4244->169), mult. (8732->180), div. (0->0), fcn. (5591->6), ass. (0->107)
t1024 = qJD(3) + qJD(4);
t1057 = qJD(4) + t1024;
t1028 = sin(qJ(4));
t1029 = sin(qJ(3));
t1031 = cos(qJ(4));
t1032 = cos(qJ(3));
t1000 = (-t1028 * t1032 - t1029 * t1031) * qJD(1);
t999 = t1000 ^ 2;
t1052 = qJD(1) * t1032;
t1053 = qJD(1) * t1029;
t1002 = -t1028 * t1053 + t1031 * t1052;
t1056 = t1002 ^ 2;
t1055 = t1024 ^ 2;
t1054 = pkin(6) + pkin(1);
t1044 = qJD(3) * t1053;
t1046 = t1032 * qJDD(1);
t1007 = -t1044 + t1046;
t1035 = qJD(1) ^ 2;
t1043 = t1029 * t1035 * t1032;
t1013 = qJDD(3) - t1043;
t1030 = sin(qJ(1));
t1033 = cos(qJ(1));
t1015 = t1030 * g(1) - t1033 * g(2);
t1037 = -t1035 * qJ(2) + qJDD(2) - t1015;
t995 = -t1054 * qJDD(1) + t1037;
t985 = t1029 * g(3) + t1032 * t995;
t976 = (-t1007 - t1044) * pkin(7) + t1013 * pkin(3) + t985;
t1042 = qJD(3) * t1052;
t1047 = t1029 * qJDD(1);
t1006 = -t1042 - t1047;
t1038 = qJD(3) * pkin(3) - pkin(7) * t1052;
t1026 = t1029 ^ 2;
t1050 = t1026 * t1035;
t986 = -t1032 * g(3) + t1029 * t995;
t977 = -pkin(3) * t1050 + t1006 * pkin(7) - qJD(3) * t1038 + t986;
t962 = t1028 * t976 + t1031 * t977;
t1051 = t1002 * t1000;
t1049 = qJD(4) - t1024;
t1027 = t1032 ^ 2;
t1048 = t1026 + t1027;
t1045 = -qJDD(3) - qJDD(4);
t961 = -t1028 * t977 + t1031 * t976;
t1041 = -t1031 * t1006 + t1028 * t1007;
t1016 = -t1033 * g(1) - t1030 * g(2);
t1040 = t1024 * pkin(4) - t1002 * qJ(5);
t983 = -t1045 + t1051;
t1039 = -t1028 * t1006 - t1031 * t1007;
t1036 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t1016;
t968 = -t1049 * t1000 + t1039;
t978 = t1006 * pkin(3) - t1038 * t1052 + (pkin(7) * t1026 + t1054) * t1035 + t1036;
t1034 = qJD(3) ^ 2;
t1018 = -t1027 * t1035 - t1034;
t1017 = -t1034 - t1050;
t1014 = -qJDD(3) - t1043;
t1012 = t1048 * t1035;
t1011 = t1030 * qJDD(1) + t1033 * t1035;
t1010 = t1033 * qJDD(1) - t1030 * t1035;
t1009 = t1048 * qJDD(1);
t1008 = -0.2e1 * t1044 + t1046;
t1005 = 0.2e1 * t1042 + t1047;
t998 = qJDD(1) * pkin(1) - t1037;
t997 = t1035 * pkin(1) + t1036;
t994 = t1054 * t1035 + t1036;
t991 = -t1055 - t1056;
t990 = t1032 * t1014 - t1029 * t1018;
t989 = -t1029 * t1013 + t1032 * t1017;
t988 = t1029 * t1014 + t1032 * t1018;
t987 = t1032 * t1013 + t1029 * t1017;
t984 = t1045 + t1051;
t982 = -t1055 - t999;
t980 = -t999 - t1056;
t979 = -t1002 * qJD(4) - t1041;
t972 = -t1029 * t985 + t1032 * t986;
t971 = t1029 * t986 + t1032 * t985;
t970 = -t1028 * t991 + t1031 * t984;
t969 = t1028 * t984 + t1031 * t991;
t967 = t1057 * t1000 - t1039;
t966 = -t1049 * t1002 - t1041;
t965 = t1057 * t1002 + t1041;
t964 = -t1028 * t983 + t1031 * t982;
t963 = t1028 * t982 + t1031 * t983;
t960 = -t1029 * t969 + t1032 * t970;
t959 = t1029 * t970 + t1032 * t969;
t958 = -t1028 * t968 + t1031 * t966;
t957 = t1028 * t966 + t1031 * t968;
t956 = t979 * pkin(4) + t999 * qJ(5) - t1002 * t1040 - qJDD(5) + t978;
t955 = -t1029 * t963 + t1032 * t964;
t954 = t1029 * t964 + t1032 * t963;
t953 = t1030 * t959 + t1033 * t967;
t952 = t1030 * t967 - t1033 * t959;
t951 = -t999 * pkin(4) + t979 * qJ(5) + 0.2e1 * qJD(5) * t1000 - t1024 * t1040 + t962;
t950 = t1030 * t954 + t1033 * t965;
t949 = t1030 * t965 - t1033 * t954;
t948 = t983 * pkin(4) + t968 * qJ(5) - 0.2e1 * qJD(5) * t1002 + t961;
t947 = -t1028 * t961 + t1031 * t962;
t946 = t1028 * t962 + t1031 * t961;
t945 = -t1029 * t957 + t1032 * t958;
t944 = t1029 * t958 + t1032 * t957;
t943 = t1030 * t944 + t1033 * t980;
t942 = t1030 * t980 - t1033 * t944;
t941 = -t1028 * t948 + t1031 * t951;
t940 = t1028 * t951 + t1031 * t948;
t939 = -t1029 * t946 + t1032 * t947;
t938 = t1029 * t947 + t1032 * t946;
t937 = -t1029 * t940 + t1032 * t941;
t936 = t1029 * t941 + t1032 * t940;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1011, -t1010, 0, -t1030 * t1015 + t1033 * t1016, 0, 0, 0, 0, 0, 0, 0, t1011, t1010, -t1030 * t998 - t1033 * t997, 0, 0, 0, 0, 0, 0, t1033 * t1005 + t1030 * t987, t1033 * t1008 + t1030 * t988, -t1030 * t1009 - t1033 * t1012, t1030 * t971 - t1033 * t994, 0, 0, 0, 0, 0, 0, t950, t953, t943, t1030 * t938 - t1033 * t978, 0, 0, 0, 0, 0, 0, t950, t953, t943, t1030 * t936 - t1033 * t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1010, -t1011, 0, t1033 * t1015 + t1030 * t1016, 0, 0, 0, 0, 0, 0, 0, -t1010, t1011, -t1030 * t997 + t1033 * t998, 0, 0, 0, 0, 0, 0, t1030 * t1005 - t1033 * t987, t1030 * t1008 - t1033 * t988, t1033 * t1009 - t1030 * t1012, -t1030 * t994 - t1033 * t971, 0, 0, 0, 0, 0, 0, t949, t952, t942, -t1030 * t978 - t1033 * t938, 0, 0, 0, 0, 0, 0, t949, t952, t942, -t1030 * t956 - t1033 * t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t989, t990, 0, t972, 0, 0, 0, 0, 0, 0, t955, t960, t945, t939, 0, 0, 0, 0, 0, 0, t955, t960, t945, t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1035, -qJDD(1), 0, t1016, 0, 0, 0, 0, 0, 0, 0, t1035, qJDD(1), -t997, 0, 0, 0, 0, 0, 0, t1005, t1008, -t1012, -t994, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t978, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1035, 0, t1015, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t1035, t998, 0, 0, 0, 0, 0, 0, -t987, -t988, t1009, -t971, 0, 0, 0, 0, 0, 0, -t954, -t959, -t944, -t938, 0, 0, 0, 0, 0, 0, -t954, -t959, -t944, -t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t989, t990, 0, t972, 0, 0, 0, 0, 0, 0, t955, t960, t945, t939, 0, 0, 0, 0, 0, 0, t955, t960, t945, t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t989, t990, 0, t972, 0, 0, 0, 0, 0, 0, t955, t960, t945, t939, 0, 0, 0, 0, 0, 0, t955, t960, t945, t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1035, -qJDD(1), t997, 0, 0, 0, 0, 0, 0, -t1005, -t1008, t1012, t994, 0, 0, 0, 0, 0, 0, -t965, -t967, -t980, t978, 0, 0, 0, 0, 0, 0, -t965, -t967, -t980, t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1035, -t998, 0, 0, 0, 0, 0, 0, t987, t988, -t1009, t971, 0, 0, 0, 0, 0, 0, t954, t959, t944, t938, 0, 0, 0, 0, 0, 0, t954, t959, t944, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1017, t1014, -t1047, t986, 0, 0, 0, 0, 0, 0, t964, t970, t958, t947, 0, 0, 0, 0, 0, 0, t964, t970, t958, t941; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1013, t1018, -t1046, t985, 0, 0, 0, 0, 0, 0, t963, t969, t957, t946, 0, 0, 0, 0, 0, 0, t963, t969, t957, t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1005, t1008, -t1012, -t994, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t978, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t982, t984, t966, t962, 0, 0, 0, 0, 0, 0, t982, t984, t966, t951; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t983, t991, t968, t961, 0, 0, 0, 0, 0, 0, t983, t991, t968, t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t978, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t982, t984, t966, t951; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t983, t991, t968, t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t965, t967, t980, -t956;];
f_new_reg = t1;