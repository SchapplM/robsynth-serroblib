% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:14
% EndTime: 2019-12-31 16:45:16
% DurationCPUTime: 1.37s
% Computational Cost: add. (1999->139), mult. (5029->164), div. (0->0), fcn. (3355->6), ass. (0->94)
t1019 = sin(qJ(1));
t1021 = cos(qJ(1));
t1020 = cos(qJ(3));
t1016 = sin(pkin(6));
t1018 = sin(qJ(3));
t1040 = t1016 * t1018;
t1017 = cos(pkin(6));
t1042 = qJD(1) * t1017;
t994 = qJD(1) * t1040 - t1020 * t1042;
t1034 = -0.2e1 * qJD(3) * t994;
t1029 = t1016 * t1020 + t1017 * t1018;
t993 = t1029 * qJDD(1);
t1028 = t993 + t1034;
t1022 = qJD(3) ^ 2;
t996 = t1029 * qJD(1);
t991 = t996 ^ 2;
t1054 = -t991 - t1022;
t1047 = t996 * t994;
t973 = qJDD(3) + t1047;
t946 = t1018 * t973 - t1020 * t1054;
t948 = t1018 * t1054 + t1020 * t973;
t940 = t1016 * t946 - t1017 * t948;
t1067 = t1019 * t940 - t1021 * t1028;
t1066 = t1019 * t1028 + t1021 * t940;
t978 = t994 ^ 2;
t1055 = -t978 - t1022;
t974 = qJDD(3) - t1047;
t1058 = -t1018 * t974 + t1020 * t1055;
t1059 = t1018 * t1055 + t1020 * t974;
t1061 = -t1016 * t1059 + t1017 * t1058;
t1010 = t1017 * qJDD(1);
t1007 = t1020 * t1010;
t1037 = t1016 * qJDD(1);
t1030 = -t1018 * t1037 + t1007;
t1043 = t996 * qJD(3);
t976 = -t1030 + 0.2e1 * t1043;
t1065 = t1019 * t1061 - t1021 * t976;
t1064 = t1019 * t976 + t1021 * t1061;
t930 = t1016 * t948 + t1017 * t946;
t1050 = t1018 * t993 + t1020 * t1030;
t1051 = t1018 * t1030 - t1020 * t993;
t1057 = -t1016 * t1051 + t1017 * t1050;
t960 = t991 + t978;
t1063 = -t1019 * t960 + t1021 * t1057;
t1062 = t1019 * t1057 + t1021 * t960;
t1060 = t1016 * t1058 + t1017 * t1059;
t1056 = t1016 * t1050 + t1017 * t1051;
t1023 = qJD(1) ^ 2;
t1012 = t1016 ^ 2;
t1013 = t1017 ^ 2;
t1038 = t1012 + t1013;
t1001 = t1038 * t1023;
t1049 = 2 * qJD(4);
t1048 = t1017 * g(3);
t1005 = -t1021 * g(1) - t1019 * g(2);
t997 = -t1023 * pkin(1) + qJDD(1) * qJ(2) + t1005;
t1032 = -0.2e1 * qJD(1) * qJD(2) - t997;
t1039 = t1017 * t1023;
t1026 = -t1048 + (pkin(2) * t1039 - pkin(5) * qJDD(1) + t1032) * t1016;
t1041 = t1013 * t1023;
t980 = -t1016 * g(3) + 0.2e1 * qJD(2) * t1042 + t1017 * t997;
t964 = -pkin(2) * t1041 + pkin(5) * t1010 + t980;
t942 = t1018 * t1026 + t1020 * t964;
t1036 = t1019 * qJDD(1);
t1035 = t1021 * qJDD(1);
t1033 = t1017 * pkin(2) + pkin(1);
t1004 = t1019 * g(1) - t1021 * g(2);
t941 = -t1018 * t964 + t1020 * t1026;
t1031 = -qJDD(2) + t1004;
t1027 = (t1038 * pkin(5) + qJ(2)) * t1023 + t1031;
t1006 = t1016 * t1039;
t1003 = -t1021 * t1023 - t1036;
t1002 = -t1019 * t1023 + t1035;
t1000 = t1038 * qJDD(1);
t999 = t1017 * t1001;
t998 = t1016 * t1001;
t990 = qJDD(1) * pkin(1) + t1023 * qJ(2) + t1031;
t979 = t1032 * t1016 - t1048;
t972 = t1033 * qJDD(1) + t1027;
t969 = t994 * pkin(3) - t996 * qJ(4);
t952 = -t1016 * t979 + t1017 * t980;
t951 = t1016 * t980 + t1017 * t979;
t936 = -qJDD(3) * pkin(3) - t1022 * qJ(4) + t996 * t969 + qJDD(4) - t941;
t935 = -t1022 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t1049 - t994 * t969 + t942;
t928 = qJ(4) * t1034 + t996 * t1049 + (t1029 * qJ(4) + t1033) * qJDD(1) + t1027 + (-t1040 * qJDD(1) + t1007 - 0.2e1 * t1043) * pkin(3);
t927 = -t1018 * t941 + t1020 * t942;
t926 = t1018 * t942 + t1020 * t941;
t925 = t1018 * t936 + t1020 * t935;
t924 = t1018 * t935 - t1020 * t936;
t923 = -t1016 * t926 + t1017 * t927;
t922 = t1016 * t927 + t1017 * t926;
t921 = -t1016 * t924 + t1017 * t925;
t920 = t1016 * t925 + t1017 * t924;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1003, -t1002, 0, -t1019 * t1004 + t1021 * t1005, 0, 0, 0, 0, 0, 0, -t1017 * t1036 - t1021 * t999, t1016 * t1036 + t1021 * t998, t1021 * t1000 - t1019 * t1001, -t1019 * t990 + t1021 * t952, 0, 0, 0, 0, 0, 0, t1064, t1066, t1063, -t1019 * t972 + t1021 * t923, 0, 0, 0, 0, 0, 0, t1064, t1063, -t1066, -t1019 * t928 + t1021 * t921; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1002, t1003, 0, t1021 * t1004 + t1019 * t1005, 0, 0, 0, 0, 0, 0, t1017 * t1035 - t1019 * t999, -t1016 * t1035 + t1019 * t998, t1019 * t1000 + t1021 * t1001, t1019 * t952 + t1021 * t990, 0, 0, 0, 0, 0, 0, t1065, t1067, t1062, t1019 * t923 + t1021 * t972, 0, 0, 0, 0, 0, 0, t1065, t1062, -t1067, t1019 * t921 + t1021 * t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t951, 0, 0, 0, 0, 0, 0, t1060, -t930, t1056, t922, 0, 0, 0, 0, 0, 0, t1060, t1056, t930, t920; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1023, -qJDD(1), 0, t1005, 0, 0, 0, 0, 0, 0, -t999, t998, t1000, t952, 0, 0, 0, 0, 0, 0, t1061, t940, t1057, t923, 0, 0, 0, 0, 0, 0, t1061, t1057, -t940, t921; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1023, 0, t1004, 0, 0, 0, 0, 0, 0, t1010, -t1037, t1001, t990, 0, 0, 0, 0, 0, 0, -t976, -t1028, t960, t972, 0, 0, 0, 0, 0, 0, -t976, t960, t1028, t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t951, 0, 0, 0, 0, 0, 0, t1060, -t930, t1056, t922, 0, 0, 0, 0, 0, 0, t1060, t1056, t930, t920; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1041, t1006, t1010, t980, 0, 0, 0, 0, 0, 0, t1058, -t948, t1050, t927, 0, 0, 0, 0, 0, 0, t1058, t1050, t948, t925; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1006, -t1012 * t1023, -t1037, t979, 0, 0, 0, 0, 0, 0, t1059, -t946, t1051, t926, 0, 0, 0, 0, 0, 0, t1059, t1051, t946, t924; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1010, t1037, -t1001, -t990, 0, 0, 0, 0, 0, 0, t976, t1028, -t960, -t972, 0, 0, 0, 0, 0, 0, t976, -t960, -t1028, -t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1055, -t973, t1030, t942, 0, 0, 0, 0, 0, 0, t1055, t1030, t973, t935; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t974, t1054, -t993, t941, 0, 0, 0, 0, 0, 0, t974, -t993, -t1054, -t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t976, t1028, -t960, -t972, 0, 0, 0, 0, 0, 0, t976, -t960, -t1028, -t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1055, t1030, t973, t935; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t976, -t960, -t1028, -t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t974, t993, t1054, t936;];
f_new_reg = t1;