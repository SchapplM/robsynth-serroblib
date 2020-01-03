% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5PPRRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:36
% EndTime: 2019-12-05 15:11:37
% DurationCPUTime: 1.66s
% Computational Cost: add. (2177->149), mult. (4094->165), div. (0->0), fcn. (2797->8), ass. (0->94)
t972 = sin(qJ(4));
t974 = cos(qJ(4));
t992 = qJD(3) ^ 2;
t955 = t972 * t992 * t974;
t951 = qJDD(4) - t955;
t963 = t972 ^ 2;
t976 = qJD(4) ^ 2;
t952 = t963 * t992 + t976;
t927 = t951 * t974 - t952 * t972;
t980 = t972 * qJDD(3);
t981 = qJD(3) * qJD(4);
t941 = 0.2e1 * t974 * t981 + t980;
t973 = sin(qJ(3));
t975 = cos(qJ(3));
t914 = t927 * t975 - t973 * t941;
t924 = t951 * t972 + t952 * t974;
t967 = sin(pkin(8));
t969 = cos(pkin(8));
t899 = t914 * t969 + t967 * t924;
t911 = t927 * t973 + t975 * t941;
t968 = sin(pkin(7));
t970 = cos(pkin(7));
t998 = t899 * t968 - t970 * t911;
t997 = t899 * t970 + t968 * t911;
t896 = t914 * t967 - t969 * t924;
t949 = -g(1) * t970 - g(2) * t968;
t984 = -g(3) + qJDD(1);
t932 = t949 * t967 - t969 * t984;
t991 = t967 * t932;
t964 = t974 ^ 2;
t983 = t963 + t964;
t944 = t983 * qJDD(3);
t947 = t983 * t992;
t921 = t944 * t975 - t947 * t973;
t919 = t969 * t921;
t945 = qJDD(3) * t973 + t975 * t992;
t990 = t969 * t945;
t946 = qJDD(3) * t975 - t973 * t992;
t989 = t969 * t946;
t988 = t970 * t945;
t987 = t970 * t946;
t933 = t949 * t969 + t967 * t984;
t948 = g(1) * t968 - g(2) * t970;
t943 = -qJDD(2) + t948;
t917 = t933 * t975 - t943 * t973;
t909 = -pkin(3) * t992 + qJDD(3) * pkin(6) + t917;
t902 = t909 * t974 + t932 * t972;
t982 = t992 * (-pkin(4) * t974 - qJ(5) * t972);
t979 = t974 * qJDD(3);
t978 = t972 * t981;
t916 = -t973 * t933 - t943 * t975;
t908 = -qJDD(3) * pkin(3) - pkin(6) * t992 - t916;
t953 = -t964 * t992 - t976;
t950 = qJDD(4) + t955;
t942 = -0.2e1 * t978 + t979;
t937 = t967 * t945;
t936 = t967 * t946;
t930 = t974 * t932;
t926 = -t950 * t972 + t953 * t974;
t923 = t950 * t974 + t953 * t972;
t922 = t969 * t932;
t920 = t944 * t973 + t947 * t975;
t918 = t967 * t921;
t913 = t926 * t975 - t942 * t973;
t910 = t926 * t973 + t942 * t975;
t906 = t933 * t969 + t991;
t905 = t933 * t967 - t922;
t904 = t919 * t970 + t920 * t968;
t903 = t919 * t968 - t920 * t970;
t901 = -t909 * t972 + t930;
t898 = t913 * t969 + t923 * t967;
t895 = t913 * t967 - t923 * t969;
t894 = -t916 * t973 + t917 * t975;
t893 = t916 * t975 + t917 * t973;
t892 = -(-t978 + t979) * pkin(4) + (pkin(4) * qJD(4) - (2 * qJD(5))) * t972 * qJD(3) + t908 - t941 * qJ(5);
t891 = qJDD(5) - t930 - t976 * qJ(5) - qJDD(4) * pkin(4) + (t909 + t982) * t972;
t890 = -pkin(4) * t976 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t974 * t982 + t902;
t889 = t894 * t969 + t991;
t888 = t894 * t967 - t922;
t887 = t898 * t970 + t910 * t968;
t886 = t898 * t968 - t910 * t970;
t885 = -t901 * t972 + t902 * t974;
t884 = t901 * t974 + t902 * t972;
t883 = t890 * t974 + t891 * t972;
t882 = t890 * t972 - t891 * t974;
t881 = t885 * t975 + t908 * t973;
t880 = t885 * t973 - t908 * t975;
t879 = t883 * t975 + t892 * t973;
t878 = t883 * t973 - t892 * t975;
t877 = t881 * t969 + t884 * t967;
t876 = t881 * t967 - t884 * t969;
t875 = t879 * t969 + t882 * t967;
t874 = t879 * t967 - t882 * t969;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t948 * t968 + t949 * t970, 0, 0, 0, 0, 0, 0, 0, 0, 0, t906 * t970 - t943 * t968, 0, 0, 0, 0, 0, 0, t946 * t968 - t969 * t988, -t945 * t968 - t969 * t987, 0, t889 * t970 + t893 * t968, 0, 0, 0, 0, 0, 0, t887, -t997, t904, t877 * t970 + t880 * t968, 0, 0, 0, 0, 0, 0, t887, t904, t997, t875 * t970 + t878 * t968; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t948 * t970 + t949 * t968, 0, 0, 0, 0, 0, 0, 0, 0, 0, t906 * t968 + t943 * t970, 0, 0, 0, 0, 0, 0, -t968 * t990 - t987, -t968 * t989 + t988, 0, t889 * t968 - t893 * t970, 0, 0, 0, 0, 0, 0, t886, -t998, t903, t877 * t968 - t880 * t970, 0, 0, 0, 0, 0, 0, t886, t903, t998, t875 * t968 - t878 * t970; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t984, 0, 0, 0, 0, 0, 0, 0, 0, 0, t905, 0, 0, 0, 0, 0, 0, -t937, -t936, 0, t888, 0, 0, 0, 0, 0, 0, t895, -t896, t918, t876, 0, 0, 0, 0, 0, 0, t895, t918, t896, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t949, 0, 0, 0, 0, 0, 0, 0, 0, 0, t906, 0, 0, 0, 0, 0, 0, -t990, -t989, 0, t889, 0, 0, 0, 0, 0, 0, t898, -t899, t919, t877, 0, 0, 0, 0, 0, 0, t898, t919, t899, t875; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t948, 0, 0, 0, 0, 0, 0, 0, 0, 0, t943, 0, 0, 0, 0, 0, 0, -t946, t945, 0, -t893, 0, 0, 0, 0, 0, 0, -t910, t911, -t920, -t880, 0, 0, 0, 0, 0, 0, -t910, -t920, -t911, -t878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t984, 0, 0, 0, 0, 0, 0, 0, 0, 0, t905, 0, 0, 0, 0, 0, 0, -t937, -t936, 0, t888, 0, 0, 0, 0, 0, 0, t895, -t896, t918, t876, 0, 0, 0, 0, 0, 0, t895, t918, t896, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, 0, 0, 0, 0, 0, 0, -t945, -t946, 0, t894, 0, 0, 0, 0, 0, 0, t913, -t914, t921, t881, 0, 0, 0, 0, 0, 0, t913, t921, t914, t879; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t932, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t932, 0, 0, 0, 0, 0, 0, -t923, t924, 0, -t884, 0, 0, 0, 0, 0, 0, -t923, 0, -t924, -t882; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t943, 0, 0, 0, 0, 0, 0, t946, -t945, 0, t893, 0, 0, 0, 0, 0, 0, t910, -t911, t920, t880, 0, 0, 0, 0, 0, 0, t910, t920, t911, t878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t992, -qJDD(3), 0, t917, 0, 0, 0, 0, 0, 0, t926, -t927, t944, t885, 0, 0, 0, 0, 0, 0, t926, t944, t927, t883; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t992, 0, t916, 0, 0, 0, 0, 0, 0, t942, -t941, t947, -t908, 0, 0, 0, 0, 0, 0, t942, t947, t941, -t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t932, 0, 0, 0, 0, 0, 0, t923, -t924, 0, t884, 0, 0, 0, 0, 0, 0, t923, 0, t924, t882; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t953, -t951, t979, t902, 0, 0, 0, 0, 0, 0, t953, t979, t951, t890; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t950, -t952, -t980, t901, 0, 0, 0, 0, 0, 0, t950, -t980, t952, -t891; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t942, t941, -t947, t908, 0, 0, 0, 0, 0, 0, -t942, -t947, -t941, t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t953, t979, t951, t890; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t942, -t947, -t941, t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t950, t980, -t952, t891;];
f_new_reg = t1;