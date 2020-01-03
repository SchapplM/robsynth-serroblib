% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR10_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:08
% EndTime: 2019-12-31 17:12:09
% DurationCPUTime: 1.07s
% Computational Cost: add. (1918->165), mult. (4138->171), div. (0->0), fcn. (2331->6), ass. (0->109)
t912 = sin(qJ(2));
t915 = cos(qJ(2));
t918 = qJD(1) ^ 2;
t943 = t915 * t918;
t896 = t912 * t943;
t888 = qJDD(2) + t896;
t917 = qJD(2) ^ 2;
t908 = t915 ^ 2;
t946 = t908 * t918;
t894 = -t917 - t946;
t868 = t912 * t888 - t915 * t894;
t935 = qJD(1) * qJD(2);
t902 = t912 * t935;
t933 = t915 * qJDD(1);
t883 = -0.2e1 * t902 + t933;
t913 = sin(qJ(1));
t916 = cos(qJ(1));
t958 = t913 * t868 - t916 * t883;
t957 = t916 * t868 + t913 * t883;
t889 = qJDD(2) - t896;
t907 = t912 ^ 2;
t893 = -t907 * t918 - t917;
t869 = t915 * t889 + t912 * t893;
t931 = t915 * t935;
t934 = t912 * qJDD(1);
t880 = 0.2e1 * t931 + t934;
t956 = t913 * t869 + t916 * t880;
t955 = t916 * t869 - t913 * t880;
t881 = t931 + t934;
t891 = t913 * g(1) - t916 * g(2);
t924 = -qJDD(1) * pkin(1) - t891;
t938 = t912 * qJD(1);
t950 = 2 * qJD(3);
t954 = pkin(2) * t902 - (t881 + t931) * qJ(3) - t938 * t950 + t924;
t911 = sin(qJ(4));
t914 = cos(qJ(4));
t939 = qJD(1) * t915;
t876 = t911 * qJD(2) + t914 * t939;
t953 = t876 ^ 2;
t878 = t914 * qJD(2) - t911 * t939;
t952 = t878 ^ 2;
t898 = qJD(4) + t938;
t951 = t898 ^ 2;
t949 = t915 * g(3);
t948 = t918 * pkin(5);
t947 = t878 * t876;
t940 = t907 + t908;
t937 = qJD(4) - t898;
t936 = qJD(4) + t898;
t932 = -t951 - t952;
t892 = -t916 * g(1) - t913 * g(2);
t874 = -t918 * pkin(1) + qJDD(1) * pkin(5) + t892;
t930 = t918 * (-pkin(2) * t915 - qJ(3) * t912) + t874;
t882 = -t902 + t933;
t929 = -t911 * qJDD(2) - t914 * t882;
t890 = pkin(3) * t938 - qJD(2) * pkin(6);
t840 = -t890 * t938 + (-pkin(3) * t908 - pkin(5)) * t918 + (-pkin(2) - pkin(6)) * t882 + t954;
t925 = -qJDD(2) * pkin(2) - t917 * qJ(3) + qJDD(3) + t949;
t844 = -qJDD(2) * pkin(6) + (t881 - t931) * pkin(3) + (-pkin(6) * t943 + t930) * t912 + t925;
t927 = -t911 * t840 + t914 * t844;
t863 = t915 * t888 + t912 * t894;
t865 = t912 * t889 - t915 * t893;
t926 = -t914 * qJDD(2) + t911 * t882;
t923 = -qJDD(4) - t881;
t922 = t937 * t876 + t926;
t921 = -t923 - t947;
t904 = t912 * g(3);
t920 = -t917 * pkin(2) + qJDD(2) * qJ(3) + t930 * t915 - t904;
t887 = t940 * t918;
t886 = -t913 * qJDD(1) - t916 * t918;
t885 = t916 * qJDD(1) - t913 * t918;
t884 = t940 * qJDD(1);
t873 = -t924 + t948;
t872 = t915 * t874 - t904;
t871 = -t912 * t874 - t949;
t862 = t916 * t884 - t913 * t887;
t861 = t913 * t884 + t916 * t887;
t859 = -t951 - t953;
t858 = t923 - t947;
t857 = -t952 - t953;
t855 = t930 * t912 + t925;
t854 = qJD(2) * t950 + t920;
t853 = -t912 * t871 + t915 * t872;
t852 = t915 * t871 + t912 * t872;
t851 = -t936 * t876 - t926;
t850 = -t937 * t878 + t929;
t849 = t936 * t878 - t929;
t847 = t882 * pkin(2) + t948 - t954;
t846 = t914 * t858 - t911 * t932;
t845 = t911 * t858 + t914 * t932;
t843 = -pkin(6) * t946 + t882 * pkin(3) + (t950 + t890) * qJD(2) + t920;
t842 = t914 * t859 - t911 * t921;
t841 = t911 * t859 + t914 * t921;
t839 = t915 * t854 + t912 * t855;
t838 = t912 * t854 - t915 * t855;
t837 = t914 * t850 - t911 * t922;
t836 = t911 * t850 + t914 * t922;
t835 = t912 * t845 + t915 * t851;
t834 = -t915 * t845 + t912 * t851;
t833 = t912 * t841 + t915 * t849;
t832 = -t915 * t841 + t912 * t849;
t831 = t912 * t836 + t915 * t857;
t830 = -t915 * t836 + t912 * t857;
t829 = t914 * t840 + t911 * t844;
t827 = t914 * t829 - t911 * t927;
t826 = t911 * t829 + t914 * t927;
t825 = t912 * t826 + t915 * t843;
t824 = -t915 * t826 + t912 * t843;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t886, -t885, 0, -t913 * t891 + t916 * t892, 0, 0, 0, 0, 0, 0, -t957, -t955, t862, t916 * t853 - t913 * t873, 0, 0, 0, 0, 0, 0, t862, t957, t955, t916 * t839 - t913 * t847, 0, 0, 0, 0, 0, 0, t916 * t833 + t913 * t842, t916 * t835 + t913 * t846, t916 * t831 + t913 * t837, t916 * t825 + t913 * t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t885, t886, 0, t916 * t891 + t913 * t892, 0, 0, 0, 0, 0, 0, -t958, -t956, t861, t913 * t853 + t916 * t873, 0, 0, 0, 0, 0, 0, t861, t958, t956, t913 * t839 + t916 * t847, 0, 0, 0, 0, 0, 0, t913 * t833 - t916 * t842, t913 * t835 - t916 * t846, t913 * t831 - t916 * t837, t913 * t825 - t916 * t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t863, -t865, 0, t852, 0, 0, 0, 0, 0, 0, 0, -t863, t865, t838, 0, 0, 0, 0, 0, 0, t832, t834, t830, t824; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t918, -qJDD(1), 0, t892, 0, 0, 0, 0, 0, 0, -t868, -t869, t884, t853, 0, 0, 0, 0, 0, 0, t884, t868, t869, t839, 0, 0, 0, 0, 0, 0, t833, t835, t831, t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t918, 0, t891, 0, 0, 0, 0, 0, 0, t883, -t880, t887, t873, 0, 0, 0, 0, 0, 0, t887, -t883, t880, t847, 0, 0, 0, 0, 0, 0, -t842, -t846, -t837, -t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t863, -t865, 0, t852, 0, 0, 0, 0, 0, 0, 0, -t863, t865, t838, 0, 0, 0, 0, 0, 0, t832, t834, t830, t824; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t894, -t889, t933, t872, 0, 0, 0, 0, 0, 0, t933, -t894, t889, t854, 0, 0, 0, 0, 0, 0, t849, t851, t857, t843; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t888, t893, -t934, t871, 0, 0, 0, 0, 0, 0, -t934, -t888, -t893, -t855, 0, 0, 0, 0, 0, 0, -t841, -t845, -t836, -t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t883, t880, -t887, -t873, 0, 0, 0, 0, 0, 0, -t887, t883, -t880, -t847, 0, 0, 0, 0, 0, 0, t842, t846, t837, t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t887, t883, -t880, -t847, 0, 0, 0, 0, 0, 0, t842, t846, t837, t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t933, t894, -t889, -t854, 0, 0, 0, 0, 0, 0, -t849, -t851, -t857, -t843; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, t888, t893, t855, 0, 0, 0, 0, 0, 0, t841, t845, t836, t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t859, t858, t850, t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t921, t932, t922, t927; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t851, t857, t843;];
f_new_reg = t1;
