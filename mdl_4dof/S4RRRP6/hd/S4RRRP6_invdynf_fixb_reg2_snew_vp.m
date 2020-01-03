% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRP6
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRP6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:23
% EndTime: 2019-12-31 17:19:24
% DurationCPUTime: 0.98s
% Computational Cost: add. (2697->137), mult. (5537->172), div. (0->0), fcn. (3584->6), ass. (0->103)
t910 = cos(qJ(2));
t921 = t910 * qJD(1);
t894 = -qJD(3) + t921;
t931 = qJD(3) - t894;
t906 = sin(qJ(3));
t909 = cos(qJ(3));
t907 = sin(qJ(2));
t922 = qJD(1) * t907;
t876 = -t909 * qJD(2) + t906 * t922;
t930 = t876 ^ 2;
t878 = t906 * qJD(2) + t909 * t922;
t929 = t878 ^ 2;
t928 = t894 ^ 2;
t927 = qJD(2) ^ 2;
t926 = -2 * qJD(4);
t925 = t910 * g(3);
t924 = t878 * t876;
t908 = sin(qJ(1));
t911 = cos(qJ(1));
t889 = t908 * g(1) - t911 * g(2);
t912 = qJD(1) ^ 2;
t873 = qJDD(1) * pkin(1) + t912 * pkin(5) + t889;
t919 = qJD(1) * qJD(2);
t917 = t910 * t919;
t918 = t907 * qJDD(1);
t881 = t917 + t918;
t897 = t907 * t919;
t899 = t910 * qJDD(1);
t882 = t899 - 0.2e1 * t897;
t848 = (-t881 - t917) * pkin(6) - t882 * pkin(2) - t873;
t890 = -t911 * g(1) - t908 * g(2);
t874 = -t912 * pkin(1) + qJDD(1) * pkin(5) + t890;
t870 = -t907 * g(3) + t910 * t874;
t879 = (-pkin(2) * t910 - pkin(6) * t907) * qJD(1);
t857 = -pkin(2) * t927 + qJDD(2) * pkin(6) + t879 * t921 + t870;
t841 = t906 * t848 + t909 * t857;
t902 = t907 ^ 2;
t903 = t910 ^ 2;
t923 = t902 + t903;
t920 = qJD(3) + t894;
t916 = t899 - qJDD(3) - t897;
t840 = t909 * t848 - t906 * t857;
t915 = -t909 * qJDD(2) + t906 * t881;
t914 = -t906 * qJDD(2) - t909 * t881;
t860 = -t916 - t924;
t913 = -t878 * qJD(3) - t915;
t852 = t876 * t920 + t914;
t856 = t925 - qJDD(2) * pkin(2) - t927 * pkin(6) + (qJD(1) * t879 + t874) * t907;
t893 = t910 * t912 * t907;
t892 = -t903 * t912 - t927;
t891 = -t902 * t912 - t927;
t888 = -qJDD(2) + t893;
t887 = qJDD(2) + t893;
t886 = t923 * t912;
t885 = -t908 * qJDD(1) - t911 * t912;
t884 = t911 * qJDD(1) - t908 * t912;
t883 = t923 * qJDD(1);
t880 = 0.2e1 * t917 + t918;
t869 = -t907 * t874 - t925;
t868 = -t894 * pkin(3) - t878 * qJ(4);
t867 = t910 * t888 - t907 * t891;
t866 = -t907 * t887 + t910 * t892;
t865 = t907 * t888 + t910 * t891;
t864 = t910 * t887 + t907 * t892;
t863 = -t928 - t929;
t862 = -t928 - t930;
t859 = t916 - t924;
t858 = -t929 - t930;
t854 = -t907 * t869 + t910 * t870;
t853 = t910 * t869 + t907 * t870;
t851 = -t931 * t876 - t914;
t850 = -t878 * t920 - t915;
t849 = t931 * t878 + t915;
t845 = t909 * t859 - t906 * t863;
t844 = t906 * t859 + t909 * t863;
t843 = -t906 * t860 + t909 * t862;
t842 = t909 * t860 + t906 * t862;
t839 = t909 * t850 - t906 * t852;
t838 = t906 * t850 + t909 * t852;
t837 = -pkin(3) * t913 - qJ(4) * t930 + t878 * t868 + qJDD(4) + t856;
t836 = t910 * t845 + t907 * t851;
t835 = t907 * t845 - t910 * t851;
t834 = t910 * t843 + t907 * t849;
t833 = t907 * t843 - t910 * t849;
t832 = t910 * t839 + t907 * t858;
t831 = t907 * t839 - t910 * t858;
t830 = t894 * t868 + t913 * qJ(4) + (-pkin(3) * t876 + t926) * t876 + t841;
t829 = t911 * t836 + t908 * t844;
t828 = t908 * t836 - t911 * t844;
t827 = pkin(3) * t860 + qJ(4) * t852 + t878 * t926 + t840;
t826 = -t906 * t840 + t909 * t841;
t825 = t909 * t840 + t906 * t841;
t824 = t911 * t834 + t908 * t842;
t823 = t908 * t834 - t911 * t842;
t822 = t910 * t826 + t907 * t856;
t821 = t907 * t826 - t910 * t856;
t820 = t911 * t832 + t908 * t838;
t819 = t908 * t832 - t911 * t838;
t818 = -t906 * t827 + t909 * t830;
t817 = t909 * t827 + t906 * t830;
t816 = t910 * t818 + t907 * t837;
t815 = t907 * t818 - t910 * t837;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t885, -t884, 0, -t908 * t889 + t911 * t890, 0, 0, 0, 0, 0, 0, t911 * t866 - t908 * t882, t911 * t867 + t908 * t880, t911 * t883 - t908 * t886, t911 * t854 - t908 * t873, 0, 0, 0, 0, 0, 0, t824, t829, t820, t911 * t822 + t908 * t825, 0, 0, 0, 0, 0, 0, t824, t829, t820, t911 * t816 + t908 * t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t884, t885, 0, t911 * t889 + t908 * t890, 0, 0, 0, 0, 0, 0, t908 * t866 + t911 * t882, t908 * t867 - t911 * t880, t908 * t883 + t911 * t886, t908 * t854 + t911 * t873, 0, 0, 0, 0, 0, 0, t823, t828, t819, t908 * t822 - t911 * t825, 0, 0, 0, 0, 0, 0, t823, t828, t819, t908 * t816 - t911 * t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t864, t865, 0, t853, 0, 0, 0, 0, 0, 0, t833, t835, t831, t821, 0, 0, 0, 0, 0, 0, t833, t835, t831, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t912, -qJDD(1), 0, t890, 0, 0, 0, 0, 0, 0, t866, t867, t883, t854, 0, 0, 0, 0, 0, 0, t834, t836, t832, t822, 0, 0, 0, 0, 0, 0, t834, t836, t832, t816; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t912, 0, t889, 0, 0, 0, 0, 0, 0, t882, -t880, t886, t873, 0, 0, 0, 0, 0, 0, -t842, -t844, -t838, -t825, 0, 0, 0, 0, 0, 0, -t842, -t844, -t838, -t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t864, t865, 0, t853, 0, 0, 0, 0, 0, 0, t833, t835, t831, t821, 0, 0, 0, 0, 0, 0, t833, t835, t831, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t892, t888, t899, t870, 0, 0, 0, 0, 0, 0, t843, t845, t839, t826, 0, 0, 0, 0, 0, 0, t843, t845, t839, t818; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t887, t891, -t918, t869, 0, 0, 0, 0, 0, 0, -t849, -t851, -t858, -t856, 0, 0, 0, 0, 0, 0, -t849, -t851, -t858, -t837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t882, t880, -t886, -t873, 0, 0, 0, 0, 0, 0, t842, t844, t838, t825, 0, 0, 0, 0, 0, 0, t842, t844, t838, t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t862, t859, t850, t841, 0, 0, 0, 0, 0, 0, t862, t859, t850, t830; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t860, t863, t852, t840, 0, 0, 0, 0, 0, 0, t860, t863, t852, t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t851, t858, t856, 0, 0, 0, 0, 0, 0, t849, t851, t858, t837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t862, t859, t850, t830; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t860, t863, t852, t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t851, t858, t837;];
f_new_reg = t1;
