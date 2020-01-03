% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRP4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:53
% EndTime: 2019-12-31 17:15:54
% DurationCPUTime: 1.00s
% Computational Cost: add. (3011->135), mult. (6602->176), div. (0->0), fcn. (4357->6), ass. (0->106)
t890 = qJD(2) + qJD(3);
t921 = qJD(3) + t890;
t893 = sin(qJ(3));
t896 = cos(qJ(3));
t897 = cos(qJ(2));
t894 = sin(qJ(2));
t912 = qJD(1) * t894;
t864 = -t896 * t897 * qJD(1) + t893 * t912;
t863 = t864 ^ 2;
t866 = (t893 * t897 + t894 * t896) * qJD(1);
t920 = t866 ^ 2;
t919 = t890 ^ 2;
t918 = -2 * qJD(4);
t917 = t866 * t864;
t892 = t897 ^ 2;
t900 = qJD(1) ^ 2;
t916 = t892 * t900;
t895 = sin(qJ(1));
t898 = cos(qJ(1));
t882 = -g(1) * t898 - g(2) * t895;
t868 = -pkin(1) * t900 + qJDD(1) * pkin(5) + t882;
t915 = t894 * t868;
t914 = t894 * t900;
t910 = qJD(1) * qJD(2);
t906 = t897 * t910;
t909 = t894 * qJDD(1);
t872 = t906 + t909;
t843 = qJDD(2) * pkin(2) - t872 * pkin(6) - t915 + (pkin(2) * t914 + pkin(6) * t910 - g(3)) * t897;
t860 = -t894 * g(3) + t897 * t868;
t888 = t897 * qJDD(1);
t907 = t894 * t910;
t873 = t888 - t907;
t902 = qJD(2) * pkin(2) - pkin(6) * t912;
t846 = -pkin(2) * t916 + t873 * pkin(6) - qJD(2) * t902 + t860;
t831 = t893 * t843 + t896 * t846;
t891 = t894 ^ 2;
t913 = t891 + t892;
t911 = qJD(3) - t890;
t908 = -qJDD(2) - qJDD(3);
t881 = t895 * g(1) - t898 * g(2);
t830 = t896 * t843 - t893 * t846;
t905 = t893 * t872 - t896 * t873;
t904 = pkin(3) * t890 - qJ(4) * t866;
t903 = -t896 * t872 - t893 * t873;
t852 = -t908 - t917;
t901 = qJDD(1) * pkin(1) + t881;
t837 = t864 * t911 + t903;
t849 = t873 * pkin(2) - t902 * t912 + (pkin(6) * t892 + pkin(5)) * t900 + t901;
t899 = qJD(2) ^ 2;
t886 = t897 * t914;
t884 = -t899 - t916;
t883 = -t891 * t900 - t899;
t880 = -qJDD(2) + t886;
t879 = qJDD(2) + t886;
t878 = t913 * t900;
t877 = -qJDD(1) * t895 - t898 * t900;
t876 = qJDD(1) * t898 - t895 * t900;
t875 = t913 * qJDD(1);
t874 = t888 - 0.2e1 * t907;
t871 = 0.2e1 * t906 + t909;
t867 = pkin(5) * t900 + t901;
t859 = -t897 * g(3) - t915;
t858 = -t919 - t920;
t857 = t880 * t897 - t883 * t894;
t856 = -t879 * t894 + t884 * t897;
t855 = t880 * t894 + t883 * t897;
t854 = t879 * t897 + t884 * t894;
t853 = t908 - t917;
t851 = -t919 - t863;
t848 = -t863 - t920;
t847 = -qJD(3) * t866 - t905;
t845 = -t859 * t894 + t860 * t897;
t844 = t859 * t897 + t860 * t894;
t839 = t853 * t896 - t858 * t893;
t838 = t853 * t893 + t858 * t896;
t836 = -t864 * t921 - t903;
t835 = -t866 * t911 - t905;
t834 = t866 * t921 + t905;
t833 = t851 * t896 - t852 * t893;
t832 = t851 * t893 + t852 * t896;
t829 = -t838 * t894 + t839 * t897;
t828 = t838 * t897 + t839 * t894;
t827 = t835 * t896 - t837 * t893;
t826 = t835 * t893 + t837 * t896;
t825 = t847 * pkin(3) + t863 * qJ(4) - t866 * t904 - qJDD(4) + t849;
t824 = -t832 * t894 + t833 * t897;
t823 = t832 * t897 + t833 * t894;
t822 = t829 * t898 + t836 * t895;
t821 = t829 * t895 - t836 * t898;
t820 = -t863 * pkin(3) + t847 * qJ(4) + t864 * t918 - t890 * t904 + t831;
t819 = -t830 * t893 + t831 * t896;
t818 = t830 * t896 + t831 * t893;
t817 = t824 * t898 + t834 * t895;
t816 = t824 * t895 - t834 * t898;
t815 = pkin(3) * t852 + qJ(4) * t837 + t866 * t918 + t830;
t814 = -t826 * t894 + t827 * t897;
t813 = t826 * t897 + t827 * t894;
t812 = t814 * t898 + t848 * t895;
t811 = t814 * t895 - t848 * t898;
t810 = -t815 * t893 + t820 * t896;
t809 = t815 * t896 + t820 * t893;
t808 = -t818 * t894 + t819 * t897;
t807 = t818 * t897 + t819 * t894;
t806 = -t809 * t894 + t810 * t897;
t805 = t809 * t897 + t810 * t894;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t877, -t876, 0, -t881 * t895 + t882 * t898, 0, 0, 0, 0, 0, 0, t856 * t898 - t874 * t895, t857 * t898 + t871 * t895, t875 * t898 - t878 * t895, t845 * t898 - t867 * t895, 0, 0, 0, 0, 0, 0, t817, t822, t812, t808 * t898 - t849 * t895, 0, 0, 0, 0, 0, 0, t817, t822, t812, t806 * t898 - t825 * t895; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t876, t877, 0, t881 * t898 + t882 * t895, 0, 0, 0, 0, 0, 0, t856 * t895 + t874 * t898, t857 * t895 - t871 * t898, t875 * t895 + t878 * t898, t845 * t895 + t867 * t898, 0, 0, 0, 0, 0, 0, t816, t821, t811, t808 * t895 + t849 * t898, 0, 0, 0, 0, 0, 0, t816, t821, t811, t806 * t895 + t825 * t898; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t854, t855, 0, t844, 0, 0, 0, 0, 0, 0, t823, t828, t813, t807, 0, 0, 0, 0, 0, 0, t823, t828, t813, t805; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t900, -qJDD(1), 0, t882, 0, 0, 0, 0, 0, 0, t856, t857, t875, t845, 0, 0, 0, 0, 0, 0, t824, t829, t814, t808, 0, 0, 0, 0, 0, 0, t824, t829, t814, t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t900, 0, t881, 0, 0, 0, 0, 0, 0, t874, -t871, t878, t867, 0, 0, 0, 0, 0, 0, -t834, -t836, -t848, t849, 0, 0, 0, 0, 0, 0, -t834, -t836, -t848, t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t854, t855, 0, t844, 0, 0, 0, 0, 0, 0, t823, t828, t813, t807, 0, 0, 0, 0, 0, 0, t823, t828, t813, t805; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t884, t880, t888, t860, 0, 0, 0, 0, 0, 0, t833, t839, t827, t819, 0, 0, 0, 0, 0, 0, t833, t839, t827, t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t879, t883, -t909, t859, 0, 0, 0, 0, 0, 0, t832, t838, t826, t818, 0, 0, 0, 0, 0, 0, t832, t838, t826, t809; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t874, t871, -t878, -t867, 0, 0, 0, 0, 0, 0, t834, t836, t848, -t849, 0, 0, 0, 0, 0, 0, t834, t836, t848, -t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, t853, t835, t831, 0, 0, 0, 0, 0, 0, t851, t853, t835, t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t852, t858, t837, t830, 0, 0, 0, 0, 0, 0, t852, t858, t837, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t834, t836, t848, -t849, 0, 0, 0, 0, 0, 0, t834, t836, t848, -t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, t853, t835, t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t852, t858, t837, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t834, t836, t848, -t825;];
f_new_reg = t1;
