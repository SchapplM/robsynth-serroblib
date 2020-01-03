% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRPP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:54
% EndTime: 2019-12-31 18:16:55
% DurationCPUTime: 1.35s
% Computational Cost: add. (1190->173), mult. (2548->124), div. (0->0), fcn. (1134->4), ass. (0->78)
t843 = cos(qJ(3));
t863 = qJD(1) * qJD(3);
t855 = t843 * t863;
t841 = sin(qJ(3));
t862 = t841 * qJDD(1);
t877 = t855 + t862;
t846 = qJD(1) ^ 2;
t869 = t841 * t846;
t858 = t843 * t869;
t820 = qJDD(3) + t858;
t840 = t843 ^ 2;
t845 = qJD(3) ^ 2;
t826 = t840 * t846 + t845;
t799 = t820 * t841 + t826 * t843;
t856 = t841 * t863;
t861 = t843 * qJDD(1);
t814 = -0.2e1 * t856 + t861;
t842 = sin(qJ(1));
t844 = cos(qJ(1));
t793 = t799 * t842 - t844 * t814;
t876 = t799 * t844 + t842 * t814;
t875 = t841 ^ 2;
t874 = -2 * qJD(2);
t873 = 2 * qJD(4);
t872 = 2 * qJD(5);
t871 = -pkin(6) - pkin(1);
t870 = t841 * g(3);
t868 = t875 * t846;
t822 = t842 * g(1) - t844 * g(2);
t851 = -t846 * qJ(2) + qJDD(2) - t822;
t806 = qJDD(1) * t871 + t851;
t866 = t843 * t806;
t864 = qJD(1) * t843;
t860 = t840 + t875;
t859 = qJD(1) * t874;
t857 = qJ(5) * t864;
t795 = -t843 * g(3) + t806 * t841;
t854 = -pkin(3) * qJD(3) + t873;
t823 = -t844 * g(1) - t842 * g(2);
t815 = t860 * qJDD(1);
t818 = t860 * t846;
t796 = t815 * t844 - t818 * t842;
t853 = t815 * t842 + t818 * t844;
t802 = t820 * t843 - t826 * t841;
t852 = -qJDD(1) * qJ(2) - t823;
t811 = (pkin(3) * t841 - qJ(4) * t843) * qJD(1);
t850 = t845 * qJ(4) - t811 * t864 - qJDD(4) + t870;
t849 = -t845 * pkin(3) + qJDD(3) * qJ(4) + (qJD(3) * t873) + t795;
t848 = -t846 * t871 + t852;
t813 = -t856 + t861;
t847 = t848 + (t813 - t856) * qJ(4) - t877 * pkin(3);
t825 = -t845 - t868;
t821 = qJDD(3) - t858;
t819 = -(qJD(3) * pkin(4)) - t857;
t817 = qJDD(1) * t842 + t844 * t846;
t816 = qJDD(1) * t844 - t842 * t846;
t812 = 0.2e1 * t855 + t862;
t808 = qJDD(1) * pkin(1) - t851;
t807 = t846 * pkin(1) + t852 + t859;
t805 = t848 + t859;
t801 = -t821 * t841 + t825 * t843;
t798 = t821 * t843 + t825 * t841;
t794 = t866 + t870;
t792 = t798 * t842 + t812 * t844;
t790 = -t798 * t844 + t812 * t842;
t789 = qJDD(3) * pkin(3) + t850 + t866;
t788 = -t841 * qJD(1) * t811 + t849;
t787 = -t794 * t841 + t795 * t843;
t786 = t794 * t843 + t795 * t841;
t785 = (t843 * t854 + t874) * qJD(1) + t847;
t784 = -pkin(4) * t858 + (qJD(1) * t872 + t806) * t843 + (pkin(3) + pkin(4)) * qJDD(3) + (t813 + t856) * qJ(5) + t850;
t783 = (t819 + t857) * qJD(3) + (-pkin(4) * t869 + qJ(5) * qJDD(1) + (t872 - t811) * qJD(1)) * t841 + t849;
t782 = qJDD(5) - t877 * pkin(4) - qJ(5) * t868 + (t874 + (t819 + t854) * t843) * qJD(1) + t847;
t781 = t788 * t843 - t789 * t841;
t780 = t788 * t841 + t789 * t843;
t779 = t783 * t843 - t784 * t841;
t778 = t783 * t841 + t784 * t843;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t817, -t816, 0, -t822 * t842 + t823 * t844, 0, 0, 0, 0, 0, 0, 0, t817, t816, -t807 * t844 - t808 * t842, 0, 0, 0, 0, 0, 0, t792, -t793, -t853, t786 * t842 - t805 * t844, 0, 0, 0, 0, 0, 0, t792, -t853, t793, t780 * t842 - t785 * t844, 0, 0, 0, 0, 0, 0, t792, t793, t853, t778 * t842 - t782 * t844; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t816, -t817, 0, t822 * t844 + t823 * t842, 0, 0, 0, 0, 0, 0, 0, -t816, t817, -t807 * t842 + t808 * t844, 0, 0, 0, 0, 0, 0, t790, t876, t796, -t786 * t844 - t805 * t842, 0, 0, 0, 0, 0, 0, t790, t796, -t876, -t780 * t844 - t785 * t842, 0, 0, 0, 0, 0, 0, t790, -t876, -t796, -t778 * t844 - t782 * t842; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t801, -t802, 0, t787, 0, 0, 0, 0, 0, 0, t801, 0, t802, t781, 0, 0, 0, 0, 0, 0, t801, t802, 0, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t846, -qJDD(1), 0, t823, 0, 0, 0, 0, 0, 0, 0, t846, qJDD(1), -t807, 0, 0, 0, 0, 0, 0, t812, t814, -t818, -t805, 0, 0, 0, 0, 0, 0, t812, -t818, -t814, -t785, 0, 0, 0, 0, 0, 0, t812, -t814, t818, -t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t846, 0, t822, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t846, t808, 0, 0, 0, 0, 0, 0, -t798, t799, t815, -t786, 0, 0, 0, 0, 0, 0, -t798, t815, -t799, -t780, 0, 0, 0, 0, 0, 0, -t798, -t799, -t815, -t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t801, -t802, 0, t787, 0, 0, 0, 0, 0, 0, t801, 0, t802, t781, 0, 0, 0, 0, 0, 0, t801, t802, 0, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t801, -t802, 0, t787, 0, 0, 0, 0, 0, 0, t801, 0, t802, t781, 0, 0, 0, 0, 0, 0, t801, t802, 0, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t846, -qJDD(1), t807, 0, 0, 0, 0, 0, 0, -t812, -t814, t818, t805, 0, 0, 0, 0, 0, 0, -t812, t818, t814, t785, 0, 0, 0, 0, 0, 0, -t812, t814, -t818, t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t846, -t808, 0, 0, 0, 0, 0, 0, t798, -t799, -t815, t786, 0, 0, 0, 0, 0, 0, t798, -t815, t799, t780, 0, 0, 0, 0, 0, 0, t798, t799, t815, t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, -t820, -t862, t795, 0, 0, 0, 0, 0, 0, t825, -t862, t820, t788, 0, 0, 0, 0, 0, 0, t825, t820, t862, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t821, -t826, -t861, t794, 0, 0, 0, 0, 0, 0, t821, -t861, t826, t789, 0, 0, 0, 0, 0, 0, t821, t826, t861, t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, t814, -t818, -t805, 0, 0, 0, 0, 0, 0, t812, -t818, -t814, -t785, 0, 0, 0, 0, 0, 0, t812, -t814, t818, -t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, -t862, t820, t788, 0, 0, 0, 0, 0, 0, t825, t820, t862, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, -t818, -t814, -t785, 0, 0, 0, 0, 0, 0, t812, -t814, t818, -t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t821, t861, -t826, -t789, 0, 0, 0, 0, 0, 0, -t821, -t826, -t861, -t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, t820, t862, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t821, -t826, -t861, -t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t812, t814, -t818, t782;];
f_new_reg = t1;
