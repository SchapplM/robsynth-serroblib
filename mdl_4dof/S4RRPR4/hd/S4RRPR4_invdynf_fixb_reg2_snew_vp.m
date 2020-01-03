% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:45
% EndTime: 2019-12-31 17:02:46
% DurationCPUTime: 1.10s
% Computational Cost: add. (3582->129), mult. (5218->186), div. (0->0), fcn. (3397->8), ass. (0->105)
t878 = 2 * qJD(4);
t844 = sin(pkin(7));
t845 = cos(pkin(7));
t846 = sin(qJ(4));
t849 = cos(qJ(4));
t877 = -t844 * t846 + t845 * t849;
t843 = qJD(1) + qJD(2);
t839 = t843 ^ 2;
t847 = sin(qJ(2));
t840 = qJDD(1) + qJDD(2);
t850 = cos(qJ(2));
t862 = t850 * t840;
t821 = t847 * t839 - t862;
t848 = sin(qJ(1));
t851 = cos(qJ(1));
t863 = t847 * t840;
t858 = -t850 * t839 - t863;
t876 = t848 * t821 + t851 * t858;
t875 = t851 * t821 - t848 * t858;
t857 = t844 * t849 + t845 * t846;
t808 = t857 * t840;
t841 = t844 ^ 2;
t842 = t845 ^ 2;
t861 = t841 + t842;
t818 = t861 * t839;
t809 = t877 * t843;
t872 = t809 ^ 2;
t811 = t857 * t843;
t871 = t811 ^ 2;
t870 = t811 * t809;
t869 = t839 * t845;
t868 = t841 * t839;
t867 = t842 * t839;
t866 = t844 * t840;
t836 = t845 * t840;
t832 = -t851 * g(1) - t848 * g(2);
t853 = qJD(1) ^ 2;
t823 = -t853 * pkin(1) + t832;
t831 = t848 * g(1) - t851 * g(2);
t856 = qJDD(1) * pkin(1) + t831;
t804 = t850 * t823 + t847 * t856;
t860 = qJD(3) * t843;
t859 = -t845 * g(3) - 0.2e1 * t844 * t860;
t803 = -t847 * t823 + t850 * t856;
t796 = -t839 * pkin(2) + t840 * qJ(3) + t804;
t788 = -t844 * g(3) + (t796 + 0.2e1 * t860) * t845;
t785 = t877 * t840;
t790 = -t840 * pkin(2) - t839 * qJ(3) + qJDD(3) - t803;
t852 = qJD(4) ^ 2;
t826 = t844 * t869;
t825 = -t848 * qJDD(1) - t851 * t853;
t824 = t851 * qJDD(1) - t848 * t853;
t814 = t861 * t840;
t813 = t845 * t818;
t812 = t844 * t818;
t805 = -t852 - t871;
t802 = -t850 * t813 - t845 * t863;
t801 = t850 * t812 + t844 * t863;
t800 = -t847 * t813 + t845 * t862;
t799 = t847 * t812 - t844 * t862;
t798 = t850 * t814 - t847 * t818;
t797 = t847 * t814 + t850 * t818;
t795 = t809 * t878 + t808;
t794 = t811 * t878 - t785;
t793 = -qJDD(4) + t870;
t792 = qJDD(4) + t870;
t791 = -t852 - t872;
t787 = -t844 * t796 + t859;
t786 = -t871 - t872;
t784 = -pkin(3) * t836 + t790 + (-t867 - t868) * pkin(6);
t783 = -t847 * t803 + t850 * t804;
t782 = t850 * t803 + t847 * t804;
t781 = t849 * t793 - t846 * t805;
t780 = t846 * t793 + t849 * t805;
t779 = -pkin(3) * t867 + pkin(6) * t836 + t788;
t778 = (pkin(3) * t869 - pkin(6) * t840 - t796) * t844 + t859;
t777 = t849 * t785 + t846 * t808;
t776 = t846 * t785 - t849 * t808;
t775 = t849 * t791 - t846 * t792;
t774 = t846 * t791 + t849 * t792;
t773 = -t844 * t787 + t845 * t788;
t772 = t845 * t787 + t844 * t788;
t771 = -t844 * t780 + t845 * t781;
t770 = t845 * t780 + t844 * t781;
t769 = t846 * t778 + t849 * t779;
t768 = t849 * t778 - t846 * t779;
t767 = t850 * t773 + t847 * t790;
t766 = t847 * t773 - t850 * t790;
t765 = -t844 * t776 + t845 * t777;
t764 = t845 * t776 + t844 * t777;
t763 = -t844 * t774 + t845 * t775;
t762 = t845 * t774 + t844 * t775;
t761 = t850 * t771 + t847 * t795;
t760 = t847 * t771 - t850 * t795;
t759 = t850 * t763 + t847 * t794;
t758 = t847 * t763 - t850 * t794;
t757 = t850 * t765 + t847 * t786;
t756 = t847 * t765 - t850 * t786;
t755 = -t846 * t768 + t849 * t769;
t754 = t849 * t768 + t846 * t769;
t753 = -t844 * t754 + t845 * t755;
t752 = t845 * t754 + t844 * t755;
t751 = t850 * t753 + t847 * t784;
t750 = t847 * t753 - t850 * t784;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t825, -t824, 0, -t848 * t831 + t851 * t832, 0, 0, 0, 0, 0, 0, t876, t875, 0, -t848 * t782 + t851 * t783, 0, 0, 0, 0, 0, 0, -t848 * t800 + t851 * t802, -t848 * t799 + t851 * t801, -t848 * t797 + t851 * t798, -t848 * t766 + t851 * t767, 0, 0, 0, 0, 0, 0, -t848 * t758 + t851 * t759, -t848 * t760 + t851 * t761, -t848 * t756 + t851 * t757, -t848 * t750 + t851 * t751; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t824, t825, 0, t851 * t831 + t848 * t832, 0, 0, 0, 0, 0, 0, -t875, t876, 0, t851 * t782 + t848 * t783, 0, 0, 0, 0, 0, 0, t851 * t800 + t848 * t802, t851 * t799 + t848 * t801, t851 * t797 + t848 * t798, t851 * t766 + t848 * t767, 0, 0, 0, 0, 0, 0, t851 * t758 + t848 * t759, t851 * t760 + t848 * t761, t851 * t756 + t848 * t757, t851 * t750 + t848 * t751; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t772, 0, 0, 0, 0, 0, 0, t762, t770, t764, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t853, -qJDD(1), 0, t832, 0, 0, 0, 0, 0, 0, t858, t821, 0, t783, 0, 0, 0, 0, 0, 0, t802, t801, t798, t767, 0, 0, 0, 0, 0, 0, t759, t761, t757, t751; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t853, 0, t831, 0, 0, 0, 0, 0, 0, -t821, t858, 0, t782, 0, 0, 0, 0, 0, 0, t800, t799, t797, t766, 0, 0, 0, 0, 0, 0, t758, t760, t756, t750; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t772, 0, 0, 0, 0, 0, 0, t762, t770, t764, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t839, -t840, 0, t804, 0, 0, 0, 0, 0, 0, -t813, t812, t814, t773, 0, 0, 0, 0, 0, 0, t763, t771, t765, t753; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t840, -t839, 0, t803, 0, 0, 0, 0, 0, 0, t836, -t866, t818, -t790, 0, 0, 0, 0, 0, 0, -t794, -t795, -t786, -t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t772, 0, 0, 0, 0, 0, 0, t762, t770, t764, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t867, t826, t836, t788, 0, 0, 0, 0, 0, 0, t775, t781, t777, t755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t826, -t868, -t866, t787, 0, 0, 0, 0, 0, 0, t774, t780, t776, t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t836, t866, -t818, t790, 0, 0, 0, 0, 0, 0, t794, t795, t786, t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t791, t793, t785, t769; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t792, t805, -t808, t768; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t794, t795, t786, t784;];
f_new_reg = t1;
