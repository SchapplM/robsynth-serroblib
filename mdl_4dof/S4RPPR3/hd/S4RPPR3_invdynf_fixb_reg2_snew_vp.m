% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:06
% EndTime: 2019-12-31 16:38:07
% DurationCPUTime: 0.94s
% Computational Cost: add. (2370->126), mult. (5218->187), div. (0->0), fcn. (3397->8), ass. (0->101)
t811 = qJD(1) ^ 2;
t802 = sin(pkin(7));
t797 = t802 ^ 2;
t804 = cos(pkin(7));
t798 = t804 ^ 2;
t822 = t797 + t798;
t780 = t822 * t811;
t806 = sin(qJ(4));
t808 = cos(qJ(4));
t815 = t802 * t808 + t804 * t806;
t766 = t815 * qJDD(1);
t767 = (t802 * t806 - t804 * t808) * qJD(1);
t829 = t767 ^ 2;
t769 = t815 * qJD(1);
t828 = t769 ^ 2;
t827 = t769 * t767;
t826 = t797 * t811;
t825 = t798 * t811;
t824 = t804 * t811;
t807 = sin(qJ(1));
t809 = cos(qJ(1));
t784 = -t809 * g(1) - t807 * g(2);
t776 = -t811 * pkin(1) + t784;
t803 = sin(pkin(6));
t805 = cos(pkin(6));
t783 = t807 * g(1) - t809 * g(2);
t814 = qJDD(1) * pkin(1) + t783;
t756 = t805 * t776 + t803 * t814;
t800 = -g(3) + qJDD(2);
t821 = qJD(1) * qJD(3);
t823 = t804 * t800 - 0.2e1 * t802 * t821;
t820 = t802 * qJDD(1);
t819 = t803 * qJDD(1);
t795 = t804 * qJDD(1);
t818 = t805 * qJDD(1);
t749 = -t811 * pkin(2) + qJDD(1) * qJ(3) + t756;
t744 = t802 * t800 + (t749 + 0.2e1 * t821) * t804;
t755 = -t803 * t776 + t805 * t814;
t778 = -t805 * t811 - t819;
t779 = -t803 * t811 + t818;
t817 = t809 * t778 - t807 * t779;
t745 = t808 * t795 - t806 * t820;
t816 = t807 * t778 + t809 * t779;
t748 = -qJDD(1) * pkin(2) - t811 * qJ(3) + qJDD(3) - t755;
t810 = qJD(4) ^ 2;
t785 = t802 * t824;
t782 = -t807 * qJDD(1) - t809 * t811;
t781 = t809 * qJDD(1) - t807 * t811;
t777 = t822 * qJDD(1);
t775 = t804 * t780;
t774 = t802 * t780;
t763 = -t810 - t828;
t762 = -t805 * t775 - t803 * t795;
t761 = t805 * t774 + t802 * t819;
t760 = -t803 * t775 + t804 * t818;
t759 = t803 * t774 - t802 * t818;
t758 = t805 * t777 - t803 * t780;
t757 = t803 * t777 + t805 * t780;
t754 = -0.2e1 * t767 * qJD(4) + t766;
t753 = 0.2e1 * t769 * qJD(4) - t745;
t752 = -qJDD(4) - t827;
t751 = qJDD(4) - t827;
t750 = -t810 - t829;
t746 = -t828 - t829;
t743 = -t802 * t749 + t823;
t742 = -pkin(3) * t795 + t748 + (-t825 - t826) * pkin(5);
t741 = t808 * t752 - t806 * t763;
t740 = t806 * t752 + t808 * t763;
t739 = -pkin(3) * t825 + pkin(5) * t795 + t744;
t738 = (pkin(3) * t824 - pkin(5) * qJDD(1) - t749) * t802 + t823;
t737 = -t803 * t755 + t805 * t756;
t736 = t805 * t755 + t803 * t756;
t735 = t808 * t745 + t806 * t766;
t734 = t806 * t745 - t808 * t766;
t733 = t808 * t750 - t806 * t751;
t732 = t806 * t750 + t808 * t751;
t731 = -t802 * t743 + t804 * t744;
t730 = t804 * t743 + t802 * t744;
t729 = -t802 * t740 + t804 * t741;
t728 = t804 * t740 + t802 * t741;
t727 = t806 * t738 + t808 * t739;
t726 = t808 * t738 - t806 * t739;
t725 = -t802 * t734 + t804 * t735;
t724 = t804 * t734 + t802 * t735;
t723 = -t802 * t732 + t804 * t733;
t722 = t804 * t732 + t802 * t733;
t721 = t805 * t731 + t803 * t748;
t720 = t803 * t731 - t805 * t748;
t719 = t805 * t729 + t803 * t754;
t718 = t803 * t729 - t805 * t754;
t717 = t805 * t723 + t803 * t753;
t716 = t803 * t723 - t805 * t753;
t715 = t805 * t725 + t803 * t746;
t714 = t803 * t725 - t805 * t746;
t713 = -t806 * t726 + t808 * t727;
t712 = t808 * t726 + t806 * t727;
t711 = -t802 * t712 + t804 * t713;
t710 = t804 * t712 + t802 * t713;
t709 = t805 * t711 + t803 * t742;
t708 = t803 * t711 - t805 * t742;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t782, -t781, 0, -t807 * t783 + t809 * t784, 0, 0, 0, 0, 0, 0, t817, -t816, 0, -t807 * t736 + t809 * t737, 0, 0, 0, 0, 0, 0, -t807 * t760 + t809 * t762, -t807 * t759 + t809 * t761, -t807 * t757 + t809 * t758, -t807 * t720 + t809 * t721, 0, 0, 0, 0, 0, 0, -t807 * t716 + t809 * t717, -t807 * t718 + t809 * t719, -t807 * t714 + t809 * t715, -t807 * t708 + t809 * t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t781, t782, 0, t809 * t783 + t807 * t784, 0, 0, 0, 0, 0, 0, t816, t817, 0, t809 * t736 + t807 * t737, 0, 0, 0, 0, 0, 0, t809 * t760 + t807 * t762, t809 * t759 + t807 * t761, t809 * t757 + t807 * t758, t809 * t720 + t807 * t721, 0, 0, 0, 0, 0, 0, t809 * t716 + t807 * t717, t809 * t718 + t807 * t719, t809 * t714 + t807 * t715, t809 * t708 + t807 * t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t800, 0, 0, 0, 0, 0, 0, 0, 0, 0, t730, 0, 0, 0, 0, 0, 0, t722, t728, t724, t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t811, -qJDD(1), 0, t784, 0, 0, 0, 0, 0, 0, t778, -t779, 0, t737, 0, 0, 0, 0, 0, 0, t762, t761, t758, t721, 0, 0, 0, 0, 0, 0, t717, t719, t715, t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t811, 0, t783, 0, 0, 0, 0, 0, 0, t779, t778, 0, t736, 0, 0, 0, 0, 0, 0, t760, t759, t757, t720, 0, 0, 0, 0, 0, 0, t716, t718, t714, t708; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t800, 0, 0, 0, 0, 0, 0, 0, 0, 0, t730, 0, 0, 0, 0, 0, 0, t722, t728, t724, t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t811, -qJDD(1), 0, t756, 0, 0, 0, 0, 0, 0, -t775, t774, t777, t731, 0, 0, 0, 0, 0, 0, t723, t729, t725, t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t811, 0, t755, 0, 0, 0, 0, 0, 0, t795, -t820, t780, -t748, 0, 0, 0, 0, 0, 0, -t753, -t754, -t746, -t742; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t800, 0, 0, 0, 0, 0, 0, 0, 0, 0, t730, 0, 0, 0, 0, 0, 0, t722, t728, t724, t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t825, t785, t795, t744, 0, 0, 0, 0, 0, 0, t733, t741, t735, t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t785, -t826, -t820, t743, 0, 0, 0, 0, 0, 0, t732, t740, t734, t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t795, t820, -t780, t748, 0, 0, 0, 0, 0, 0, t753, t754, t746, t742; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t750, t752, t745, t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t751, t763, -t766, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t753, t754, t746, t742;];
f_new_reg = t1;
