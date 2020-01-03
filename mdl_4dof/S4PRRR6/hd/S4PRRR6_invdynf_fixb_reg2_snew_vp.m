% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:17
% EndTime: 2019-12-31 16:35:18
% DurationCPUTime: 0.96s
% Computational Cost: add. (2672->135), mult. (5495->189), div. (0->0), fcn. (3784->8), ass. (0->112)
t804 = qJD(3) + qJD(4);
t837 = qJD(4) + t804;
t811 = sin(qJ(4));
t814 = cos(qJ(4));
t815 = cos(qJ(3));
t812 = sin(qJ(3));
t828 = qJD(2) * t812;
t772 = -t814 * t815 * qJD(2) + t811 * t828;
t836 = t772 ^ 2;
t774 = (t811 * t815 + t812 * t814) * qJD(2);
t835 = t774 ^ 2;
t834 = t804 ^ 2;
t833 = t774 * t772;
t806 = t815 ^ 2;
t818 = qJD(2) ^ 2;
t832 = t806 * t818;
t808 = sin(pkin(7));
t809 = cos(pkin(7));
t788 = t808 * g(1) - t809 * g(2);
t831 = t808 * t788;
t830 = -g(3) + qJDD(1);
t789 = -t809 * g(1) - t808 * g(2);
t813 = sin(qJ(2));
t816 = cos(qJ(2));
t770 = t816 * t789 + t813 * t830;
t763 = -t818 * pkin(2) + qJDD(2) * pkin(5) + t770;
t753 = t815 * t763 - t812 * t788;
t805 = t812 ^ 2;
t829 = t805 + t806;
t827 = qJD(4) - t804;
t826 = qJD(2) * qJD(3);
t825 = t812 * qJDD(2);
t795 = t812 * t818 * t815;
t790 = qJDD(3) + t795;
t824 = -qJDD(3) - qJDD(4);
t823 = t812 * t826;
t822 = t815 * t826;
t750 = -t812 * t763 - t815 * t788;
t782 = t822 + t825;
t801 = t815 * qJDD(2);
t820 = -t801 + t823;
t821 = -t811 * t782 - t814 * t820;
t769 = -t813 * t789 + t816 * t830;
t762 = -qJDD(2) * pkin(2) - t818 * pkin(5) - t769;
t819 = -t814 * t782 + t811 * t820;
t817 = qJD(3) ^ 2;
t794 = -t817 - t832;
t793 = -t805 * t818 - t817;
t792 = qJD(3) * pkin(3) - pkin(6) * t828;
t791 = -qJDD(3) + t795;
t787 = t829 * t818;
t786 = t816 * qJDD(2) - t813 * t818;
t785 = -t813 * qJDD(2) - t816 * t818;
t784 = t829 * qJDD(2);
t783 = t801 - 0.2e1 * t823;
t781 = 0.2e1 * t822 + t825;
t776 = t809 * t788;
t768 = -t834 - t835;
t767 = t815 * t791 - t812 * t793;
t766 = -t812 * t790 + t815 * t794;
t765 = t812 * t791 + t815 * t793;
t764 = t815 * t790 + t812 * t794;
t761 = t816 * t784 - t813 * t787;
t760 = t813 * t784 + t816 * t787;
t758 = t824 - t833;
t757 = -t824 - t833;
t756 = -t834 - t836;
t755 = t816 * t767 + t813 * t781;
t754 = t816 * t766 - t813 * t783;
t752 = t813 * t767 - t816 * t781;
t751 = t813 * t766 + t816 * t783;
t749 = -t813 * t769 + t816 * t770;
t748 = t816 * t769 + t813 * t770;
t747 = -t835 - t836;
t746 = pkin(3) * t820 - pkin(6) * t832 + t792 * t828 + t762;
t745 = t814 * t758 - t811 * t768;
t744 = t811 * t758 + t814 * t768;
t743 = t827 * t772 + t819;
t742 = -t837 * t772 - t819;
t741 = -t827 * t774 + t821;
t740 = t837 * t774 - t821;
t739 = -pkin(3) * t832 - pkin(6) * t820 - qJD(3) * t792 + t753;
t738 = (-t782 + t822) * pkin(6) + t790 * pkin(3) + t750;
t737 = t814 * t756 - t811 * t757;
t736 = t811 * t756 + t814 * t757;
t735 = -t812 * t750 + t815 * t753;
t734 = t815 * t750 + t812 * t753;
t733 = t816 * t735 + t813 * t762;
t732 = t813 * t735 - t816 * t762;
t731 = -t812 * t744 + t815 * t745;
t730 = t815 * t744 + t812 * t745;
t729 = t814 * t741 - t811 * t743;
t728 = t811 * t741 + t814 * t743;
t727 = t811 * t738 + t814 * t739;
t726 = t814 * t738 - t811 * t739;
t725 = -t812 * t736 + t815 * t737;
t724 = t815 * t736 + t812 * t737;
t723 = t816 * t731 + t813 * t742;
t722 = t813 * t731 - t816 * t742;
t721 = t816 * t725 + t813 * t740;
t720 = t813 * t725 - t816 * t740;
t719 = -t812 * t728 + t815 * t729;
t718 = t815 * t728 + t812 * t729;
t717 = -t811 * t726 + t814 * t727;
t716 = t814 * t726 + t811 * t727;
t715 = t816 * t719 + t813 * t747;
t714 = t813 * t719 - t816 * t747;
t713 = -t812 * t716 + t815 * t717;
t712 = t815 * t716 + t812 * t717;
t711 = t816 * t713 + t813 * t746;
t710 = t813 * t713 - t816 * t746;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t809 * t789 - t831, 0, 0, 0, 0, 0, 0, t809 * t785, -t809 * t786, 0, t809 * t749 - t831, 0, 0, 0, 0, 0, 0, t809 * t754 + t808 * t764, t809 * t755 + t808 * t765, t809 * t761, t809 * t733 + t808 * t734, 0, 0, 0, 0, 0, 0, t809 * t721 + t808 * t724, t809 * t723 + t808 * t730, t809 * t715 + t808 * t718, t809 * t711 + t808 * t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t808 * t789 + t776, 0, 0, 0, 0, 0, 0, t808 * t785, -t808 * t786, 0, t808 * t749 + t776, 0, 0, 0, 0, 0, 0, t808 * t754 - t809 * t764, t808 * t755 - t809 * t765, t808 * t761, t808 * t733 - t809 * t734, 0, 0, 0, 0, 0, 0, t808 * t721 - t809 * t724, t808 * t723 - t809 * t730, t808 * t715 - t809 * t718, t808 * t711 - t809 * t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t830, 0, 0, 0, 0, 0, 0, t786, t785, 0, t748, 0, 0, 0, 0, 0, 0, t751, t752, t760, t732, 0, 0, 0, 0, 0, 0, t720, t722, t714, t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t789, 0, 0, 0, 0, 0, 0, t785, -t786, 0, t749, 0, 0, 0, 0, 0, 0, t754, t755, t761, t733, 0, 0, 0, 0, 0, 0, t721, t723, t715, t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t788, 0, 0, 0, 0, 0, 0, 0, 0, 0, t788, 0, 0, 0, 0, 0, 0, -t764, -t765, 0, -t734, 0, 0, 0, 0, 0, 0, -t724, -t730, -t718, -t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t830, 0, 0, 0, 0, 0, 0, t786, t785, 0, t748, 0, 0, 0, 0, 0, 0, t751, t752, t760, t732, 0, 0, 0, 0, 0, 0, t720, t722, t714, t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t818, -qJDD(2), 0, t770, 0, 0, 0, 0, 0, 0, t766, t767, t784, t735, 0, 0, 0, 0, 0, 0, t725, t731, t719, t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t818, 0, t769, 0, 0, 0, 0, 0, 0, t783, -t781, t787, -t762, 0, 0, 0, 0, 0, 0, -t740, -t742, -t747, -t746; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t788, 0, 0, 0, 0, 0, 0, t764, t765, 0, t734, 0, 0, 0, 0, 0, 0, t724, t730, t718, t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t794, t791, t801, t753, 0, 0, 0, 0, 0, 0, t737, t745, t729, t717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t790, t793, -t825, t750, 0, 0, 0, 0, 0, 0, t736, t744, t728, t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, t781, -t787, t762, 0, 0, 0, 0, 0, 0, t740, t742, t747, t746; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t756, t758, t741, t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t757, t768, t743, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t740, t742, t747, t746;];
f_new_reg = t1;
