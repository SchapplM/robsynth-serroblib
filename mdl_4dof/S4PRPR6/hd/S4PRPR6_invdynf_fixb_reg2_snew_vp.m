% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:50
% EndTime: 2019-12-31 16:24:51
% DurationCPUTime: 0.89s
% Computational Cost: add. (2018->115), mult. (4587->173), div. (0->0), fcn. (3244->8), ass. (0->97)
t744 = qJD(2) ^ 2;
t735 = sin(pkin(7));
t731 = t735 ^ 2;
t737 = cos(pkin(7));
t732 = t737 ^ 2;
t752 = t731 + t732;
t713 = t752 * t744;
t739 = sin(qJ(4));
t741 = cos(qJ(4));
t747 = t735 * t741 + t737 * t739;
t702 = t747 * qJDD(2);
t703 = (t735 * t739 - t737 * t741) * qJD(2);
t761 = t703 ^ 2;
t705 = t747 * qJD(2);
t760 = t705 ^ 2;
t759 = t705 * t703;
t758 = t731 * t744;
t757 = t732 * t744;
t736 = sin(pkin(6));
t738 = cos(pkin(6));
t716 = t736 * g(1) - t738 * g(2);
t756 = t736 * t716;
t755 = t737 * t744;
t754 = -g(3) + qJDD(1);
t717 = -t738 * g(1) - t736 * g(2);
t740 = sin(qJ(2));
t742 = cos(qJ(2));
t698 = t742 * t717 + t740 * t754;
t751 = qJD(2) * qJD(3);
t753 = -t737 * t716 - 0.2e1 * t735 * t751;
t750 = t735 * qJDD(2);
t728 = t737 * qJDD(2);
t749 = t740 * qJDD(2);
t748 = t742 * qJDD(2);
t692 = -t744 * pkin(2) + qJDD(2) * qJ(3) + t698;
t678 = -t735 * t716 + (t692 + 0.2e1 * t751) * t737;
t697 = -t740 * t717 + t742 * t754;
t676 = t741 * t728 - t739 * t750;
t691 = -qJDD(2) * pkin(2) - t744 * qJ(3) + qJDD(3) - t697;
t743 = qJD(4) ^ 2;
t718 = t735 * t755;
t715 = -t740 * t744 + t748;
t714 = -t742 * t744 - t749;
t712 = t752 * qJDD(2);
t710 = t738 * t716;
t707 = t737 * t713;
t706 = t735 * t713;
t699 = -t743 - t760;
t696 = -t742 * t707 - t737 * t749;
t695 = t742 * t706 + t735 * t749;
t694 = -t740 * t707 + t737 * t748;
t693 = t740 * t706 - t735 * t748;
t690 = t742 * t712 - t740 * t713;
t689 = t740 * t712 + t742 * t713;
t687 = -0.2e1 * t703 * qJD(4) + t702;
t686 = 0.2e1 * t705 * qJD(4) - t676;
t685 = -qJDD(4) - t759;
t684 = qJDD(4) - t759;
t683 = -t743 - t761;
t682 = -pkin(3) * t728 + t691 + (-t757 - t758) * pkin(5);
t681 = -t740 * t697 + t742 * t698;
t680 = t742 * t697 + t740 * t698;
t679 = -t760 - t761;
t677 = -t735 * t692 + t753;
t675 = t741 * t685 - t739 * t699;
t674 = t739 * t685 + t741 * t699;
t673 = -pkin(3) * t757 + pkin(5) * t728 + t678;
t672 = (pkin(3) * t755 - pkin(5) * qJDD(2) - t692) * t735 + t753;
t671 = t741 * t676 + t739 * t702;
t670 = t739 * t676 - t741 * t702;
t669 = t741 * t683 - t739 * t684;
t668 = t739 * t683 + t741 * t684;
t667 = -t735 * t677 + t737 * t678;
t666 = t737 * t677 + t735 * t678;
t665 = -t735 * t674 + t737 * t675;
t664 = t737 * t674 + t735 * t675;
t663 = t739 * t672 + t741 * t673;
t662 = t741 * t672 - t739 * t673;
t661 = t742 * t667 + t740 * t691;
t660 = t740 * t667 - t742 * t691;
t659 = -t735 * t670 + t737 * t671;
t658 = t737 * t670 + t735 * t671;
t657 = -t735 * t668 + t737 * t669;
t656 = t737 * t668 + t735 * t669;
t655 = t742 * t665 + t740 * t687;
t654 = t740 * t665 - t742 * t687;
t653 = t742 * t657 + t740 * t686;
t652 = t740 * t657 - t742 * t686;
t651 = t742 * t659 + t740 * t679;
t650 = t740 * t659 - t742 * t679;
t649 = -t739 * t662 + t741 * t663;
t648 = t741 * t662 + t739 * t663;
t647 = -t735 * t648 + t737 * t649;
t646 = t737 * t648 + t735 * t649;
t645 = t742 * t647 + t740 * t682;
t644 = t740 * t647 - t742 * t682;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t738 * t717 - t756, 0, 0, 0, 0, 0, 0, t738 * t714, -t738 * t715, 0, t738 * t681 - t756, 0, 0, 0, 0, 0, 0, t738 * t696, t738 * t695, t738 * t690, t738 * t661 + t736 * t666, 0, 0, 0, 0, 0, 0, t738 * t653 + t736 * t656, t738 * t655 + t736 * t664, t738 * t651 + t736 * t658, t738 * t645 + t736 * t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t736 * t717 + t710, 0, 0, 0, 0, 0, 0, t736 * t714, -t736 * t715, 0, t736 * t681 + t710, 0, 0, 0, 0, 0, 0, t736 * t696, t736 * t695, t736 * t690, t736 * t661 - t738 * t666, 0, 0, 0, 0, 0, 0, t736 * t653 - t738 * t656, t736 * t655 - t738 * t664, t736 * t651 - t738 * t658, t736 * t645 - t738 * t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t754, 0, 0, 0, 0, 0, 0, t715, t714, 0, t680, 0, 0, 0, 0, 0, 0, t694, t693, t689, t660, 0, 0, 0, 0, 0, 0, t652, t654, t650, t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t717, 0, 0, 0, 0, 0, 0, t714, -t715, 0, t681, 0, 0, 0, 0, 0, 0, t696, t695, t690, t661, 0, 0, 0, 0, 0, 0, t653, t655, t651, t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t666, 0, 0, 0, 0, 0, 0, -t656, -t664, -t658, -t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t754, 0, 0, 0, 0, 0, 0, t715, t714, 0, t680, 0, 0, 0, 0, 0, 0, t694, t693, t689, t660, 0, 0, 0, 0, 0, 0, t652, t654, t650, t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t744, -qJDD(2), 0, t698, 0, 0, 0, 0, 0, 0, -t707, t706, t712, t667, 0, 0, 0, 0, 0, 0, t657, t665, t659, t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t744, 0, t697, 0, 0, 0, 0, 0, 0, t728, -t750, t713, -t691, 0, 0, 0, 0, 0, 0, -t686, -t687, -t679, -t682; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t716, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, 0, 0, 0, 0, 0, 0, t656, t664, t658, t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t757, t718, t728, t678, 0, 0, 0, 0, 0, 0, t669, t675, t671, t649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t718, -t758, -t750, t677, 0, 0, 0, 0, 0, 0, t668, t674, t670, t648; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t728, t750, -t713, t691, 0, 0, 0, 0, 0, 0, t686, t687, t679, t682; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t683, t685, t676, t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t684, t699, -t702, t662; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t686, t687, t679, t682;];
f_new_reg = t1;
