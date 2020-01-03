% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPPRR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:28
% EndTime: 2019-12-31 18:01:29
% DurationCPUTime: 1.43s
% Computational Cost: add. (4171->151), mult. (6286->159), div. (0->0), fcn. (2678->8), ass. (0->89)
t769 = -qJD(1) + qJD(4);
t767 = t769 ^ 2;
t768 = qJDD(1) - qJDD(4);
t777 = sin(qJ(4));
t780 = cos(qJ(4));
t746 = t777 * t767 + t780 * t768;
t774 = sin(pkin(8));
t775 = cos(pkin(8));
t787 = -t780 * t767 + t777 * t768;
t727 = t775 * t746 - t774 * t787;
t778 = sin(qJ(1));
t781 = cos(qJ(1));
t797 = t774 * t746 + t775 * t787;
t801 = t778 * t727 - t781 * t797;
t800 = t781 * t727 + t778 * t797;
t794 = -pkin(1) - pkin(2);
t776 = sin(qJ(5));
t793 = t776 * t768;
t779 = cos(qJ(5));
t792 = t779 * t768;
t783 = qJD(1) ^ 2;
t760 = -t781 * g(1) - t778 * g(2);
t786 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t760;
t737 = t794 * t783 + t786;
t759 = t778 * g(1) - t781 * g(2);
t785 = -t783 * qJ(2) + qJDD(2) - t759;
t738 = t794 * qJDD(1) + t785;
t724 = t775 * t737 + t774 * t738;
t718 = -t783 * pkin(3) + t724;
t723 = -t774 * t737 + t775 * t738;
t784 = -qJDD(1) * pkin(3) + t723;
t706 = t780 * t718 + t777 * t784;
t771 = t776 ^ 2;
t772 = t779 ^ 2;
t791 = t771 + t772;
t790 = qJD(5) * t769;
t705 = -t777 * t718 + t780 * t784;
t752 = -t774 * qJDD(1) + t775 * t783;
t753 = t775 * qJDD(1) + t774 * t783;
t789 = -t778 * t752 + t781 * t753;
t788 = t781 * t752 + t778 * t753;
t782 = qJD(5) ^ 2;
t773 = g(3) + qJDD(3);
t758 = t779 * t767 * t776;
t757 = -t772 * t767 - t782;
t756 = -t771 * t767 - t782;
t755 = t781 * qJDD(1) - t778 * t783;
t754 = t778 * qJDD(1) + t781 * t783;
t751 = -qJDD(5) + t758;
t750 = qJDD(5) + t758;
t748 = t791 * t767;
t743 = t791 * t768;
t742 = -0.2e1 * t776 * t790 - t792;
t741 = 0.2e1 * t779 * t790 - t793;
t740 = qJDD(1) * pkin(1) - t785;
t739 = -t783 * pkin(1) + t786;
t734 = t779 * t751 - t776 * t756;
t733 = -t776 * t750 + t779 * t757;
t732 = t776 * t751 + t779 * t756;
t731 = t779 * t750 + t776 * t757;
t730 = -t780 * t743 - t777 * t748;
t729 = -t777 * t743 + t780 * t748;
t722 = t780 * t734 + t777 * t741;
t721 = t780 * t733 - t777 * t742;
t720 = t777 * t734 - t780 * t741;
t719 = t777 * t733 + t780 * t742;
t714 = -t774 * t729 + t775 * t730;
t713 = t775 * t729 + t774 * t730;
t712 = -t774 * t723 + t775 * t724;
t711 = t775 * t723 + t774 * t724;
t710 = -t774 * t720 + t775 * t722;
t709 = -t774 * t719 + t775 * t721;
t708 = t775 * t720 + t774 * t722;
t707 = t775 * t719 + t774 * t721;
t704 = -t767 * pkin(4) - t768 * pkin(7) + t706;
t703 = t768 * pkin(4) - t767 * pkin(7) - t705;
t702 = t779 * t704 + t776 * t773;
t701 = -t776 * t704 + t779 * t773;
t700 = -t777 * t705 + t780 * t706;
t699 = t780 * t705 + t777 * t706;
t698 = -t776 * t701 + t779 * t702;
t697 = t779 * t701 + t776 * t702;
t696 = t780 * t698 + t777 * t703;
t695 = t777 * t698 - t780 * t703;
t694 = -t774 * t699 + t775 * t700;
t693 = t775 * t699 + t774 * t700;
t692 = -t774 * t695 + t775 * t696;
t691 = t775 * t695 + t774 * t696;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t754, -t755, 0, -t778 * t759 + t781 * t760, 0, 0, 0, 0, 0, 0, -t754, 0, t755, t781 * t739 - t778 * t740, 0, 0, 0, 0, 0, 0, -t788, t789, 0, t778 * t711 + t781 * t712, 0, 0, 0, 0, 0, 0, -t801, t800, 0, t778 * t693 + t781 * t694, 0, 0, 0, 0, 0, 0, t778 * t707 + t781 * t709, t778 * t708 + t781 * t710, t778 * t713 + t781 * t714, t778 * t691 + t781 * t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t755, -t754, 0, t781 * t759 + t778 * t760, 0, 0, 0, 0, 0, 0, t755, 0, t754, t778 * t739 + t781 * t740, 0, 0, 0, 0, 0, 0, t789, t788, 0, -t781 * t711 + t778 * t712, 0, 0, 0, 0, 0, 0, t800, t801, 0, -t781 * t693 + t778 * t694, 0, 0, 0, 0, 0, 0, -t781 * t707 + t778 * t709, -t781 * t708 + t778 * t710, -t781 * t713 + t778 * t714, -t781 * t691 + t778 * t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, -t731, -t732, 0, -t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, -qJDD(1), 0, t760, 0, 0, 0, 0, 0, 0, -t783, 0, qJDD(1), t739, 0, 0, 0, 0, 0, 0, -t752, t753, 0, t712, 0, 0, 0, 0, 0, 0, t797, t727, 0, t694, 0, 0, 0, 0, 0, 0, t709, t710, t714, t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t783, 0, t759, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t783, t740, 0, 0, 0, 0, 0, 0, t753, t752, 0, -t711, 0, 0, 0, 0, 0, 0, t727, -t797, 0, -t693, 0, 0, 0, 0, 0, 0, -t707, -t708, -t713, -t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, -t731, -t732, 0, -t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, 0, qJDD(1), t739, 0, 0, 0, 0, 0, 0, -t752, t753, 0, t712, 0, 0, 0, 0, 0, 0, t797, t727, 0, t694, 0, 0, 0, 0, 0, 0, t709, t710, t714, t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, 0, 0, 0, 0, 0, 0, -t731, -t732, 0, -t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t783, -t740, 0, 0, 0, 0, 0, 0, -t753, -t752, 0, t711, 0, 0, 0, 0, 0, 0, -t727, t797, 0, t693, 0, 0, 0, 0, 0, 0, t707, t708, t713, t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, qJDD(1), 0, t724, 0, 0, 0, 0, 0, 0, t787, t746, 0, t700, 0, 0, 0, 0, 0, 0, t721, t722, t730, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t783, 0, t723, 0, 0, 0, 0, 0, 0, -t746, t787, 0, t699, 0, 0, 0, 0, 0, 0, t719, t720, t729, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, 0, 0, 0, 0, 0, 0, t731, t732, 0, t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t767, t768, 0, t706, 0, 0, 0, 0, 0, 0, t733, t734, -t743, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t768, -t767, 0, t705, 0, 0, 0, 0, 0, 0, t742, -t741, t748, -t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, 0, 0, 0, 0, 0, 0, t731, t732, 0, t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t757, t751, -t792, t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t750, t756, t793, t701; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t742, t741, -t748, t703;];
f_new_reg = t1;
