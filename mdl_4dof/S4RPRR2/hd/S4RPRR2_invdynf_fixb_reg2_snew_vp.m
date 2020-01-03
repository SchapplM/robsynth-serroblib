% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:21
% EndTime: 2019-12-31 16:48:22
% DurationCPUTime: 0.97s
% Computational Cost: add. (2171->105), mult. (3462->149), div. (0->0), fcn. (2110->8), ass. (0->84)
t711 = qJD(1) + qJD(3);
t709 = t711 ^ 2;
t710 = qJDD(1) + qJDD(3);
t718 = sin(qJ(3));
t721 = cos(qJ(3));
t687 = t718 * t709 - t721 * t710;
t715 = sin(pkin(7));
t716 = cos(pkin(7));
t726 = -t721 * t709 - t718 * t710;
t669 = t716 * t687 - t715 * t726;
t719 = sin(qJ(1));
t722 = cos(qJ(1));
t735 = t715 * t687 + t716 * t726;
t739 = t719 * t669 + t722 * t735;
t738 = t722 * t669 - t719 * t735;
t717 = sin(qJ(4));
t732 = t717 * t710;
t720 = cos(qJ(4));
t731 = t720 * t710;
t702 = t719 * g(1) - t722 * g(2);
t691 = qJDD(1) * pkin(1) + t702;
t703 = -t722 * g(1) - t719 * g(2);
t724 = qJD(1) ^ 2;
t692 = -t724 * pkin(1) + t703;
t675 = t715 * t691 + t716 * t692;
t673 = -t724 * pkin(2) + t675;
t674 = t716 * t691 - t715 * t692;
t725 = qJDD(1) * pkin(2) + t674;
t657 = t721 * t673 + t718 * t725;
t712 = t717 ^ 2;
t713 = t720 ^ 2;
t730 = t712 + t713;
t729 = qJD(4) * t711;
t656 = -t718 * t673 + t721 * t725;
t695 = -t715 * qJDD(1) - t716 * t724;
t696 = t716 * qJDD(1) - t715 * t724;
t728 = t722 * t695 - t719 * t696;
t727 = t719 * t695 + t722 * t696;
t723 = qJD(4) ^ 2;
t714 = -g(3) + qJDD(2);
t701 = t720 * t709 * t717;
t700 = -t713 * t709 - t723;
t699 = -t712 * t709 - t723;
t698 = -t719 * qJDD(1) - t722 * t724;
t697 = t722 * qJDD(1) - t719 * t724;
t694 = -qJDD(4) + t701;
t693 = qJDD(4) + t701;
t689 = t730 * t709;
t684 = t730 * t710;
t681 = -0.2e1 * t717 * t729 + t731;
t680 = 0.2e1 * t720 * t729 + t732;
t679 = t720 * t694 - t717 * t699;
t678 = -t717 * t693 + t720 * t700;
t677 = t717 * t694 + t720 * t699;
t676 = t720 * t693 + t717 * t700;
t672 = t721 * t684 - t718 * t689;
t671 = t718 * t684 + t721 * t689;
t663 = t721 * t679 + t718 * t680;
t662 = t721 * t678 - t718 * t681;
t661 = t718 * t679 - t721 * t680;
t660 = t718 * t678 + t721 * t681;
t659 = -t715 * t674 + t716 * t675;
t658 = t716 * t674 + t715 * t675;
t655 = -t715 * t671 + t716 * t672;
t654 = t716 * t671 + t715 * t672;
t653 = -t709 * pkin(3) + t710 * pkin(6) + t657;
t652 = -t710 * pkin(3) - t709 * pkin(6) - t656;
t651 = t720 * t653 + t717 * t714;
t650 = -t717 * t653 + t720 * t714;
t649 = -t715 * t661 + t716 * t663;
t648 = -t715 * t660 + t716 * t662;
t647 = t716 * t661 + t715 * t663;
t646 = t716 * t660 + t715 * t662;
t645 = -t718 * t656 + t721 * t657;
t644 = t721 * t656 + t718 * t657;
t643 = -t717 * t650 + t720 * t651;
t642 = t720 * t650 + t717 * t651;
t641 = t721 * t643 + t718 * t652;
t640 = t718 * t643 - t721 * t652;
t639 = -t715 * t644 + t716 * t645;
t638 = t716 * t644 + t715 * t645;
t637 = -t715 * t640 + t716 * t641;
t636 = t716 * t640 + t715 * t641;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t698, -t697, 0, -t719 * t702 + t722 * t703, 0, 0, 0, 0, 0, 0, t728, -t727, 0, -t719 * t658 + t722 * t659, 0, 0, 0, 0, 0, 0, t739, t738, 0, -t719 * t638 + t722 * t639, 0, 0, 0, 0, 0, 0, -t719 * t646 + t722 * t648, -t719 * t647 + t722 * t649, -t719 * t654 + t722 * t655, -t719 * t636 + t722 * t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t697, t698, 0, t722 * t702 + t719 * t703, 0, 0, 0, 0, 0, 0, t727, t728, 0, t722 * t658 + t719 * t659, 0, 0, 0, 0, 0, 0, -t738, t739, 0, t722 * t638 + t719 * t639, 0, 0, 0, 0, 0, 0, t722 * t646 + t719 * t648, t722 * t647 + t719 * t649, t722 * t654 + t719 * t655, t722 * t636 + t719 * t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t724, -qJDD(1), 0, t703, 0, 0, 0, 0, 0, 0, t695, -t696, 0, t659, 0, 0, 0, 0, 0, 0, t735, t669, 0, t639, 0, 0, 0, 0, 0, 0, t648, t649, t655, t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t724, 0, t702, 0, 0, 0, 0, 0, 0, t696, t695, 0, t658, 0, 0, 0, 0, 0, 0, -t669, t735, 0, t638, 0, 0, 0, 0, 0, 0, t646, t647, t654, t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t724, -qJDD(1), 0, t675, 0, 0, 0, 0, 0, 0, t726, t687, 0, t645, 0, 0, 0, 0, 0, 0, t662, t663, t672, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t724, 0, t674, 0, 0, 0, 0, 0, 0, -t687, t726, 0, t644, 0, 0, 0, 0, 0, 0, t660, t661, t671, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t709, -t710, 0, t657, 0, 0, 0, 0, 0, 0, t678, t679, t684, t643; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, -t709, 0, t656, 0, 0, 0, 0, 0, 0, t681, -t680, t689, -t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t700, t694, t731, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t693, t699, -t732, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t681, t680, -t689, t652;];
f_new_reg = t1;
