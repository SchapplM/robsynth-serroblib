% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:43
% EndTime: 2019-12-31 17:01:44
% DurationCPUTime: 0.97s
% Computational Cost: add. (2370->108), mult. (3462->149), div. (0->0), fcn. (2110->8), ass. (0->84)
t711 = qJD(1) + qJD(2);
t709 = t711 ^ 2;
t710 = qJDD(1) + qJDD(2);
t715 = sin(pkin(7));
t716 = cos(pkin(7));
t685 = t715 * t709 - t716 * t710;
t718 = sin(qJ(2));
t721 = cos(qJ(2));
t727 = -t716 * t709 - t715 * t710;
t669 = t721 * t685 - t718 * t727;
t719 = sin(qJ(1));
t722 = cos(qJ(1));
t738 = t718 * t685 + t721 * t727;
t742 = t719 * t669 + t722 * t738;
t741 = t722 * t669 - t719 * t738;
t691 = t718 * t709 - t721 * t710;
t726 = -t721 * t709 - t718 * t710;
t737 = t719 * t691 + t722 * t726;
t736 = t722 * t691 - t719 * t726;
t717 = sin(qJ(4));
t731 = t717 * t710;
t720 = cos(qJ(4));
t730 = t720 * t710;
t703 = t719 * g(1) - t722 * g(2);
t694 = qJDD(1) * pkin(1) + t703;
t704 = -t722 * g(1) - t719 * g(2);
t724 = qJD(1) ^ 2;
t695 = -t724 * pkin(1) + t704;
t675 = t718 * t694 + t721 * t695;
t673 = -t709 * pkin(2) + t675;
t674 = t721 * t694 - t718 * t695;
t725 = t710 * pkin(2) + t674;
t657 = t716 * t673 + t715 * t725;
t712 = t717 ^ 2;
t713 = t720 ^ 2;
t729 = t712 + t713;
t728 = qJD(4) * t711;
t656 = -t715 * t673 + t716 * t725;
t723 = qJD(4) ^ 2;
t714 = -g(3) + qJDD(3);
t702 = t720 * t709 * t717;
t701 = -t713 * t709 - t723;
t700 = -t712 * t709 - t723;
t699 = -t719 * qJDD(1) - t722 * t724;
t698 = t722 * qJDD(1) - t719 * t724;
t697 = -qJDD(4) + t702;
t696 = qJDD(4) + t702;
t693 = t729 * t709;
t688 = t729 * t710;
t681 = -0.2e1 * t717 * t728 + t730;
t680 = 0.2e1 * t720 * t728 + t731;
t679 = t720 * t697 - t717 * t700;
t678 = -t717 * t696 + t720 * t701;
t677 = t717 * t697 + t720 * t700;
t676 = t720 * t696 + t717 * t701;
t672 = t716 * t688 - t715 * t693;
t671 = t715 * t688 + t716 * t693;
t663 = t716 * t679 + t715 * t680;
t662 = t716 * t678 - t715 * t681;
t661 = t715 * t679 - t716 * t680;
t660 = t715 * t678 + t716 * t681;
t659 = -t718 * t674 + t721 * t675;
t658 = t721 * t674 + t718 * t675;
t655 = -t718 * t671 + t721 * t672;
t654 = t721 * t671 + t718 * t672;
t653 = -t709 * pkin(3) + t710 * pkin(6) + t657;
t652 = -t710 * pkin(3) - t709 * pkin(6) - t656;
t651 = t720 * t653 + t717 * t714;
t650 = -t717 * t653 + t720 * t714;
t649 = -t718 * t661 + t721 * t663;
t648 = -t718 * t660 + t721 * t662;
t647 = t721 * t661 + t718 * t663;
t646 = t721 * t660 + t718 * t662;
t645 = -t715 * t656 + t716 * t657;
t644 = t716 * t656 + t715 * t657;
t643 = -t717 * t650 + t720 * t651;
t642 = t720 * t650 + t717 * t651;
t641 = t716 * t643 + t715 * t652;
t640 = t715 * t643 - t716 * t652;
t639 = -t718 * t644 + t721 * t645;
t638 = t721 * t644 + t718 * t645;
t637 = -t718 * t640 + t721 * t641;
t636 = t721 * t640 + t718 * t641;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t699, -t698, 0, -t719 * t703 + t722 * t704, 0, 0, 0, 0, 0, 0, t737, t736, 0, -t719 * t658 + t722 * t659, 0, 0, 0, 0, 0, 0, t742, t741, 0, -t719 * t638 + t722 * t639, 0, 0, 0, 0, 0, 0, -t719 * t646 + t722 * t648, -t719 * t647 + t722 * t649, -t719 * t654 + t722 * t655, -t719 * t636 + t722 * t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t698, t699, 0, t722 * t703 + t719 * t704, 0, 0, 0, 0, 0, 0, -t736, t737, 0, t722 * t658 + t719 * t659, 0, 0, 0, 0, 0, 0, -t741, t742, 0, t722 * t638 + t719 * t639, 0, 0, 0, 0, 0, 0, t722 * t646 + t719 * t648, t722 * t647 + t719 * t649, t722 * t654 + t719 * t655, t722 * t636 + t719 * t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t724, -qJDD(1), 0, t704, 0, 0, 0, 0, 0, 0, t726, t691, 0, t659, 0, 0, 0, 0, 0, 0, t738, t669, 0, t639, 0, 0, 0, 0, 0, 0, t648, t649, t655, t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t724, 0, t703, 0, 0, 0, 0, 0, 0, -t691, t726, 0, t658, 0, 0, 0, 0, 0, 0, -t669, t738, 0, t638, 0, 0, 0, 0, 0, 0, t646, t647, t654, t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t709, -t710, 0, t675, 0, 0, 0, 0, 0, 0, t727, t685, 0, t645, 0, 0, 0, 0, 0, 0, t662, t663, t672, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, -t709, 0, t674, 0, 0, 0, 0, 0, 0, -t685, t727, 0, t644, 0, 0, 0, 0, 0, 0, t660, t661, t671, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t709, -t710, 0, t657, 0, 0, 0, 0, 0, 0, t678, t679, t688, t643; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, -t709, 0, t656, 0, 0, 0, 0, 0, 0, t681, -t680, t693, -t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, 0, t676, t677, 0, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t701, t697, t730, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t696, t700, -t731, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t681, t680, -t693, t652;];
f_new_reg = t1;
