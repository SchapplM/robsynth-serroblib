% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:23
% EndTime: 2019-12-31 17:22:24
% DurationCPUTime: 1.00s
% Computational Cost: add. (2710->113), mult. (3462->150), div. (0->0), fcn. (2110->8), ass. (0->86)
t722 = qJD(1) + qJD(2);
t718 = qJD(3) + t722;
t716 = t718 ^ 2;
t721 = qJDD(1) + qJDD(2);
t717 = qJDD(3) + t721;
t726 = sin(qJ(3));
t730 = cos(qJ(3));
t693 = t726 * t716 - t730 * t717;
t727 = sin(qJ(2));
t731 = cos(qJ(2));
t737 = -t730 * t716 - t726 * t717;
t676 = t731 * t693 - t727 * t737;
t728 = sin(qJ(1));
t732 = cos(qJ(1));
t748 = t727 * t693 + t731 * t737;
t752 = t728 * t676 + t732 * t748;
t751 = t732 * t676 - t728 * t748;
t720 = t722 ^ 2;
t700 = t727 * t720 - t731 * t721;
t736 = -t731 * t720 - t727 * t721;
t747 = t728 * t700 + t732 * t736;
t746 = t732 * t700 - t728 * t736;
t725 = sin(qJ(4));
t741 = t725 * t717;
t729 = cos(qJ(4));
t740 = t729 * t717;
t712 = t728 * g(1) - t732 * g(2);
t704 = qJDD(1) * pkin(1) + t712;
t713 = -t732 * g(1) - t728 * g(2);
t734 = qJD(1) ^ 2;
t705 = -t734 * pkin(1) + t713;
t683 = t727 * t704 + t731 * t705;
t681 = -t720 * pkin(2) + t683;
t682 = t731 * t704 - t727 * t705;
t735 = t721 * pkin(2) + t682;
t665 = t730 * t681 + t726 * t735;
t723 = t725 ^ 2;
t724 = t729 ^ 2;
t739 = t723 + t724;
t738 = qJD(4) * t718;
t664 = -t726 * t681 + t730 * t735;
t733 = qJD(4) ^ 2;
t710 = -t728 * qJDD(1) - t732 * t734;
t709 = t732 * qJDD(1) - t728 * t734;
t708 = t729 * t716 * t725;
t707 = -t724 * t716 - t733;
t706 = -t723 * t716 - t733;
t703 = -qJDD(4) + t708;
t702 = qJDD(4) + t708;
t695 = t739 * t716;
t690 = t739 * t717;
t689 = -0.2e1 * t725 * t738 + t740;
t688 = 0.2e1 * t729 * t738 + t741;
t687 = t729 * t703 - t725 * t706;
t686 = -t725 * t702 + t729 * t707;
t685 = t725 * t703 + t729 * t706;
t684 = t729 * t702 + t725 * t707;
t675 = t730 * t690 - t726 * t695;
t672 = t726 * t690 + t730 * t695;
t671 = t730 * t687 + t726 * t688;
t670 = t730 * t686 - t726 * t689;
t669 = t726 * t687 - t730 * t688;
t668 = t726 * t686 + t730 * t689;
t667 = -t727 * t682 + t731 * t683;
t666 = t731 * t682 + t727 * t683;
t663 = -t727 * t672 + t731 * t675;
t662 = t731 * t672 + t727 * t675;
t661 = -t716 * pkin(3) + t717 * pkin(7) + t665;
t660 = -t717 * pkin(3) - t716 * pkin(7) - t664;
t659 = -t725 * g(3) + t729 * t661;
t658 = -t729 * g(3) - t725 * t661;
t657 = -t727 * t669 + t731 * t671;
t656 = -t727 * t668 + t731 * t670;
t655 = t731 * t669 + t727 * t671;
t654 = t731 * t668 + t727 * t670;
t653 = -t726 * t664 + t730 * t665;
t652 = t730 * t664 + t726 * t665;
t651 = -t725 * t658 + t729 * t659;
t650 = t729 * t658 + t725 * t659;
t649 = t730 * t651 + t726 * t660;
t648 = t726 * t651 - t730 * t660;
t647 = -t727 * t652 + t731 * t653;
t646 = t731 * t652 + t727 * t653;
t645 = -t727 * t648 + t731 * t649;
t644 = t731 * t648 + t727 * t649;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t710, -t709, 0, -t728 * t712 + t732 * t713, 0, 0, 0, 0, 0, 0, t747, t746, 0, -t728 * t666 + t732 * t667, 0, 0, 0, 0, 0, 0, t752, t751, 0, -t728 * t646 + t732 * t647, 0, 0, 0, 0, 0, 0, -t728 * t654 + t732 * t656, -t728 * t655 + t732 * t657, -t728 * t662 + t732 * t663, -t728 * t644 + t732 * t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t709, t710, 0, t732 * t712 + t728 * t713, 0, 0, 0, 0, 0, 0, -t746, t747, 0, t732 * t666 + t728 * t667, 0, 0, 0, 0, 0, 0, -t751, t752, 0, t732 * t646 + t728 * t647, 0, 0, 0, 0, 0, 0, t732 * t654 + t728 * t656, t732 * t655 + t728 * t657, t732 * t662 + t728 * t663, t732 * t644 + t728 * t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t684, t685, 0, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t734, -qJDD(1), 0, t713, 0, 0, 0, 0, 0, 0, t736, t700, 0, t667, 0, 0, 0, 0, 0, 0, t748, t676, 0, t647, 0, 0, 0, 0, 0, 0, t656, t657, t663, t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t734, 0, t712, 0, 0, 0, 0, 0, 0, -t700, t736, 0, t666, 0, 0, 0, 0, 0, 0, -t676, t748, 0, t646, 0, 0, 0, 0, 0, 0, t654, t655, t662, t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t684, t685, 0, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t720, -t721, 0, t683, 0, 0, 0, 0, 0, 0, t737, t693, 0, t653, 0, 0, 0, 0, 0, 0, t670, t671, t675, t649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t721, -t720, 0, t682, 0, 0, 0, 0, 0, 0, -t693, t737, 0, t652, 0, 0, 0, 0, 0, 0, t668, t669, t672, t648; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t684, t685, 0, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t716, -t717, 0, t665, 0, 0, 0, 0, 0, 0, t686, t687, t690, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t717, -t716, 0, t664, 0, 0, 0, 0, 0, 0, t689, -t688, t695, -t660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t684, t685, 0, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, t703, t740, t659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t702, t706, -t741, t658; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t689, t688, -t695, t660;];
f_new_reg = t1;
