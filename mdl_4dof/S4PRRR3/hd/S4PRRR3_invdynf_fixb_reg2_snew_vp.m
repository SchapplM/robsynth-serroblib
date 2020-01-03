% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:47
% EndTime: 2019-12-31 16:31:48
% DurationCPUTime: 0.88s
% Computational Cost: add. (1792->95), mult. (2890->143), div. (0->0), fcn. (2102->8), ass. (0->80)
t685 = qJD(2) + qJD(3);
t683 = t685 ^ 2;
t684 = qJDD(2) + qJDD(3);
t692 = sin(qJ(3));
t695 = cos(qJ(3));
t664 = t692 * t683 - t695 * t684;
t693 = sin(qJ(2));
t696 = cos(qJ(2));
t700 = -t695 * t683 - t692 * t684;
t647 = t696 * t664 - t693 * t700;
t689 = sin(pkin(7));
t690 = cos(pkin(7));
t709 = t693 * t664 + t696 * t700;
t713 = t689 * t647 + t690 * t709;
t712 = t690 * t647 - t689 * t709;
t691 = sin(qJ(4));
t706 = t691 * t684;
t694 = cos(qJ(4));
t705 = t694 * t684;
t674 = t689 * g(1) - t690 * g(2);
t675 = -t690 * g(1) - t689 * g(2);
t658 = t693 * t674 + t696 * t675;
t698 = qJD(2) ^ 2;
t652 = -t698 * pkin(2) + t658;
t657 = t696 * t674 - t693 * t675;
t699 = qJDD(2) * pkin(2) + t657;
t636 = t695 * t652 + t692 * t699;
t686 = t691 ^ 2;
t687 = t694 ^ 2;
t704 = t686 + t687;
t703 = qJD(4) * t685;
t635 = -t692 * t652 + t695 * t699;
t672 = t696 * qJDD(2) - t693 * t698;
t673 = -t693 * qJDD(2) - t696 * t698;
t702 = -t689 * t672 + t690 * t673;
t701 = t690 * t672 + t689 * t673;
t697 = qJD(4) ^ 2;
t688 = -g(3) + qJDD(1);
t678 = t691 * t683 * t694;
t677 = -t687 * t683 - t697;
t676 = -t686 * t683 - t697;
t671 = -qJDD(4) + t678;
t670 = qJDD(4) + t678;
t666 = t704 * t683;
t661 = t704 * t684;
t660 = -0.2e1 * t691 * t703 + t705;
t659 = 0.2e1 * t694 * t703 + t706;
t656 = t694 * t671 - t691 * t676;
t655 = -t691 * t670 + t694 * t677;
t654 = t691 * t671 + t694 * t676;
t653 = t694 * t670 + t691 * t677;
t646 = t695 * t661 - t692 * t666;
t643 = t692 * t661 + t695 * t666;
t642 = t695 * t656 + t692 * t659;
t641 = t695 * t655 - t692 * t660;
t640 = t692 * t656 - t695 * t659;
t639 = t692 * t655 + t695 * t660;
t638 = -t693 * t657 + t696 * t658;
t637 = t696 * t657 + t693 * t658;
t634 = -t683 * pkin(3) + t684 * pkin(6) + t636;
t633 = -t684 * pkin(3) - t683 * pkin(6) - t635;
t632 = -t693 * t643 + t696 * t646;
t631 = t696 * t643 + t693 * t646;
t630 = t694 * t634 + t691 * t688;
t629 = -t691 * t634 + t694 * t688;
t628 = -t693 * t640 + t696 * t642;
t627 = -t693 * t639 + t696 * t641;
t626 = t696 * t640 + t693 * t642;
t625 = t696 * t639 + t693 * t641;
t624 = -t692 * t635 + t695 * t636;
t623 = t695 * t635 + t692 * t636;
t622 = -t691 * t629 + t694 * t630;
t621 = t694 * t629 + t691 * t630;
t620 = t695 * t622 + t692 * t633;
t619 = t692 * t622 - t695 * t633;
t618 = -t693 * t623 + t696 * t624;
t617 = t696 * t623 + t693 * t624;
t616 = -t693 * t619 + t696 * t620;
t615 = t696 * t619 + t693 * t620;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t689 * t674 + t690 * t675, 0, 0, 0, 0, 0, 0, t702, -t701, 0, -t689 * t637 + t690 * t638, 0, 0, 0, 0, 0, 0, t713, t712, 0, -t689 * t617 + t690 * t618, 0, 0, 0, 0, 0, 0, -t689 * t625 + t690 * t627, -t689 * t626 + t690 * t628, -t689 * t631 + t690 * t632, -t689 * t615 + t690 * t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t690 * t674 + t689 * t675, 0, 0, 0, 0, 0, 0, t701, t702, 0, t690 * t637 + t689 * t638, 0, 0, 0, 0, 0, 0, -t712, t713, 0, t690 * t617 + t689 * t618, 0, 0, 0, 0, 0, 0, t690 * t625 + t689 * t627, t690 * t626 + t689 * t628, t690 * t631 + t689 * t632, t690 * t615 + t689 * t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, t653, t654, 0, t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t675, 0, 0, 0, 0, 0, 0, t673, -t672, 0, t638, 0, 0, 0, 0, 0, 0, t709, t647, 0, t618, 0, 0, 0, 0, 0, 0, t627, t628, t632, t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t674, 0, 0, 0, 0, 0, 0, t672, t673, 0, t637, 0, 0, 0, 0, 0, 0, -t647, t709, 0, t617, 0, 0, 0, 0, 0, 0, t625, t626, t631, t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, t653, t654, 0, t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t698, -qJDD(2), 0, t658, 0, 0, 0, 0, 0, 0, t700, t664, 0, t624, 0, 0, 0, 0, 0, 0, t641, t642, t646, t620; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t698, 0, t657, 0, 0, 0, 0, 0, 0, -t664, t700, 0, t623, 0, 0, 0, 0, 0, 0, t639, t640, t643, t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, t653, t654, 0, t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t683, -t684, 0, t636, 0, 0, 0, 0, 0, 0, t655, t656, t661, t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t684, -t683, 0, t635, 0, 0, 0, 0, 0, 0, t660, -t659, t666, -t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t688, 0, 0, 0, 0, 0, 0, t653, t654, 0, t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t677, t671, t705, t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t670, t676, -t706, t629; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t660, t659, -t666, t633;];
f_new_reg = t1;
