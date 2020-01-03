% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR8
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:24
% EndTime: 2019-12-31 16:55:25
% DurationCPUTime: 0.80s
% Computational Cost: add. (1874->139), mult. (3779->154), div. (0->0), fcn. (2352->6), ass. (0->92)
t681 = sin(qJ(4));
t682 = sin(qJ(3));
t684 = cos(qJ(4));
t685 = cos(qJ(3));
t653 = (-t681 * t685 - t682 * t684) * qJD(1);
t710 = t653 ^ 2;
t703 = qJD(1) * t685;
t655 = -t681 * t682 * qJD(1) + t684 * t703;
t709 = t655 ^ 2;
t677 = qJD(3) + qJD(4);
t708 = t677 ^ 2;
t707 = pkin(5) + pkin(1);
t706 = t655 * t653;
t679 = t682 ^ 2;
t688 = qJD(1) ^ 2;
t705 = t679 * t688;
t683 = sin(qJ(1));
t686 = cos(qJ(1));
t668 = t683 * g(1) - t686 * g(2);
t690 = -t688 * qJ(2) + qJDD(2) - t668;
t650 = -t707 * qJDD(1) + t690;
t640 = t682 * g(3) + t685 * t650;
t680 = t685 ^ 2;
t704 = t679 + t680;
t702 = qJD(4) - t677;
t701 = qJD(4) + t677;
t700 = qJD(1) * qJD(3);
t699 = t682 * qJDD(1);
t698 = t685 * qJDD(1);
t697 = -qJDD(3) - qJDD(4);
t696 = t682 * t688 * t685;
t695 = t682 * t700;
t694 = t685 * t700;
t641 = -t685 * g(3) + t682 * t650;
t659 = -t694 - t699;
t660 = -t695 + t698;
t693 = t684 * t659 - t681 * t660;
t669 = -t686 * g(1) - t683 * g(2);
t666 = qJDD(3) - t696;
t692 = -t681 * t659 - t684 * t660;
t691 = -qJD(3) * pkin(3) + pkin(6) * t703;
t689 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t669;
t687 = qJD(3) ^ 2;
t671 = -t680 * t688 - t687;
t670 = -t687 - t705;
t667 = -qJDD(3) - t696;
t665 = t704 * t688;
t664 = t683 * qJDD(1) + t686 * t688;
t663 = t686 * qJDD(1) - t683 * t688;
t662 = t704 * qJDD(1);
t661 = -0.2e1 * t695 + t698;
t658 = 0.2e1 * t694 + t699;
t652 = qJDD(1) * pkin(1) - t690;
t651 = t688 * pkin(1) + t689;
t649 = t707 * t688 + t689;
t646 = -t708 - t709;
t645 = t685 * t667 - t682 * t671;
t644 = -t682 * t666 + t685 * t670;
t643 = t682 * t667 + t685 * t671;
t642 = t685 * t666 + t682 * t670;
t639 = t697 + t706;
t638 = -t697 + t706;
t637 = -t708 - t710;
t636 = -t709 - t710;
t635 = t659 * pkin(3) + t691 * t703 + (pkin(6) * t679 + t707) * t688 + t689;
t634 = -pkin(3) * t705 + t659 * pkin(6) + qJD(3) * t691 + t641;
t633 = (-t660 - t695) * pkin(6) + t666 * pkin(3) + t640;
t632 = -t682 * t640 + t685 * t641;
t631 = t685 * t640 + t682 * t641;
t630 = t684 * t639 - t681 * t646;
t629 = t681 * t639 + t684 * t646;
t628 = -t702 * t653 + t692;
t627 = t701 * t653 - t692;
t626 = -t702 * t655 + t693;
t625 = t701 * t655 - t693;
t624 = t684 * t637 - t681 * t638;
t623 = t681 * t637 + t684 * t638;
t622 = t681 * t633 + t684 * t634;
t621 = t684 * t633 - t681 * t634;
t620 = -t682 * t629 + t685 * t630;
t619 = t685 * t629 + t682 * t630;
t618 = t684 * t626 - t681 * t628;
t617 = t681 * t626 + t684 * t628;
t616 = -t682 * t623 + t685 * t624;
t615 = t685 * t623 + t682 * t624;
t614 = -t681 * t621 + t684 * t622;
t613 = t684 * t621 + t681 * t622;
t612 = -t682 * t617 + t685 * t618;
t611 = t685 * t617 + t682 * t618;
t610 = -t682 * t613 + t685 * t614;
t609 = t685 * t613 + t682 * t614;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t664, -t663, 0, -t683 * t668 + t686 * t669, 0, 0, 0, 0, 0, 0, 0, t664, t663, -t686 * t651 - t683 * t652, 0, 0, 0, 0, 0, 0, t683 * t642 + t686 * t658, t683 * t643 + t686 * t661, -t683 * t662 - t686 * t665, t683 * t631 - t686 * t649, 0, 0, 0, 0, 0, 0, t683 * t615 + t686 * t625, t683 * t619 + t686 * t627, t683 * t611 + t686 * t636, t683 * t609 - t686 * t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t663, -t664, 0, t686 * t668 + t683 * t669, 0, 0, 0, 0, 0, 0, 0, -t663, t664, -t683 * t651 + t686 * t652, 0, 0, 0, 0, 0, 0, -t686 * t642 + t683 * t658, -t686 * t643 + t683 * t661, t686 * t662 - t683 * t665, -t686 * t631 - t683 * t649, 0, 0, 0, 0, 0, 0, -t686 * t615 + t683 * t625, -t686 * t619 + t683 * t627, -t686 * t611 + t683 * t636, -t686 * t609 - t683 * t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t644, t645, 0, t632, 0, 0, 0, 0, 0, 0, t616, t620, t612, t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t688, -qJDD(1), 0, t669, 0, 0, 0, 0, 0, 0, 0, t688, qJDD(1), -t651, 0, 0, 0, 0, 0, 0, t658, t661, -t665, -t649, 0, 0, 0, 0, 0, 0, t625, t627, t636, -t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t688, 0, t668, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t688, t652, 0, 0, 0, 0, 0, 0, -t642, -t643, t662, -t631, 0, 0, 0, 0, 0, 0, -t615, -t619, -t611, -t609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t644, t645, 0, t632, 0, 0, 0, 0, 0, 0, t616, t620, t612, t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t644, t645, 0, t632, 0, 0, 0, 0, 0, 0, t616, t620, t612, t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t688, -qJDD(1), t651, 0, 0, 0, 0, 0, 0, -t658, -t661, t665, t649, 0, 0, 0, 0, 0, 0, -t625, -t627, -t636, t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t688, -t652, 0, 0, 0, 0, 0, 0, t642, t643, -t662, t631, 0, 0, 0, 0, 0, 0, t615, t619, t611, t609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t670, t667, -t699, t641, 0, 0, 0, 0, 0, 0, t624, t630, t618, t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, t671, -t698, t640, 0, 0, 0, 0, 0, 0, t623, t629, t617, t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t658, t661, -t665, -t649, 0, 0, 0, 0, 0, 0, t625, t627, t636, -t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t637, t639, t626, t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t638, t646, t628, t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t625, t627, t636, -t635;];
f_new_reg = t1;
