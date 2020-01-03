% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR9_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:37
% EndTime: 2019-12-31 16:56:38
% DurationCPUTime: 0.84s
% Computational Cost: add. (1672->140), mult. (3281->154), div. (0->0), fcn. (1966->6), ass. (0->92)
t688 = sin(qJ(3));
t707 = t688 * qJD(1);
t678 = qJD(4) + t707;
t716 = qJD(4) + t678;
t694 = qJD(1) ^ 2;
t689 = sin(qJ(1));
t692 = cos(qJ(1));
t675 = -t692 * g(1) - t689 * g(2);
t695 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t675;
t712 = pkin(5) + pkin(1);
t656 = t712 * t694 + t695;
t687 = sin(qJ(4));
t690 = cos(qJ(4));
t691 = cos(qJ(3));
t708 = qJD(1) * t691;
t661 = -t690 * qJD(3) + t687 * t708;
t715 = t661 ^ 2;
t663 = t687 * qJD(3) + t690 * t708;
t714 = t663 ^ 2;
t713 = t678 ^ 2;
t711 = t688 * g(3);
t710 = t663 * t661;
t685 = t688 ^ 2;
t686 = t691 ^ 2;
t709 = t685 + t686;
t706 = qJD(4) - t678;
t705 = qJD(1) * qJD(3);
t704 = t688 * qJDD(1);
t703 = t691 * qJDD(1);
t702 = t688 * t694 * t691;
t701 = t688 * t705;
t700 = t691 * t705;
t674 = t689 * g(1) - t692 * g(2);
t697 = -t694 * qJ(2) + qJDD(2) - t674;
t657 = -t712 * qJDD(1) + t697;
t649 = -t691 * g(3) + t688 * t657;
t666 = -t701 + t703;
t699 = t690 * qJDD(3) - t687 * t666;
t698 = -t687 * qJDD(3) - t690 * t666;
t665 = 0.2e1 * t700 + t704;
t696 = -qJDD(4) - t700 - t704;
t693 = qJD(3) ^ 2;
t677 = -t686 * t694 - t693;
t676 = -t685 * t694 - t693;
t673 = -qJDD(3) - t702;
t672 = qJDD(3) - t702;
t671 = t709 * t694;
t670 = t689 * qJDD(1) + t692 * t694;
t669 = t692 * qJDD(1) - t689 * t694;
t668 = t709 * qJDD(1);
t667 = -0.2e1 * t701 + t703;
t664 = (pkin(3) * t688 - pkin(6) * t691) * qJD(1);
t660 = qJDD(1) * pkin(1) - t697;
t659 = t694 * pkin(1) + t695;
t654 = t691 * t673 - t688 * t677;
t653 = -t688 * t672 + t691 * t676;
t652 = t688 * t673 + t691 * t677;
t651 = t691 * t672 + t688 * t676;
t650 = -t713 - t714;
t648 = t691 * t657 + t711;
t647 = -t713 - t715;
t646 = t696 - t710;
t645 = -t696 - t710;
t644 = -t714 - t715;
t643 = -t693 * pkin(3) + qJDD(3) * pkin(6) - t664 * t707 + t649;
t642 = qJDD(3) * pkin(3) + t693 * pkin(6) + t711 + (-qJD(1) * t664 + t657) * t691;
t641 = t706 * t661 + t698;
t640 = -t716 * t661 - t698;
t639 = -t706 * t663 + t699;
t638 = t716 * t663 - t699;
t637 = (-t666 + t701) * pkin(6) + t665 * pkin(3) - t656;
t636 = -t688 * t648 + t691 * t649;
t635 = t691 * t648 + t688 * t649;
t634 = t690 * t646 - t687 * t650;
t633 = t687 * t646 + t690 * t650;
t632 = -t687 * t645 + t690 * t647;
t631 = t690 * t645 + t687 * t647;
t630 = t690 * t639 - t687 * t641;
t629 = t687 * t639 + t690 * t641;
t628 = t687 * t637 + t690 * t643;
t627 = t690 * t637 - t687 * t643;
t626 = t691 * t634 + t688 * t640;
t625 = t688 * t634 - t691 * t640;
t624 = t691 * t632 + t688 * t638;
t623 = t688 * t632 - t691 * t638;
t622 = t691 * t630 + t688 * t644;
t621 = t688 * t630 - t691 * t644;
t620 = -t687 * t627 + t690 * t628;
t619 = t690 * t627 + t687 * t628;
t618 = t691 * t620 - t688 * t642;
t617 = t688 * t620 + t691 * t642;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t670, -t669, 0, -t689 * t674 + t692 * t675, 0, 0, 0, 0, 0, 0, 0, t670, t669, -t692 * t659 - t689 * t660, 0, 0, 0, 0, 0, 0, t689 * t651 + t692 * t665, t689 * t652 + t692 * t667, -t689 * t668 - t692 * t671, t689 * t635 - t692 * t656, 0, 0, 0, 0, 0, 0, t689 * t623 + t692 * t631, t689 * t625 + t692 * t633, t689 * t621 + t692 * t629, t689 * t617 + t692 * t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t669, -t670, 0, t692 * t674 + t689 * t675, 0, 0, 0, 0, 0, 0, 0, -t669, t670, -t689 * t659 + t692 * t660, 0, 0, 0, 0, 0, 0, -t692 * t651 + t689 * t665, -t692 * t652 + t689 * t667, t692 * t668 - t689 * t671, -t692 * t635 - t689 * t656, 0, 0, 0, 0, 0, 0, -t692 * t623 + t689 * t631, -t692 * t625 + t689 * t633, -t692 * t621 + t689 * t629, -t692 * t617 + t689 * t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t653, t654, 0, t636, 0, 0, 0, 0, 0, 0, t624, t626, t622, t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t694, -qJDD(1), 0, t675, 0, 0, 0, 0, 0, 0, 0, t694, qJDD(1), -t659, 0, 0, 0, 0, 0, 0, t665, t667, -t671, -t656, 0, 0, 0, 0, 0, 0, t631, t633, t629, t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t694, 0, t674, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t694, t660, 0, 0, 0, 0, 0, 0, -t651, -t652, t668, -t635, 0, 0, 0, 0, 0, 0, -t623, -t625, -t621, -t617; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t653, t654, 0, t636, 0, 0, 0, 0, 0, 0, t624, t626, t622, t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t653, t654, 0, t636, 0, 0, 0, 0, 0, 0, t624, t626, t622, t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t694, -qJDD(1), t659, 0, 0, 0, 0, 0, 0, -t665, -t667, t671, t656, 0, 0, 0, 0, 0, 0, -t631, -t633, -t629, -t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t694, -t660, 0, 0, 0, 0, 0, 0, t651, t652, -t668, t635, 0, 0, 0, 0, 0, 0, t623, t625, t621, t617; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t676, t673, -t704, t649, 0, 0, 0, 0, 0, 0, t632, t634, t630, t620; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t672, t677, -t703, t648, 0, 0, 0, 0, 0, 0, -t638, -t640, -t644, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t665, t667, -t671, -t656, 0, 0, 0, 0, 0, 0, t631, t633, t629, t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t647, t646, t639, t628; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t645, t650, t641, t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t638, t640, t644, -t642;];
f_new_reg = t1;
