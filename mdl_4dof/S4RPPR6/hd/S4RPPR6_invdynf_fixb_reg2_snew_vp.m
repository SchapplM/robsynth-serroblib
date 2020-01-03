% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:52
% EndTime: 2019-12-31 16:40:53
% DurationCPUTime: 0.93s
% Computational Cost: add. (1364->135), mult. (3404->149), div. (0->0), fcn. (2123->6), ass. (0->94)
t731 = 2 * qJD(2);
t700 = sin(qJ(1));
t702 = cos(qJ(1));
t683 = -t702 * g(1) - t700 * g(2);
t704 = qJD(1) ^ 2;
t674 = -t704 * pkin(1) + qJDD(1) * qJ(2) + t683;
t730 = (qJD(1) * t731) + t674;
t697 = sin(pkin(6));
t695 = t697 ^ 2;
t698 = cos(pkin(6));
t696 = t698 ^ 2;
t720 = t695 + t696;
t679 = t720 * t704;
t699 = sin(qJ(4));
t701 = cos(qJ(4));
t712 = t697 * t699 + t698 * t701;
t648 = t712 * qJDD(1);
t671 = t712 * qJD(1);
t729 = t671 ^ 2;
t718 = qJD(1) * t698;
t719 = qJD(1) * t697;
t673 = -t699 * t718 + t701 * t719;
t728 = t673 ^ 2;
t727 = 2 * qJD(4);
t726 = t698 * g(3);
t725 = qJ(3) * t697;
t724 = t673 * t671;
t723 = t695 * t704;
t722 = t696 * t704;
t721 = t698 * t704;
t692 = t697 * qJDD(1);
t693 = t698 * qJDD(1);
t716 = t700 * qJDD(1);
t715 = t702 * qJDD(1);
t682 = t700 * g(1) - t702 * g(2);
t714 = qJDD(3) + t726;
t659 = -t697 * g(3) + t730 * t698;
t713 = -pkin(2) * t698 - t725;
t670 = -t701 * t692 + t699 * t693;
t677 = t713 * qJD(1);
t651 = t677 * t718 + t659;
t711 = t674 + (t731 + t677) * qJD(1);
t675 = t697 * t679;
t710 = -t700 * t675 + t697 * t715;
t709 = t702 * t675 + t697 * t716;
t708 = t704 * qJ(2) - qJDD(2) + t682;
t707 = 0.2e1 * qJD(3) * t719 + t708;
t703 = qJD(4) ^ 2;
t684 = t697 * t721;
t681 = -t702 * t704 - t716;
t680 = -t700 * t704 + t715;
t678 = t720 * qJDD(1);
t676 = t698 * t679;
t669 = qJDD(1) * pkin(1) + t708;
t664 = -t703 - t728;
t663 = -t702 * t676 - t698 * t716;
t662 = -t700 * t676 + t698 * t715;
t661 = t702 * t678 - t700 * t679;
t660 = t700 * t678 + t702 * t679;
t658 = -t730 * t697 - t726;
t657 = -t671 * t727 - t670;
t656 = t673 * t727 + t648;
t655 = -qJDD(4) - t724;
t654 = qJDD(4) - t724;
t653 = -t703 - t729;
t652 = (pkin(1) - t713) * qJDD(1) + t707;
t650 = -t728 - t729;
t649 = t711 * t697 + t714;
t647 = (t725 + pkin(1) + (pkin(2) + pkin(3)) * t698) * qJDD(1) + t707 + (-t723 - t722) * pkin(5);
t646 = -pkin(3) * t722 - pkin(5) * t693 + t651;
t645 = t701 * t655 - t699 * t664;
t644 = t699 * t655 + t701 * t664;
t643 = (-pkin(3) * t721 - pkin(5) * qJDD(1) + t711) * t697 + t714;
t642 = -t697 * t658 + t698 * t659;
t641 = t698 * t658 + t697 * t659;
t640 = -t701 * t648 - t699 * t670;
t639 = -t699 * t648 + t701 * t670;
t638 = t701 * t653 - t699 * t654;
t637 = t699 * t653 + t701 * t654;
t636 = t697 * t649 + t698 * t651;
t635 = -t698 * t649 + t697 * t651;
t634 = t697 * t644 + t698 * t645;
t633 = -t698 * t644 + t697 * t645;
t632 = t699 * t643 + t701 * t646;
t631 = t701 * t643 - t699 * t646;
t630 = t697 * t639 + t698 * t640;
t629 = -t698 * t639 + t697 * t640;
t628 = t697 * t637 + t698 * t638;
t627 = -t698 * t637 + t697 * t638;
t626 = -t699 * t631 + t701 * t632;
t625 = t701 * t631 + t699 * t632;
t624 = t697 * t625 + t698 * t626;
t623 = -t698 * t625 + t697 * t626;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t681, -t680, 0, -t700 * t682 + t702 * t683, 0, 0, 0, 0, 0, 0, t663, t709, t661, t702 * t642 - t700 * t669, 0, 0, 0, 0, 0, 0, t663, t661, -t709, t702 * t636 - t700 * t652, 0, 0, 0, 0, 0, 0, t702 * t628 - t700 * t656, t702 * t634 - t700 * t657, t702 * t630 - t700 * t650, t702 * t624 - t700 * t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t680, t681, 0, t702 * t682 + t700 * t683, 0, 0, 0, 0, 0, 0, t662, -t710, t660, t700 * t642 + t702 * t669, 0, 0, 0, 0, 0, 0, t662, t660, t710, t700 * t636 + t702 * t652, 0, 0, 0, 0, 0, 0, t700 * t628 + t702 * t656, t700 * t634 + t702 * t657, t700 * t630 + t702 * t650, t700 * t624 + t702 * t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t641, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, 0, 0, 0, 0, 0, 0, t627, t633, t629, t623; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t704, -qJDD(1), 0, t683, 0, 0, 0, 0, 0, 0, -t676, t675, t678, t642, 0, 0, 0, 0, 0, 0, -t676, t678, -t675, t636, 0, 0, 0, 0, 0, 0, t628, t634, t630, t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t704, 0, t682, 0, 0, 0, 0, 0, 0, t693, -t692, t679, t669, 0, 0, 0, 0, 0, 0, t693, t679, t692, t652, 0, 0, 0, 0, 0, 0, t656, t657, t650, t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t641, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, 0, 0, 0, 0, 0, 0, t627, t633, t629, t623; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t722, t684, t693, t659, 0, 0, 0, 0, 0, 0, -t722, t693, -t684, t651, 0, 0, 0, 0, 0, 0, t638, t645, t640, t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t684, -t723, -t692, t658, 0, 0, 0, 0, 0, 0, t684, -t692, t723, -t649, 0, 0, 0, 0, 0, 0, -t637, -t644, -t639, -t625; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t693, t692, -t679, -t669, 0, 0, 0, 0, 0, 0, -t693, -t679, -t692, -t652, 0, 0, 0, 0, 0, 0, -t656, -t657, -t650, -t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t722, t693, -t684, t651, 0, 0, 0, 0, 0, 0, t638, t645, t640, t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t693, -t679, -t692, -t652, 0, 0, 0, 0, 0, 0, -t656, -t657, -t650, -t647; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t684, t692, -t723, t649, 0, 0, 0, 0, 0, 0, t637, t644, t639, t625; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t653, t655, -t648, t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t654, t664, t670, t631; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t656, t657, t650, t647;];
f_new_reg = t1;
