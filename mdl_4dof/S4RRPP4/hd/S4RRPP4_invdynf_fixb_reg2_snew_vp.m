% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPP4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:26
% EndTime: 2019-12-31 16:59:27
% DurationCPUTime: 0.78s
% Computational Cost: add. (806->132), mult. (1894->118), div. (0->0), fcn. (891->4), ass. (0->71)
t722 = sin(qJ(2));
t724 = cos(qJ(2));
t726 = qJD(1) ^ 2;
t733 = t724 * t726 * t722;
t702 = qJDD(2) - t733;
t718 = t722 ^ 2;
t750 = qJD(2) ^ 2;
t706 = t718 * t726 + t750;
t682 = t702 * t724 - t706 * t722;
t736 = qJD(1) * qJD(2);
t731 = t724 * t736;
t735 = t722 * qJDD(1);
t693 = 0.2e1 * t731 + t735;
t723 = sin(qJ(1));
t725 = cos(qJ(1));
t673 = t682 * t723 + t725 * t693;
t675 = t682 * t725 - t723 * t693;
t751 = -t733 - qJDD(2);
t749 = 2 * qJD(3);
t748 = -2 * qJD(4);
t747 = t724 * g(3);
t746 = qJDD(1) * pkin(1);
t719 = t724 ^ 2;
t745 = t719 * t726;
t744 = t722 * qJ(3);
t704 = -g(1) * t725 - g(2) * t723;
t689 = -pkin(1) * t726 + qJDD(1) * pkin(5) + t704;
t743 = t722 * t689;
t740 = t718 + t719;
t739 = qJD(1) * t722;
t738 = qJD(1) * t724;
t737 = qJD(2) * t724;
t734 = t724 * qJDD(1);
t732 = t722 * t736;
t703 = t723 * g(1) - t725 * g(2);
t685 = -t722 * g(3) + t724 * t689;
t696 = t740 * qJDD(1);
t699 = t740 * t726;
t677 = t696 * t725 - t699 * t723;
t676 = t696 * t723 + t699 * t725;
t679 = t702 * t722 + t706 * t724;
t730 = t726 * pkin(5) + t703;
t729 = -qJDD(2) * pkin(2) - t750 * qJ(3) + qJDD(3) + t747;
t692 = (-pkin(2) * t724 - t744) * qJD(1);
t728 = qJDD(2) * qJ(3) + qJD(2) * t749 + t692 * t738 + t685;
t694 = -t732 + t734;
t727 = t730 + (t694 - t732) * pkin(2);
t707 = -t745 - t750;
t700 = -qJD(2) * pkin(3) - qJ(4) * t739;
t698 = -qJDD(1) * t723 - t725 * t726;
t697 = qJDD(1) * t725 - t723 * t726;
t695 = -0.2e1 * t732 + t734;
t688 = t730 + t746;
t684 = -t743 - t747;
t681 = t707 * t724 + t722 * t751;
t678 = t707 * t722 - t724 * t751;
t674 = t681 * t725 - t695 * t723;
t672 = t681 * t723 + t695 * t725;
t671 = (qJD(1) * t692 + t689) * t722 + t729;
t670 = -pkin(2) * t750 + t728;
t669 = -t684 * t722 + t685 * t724;
t668 = t684 * t724 + t685 * t722;
t667 = qJ(3) * t693 + t739 * t749 + t727 + t746;
t666 = t743 - (t731 + t735) * qJ(4) + (qJ(4) * t737 + (t748 + t692) * t722) * qJD(1) + t729 + t751 * pkin(3);
t665 = -pkin(3) * t745 + t738 * t748 - t694 * qJ(4) + (-pkin(2) * qJD(2) + t700) * qJD(2) + t728;
t664 = qJDD(4) + t694 * pkin(3) - qJ(4) * t745 + (pkin(1) + t744) * qJDD(1) + (0.2e1 * qJ(3) * t737 + (t749 + t700) * t722) * qJD(1) + t727;
t663 = t670 * t724 + t671 * t722;
t662 = t670 * t722 - t671 * t724;
t661 = t665 * t724 + t666 * t722;
t660 = t665 * t722 - t666 * t724;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t698, -t697, 0, -t703 * t723 + t704 * t725, 0, 0, 0, 0, 0, 0, t674, -t675, t677, t669 * t725 - t688 * t723, 0, 0, 0, 0, 0, 0, t674, t677, t675, t663 * t725 - t667 * t723, 0, 0, 0, 0, 0, 0, t674, t675, -t677, t661 * t725 - t664 * t723; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t697, t698, 0, t703 * t725 + t704 * t723, 0, 0, 0, 0, 0, 0, t672, -t673, t676, t669 * t723 + t688 * t725, 0, 0, 0, 0, 0, 0, t672, t676, t673, t663 * t723 + t667 * t725, 0, 0, 0, 0, 0, 0, t672, t673, -t676, t661 * t723 + t664 * t725; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t678, -t679, 0, t668, 0, 0, 0, 0, 0, 0, t678, 0, t679, t662, 0, 0, 0, 0, 0, 0, t678, t679, 0, t660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t726, -qJDD(1), 0, t704, 0, 0, 0, 0, 0, 0, t681, -t682, t696, t669, 0, 0, 0, 0, 0, 0, t681, t696, t682, t663, 0, 0, 0, 0, 0, 0, t681, t682, -t696, t661; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t726, 0, t703, 0, 0, 0, 0, 0, 0, t695, -t693, t699, t688, 0, 0, 0, 0, 0, 0, t695, t699, t693, t667, 0, 0, 0, 0, 0, 0, t695, t693, -t699, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t678, -t679, 0, t668, 0, 0, 0, 0, 0, 0, t678, 0, t679, t662, 0, 0, 0, 0, 0, 0, t678, t679, 0, t660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, -t702, t734, t685, 0, 0, 0, 0, 0, 0, t707, t734, t702, t670, 0, 0, 0, 0, 0, 0, t707, t702, -t734, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t751, -t706, -t735, t684, 0, 0, 0, 0, 0, 0, -t751, -t735, t706, -t671, 0, 0, 0, 0, 0, 0, -t751, t706, t735, -t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t695, t693, -t699, -t688, 0, 0, 0, 0, 0, 0, -t695, -t699, -t693, -t667, 0, 0, 0, 0, 0, 0, -t695, -t693, t699, -t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, t734, t702, t670, 0, 0, 0, 0, 0, 0, t707, t702, -t734, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t695, -t699, -t693, -t667, 0, 0, 0, 0, 0, 0, -t695, -t693, t699, -t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t751, t735, -t706, t671, 0, 0, 0, 0, 0, 0, t751, -t706, -t735, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, t702, -t734, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t751, -t706, -t735, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t695, t693, -t699, t664;];
f_new_reg = t1;
