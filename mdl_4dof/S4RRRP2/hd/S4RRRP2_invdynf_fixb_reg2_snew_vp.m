% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:17
% EndTime: 2019-12-31 17:13:18
% DurationCPUTime: 0.83s
% Computational Cost: add. (1635->106), mult. (2413->131), div. (0->0), fcn. (1325->6), ass. (0->79)
t728 = qJD(1) + qJD(2);
t726 = t728 ^ 2;
t727 = qJDD(1) + qJDD(2);
t732 = sin(qJ(2));
t735 = cos(qJ(2));
t707 = t732 * t726 - t735 * t727;
t733 = sin(qJ(1));
t736 = cos(qJ(1));
t741 = -t735 * t726 - t732 * t727;
t754 = t733 * t707 + t736 * t741;
t753 = t736 * t707 - t733 * t741;
t734 = cos(qJ(3));
t750 = t734 * g(3);
t749 = t726 * t734;
t731 = sin(qJ(3));
t748 = t728 * t731;
t730 = t734 ^ 2;
t747 = t730 * t726;
t746 = t731 * t727;
t745 = t734 * t727;
t720 = -t736 * g(1) - t733 * g(2);
t738 = qJD(1) ^ 2;
t710 = -t738 * pkin(1) + t720;
t719 = t733 * g(1) - t736 * g(2);
t739 = qJDD(1) * pkin(1) + t719;
t693 = t735 * t710 + t732 * t739;
t729 = t731 ^ 2;
t744 = t729 + t730;
t743 = 0.2e1 * t728 * t734;
t742 = qJD(3) * t748;
t689 = -t726 * pkin(2) + t727 * pkin(6) + t693;
t686 = -t731 * g(3) + t734 * t689;
t692 = -t732 * t710 + t735 * t739;
t688 = -t727 * pkin(2) - t726 * pkin(6) - t692;
t740 = -t742 + t745;
t737 = qJD(3) ^ 2;
t718 = t731 * t749;
t717 = -t737 - t747;
t716 = -t729 * t726 - t737;
t715 = -t733 * qJDD(1) - t736 * t738;
t714 = t736 * qJDD(1) - t733 * t738;
t713 = -qJDD(3) + t718;
t712 = qJDD(3) + t718;
t711 = qJD(3) * pkin(3) - qJ(4) * t748;
t709 = t744 * t726;
t704 = t744 * t727;
t700 = -0.2e1 * t742 + t745;
t699 = qJD(3) * t743 + t746;
t697 = t734 * t713 - t731 * t716;
t696 = -t731 * t712 + t734 * t717;
t695 = t731 * t713 + t734 * t716;
t694 = t734 * t712 + t731 * t717;
t691 = t735 * t704 - t732 * t709;
t690 = t732 * t704 + t735 * t709;
t685 = -t731 * t689 - t750;
t684 = t735 * t697 + t732 * t699;
t683 = t735 * t696 - t732 * t700;
t682 = t732 * t697 - t735 * t699;
t681 = t732 * t696 + t735 * t700;
t680 = -t732 * t692 + t735 * t693;
t679 = t735 * t692 + t732 * t693;
t678 = -t733 * t690 + t736 * t691;
t677 = t736 * t690 + t733 * t691;
t676 = -t740 * pkin(3) - qJ(4) * t747 + t711 * t748 + qJDD(4) + t688;
t675 = -pkin(3) * t747 + t740 * qJ(4) - qJD(3) * t711 + qJD(4) * t743 + t686;
t674 = qJDD(3) * pkin(3) - t750 + (pkin(3) * t749 - t727 * qJ(4) - 0.2e1 * qJD(4) * t728 - t689) * t731;
t673 = -t731 * t685 + t734 * t686;
t672 = t734 * t685 + t731 * t686;
t671 = -t733 * t682 + t736 * t684;
t670 = -t733 * t681 + t736 * t683;
t669 = t736 * t682 + t733 * t684;
t668 = t736 * t681 + t733 * t683;
t667 = t735 * t673 + t732 * t688;
t666 = t732 * t673 - t735 * t688;
t665 = -t731 * t674 + t734 * t675;
t664 = t734 * t674 + t731 * t675;
t663 = t735 * t665 + t732 * t676;
t662 = t732 * t665 - t735 * t676;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t715, -t714, 0, -t733 * t719 + t736 * t720, 0, 0, 0, 0, 0, 0, t754, t753, 0, -t733 * t679 + t736 * t680, 0, 0, 0, 0, 0, 0, t670, t671, t678, -t733 * t666 + t736 * t667, 0, 0, 0, 0, 0, 0, t670, t671, t678, -t733 * t662 + t736 * t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t714, t715, 0, t736 * t719 + t733 * t720, 0, 0, 0, 0, 0, 0, -t753, t754, 0, t736 * t679 + t733 * t680, 0, 0, 0, 0, 0, 0, t668, t669, t677, t736 * t666 + t733 * t667, 0, 0, 0, 0, 0, 0, t668, t669, t677, t736 * t662 + t733 * t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t694, t695, 0, t672, 0, 0, 0, 0, 0, 0, t694, t695, 0, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t738, -qJDD(1), 0, t720, 0, 0, 0, 0, 0, 0, t741, t707, 0, t680, 0, 0, 0, 0, 0, 0, t683, t684, t691, t667, 0, 0, 0, 0, 0, 0, t683, t684, t691, t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t738, 0, t719, 0, 0, 0, 0, 0, 0, -t707, t741, 0, t679, 0, 0, 0, 0, 0, 0, t681, t682, t690, t666, 0, 0, 0, 0, 0, 0, t681, t682, t690, t662; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t694, t695, 0, t672, 0, 0, 0, 0, 0, 0, t694, t695, 0, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t726, -t727, 0, t693, 0, 0, 0, 0, 0, 0, t696, t697, t704, t673, 0, 0, 0, 0, 0, 0, t696, t697, t704, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t727, -t726, 0, t692, 0, 0, 0, 0, 0, 0, t700, -t699, t709, -t688, 0, 0, 0, 0, 0, 0, t700, -t699, t709, -t676; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t694, t695, 0, t672, 0, 0, 0, 0, 0, 0, t694, t695, 0, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t717, t713, t745, t686, 0, 0, 0, 0, 0, 0, t717, t713, t745, t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t712, t716, -t746, t685, 0, 0, 0, 0, 0, 0, t712, t716, -t746, t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t700, t699, -t709, t688, 0, 0, 0, 0, 0, 0, -t700, t699, -t709, t676; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t717, t713, t745, t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t712, t716, -t746, t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t700, t699, -t709, t676;];
f_new_reg = t1;
