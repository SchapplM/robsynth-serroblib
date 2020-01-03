% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:45
% EndTime: 2019-12-31 17:00:46
% DurationCPUTime: 0.84s
% Computational Cost: add. (804->140), mult. (1880->113), div. (0->0), fcn. (881->4), ass. (0->67)
t734 = sin(qJ(2));
t736 = cos(qJ(2));
t738 = qJD(1) ^ 2;
t745 = t736 * t738 * t734;
t710 = qJDD(2) + t745;
t731 = t736 ^ 2;
t756 = t731 * t738;
t760 = qJD(2) ^ 2;
t716 = -t756 - t760;
t692 = t734 * t710 - t736 * t716;
t748 = qJD(1) * qJD(2);
t744 = t734 * t748;
t746 = t736 * qJDD(1);
t704 = -0.2e1 * t744 + t746;
t735 = sin(qJ(1));
t737 = cos(qJ(1));
t762 = t735 * t692 - t737 * t704;
t761 = t737 * t692 + t735 * t704;
t711 = qJDD(2) - t745;
t730 = t734 ^ 2;
t715 = -t730 * t738 - t760;
t693 = t736 * t711 + t734 * t715;
t723 = t736 * t748;
t747 = t734 * qJDD(1);
t702 = 0.2e1 * t723 + t747;
t682 = t735 * t693 + t737 * t702;
t684 = t737 * t693 - t735 * t702;
t759 = 2 * qJD(3);
t758 = t736 * g(3);
t757 = t738 * pkin(5);
t751 = t730 + t731;
t750 = qJD(1) * t734;
t749 = pkin(3) * t750 - qJD(2) * qJ(4) + t759;
t712 = t735 * g(1) - t737 * g(2);
t713 = -t737 * g(1) - t735 * g(2);
t698 = -t738 * pkin(1) + qJDD(1) * pkin(5) + t713;
t743 = t738 * (-pkin(2) * t736 - t734 * qJ(3)) + t698;
t742 = t723 + t747;
t687 = t736 * t710 + t734 * t716;
t689 = t734 * t711 - t736 * t715;
t741 = qJDD(1) * pkin(1) + t712;
t703 = -t744 + t746;
t740 = t741 + (t703 - t744) * pkin(2);
t726 = t734 * g(3);
t739 = -t760 * pkin(2) + qJDD(2) * qJ(3) + t743 * t736 - t726;
t680 = -qJDD(2) * pkin(2) - t760 * qJ(3) + t743 * t734 + qJDD(3) + t758;
t708 = t751 * t738;
t707 = -t735 * qJDD(1) - t737 * t738;
t706 = t737 * qJDD(1) - t735 * t738;
t705 = t751 * qJDD(1);
t697 = t741 + t757;
t696 = t736 * t698 - t726;
t695 = -t734 * t698 - t758;
t686 = t737 * t705 - t735 * t708;
t685 = t735 * t705 + t737 * t708;
t679 = qJD(2) * t759 + t739;
t678 = -t734 * t695 + t736 * t696;
t677 = t736 * t695 + t734 * t696;
t676 = t750 * t759 + t757 + (t742 + t723) * qJ(3) + t740;
t675 = t703 * pkin(3) - qJ(4) * t756 + t749 * qJD(2) + qJDD(4) + t739;
t674 = -0.2e1 * qJD(4) * qJD(2) + t680 - t710 * qJ(4) + (-t723 + t742) * pkin(3);
t673 = t742 * qJ(3) + t703 * qJ(4) + (t731 * pkin(3) + pkin(5)) * t738 + ((qJ(3) * qJD(2) + 0.2e1 * qJD(4)) * t736 + t749 * t734) * qJD(1) + t740;
t672 = t736 * t679 + t734 * t680;
t671 = t734 * t679 - t736 * t680;
t670 = t734 * t674 + t736 * t675;
t669 = -t736 * t674 + t734 * t675;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t707, -t706, 0, -t735 * t712 + t737 * t713, 0, 0, 0, 0, 0, 0, -t761, -t684, t686, t737 * t678 - t735 * t697, 0, 0, 0, 0, 0, 0, t686, t761, t684, t737 * t672 - t735 * t676, 0, 0, 0, 0, 0, 0, t686, t684, -t761, t737 * t670 - t735 * t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t706, t707, 0, t737 * t712 + t735 * t713, 0, 0, 0, 0, 0, 0, -t762, -t682, t685, t735 * t678 + t737 * t697, 0, 0, 0, 0, 0, 0, t685, t762, t682, t735 * t672 + t737 * t676, 0, 0, 0, 0, 0, 0, t685, t682, -t762, t735 * t670 + t737 * t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t687, -t689, 0, t677, 0, 0, 0, 0, 0, 0, 0, -t687, t689, t671, 0, 0, 0, 0, 0, 0, 0, t689, t687, t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t738, -qJDD(1), 0, t713, 0, 0, 0, 0, 0, 0, -t692, -t693, t705, t678, 0, 0, 0, 0, 0, 0, t705, t692, t693, t672, 0, 0, 0, 0, 0, 0, t705, t693, -t692, t670; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t738, 0, t712, 0, 0, 0, 0, 0, 0, t704, -t702, t708, t697, 0, 0, 0, 0, 0, 0, t708, -t704, t702, t676, 0, 0, 0, 0, 0, 0, t708, t702, t704, t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t687, -t689, 0, t677, 0, 0, 0, 0, 0, 0, 0, -t687, t689, t671, 0, 0, 0, 0, 0, 0, 0, t689, t687, t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, -t711, t746, t696, 0, 0, 0, 0, 0, 0, t746, -t716, t711, t679, 0, 0, 0, 0, 0, 0, t746, t711, t716, t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, t715, -t747, t695, 0, 0, 0, 0, 0, 0, -t747, -t710, -t715, -t680, 0, 0, 0, 0, 0, 0, -t747, -t715, t710, -t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t704, t702, -t708, -t697, 0, 0, 0, 0, 0, 0, -t708, t704, -t702, -t676, 0, 0, 0, 0, 0, 0, -t708, -t702, -t704, -t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t708, t704, -t702, -t676, 0, 0, 0, 0, 0, 0, -t708, -t702, -t704, -t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t746, t716, -t711, -t679, 0, 0, 0, 0, 0, 0, -t746, -t711, -t716, -t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t747, t710, t715, t680, 0, 0, 0, 0, 0, 0, t747, t715, -t710, t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t708, -t702, -t704, -t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t747, t715, -t710, t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t746, t711, t716, t675;];
f_new_reg = t1;
