% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:54
% EndTime: 2019-12-31 16:42:55
% DurationCPUTime: 0.67s
% Computational Cost: add. (1135->102), mult. (2413->130), div. (0->0), fcn. (1325->6), ass. (0->77)
t712 = cos(qJ(3));
t704 = t712 ^ 2;
t715 = qJD(1) ^ 2;
t727 = t704 * t715;
t726 = t712 * t715;
t711 = sin(qJ(1));
t713 = cos(qJ(1));
t693 = -t713 * g(1) - t711 * g(2);
t680 = -t715 * pkin(1) + t693;
t707 = sin(pkin(6));
t708 = cos(pkin(6));
t692 = t711 * g(1) - t713 * g(2);
t716 = qJDD(1) * pkin(1) + t692;
t668 = t708 * t680 + t707 * t716;
t666 = -t715 * pkin(2) + qJDD(1) * pkin(5) + t668;
t705 = -g(3) + qJDD(2);
t710 = sin(qJ(3));
t659 = t712 * t666 + t710 * t705;
t703 = t710 ^ 2;
t725 = t703 + t704;
t724 = qJD(1) * t710;
t723 = t710 * qJDD(1);
t722 = t712 * qJDD(1);
t721 = 0.2e1 * qJD(1) * t712;
t720 = qJD(3) * t724;
t667 = -t707 * t680 + t708 * t716;
t683 = -t707 * qJDD(1) - t708 * t715;
t684 = t708 * qJDD(1) - t707 * t715;
t719 = t713 * t683 - t711 * t684;
t718 = t711 * t683 + t713 * t684;
t665 = -qJDD(1) * pkin(2) - t715 * pkin(5) - t667;
t717 = -t720 + t722;
t714 = qJD(3) ^ 2;
t700 = t712 * t705;
t696 = t710 * t726;
t695 = -t714 - t727;
t694 = -t703 * t715 - t714;
t691 = -qJDD(3) + t696;
t690 = qJDD(3) + t696;
t689 = qJD(3) * pkin(3) - qJ(4) * t724;
t688 = t725 * t715;
t687 = -t711 * qJDD(1) - t713 * t715;
t686 = t713 * qJDD(1) - t711 * t715;
t685 = t725 * qJDD(1);
t682 = -0.2e1 * t720 + t722;
t681 = qJD(3) * t721 + t723;
t674 = t712 * t691 - t710 * t694;
t673 = -t710 * t690 + t712 * t695;
t672 = t710 * t691 + t712 * t694;
t671 = t712 * t690 + t710 * t695;
t670 = t708 * t685 - t707 * t688;
t669 = t707 * t685 + t708 * t688;
t663 = t708 * t674 + t707 * t681;
t662 = t708 * t673 - t707 * t682;
t661 = t707 * t674 - t708 * t681;
t660 = t707 * t673 + t708 * t682;
t658 = -t710 * t666 + t700;
t657 = -t711 * t669 + t713 * t670;
t656 = t713 * t669 + t711 * t670;
t655 = -t707 * t667 + t708 * t668;
t654 = t708 * t667 + t707 * t668;
t653 = -t717 * pkin(3) - qJ(4) * t727 + t689 * t724 + qJDD(4) + t665;
t652 = -pkin(3) * t727 + t717 * qJ(4) - qJD(3) * t689 + qJD(4) * t721 + t659;
t651 = qJDD(3) * pkin(3) + t700 + (pkin(3) * t726 - qJDD(1) * qJ(4) - 0.2e1 * qJD(1) * qJD(4) - t666) * t710;
t650 = -t711 * t661 + t713 * t663;
t649 = -t711 * t660 + t713 * t662;
t648 = t713 * t661 + t711 * t663;
t647 = t713 * t660 + t711 * t662;
t646 = -t710 * t658 + t712 * t659;
t645 = t712 * t658 + t710 * t659;
t644 = t708 * t646 + t707 * t665;
t643 = t707 * t646 - t708 * t665;
t642 = -t710 * t651 + t712 * t652;
t641 = t712 * t651 + t710 * t652;
t640 = t708 * t642 + t707 * t653;
t639 = t707 * t642 - t708 * t653;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t687, -t686, 0, -t711 * t692 + t713 * t693, 0, 0, 0, 0, 0, 0, t719, -t718, 0, -t711 * t654 + t713 * t655, 0, 0, 0, 0, 0, 0, t649, t650, t657, -t711 * t643 + t713 * t644, 0, 0, 0, 0, 0, 0, t649, t650, t657, -t711 * t639 + t713 * t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t686, t687, 0, t713 * t692 + t711 * t693, 0, 0, 0, 0, 0, 0, t718, t719, 0, t713 * t654 + t711 * t655, 0, 0, 0, 0, 0, 0, t647, t648, t656, t713 * t643 + t711 * t644, 0, 0, 0, 0, 0, 0, t647, t648, t656, t713 * t639 + t711 * t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t705, 0, 0, 0, 0, 0, 0, t671, t672, 0, t645, 0, 0, 0, 0, 0, 0, t671, t672, 0, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t715, -qJDD(1), 0, t693, 0, 0, 0, 0, 0, 0, t683, -t684, 0, t655, 0, 0, 0, 0, 0, 0, t662, t663, t670, t644, 0, 0, 0, 0, 0, 0, t662, t663, t670, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t715, 0, t692, 0, 0, 0, 0, 0, 0, t684, t683, 0, t654, 0, 0, 0, 0, 0, 0, t660, t661, t669, t643, 0, 0, 0, 0, 0, 0, t660, t661, t669, t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t705, 0, 0, 0, 0, 0, 0, t671, t672, 0, t645, 0, 0, 0, 0, 0, 0, t671, t672, 0, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t715, -qJDD(1), 0, t668, 0, 0, 0, 0, 0, 0, t673, t674, t685, t646, 0, 0, 0, 0, 0, 0, t673, t674, t685, t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t715, 0, t667, 0, 0, 0, 0, 0, 0, t682, -t681, t688, -t665, 0, 0, 0, 0, 0, 0, t682, -t681, t688, -t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t705, 0, 0, 0, 0, 0, 0, t671, t672, 0, t645, 0, 0, 0, 0, 0, 0, t671, t672, 0, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t695, t691, t722, t659, 0, 0, 0, 0, 0, 0, t695, t691, t722, t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t690, t694, -t723, t658, 0, 0, 0, 0, 0, 0, t690, t694, -t723, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t682, t681, -t688, t665, 0, 0, 0, 0, 0, 0, -t682, t681, -t688, t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t695, t691, t722, t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t690, t694, -t723, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t682, t681, -t688, t653;];
f_new_reg = t1;
