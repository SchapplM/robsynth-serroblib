% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:02
% EndTime: 2019-12-31 16:27:03
% DurationCPUTime: 0.64s
% Computational Cost: add. (974->92), mult. (2168->124), div. (0->0), fcn. (1317->6), ass. (0->73)
t694 = cos(qJ(3));
t686 = t694 ^ 2;
t697 = qJD(2) ^ 2;
t709 = t686 * t697;
t708 = t694 * t697;
t689 = sin(pkin(6));
t690 = cos(pkin(6));
t672 = -t690 * g(1) - t689 * g(2);
t693 = sin(qJ(2));
t695 = cos(qJ(2));
t700 = t689 * g(1) - t690 * g(2);
t656 = t695 * t672 + t693 * t700;
t652 = -t697 * pkin(2) + qJDD(2) * pkin(5) + t656;
t687 = -g(3) + qJDD(1);
t692 = sin(qJ(3));
t647 = t694 * t652 + t692 * t687;
t685 = t692 ^ 2;
t707 = t685 + t686;
t706 = qJD(2) * t692;
t705 = t692 * qJDD(2);
t704 = t694 * qJDD(2);
t703 = 0.2e1 * qJD(2) * t694;
t702 = qJD(3) * t706;
t669 = t695 * qJDD(2) - t693 * t697;
t670 = -t693 * qJDD(2) - t695 * t697;
t701 = -t689 * t669 + t690 * t670;
t655 = -t693 * t672 + t695 * t700;
t699 = t690 * t669 + t689 * t670;
t651 = -qJDD(2) * pkin(2) - t697 * pkin(5) - t655;
t698 = -t702 + t704;
t696 = qJD(3) ^ 2;
t682 = t694 * t687;
t678 = t692 * t708;
t677 = -t696 - t709;
t676 = -t685 * t697 - t696;
t675 = -qJDD(3) + t678;
t674 = qJDD(3) + t678;
t673 = qJD(3) * pkin(3) - qJ(4) * t706;
t671 = t707 * t697;
t668 = t707 * qJDD(2);
t667 = -0.2e1 * t702 + t704;
t666 = qJD(3) * t703 + t705;
t660 = t694 * t675 - t692 * t676;
t659 = -t692 * t674 + t694 * t677;
t658 = t692 * t675 + t694 * t676;
t657 = t694 * t674 + t692 * t677;
t654 = t695 * t668 - t693 * t671;
t653 = t693 * t668 + t695 * t671;
t649 = t695 * t660 + t693 * t666;
t648 = t695 * t659 - t693 * t667;
t646 = t693 * t660 - t695 * t666;
t645 = t693 * t659 + t695 * t667;
t644 = -t692 * t652 + t682;
t643 = -t693 * t655 + t695 * t656;
t642 = t695 * t655 + t693 * t656;
t641 = -t689 * t653 + t690 * t654;
t640 = t690 * t653 + t689 * t654;
t639 = -t698 * pkin(3) - qJ(4) * t709 + t673 * t706 + qJDD(4) + t651;
t638 = -pkin(3) * t709 + t698 * qJ(4) - qJD(3) * t673 + qJD(4) * t703 + t647;
t637 = qJDD(3) * pkin(3) + t682 + (pkin(3) * t708 - qJDD(2) * qJ(4) - 0.2e1 * qJD(2) * qJD(4) - t652) * t692;
t636 = -t692 * t644 + t694 * t647;
t635 = t694 * t644 + t692 * t647;
t634 = -t689 * t646 + t690 * t649;
t633 = -t689 * t645 + t690 * t648;
t632 = t690 * t646 + t689 * t649;
t631 = t690 * t645 + t689 * t648;
t630 = t695 * t636 + t693 * t651;
t629 = t693 * t636 - t695 * t651;
t628 = -t692 * t637 + t694 * t638;
t627 = t694 * t637 + t692 * t638;
t626 = t695 * t628 + t693 * t639;
t625 = t693 * t628 - t695 * t639;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t690 * t672 - t689 * t700, 0, 0, 0, 0, 0, 0, t701, -t699, 0, -t689 * t642 + t690 * t643, 0, 0, 0, 0, 0, 0, t633, t634, t641, -t689 * t629 + t690 * t630, 0, 0, 0, 0, 0, 0, t633, t634, t641, -t689 * t625 + t690 * t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t689 * t672 + t690 * t700, 0, 0, 0, 0, 0, 0, t699, t701, 0, t690 * t642 + t689 * t643, 0, 0, 0, 0, 0, 0, t631, t632, t640, t690 * t629 + t689 * t630, 0, 0, 0, 0, 0, 0, t631, t632, t640, t690 * t625 + t689 * t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, 0, 0, 0, 0, 0, 0, t657, t658, 0, t635, 0, 0, 0, 0, 0, 0, t657, t658, 0, t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t672, 0, 0, 0, 0, 0, 0, t670, -t669, 0, t643, 0, 0, 0, 0, 0, 0, t648, t649, t654, t630, 0, 0, 0, 0, 0, 0, t648, t649, t654, t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t700, 0, 0, 0, 0, 0, 0, t669, t670, 0, t642, 0, 0, 0, 0, 0, 0, t645, t646, t653, t629, 0, 0, 0, 0, 0, 0, t645, t646, t653, t625; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, 0, 0, 0, 0, 0, 0, t657, t658, 0, t635, 0, 0, 0, 0, 0, 0, t657, t658, 0, t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t697, -qJDD(2), 0, t656, 0, 0, 0, 0, 0, 0, t659, t660, t668, t636, 0, 0, 0, 0, 0, 0, t659, t660, t668, t628; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t697, 0, t655, 0, 0, 0, 0, 0, 0, t667, -t666, t671, -t651, 0, 0, 0, 0, 0, 0, t667, -t666, t671, -t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, 0, 0, 0, 0, 0, 0, t657, t658, 0, t635, 0, 0, 0, 0, 0, 0, t657, t658, 0, t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t677, t675, t704, t647, 0, 0, 0, 0, 0, 0, t677, t675, t704, t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t674, t676, -t705, t644, 0, 0, 0, 0, 0, 0, t674, t676, -t705, t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t667, t666, -t671, t651, 0, 0, 0, 0, 0, 0, -t667, t666, -t671, t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t677, t675, t704, t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t674, t676, -t705, t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t667, t666, -t671, t639;];
f_new_reg = t1;
