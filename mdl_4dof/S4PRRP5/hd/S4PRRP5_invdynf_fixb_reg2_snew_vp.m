% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:28
% EndTime: 2019-12-31 16:29:29
% DurationCPUTime: 0.64s
% Computational Cost: add. (926->96), mult. (2011->120), div. (0->0), fcn. (1196->6), ass. (0->73)
t652 = cos(qJ(3));
t645 = t652 ^ 2;
t655 = qJD(2) ^ 2;
t666 = t645 * t655;
t647 = sin(pkin(6));
t648 = cos(pkin(6));
t629 = t647 * g(1) - t648 * g(2);
t665 = t647 * t629;
t664 = t652 * t655;
t663 = -g(3) + qJDD(1);
t630 = -t648 * g(1) - t647 * g(2);
t651 = sin(qJ(2));
t653 = cos(qJ(2));
t617 = t653 * t630 + t651 * t663;
t611 = -t655 * pkin(2) + qJDD(2) * pkin(5) + t617;
t650 = sin(qJ(3));
t602 = t652 * t611 - t650 * t629;
t644 = t650 ^ 2;
t662 = t644 + t645;
t661 = qJD(2) * t650;
t660 = t650 * qJDD(2);
t659 = t652 * qJDD(2);
t658 = 0.2e1 * qJD(2) * t652;
t657 = qJD(3) * t661;
t616 = -t651 * t630 + t653 * t663;
t610 = -qJDD(2) * pkin(2) - t655 * pkin(5) - t616;
t656 = -t657 + t659;
t654 = qJD(3) ^ 2;
t636 = t650 * t664;
t635 = -t654 - t666;
t634 = -t644 * t655 - t654;
t633 = -qJDD(3) + t636;
t632 = qJDD(3) + t636;
t631 = qJD(3) * pkin(3) - qJ(4) * t661;
t628 = t662 * t655;
t627 = t653 * qJDD(2) - t651 * t655;
t626 = -t651 * qJDD(2) - t653 * t655;
t625 = t662 * qJDD(2);
t624 = -0.2e1 * t657 + t659;
t623 = qJD(3) * t658 + t660;
t621 = t652 * t629;
t618 = t648 * t629;
t615 = t652 * t633 - t650 * t634;
t614 = -t650 * t632 + t652 * t635;
t613 = t650 * t633 + t652 * t634;
t612 = t652 * t632 + t650 * t635;
t609 = t653 * t625 - t651 * t628;
t608 = t651 * t625 + t653 * t628;
t606 = t648 * t609;
t605 = t647 * t609;
t604 = t653 * t615 + t651 * t623;
t603 = t653 * t614 - t651 * t624;
t601 = t651 * t615 - t653 * t623;
t600 = t651 * t614 + t653 * t624;
t599 = -t650 * t611 - t621;
t598 = -t651 * t616 + t653 * t617;
t597 = t653 * t616 + t651 * t617;
t596 = -t656 * pkin(3) - qJ(4) * t666 + t631 * t661 + qJDD(4) + t610;
t595 = t648 * t604 + t647 * t613;
t594 = t648 * t603 + t647 * t612;
t593 = t647 * t604 - t648 * t613;
t592 = t647 * t603 - t648 * t612;
t591 = -pkin(3) * t666 + t656 * qJ(4) - qJD(3) * t631 + qJD(4) * t658 + t602;
t590 = qJDD(3) * pkin(3) - t621 + (pkin(3) * t664 - qJDD(2) * qJ(4) - 0.2e1 * qJD(2) * qJD(4) - t611) * t650;
t589 = -t650 * t599 + t652 * t602;
t588 = t652 * t599 + t650 * t602;
t587 = t653 * t589 + t651 * t610;
t586 = t651 * t589 - t653 * t610;
t585 = -t650 * t590 + t652 * t591;
t584 = t652 * t590 + t650 * t591;
t583 = t653 * t585 + t651 * t596;
t582 = t651 * t585 - t653 * t596;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t648 * t630 - t665, 0, 0, 0, 0, 0, 0, t648 * t626, -t648 * t627, 0, t648 * t598 - t665, 0, 0, 0, 0, 0, 0, t594, t595, t606, t648 * t587 + t647 * t588, 0, 0, 0, 0, 0, 0, t594, t595, t606, t648 * t583 + t647 * t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t647 * t630 + t618, 0, 0, 0, 0, 0, 0, t647 * t626, -t647 * t627, 0, t647 * t598 + t618, 0, 0, 0, 0, 0, 0, t592, t593, t605, t647 * t587 - t648 * t588, 0, 0, 0, 0, 0, 0, t592, t593, t605, t647 * t583 - t648 * t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t663, 0, 0, 0, 0, 0, 0, t627, t626, 0, t597, 0, 0, 0, 0, 0, 0, t600, t601, t608, t586, 0, 0, 0, 0, 0, 0, t600, t601, t608, t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t630, 0, 0, 0, 0, 0, 0, t626, -t627, 0, t598, 0, 0, 0, 0, 0, 0, t603, t604, t609, t587, 0, 0, 0, 0, 0, 0, t603, t604, t609, t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t629, 0, 0, 0, 0, 0, 0, 0, 0, 0, t629, 0, 0, 0, 0, 0, 0, -t612, -t613, 0, -t588, 0, 0, 0, 0, 0, 0, -t612, -t613, 0, -t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t663, 0, 0, 0, 0, 0, 0, t627, t626, 0, t597, 0, 0, 0, 0, 0, 0, t600, t601, t608, t586, 0, 0, 0, 0, 0, 0, t600, t601, t608, t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t655, -qJDD(2), 0, t617, 0, 0, 0, 0, 0, 0, t614, t615, t625, t589, 0, 0, 0, 0, 0, 0, t614, t615, t625, t585; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t655, 0, t616, 0, 0, 0, 0, 0, 0, t624, -t623, t628, -t610, 0, 0, 0, 0, 0, 0, t624, -t623, t628, -t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t629, 0, 0, 0, 0, 0, 0, t612, t613, 0, t588, 0, 0, 0, 0, 0, 0, t612, t613, 0, t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, t633, t659, t602, 0, 0, 0, 0, 0, 0, t635, t633, t659, t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t632, t634, -t660, t599, 0, 0, 0, 0, 0, 0, t632, t634, -t660, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t624, t623, -t628, t610, 0, 0, 0, 0, 0, 0, -t624, t623, -t628, t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, t633, t659, t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t632, t634, -t660, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t624, t623, -t628, t596;];
f_new_reg = t1;
