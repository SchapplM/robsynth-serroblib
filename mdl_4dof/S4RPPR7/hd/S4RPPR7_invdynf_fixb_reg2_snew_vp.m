% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR7_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:49
% EndTime: 2019-12-31 16:41:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (1396->123), mult. (3109->138), div. (0->0), fcn. (1956->6), ass. (0->77)
t617 = sin(qJ(1));
t619 = cos(qJ(1));
t598 = t617 * g(1) - t619 * g(2);
t621 = qJD(1) ^ 2;
t626 = -t621 * qJ(2) + qJDD(2) - t598;
t636 = -qJ(3) - pkin(1);
t638 = -2 * qJD(1);
t643 = (qJD(3) * t638) + t636 * qJDD(1) + t626;
t614 = sin(pkin(6));
t611 = t614 ^ 2;
t615 = cos(pkin(6));
t612 = t615 ^ 2;
t632 = t611 + t612;
t642 = t632 * t621;
t616 = sin(qJ(4));
t618 = cos(qJ(4));
t627 = t614 * t618 + t615 * t616;
t569 = t627 * qJDD(1);
t589 = t627 * qJD(1);
t640 = t589 ^ 2;
t591 = (-t614 * t616 + t615 * t618) * qJD(1);
t639 = t591 ^ 2;
t637 = 2 * qJD(4);
t635 = t591 * t589;
t634 = t611 * t621;
t633 = t614 * t621;
t605 = t614 * qJDD(1);
t606 = t615 * qJDD(1);
t631 = t617 * qJDD(1);
t630 = t619 * qJDD(1);
t628 = t615 * t633;
t574 = t614 * g(3) + t643 * t615;
t575 = -t615 * g(3) + t643 * t614;
t599 = -t619 * g(1) - t617 * g(2);
t588 = t616 * t605 - t618 * t606;
t625 = -qJDD(1) * qJ(2) + (qJD(2) * t638) - t599;
t624 = -qJDD(3) + t625;
t620 = qJD(4) ^ 2;
t597 = t619 * t621 + t631;
t596 = -t617 * t621 + t630;
t594 = t632 * qJDD(1);
t593 = t614 * t642;
t592 = t615 * t642;
t587 = qJDD(1) * pkin(1) - t626;
t586 = t621 * pkin(1) + t625;
t585 = -t620 - t639;
t583 = -t636 * t621 + t624;
t580 = -t589 * t637 - t588;
t579 = t591 * t637 + t569;
t578 = -qJDD(4) - t635;
t577 = qJDD(4) - t635;
t576 = -t620 - t640;
t573 = -pkin(3) * t605 + (t632 * pkin(5) - t636) * t621 + t624;
t572 = -pkin(3) * t634 - pkin(5) * t605 + t575;
t571 = -t639 - t640;
t570 = (-pkin(3) * t633 - pkin(5) * qJDD(1)) * t615 + t574;
t568 = t618 * t578 - t616 * t585;
t567 = t616 * t578 + t618 * t585;
t566 = -t618 * t569 - t616 * t588;
t565 = -t616 * t569 + t618 * t588;
t564 = t618 * t576 - t616 * t577;
t563 = t616 * t576 + t618 * t577;
t562 = -t614 * t574 + t615 * t575;
t561 = t615 * t574 + t614 * t575;
t560 = t616 * t570 + t618 * t572;
t559 = t618 * t570 - t616 * t572;
t558 = -t614 * t567 + t615 * t568;
t557 = t615 * t567 + t614 * t568;
t556 = -t614 * t565 + t615 * t566;
t555 = t615 * t565 + t614 * t566;
t554 = -t614 * t563 + t615 * t564;
t553 = t615 * t563 + t614 * t564;
t552 = -t616 * t559 + t618 * t560;
t551 = t618 * t559 + t616 * t560;
t550 = -t614 * t551 + t615 * t552;
t549 = t615 * t551 + t614 * t552;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t597, -t596, 0, -t617 * t598 + t619 * t599, 0, 0, 0, 0, 0, 0, 0, t597, t596, -t619 * t586 - t617 * t587, 0, 0, 0, 0, 0, 0, -t617 * t593 + t614 * t630, -t617 * t592 + t615 * t630, -t617 * t594 - t619 * t642, t617 * t561 - t619 * t583, 0, 0, 0, 0, 0, 0, t617 * t553 + t619 * t579, t617 * t557 + t619 * t580, t617 * t555 + t619 * t571, t617 * t549 - t619 * t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t596, -t597, 0, t619 * t598 + t617 * t599, 0, 0, 0, 0, 0, 0, 0, -t596, t597, -t617 * t586 + t619 * t587, 0, 0, 0, 0, 0, 0, t619 * t593 + t614 * t631, t619 * t592 + t615 * t631, t619 * t594 - t617 * t642, -t619 * t561 - t617 * t583, 0, 0, 0, 0, 0, 0, -t619 * t553 + t617 * t579, -t619 * t557 + t617 * t580, -t619 * t555 + t617 * t571, -t619 * t549 - t617 * t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, 0, 0, 0, 0, 0, 0, t554, t558, t556, t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t621, -qJDD(1), 0, t599, 0, 0, 0, 0, 0, 0, 0, t621, qJDD(1), -t586, 0, 0, 0, 0, 0, 0, t605, t606, -t642, -t583, 0, 0, 0, 0, 0, 0, t579, t580, t571, -t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t621, 0, t598, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t621, t587, 0, 0, 0, 0, 0, 0, t593, t592, t594, -t561, 0, 0, 0, 0, 0, 0, -t553, -t557, -t555, -t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, 0, 0, 0, 0, 0, 0, t554, t558, t556, t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, 0, 0, 0, 0, 0, 0, t554, t558, t556, t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t621, -qJDD(1), t586, 0, 0, 0, 0, 0, 0, -t605, -t606, t642, t583, 0, 0, 0, 0, 0, 0, -t579, -t580, -t571, t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t621, -t587, 0, 0, 0, 0, 0, 0, -t593, -t592, -t594, t561, 0, 0, 0, 0, 0, 0, t553, t557, t555, t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t634, -t628, -t605, t575, 0, 0, 0, 0, 0, 0, t564, t568, t566, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t628, -t612 * t621, -t606, t574, 0, 0, 0, 0, 0, 0, t563, t567, t565, t551; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, t606, -t642, -t583, 0, 0, 0, 0, 0, 0, t579, t580, t571, -t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t576, t578, -t569, t560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t585, t588, t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, t580, t571, -t573;];
f_new_reg = t1;
