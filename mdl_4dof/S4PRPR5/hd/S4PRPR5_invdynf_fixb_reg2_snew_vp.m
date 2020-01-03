% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:26
% EndTime: 2019-12-31 16:23:27
% DurationCPUTime: 0.69s
% Computational Cost: add. (1319->96), mult. (2364->138), div. (0->0), fcn. (1646->8), ass. (0->76)
t611 = sin(pkin(6));
t613 = cos(pkin(6));
t594 = t611 * g(1) - t613 * g(2);
t627 = t611 * t594;
t595 = -t613 * g(1) - t611 * g(2);
t608 = -g(3) + qJDD(1);
t616 = sin(qJ(2));
t618 = cos(qJ(2));
t581 = t618 * t595 + t616 * t608;
t620 = qJD(2) ^ 2;
t579 = -t620 * pkin(2) + t581;
t610 = sin(pkin(7));
t612 = cos(pkin(7));
t580 = -t616 * t595 + t618 * t608;
t621 = qJDD(2) * pkin(2) + t580;
t559 = t612 * t579 + t610 * t621;
t615 = sin(qJ(4));
t606 = t615 ^ 2;
t617 = cos(qJ(4));
t607 = t617 ^ 2;
t626 = t606 + t607;
t625 = qJD(2) * qJD(4);
t624 = t615 * qJDD(2);
t623 = t617 * qJDD(2);
t558 = -t610 * t579 + t612 * t621;
t588 = t612 * qJDD(2) - t610 * t620;
t589 = -t610 * qJDD(2) - t612 * t620;
t622 = -t616 * t588 + t618 * t589;
t567 = t618 * t588 + t616 * t589;
t619 = qJD(4) ^ 2;
t600 = t615 * t620 * t617;
t599 = -t607 * t620 - t619;
t598 = -t606 * t620 - t619;
t597 = -qJDD(4) + t600;
t596 = qJDD(4) + t600;
t593 = t626 * t620;
t592 = t618 * qJDD(2) - t616 * t620;
t591 = -t616 * qJDD(2) - t618 * t620;
t590 = t626 * qJDD(2);
t587 = -qJDD(3) + t594;
t586 = -0.2e1 * t615 * t625 + t623;
t585 = 0.2e1 * t617 * t625 + t624;
t583 = t613 * t594;
t578 = t617 * t597 - t615 * t598;
t577 = -t615 * t596 + t617 * t599;
t576 = t615 * t597 + t617 * t598;
t575 = t617 * t596 + t615 * t599;
t571 = t612 * t590 - t610 * t593;
t570 = t610 * t590 + t612 * t593;
t565 = t612 * t578 + t610 * t585;
t564 = t612 * t577 - t610 * t586;
t563 = t610 * t578 - t612 * t585;
t562 = t610 * t577 + t612 * t586;
t561 = -t616 * t580 + t618 * t581;
t560 = t618 * t580 + t616 * t581;
t557 = -t620 * pkin(3) + qJDD(2) * pkin(5) + t559;
t556 = -qJDD(2) * pkin(3) - t620 * pkin(5) - t558;
t555 = -t616 * t570 + t618 * t571;
t554 = t618 * t570 + t616 * t571;
t553 = t617 * t557 - t615 * t587;
t552 = -t615 * t557 - t617 * t587;
t551 = -t616 * t563 + t618 * t565;
t550 = -t616 * t562 + t618 * t564;
t549 = t618 * t563 + t616 * t565;
t548 = t618 * t562 + t616 * t564;
t547 = -t610 * t558 + t612 * t559;
t546 = t612 * t558 + t610 * t559;
t545 = -t615 * t552 + t617 * t553;
t544 = t617 * t552 + t615 * t553;
t543 = t612 * t545 + t610 * t556;
t542 = t610 * t545 - t612 * t556;
t541 = -t616 * t546 + t618 * t547;
t540 = t618 * t546 + t616 * t547;
t539 = -t616 * t542 + t618 * t543;
t538 = t618 * t542 + t616 * t543;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t613 * t595 - t627, 0, 0, 0, 0, 0, 0, t613 * t591, -t613 * t592, 0, t613 * t561 - t627, 0, 0, 0, 0, 0, 0, t613 * t622, -t613 * t567, 0, t613 * t541 - t611 * t587, 0, 0, 0, 0, 0, 0, t613 * t550 + t611 * t575, t613 * t551 + t611 * t576, t613 * t555, t613 * t539 + t611 * t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t611 * t595 + t583, 0, 0, 0, 0, 0, 0, t611 * t591, -t611 * t592, 0, t611 * t561 + t583, 0, 0, 0, 0, 0, 0, t611 * t622, -t611 * t567, 0, t611 * t541 + t613 * t587, 0, 0, 0, 0, 0, 0, t611 * t550 - t613 * t575, t611 * t551 - t613 * t576, t611 * t555, t611 * t539 - t613 * t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t608, 0, 0, 0, 0, 0, 0, t592, t591, 0, t560, 0, 0, 0, 0, 0, 0, t567, t622, 0, t540, 0, 0, 0, 0, 0, 0, t548, t549, t554, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t595, 0, 0, 0, 0, 0, 0, t591, -t592, 0, t561, 0, 0, 0, 0, 0, 0, t622, -t567, 0, t541, 0, 0, 0, 0, 0, 0, t550, t551, t555, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t594, 0, 0, 0, 0, 0, 0, 0, 0, 0, t594, 0, 0, 0, 0, 0, 0, 0, 0, 0, t587, 0, 0, 0, 0, 0, 0, -t575, -t576, 0, -t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t608, 0, 0, 0, 0, 0, 0, t592, t591, 0, t560, 0, 0, 0, 0, 0, 0, t567, t622, 0, t540, 0, 0, 0, 0, 0, 0, t548, t549, t554, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, -qJDD(2), 0, t581, 0, 0, 0, 0, 0, 0, t589, -t588, 0, t547, 0, 0, 0, 0, 0, 0, t564, t565, t571, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t620, 0, t580, 0, 0, 0, 0, 0, 0, t588, t589, 0, t546, 0, 0, 0, 0, 0, 0, t562, t563, t570, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t594, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, 0, 0, 0, 0, 0, 0, t575, t576, 0, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, -qJDD(2), 0, t559, 0, 0, 0, 0, 0, 0, t577, t578, t590, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t620, 0, t558, 0, 0, 0, 0, 0, 0, t586, -t585, t593, -t556; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, 0, 0, 0, 0, 0, 0, t575, t576, 0, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599, t597, t623, t553; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, t598, -t624, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t586, t585, -t593, t556;];
f_new_reg = t1;
