% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:07
% EndTime: 2019-12-31 16:22:08
% DurationCPUTime: 0.62s
% Computational Cost: add. (667->89), mult. (1293->108), div. (0->0), fcn. (840->6), ass. (0->57)
t559 = sin(qJ(2));
t561 = cos(qJ(2));
t563 = qJD(2) ^ 2;
t539 = t561 * qJDD(2) - t559 * t563;
t540 = t559 * qJDD(2) + t561 * t563;
t555 = sin(pkin(6));
t556 = cos(pkin(6));
t565 = t555 * t539 + t556 * t540;
t567 = t556 * t539 - t555 * t540;
t558 = sin(qJ(4));
t551 = t558 ^ 2;
t560 = cos(qJ(4));
t552 = t560 ^ 2;
t572 = t551 + t552;
t571 = qJD(2) * qJD(4);
t570 = t558 * qJDD(2);
t569 = t560 * qJDD(2);
t568 = t560 * t563 * t558;
t542 = t555 * g(1) - t556 * g(2);
t543 = -t556 * g(1) - t555 * g(2);
t525 = t561 * t542 - t559 * t543;
t526 = t559 * t542 + t561 * t543;
t522 = -qJDD(2) * pkin(2) - t563 * qJ(3) + qJDD(3) - t525;
t521 = -t563 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t526;
t562 = qJD(4) ^ 2;
t553 = -g(3) + qJDD(1);
t547 = -t552 * t563 - t562;
t546 = -t551 * t563 - t562;
t545 = -qJDD(4) - t568;
t544 = qJDD(4) - t568;
t541 = t572 * t563;
t538 = t572 * qJDD(2);
t537 = -0.2e1 * t558 * t571 + t569;
t536 = 0.2e1 * t560 * t571 + t570;
t530 = t560 * t545 - t558 * t547;
t529 = -t558 * t544 + t560 * t546;
t528 = t558 * t545 + t560 * t547;
t527 = t560 * t544 + t558 * t546;
t524 = -t559 * t538 - t561 * t541;
t523 = t561 * t538 - t559 * t541;
t520 = -qJDD(2) * pkin(5) + t522;
t519 = -t563 * pkin(5) + t521;
t518 = t559 * t528 + t561 * t537;
t517 = t559 * t527 + t561 * t536;
t516 = -t561 * t528 + t559 * t537;
t515 = -t561 * t527 + t559 * t536;
t514 = t558 * t520 + t560 * t553;
t513 = t560 * t520 - t558 * t553;
t512 = -t559 * t525 + t561 * t526;
t511 = t561 * t525 + t559 * t526;
t510 = t561 * t521 + t559 * t522;
t509 = t559 * t521 - t561 * t522;
t508 = -t558 * t513 + t560 * t514;
t507 = t560 * t513 + t558 * t514;
t506 = t559 * t507 + t561 * t519;
t505 = -t561 * t507 + t559 * t519;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t555 * t542 + t556 * t543, 0, 0, 0, 0, 0, 0, -t565, -t567, 0, -t555 * t511 + t556 * t512, 0, 0, 0, 0, 0, 0, 0, t565, t567, -t555 * t509 + t556 * t510, 0, 0, 0, 0, 0, 0, -t555 * t515 + t556 * t517, -t555 * t516 + t556 * t518, -t555 * t523 + t556 * t524, -t555 * t505 + t556 * t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t556 * t542 + t555 * t543, 0, 0, 0, 0, 0, 0, t567, -t565, 0, t556 * t511 + t555 * t512, 0, 0, 0, 0, 0, 0, 0, -t567, t565, t556 * t509 + t555 * t510, 0, 0, 0, 0, 0, 0, t556 * t515 + t555 * t517, t556 * t516 + t555 * t518, t556 * t523 + t555 * t524, t556 * t505 + t555 * t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, t529, t530, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t543, 0, 0, 0, 0, 0, 0, -t540, -t539, 0, t512, 0, 0, 0, 0, 0, 0, 0, t540, t539, t510, 0, 0, 0, 0, 0, 0, t517, t518, t524, t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t542, 0, 0, 0, 0, 0, 0, t539, -t540, 0, t511, 0, 0, 0, 0, 0, 0, 0, -t539, t540, t509, 0, 0, 0, 0, 0, 0, t515, t516, t523, t505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, t529, t530, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t563, -qJDD(2), 0, t526, 0, 0, 0, 0, 0, 0, 0, t563, qJDD(2), t521, 0, 0, 0, 0, 0, 0, t536, t537, -t541, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t563, 0, t525, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t563, -t522, 0, 0, 0, 0, 0, 0, -t527, -t528, t538, -t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, t529, t530, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, 0, 0, 0, 0, 0, 0, t529, t530, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t563, -qJDD(2), -t521, 0, 0, 0, 0, 0, 0, -t536, -t537, t541, -t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t563, t522, 0, 0, 0, 0, 0, 0, t527, t528, -t538, t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t546, t545, -t570, t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t544, t547, -t569, t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t536, t537, -t541, t519;];
f_new_reg = t1;
