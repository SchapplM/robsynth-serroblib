% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:04:44
% EndTime: 2019-05-04 19:04:45
% DurationCPUTime: 0.83s
% Computational Cost: add. (2324->76), mult. (3626->99), div. (0->0), fcn. (2872->8), ass. (0->59)
t567 = qJD(2) + qJD(3);
t563 = qJD(4) + t567;
t561 = t563 ^ 2;
t566 = qJDD(2) + qJDD(3);
t562 = qJDD(4) + t566;
t571 = sin(qJ(4));
t574 = cos(qJ(4));
t547 = t571 * t561 - t574 * t562;
t572 = sin(qJ(3));
t575 = cos(qJ(3));
t579 = -t574 * t561 - t571 * t562;
t532 = t575 * t547 - t572 * t579;
t573 = sin(qJ(2));
t576 = cos(qJ(2));
t587 = t572 * t547 + t575 * t579;
t522 = t576 * t532 - t573 * t587;
t569 = sin(pkin(7));
t570 = cos(pkin(7));
t594 = t573 * t532 + t576 * t587;
t598 = t569 * t522 + t570 * t594;
t597 = t570 * t522 - t569 * t594;
t565 = t567 ^ 2;
t551 = t572 * t565 - t575 * t566;
t578 = -t575 * t565 - t572 * t566;
t536 = t576 * t551 - t573 * t578;
t586 = t573 * t551 + t576 * t578;
t593 = t569 * t536 + t570 * t586;
t592 = t570 * t536 - t569 * t586;
t559 = t569 * g(1) - t570 * g(2);
t560 = -t570 * g(1) - t569 * g(2);
t543 = t576 * t559 - t573 * t560;
t541 = qJDD(2) * pkin(2) + t543;
t544 = t573 * t559 + t576 * t560;
t577 = qJD(2) ^ 2;
t542 = -t577 * pkin(2) + t544;
t527 = t572 * t541 + t575 * t542;
t526 = t575 * t541 - t572 * t542;
t557 = t576 * qJDD(2) - t573 * t577;
t558 = -t573 * qJDD(2) - t576 * t577;
t581 = -t569 * t557 + t570 * t558;
t580 = t570 * t557 + t569 * t558;
t568 = -g(3) + qJDD(1);
t529 = -t573 * t543 + t576 * t544;
t528 = t576 * t543 + t573 * t544;
t525 = -t565 * pkin(3) + t527;
t524 = t566 * pkin(3) + t526;
t519 = -t572 * t526 + t575 * t527;
t518 = t575 * t526 + t572 * t527;
t517 = t571 * t524 + t574 * t525;
t516 = t574 * t524 - t571 * t525;
t515 = -t573 * t518 + t576 * t519;
t514 = t576 * t518 + t573 * t519;
t513 = -t571 * t516 + t574 * t517;
t512 = t574 * t516 + t571 * t517;
t511 = -t572 * t512 + t575 * t513;
t510 = t575 * t512 + t572 * t513;
t509 = -t573 * t510 + t576 * t511;
t508 = t576 * t510 + t573 * t511;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t569 * t559 + t570 * t560, 0, 0, 0, 0, 0, 0, t581, -t580, 0, -t569 * t528 + t570 * t529, 0, 0, 0, 0, 0, 0, t593, t592, 0, -t569 * t514 + t570 * t515, 0, 0, 0, 0, 0, 0, t598, t597, 0, -t569 * t508 + t570 * t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t570 * t559 + t569 * t560, 0, 0, 0, 0, 0, 0, t580, t581, 0, t570 * t528 + t569 * t529, 0, 0, 0, 0, 0, 0, -t592, t593, 0, t570 * t514 + t569 * t515, 0, 0, 0, 0, 0, 0, -t597, t598, 0, t570 * t508 + t569 * t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t560, 0, 0, 0, 0, 0, 0, t558, -t557, 0, t529, 0, 0, 0, 0, 0, 0, t586, t536, 0, t515, 0, 0, 0, 0, 0, 0, t594, t522, 0, t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t559, 0, 0, 0, 0, 0, 0, t557, t558, 0, t528, 0, 0, 0, 0, 0, 0, -t536, t586, 0, t514, 0, 0, 0, 0, 0, 0, -t522, t594, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, -qJDD(2), 0, t544, 0, 0, 0, 0, 0, 0, t578, t551, 0, t519, 0, 0, 0, 0, 0, 0, t587, t532, 0, t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t577, 0, t543, 0, 0, 0, 0, 0, 0, -t551, t578, 0, t518, 0, 0, 0, 0, 0, 0, -t532, t587, 0, t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t565, -t566, 0, t527, 0, 0, 0, 0, 0, 0, t579, t547, 0, t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t566, -t565, 0, t526, 0, 0, 0, 0, 0, 0, -t547, t579, 0, t512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t561, -t562, 0, t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, -t561, 0, t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568;];
f_new_reg  = t1;
