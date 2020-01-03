% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:47
% EndTime: 2019-12-31 16:18:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (1154->86), mult. (2094->128), div. (0->0), fcn. (1622->8), ass. (0->72)
t565 = sin(pkin(6));
t567 = cos(pkin(6));
t550 = t565 * g(1) - t567 * g(2);
t545 = -qJDD(2) + t550;
t580 = t565 * t545;
t551 = -t567 * g(1) - t565 * g(2);
t562 = -g(3) + qJDD(1);
t564 = sin(pkin(7));
t566 = cos(pkin(7));
t539 = -t564 * t551 + t566 * t562;
t540 = t566 * t551 + t564 * t562;
t570 = sin(qJ(3));
t572 = cos(qJ(3));
t521 = t570 * t539 + t572 * t540;
t569 = sin(qJ(4));
t560 = t569 ^ 2;
t571 = cos(qJ(4));
t561 = t571 ^ 2;
t579 = t560 + t561;
t578 = qJD(3) * qJD(4);
t577 = t569 * qJDD(3);
t576 = t571 * qJDD(3);
t520 = t572 * t539 - t570 * t540;
t574 = qJD(3) ^ 2;
t547 = t572 * qJDD(3) - t570 * t574;
t548 = -t570 * qJDD(3) - t572 * t574;
t575 = -t564 * t547 + t566 * t548;
t527 = t566 * t547 + t564 * t548;
t573 = qJD(4) ^ 2;
t556 = t569 * t574 * t571;
t555 = -t561 * t574 - t573;
t554 = -t560 * t574 - t573;
t553 = -qJDD(4) + t556;
t552 = qJDD(4) + t556;
t549 = t579 * t574;
t546 = t579 * qJDD(3);
t544 = -0.2e1 * t569 * t578 + t576;
t543 = 0.2e1 * t571 * t578 + t577;
t541 = t567 * t545;
t535 = t571 * t553 - t569 * t554;
t534 = -t569 * t552 + t571 * t555;
t533 = t569 * t553 + t571 * t554;
t532 = t571 * t552 + t569 * t555;
t531 = t572 * t546 - t570 * t549;
t530 = t570 * t546 + t572 * t549;
t525 = t572 * t535 + t570 * t543;
t524 = t572 * t534 - t570 * t544;
t523 = t570 * t535 - t572 * t543;
t522 = t570 * t534 + t572 * t544;
t519 = -t564 * t539 + t566 * t540;
t518 = t566 * t539 + t564 * t540;
t517 = -t574 * pkin(3) + qJDD(3) * pkin(5) + t521;
t516 = -qJDD(3) * pkin(3) - t574 * pkin(5) - t520;
t515 = -t564 * t530 + t566 * t531;
t514 = t566 * t530 + t564 * t531;
t513 = t571 * t517 - t569 * t545;
t512 = -t569 * t517 - t571 * t545;
t511 = -t564 * t523 + t566 * t525;
t510 = -t564 * t522 + t566 * t524;
t509 = t566 * t523 + t564 * t525;
t508 = t566 * t522 + t564 * t524;
t507 = -t570 * t520 + t572 * t521;
t506 = t572 * t520 + t570 * t521;
t505 = -t569 * t512 + t571 * t513;
t504 = t571 * t512 + t569 * t513;
t503 = t572 * t505 + t570 * t516;
t502 = t570 * t505 - t572 * t516;
t501 = -t564 * t506 + t566 * t507;
t500 = t566 * t506 + t564 * t507;
t499 = -t564 * t502 + t566 * t503;
t498 = t566 * t502 + t564 * t503;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t565 * t550 + t567 * t551, 0, 0, 0, 0, 0, 0, 0, 0, 0, t567 * t519 - t580, 0, 0, 0, 0, 0, 0, t567 * t575, -t567 * t527, 0, t567 * t501 - t580, 0, 0, 0, 0, 0, 0, t567 * t510 + t565 * t532, t567 * t511 + t565 * t533, t567 * t515, t567 * t499 + t565 * t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t567 * t550 + t565 * t551, 0, 0, 0, 0, 0, 0, 0, 0, 0, t565 * t519 + t541, 0, 0, 0, 0, 0, 0, t565 * t575, -t565 * t527, 0, t565 * t501 + t541, 0, 0, 0, 0, 0, 0, t565 * t510 - t567 * t532, t565 * t511 - t567 * t533, t565 * t515, t565 * t499 - t567 * t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, 0, 0, 0, 0, 0, 0, 0, 0, 0, t518, 0, 0, 0, 0, 0, 0, t527, t575, 0, t500, 0, 0, 0, 0, 0, 0, t508, t509, t514, t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t551, 0, 0, 0, 0, 0, 0, 0, 0, 0, t519, 0, 0, 0, 0, 0, 0, t575, -t527, 0, t501, 0, 0, 0, 0, 0, 0, t510, t511, t515, t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t550, 0, 0, 0, 0, 0, 0, 0, 0, 0, t545, 0, 0, 0, 0, 0, 0, 0, 0, 0, t545, 0, 0, 0, 0, 0, 0, -t532, -t533, 0, -t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, 0, 0, 0, 0, 0, 0, 0, 0, 0, t518, 0, 0, 0, 0, 0, 0, t527, t575, 0, t500, 0, 0, 0, 0, 0, 0, t508, t509, t514, t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t540, 0, 0, 0, 0, 0, 0, t548, -t547, 0, t507, 0, 0, 0, 0, 0, 0, t524, t525, t531, t503; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t539, 0, 0, 0, 0, 0, 0, t547, t548, 0, t506, 0, 0, 0, 0, 0, 0, t522, t523, t530, t502; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t545, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t545, 0, 0, 0, 0, 0, 0, t532, t533, 0, t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t574, -qJDD(3), 0, t521, 0, 0, 0, 0, 0, 0, t534, t535, t546, t505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t574, 0, t520, 0, 0, 0, 0, 0, 0, t544, -t543, t549, -t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t545, 0, 0, 0, 0, 0, 0, t532, t533, 0, t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t555, t553, t576, t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t552, t554, -t577, t512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t544, t543, -t549, t516;];
f_new_reg = t1;
