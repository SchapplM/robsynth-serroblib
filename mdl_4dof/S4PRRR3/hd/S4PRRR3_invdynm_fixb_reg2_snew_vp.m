% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRRR3
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
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRRR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:43
% EndTime: 2019-12-31 16:31:47
% DurationCPUTime: 3.65s
% Computational Cost: add. (9540->286), mult. (15198->391), div. (0->0), fcn. (10702->8), ass. (0->202)
t529 = qJD(2) + qJD(3);
t527 = t529 ^ 2;
t539 = cos(qJ(3));
t528 = qJDD(2) + qJDD(3);
t536 = sin(qJ(3));
t572 = t536 * t528;
t497 = t539 * t527 + t572;
t568 = t539 * t528;
t500 = t536 * t527 - t568;
t537 = sin(qJ(2));
t540 = cos(qJ(2));
t445 = t540 * t497 - t537 * t500;
t532 = g(3) - qJDD(1);
t478 = pkin(5) * t497 - t539 * t532;
t598 = pkin(5) * t500 - t536 * t532;
t394 = pkin(4) * t445 + t540 * t478 - t537 * t598;
t533 = sin(pkin(7));
t534 = cos(pkin(7));
t450 = t537 * t497 + t540 * t500;
t597 = t534 * t445 - t533 * t450;
t604 = pkin(4) * t450 + t537 * t478 + t540 * t598;
t614 = qJ(1) * t597 + t534 * t394 - t533 * t604;
t605 = t533 * t445 + t534 * t450;
t613 = qJ(1) * t605 + t533 * t394 + t534 * t604;
t511 = t533 * g(1) - t534 * g(2);
t512 = t534 * g(1) + t533 * g(2);
t472 = t537 * t511 - t540 * t512;
t542 = qJD(2) ^ 2;
t458 = -t542 * pkin(2) + t472;
t554 = t540 * t511 + t537 * t512;
t544 = qJDD(2) * pkin(2) + t554;
t414 = t536 * t458 - t539 * t544;
t415 = t539 * t458 + t536 * t544;
t558 = t536 * t414 + t539 * t415;
t373 = t539 * t414 - t536 * t415;
t566 = t540 * t373;
t358 = -t537 * t558 + t566;
t570 = t537 * t373;
t601 = t540 * t558 + t570;
t611 = t533 * t358 + t534 * t601;
t610 = t534 * t358 - t533 * t601;
t557 = t540 * t472 - t537 * t554;
t421 = -t537 * t472 - t540 * t554;
t577 = t534 * t421;
t600 = -t533 * t557 + t577;
t581 = t533 * t421;
t599 = t534 * t557 + t581;
t509 = t537 * qJDD(2) + t540 * t542;
t482 = pkin(4) * t509 - t540 * t532;
t510 = t540 * qJDD(2) - t537 * t542;
t547 = -pkin(4) * t510 - t537 * t532;
t586 = t534 * t509 + t533 * t510;
t596 = qJ(1) * t586 + t534 * t482 - t533 * t547;
t459 = -t533 * t509 + t534 * t510;
t595 = -qJ(1) * t459 + t533 * t482 + t534 * t547;
t411 = -t527 * pkin(3) + t528 * pkin(6) + t415;
t535 = sin(qJ(4));
t538 = cos(qJ(4));
t401 = t535 * t411 + t538 * t532;
t402 = t538 * t411 - t535 * t532;
t368 = t535 * t401 + t538 * t402;
t584 = pkin(1) * t532;
t530 = t535 ^ 2;
t582 = t530 * t527;
t578 = t533 * t532;
t576 = t534 * t532;
t410 = -t528 * pkin(3) - t527 * pkin(6) + t414;
t405 = t535 * t410;
t517 = t535 * t527 * t538;
t507 = qJDD(4) + t517;
t575 = t535 * t507;
t508 = qJDD(4) - t517;
t574 = t535 * t508;
t573 = t535 * t528;
t406 = t538 * t410;
t569 = t538 * t508;
t520 = t538 * t528;
t565 = -pkin(3) * t410 + pkin(6) * t368;
t531 = t538 ^ 2;
t564 = t530 + t531;
t563 = qJD(4) * t529;
t541 = qJD(4) ^ 2;
t514 = -t541 - t582;
t468 = -t535 * t514 - t569;
t519 = t538 * t563;
t490 = 0.2e1 * t519 + t573;
t562 = -pkin(3) * t490 + pkin(6) * t468 + t405;
t521 = t531 * t527;
t516 = -t521 - t541;
t466 = t538 * t516 - t575;
t559 = t535 * t563;
t493 = t520 - 0.2e1 * t559;
t561 = pkin(3) * t493 + pkin(6) * t466 - t406;
t361 = t536 * t368 - t539 * t410;
t560 = pkin(2) * t361 + t565;
t555 = -t533 * t511 - t534 * t512;
t553 = t536 * t517;
t552 = t539 * t517;
t495 = t564 * t528;
t501 = t521 + t582;
t551 = pkin(3) * t501 + pkin(6) * t495 + t368;
t428 = t536 * t468 - t539 * t490;
t550 = pkin(2) * t428 + t562;
t427 = t536 * t466 + t539 * t493;
t549 = pkin(2) * t427 + t561;
t548 = -pkin(2) * t500 - t414;
t444 = t536 * t495 + t539 * t501;
t546 = pkin(2) * t444 + t551;
t367 = t538 * t401 - t535 * t402;
t545 = t534 * t511 - t533 * t512;
t543 = -pkin(2) * t497 - t415;
t515 = t521 - t541;
t513 = t541 - t582;
t502 = -t521 + t582;
t496 = t538 * t507;
t492 = t520 - t559;
t491 = t519 + t573;
t486 = t564 * t563;
t474 = t536 * qJDD(4) + t539 * t486;
t473 = -t539 * qJDD(4) + t536 * t486;
t470 = t538 * t491 - t530 * t563;
t469 = -t535 * t492 - t531 * t563;
t467 = -t535 * t513 + t496;
t465 = t538 * t515 - t574;
t463 = t538 * t514 - t574;
t462 = t538 * t513 + t575;
t461 = t535 * t516 + t496;
t460 = t535 * t515 + t569;
t454 = (t491 + t519) * t535;
t453 = (t492 - t559) * t538;
t448 = t539 * t495 - t536 * t501;
t442 = -pkin(1) * t509 - t472;
t441 = pkin(1) * t510 + t554;
t440 = -t535 * t490 + t538 * t493;
t439 = t538 * t490 + t535 * t493;
t438 = t539 * t467 + t535 * t572;
t437 = t539 * t465 + t536 * t520;
t436 = t536 * t467 - t535 * t568;
t435 = t536 * t465 - t538 * t568;
t434 = t539 * t470 - t553;
t433 = t539 * t469 + t553;
t432 = t536 * t470 + t552;
t431 = t536 * t469 - t552;
t430 = t539 * t468 + t536 * t490;
t429 = t539 * t466 - t536 * t493;
t424 = -t537 * t473 + t540 * t474;
t423 = t540 * t473 + t537 * t474;
t420 = t539 * t440 + t536 * t502;
t417 = t536 * t440 - t539 * t502;
t416 = pkin(1) * t421;
t413 = pkin(4) * t557 + t584;
t408 = -t537 * t444 + t540 * t448;
t407 = t540 * t444 + t537 * t448;
t398 = -t537 * t436 + t540 * t438;
t397 = -t537 * t435 + t540 * t437;
t396 = t540 * t436 + t537 * t438;
t395 = t540 * t435 + t537 * t437;
t390 = -t537 * t432 + t540 * t434;
t389 = -t537 * t431 + t540 * t433;
t388 = t540 * t432 + t537 * t434;
t387 = t540 * t431 + t537 * t433;
t386 = -t537 * t428 + t540 * t430;
t385 = -t537 * t427 + t540 * t429;
t384 = t540 * t428 + t537 * t430;
t383 = t540 * t427 + t537 * t429;
t382 = -pkin(6) * t463 + t406;
t381 = -pkin(6) * t461 + t405;
t380 = -pkin(3) * t463 + t402;
t379 = -pkin(3) * t461 + t401;
t378 = -pkin(1) * t450 + t548;
t377 = -pkin(1) * t445 + t543;
t376 = -t537 * t417 + t540 * t420;
t375 = t540 * t417 + t537 * t420;
t370 = pkin(2) * t373;
t369 = pkin(2) * t532 + pkin(5) * t558;
t364 = -pkin(5) * t444 + t539 * t367;
t363 = pkin(5) * t448 + t536 * t367;
t362 = t539 * t368 + t536 * t410;
t355 = -pkin(5) * t428 - t536 * t380 + t539 * t382;
t354 = -pkin(5) * t427 - t536 * t379 + t539 * t381;
t353 = pkin(1) * t384 + t550;
t352 = pkin(1) * t383 + t549;
t351 = -pkin(2) * t463 + pkin(5) * t430 + t539 * t380 + t536 * t382;
t350 = -pkin(2) * t461 + pkin(5) * t429 + t539 * t379 + t536 * t381;
t349 = pkin(1) * t407 + t546;
t348 = -pkin(1) * t358 - t370;
t347 = -pkin(4) * t407 - t537 * t363 + t540 * t364;
t346 = pkin(4) * t408 + t540 * t363 + t537 * t364;
t345 = -t537 * t361 + t540 * t362;
t344 = t540 * t361 + t537 * t362;
t343 = pkin(4) * t358 + pkin(5) * t566 - t537 * t369;
t342 = pkin(4) * t601 + pkin(5) * t570 + t540 * t369 + t584;
t341 = -pkin(5) * t361 - (pkin(3) * t536 - pkin(6) * t539) * t367;
t340 = -pkin(4) * t384 - t537 * t351 + t540 * t355;
t339 = -pkin(4) * t383 - t537 * t350 + t540 * t354;
t338 = -pkin(1) * t463 + pkin(4) * t386 + t540 * t351 + t537 * t355;
t337 = -pkin(1) * t461 + pkin(4) * t385 + t540 * t350 + t537 * t354;
t336 = pkin(5) * t362 - (-pkin(3) * t539 - pkin(6) * t536 - pkin(2)) * t367;
t335 = pkin(1) * t344 + t560;
t334 = -pkin(4) * t344 - t537 * t336 + t540 * t341;
t333 = pkin(1) * t367 + pkin(4) * t345 + t540 * t336 + t537 * t341;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t578, -t576, -t545, -qJ(1) * t545, 0, 0, t459, 0, -t586, 0, t595, t596, t600, pkin(4) * t577 + qJ(1) * t600 - t533 * t413, 0, 0, -t605, 0, -t597, 0, t613, t614, t610, qJ(1) * t610 - t533 * t342 + t534 * t343, -t533 * t388 + t534 * t390, -t533 * t375 + t534 * t376, -t533 * t396 + t534 * t398, -t533 * t387 + t534 * t389, -t533 * t395 + t534 * t397, -t533 * t423 + t534 * t424, t534 * t339 - t533 * t337 - qJ(1) * (t534 * t383 + t533 * t385), t534 * t340 - t533 * t338 - qJ(1) * (t534 * t384 + t533 * t386), t534 * t347 - t533 * t346 - qJ(1) * (t534 * t407 + t533 * t408), t534 * t334 - t533 * t333 - qJ(1) * (t534 * t344 + t533 * t345); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t576, -t578, t555, qJ(1) * t555, 0, 0, t586, 0, t459, 0, -t596, t595, t599, pkin(4) * t581 + qJ(1) * t599 + t534 * t413, 0, 0, t597, 0, -t605, 0, -t614, t613, t611, qJ(1) * t611 + t534 * t342 + t533 * t343, t534 * t388 + t533 * t390, t534 * t375 + t533 * t376, t534 * t396 + t533 * t398, t534 * t387 + t533 * t389, t534 * t395 + t533 * t397, t534 * t423 + t533 * t424, t533 * t339 + t534 * t337 + qJ(1) * (-t533 * t383 + t534 * t385), t533 * t340 + t534 * t338 + qJ(1) * (-t533 * t384 + t534 * t386), t533 * t347 + t534 * t346 + qJ(1) * (-t533 * t407 + t534 * t408), t533 * t334 + t534 * t333 + qJ(1) * (-t533 * t344 + t534 * t345); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t511, t512, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t441, t442, 0, -t416, 0, 0, 0, 0, 0, t528, t378, t377, 0, t348, t454, t439, t462, t453, t460, 0, t352, t353, t349, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t532, -t511, 0, 0, 0, t510, 0, -t509, 0, t547, t482, t421, pkin(4) * t421, 0, 0, -t450, 0, -t445, 0, t604, t394, t358, t343, t390, t376, t398, t389, t397, t424, t339, t340, t347, t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t532, 0, -t512, 0, 0, 0, t509, 0, t510, 0, -t482, t547, t557, t413, 0, 0, t445, 0, -t450, 0, -t394, t604, t601, t342, t388, t375, t396, t387, t395, t423, t337, t338, t346, t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t511, t512, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t441, t442, 0, -t416, 0, 0, 0, 0, 0, t528, t378, t377, 0, t348, t454, t439, t462, t453, t460, 0, t352, t353, t349, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t542, 0, 0, -t532, -t554, 0, 0, 0, -t500, 0, -t497, 0, t598, t478, t373, pkin(5) * t373, t434, t420, t438, t433, t437, t474, t354, t355, t364, t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t542, 0, qJDD(2), 0, t532, 0, t472, 0, 0, 0, t497, 0, -t500, 0, -t478, t598, t558, t369, t432, t417, t436, t431, t435, t473, t350, t351, t363, t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t554, -t472, 0, 0, 0, 0, 0, 0, 0, t528, t548, t543, 0, -t370, t454, t439, t462, t453, t460, 0, t549, t550, t546, t560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t528, 0, -t527, 0, 0, -t532, t414, 0, t470, t440, t467, t469, t465, t486, t381, t382, t367, pkin(6) * t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t527, 0, t528, 0, t532, 0, t415, 0, t517, -t502, -t573, -t517, -t520, -qJDD(4), t379, t380, 0, pkin(3) * t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t528, -t414, -t415, 0, 0, t454, t439, t462, t453, t460, 0, t561, t562, t551, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t491, t493, t507, -t519, t515, t519, 0, t410, t401, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t559, t490, t513, t492, t508, -t559, -t410, 0, t402, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t517, t502, t573, t517, t520, qJDD(4), -t401, -t402, 0, 0;];
m_new_reg = t1;