% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:49
% EndTime: 2019-03-09 02:54:59
% DurationCPUTime: 7.24s
% Computational Cost: add. (4432->472), mult. (8939->617), div. (0->0), fcn. (6324->14), ass. (0->204)
t469 = sin(qJ(3));
t472 = cos(qJ(3));
t561 = sin(pkin(9));
t562 = cos(pkin(9));
t483 = t469 * t562 + t472 * t561;
t571 = t483 * qJD(1);
t398 = qJD(6) + t571;
t411 = -t561 * t469 + t562 * t472;
t404 = t411 * qJD(1);
t465 = sin(pkin(10));
t466 = cos(pkin(10));
t379 = -t466 * qJD(3) + t404 * t465;
t471 = cos(qJ(6));
t381 = qJD(3) * t465 + t404 * t466;
t468 = sin(qJ(6));
t553 = t381 * t468;
t578 = -t471 * t379 - t553;
t579 = t578 * t398;
t474 = -pkin(1) - pkin(7);
t422 = qJDD(1) * t474 + qJDD(2);
t415 = t472 * t422;
t423 = qJD(1) * t474 + qJD(2);
t523 = qJDD(1) * t472;
t525 = qJD(1) * qJD(4);
t526 = qJD(1) * qJD(3);
t532 = qJD(3) * t469;
t344 = -t472 * t525 - t423 * t532 + qJDD(3) * pkin(3) + t415 + (t469 * t526 - t523) * qJ(4);
t531 = qJD(3) * t472;
t359 = (-qJ(4) * qJD(1) + t423) * t531 + (-qJ(4) * qJDD(1) + t422 - t525) * t469;
t316 = t344 * t562 - t561 * t359;
t313 = -qJDD(3) * pkin(4) + qJDD(5) - t316;
t459 = qJ(3) + pkin(9);
t448 = sin(t459);
t450 = cos(t459);
t470 = sin(qJ(1));
t473 = cos(qJ(1));
t572 = g(1) * t470 - g(2) * t473;
t482 = -g(3) * t448 + t450 * t572;
t574 = t313 + t482;
t493 = t379 * t468 - t381 * t471;
t577 = t398 * t493;
t534 = qJD(1) * t469;
t396 = -qJ(4) * t534 + t423 * t469;
t389 = t561 * t396;
t533 = qJD(1) * t472;
t397 = -qJ(4) * t533 + t472 * t423;
t392 = qJD(3) * pkin(3) + t397;
t350 = t392 * t562 - t389;
t336 = -qJD(3) * pkin(4) + qJD(5) - t350;
t513 = qJD(3) * t561;
t514 = qJD(3) * t562;
t403 = -t469 * t514 - t472 * t513;
t479 = t313 * t411 + t336 * t403 + t572;
t476 = qJD(1) ^ 2;
t484 = -qJ(2) * t476 - t572;
t505 = g(1) * t473 + g(2) * t470;
t575 = t448 * t505;
t509 = qJDD(1) * t561;
t510 = qJDD(1) * t562;
t371 = -qJD(3) * t404 - t469 * t510 - t472 * t509;
t369 = -qJDD(6) + t371;
t414 = t465 * t471 + t466 * t468;
t407 = t414 * qJD(6);
t413 = t465 * t468 - t471 * t466;
t573 = t413 * t369 - t407 * t398;
t460 = qJDD(1) * qJ(2);
t461 = qJD(1) * qJD(2);
t521 = 0.2e1 * t461;
t570 = 0.2e1 * t460 + t521 - t505;
t406 = t413 * qJD(6);
t541 = -t413 * t571 - t406;
t556 = t369 * t414;
t569 = -t398 * t541 + t556;
t567 = pkin(8) * t466;
t565 = g(3) * t450;
t564 = g(3) * t469;
t454 = t469 * pkin(3);
t434 = pkin(3) * t561 + qJ(5);
t563 = pkin(8) + t434;
t560 = pkin(1) * qJDD(1);
t558 = t578 * t404;
t557 = t493 * t404;
t555 = t371 * t465;
t554 = t371 * t466;
t402 = t469 * t513 - t472 * t514;
t552 = t398 * t402;
t551 = t571 * t465;
t550 = t403 * t465;
t549 = t411 * t465;
t548 = t411 * t466;
t458 = pkin(10) + qJ(6);
t447 = sin(t458);
t547 = t447 * t470;
t546 = t447 * t473;
t449 = cos(t458);
t545 = t449 * t470;
t544 = t449 * t473;
t543 = qJ(2) + t454;
t542 = qJ(4) - t474;
t372 = qJD(3) * t571 + t469 * t509 - t472 * t510;
t519 = t472 * t526;
t524 = qJDD(1) * t469;
t492 = qJDD(4) + t460 + t461 + (t519 + t524) * pkin(3);
t307 = -pkin(4) * t371 + qJ(5) * t372 - qJD(5) * t404 + t492;
t317 = t561 * t344 + t562 * t359;
t311 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t317;
t295 = t465 * t307 + t466 * t311;
t517 = t562 * t396;
t351 = t561 * t392 + t517;
t337 = qJD(3) * qJ(5) + t351;
t417 = pkin(3) * t534 + qJD(1) * qJ(2) + qJD(4);
t352 = pkin(4) * t571 - qJ(5) * t404 + t417;
t315 = t466 * t337 + t465 * t352;
t528 = pkin(3) * t531 + qJD(2);
t330 = -pkin(4) * t402 - qJ(5) * t403 - qJD(5) * t411 + t528;
t394 = -qJD(4) * t472 + t532 * t542;
t419 = t542 * t472;
t395 = -qJD(3) * t419 - qJD(4) * t469;
t358 = t394 * t561 + t395 * t562;
t309 = t465 * t330 + t466 * t358;
t363 = t397 * t562 - t389;
t364 = pkin(3) * t533 + pkin(4) * t404 + qJ(5) * t571;
t319 = t466 * t363 + t465 * t364;
t370 = pkin(4) * t483 - qJ(5) * t411 + t543;
t418 = t542 * t469;
t377 = -t418 * t562 - t419 * t561;
t325 = t465 * t370 + t466 * t377;
t355 = t414 * t571;
t540 = t407 + t355;
t539 = t473 * pkin(1) + t470 * qJ(2);
t464 = t472 ^ 2;
t537 = t469 ^ 2 - t464;
t475 = qJD(3) ^ 2;
t536 = -t475 - t476;
t535 = qJD(1) * t417;
t530 = qJD(6) * t468;
t529 = qJD(6) * t471;
t527 = -qJD(5) + t336;
t522 = qJDD(3) * t469;
t360 = -t466 * qJDD(3) - t372 * t465;
t361 = qJDD(3) * t465 - t372 * t466;
t520 = -t468 * t360 + t471 * t361 - t379 * t529;
t518 = -pkin(1) * t470 + t473 * qJ(2);
t294 = t466 * t307 - t311 * t465;
t290 = -pkin(5) * t371 - pkin(8) * t361 + t294;
t293 = -pkin(8) * t360 + t295;
t512 = t471 * t290 - t293 * t468;
t308 = t466 * t330 - t358 * t465;
t314 = -t337 * t465 + t466 * t352;
t511 = t471 * t360 + t468 * t361;
t318 = -t363 * t465 + t466 * t364;
t324 = t466 * t370 - t377 * t465;
t508 = qJDD(2) - t560;
t439 = -pkin(3) * t562 - pkin(4);
t503 = -t355 * t398 + t573;
t502 = pkin(4) * t448 - qJ(5) * t450;
t357 = -t562 * t394 + t395 * t561;
t362 = t397 * t561 + t517;
t376 = -t418 * t561 + t562 * t419;
t501 = t290 * t468 + t293 * t471;
t500 = -t294 * t466 - t295 * t465;
t499 = t294 * t483 - t314 * t402;
t498 = -t295 * t483 + t315 * t402;
t300 = pkin(5) * t571 - pkin(8) * t381 + t314;
t303 = -pkin(8) * t379 + t315;
t291 = t300 * t471 - t303 * t468;
t292 = t300 * t468 + t303 * t471;
t312 = pkin(5) * t483 - pkin(8) * t548 + t324;
t320 = -pkin(8) * t549 + t325;
t497 = t312 * t471 - t320 * t468;
t496 = t312 * t468 + t320 * t471;
t494 = -t314 * t465 + t315 * t466;
t467 = -qJ(4) - pkin(7);
t491 = t473 * t454 + t470 * t467 + t518;
t490 = t470 * t454 - t467 * t473 + t539;
t409 = t563 * t466;
t489 = pkin(5) * t404 + qJD(5) * t465 + qJD(6) * t409 + t567 * t571 + t318;
t408 = t563 * t465;
t488 = pkin(8) * t551 - qJD(5) * t466 + qJD(6) * t408 + t319;
t487 = t398 * t414;
t486 = 0.2e1 * qJ(2) * t526 + qJDD(3) * t474;
t296 = -t381 * t530 + t520;
t297 = -qJD(6) * t493 + t511;
t478 = t316 * t411 + t317 * t483 + t350 * t403 - t351 * t402 - t572;
t477 = -t474 * t475 + t570;
t451 = qJDD(3) * t472;
t420 = -t466 * pkin(5) + t439;
t399 = t571 ^ 2;
t388 = t448 * t544 - t547;
t387 = t448 * t546 + t545;
t386 = t448 * t545 + t546;
t385 = -t448 * t547 + t544;
t366 = t413 * t411;
t365 = t414 * t411;
t345 = pkin(5) * t549 + t376;
t327 = -pkin(5) * t551 + t362;
t326 = pkin(5) * t550 + t357;
t323 = t379 * pkin(5) + t336;
t322 = t403 * t414 + t529 * t548 - t530 * t549;
t321 = -t403 * t413 - t407 * t411;
t301 = -pkin(8) * t550 + t309;
t299 = -pkin(5) * t402 - t403 * t567 + t308;
t298 = t360 * pkin(5) + t313;
t1 = [(-t296 * t365 + t297 * t366 + t321 * t578 + t322 * t493) * MDP(21) + ((t299 * t471 - t301 * t468) * t398 - t497 * t369 + t512 * t483 - t291 * t402 - t326 * t578 + t345 * t297 + t298 * t365 + t323 * t322 - g(1) * t388 - g(2) * t386 + (-t292 * t483 - t398 * t496) * qJD(6)) * MDP(25) + (-t297 * t483 - t322 * t398 + t365 * t369 - t402 * t578) * MDP(23) + t572 * MDP(2) + (qJDD(2) - t572 - 0.2e1 * t560) * MDP(4) + (-t296 * t366 - t321 * t493) * MDP(20) + (t357 * t404 - t358 * t571 + t371 * t377 - t372 * t376 - t478) * MDP(14) + (t308 * t571 - t324 * t371 + t357 * t379 + t376 * t360 + t465 * t479 - t466 * t575 + t499) * MDP(16) + (-t309 * t571 + t325 * t371 + t357 * t381 + t376 * t361 + t465 * t575 + t466 * t479 + t498) * MDP(17) + (-t369 * t483 - t552) * MDP(24) + (-(t299 * t468 + t301 * t471) * t398 + t496 * t369 - t501 * t483 + t292 * t402 - t326 * t493 + t345 * t296 - t298 * t366 + t323 * t321 + g(1) * t387 - g(2) * t385 + (-t291 * t483 - t398 * t497) * qJD(6)) * MDP(26) + (t296 * t483 + t321 * t398 + t366 * t369 + t402 * t493) * MDP(22) + t570 * MDP(5) + qJDD(1) * MDP(1) + (t469 * t477 + t472 * t486) * MDP(12) + (-t469 * t486 + t472 * t477) * MDP(13) + (t295 * t325 + t315 * t309 + t294 * t324 + t314 * t308 + t313 * t376 + t336 * t357 - g(1) * (t473 * t502 + t491) - g(2) * (t470 * t502 + t490)) * MDP(19) + (-t308 * t381 - t309 * t379 - t324 * t361 - t325 * t360 + t505 * t450 + t500 * t411 + (-t314 * t466 - t315 * t465) * t403) * MDP(18) + t505 * MDP(3) + (qJDD(1) * t464 - 0.2e1 * t469 * t519) * MDP(7) + (-t472 * t475 - t522) * MDP(10) + 0.2e1 * (-t469 * t523 + t526 * t537) * MDP(8) + (-t508 * pkin(1) - g(1) * t518 - g(2) * t539 + (t521 + t460) * qJ(2)) * MDP(6) + (-g(1) * t491 - g(2) * t490 - t316 * t376 + t317 * t377 - t350 * t357 + t351 * t358 + t417 * t528 + t492 * t543) * MDP(15) + (-t469 * t475 + t451) * MDP(9); qJDD(1) * MDP(4) - t476 * MDP(5) + (t508 + t484) * MDP(6) + (t469 * t536 + t451) * MDP(12) + (t472 * t536 - t522) * MDP(13) + (t371 * t483 + t372 * t411 + t402 * t571 - t403 * t404) * MDP(14) + (t478 - t535) * MDP(15) + (t483 * t555 - t360 * t411 - t379 * t403 + (-qJD(1) * t466 + t402 * t465) * t571) * MDP(16) + (t483 * t554 - t361 * t411 - t381 * t403 + (qJD(1) * t465 + t402 * t466) * t571) * MDP(17) + ((qJD(1) * t381 - t360 * t483 + t379 * t402) * t466 + (qJD(1) * t379 + t361 * t483 - t381 * t402) * t465) * MDP(18) + ((-qJD(1) * t314 - t498) * t466 + (-qJD(1) * t315 - t499) * t465 - t479) * MDP(19) + (-t411 * t297 + t403 * t578 + t402 * t487 + t413 * t398 * qJD(1) - (-t398 * t406 - t556) * t483) * MDP(25) + (qJD(1) * t487 - t411 * t296 + t403 * t493 - t413 * t552 - t483 * t573) * MDP(26); MDP(9) * t523 - MDP(10) * t524 + qJDD(3) * MDP(11) + (t472 * t484 + t415 + t564) * MDP(12) + (g(3) * t472 + (-t422 - t484) * t469) * MDP(13) + ((t351 - t362) * t404 + (-t350 + t363) * t571 + (t371 * t561 + t372 * t562) * pkin(3)) * MDP(14) + (t350 * t362 - t351 * t363 + (t561 * t317 + t562 * t316 + t564 + (-t572 - t535) * t472) * pkin(3)) * MDP(15) + (t434 * t555 - t314 * t404 + t360 * t439 - t362 * t379 + (t465 * t527 - t318) * t571 - t574 * t466) * MDP(16) + (t434 * t554 + t315 * t404 + t361 * t439 - t362 * t381 + (t466 * t527 + t319) * t571 + t574 * t465) * MDP(17) + (-t565 + t318 * t381 + t319 * t379 - t572 * t448 + (-qJD(5) * t379 - t314 * t571 - t360 * t434 + t295) * t466 + (qJD(5) * t381 - t315 * t571 + t361 * t434 - t294) * t465) * MDP(18) + (t313 * t439 - t315 * t319 - t314 * t318 - t336 * t362 - g(3) * (-t454 - t502) + (-t294 * t465 + t295 * t466) * t434 + t494 * qJD(5) - t572 * (pkin(3) * t472 + pkin(4) * t450 + qJ(5) * t448)) * MDP(19) + (t296 * t414 - t493 * t541) * MDP(20) + (-t296 * t413 - t297 * t414 + t493 * t540 + t541 * t578) * MDP(21) + (t557 - t569) * MDP(22) + (t503 - t558) * MDP(23) - t398 * t404 * MDP(24) + (-(-t408 * t471 - t409 * t468) * t369 + t420 * t297 + t298 * t413 - t291 * t404 + t327 * t578 + (t468 * t488 - t471 * t489) * t398 + t540 * t323 - t482 * t449) * MDP(25) + ((-t408 * t468 + t409 * t471) * t369 + t420 * t296 + t298 * t414 + t292 * t404 + t327 * t493 + (t468 * t489 + t471 * t488) * t398 + t541 * t323 + t482 * t447) * MDP(26) + (MDP(7) * t469 * t472 - MDP(8) * t537) * t476; (-t404 ^ 2 - t399) * MDP(14) + (t350 * t404 + t492 - t505) * MDP(15) + (-t379 * t404 - t554) * MDP(16) + (-t381 * t404 - t399 * t466 + t555) * MDP(17) + (-t360 * t465 - t361 * t466) * MDP(18) + (-t336 * t404 - t500 - t505) * MDP(19) + (t503 + t558) * MDP(25) + (t557 + t569) * MDP(26) + (t351 * MDP(15) + (-t379 * t466 + t381 * t465) * MDP(18) + t494 * MDP(19) - MDP(16) * t551) * t571; (t381 * t571 + t360) * MDP(16) + (-t379 * t571 + t361) * MDP(17) + (-t379 ^ 2 - t381 ^ 2) * MDP(18) + (t314 * t381 + t315 * t379 + t574) * MDP(19) + (t297 - t577) * MDP(25) + (t296 + t579) * MDP(26); t493 * t578 * MDP(20) + (t493 ^ 2 - t578 ^ 2) * MDP(21) + (t520 - t579) * MDP(22) + (-t511 - t577) * MDP(23) - t369 * MDP(24) + (-g(1) * t385 - g(2) * t387 + t292 * t398 + t323 * t493 + t447 * t565 + t512) * MDP(25) + (g(1) * t386 - g(2) * t388 + t291 * t398 - t323 * t578 + t449 * t565 - t501) * MDP(26) + (-MDP(22) * t553 + MDP(23) * t493 - MDP(25) * t292 - MDP(26) * t291) * qJD(6);];
tau  = t1;
