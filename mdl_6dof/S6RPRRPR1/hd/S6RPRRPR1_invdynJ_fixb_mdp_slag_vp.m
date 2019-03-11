% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:56
% EndTime: 2019-03-09 04:59:03
% DurationCPUTime: 4.78s
% Computational Cost: add. (6079->424), mult. (13505->556), div. (0->0), fcn. (9889->18), ass. (0->213)
t494 = qJ(3) + qJ(4);
t485 = sin(t494);
t486 = cos(t494);
t491 = qJ(1) + pkin(10);
t480 = sin(t491);
t481 = cos(t491);
t543 = g(1) * t481 + g(2) * t480;
t624 = -g(3) * t486 + t485 * t543;
t496 = sin(pkin(10));
t470 = pkin(1) * t496 + pkin(7);
t604 = pkin(8) + t470;
t490 = qJD(3) + qJD(4);
t498 = sin(qJ(6));
t502 = cos(qJ(6));
t503 = cos(qJ(4));
t504 = cos(qJ(3));
t577 = qJD(1) * t504;
t563 = t503 * t577;
t499 = sin(qJ(4));
t500 = sin(qJ(3));
t578 = qJD(1) * t500;
t564 = t499 * t578;
t428 = -t563 + t564;
t430 = -t499 * t577 - t503 * t578;
t495 = sin(pkin(11));
t601 = cos(pkin(11));
t523 = -t495 * t428 - t430 * t601;
t382 = -t502 * t490 + t498 * t523;
t552 = -t601 * t428 + t430 * t495;
t615 = qJD(6) - t552;
t623 = t382 * t615;
t384 = t490 * t498 + t502 * t523;
t622 = t384 * t615;
t570 = qJD(1) * qJD(3);
t562 = t504 * t570;
t568 = qJDD(1) * t504;
t569 = qJDD(1) * t500;
t378 = qJD(4) * t563 - t490 * t564 + t499 * t568 + (t562 + t569) * t503;
t488 = qJDD(3) + qJDD(4);
t557 = t604 * qJD(1);
t403 = qJD(2) * t500 + t504 * t557;
t482 = t504 * qJDD(2);
t448 = t470 * qJDD(1);
t556 = pkin(8) * qJDD(1) + t448;
t373 = qJDD(3) * pkin(3) - qJD(3) * t403 - t556 * t500 + t482;
t402 = t504 * qJD(2) - t557 * t500;
t376 = qJD(3) * t402 + t500 * qJDD(2) + t556 * t504;
t400 = t503 * t403;
t602 = qJD(3) * pkin(3);
t401 = t402 + t602;
t534 = -t401 * t499 - t400;
t516 = qJD(4) * t534 + t503 * t373 - t499 * t376;
t328 = pkin(4) * t488 - qJ(5) * t378 + qJD(5) * t430 + t516;
t437 = t499 * t504 + t500 * t503;
t396 = t490 * t437;
t540 = t499 * t569 - t503 * t568;
t379 = qJD(1) * t396 + t540;
t575 = qJD(4) * t499;
t613 = t503 * (qJD(4) * t401 + t376) + t499 * t373 - t403 * t575;
t330 = -qJ(5) * t379 - qJD(5) * t428 + t613;
t316 = t328 * t601 - t495 * t330;
t314 = -t488 * pkin(5) - t316;
t483 = pkin(11) + t494;
t468 = cos(t483);
t620 = g(3) * t468 + t314;
t478 = pkin(3) * t578;
t497 = cos(pkin(10));
t472 = -pkin(1) * t497 - pkin(2);
t487 = t504 * pkin(3);
t616 = t472 - t487;
t408 = qJD(3) * t478 + qJDD(1) * t616;
t513 = pkin(4) * t379 + qJDD(5) + t408;
t550 = t615 * t502;
t352 = -t378 * t495 - t601 * t379;
t351 = qJDD(6) - t352;
t587 = t498 * t351;
t619 = -t615 * t550 - t587;
t610 = pkin(4) * t430;
t357 = pkin(5) * t523 - pkin(9) * t552 - t610;
t469 = pkin(4) * t495 + pkin(9);
t618 = t615 * (qJD(6) * t469 + t357);
t434 = t604 * t500;
t435 = t604 * t504;
t581 = -t499 * t434 + t503 * t435;
t424 = t430 * qJ(5);
t398 = t499 * t403;
t554 = t503 * t401 - t398;
t365 = t424 + t554;
t436 = t499 * t500 - t503 * t504;
t559 = qJD(3) * t604;
t425 = t500 * t559;
t426 = t504 * t559;
t574 = qJD(4) * t503;
t522 = -t503 * t425 - t499 * t426 - t434 * t574 - t435 * t575;
t346 = -qJ(5) * t396 - qJD(5) * t436 + t522;
t395 = t490 * t436;
t515 = -qJD(4) * t581 + t425 * t499 - t503 * t426;
t510 = qJ(5) * t395 - qJD(5) * t437 + t515;
t326 = t346 * t601 + t495 * t510;
t361 = pkin(4) * t490 + t365;
t600 = qJ(5) * t428;
t366 = -t534 - t600;
t589 = t495 * t366;
t339 = t361 * t601 - t589;
t337 = -t490 * pkin(5) - t339;
t380 = -qJ(5) * t436 + t581;
t551 = -t503 * t434 - t435 * t499;
t527 = -qJ(5) * t437 + t551;
t355 = t380 * t601 + t495 * t527;
t392 = t436 * t601 + t437 * t495;
t393 = -t495 * t436 + t437 * t601;
t524 = pkin(4) * t436 + t616;
t358 = pkin(5) * t392 - pkin(9) * t393 + t524;
t371 = -t395 * t601 - t495 * t396;
t317 = t495 * t328 + t601 * t330;
t315 = pkin(9) * t488 + t317;
t431 = t616 * qJD(1);
t391 = pkin(4) * t428 + qJD(5) + t431;
t348 = -pkin(5) * t552 - pkin(9) * t523 + t391;
t548 = qJD(6) * t348 + t315;
t612 = t314 * t393 + t337 * t371 - t355 * t351 - (qJD(6) * t358 + t326) * t615 - t392 * t548;
t467 = sin(t483);
t607 = g(3) * t467;
t603 = pkin(3) * qJD(4);
t353 = t378 * t601 - t495 * t379;
t572 = qJD(6) * t502;
t565 = t502 * t353 + t498 * t488 + t490 * t572;
t573 = qJD(6) * t498;
t335 = -t523 * t573 + t565;
t599 = t335 * t498;
t598 = t337 * t552;
t597 = t337 * t393;
t596 = t358 * t351;
t595 = t382 * t523;
t594 = t384 * t523;
t593 = t480 * t498;
t592 = t480 * t502;
t591 = t481 * t498;
t590 = t481 * t502;
t588 = t495 * t499;
t349 = t502 * t351;
t586 = qJDD(2) - g(3);
t370 = -t395 * t495 + t396 * t601;
t585 = t335 * t392 + t384 * t370;
t362 = t601 * t366;
t340 = t495 * t361 + t362;
t584 = t503 * t402 - t398;
t367 = t424 + t584;
t553 = -t402 * t499 - t400;
t528 = t553 + t600;
t558 = t601 * t499;
t583 = -t367 * t495 + t601 * t528 + (t495 * t503 + t558) * t603;
t582 = -t601 * t367 - t495 * t528 + (t503 * t601 - t588) * t603;
t477 = pkin(3) * t503 + pkin(4);
t423 = pkin(3) * t558 + t495 * t477;
t580 = pkin(4) * t486 + t487;
t492 = t500 ^ 2;
t579 = -t504 ^ 2 + t492;
t451 = qJD(1) * t472;
t479 = t500 * t602;
t567 = t393 * t587;
t566 = t393 * t349;
t561 = pkin(4) * t396 + t479;
t322 = -pkin(5) * t352 - pkin(9) * t353 + t513;
t338 = pkin(9) * t490 + t340;
t547 = qJD(6) * t338 - t322;
t416 = pkin(9) + t423;
t545 = qJD(6) * t416 + t357 + t478;
t542 = g(1) * t480 - g(2) * t481;
t501 = sin(qJ(1));
t505 = cos(qJ(1));
t541 = g(1) * t501 - g(2) * t505;
t464 = t502 * t488;
t336 = qJD(6) * t384 + t353 * t498 - t464;
t538 = -t336 * t392 - t370 * t382;
t537 = -t351 * t416 - t598;
t536 = t339 * t552 + t340 * t523;
t535 = -t395 * t490 + t437 * t488;
t324 = t338 * t502 + t348 * t498;
t532 = t324 * t523 + t337 * t572 + t498 * t620;
t323 = -t338 * t498 + t348 * t502;
t531 = -t323 * t523 + t337 * t573 + (g(1) * t590 + g(2) * t592) * t467;
t530 = t349 + (t498 * t552 - t573) * t615;
t529 = t543 * t467;
t526 = -t371 * t498 - t393 * t572;
t525 = -t371 * t502 + t393 * t573;
t422 = -pkin(3) * t588 + t477 * t601;
t342 = t365 * t601 - t589;
t519 = t342 * t615 - t351 * t469 - t598;
t518 = -qJD(1) * t451 - t448 + t543;
t517 = 0.2e1 * qJD(3) * t451 - qJDD(3) * t470;
t506 = qJD(3) ^ 2;
t514 = -0.2e1 * qJDD(1) * t472 - t470 * t506 + t542;
t512 = g(3) * t485 + t431 * t428 + t486 * t543 - t613;
t511 = -t430 * t428 * MDP(12) - t615 * t523 * MDP(25) + ((t335 - t623) * t502 + (-t336 - t622) * t498) * MDP(22) + (t530 + t595) * MDP(24) + (-t594 - t619) * MDP(23) + (t384 * t550 + t599) * MDP(21) + (t428 * t490 + t378) * MDP(14) + (-t540 + (-qJD(1) * t437 - t430) * t490) * MDP(15) + (-t428 ^ 2 + t430 ^ 2) * MDP(13) + t488 * MDP(16);
t509 = t431 * t430 + t516 + t624;
t489 = -qJ(5) - pkin(8) - pkin(7);
t471 = -pkin(4) * t601 - pkin(5);
t446 = qJDD(3) * t504 - t500 * t506;
t445 = qJDD(3) * t500 + t504 * t506;
t438 = pkin(2) + t580;
t415 = -pkin(5) - t422;
t407 = t468 * t590 + t593;
t406 = -t468 * t591 + t592;
t405 = -t468 * t592 + t591;
t404 = t468 * t593 + t590;
t377 = -t396 * t490 - t436 * t488;
t354 = t380 * t495 - t527 * t601;
t341 = t365 * t495 + t362;
t332 = pkin(5) * t370 - pkin(9) * t371 + t561;
t325 = t346 * t495 - t510 * t601;
t321 = t502 * t322;
t1 = [(-g(1) * t404 - g(2) * t406 - t324 * t370 + t325 * t384 + t354 * t335 + (-(-qJD(6) * t355 + t332) * t615 - t596 + t547 * t392 - qJD(6) * t597) * t498 + t612 * t502) * MDP(27) + (-g(1) * t405 - g(2) * t407 + t321 * t392 + t323 * t370 + t325 * t382 + t354 * t336 + (t332 * t615 + t596 + (-t338 * t392 - t355 * t615 + t597) * qJD(6)) * t502 + t612 * t498) * MDP(26) + (-t525 * t615 + t566 + t585) * MDP(23) + (t526 * t615 + t538 - t567) * MDP(24) + (t351 * t392 + t370 * t615) * MDP(25) + (t379 * t616 + t431 * t396 + t408 * t436 + t428 * t479 + t486 * t542 + t488 * t551 + t490 * t515) * MDP(17) + (t378 * t616 - t431 * t395 + t408 * t437 - t430 * t479 - t485 * t542 - t488 * t581 - t490 * t522) * MDP(18) + (t378 * t437 + t395 * t430) * MDP(12) + (-t378 * t436 - t379 * t437 + t395 * t428 + t396 * t430) * MDP(13) + t377 * MDP(15) + (-t316 * t393 - t317 * t392 + t325 * t523 + t326 * t552 - t339 * t371 - t340 * t370 + t352 * t355 + t353 * t354 - t543) * MDP(19) + (g(1) * t505 + g(2) * t501) * MDP(3) + ((-t382 * t502 - t384 * t498) * t371 + (-t599 - t336 * t502 + (t382 * t498 - t384 * t502) * qJD(6)) * t393) * MDP(22) + (t541 + (t496 ^ 2 + t497 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + 0.2e1 * (t500 * t568 - t570 * t579) * MDP(6) + (t317 * t355 + t340 * t326 - t316 * t354 - t339 * t325 + t513 * t524 + t391 * t561 - g(1) * (-pkin(1) * t501 - t438 * t480 - t481 * t489) - g(2) * (pkin(1) * t505 + t438 * t481 - t480 * t489)) * MDP(20) + (qJDD(1) * t492 + 0.2e1 * t500 * t562) * MDP(5) + t541 * MDP(2) + t535 * MDP(14) + qJDD(1) * MDP(1) + (t500 * t517 + t504 * t514) * MDP(10) + (-t500 * t514 + t504 * t517) * MDP(11) + (t335 * t393 * t502 - t384 * t525) * MDP(21) + t445 * MDP(7) + t446 * MDP(8); t586 * MDP(4) + t446 * MDP(10) - t445 * MDP(11) + t377 * MDP(17) - t535 * MDP(18) + (t352 * t393 + t353 * t392 + t370 * t523 + t371 * t552) * MDP(19) + (-t316 * t392 + t317 * t393 - t339 * t370 + t340 * t371 - g(3)) * MDP(20) + (-t538 - t567) * MDP(26) + (-t566 + t585) * MDP(27) + (MDP(26) * t526 + MDP(27) * t525) * t615; (t317 * t423 + t316 * t422 - t391 * (t478 - t610) - g(3) * t580 - t543 * (-pkin(3) * t500 - pkin(4) * t485) + t582 * t340 - t583 * t339) * MDP(20) + (t584 * t490 + (t430 * t578 - t488 * t499 - t490 * t574) * pkin(3) + t512) * MDP(18) + (t352 * t423 - t353 * t422 + t523 * t583 + t552 * t582 + t536) * MDP(19) + (t415 * t336 - t620 * t502 + t537 * t498 + t583 * t382 + (-t498 * t582 - t502 * t545) * t615 + t531) * MDP(26) + t511 + (-g(3) * t504 + t500 * t518 + t482) * MDP(10) + (-t500 * t586 + t504 * t518) * MDP(11) + (-t553 * t490 + (-t428 * t578 + t488 * t503 - t490 * t575) * pkin(3) + t509) * MDP(17) + qJDD(3) * MDP(9) + MDP(7) * t569 + MDP(8) * t568 + (t415 * t335 + t537 * t502 - t498 * t529 + t583 * t384 + (t498 * t545 - t502 * t582) * t615 + t532) * MDP(27) + (-MDP(5) * t500 * t504 + MDP(6) * t579) * qJD(1) ^ 2; (-t341 * t523 - t342 * t552 + (t352 * t495 - t353 * t601) * pkin(4) + t536) * MDP(19) + (t471 * t336 - t341 * t382 + t519 * t498 + (-t620 - t618) * t502 + t531) * MDP(26) + t511 + (t471 * t335 - t341 * t384 + t519 * t502 + (-t529 + t618) * t498 + t532) * MDP(27) + (t339 * t341 - t340 * t342 + (t316 * t601 + t317 * t495 + t391 * t430 + t624) * pkin(4)) * MDP(20) + (t490 * t554 + t512) * MDP(18) + (-t490 * t534 + t509) * MDP(17); (-t523 ^ 2 - t552 ^ 2) * MDP(19) + (t530 - t595) * MDP(26) + (-t594 + t619) * MDP(27) + (t339 * t523 - t340 * t552 + t513 - t542) * MDP(20); t384 * t382 * MDP(21) + (-t382 ^ 2 + t384 ^ 2) * MDP(22) + (t565 + t623) * MDP(23) + (t464 + t622) * MDP(24) + t351 * MDP(25) + (-g(1) * t406 + g(2) * t404 + t324 * t615 - t337 * t384 + t321) * MDP(26) + (g(1) * t407 - g(2) * t405 + t323 * t615 + t337 * t382) * MDP(27) + ((-t315 + t607) * MDP(27) + (-MDP(24) * t523 - MDP(26) * t338 - MDP(27) * t348) * qJD(6)) * t502 + (-qJD(6) * t523 * MDP(23) + (-qJD(6) * t490 - t353) * MDP(24) + (-t548 + t607) * MDP(26) + t547 * MDP(27)) * t498;];
tau  = t1;
