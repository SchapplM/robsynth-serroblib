% Calculate vector of inverse dynamics joint torques for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:12
% EndTime: 2019-12-31 21:02:20
% DurationCPUTime: 6.16s
% Computational Cost: add. (4099->483), mult. (9152->614), div. (0->0), fcn. (6016->10), ass. (0->208)
t479 = sin(qJ(3));
t538 = qJD(3) * t479;
t483 = cos(qJ(2));
t544 = qJD(1) * t483;
t599 = -t479 * t544 + t538;
t480 = sin(qJ(2));
t545 = qJD(1) * t480;
t525 = t479 * t545;
t482 = cos(qJ(3));
t534 = t482 * qJD(2);
t419 = t525 - t534;
t541 = qJD(2) * t479;
t421 = t482 * t545 + t541;
t476 = sin(pkin(8));
t477 = cos(pkin(8));
t367 = t477 * t419 + t421 * t476;
t445 = -qJD(3) + t544;
t598 = t367 * t445;
t415 = t476 * t482 + t477 * t479;
t400 = t415 * qJD(3);
t554 = t415 * t544 - t400;
t536 = qJD(3) * t482;
t572 = t477 * t482;
t553 = t599 * t476 - t477 * t536 + t544 * t572;
t537 = qJD(3) * t480;
t597 = -qJD(1) * t537 + qJDD(2);
t531 = qJDD(1) * t480;
t362 = t479 * ((qJD(3) + t544) * qJD(2) + t531) - t597 * t482;
t503 = -t419 * t476 + t477 * t421;
t596 = t503 ^ 2;
t481 = sin(qJ(1));
t484 = cos(qJ(1));
t508 = g(1) * t484 + g(2) * t481;
t595 = t508 * t480;
t562 = t482 * t483;
t587 = pkin(3) * t480;
t501 = -qJ(4) * t562 + t587;
t509 = pkin(2) * t480 - pkin(7) * t483;
t422 = t509 * qJD(1);
t550 = pkin(6) * t525 + t482 * t422;
t360 = t501 * qJD(1) + t550;
t406 = t479 * t422;
t566 = t480 * t482;
t569 = t479 * t483;
t371 = t406 + (-pkin(6) * t566 - qJ(4) * t569) * qJD(1);
t579 = qJ(4) + pkin(7);
t513 = qJD(3) * t579;
t535 = qJD(4) * t482;
t395 = -t479 * t513 + t535;
t491 = -qJD(4) * t479 - t482 * t513;
t557 = (-t360 + t491) * t477 + (t371 - t395) * t476;
t466 = t483 * qJDD(1);
t532 = qJD(1) * qJD(2);
t594 = -t480 * t532 + t466;
t460 = pkin(6) * t544;
t593 = t599 * pkin(3) - t460;
t458 = pkin(6) * t531;
t516 = t483 * t532;
t397 = -qJDD(2) * pkin(2) + pkin(6) * t516 + t458;
t472 = g(3) * t483;
t489 = -t472 + t595;
t592 = qJD(3) * pkin(7) * t445 - t397 + t489;
t429 = -qJD(2) * pkin(2) + pkin(6) * t545;
t382 = pkin(3) * t419 + qJD(4) + t429;
t327 = pkin(4) * t367 - qJ(5) * t503 + t382;
t473 = qJ(3) + pkin(8);
t462 = sin(t473);
t463 = cos(t473);
t563 = t481 * t483;
t387 = t462 * t563 + t463 * t484;
t559 = t484 * t462;
t389 = -t481 * t463 + t483 * t559;
t361 = qJD(3) * t534 + (t516 + t531) * t482 + t597 * t479;
t423 = t509 * qJD(2);
t427 = -pkin(2) * t483 - pkin(7) * t480 - pkin(1);
t378 = qJD(1) * t423 + t427 * qJDD(1);
t373 = t482 * t378;
t411 = t427 * qJD(1);
t430 = qJD(2) * pkin(7) + t460;
t377 = t411 * t479 + t430 * t482;
t396 = pkin(6) * t594 + qJDD(2) * pkin(7);
t413 = qJDD(3) - t594;
t315 = pkin(3) * t413 - qJ(4) * t361 - t377 * qJD(3) - qJD(4) * t421 - t396 * t479 + t373;
t495 = t479 * t378 + t482 * t396 + t411 * t536 - t430 * t538;
t318 = -qJ(4) * t362 - qJD(4) * t419 + t495;
t306 = t477 * t315 - t476 * t318;
t514 = -qJDD(5) + t306;
t580 = g(3) * t480;
t591 = g(1) * t389 + g(2) * t387 - t327 * t503 + t462 * t580 + t514;
t590 = t508 * t483 + t580;
t588 = -0.2e1 * pkin(1);
t586 = pkin(4) * t413;
t585 = pkin(6) * t479;
t584 = g(1) * t481;
t581 = g(2) * t484;
t353 = -qJ(4) * t419 + t377;
t350 = t477 * t353;
t376 = t482 * t411 - t430 * t479;
t352 = -qJ(4) * t421 + t376;
t325 = t352 * t476 + t350;
t578 = t325 * t503;
t577 = t353 * t476;
t576 = t361 * t479;
t575 = t419 * t445;
t574 = t421 * t445;
t573 = t421 * t482;
t571 = t479 * t480;
t570 = t479 * t481;
t568 = t479 * t484;
t567 = t480 * t481;
t565 = t480 * t484;
t564 = t481 * t482;
t561 = t482 * t484;
t457 = pkin(3) * t482 + pkin(2);
t431 = t483 * t457;
t560 = t483 * t484;
t307 = t476 * t315 + t477 * t318;
t450 = pkin(6) * t562;
t540 = qJD(2) * t480;
t551 = t482 * t423 + t540 * t585;
t336 = -t480 * t535 + t501 * qJD(2) + (-t450 + (qJ(4) * t480 - t427) * t479) * qJD(3) + t551;
t552 = t479 * t423 + t427 * t536;
t341 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t566 + (-qJD(4) * t480 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t483) * t479 + t552;
t314 = t476 * t336 + t477 * t341;
t558 = -pkin(4) * t554 + qJ(5) * t553 - qJD(5) * t415 + t593;
t348 = -pkin(3) * t445 + t352;
t324 = t476 * t348 + t350;
t556 = pkin(4) * t545 - t557;
t335 = t476 * t360 + t477 * t371;
t329 = qJ(5) * t545 + t335;
t357 = t477 * t395 + t476 * t491;
t555 = t357 - t329;
t417 = t482 * t427;
t374 = -qJ(4) * t566 + t417 + (-pkin(3) - t585) * t483;
t548 = t479 * t427 + t450;
t379 = -qJ(4) * t571 + t548;
t343 = t476 * t374 + t477 * t379;
t447 = pkin(3) * t571;
t547 = t480 * pkin(6) + t447;
t474 = t480 ^ 2;
t546 = -t483 ^ 2 + t474;
t543 = qJD(2) * t419;
t542 = qJD(2) * t421;
t539 = qJD(2) * t483;
t326 = t352 * t477 - t577;
t533 = qJD(5) - t326;
t528 = t479 * t560;
t527 = t413 * qJ(5) + t307;
t522 = t479 * t539;
t526 = pkin(3) * t522 + pkin(6) * t539 + t536 * t587;
t523 = t445 * t534;
t521 = t483 * t534;
t520 = t445 * t538;
t519 = t445 * t536;
t518 = t579 * t479;
t331 = t361 * t476 + t477 * t362;
t512 = -qJD(3) * t411 - t396;
t452 = g(1) * t567;
t510 = -g(2) * t565 + t452;
t507 = t430 * t536 - t373;
t506 = pkin(4) * t463 + qJ(5) * t462;
t505 = -pkin(7) * t413 + qJD(3) * t429;
t313 = t336 * t477 - t341 * t476;
t323 = t348 * t477 - t577;
t332 = t361 * t477 - t362 * t476;
t342 = t374 * t477 - t379 * t476;
t414 = t476 * t479 - t572;
t498 = -pkin(6) * qJDD(2) + t532 * t588;
t402 = t479 * t563 + t561;
t497 = t413 * t479 - t519;
t496 = t413 * t482 + t520;
t486 = qJD(1) ^ 2;
t494 = pkin(1) * t486 + t508;
t485 = qJD(2) ^ 2;
t493 = pkin(6) * t485 + qJDD(1) * t588 + t581;
t492 = t484 * pkin(1) + pkin(3) * t570 + t481 * pkin(6) + t457 * t560 + t565 * t579;
t347 = pkin(3) * t362 + qJDD(4) + t397;
t490 = -t579 * t567 + pkin(3) * t568 + t484 * pkin(6) + (-pkin(1) - t431) * t481;
t428 = t579 * t482;
t380 = t428 * t476 + t477 * t518;
t381 = t477 * t428 - t476 * t518;
t487 = -t381 * t331 + t332 * t380 - t357 * t367 - t590;
t308 = pkin(4) * t331 - qJ(5) * t332 - qJD(5) * t503 + t347;
t455 = -pkin(3) * t477 - pkin(4);
t453 = pkin(3) * t476 + qJ(5);
t449 = pkin(3) * t564;
t405 = t482 * t560 + t570;
t404 = -t528 + t564;
t403 = -t481 * t562 + t568;
t393 = -t476 * t571 + t477 * t566;
t392 = t415 * t480;
t390 = t462 * t481 + t463 * t560;
t388 = t463 * t563 - t559;
t364 = pkin(4) * t414 - qJ(5) * t415 - t457;
t355 = t480 * t400 + t476 * t522 - t477 * t521;
t354 = t414 * t537 - t415 * t539;
t349 = pkin(4) * t392 - qJ(5) * t393 + t547;
t340 = pkin(4) * t483 - t342;
t339 = -qJ(5) * t483 + t343;
t328 = pkin(3) * t421 + pkin(4) * t503 + qJ(5) * t367;
t321 = -qJ(5) * t445 + t324;
t320 = pkin(4) * t445 + qJD(5) - t323;
t319 = -pkin(4) * t354 + qJ(5) * t355 - qJD(5) * t393 + t526;
t310 = -pkin(4) * t540 - t313;
t309 = qJ(5) * t540 - qJD(5) * t483 + t314;
t305 = -t514 - t586;
t304 = -qJD(5) * t445 + t527;
t1 = [qJDD(1) * MDP(1) + (-t581 + t584) * MDP(2) + t508 * MDP(3) + (qJDD(1) * t474 + 0.2e1 * t480 * t516) * MDP(4) + 0.2e1 * (t480 * t466 - t546 * t532) * MDP(5) + (qJDD(2) * t480 + t483 * t485) * MDP(6) + (qJDD(2) * t483 - t480 * t485) * MDP(7) + (t498 * t480 + (-t493 + t584) * t483) * MDP(9) + (t493 * t480 + t498 * t483 - t452) * MDP(10) + (t361 * t566 + (-t479 * t537 + t521) * t421) * MDP(11) + ((-t419 * t482 - t421 * t479) * t539 + (-t576 - t362 * t482 + (t419 * t479 - t573) * qJD(3)) * t480) * MDP(12) + ((-t361 - t523) * t483 + (t496 + t542) * t480) * MDP(13) + ((t445 * t541 + t362) * t483 + (-t497 - t543) * t480) * MDP(14) + (-t413 * t483 - t445 * t540) * MDP(15) + (-(-t427 * t538 + t551) * t445 + t417 * t413 - g(1) * t403 - g(2) * t405 + ((t519 + t543) * pkin(6) + (-pkin(6) * t413 + qJD(2) * t429 - t512) * t479 + t507) * t483 + (pkin(6) * t362 + qJD(2) * t376 + t397 * t479 + t429 * t536) * t480) * MDP(16) + (t552 * t445 - t548 * t413 - g(1) * t402 - g(2) * t404 + (t429 * t534 + (-t520 + t542) * pkin(6) + t495) * t483 + (-t429 * t538 - t377 * qJD(2) + t397 * t482 + (t361 - t523) * pkin(6)) * t480) * MDP(17) + (-t306 * t393 - t307 * t392 - t313 * t503 - t314 * t367 + t323 * t355 + t324 * t354 - t331 * t343 - t332 * t342 + t510) * MDP(18) + (-g(1) * t490 - g(2) * t492 + t306 * t342 + t307 * t343 + t323 * t313 + t324 * t314 + t347 * t547 + t382 * t526) * MDP(19) + (g(1) * t388 - g(2) * t390 + t305 * t483 + t308 * t392 + t310 * t445 + t319 * t367 - t320 * t540 - t327 * t354 + t331 * t349 - t340 * t413) * MDP(20) + (-t304 * t392 + t305 * t393 - t309 * t367 + t310 * t503 - t320 * t355 + t321 * t354 - t331 * t339 + t332 * t340 + t510) * MDP(21) + (g(1) * t387 - g(2) * t389 - t304 * t483 - t308 * t393 - t309 * t445 - t319 * t503 + t321 * t540 + t327 * t355 - t332 * t349 + t339 * t413) * MDP(22) + (t304 * t339 + t321 * t309 + t308 * t349 + t327 * t319 + t305 * t340 + t320 * t310 - g(1) * (-pkin(4) * t388 - qJ(5) * t387 + t490) - g(2) * (pkin(4) * t390 + qJ(5) * t389 + t492)) * MDP(23); MDP(6) * t531 + MDP(7) * t466 + qJDD(2) * MDP(8) + (t494 * t480 - t458 - t472) * MDP(9) + (t580 + (-pkin(6) * qJDD(1) + t494) * t483) * MDP(10) + (-t445 * t573 + t576) * MDP(11) + ((t361 + t575) * t482 + (-t362 + t574) * t479) * MDP(12) + ((-t421 * t480 + t445 * t562) * qJD(1) + t497) * MDP(13) + ((t419 * t480 - t445 * t569) * qJD(1) + t496) * MDP(14) + t445 * MDP(15) * t545 + (-pkin(2) * t362 + t550 * t445 + t505 * t479 + (-t376 * t480 + (-pkin(6) * t419 - t429 * t479) * t483) * qJD(1) + t592 * t482) * MDP(16) + (-pkin(2) * t361 - t406 * t445 + t505 * t482 + (-t429 * t562 + t377 * t480 + (-t421 * t483 + t445 * t566) * pkin(6)) * qJD(1) - t592 * t479) * MDP(17) + (-t306 * t415 - t307 * t414 + t553 * t323 + t554 * t324 + t335 * t367 - t503 * t557 + t487) * MDP(18) + (t307 * t381 - t306 * t380 - t347 * t457 - g(3) * (t480 * t579 + t431) + t593 * t382 + (t357 - t335) * t324 + t557 * t323 + t508 * (t457 * t480 - t483 * t579)) * MDP(19) + (t308 * t414 + t320 * t545 - t554 * t327 + t331 * t364 + t558 * t367 - t380 * t413 + t556 * t445 + t489 * t463) * MDP(20) + (-t304 * t414 + t305 * t415 - t553 * t320 + t554 * t321 + t329 * t367 + t503 * t556 + t487) * MDP(21) + (-t308 * t415 - t321 * t545 + t553 * t327 - t332 * t364 + t381 * t413 - t555 * t445 + t489 * t462 - t503 * t558) * MDP(22) + (-g(3) * t431 + t304 * t381 + t305 * t380 + t308 * t364 + t558 * t327 + t555 * t321 + t556 * t320 + (-g(3) * t506 - t508 * t579) * t483 + (-g(3) * t579 + t508 * (t457 + t506)) * t480) * MDP(23) + (-MDP(4) * t480 * t483 + MDP(5) * t546) * t486; t421 * t419 * MDP(11) + (-t419 ^ 2 + t421 ^ 2) * MDP(12) + (t361 - t575) * MDP(13) + (-t362 - t574) * MDP(14) + t413 * MDP(15) + (-g(1) * t404 + g(2) * t402 - t377 * t445 - t421 * t429 + (t512 + t580) * t479 - t507) * MDP(16) + (g(1) * t405 - g(2) * t403 + g(3) * t566 - t376 * t445 + t419 * t429 - t495) * MDP(17) + (t324 * t503 - t578 + (-t331 * t476 - t332 * t477) * pkin(3) + (-t323 + t326) * t367) * MDP(18) + (-g(1) * t449 + t323 * t325 - t324 * t326 + (g(2) * t561 + t306 * t477 + t307 * t476 - t382 * t421 + t479 * t590) * pkin(3)) * MDP(19) + (-t325 * t445 - t328 * t367 + (pkin(4) - t455) * t413 + t591) * MDP(20) + (t321 * t503 - t331 * t453 + t332 * t455 - t578 + (t320 - t533) * t367) * MDP(21) + (-t463 * t580 - g(1) * t390 - g(2) * t388 - t327 * t367 + t328 * t503 + t413 * t453 + (-0.2e1 * qJD(5) + t326) * t445 + t527) * MDP(22) + (t304 * t453 + t305 * t455 - t327 * t328 - t320 * t325 - g(1) * (-pkin(3) * t528 - pkin(4) * t389 + qJ(5) * t390 + t449) - g(2) * (-t402 * pkin(3) - pkin(4) * t387 + qJ(5) * t388) - g(3) * (-t447 + (-pkin(4) * t462 + qJ(5) * t463) * t480) + t533 * t321) * MDP(23); (t323 * t503 + t324 * t367 + t347 + t472) * MDP(19) + (-t445 * t503 + t331) * MDP(20) + (-t332 - t598) * MDP(22) + (-t320 * t503 + t321 * t367 + t308 + t472) * MDP(23) + (-MDP(19) - MDP(23)) * t595 + (MDP(18) + MDP(21)) * (-t367 ^ 2 - t596); (t367 * t503 - t413) * MDP(20) + (t332 - t598) * MDP(21) + (-t445 ^ 2 - t596) * MDP(22) + (t321 * t445 - t586 - t591) * MDP(23);];
tau = t1;
