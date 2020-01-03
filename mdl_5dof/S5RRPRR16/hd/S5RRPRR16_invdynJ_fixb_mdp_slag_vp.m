% Calculate vector of inverse dynamics joint torques for
% S5RRPRR16
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR16_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR16_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR16_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:40
% EndTime: 2019-12-31 20:47:52
% DurationCPUTime: 7.30s
% Computational Cost: add. (3553->526), mult. (8721->705), div. (0->0), fcn. (6627->10), ass. (0->238)
t467 = sin(pkin(5));
t471 = sin(qJ(2));
t561 = qJD(1) * qJD(2);
t544 = t471 * t561;
t529 = t467 * t544;
t475 = cos(qJ(2));
t558 = qJDD(1) * t475;
t628 = -t467 * t558 + t529;
t574 = qJD(1) * t471;
t548 = t467 * t574;
t468 = cos(pkin(5));
t575 = qJD(1) * t468;
t555 = pkin(1) * t575;
t407 = pkin(7) * t548 - t475 * t555;
t562 = qJD(3) + t407;
t455 = qJD(2) + t575;
t470 = sin(qJ(4));
t474 = cos(qJ(4));
t573 = qJD(1) * t475;
t547 = t467 * t573;
t398 = t455 * t470 + t474 * t547;
t394 = qJD(5) + t398;
t434 = qJD(4) + t548;
t473 = cos(qJ(5));
t530 = t470 * t547;
t400 = t455 * t474 - t530;
t469 = sin(qJ(5));
t607 = t400 * t469;
t356 = -t473 * t434 + t607;
t627 = t356 * t434;
t599 = t467 * t471;
t457 = pkin(7) * t599;
t549 = -pkin(1) * t475 - pkin(2);
t367 = pkin(3) * t599 + t457 + (-pkin(8) + t549) * t468;
t477 = -pkin(2) - pkin(8);
t614 = qJ(3) * t471;
t494 = t475 * t477 - pkin(1) - t614;
t383 = t494 * t467;
t626 = t470 * t367 + t474 * t383;
t476 = cos(qJ(1));
t584 = t475 * t476;
t472 = sin(qJ(1));
t591 = t471 * t472;
t412 = -t468 * t584 + t591;
t588 = t472 * t475;
t589 = t471 * t476;
t414 = t468 * t588 + t589;
t597 = t467 * t475;
t625 = g(1) * t414 + g(2) * t412 - g(3) * t597;
t563 = pkin(3) * t548 + t562;
t521 = pkin(2) * t475 + t614;
t509 = -pkin(1) - t521;
t402 = t509 * t467;
t390 = qJD(1) * t402;
t624 = t390 * t548 + qJDD(3);
t543 = t475 * t561;
t559 = qJDD(1) * t471;
t497 = t543 + t559;
t489 = t497 * t467;
t404 = qJDD(4) + t489;
t355 = t455 * t477 + t563;
t364 = qJD(1) * t383;
t328 = t355 * t470 + t364 * t474;
t560 = qJDD(1) * t468;
t454 = qJDD(2) + t560;
t449 = pkin(7) * t547;
t554 = pkin(1) * qJD(2) * t468;
t534 = qJD(1) * t554;
t553 = pkin(1) * t560;
t531 = qJD(2) * t449 + qJDD(1) * t457 + t471 * t534 - t475 * t553;
t512 = qJDD(3) + t531;
t332 = pkin(3) * t489 + t454 * t477 + t512;
t429 = pkin(2) * t529;
t613 = qJ(3) * t475;
t520 = pkin(8) * t471 - t613;
t570 = qJD(3) * t471;
t483 = qJD(2) * t520 - t570;
t339 = t429 + (qJD(1) * t483 + qJDD(1) * t494) * t467;
t482 = -qJD(4) * t328 + t474 * t332 - t339 * t470;
t315 = -pkin(4) * t404 - t482;
t598 = t467 * t472;
t376 = t414 * t474 - t470 * t598;
t410 = t468 * t470 + t474 * t597;
t592 = t470 * t476;
t606 = t412 * t474;
t493 = g(1) * t376 + g(2) * (t467 * t592 + t606) - g(3) * t410;
t623 = t394 * (pkin(4) * t400 + t394 * pkin(9)) + t315 + t493;
t346 = -qJD(4) * t530 + t454 * t470 + (qJD(4) * t455 - t628) * t474;
t344 = qJDD(5) + t346;
t423 = pkin(4) * t470 - pkin(9) * t474 + qJ(3);
t528 = pkin(4) * t474 + pkin(9) * t470;
t622 = t394 * ((-pkin(3) - t528) * t548 - qJD(4) * t528 - t562) - t423 * t344;
t413 = t468 * t589 + t588;
t596 = t467 * t476;
t507 = -t412 * t470 + t474 * t596;
t621 = t413 * t473 + t469 * t507;
t620 = -t413 * t469 + t473 * t507;
t572 = qJD(2) * t471;
t546 = t467 * t572;
t448 = pkin(2) * t546;
t363 = t467 * t483 + t448;
t617 = pkin(1) * t471;
t460 = t468 * t617;
t618 = pkin(3) + pkin(7);
t388 = (t597 * t618 + t460) * qJD(2);
t619 = -qJD(4) * t626 - t363 * t470 + t388 * t474;
t616 = pkin(2) * t454;
t615 = MDP(7) * t467;
t345 = -qJD(4) * t398 + t474 * t454 + t628 * t470;
t565 = qJD(5) * t473;
t550 = t473 * t345 + t469 * t404 + t434 * t565;
t566 = qJD(5) * t469;
t322 = -t400 * t566 + t550;
t612 = t322 * t469;
t611 = t356 * t394;
t358 = t400 * t473 + t434 * t469;
t610 = t358 * t394;
t609 = t398 * t434;
t608 = t400 * t434;
t602 = t434 * t474;
t435 = t454 * qJ(3);
t601 = t454 * t468;
t464 = t467 ^ 2;
t600 = t464 * qJD(1) ^ 2;
t595 = t469 * t344;
t594 = t469 * t477;
t593 = t470 * t404;
t590 = t471 * t473;
t587 = t473 * t344;
t586 = t473 * t477;
t585 = t474 * t322;
t583 = t477 * t404;
t453 = pkin(2) * t548;
t576 = qJD(1) * t467;
t385 = t520 * t576 + t453;
t408 = t471 * t555 + t449;
t387 = pkin(3) * t547 + t408;
t580 = t474 * t385 + t470 * t387;
t579 = pkin(7) * t597 + t460;
t465 = t471 ^ 2;
t578 = -t475 ^ 2 + t465;
t577 = MDP(11) * t467;
t571 = qJD(2) * t475;
t569 = qJD(4) * t470;
t568 = qJD(4) * t474;
t567 = qJD(4) * t477;
t564 = qJD(2) - t455;
t557 = 0.2e1 * t464;
t556 = g(3) * t599;
t552 = t475 * t600;
t551 = t470 * t597;
t401 = -t468 * qJ(3) - t579;
t545 = t467 * t571;
t496 = t470 * t332 + t474 * t339 + t355 * t568 - t364 * t569;
t314 = pkin(9) * t404 + t496;
t436 = t455 * qJD(3);
t532 = t628 * pkin(7) - t471 * t553 - t475 * t534;
t347 = -t435 - t436 + t532;
t334 = -t628 * pkin(3) - t347;
t317 = pkin(4) * t346 - pkin(9) * t345 + t334;
t540 = -t469 * t314 + t473 * t317;
t538 = t345 * t469 - t473 * t404;
t537 = t394 * t473;
t536 = t455 + t575;
t535 = t454 + t560;
t533 = t471 * t552;
t437 = t455 * qJ(3);
t361 = t437 + t387;
t382 = pkin(3) * t597 - t401;
t526 = g(1) * t412 - g(2) * t414;
t415 = -t468 * t591 + t584;
t525 = -g(1) * t415 - g(2) * t413;
t524 = g(1) * t476 + g(2) * t472;
t392 = (t469 * t475 + t470 * t590) * t576;
t522 = -t473 * t569 - t392;
t519 = t473 * t314 + t469 * t317;
t325 = pkin(9) * t434 + t328;
t329 = pkin(4) * t398 - pkin(9) * t400 + t361;
t319 = t325 * t473 + t329 * t469;
t518 = t325 * t469 - t329 * t473;
t338 = pkin(9) * t599 + t626;
t411 = t468 * t474 - t551;
t342 = pkin(4) * t410 - pkin(9) * t411 + t382;
t517 = t338 * t473 + t342 * t469;
t516 = -t338 * t469 + t342 * t473;
t327 = t355 * t474 - t364 * t470;
t514 = t367 * t474 - t383 * t470;
t511 = -0.2e1 * qJD(2) * t390;
t452 = t475 * t554;
t510 = -pkin(7) * t546 + t452;
t508 = -t411 * t469 + t467 * t590;
t375 = t411 * t473 + t469 * t599;
t506 = -t394 * t565 - t595;
t505 = -t394 * t566 + t587;
t502 = t434 * t358;
t501 = -qJ(3) * t571 - t570;
t348 = t429 + (qJD(1) * t501 + qJDD(1) * t509) * t467;
t389 = t467 * t501 + t448;
t498 = qJD(1) * t389 + qJDD(1) * t402 + t348;
t495 = t474 * t363 + t367 * t568 - t383 * t569 + t470 * t388;
t409 = t579 * qJD(2);
t492 = -g(1) * t413 + g(2) * t415 + t409 * t455;
t462 = t468 * qJD(3);
t366 = -t546 * t618 + t452 + t462;
t491 = -t525 + t556;
t488 = t525 - t532;
t487 = t334 - t491;
t324 = -pkin(4) * t434 - t327;
t486 = -pkin(9) * t344 + (t324 + t327) * t394;
t485 = -t531 + t625;
t484 = t467 * (t564 * t573 + t559);
t481 = qJD(5) * t394 * t477 + t491;
t480 = t408 * t455 + t485;
t479 = (pkin(9) * t547 - qJD(5) * t423 + t580) * t394 + t625;
t406 = -qJ(3) * t547 + t453;
t403 = t468 * t549 + t457;
t397 = t474 * t404;
t393 = -t462 - t510;
t391 = t469 * t470 * t548 - t473 * t547;
t384 = -t437 - t408;
t381 = -pkin(2) * t455 + t562;
t377 = t414 * t470 + t474 * t598;
t373 = -qJD(4) * t551 + t468 * t568 - t474 * t546;
t372 = -qJD(4) * t410 + t470 * t546;
t352 = t512 - t616;
t350 = t377 * t473 + t415 * t469;
t349 = -t377 * t469 + t415 * t473;
t340 = -pkin(4) * t547 + t385 * t470 - t387 * t474;
t337 = -pkin(4) * t599 - t514;
t336 = qJD(5) * t508 + t372 * t473 + t469 * t545;
t335 = qJD(5) * t375 + t372 * t469 - t473 * t545;
t326 = pkin(4) * t373 - pkin(9) * t372 + t366;
t323 = qJD(5) * t358 + t538;
t321 = -pkin(4) * t545 - t619;
t320 = pkin(9) * t545 + t495;
t313 = -t319 * qJD(5) + t540;
t312 = -t518 * qJD(5) + t519;
t1 = [t524 * MDP(3) + (-t347 * t468 - t393 * t455 - t401 * t454 + t526) * MDP(13) + (t352 * t468 + t403 * t454 + t492) * MDP(12) + ((-t471 * t498 + t475 * t511) * MDP(13) + (g(1) * t592 - t328 * t571 - t471 * t496) * MDP(21) + (t471 * t511 + t475 * t498) * MDP(12) + (-t346 * t471 - t398 * t571) * MDP(18) + (t327 * t571 + t471 * t482) * MDP(20) + (t345 * t471 + t400 * t571) * MDP(17) + (t471 * t535 + t536 * t571) * MDP(6) + (t404 * t471 + t434 * t571) * MDP(19)) * t467 + (-t373 * t434 - t404 * t410) * MDP(18) + (t372 * t434 + t404 * t411) * MDP(17) + ((qJDD(1) * t465 + 0.2e1 * t471 * t543) * MDP(4) + 0.2e1 * (t471 * t558 - t561 * t578) * MDP(5)) * t464 + ((-qJD(5) * t517 - t320 * t469 + t326 * t473) * t394 + t516 * t344 + t313 * t410 - t518 * t373 + t321 * t356 + t337 * t323 - t315 * t508 + t324 * t335 - g(1) * t620 - g(2) * t350) * MDP(27) + (t322 * t508 - t323 * t375 - t335 * t358 - t336 * t356) * MDP(23) + (-t323 * t410 - t335 * t394 + t344 * t508 - t356 * t373) * MDP(25) + (g(1) * t606 - g(2) * t376 + t334 * t411 + t382 * t345 + t361 * t372 + t366 * t400 - t404 * t626 - t495 * t434) * MDP(21) + (t348 * t402 + t390 * t389 + t347 * t401 + t384 * t393 + t352 * t403 + t381 * t409 - g(1) * (-pkin(1) * t472 - pkin(2) * t413 + pkin(7) * t596 - qJ(3) * t412) - g(2) * (pkin(1) * t476 + pkin(2) * t415 + pkin(7) * t598 + qJ(3) * t414)) * MDP(14) + (g(1) * t472 - g(2) * t476) * MDP(2) + (-pkin(1) * t497 * t557 - t454 * t579 - t455 * t510 + t468 * t532 - t526) * MDP(10) + (-(qJD(5) * t516 + t320 * t473 + t326 * t469) * t394 - t517 * t344 - t312 * t410 - t319 * t373 + t321 * t358 + t337 * t322 + t315 * t375 + t324 * t336 + g(1) * t621 - g(2) * t349) * MDP(28) + ((qJD(2) * t381 - qJDD(1) * t401 - t347 + (qJD(2) * t403 - t393) * qJD(1)) * t475 + (qJD(2) * t384 + qJDD(1) * t403 + t352 + (qJD(2) * t401 + t409) * qJD(1)) * t471 - t524) * t577 + (t322 * t375 + t336 * t358) * MDP(22) + (-t457 * t454 - t531 * t468 + (t475 * t601 + (-t544 + t558) * t557) * pkin(1) - t492) * MDP(9) + (-g(1) * t507 - g(2) * t377 + t334 * t410 + t382 * t346 + t361 * t373 + t366 * t398 + t514 * t404 + t619 * t434) * MDP(20) + (t322 * t410 + t336 * t394 + t344 * t375 + t358 * t373) * MDP(24) + (t344 * t410 + t373 * t394) * MDP(26) + (t345 * t411 + t372 * t400) * MDP(15) + (-t345 * t410 - t346 * t411 - t372 * t398 - t373 * t400) * MDP(16) + MDP(8) * t601 + qJDD(1) * MDP(1) + (t475 * t535 - t536 * t572) * t615; -MDP(4) * t533 + (t322 * t470 + t522 * t394 + (t502 + t505) * t474) * MDP(24) + (-t327 * t547 + qJ(3) * t346 + t563 * t398 + (t583 + (t361 - t387) * t434) * t474 + ((t385 - t567) * t434 + t487) * t470) * MDP(20) + (t600 * t617 + t480) * MDP(9) + (-t324 * t391 - t340 * t356 - t622 * t473 + t479 * t469 + (-t344 * t594 + t313 + (-t324 * t469 + t356 * t477) * qJD(4) - t481 * t473) * t470 + (-t518 * t548 + t324 * t565 + t315 * t469 - t477 * t323 + (-t394 * t594 - t518) * qJD(4)) * t474) * MDP(27) + ((-t346 - t608) * t474 + (-t345 + t609) * t470) * MDP(16) + (t345 * t474 - t470 * t608) * MDP(15) + t454 * MDP(8) + (-t347 * qJ(3) - t352 * pkin(2) - t390 * t406 - t381 * t408 - g(1) * (-pkin(2) * t414 + qJ(3) * t415) - g(2) * (-pkin(2) * t412 + qJ(3) * t413) - g(3) * t521 * t467 - t562 * t384) * MDP(14) + (qJ(3) * t345 + t580 * t434 + t328 * t547 + t563 * t400 + (-t361 * t434 - t583) * t470 + (-t434 * t567 + t487) * t474) * MDP(21) + (t473 * t585 + (-t474 * t566 + t522) * t358) * MDP(22) + (-t323 * t470 + (t469 * t569 + t391) * t394 + (t506 - t627) * t474) * MDP(25) + (-t406 * t547 - t480 - 0.2e1 * t616 + t624) * MDP(12) + (-t324 * t392 - t340 * t358 + t622 * t469 + t479 * t473 + (-t344 * t586 - t312 + (-t324 * t473 + t358 * t477) * qJD(4) + t481 * t469) * t470 + (-t319 * t548 - t324 * t566 + t315 * t473 - t477 * t322 + (-t394 * t586 - t319) * qJD(4)) * t474) * MDP(28) + MDP(6) * t484 - t434 * MDP(19) * t547 + (pkin(1) * t552 - t407 * t455 - t488 + t556) * MDP(10) + (0.2e1 * t435 + t436 + t562 * t455 + (-g(3) * t471 + (t390 * t475 + t406 * t471) * qJD(1)) * t467 + t488) * MDP(13) + (-t434 * t568 - t593 + (t398 * t475 - t471 * t602) * t576) * MDP(18) + (t344 * t470 + t394 * t602) * MDP(26) + (t356 * t392 + t358 * t391 + (t356 * t473 + t358 * t469) * t569 + (-t612 - t323 * t473 + (t356 * t469 - t358 * t473) * qJD(5)) * t474) * MDP(23) + ((-pkin(2) * t471 + t613) * qJDD(1) + ((-qJ(3) * qJD(2) - t384 - t408) * t471 + (-pkin(2) * qJD(2) - t381 + t562) * t475) * qJD(1)) * t577 + (-t434 * t569 + t397 + (-t434 * t470 * t471 - t400 * t475) * t576) * MDP(17) + t578 * MDP(5) * t600 + (-t564 * t574 + t558) * t615; MDP(11) * t484 + (t454 + t533) * MDP(12) + (-t455 ^ 2 - t465 * t600) * MDP(13) + (t384 * t455 - t485 - t616 + t624) * MDP(14) + (-t434 ^ 2 * t470 - t398 * t455 + t397) * MDP(20) + (-t400 * t455 - t434 * t602 - t593) * MDP(21) + (-t474 * t323 + (-t455 * t473 - t469 * t602) * t394 + (t506 + t627) * t470) * MDP(27) + (-t585 + (t455 * t469 - t473 * t602) * t394 + (t502 - t505) * t470) * MDP(28); -t398 ^ 2 * MDP(16) + (t345 + t609) * MDP(17) + (-t346 + t608) * MDP(18) + t404 * MDP(19) + (t328 * t434 + t482 - t493) * MDP(20) + (g(1) * t377 - g(2) * t507 + g(3) * t411 + t327 * t434 + t361 * t398 - t496) * MDP(21) + (t358 * t537 + t612) * MDP(22) + ((t322 - t611) * t473 + (-t323 - t610) * t469) * MDP(23) + (t394 * t537 + t595) * MDP(24) + (-t394 ^ 2 * t469 + t587) * MDP(25) + (-pkin(4) * t323 - t328 * t356 + t486 * t469 - t623 * t473) * MDP(27) + (-pkin(4) * t322 - t328 * t358 + t623 * t469 + t486 * t473) * MDP(28) + (MDP(15) * t398 + t400 * MDP(16) - t361 * MDP(20) - t358 * MDP(24) + t356 * MDP(25) - t394 * MDP(26) + MDP(27) * t518 + t319 * MDP(28)) * t400; t358 * t356 * MDP(22) + (-t356 ^ 2 + t358 ^ 2) * MDP(23) + (t550 + t611) * MDP(24) + (-t538 + t610) * MDP(25) + t344 * MDP(26) + (-g(1) * t349 - g(2) * t621 - g(3) * t508 + t319 * t394 - t324 * t358 + t540) * MDP(27) + (g(1) * t350 - g(2) * t620 + g(3) * t375 + t324 * t356 - t394 * t518 - t519) * MDP(28) + (-MDP(24) * t607 - MDP(25) * t358 - MDP(27) * t319 + MDP(28) * t518) * qJD(5);];
tau = t1;
