% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:34
% EndTime: 2021-01-16 02:21:49
% DurationCPUTime: 8.15s
% Computational Cost: add. (3328->511), mult. (7837->656), div. (0->0), fcn. (6185->14), ass. (0->229)
t464 = sin(pkin(11));
t471 = sin(qJ(3));
t565 = qJD(2) * t471;
t467 = cos(pkin(11));
t474 = cos(qJ(3));
t579 = t467 * t474;
t415 = -qJD(2) * t579 + t464 * t565;
t470 = sin(qJ(6));
t473 = cos(qJ(6));
t391 = qJD(3) * t470 - t473 * t415;
t580 = t467 * t471;
t588 = t464 * t474;
t424 = t580 + t588;
t418 = t424 * qJD(2);
t635 = qJD(6) + t418;
t638 = t391 * t635;
t469 = qJ(4) + pkin(8);
t433 = t469 * t471;
t434 = t469 * t474;
t390 = -t433 * t464 + t434 * t467;
t460 = qJ(3) + pkin(11);
t457 = sin(t460);
t468 = cos(pkin(6));
t475 = cos(qJ(2));
t603 = cos(pkin(10));
t533 = t603 * t475;
t465 = sin(pkin(10));
t472 = sin(qJ(2));
t586 = t465 * t472;
t412 = -t468 * t533 + t586;
t534 = t603 * t472;
t585 = t465 * t475;
t414 = t468 * t585 + t534;
t522 = g(1) * t414 + g(2) * t412;
t466 = sin(pkin(6));
t582 = t466 * t475;
t496 = g(3) * t582 - t522;
t490 = t496 * t457;
t637 = -qJDD(3) * t390 + t490;
t567 = qJD(1) * t466;
t546 = t472 * t567;
t563 = qJD(3) * t471;
t636 = pkin(3) * t563 - t546;
t584 = t466 * t472;
t421 = t468 * t474 - t471 * t584;
t528 = t469 * qJD(2) + t546;
t566 = qJD(1) * t468;
t386 = t471 * t566 + t474 * t528;
t379 = t464 * t386;
t385 = -t471 * t528 + t474 * t566;
t352 = t385 * t467 - t379;
t624 = -qJD(5) + t352;
t551 = MDP(12) - MDP(17);
t550 = MDP(13) - MDP(18);
t634 = MDP(14) + MDP(16);
t413 = t468 * t534 + t585;
t458 = cos(t460);
t535 = t466 * t603;
t377 = t413 * t458 - t457 * t535;
t401 = t457 * t468 + t458 * t584;
t411 = t468 * t586 - t533;
t587 = t465 * t466;
t506 = -t411 * t458 + t457 * t587;
t633 = -g(1) * t506 - g(2) * t377 - g(3) * t401;
t632 = t470 * t635;
t536 = qJD(3) * t469;
t408 = qJD(4) * t474 - t471 * t536;
t409 = -qJD(4) * t471 - t474 * t536;
t545 = t475 * t567;
t572 = t408 * t464 - t467 * t409 - t424 * t545;
t625 = -t464 * t471 + t579;
t571 = t408 * t467 + t409 * t464 - t625 * t545;
t507 = t411 * t457 + t458 * t587;
t583 = t466 * t474;
t422 = t468 * t471 + t472 * t583;
t364 = -t467 * t421 + t422 * t464;
t577 = t470 * t475;
t439 = t466 * t577;
t630 = t364 * t473 + t439;
t400 = t457 * t584 - t468 * t458;
t498 = -t413 * t457 - t458 * t535;
t552 = t468 * qJDD(1);
t442 = t474 * t552;
t558 = qJD(1) * qJD(2);
t397 = qJDD(2) * pkin(8) + (qJDD(1) * t472 + t475 * t558) * t466;
t486 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t566 + t397;
t515 = t528 * qJD(3);
t339 = qJDD(3) * pkin(3) - t471 * t486 - t474 * t515 + t442;
t340 = (-t515 + t552) * t471 + t486 * t474;
t575 = -t467 * t339 + t464 * t340;
t628 = -g(1) * t507 - g(2) * t498 + g(3) * t400 - t575;
t420 = qJD(3) * t579 - t464 * t563;
t627 = qJ(5) * t420 + qJD(5) * t424 - t636;
t557 = qJD(2) * qJD(3);
t540 = t471 * t557;
t436 = t464 * t540;
t539 = t474 * t557;
t374 = qJDD(2) * t424 + t467 * t539 - t436;
t626 = -qJ(5) * t374 - qJD(5) * t418;
t607 = pkin(5) * t418;
t623 = t607 - t624;
t622 = pkin(3) * t540 + qJDD(4);
t454 = pkin(3) * t474 + pkin(2);
t405 = -qJD(2) * t454 + qJD(4) - t545;
t482 = -qJ(5) * t418 + t405;
t357 = pkin(4) * t415 + t482;
t621 = t357 * t418 + qJDD(5);
t383 = qJD(3) * pkin(3) + t385;
t581 = t467 * t386;
t348 = t464 * t383 + t581;
t343 = -qJD(3) * qJ(5) - t348;
t608 = pkin(5) * t415;
t332 = -t343 - t608;
t350 = t385 * t464 + t581;
t620 = (t350 - t608) * t635 - t332 * t418;
t389 = t467 * t433 + t434 * t464;
t619 = -qJDD(3) * t389 - t496 * t458;
t476 = qJD(3) ^ 2;
t541 = t472 * t558;
t437 = t466 * t541;
t538 = qJDD(1) * t582;
t517 = t437 - t538;
t618 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t476 + t466 * (-g(3) * t475 + t541) - t517 + t522;
t328 = t464 * t339 + t467 * t340;
t461 = qJDD(3) * qJ(5);
t323 = -qJD(3) * qJD(5) - t328 - t461;
t417 = t424 * qJD(3);
t555 = qJDD(2) * t474;
t441 = t467 * t555;
t556 = qJDD(2) * t471;
t373 = qJD(2) * t417 + t464 * t556 - t441;
t322 = -pkin(5) * t373 - t323;
t611 = pkin(3) * t467;
t453 = -pkin(4) - t611;
t448 = -pkin(9) + t453;
t455 = pkin(3) * t565;
t532 = qJ(5) * t415 + t455;
t613 = pkin(4) + pkin(9);
t615 = t322 + t633 + (-qJD(6) * t448 + t418 * t613 + t532) * t635;
t410 = t418 ^ 2;
t612 = pkin(3) * t464;
t610 = pkin(3) * t471;
t609 = pkin(4) * t373;
t604 = qJD(2) * pkin(2);
t600 = qJDD(3) * pkin(4);
t560 = qJD(6) * t473;
t561 = qJD(6) * t470;
t345 = -qJD(3) * t561 + t473 * qJDD(3) + t470 * t373 + t415 * t560;
t598 = t345 * t473;
t511 = -qJ(5) * t424 - t454;
t358 = -t613 * t625 + t511;
t369 = qJDD(6) + t374;
t596 = t358 * t369;
t594 = t369 * t470;
t593 = t391 * t415;
t393 = qJD(3) * t473 + t415 * t470;
t592 = t393 * t415;
t590 = t625 * t470;
t366 = t473 * t369;
t576 = qJDD(1) - g(3);
t574 = pkin(5) * t420 + t572;
t573 = -pkin(5) * t417 + t571;
t570 = -t411 * t469 - t414 * t454;
t569 = t454 * t582 + t469 * t584;
t462 = t471 ^ 2;
t568 = -t474 ^ 2 + t462;
t564 = qJD(2) * t472;
t549 = g(3) * t584;
t547 = t473 * t582;
t544 = t466 * t564;
t543 = qJD(2) * t582;
t537 = qJDD(5) + t575;
t347 = t383 * t467 - t379;
t531 = -t412 * t454 + t413 * t469;
t530 = t473 * t635;
t529 = t635 * t393;
t527 = MDP(25) * t635;
t526 = qJDD(3) * t470 - t473 * t373;
t525 = t471 * t543;
t524 = t465 * pkin(3) * t583 + t411 * t610;
t523 = g(1) * t411 - g(2) * t413;
t521 = -t417 * t613 + t627;
t520 = -pkin(4) * t417 + t627;
t519 = qJD(5) - t347;
t518 = pkin(4) * t458 + qJ(5) * t457;
t331 = -qJD(3) * t613 + t519 + t607;
t344 = t415 * t613 + t482;
t325 = t331 * t470 + t344 * t473;
t516 = t421 * pkin(3);
t477 = qJD(2) ^ 2;
t514 = qJDD(2) * t475 - t472 * t477;
t510 = -g(1) * t465 + g(2) * t603;
t508 = -t364 * t470 + t547;
t503 = t417 * t470 - t560 * t625;
t500 = t514 * t466;
t497 = t422 * qJD(3);
t431 = -t545 - t604;
t495 = -qJD(2) * t431 - t397 - t523;
t493 = (-t413 * t471 - t474 * t535) * pkin(3);
t492 = t471 * t551 + t474 * t550;
t491 = -t328 - t633;
t489 = t538 - t496;
t485 = qJD(3) * t350 + t628;
t484 = -pkin(8) * qJDD(3) + (t431 + t545 - t604) * qJD(3);
t362 = pkin(5) * t424 + t389;
t483 = t322 * t625 - t332 * t417 + t362 * t369 + t523;
t372 = -qJDD(2) * t454 + t517 + t622;
t481 = t437 - t489 + t622;
t480 = -t497 - t525;
t479 = t372 + t626;
t478 = -t373 * t390 + t374 * t389 - t415 * t571 + t418 * t572 + t523 - t549;
t450 = qJ(5) + t612;
t406 = qJD(3) * t415;
t384 = qJD(3) * t421 + t474 * t543;
t370 = -pkin(4) * t625 + t511;
t365 = t421 * t464 + t422 * t467;
t363 = pkin(5) * t625 + t390;
t359 = pkin(4) * t418 + t532;
t351 = t467 * t384 + t464 * t480;
t349 = t384 * t464 - t467 * t480;
t346 = t393 * qJD(6) + t526;
t342 = -qJD(3) * pkin(4) + t519;
t330 = t479 + t609;
t329 = t373 * t613 + t479;
t326 = t537 - t600;
t324 = t331 * t473 - t344 * t470;
t321 = pkin(5) * t374 - qJDD(3) * t613 + t537;
t320 = t473 * t321;
t1 = [t576 * MDP(1) + MDP(3) * t500 + (-qJDD(2) * t472 - t475 * t477) * t466 * MDP(4) + (t421 * qJDD(3) + t474 * t500 + (-t497 - 0.2e1 * t525) * qJD(3)) * MDP(10) + (-qJD(3) * t384 - qJDD(3) * t422 + (-t471 * t514 - t475 * t539) * t466) * MDP(11) + (t575 * t364 + t328 * t365 - t347 * t349 + t348 * t351 - g(3) + (-t372 * t475 + t405 * t564) * t466) * MDP(15) + (-t323 * t365 + t326 * t364 + t342 * t349 - t343 * t351 - g(3) + (-t330 * t475 + t357 * t564) * t466) * MDP(19) + ((qJD(6) * t508 + t349 * t473 - t470 * t544) * t635 + t630 * t369 + t351 * t391 + t365 * t346) * MDP(25) + (-(qJD(6) * t630 + t349 * t470 + t473 * t544) * t635 + t508 * t369 + t351 * t393 + t365 * t345) * MDP(26) + t634 * (t349 * t418 - t351 * t415 + t364 * t374 - t365 * t373) - t550 * (t466 * (t374 * t475 - t418 * t564) + qJD(3) * t351 + qJDD(3) * t365) + t551 * (t466 * (-t373 * t475 + t415 * t564) - qJD(3) * t349 - qJDD(3) * t364); qJDD(2) * MDP(2) + t489 * MDP(3) + (-t576 * t584 - t523) * MDP(4) + (qJDD(2) * t462 + 0.2e1 * t471 * t539) * MDP(5) + 0.2e1 * (t471 * t555 - t557 * t568) * MDP(6) + (qJDD(3) * t471 + t474 * t476) * MDP(7) + (qJDD(3) * t474 - t471 * t476) * MDP(8) + (t484 * t471 + t474 * t618) * MDP(10) + (-t471 * t618 + t484 * t474) * MDP(11) + (-t415 * t546 - t372 * t625 - t373 * t454 + t405 * t417 + (t415 * t610 - t572) * qJD(3) + t619) * MDP(12) + (-t418 * t546 + t372 * t424 - t374 * t454 + t405 * t420 + (t418 * t610 - t571) * qJD(3) + t637) * MDP(13) + (t328 * t625 - t347 * t420 - t348 * t417 + t424 * t575 + t478) * MDP(14) + (-g(1) * t570 - g(2) * t531 - g(3) * t569 + t328 * t390 - t572 * t347 + t571 * t348 - t372 * t454 + t575 * t389 + t636 * t405) * MDP(15) + (-t323 * t625 + t326 * t424 + t342 * t420 + t343 * t417 + t478) * MDP(16) + (t572 * qJD(3) + t330 * t625 - t357 * t417 - t370 * t373 + t415 * t520 - t619) * MDP(17) + (t571 * qJD(3) - t330 * t424 - t357 * t420 - t370 * t374 + t418 * t520 - t637) * MDP(18) + (t330 * t370 - t323 * t390 + t326 * t389 - g(1) * (-t414 * t518 + t570) - g(2) * (-t412 * t518 + t531) - g(3) * (t518 * t582 + t569) - t520 * t357 - t571 * t343 + t572 * t342) * MDP(19) + (-t345 * t590 + t393 * t503) * MDP(20) + ((-t391 * t470 + t393 * t473) * t417 - (t598 - t346 * t470 + (-t391 * t473 - t393 * t470) * qJD(6)) * t625) * MDP(21) + (t345 * t424 - t369 * t590 + t393 * t420 + t503 * t635) * MDP(22) + (-t625 * t366 - t346 * t424 - t391 * t420 + (t417 * t473 + t561 * t625) * t635) * MDP(23) + (t369 * t424 + t420 * t635) * MDP(24) + (t320 * t424 + t324 * t420 + t363 * t346 + (-t329 * t424 + t457 * t522 - t596) * t470 + t483 * t473 - g(3) * (t457 * t577 + t472 * t473) * t466 + (t521 * t470 + t574 * t473) * t635 + t573 * t391 + ((-t358 * t473 - t362 * t470) * t635 - t325 * t424 - t332 * t590) * qJD(6)) * MDP(25) + (-t325 * t420 + t363 * t345 + t573 * t393 + (-t596 - (qJD(6) * t331 + t329) * t424 - t332 * qJD(6) * t625 + (-qJD(6) * t362 + t521) * t635 - t490) * t473 + (-(-qJD(6) * t344 + t321) * t424 + t549 + (qJD(6) * t358 - t574) * t635 - t483) * t470) * MDP(26); MDP(7) * t556 + MDP(8) * t555 + qJDD(3) * MDP(9) + (-g(3) * t421 + t471 * t495 + t510 * t583 + t442) * MDP(10) + (g(3) * t422 + (-t466 * t510 - t552) * t471 + t495 * t474) * MDP(11) + (-t405 * t418 + (qJDD(3) * t467 - t415 * t565) * pkin(3) + t485) * MDP(12) + (qJD(3) * t352 + t405 * t415 + (-qJDD(3) * t464 - t418 * t565) * pkin(3) + t491) * MDP(13) + ((t348 - t350) * t418 + (-t347 + t352) * t415 + (-t373 * t464 - t374 * t467) * pkin(3)) * MDP(14) + (-g(1) * t524 - g(2) * t493 - g(3) * t516 + t328 * t612 + t347 * t350 - t348 * t352 - t405 * t455 - t575 * t611) * MDP(15) + (-t373 * t450 + t374 * t453 + (-t343 - t350) * t418 + (t342 + t624) * t415) * MDP(16) + (t359 * t415 + (-pkin(4) + t453) * qJDD(3) - t485 + t621) * MDP(17) + (qJDD(3) * t450 - t357 * t415 + t359 * t418 + t461 + (0.2e1 * qJD(5) - t352) * qJD(3) - t491) * MDP(18) + (-t323 * t450 + t326 * t453 - t357 * t359 - t342 * t350 - g(1) * (pkin(4) * t507 + qJ(5) * t506 + t524) - g(2) * (pkin(4) * t498 + t377 * qJ(5) + t493) - g(3) * (-pkin(4) * t400 + qJ(5) * t401 + t516) + t624 * t343) * MDP(19) + (-t470 * t529 + t598) * MDP(20) + ((-t346 - t529) * t473 + (-t345 + t638) * t470) * MDP(21) + (-t632 * t635 + t366 + t592) * MDP(22) + (-t530 * t635 - t593 - t594) * MDP(23) + t635 * t415 * MDP(24) + (t324 * t415 + t332 * t560 + t450 * t346 + t448 * t366 + t391 * t623 + t470 * t615 - t473 * t620) * MDP(25) + (-t325 * t415 - t332 * t561 + t450 * t345 + t393 * t623 - t448 * t594 + t470 * t620 + t473 * t615) * MDP(26) + (-MDP(5) * t471 * t474 + MDP(6) * t568) * t477; -t436 * MDP(13) + (t347 * t418 + t348 * t415 + t481) * MDP(15) + (t406 + t436) * MDP(18) + (-t342 * t418 - t343 * t415 + t481 + t609 + t626) * MDP(19) + (t593 - t594) * MDP(25) + (-t366 + t592) * MDP(26) - t551 * t441 + (MDP(26) * t632 - t473 * t527) * t635 + (t464 * t492 + t550 * t580 + (-MDP(15) - MDP(19)) * t454) * qJDD(2) + (-t415 * MDP(13) + t551 * t418 + (t467 * t492 + t551 * t588) * qJD(2)) * qJD(3) + t634 * (-t415 ^ 2 - t410); (t406 + t374) * MDP(16) + (-t415 * t418 + qJDD(3)) * MDP(17) + (-t410 - t476) * MDP(18) + (t343 * qJD(3) - t600 + t621 - t628) * MDP(19) + (-qJD(3) * t391 + t366) * MDP(25) + (-qJD(3) * t393 - t594) * MDP(26) + (-MDP(26) * t530 - t470 * t527) * t635; -t391 ^ 2 * MDP(21) + (t345 + t638) * MDP(22) - t526 * MDP(23) + t369 * MDP(24) + (-t344 * t560 - t470 * t329 - t331 * t561 + t320 + t325 * t635 - g(1) * (-t414 * t470 - t473 * t507) - g(2) * (-t412 * t470 - t473 * t498) - g(3) * (t400 * t473 + t439)) * MDP(25) + (t344 * t561 - t473 * t329 - t331 * t560 - t470 * t321 + t324 * t635 + t332 * t391 - g(1) * (-t414 * t473 + t470 * t507) - g(2) * (-t412 * t473 + t470 * t498) - g(3) * (-t400 * t470 + t547)) * MDP(26) + (t391 * MDP(20) + (-qJD(6) + t635) * MDP(23) - t332 * MDP(25) + MDP(21) * t393) * t393;];
tau = t1;
