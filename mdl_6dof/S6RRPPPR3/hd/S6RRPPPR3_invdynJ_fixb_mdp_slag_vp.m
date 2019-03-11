% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:51
% EndTime: 2019-03-09 08:15:59
% DurationCPUTime: 7.66s
% Computational Cost: add. (2915->568), mult. (5831->679), div. (0->0), fcn. (3450->10), ass. (0->247)
t501 = cos(qJ(2));
t593 = qJD(4) * t501;
t498 = sin(qJ(2));
t595 = qJD(2) * t498;
t529 = pkin(7) * t595 + t593;
t588 = qJD(1) * qJD(2);
t570 = t498 * t588;
t435 = qJ(4) * t570;
t585 = qJDD(1) * t501;
t457 = pkin(7) * t585;
t485 = qJDD(2) * qJ(3);
t486 = qJD(2) * qJD(3);
t578 = t457 + t485 + t486;
t557 = -t435 - t578;
t359 = qJ(4) * t585 + qJD(1) * t529 + t557;
t628 = qJDD(2) * pkin(4) + qJDD(5);
t351 = -t359 + t628;
t621 = g(3) * t498;
t499 = sin(qJ(1));
t480 = g(2) * t499;
t502 = cos(qJ(1));
t629 = g(1) * t502 + t480;
t510 = t501 * t629 + t621;
t637 = t351 - t510;
t598 = qJD(1) * t498;
t441 = qJD(6) + t598;
t492 = sin(pkin(9));
t597 = qJD(1) * t501;
t572 = t492 * t597;
t493 = cos(pkin(9));
t592 = t493 * qJD(2);
t414 = t572 - t592;
t500 = cos(qJ(6));
t406 = t500 * t414;
t571 = t493 * t597;
t596 = qJD(2) * t492;
t415 = t571 + t596;
t497 = sin(qJ(6));
t613 = t415 * t497;
t562 = t406 + t613;
t636 = t441 * t562;
t626 = pkin(2) + pkin(3);
t635 = t498 * t626;
t417 = t492 * t497 - t493 * t500;
t408 = t417 * qJD(6);
t524 = t498 * t417;
t634 = qJD(1) * t524 + t408;
t569 = t501 * t588;
t586 = qJDD(1) * t498;
t523 = t569 + t586;
t416 = qJDD(6) + t523;
t418 = t492 * t500 + t493 * t497;
t409 = t418 * qJD(6);
t605 = t418 * t598 + t409;
t633 = t416 * t417 + t441 * t605;
t574 = t626 * qJD(2);
t632 = t626 * qJDD(2);
t460 = pkin(7) * t597;
t573 = qJ(4) * t597;
t424 = t460 - t573;
t487 = qJD(2) * qJ(3);
t411 = -t424 - t487;
t473 = t498 * pkin(4);
t546 = qJ(5) * t501 + t473;
t481 = g(1) * t499;
t622 = g(2) * t502;
t630 = t481 - t622;
t412 = -qJD(1) * pkin(1) - pkin(2) * t597 - qJ(3) * t598;
t468 = t498 * qJ(3);
t476 = t501 * pkin(2);
t601 = t476 + t468;
t426 = -pkin(1) - t601;
t617 = pkin(7) * qJDD(2);
t627 = (qJD(1) * t426 + t412) * qJD(2) - t617;
t625 = pkin(8) * t498;
t624 = pkin(8) * t501;
t623 = g(2) * qJ(4);
t479 = g(3) * t501;
t475 = t501 * pkin(3);
t620 = pkin(7) - qJ(4);
t584 = qJ(5) + t626;
t619 = pkin(8) + t584;
t496 = qJ(3) + pkin(4);
t618 = qJD(2) * pkin(2);
t616 = qJ(3) * t501;
t491 = qJDD(1) * pkin(1);
t615 = qJDD(2) * pkin(2);
t367 = t414 * t497 - t415 * t500;
t614 = t367 * t441;
t488 = t498 ^ 2;
t505 = qJD(1) ^ 2;
t611 = t488 * t505;
t610 = t498 * t499;
t609 = t498 * t502;
t608 = t498 * t505;
t607 = t499 * t501;
t606 = t501 * t502;
t522 = pkin(4) * t501 - t498 * t584;
t508 = qJD(2) * t522 + qJD(5) * t501;
t467 = t498 * qJD(3);
t545 = pkin(2) * t585 + qJ(3) * t523 + qJD(1) * t467 + t491;
t519 = pkin(3) * t585 + qJDD(4) + t545;
t333 = qJD(1) * t508 + qJDD(1) * t546 + t519;
t440 = pkin(7) * t569;
t456 = pkin(7) * t586;
t567 = qJDD(3) + t440 + t456;
t587 = qJD(1) * qJD(4);
t506 = -qJ(4) * t523 - t498 * t587 + t567;
t349 = -qJD(2) * qJD(5) - qJDD(2) * t584 + t506;
t327 = t333 * t492 + t349 * t493;
t386 = pkin(3) * t597 + qJD(4) - t412;
t368 = qJD(1) * t546 + t386;
t451 = qJ(4) * t598;
t582 = pkin(7) * t598;
t551 = qJD(3) + t582;
t531 = -t451 + t551;
t382 = -qJD(2) * t584 + t531;
t339 = t368 * t492 + t382 * t493;
t594 = qJD(2) * t501;
t602 = qJ(3) * t594 + t467;
t360 = t508 + t602;
t404 = -qJD(4) * t498 + t594 * t620;
t342 = t360 * t492 + t404 * t493;
t453 = qJ(3) * t597;
t371 = qJD(1) * t522 + t453;
t348 = t371 * t492 + t424 * t493;
t577 = t475 + t601;
t543 = t577 + t546;
t381 = pkin(1) + t543;
t429 = t620 * t498;
t358 = t381 * t492 + t429 * t493;
t603 = -qJDD(2) * t492 - t493 * t585;
t489 = t501 ^ 2;
t599 = t488 - t489;
t575 = pkin(5) * t492 - pkin(7);
t553 = t575 * t498;
t591 = -qJD(1) * t553 + qJD(3) - t451;
t422 = -t451 + t582;
t590 = qJD(3) + t422;
t589 = -qJD(4) - t386;
t583 = t492 * t625;
t581 = t501 * t608;
t558 = qJDD(2) * t493 - t492 * t585;
t383 = t492 * t570 + t558;
t384 = t493 * t570 + t603;
t580 = qJD(6) * t406 - t497 * t383 + t384 * t500;
t579 = t457 + 0.2e1 * t485 + 0.2e1 * t486;
t576 = -g(1) * t609 - g(2) * t610 + t479;
t568 = t584 * t586;
t566 = -pkin(1) - t468;
t326 = t333 * t493 - t349 * t492;
t322 = pkin(5) * t523 - pkin(8) * t384 + t326;
t325 = -pkin(8) * t383 + t327;
t564 = t322 * t500 - t325 * t497;
t341 = t360 * t493 - t404 * t492;
t338 = t368 * t493 - t382 * t492;
t347 = t371 * t493 - t424 * t492;
t357 = t381 * t493 - t429 * t492;
t563 = t500 * t383 + t497 * t384;
t561 = g(1) * t584;
t413 = pkin(1) + t577;
t560 = qJD(1) * t413 + t386;
t556 = pkin(1) * t502 + pkin(2) * t606 + pkin(7) * t499 + qJ(3) * t609;
t555 = t456 + t576;
t554 = t498 * t574;
t504 = qJD(2) ^ 2;
t552 = pkin(7) * t504 + t622;
t548 = -t418 * t416 + t441 * t634;
t477 = t502 * pkin(7);
t547 = g(1) * (-qJ(4) * t502 + t477);
t544 = pkin(3) * t606 + t556;
t542 = t322 * t497 + t325 * t500;
t541 = -t326 * t492 + t327 * t493;
t332 = pkin(5) * t598 + pkin(8) * t415 + t338;
t334 = pkin(8) * t414 + t339;
t323 = t332 * t500 - t334 * t497;
t324 = t332 * t497 + t334 * t500;
t540 = -t338 * t493 - t339 * t492;
t539 = t338 * t492 - t339 * t493;
t344 = pkin(5) * t498 + t493 * t624 + t357;
t350 = t492 * t624 + t358;
t538 = t344 * t500 - t350 * t497;
t537 = t344 * t497 + t350 * t500;
t425 = t551 - t618;
t427 = t460 + t487;
t536 = t425 * t501 - t427 * t498;
t535 = -qJDD(3) - t555;
t534 = -MDP(19) * t493 + MDP(20) * t492;
t533 = pkin(5) * t501 - t493 * t625;
t532 = t566 - t476;
t530 = -0.2e1 * pkin(1) * t588 - t617;
t528 = t440 - t535;
t421 = t619 * t493;
t527 = -qJD(1) * t533 + qJD(5) * t492 + qJD(6) * t421 - t347;
t420 = t619 * t492;
t526 = -qJD(1) * t583 + qJD(5) * t493 - qJD(6) * t420 + t348;
t525 = -MDP(19) * t492 - MDP(20) * t493 - MDP(17);
t328 = qJD(6) * t613 + t580;
t521 = -t552 + 0.2e1 * t491;
t520 = -t351 * t501 + t629;
t518 = -qJ(4) * qJDD(1) - t629;
t517 = t326 * t493 + t327 * t492 - t622;
t392 = qJD(2) * pkin(4) + qJD(5) - t411;
t515 = t584 * t594 + (qJD(5) - t392) * t498;
t514 = t392 * MDP(22) - MDP(28) * t562 + t367 * MDP(29);
t346 = -qJD(1) * t554 + t519;
t385 = -t554 + t602;
t513 = -qJD(1) * t385 - qJDD(1) * t413 - t346 + t622;
t512 = t528 - t615;
t329 = qJD(6) * t367 + t563;
t363 = pkin(2) * t570 - t545;
t405 = pkin(2) * t595 - t602;
t509 = -qJD(1) * t405 - qJDD(1) * t426 - t363 - t552;
t391 = -pkin(7) * t570 + t578;
t400 = t567 - t615;
t507 = qJD(2) * t536 + t391 * t501 + t400 * t498;
t484 = pkin(9) + qJ(6);
t470 = t501 * qJ(4);
t466 = cos(t484);
t465 = sin(t484);
t454 = qJ(4) * t595;
t446 = g(1) * t607;
t445 = g(1) * t610;
t439 = qJ(3) * t606;
t437 = qJ(3) * t607;
t434 = pkin(5) * t493 + t496;
t430 = pkin(7) * t501 - t470;
t423 = pkin(2) * t598 - t453;
t407 = -t501 * t575 - t470;
t403 = -t454 + t529;
t402 = -t598 * t626 + t453;
t399 = t417 * t501;
t398 = t418 * t501;
t397 = -t574 + t531;
t396 = -t465 * t499 + t466 * t609;
t395 = -t465 * t609 - t466 * t499;
t394 = -t465 * t502 - t466 * t610;
t393 = t465 * t610 - t466 * t502;
t379 = qJD(2) * t553 + t454 - t593;
t364 = -pkin(5) * t414 + t392;
t356 = t506 - t632;
t355 = -t408 * t501 - t418 * t595;
t354 = -qJD(2) * t524 + t409 * t501;
t337 = pkin(5) * t383 + t351;
t336 = -qJD(2) * t583 + t342;
t335 = qJD(2) * t533 + t341;
t1 = [(t498 * t530 + t501 * t521 + t446) * MDP(9) + (qJDD(1) * t488 + 0.2e1 * t498 * t569) * MDP(4) + (t328 * t398 - t329 * t399 + t354 * t562 + t355 * t367) * MDP(24) + (-t329 * t498 + t355 * t441 + t398 * t416 + t562 * t594) * MDP(26) + ((t335 * t500 - t336 * t497) * t441 + t538 * t416 + t564 * t498 + t323 * t594 - t379 * t562 + t407 * t329 - t337 * t398 - t364 * t355 - g(1) * t394 - g(2) * t396 + (-t324 * t498 - t441 * t537) * qJD(6)) * MDP(28) + (t430 * t384 + t403 * t415 + (-qJD(1) * t358 - t339) * t594 + t520 * t493 + (-qJD(1) * t342 - qJDD(1) * t358 + t392 * t592 - t492 * t630 - t327) * t498) * MDP(20) + (t430 * t383 + t403 * t414 + (qJD(1) * t357 + t338) * t594 + t520 * t492 + (qJD(1) * t341 + qJDD(1) * t357 + t392 * t596 + t493 * t630 + t326) * t498) * MDP(19) + t630 * MDP(2) + t629 * MDP(3) + ((-qJD(2) * t397 - qJDD(1) * t430 + t359 + (-qJD(2) * t429 + t403) * qJD(1)) * t501 + (-qJD(2) * t411 - qJDD(1) * t429 - t356 + (qJD(2) * t430 - t404) * qJD(1)) * t498 + t629) * MDP(17) + ((t488 + t489) * qJDD(1) * pkin(7) + t507 - t629) * MDP(12) + (pkin(7) * t507 - g(1) * t477 - g(2) * t556 + t363 * t426 + t412 * t405 - t481 * t532) * MDP(14) + (qJDD(2) * t429 - t446 + (t498 * t560 + t404) * qJD(2) + t513 * t501) * MDP(16) + (qJDD(2) * t430 + t445 + (t501 * t560 - t403) * qJD(2) - t513 * t498) * MDP(15) + (t356 * t429 + t397 * t404 - t359 * t430 + t411 * t403 + t346 * t413 + t386 * t385 - t547 - g(2) * t544 + (-g(1) * (t532 - t475) + t623) * t499) * MDP(18) + qJDD(1) * MDP(1) + (t328 * t399 + t354 * t367) * MDP(23) + (t327 * t358 + t339 * t342 + t326 * t357 + t338 * t341 + t351 * t430 - t392 * t403 - t547 - g(2) * (pkin(4) * t609 + qJ(5) * t606 + t544) + (-g(1) * (t566 - t473) + t623 + t501 * t561) * t499) * MDP(22) + (t498 * t627 + t509 * t501 + t446) * MDP(11) + (t509 * t498 - t501 * t627 + t445) * MDP(13) + (qJDD(2) * t498 + t501 * t504) * MDP(6) + (qJDD(2) * t501 - t498 * t504) * MDP(7) + (-t498 * t521 + t501 * t530 - t445) * MDP(10) + 0.2e1 * (t498 * t585 - t588 * t599) * MDP(5) + (-(t335 * t497 + t336 * t500) * t441 - t537 * t416 - t542 * t498 - t324 * t594 + t379 * t367 + t407 * t328 + t337 * t399 + t364 * t354 - g(1) * t393 - g(2) * t395 + (-t323 * t498 - t441 * t538) * qJD(6)) * MDP(29) + (t416 * t498 + t441 * t594) * MDP(27) + (t328 * t498 + t354 * t441 + t367 * t594 + t399 * t416) * MDP(25) + (t341 * t415 + t342 * t414 - t357 * t384 - t358 * t383 + t501 * t517 + t540 * t595 + t446) * MDP(21); (qJD(2) * t422 + t435 + (-g(3) + (-pkin(7) * qJD(2) - t402) * qJD(1)) * t498 + (qJD(1) * t589 + t518) * t501 + t579) * MDP(15) - MDP(4) * t581 + (t351 * t496 - t339 * t348 - t338 * t347 - g(1) * (pkin(4) * t606 + t439) - g(2) * (pkin(4) * t607 + t437) - g(3) * t543 - t541 * t584 + t590 * t392 + t539 * qJD(5) + (t480 * t584 + t502 * t561) * t498) * MDP(22) + (-t347 * t415 - t348 * t414 + (-qJD(5) * t414 + t338 * t598 + t383 * t584 - t327) * t493 + (qJD(5) * t415 + t339 * t598 - t384 * t584 + t326) * t492 - t576) * MDP(21) + ((qJD(1) * t423 - g(3)) * t498 + (qJD(1) * t412 - t629) * t501 + t579) * MDP(13) + (t621 - t457 + (pkin(1) * t505 + t629) * t501) * MDP(10) + (t492 * t568 + t383 * t496 - t590 * t414 + t637 * t493 + (-t338 * t501 - t347 * t498 + t492 * t515) * qJD(1)) * MDP(19) + (t493 * t568 + t384 * t496 - t590 * t415 - t637 * t492 + (t339 * t501 + t348 * t498 + t493 * t515) * qJD(1)) * MDP(20) + (0.2e1 * t615 + (-t412 * t498 + t423 * t501) * qJD(1) + t535) * MDP(11) + ((-pkin(2) * t498 + t616) * qJDD(1) + ((t427 - t487) * t498 + (qJD(3) - t425 - t618) * t501) * qJD(1)) * MDP(12) + MDP(6) * t586 + qJDD(2) * MDP(8) + t599 * MDP(5) * t505 + (t328 * t417 + t329 * t418 + t367 * t605 + t562 * t634) * MDP(24) + (-(t420 * t497 - t421 * t500) * t416 + t434 * t328 - t337 * t418 + t324 * t597 + (-t497 * t527 + t500 * t526) * t441 + t591 * t367 + t634 * t364 + t510 * t465) * MDP(29) + (-t328 * t418 + t367 * t634) * MDP(23) - t441 * MDP(27) * t597 + MDP(7) * t585 + ((t420 * t500 + t421 * t497) * t416 + t434 * t329 - t337 * t417 - t323 * t597 + (t497 * t526 + t500 * t527) * t441 - t591 * t562 - t605 * t364 - t510 * t466) * MDP(28) + (t391 * qJ(3) + t427 * qJD(3) - t400 * pkin(2) - t412 * t423 - g(1) * (-pkin(2) * t609 + t439) - g(2) * (-pkin(2) * t610 + t437) - g(3) * t601 - t536 * qJD(1) * pkin(7)) * MDP(14) + (pkin(1) * t608 - t555) * MDP(9) + (-qJ(4) * t586 - qJD(2) * t424 - 0.2e1 * t632 + ((-qJ(4) * qJD(2) + t402) * t501 + t589 * t498) * qJD(1) + t528) * MDP(16) + (-t562 * t597 + t633) * MDP(26) + (-g(1) * t439 - g(2) * t437 - g(3) * t577 - t359 * qJ(3) - t356 * t626 - t386 * t402 - t397 * t424 - t590 * t411 + t629 * t635) * MDP(18) + ((-t616 + t635) * qJDD(1) + (t397 - t590 + t574) * t597) * MDP(17) + (-t367 * t597 + t548) * MDP(25); t512 * MDP(14) + (-qJDD(2) * pkin(3) + t512) * MDP(18) + (-t383 * t493 + t384 * t492) * MDP(21) + (t541 + t576) * MDP(22) + t548 * MDP(28) + t633 * MDP(29) + t534 * t611 + (MDP(13) + MDP(15)) * (-t504 - t611) + (-MDP(11) + MDP(16)) * (qJDD(2) + t581) + (-t427 * MDP(14) + (t411 - t573) * MDP(18) + (t414 - t572) * MDP(19) + (t415 - t571) * MDP(20) - t514) * qJD(2) + ((-MDP(18) * qJ(4) + MDP(12) + t525) * qJDD(1) + (t412 * MDP(14) + t589 * MDP(18) + (-t414 * t492 - t415 * t493) * MDP(21) + t540 * MDP(22)) * qJD(1)) * t498; (t519 + t630) * MDP(18) + (-t383 * t492 - t384 * t493) * MDP(21) + (t481 + t517) * MDP(22) + (-MDP(28) * t417 - MDP(29) * t418) * t416 + (-MDP(28) * t605 + MDP(29) * t634) * t441 + (-MDP(17) * t489 + t488 * t525) * t505 + (-t501 * MDP(16) + (MDP(15) - t534) * t498) * qJDD(1) + ((t397 * MDP(18) + (t414 * t493 - t415 * t492) * MDP(21) - t539 * MDP(22) + (-MDP(18) * t626 + 0.2e1 * MDP(16)) * qJD(2)) * t498 + (0.2e1 * qJD(2) * MDP(15) - t411 * MDP(18) + (-t414 + t592) * MDP(19) + (-t415 - t596) * MDP(20) + t514) * t501) * qJD(1); ((-t415 + t596) * t598 + t558) * MDP(19) + ((t414 + t592) * t598 + t603) * MDP(20) + (-t414 ^ 2 - t415 ^ 2) * MDP(21) + (-t338 * t415 - t339 * t414 + (-pkin(7) * t588 - g(3)) * t498 + (t518 - t587) * t501 - t557 + t628) * MDP(22) + (t329 + t614) * MDP(28) + (t328 + t636) * MDP(29); -t367 * t562 * MDP(23) + (t367 ^ 2 - t562 ^ 2) * MDP(24) + (t580 - t636) * MDP(25) + (-t563 + t614) * MDP(26) + t416 * MDP(27) + (-g(1) * t395 + g(2) * t393 + t324 * t441 - t364 * t367 - t465 * t479 + t564) * MDP(28) + (g(1) * t396 - g(2) * t394 + t323 * t441 - t364 * t562 - t466 * t479 - t542) * MDP(29) + (MDP(25) * t613 - MDP(26) * t367 - MDP(28) * t324 - MDP(29) * t323) * qJD(6);];
tau  = t1;
