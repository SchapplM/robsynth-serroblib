% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:44:08
% EndTime: 2019-03-08 21:44:16
% DurationCPUTime: 5.64s
% Computational Cost: add. (2751->474), mult. (6055->612), div. (0->0), fcn. (4346->10), ass. (0->227)
t464 = sin(pkin(6));
t469 = sin(qJ(2));
t472 = cos(qJ(2));
t549 = qJD(1) * qJD(2);
t529 = t472 * t549;
t465 = cos(pkin(6));
t562 = qJD(3) * t465;
t599 = qJDD(2) * pkin(8);
t634 = t599 + (qJDD(1) * t469 + t529) * t464 + qJD(1) * t562;
t468 = sin(qJ(3));
t548 = qJD(2) * qJD(3);
t527 = t468 * t548;
t471 = cos(qJ(3));
t544 = qJDD(2) * t471;
t629 = t527 - t544;
t611 = pkin(4) + pkin(8);
t610 = pkin(3) * t471;
t633 = pkin(2) + t610;
t467 = sin(qJ(5));
t470 = cos(qJ(5));
t561 = qJD(3) * t467;
t566 = qJD(2) * t471;
t412 = t470 * t566 + t561;
t556 = qJD(5) * t412;
t348 = -t470 * qJDD(3) - t467 * t629 + t556;
t526 = t471 * t548;
t545 = qJDD(2) * t468;
t497 = t526 + t545;
t411 = qJDD(5) + t497;
t531 = t467 * t566;
t559 = qJD(3) * t470;
t414 = -t531 + t559;
t473 = -pkin(3) - pkin(9);
t547 = qJDD(1) * t465;
t571 = qJD(1) * t464;
t537 = t469 * t571;
t421 = qJD(2) * pkin(8) + t537;
t558 = qJD(3) * t471;
t616 = t421 * t558 + t634 * t468;
t491 = -t471 * t547 + t616;
t486 = qJDD(4) + t491;
t330 = pkin(4) * t497 + qJDD(3) * t473 + t486;
t329 = t470 * t330;
t567 = qJD(2) * t468;
t570 = qJD(1) * t471;
t375 = t468 * t421 - t465 * t570;
t628 = qJD(4) + t375;
t620 = pkin(4) * t567 + t628;
t350 = qJD(3) * t473 + t620;
t602 = qJ(4) * t468;
t525 = -pkin(2) - t602;
t408 = t471 * t473 + t525;
t536 = t472 * t571;
t366 = qJD(2) * t408 - t536;
t335 = t350 * t467 + t366 * t470;
t601 = qJ(4) * t471;
t514 = pkin(9) * t468 - t601;
t557 = qJD(4) * t468;
t485 = qJD(3) * t514 - t557;
t530 = t469 * t549;
t589 = t464 * t472;
t513 = -qJDD(1) * t589 + t464 * t530;
t503 = pkin(3) * t527 + t513;
t338 = qJD(2) * t485 + qJDD(2) * t408 + t503;
t483 = -t335 * qJD(5) - t467 * t338 + t329;
t320 = pkin(5) * t411 + qJ(6) * t348 - qJD(6) * t414 + t483;
t327 = -qJ(6) * t412 + t335;
t448 = qJD(5) + t567;
t632 = t327 * t448 + t320;
t349 = -qJD(5) * t531 + qJDD(3) * t467 + (qJD(3) * qJD(5) - t629) * t470;
t554 = qJD(5) * t470;
t540 = -t467 * t330 - t470 * t338 - t350 * t554;
t555 = qJD(5) * t467;
t495 = t366 * t555 + t540;
t321 = -qJ(6) * t349 - qJD(6) * t412 - t495;
t334 = t470 * t350 - t366 * t467;
t326 = -qJ(6) * t414 + t334;
t325 = pkin(5) * t448 + t326;
t631 = -t325 * t448 + t321;
t595 = t414 * t448;
t630 = -t349 + t595;
t543 = MDP(12) * qJDD(2);
t627 = -t602 - t633;
t420 = t611 * t558;
t584 = t470 * t472;
t490 = (-t467 * t469 + t468 * t584) * t464;
t624 = -qJD(1) * t490 + t470 * t420;
t560 = qJD(3) * t468;
t453 = pkin(3) * t560;
t378 = t453 + t485;
t433 = t611 * t468;
t586 = t468 * t472;
t489 = (t467 * t586 + t469 * t470) * t464;
t623 = qJD(1) * t489 - t470 * t378 + t408 * t555 - t467 * t420 - t433 * t554;
t574 = t470 * t408 + t467 * t433;
t588 = t465 * t468;
t376 = qJD(1) * t588 + t471 * t421;
t460 = qJD(3) * qJ(4);
t369 = -t460 - t376;
t604 = cos(pkin(10));
t522 = t604 * t469;
t463 = sin(pkin(10));
t592 = t463 * t472;
t397 = t465 * t522 + t592;
t523 = t464 * t604;
t356 = t397 * t468 + t471 * t523;
t521 = t604 * t472;
t593 = t463 * t469;
t399 = -t465 * t593 + t521;
t590 = t464 * t471;
t358 = t399 * t468 - t463 * t590;
t591 = t464 * t469;
t542 = t468 * t591;
t402 = -t465 * t471 + t542;
t362 = t402 * t470 + t467 * t589;
t396 = -t465 * t521 + t593;
t398 = t465 * t592 + t522;
t622 = -g(1) * (t358 * t470 - t398 * t467) - g(2) * (t356 * t470 - t396 * t467) - g(3) * t362;
t368 = -qJD(3) * pkin(3) + t628;
t395 = t470 * t411;
t621 = -t448 * t555 + t395;
t619 = -MDP(10) + MDP(13);
t618 = MDP(11) - MDP(14);
t425 = t525 - t610;
t569 = qJD(2) * t425;
t377 = -t536 + t569;
t617 = t377 * t567 + qJDD(4);
t365 = pkin(4) * t566 + t376;
t351 = t460 + t365;
t615 = t351 * t448 + t473 * t411;
t517 = g(1) * t398 + g(2) * t396;
t474 = qJD(3) ^ 2;
t609 = pkin(8) * t474;
t614 = 0.2e1 * qJDD(2) * pkin(2) + t464 * (-g(3) * t472 + t530) - t513 + t517 - t609;
t499 = -qJ(4) * t558 - t557;
t546 = qJDD(2) * t425;
t342 = qJD(2) * t499 + t503 + t546;
t394 = t453 + t499;
t488 = -g(3) * t589 + t517;
t613 = qJD(2) * (-t394 + t537) - t342 + t488 - t546 - t609;
t612 = t414 ^ 2;
t607 = g(3) * t469;
t450 = pkin(5) * t470 + pkin(4);
t606 = pkin(8) + t450;
t605 = qJD(2) * pkin(2);
t603 = pkin(8) * qJDD(3);
t598 = qJDD(3) * pkin(3);
t597 = t348 * t470;
t596 = t412 * t448;
t594 = t448 * t470;
t587 = t467 * t468;
t585 = t470 * t471;
t582 = qJ(6) - t473;
t581 = qJDD(1) - g(3);
t519 = qJ(6) * t471 - t408;
t551 = qJD(6) * t471;
t580 = pkin(5) * t558 + t519 * t554 + (-qJ(6) * t560 - qJD(5) * t433 - t378 + t551) * t467 + t624;
t553 = qJD(5) * t471;
t532 = t467 * t553;
t579 = -t470 * t551 + (t468 * t559 + t532) * qJ(6) - t623;
t578 = -t326 + t325;
t454 = pkin(3) * t567;
t389 = qJD(2) * t514 + t454;
t577 = t467 * t365 + t470 * t389;
t353 = t470 * t365;
t576 = -qJD(6) * t470 + t555 * t582 + t389 * t467 - t353 - (pkin(5) * t471 - qJ(6) * t587) * qJD(2);
t424 = t582 * t470;
t575 = -qJ(6) * t470 * t567 - qJD(5) * t424 - qJD(6) * t467 - t577;
t434 = t611 * t471;
t461 = t468 ^ 2;
t462 = t471 ^ 2;
t573 = t461 - t462;
t572 = t461 + t462;
t568 = qJD(2) * t464;
t565 = qJD(3) * t376;
t564 = qJD(3) * t412;
t563 = qJD(3) * t414;
t552 = qJD(5) * t473;
t475 = qJD(2) ^ 2;
t541 = t468 * t471 * t475;
t539 = t627 * t396;
t538 = t627 * t398;
t535 = t469 * t568;
t534 = t472 * t568;
t524 = pkin(5) * t467 + qJ(4);
t518 = t421 * t560 - t468 * t547 - t634 * t471;
t516 = g(1) * t399 + g(2) * t397;
t511 = g(3) * (t464 * qJ(4) * t586 + pkin(8) * t591 + t589 * t633);
t510 = t448 ^ 2;
t466 = -qJ(6) - pkin(9);
t504 = pkin(5) * t587 - t466 * t471;
t458 = qJDD(3) * qJ(4);
t459 = qJD(3) * qJD(4);
t339 = -t458 - t459 + t518;
t363 = -t402 * t467 + t464 * t584;
t403 = t469 * t590 + t588;
t500 = -t411 * t467 - t448 * t554;
t493 = g(1) * t358 + g(2) * t356 + g(3) * t402;
t357 = t397 * t471 - t468 * t523;
t359 = t463 * t464 * t468 + t399 * t471;
t492 = -g(1) * t359 - g(2) * t357 - g(3) * t403;
t331 = -pkin(4) * t629 - t339;
t487 = t331 + t492;
t422 = -t536 - t605;
t482 = -t603 + (t422 + t536 - t605) * qJD(3);
t481 = t603 + (-t377 - t536 - t569) * qJD(3);
t480 = qJD(3) * t375 + t492 - t518;
t479 = -t491 + t493;
t322 = pkin(5) * t349 + qJDD(6) + t331;
t340 = t486 - t598;
t476 = -t339 * t471 + t340 * t468 + (t368 * t471 + t369 * t468) * qJD(3) - t516;
t423 = t582 * t467;
t419 = t611 * t560;
t418 = -qJ(4) * t566 + t454;
t417 = t470 * t433;
t410 = t412 ^ 2;
t393 = t402 * pkin(3);
t361 = qJD(3) * t403 + t468 * t534;
t360 = -qJD(3) * t542 + (t534 + t562) * t471;
t355 = t358 * pkin(3);
t354 = t356 * pkin(3);
t347 = -qJ(6) * t585 + t574;
t344 = pkin(5) * t468 + t467 * t519 + t417;
t343 = pkin(5) * t412 + qJD(6) + t351;
t333 = qJD(5) * t362 + t361 * t467 + t470 * t535;
t332 = qJD(5) * t363 + t361 * t470 - t467 * t535;
t1 = [t581 * MDP(1) + (-t339 * t403 + t340 * t402 - t360 * t369 + t361 * t368 - g(3)) * MDP(15) + (t332 * t448 + t349 * t403 + t360 * t412 + t362 * t411) * MDP(21) + (-t333 * t448 - t348 * t403 + t360 * t414 + t363 * t411) * MDP(22) + (-t332 * t414 - t333 * t412 + t348 * t362 + t349 * t363) * MDP(23) + (t320 * t362 - t321 * t363 + t322 * t403 + t325 * t332 + t327 * t333 + t343 * t360 - g(3)) * MDP(24) + (t402 * t468 + t403 * t471) * t543 + (t360 * t471 + t361 * t468 + (t402 * t471 - t403 * t468) * qJD(3)) * MDP(12) * qJD(2) + ((qJD(2) * t377 * MDP(15) - qJDD(2) * MDP(4) + (t468 * t618 + t471 * t619 - MDP(3)) * t475) * t469 + (-t342 * MDP(15) + qJDD(2) * MDP(3) - t475 * MDP(4) - t618 * t497 + t619 * t629) * t472) * t464 - t618 * (qJD(3) * t360 + qJDD(3) * t403) + t619 * (qJD(3) * t361 + qJDD(3) * t402); qJDD(2) * MDP(2) + (t581 * t589 + t517) * MDP(3) + (-t581 * t591 + t516) * MDP(4) + (qJDD(2) * t461 + 0.2e1 * t468 * t526) * MDP(5) + 0.2e1 * (t468 * t544 - t548 * t573) * MDP(6) + (qJDD(3) * t468 + t471 * t474) * MDP(7) + (qJDD(3) * t471 - t468 * t474) * MDP(8) + (t482 * t468 + t471 * t614) * MDP(10) + (-t468 * t614 + t482 * t471) * MDP(11) + (t572 * t599 + (-t529 * t572 - t607) * t464 + t476) * MDP(12) + (t481 * t468 - t471 * t613) * MDP(13) + (t468 * t613 + t481 * t471) * MDP(14) + (t342 * t425 + t377 * t394 - g(1) * t538 - g(2) * t539 - t511 + (-t377 * t469 + (-t368 * t468 + t369 * t471) * t472) * t571 + t476 * pkin(8)) * MDP(15) + (t348 * t467 * t471 + (t467 * t560 - t470 * t553) * t414) * MDP(16) + ((-t412 * t467 + t414 * t470) * t560 + (t597 + t349 * t467 + (t412 * t470 + t414 * t467) * qJD(5)) * t471) * MDP(17) + ((t448 * t561 - t348) * t468 + (t500 + t563) * t471) * MDP(18) + ((t448 * t559 - t349) * t468 + (-t564 - t621) * t471) * MDP(19) + (t411 * t468 + t448 * t558) * MDP(20) + ((-t408 * t467 + t417) * t411 - t419 * t412 + t434 * t349 - t516 * t470 + (-t351 * t559 + t329 + (-t338 + t517) * t467) * t468 - g(3) * t489 + (-t378 * t467 + t624) * t448 + (-t335 * t468 - t448 * t574) * qJD(5) + (qJD(3) * t334 + t331 * t470 - t351 * t555 - t412 * t536) * t471) * MDP(21) + (-t574 * t411 - t419 * t414 - t434 * t348 + t516 * t467 + (t517 * t470 + (qJD(3) * t351 + qJD(5) * t366) * t467 + t540) * t468 - g(3) * t490 + t623 * t448 + (-qJD(3) * t335 - t331 * t467 - t351 * t554 - t414 * t536) * t471) * MDP(22) + (t344 * t348 - t347 * t349 - t580 * t414 - t579 * t412 + (-t325 * t467 + t327 * t470) * t560 + (t320 * t467 - t321 * t470 + (t325 * t470 + t327 * t467) * qJD(5) + t488) * t471) * MDP(23) + (t321 * t347 + t320 * t344 + t322 * (pkin(5) * t585 + t434) - t343 * pkin(5) * t532 - g(1) * (-t398 * t504 + t399 * t606 + t538) - g(2) * (-t396 * t504 + t397 * t606 + t539) - t511 + t579 * t327 + t580 * t325 - t343 * t606 * t560 + (-t450 * t607 + (-g(3) * t504 - t343 * t570) * t472) * t464) * MDP(24); -MDP(5) * t541 + t573 * t475 * MDP(6) + MDP(7) * t545 + MDP(8) * t544 + qJDD(3) * MDP(9) + (-t422 * t567 + t479 + t565) * MDP(10) + (-t422 * t566 - t480) * MDP(11) + (-pkin(3) * t468 + t601) * t543 + (-0.2e1 * t598 - t565 + (-qJD(2) * t418 - t547) * t471 - t493 + t616 + t617) * MDP(13) + (0.2e1 * t458 + 0.2e1 * t459 + (t377 * t471 + t418 * t468) * qJD(2) + t480) * MDP(14) + (-t339 * qJ(4) - t340 * pkin(3) - t377 * t418 - t368 * t376 - g(1) * (qJ(4) * t359 - t355) - g(2) * (qJ(4) * t357 - t354) - g(3) * (qJ(4) * t403 - t393) - t628 * t369) * MDP(15) + (-t467 * t595 - t597) * MDP(16) + ((-t349 - t595) * t470 + (t348 + t596) * t467) * MDP(17) + ((-t414 * t471 - t448 * t587) * qJD(2) + t621) * MDP(18) + ((t412 * t471 - t468 * t594) * qJD(2) + t500) * MDP(19) - t448 * MDP(20) * t566 + (-t334 * t566 + qJ(4) * t349 - t353 * t448 + t620 * t412 + t615 * t470 + ((t389 - t552) * t448 + t487) * t467) * MDP(21) + (-qJ(4) * t348 + t577 * t448 + t335 * t566 + t620 * t414 - t615 * t467 + (-t448 * t552 + t487) * t470) * MDP(22) + (-t348 * t424 + t349 * t423 - t575 * t412 - t576 * t414 - t467 * t631 - t470 * t632 + t493) * MDP(23) + (-t321 * t423 - t320 * t424 + t322 * t524 - g(1) * (t358 * t466 + t359 * t524 - t355) - g(2) * (t356 * t466 + t357 * t524 - t354) - g(3) * (t402 * t466 + t403 * t524 - t393) + (pkin(5) * t594 + t620) * t343 + t575 * t327 + t576 * t325) * MDP(24); t468 * t543 + (qJDD(3) + t541) * MDP(13) + (-t461 * t475 - t474) * MDP(14) + (qJD(3) * t369 - t479 - t598 + t617) * MDP(15) + (t395 - t564) * MDP(21) - MDP(22) * t563 + (-qJD(3) * t343 - t493) * MDP(24) + ((-t412 * t567 + t348 - t556) * MDP(23) + t632 * MDP(24) - MDP(22) * t510) * t470 + (-MDP(21) * t510 - t411 * MDP(22) + MDP(23) * t630 + MDP(24) * t631) * t467; t414 * t412 * MDP(16) + (-t410 + t612) * MDP(17) + (t596 - t348) * MDP(18) + t630 * MDP(19) + t411 * MDP(20) + (t335 * t448 - t351 * t414 + t483 + t622) * MDP(21) + (t334 * t448 + t351 * t412 - g(1) * (-t358 * t467 - t398 * t470) - g(2) * (-t356 * t467 - t396 * t470) - g(3) * t363 + t495) * MDP(22) + (pkin(5) * t348 - t412 * t578) * MDP(23) + (t578 * t327 + (-t343 * t414 + t320 + t622) * pkin(5)) * MDP(24); (-t410 - t612) * MDP(23) + (t325 * t414 + t327 * t412 + t322 + t492) * MDP(24);];
tau  = t1;
