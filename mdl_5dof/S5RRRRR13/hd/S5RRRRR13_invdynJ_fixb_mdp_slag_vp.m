% Calculate vector of inverse dynamics joint torques for
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:19
% EndTime: 2024-09-27 17:32:24
% DurationCPUTime: 2.23s
% Computational Cost: add. (4486->378), mult. (6854->503), div. (0->0), fcn. (4484->22), ass. (0->228)
t485 = sin(pkin(5));
t492 = cos(qJ(5));
t493 = cos(qJ(4));
t596 = t492 * t493;
t562 = t485 * t596;
t487 = sin(qJ(5));
t488 = sin(qJ(4));
t602 = t487 * t488;
t641 = -t485 * t602 + t562;
t519 = t487 * t493 + t488 * t492;
t640 = t519 * t485;
t478 = qJDD(1) + qJDD(2);
t468 = qJDD(3) + t478;
t480 = qJD(1) + qJD(2);
t469 = qJD(3) + t480;
t627 = qJD(4) + qJD(5);
t501 = t627 * t519;
t337 = -t468 * t562 + t485 * (t468 * t602 + t501 * t469);
t639 = t641 * t469;
t495 = cos(qJ(2));
t625 = pkin(1) * t495;
t465 = qJDD(1) * t625;
t490 = sin(qJ(2));
t617 = pkin(1) * qJD(1);
t571 = t490 * t617;
t404 = pkin(2) * t478 - qJD(2) * t571 + t465;
t434 = pkin(2) * t480 + t495 * t617;
t489 = sin(qJ(3));
t581 = qJD(3) * t490;
t554 = qJD(1) * t581;
t539 = pkin(1) * t554;
t437 = t489 * t539;
t494 = cos(qJ(3));
t584 = qJD(2) * t495;
t555 = qJD(1) * t584;
t574 = qJDD(1) * t490;
t508 = (t555 + t574) * pkin(1);
t637 = -t489 * t404 - (qJD(3) * t434 + t508) * t494 + t437;
t597 = t490 * t494;
t517 = -t489 * t495 - t597;
t411 = t517 * t617;
t582 = qJD(3) * t489;
t570 = pkin(2) * t582;
t534 = t411 + t570;
t395 = -t434 * t489 - t494 * t571;
t579 = qJD(4) * t488;
t633 = t485 * (pkin(4) * t579 + t395);
t394 = t434 * t494 - t489 * t571;
t486 = cos(pkin(5));
t578 = qJD(4) * t493;
t556 = t486 * t578;
t604 = t486 * t488;
t632 = -pkin(3) * t556 + t394 * t493 + t395 * t604;
t598 = t489 * t490;
t516 = t494 * t495 - t598;
t412 = t516 * t617;
t603 = t486 * t493;
t631 = -t411 * t603 + t412 * t488 + (-t488 * t494 - t489 * t603) * qJD(3) * pkin(2);
t463 = pkin(2) * t494 + pkin(3);
t580 = qJD(3) * t494;
t630 = t411 * t604 - t463 * t556 + (-pkin(2) * t580 + t412) * t493;
t464 = pkin(2) + t625;
t440 = t494 * t464;
t409 = -pkin(1) * t598 + pkin(3) + t440;
t392 = t409 * t604;
t475 = t485 * pkin(9);
t591 = pkin(1) * t597 + t464 * t489;
t402 = t475 + t591;
t595 = t402 * t493 + t392;
t432 = t463 * t604;
t438 = pkin(2) * t489 + t475;
t593 = t438 * t493 + t432;
t484 = qJ(1) + qJ(2);
t477 = qJ(3) + t484;
t461 = sin(t477);
t462 = cos(t477);
t629 = g(1) * t462 + g(2) * t461;
t628 = g(1) * t461 - g(2) * t462;
t610 = t469 * t485;
t373 = pkin(9) * t610 - t395;
t609 = t469 * t486;
t436 = qJD(4) + t609;
t575 = qJD(4) - t436;
t626 = g(3) * t485 + t373 * t575;
t536 = pkin(10) * t610 + t373;
t622 = pkin(3) * t469;
t381 = t394 + t622;
t568 = t381 * t604;
t343 = t493 * t536 + t568;
t624 = pkin(2) * t468;
t623 = pkin(3) * t468;
t621 = pkin(4) * t493;
t620 = pkin(10) * t485;
t616 = t343 * t492;
t559 = t490 * t580;
t372 = -t464 * t582 + (qJD(2) * t517 - t559) * pkin(1);
t615 = t372 * t469;
t614 = t395 * t469;
t613 = t468 * t485;
t612 = t468 * t486;
t611 = t468 * t493;
t608 = t469 * t493;
t479 = t485 ^ 2;
t607 = t479 * t493;
t606 = t485 * t488;
t605 = t485 * t493;
t599 = t488 * t493;
t483 = qJ(4) + qJ(5);
t454 = pkin(3) * t604;
t592 = pkin(9) * t605 + t454;
t472 = sin(t484);
t474 = cos(t484);
t590 = g(1) * t474 + g(2) * t472;
t481 = t488 ^ 2;
t589 = -t493 ^ 2 + t481;
t588 = MDP(10) * t479;
t587 = MDP(11) * t479;
t586 = MDP(12) * t485;
t585 = MDP(13) * t485;
t577 = qJD(5) * t487;
t435 = qJDD(4) + t612;
t428 = qJDD(5) + t435;
t416 = t428 * MDP(21);
t576 = t435 * MDP(14);
t470 = pkin(5) + t483;
t447 = cos(t470) / 0.2e1;
t553 = pkin(5) - t483;
t457 = cos(t553);
t573 = t457 / 0.2e1 + t447;
t572 = sin(t470) / 0.2e1;
t569 = t485 * (-pkin(9) - pkin(10));
t566 = t469 * t606;
t348 = pkin(9) * t613 - t637;
t397 = t494 * t404;
t531 = -t434 * t582 + t397;
t498 = (-t489 * t574 + (-t489 * t584 - t559) * qJD(1)) * pkin(1) + t531;
t349 = t498 + t623;
t561 = t348 * t493 + t349 * t604 + t381 * t556;
t371 = t464 * t580 + (qJD(2) * t516 - t489 * t581) * pkin(1);
t560 = t371 * t493 + t372 * t604 + t409 * t556;
t558 = t469 * t579;
t557 = t485 * t579;
t552 = -t381 - t622;
t551 = -t402 - t620;
t550 = -t438 - t620;
t549 = -t371 * t488 + t372 * t603;
t442 = pkin(4) * t557;
t547 = t485 * t534 + t442;
t546 = t435 + t612;
t545 = t436 + t609;
t324 = -t373 * t579 + t561;
t514 = t558 - t611;
t322 = -t514 * t620 + t324;
t370 = t381 * t603;
t342 = -t488 * t536 + t370;
t335 = pkin(4) * t436 + t342;
t544 = qJD(5) * t335 + t322;
t543 = qJD(1) * (-qJD(2) + t480);
t542 = qJD(2) * (-qJD(1) - t480);
t541 = t486 * t570;
t540 = t488 * t569;
t535 = g(1) * t472 - g(2) * t474 + t465;
t532 = sin(t553);
t353 = -t394 * t488 + t395 * t603;
t452 = pkin(10) * t605;
t396 = t452 + t592;
t530 = qJD(5) * t396 + t353 - (t493 * t569 - t454) * qJD(4);
t455 = pkin(3) * t603;
t476 = t486 * pkin(4);
t387 = t455 + t476 + t540;
t529 = -qJD(4) * t540 - qJD(5) * t387 + t632;
t433 = t463 * t603;
t369 = t488 * t550 + t433 + t476;
t528 = -qJD(5) * t369 - (qJD(4) * t550 - t541) * t488 + t630;
t374 = t452 + t593;
t527 = qJD(5) * t374 - (t493 * t550 - t432) * qJD(4) - t631;
t526 = t551 * t488;
t525 = -t614 + t623;
t524 = -t335 * t487 - t616;
t393 = t409 * t603;
t352 = t393 + t476 + t526;
t357 = t452 + t595;
t523 = t352 * t492 - t357 * t487;
t522 = t352 * t487 + t357 * t492;
t521 = t409 * t468 + t615;
t419 = t572 - t532 / 0.2e1;
t473 = cos(t483);
t520 = -t419 * t462 - t461 * t473;
t380 = t419 * t461 - t462 * t473;
t515 = qJD(4) * (-t409 * t469 - t381);
t346 = t349 * t603;
t399 = -t461 * t493 - t462 * t604;
t401 = -t461 * t604 + t462 * t493;
t513 = -g(1) * t399 - g(2) * t401 + (-t488 * t348 + t346 + (-t373 * t493 - t568) * qJD(4)) * t486 + t349 * t607;
t398 = t461 * t488 - t462 * t603;
t400 = -t461 * t603 - t462 * t488;
t511 = -g(1) * t398 - g(2) * t400 - t324 * t486;
t334 = (pkin(4) * t514 - t349) * t485;
t359 = (-pkin(4) * t608 - t381) * t485;
t366 = t501 * t485;
t321 = pkin(4) * t435 + t346 + (-pkin(10) * t613 - t348) * t488 - t343 * qJD(4);
t505 = qJD(5) * t524 + t321 * t492 - t487 * t322;
t510 = -g(1) * t520 + g(2) * t380 - t334 * t641 + t359 * t366 + t486 * t505;
t336 = t468 * t640 + t627 * t639;
t386 = t640 * t469;
t431 = qJD(5) + t436;
t509 = -t386 * t639 * MDP(17) + (-t431 * t639 + t336) * MDP(19) + (t386 * t431 - t337) * MDP(20) + (t386 ^ 2 - t639 ^ 2) * MDP(18) + t416;
t338 = t343 * t577;
t365 = t627 * (t596 - t602) * t485;
t471 = sin(t483);
t375 = t461 * t471 - t462 * t573;
t378 = -t461 * t573 - t462 * t471;
t506 = -g(1) * t375 - g(2) * t378 - (t487 * t321 + t544 * t492 - t338) * t486 + t334 * t640 + t359 * t365;
t504 = (t336 * t641 - t337 * t640 + t365 * t639 - t366 * t386) * MDP(18) + (t365 * t431 + t428 * t640) * MDP(19) + (-t366 * t431 + t428 * t641) * MDP(20) + (t336 * t640 + t365 * t386) * MDP(17) + (t493 * t546 - t545 * t579) * t585 + (t488 * t546 + t545 * t578) * t586 + 0.2e1 * (-qJD(4) * t469 * t589 + t468 * t599) * t587 + (t468 * t481 + 0.2e1 * t493 * t558) * t588 + t468 * MDP(7) + (t336 * MDP(19) - t337 * MDP(20) + t416 + t576) * t486;
t503 = MDP(4) * t478 + t504;
t502 = (-pkin(2) * t469 - t434) * qJD(3) - t508;
t500 = t338 - t359 * t639 + (-t343 * t431 - t321) * t487 - g(1) * t380 - g(2) * t520 - g(3) * (t447 - t457 / 0.2e1);
t499 = t629 + t637;
t497 = -g(1) * t378 + g(2) * t375 - g(3) * (t572 + t532 / 0.2e1) - t359 * t386 + t505;
t496 = cos(qJ(1));
t491 = sin(qJ(1));
t430 = (-pkin(3) - t621) * t485;
t410 = (-t463 - t621) * t485;
t383 = (-t409 - t621) * t485;
t358 = -t372 * t485 + t442;
t329 = (t493 * t551 - t392) * qJD(4) + t549;
t328 = qJD(4) * t526 + t560;
t1 = [((t478 * t495 + t490 * t542) * pkin(1) + t535) * MDP(5) + (-(-t402 * t579 + t560) * t436 - t595 * t435 + (t493 * t515 + (-t349 - t521) * t488) * t479 + t511) * MDP(16) + (((-qJDD(1) - t478) * t490 + t495 * t542) * pkin(1) + t590) * MDP(6) + (t615 + t440 * t468 + (-t494 * t554 + (-t555 + (-qJDD(1) - t468) * t490) * t489) * pkin(1) + t531 + t628) * MDP(8) + (g(1) * t491 - g(2) * t496) * MDP(2) + (g(1) * t496 + g(2) * t491) * MDP(3) + qJDD(1) * MDP(1) + ((-qJD(4) * t595 + t549) * t436 + (-t402 * t488 + t393) * t435 + (t488 * t515 + t493 * t521) * t479 + t513) * MDP(15) + ((-qJD(5) * t522 - t328 * t487 + t329 * t492) * t431 + t523 * t428 - t358 * t639 + t383 * t337 + t510) * MDP(22) + (-(qJD(5) * t523 + t328 * t492 + t329 * t487) * t431 - t522 * t428 + t358 * t386 + t383 * t336 + t506) * MDP(23) + t503 + (-t371 * t469 - t468 * t591 + t499) * MDP(9); (-(t369 * t487 + t374 * t492) * t428 + t410 * t336 + (t487 * t527 + t492 * t528) * t431 + t547 * t386 + t506) * MDP(23) + ((t369 * t492 - t374 * t487) * t428 + t410 * t337 + (t487 * t528 - t492 * t527) * t431 - t547 * t639 + t510) * MDP(22) + (-t593 * t435 + ((qJD(4) * t438 + t541) * t488 + t630) * t436 + ((-t463 * t469 - t381) * t578 + (-t463 * t468 + t469 * t534 - t349) * t488) * t479 + t511) * MDP(16) + ((t495 * t543 - t574) * pkin(1) + t590) * MDP(6) + (pkin(1) * t490 * t543 + t535) * MDP(5) + ((-t438 * t488 + t433) * t435 + (-qJD(4) * t593 + t631) * t436 + (-t381 * t579 + t463 * t611 + (-t463 * t579 - t493 * t534) * t469) * t479 + t513) * MDP(15) + t503 + (-t411 * t469 + t397 + (-t539 + t624) * t494 + t502 * t489 + t628) * MDP(8) + (t412 * t469 + t437 + (-t404 - t624) * t489 + t502 * t494 + t629) * MDP(9); ((t387 * t492 - t396 * t487) * t428 + t430 * t337 + (t487 * t529 - t492 * t530) * t431 - t639 * t633 + t510) * MDP(22) + ((-pkin(9) * t606 + t455) * t435 - t353 * t436 + t525 * t607 + (t479 * t488 * t552 - t436 * t592) * qJD(4) + t513) * MDP(15) + (-t592 * t435 + (pkin(9) * t557 + t632) * t436 + (t552 * t578 + (-t349 - t525) * t488) * t479 + t511) * MDP(16) + (-(t387 * t487 + t396 * t492) * t428 + t430 * t336 + (t487 * t530 + t492 * t529) * t431 + t386 * t633 + t506) * MDP(23) + t504 + (t498 - t614 + t628) * MDP(8) + (t394 * t469 + t499) * MDP(9); (t468 * t488 + t575 * t608) * t586 + (-t469 * t488 * t575 + t611) * t585 + t576 + (-g(1) * t400 + g(2) * t398 + t346 - t626 * t493 + (-t348 + (t469 * t479 - t486 * t575) * t381) * t488) * MDP(15) + (t381 * t469 * t607 + g(1) * t401 - g(2) * t399 + t370 * t436 + t488 * t626 - t561) * MDP(16) + (-(-t342 * t487 - t616) * t431 + (t428 * t492 - t431 * t577 + t566 * t639) * pkin(4) + t497) * MDP(22) + ((t342 * t431 - t544) * t492 + (-qJD(5) * t431 * t492 - t386 * t566 - t428 * t487) * pkin(4) + t500) * MDP(23) + t509 + (t587 * t589 - t588 * t599) * t469 ^ 2; (-t431 * t524 + t497) * MDP(22) + ((-t322 + (-qJD(5) + t431) * t335) * t492 + t500) * MDP(23) + t509;];
tau = t1;
