% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:11
% EndTime: 2019-03-09 22:04:22
% DurationCPUTime: 6.77s
% Computational Cost: add. (7904->469), mult. (19937->594), div. (0->0), fcn. (14831->8), ass. (0->211)
t505 = sin(qJ(2));
t507 = cos(qJ(2));
t619 = sin(qJ(3));
t572 = qJD(1) * t619;
t621 = cos(qJ(3));
t573 = qJD(1) * t621;
t462 = -t505 * t572 + t507 * t573;
t463 = -t505 * t573 - t507 * t572;
t504 = sin(qJ(4));
t620 = cos(qJ(4));
t544 = t462 * t504 - t463 * t620;
t436 = t462 * t620 + t463 * t504;
t611 = qJ(5) * t436;
t401 = pkin(4) * t544 - t611;
t628 = qJD(6) + t544;
t506 = cos(qJ(6));
t500 = qJD(2) + qJD(3);
t534 = t505 * t619 - t507 * t621;
t518 = t500 * t534;
t513 = t620 * t518;
t474 = t505 * t621 + t507 * t619;
t520 = t500 * t474;
t517 = t520 * qJD(1);
t571 = qJD(4) * t620;
t584 = qJD(4) * t504;
t538 = -qJD(1) * t513 + t462 * t571 + t463 * t584 - t504 * t517;
t369 = t506 * t538;
t503 = sin(qJ(6));
t639 = t628 * t503;
t644 = -t628 * t639 + t369;
t607 = t538 * t503;
t638 = t628 * t506;
t643 = -t628 * t638 - t607;
t498 = qJD(4) + t500;
t583 = qJD(6) * t503;
t514 = t620 * t520;
t516 = t518 * qJD(1);
t512 = qJD(1) * t514 - t504 * t516;
t376 = qJD(4) * t544 + t512;
t582 = qJD(6) * t506;
t596 = t376 * t503 - t436 * t582;
t348 = -t498 * t583 + t596;
t347 = t348 * t506;
t370 = t506 * t376;
t417 = -t436 * t503 + t498 * t506;
t349 = qJD(6) * t417 - t370;
t528 = -t436 * t498 + t538;
t603 = t498 * t503;
t415 = t436 * t506 + t603;
t640 = t415 * t628;
t642 = (-t462 * t584 + t463 * t571 + t498 * t544 - t512) * MDP(21) + t528 * MDP(20) + (-t417 * t639 + t347) * MDP(29) + (-t417 * t638 + (-t348 + t640) * t503 - t506 * t349) * MDP(30);
t641 = t436 * pkin(5);
t459 = t463 * pkin(9);
t622 = pkin(7) + pkin(8);
t482 = t622 * t507;
t477 = qJD(1) * t482;
t464 = t619 * t477;
t481 = t622 * t505;
t475 = qJD(1) * t481;
t612 = qJD(2) * pkin(2);
t470 = -t475 + t612;
t548 = t470 * t621 - t464;
t413 = t459 + t548;
t408 = pkin(3) * t500 + t413;
t468 = t621 * t477;
t536 = t470 * t619 + t468;
t613 = t462 * pkin(9);
t414 = t536 + t613;
t412 = t620 * t414;
t378 = t408 * t504 + t412;
t371 = -qJ(5) * t498 - t378;
t354 = -t371 + t641;
t565 = t628 * t354;
t550 = t475 * t619 - t468;
t419 = t550 - t613;
t591 = -t475 * t621 - t464;
t420 = t459 + t591;
t494 = pkin(2) * t621 + pkin(3);
t569 = t619 * qJD(3);
t570 = t621 * qJD(3);
t637 = t504 * t419 + t420 * t620 - t494 * t571 + ((qJD(4) * t619 + t569) * t504 - t620 * t570) * pkin(2);
t411 = t504 * t414;
t377 = -t408 * t620 + t411;
t580 = qJD(5) + t377;
t624 = t544 ^ 2;
t577 = qJD(1) * qJD(2);
t634 = -0.2e1 * t577;
t615 = t544 * pkin(5);
t632 = MDP(5) * (t505 ^ 2 - t507 ^ 2);
t631 = t354 * t544;
t556 = t620 * t619;
t594 = -t419 * t620 + t504 * t420 - t494 * t584 - (qJD(4) * t556 + (t504 * t621 + t556) * qJD(3)) * pkin(2);
t630 = t498 * t594;
t385 = t413 * t504 + t412;
t552 = pkin(3) * t584 - t385;
t386 = t413 * t620 - t411;
t595 = -pkin(3) * t571 - qJD(5) + t386;
t593 = -qJD(5) + t637;
t629 = -MDP(18) * t544 - MDP(33) * t628;
t581 = t615 + t580;
t574 = qJD(2) * t622;
t555 = qJD(1) * t574;
t471 = t505 * t555;
t472 = t507 * t555;
t525 = t470 * t570 - t471 * t621 - t472 * t619 - t477 * t569;
t383 = -pkin(9) * t517 + t525;
t551 = t471 * t619 - t472 * t621;
t384 = pkin(9) * t516 - t470 * t569 - t477 * t570 + t551;
t335 = t383 * t504 - t384 * t620 + t408 * t584 + t414 * t571;
t495 = -t507 * pkin(2) - pkin(1);
t480 = qJD(1) * t495;
t447 = -pkin(3) * t462 + t480;
t526 = -qJ(5) * t544 + t447;
t387 = -pkin(4) * t436 + t526;
t627 = t387 * t544 + t335;
t540 = -t447 * t544 - t335;
t623 = pkin(4) + pkin(10);
t617 = pkin(3) * t463;
t476 = t505 * t574;
t478 = t507 * t574;
t530 = -t476 * t621 - t478 * t619 - t481 * t570 - t482 * t569;
t396 = -pkin(9) * t520 + t530;
t549 = t476 * t619 - t478 * t621;
t397 = pkin(9) * t518 + t481 * t569 - t482 * t570 + t549;
t426 = -pkin(9) * t474 - t481 * t621 - t482 * t619;
t535 = t481 * t619 - t482 * t621;
t427 = -pkin(9) * t534 - t535;
t343 = -t396 * t620 - t397 * t504 - t426 * t571 + t427 * t584;
t610 = t343 * t498;
t545 = t426 * t504 + t427 * t620;
t344 = qJD(4) * t545 + t396 * t504 - t397 * t620;
t609 = t344 * t498;
t527 = t620 * t534;
t445 = t474 * t504 + t527;
t529 = t504 * t534;
t446 = t474 * t620 - t529;
t450 = pkin(3) * t534 + t495;
t524 = -t446 * qJ(5) + t450;
t366 = t445 * t623 + t524;
t608 = t366 * t538;
t606 = t378 * t498;
t604 = t445 * t503;
t509 = qJD(2) ^ 2;
t602 = t505 * t509;
t601 = t507 * t509;
t510 = qJD(1) ^ 2;
t600 = t507 * t510;
t599 = t623 * t538;
t557 = t383 * t620 + t384 * t504 + t408 * t571 - t414 * t584;
t334 = -qJD(5) * t498 - t557;
t327 = -pkin(5) * t376 - t334;
t598 = t327 * t503 + t354 * t582;
t597 = t641 + t594;
t592 = t615 - t593;
t589 = t615 - t595;
t585 = qJD(1) * t505;
t576 = t619 * pkin(2);
t497 = t505 * t612;
t496 = pkin(2) * t585;
t568 = t505 * t577;
t567 = pkin(1) * t634;
t326 = t327 * t506;
t350 = -t498 * t623 + t581;
t361 = -t436 * t623 + t526;
t341 = t350 * t503 + t361 * t506;
t566 = t341 * t436 + t326;
t393 = -t617 + t401;
t390 = t393 + t496;
t428 = t544 * pkin(10);
t457 = -t494 * t620 + t504 * t576 - pkin(4);
t453 = -pkin(10) + t457;
t560 = -qJD(6) * t453 + t390 + t428;
t493 = -pkin(3) * t620 - pkin(4);
t491 = -pkin(10) + t493;
t559 = -qJD(6) * t491 + t393 + t428;
t558 = t623 * t628 - t611;
t488 = pkin(2) * t568;
t553 = -t641 + t552;
t399 = -t426 * t620 + t427 * t504;
t340 = t350 * t506 - t361 * t503;
t367 = -pkin(4) * t498 + t580;
t547 = -t367 * t436 - t371 * t544;
t546 = -t340 * t436 + t506 * t631 + t598;
t392 = -qJD(4) * t529 + t474 * t571 - t504 * t518 + t514;
t543 = t392 * t503 + t445 * t582;
t539 = -t447 * t436 - t557;
t537 = t480 * t463 + t551;
t533 = t387 * t436 - t334;
t372 = pkin(5) * t446 + t399;
t531 = t327 * t445 + t354 * t392 - t372 * t538;
t523 = -t480 * t462 - t525;
t515 = t463 * t462 * MDP(11) + (-t417 * t436 + t644) * MDP(31) + (t415 * t436 + t643) * MDP(32) + (-t436 ^ 2 + t624) * MDP(19) + (-t462 * t500 - t516) * MDP(13) + (-t463 * t500 - t517) * MDP(14) + (-t462 ^ 2 + t463 ^ 2) * MDP(12) + t629 * t436 + t642;
t439 = pkin(3) * t520 + t497;
t421 = pkin(3) * t517 + t488;
t391 = qJD(4) * t527 + t474 * t584 + t504 * t520 + t513;
t342 = pkin(4) * t392 + t391 * qJ(5) - t446 * qJD(5) + t439;
t339 = t376 * pkin(4) - qJ(5) * t538 - qJD(5) * t544 + t421;
t492 = pkin(3) * t504 + qJ(5);
t456 = pkin(2) * t556 + t494 * t504 + qJ(5);
t448 = t496 - t617;
t398 = t445 * pkin(4) + t524;
t373 = -pkin(5) * t445 + t545;
t358 = t538 * t446;
t357 = t378 + t641;
t333 = -pkin(5) * t391 + t344;
t332 = -pkin(5) * t392 - t343;
t331 = pkin(10) * t392 + t342;
t330 = t376 * pkin(10) + t339;
t329 = pkin(5) * t538 + t335;
t328 = t506 * t329;
t1 = [(-t463 * t497 + t474 * t488 - t480 * t518 - t495 * t516) * MDP(17) + (t463 * t518 - t474 * t516) * MDP(11) + (-MDP(20) * t391 - MDP(21) * t392) * t498 + (-t518 * MDP(13) - t520 * MDP(14) + (qJD(1) * (-t474 ^ 2 + t534 ^ 2) - t462 * t534 + t463 * t474) * MDP(12) - t530 * MDP(17) + (qJD(3) * t535 + 0.2e1 * t474 * t480 + t549) * MDP(16)) * t500 + (-t339 * t446 - t342 * t544 + t387 * t391 - t398 * t538 - t610) * MDP(27) + (-t391 * t447 + t421 * t446 + t439 * t544 + t450 * t538 + t610) * MDP(24) + (-t334 * t545 + t335 * t399 + t339 * t398 + t342 * t387 + t343 * t371 + t344 * t367) * MDP(28) + (-t462 * t497 + t534 * t488) * MDP(16) + t632 * t634 + 0.2e1 * t507 * MDP(4) * t568 + (t348 * t446 - t391 * t417 + t538 * t604 + t543 * t628) * MDP(31) + (t332 * t417 + t341 * t391 + t373 * t348 + (-(qJD(6) * t372 + t331) * t628 - t608 - (qJD(6) * t350 + t330) * t446 + t354 * qJD(6) * t445) * t506 + (-(-qJD(6) * t366 + t333) * t628 - (-qJD(6) * t361 + t329) * t446 + t531) * t503) * MDP(35) + (t328 * t446 + t332 * t415 - t340 * t391 + t373 * t349 + (-t330 * t446 - t331 * t628 - t608) * t503 + (t333 * t628 - t531) * t506 + ((-t366 * t506 - t372 * t503) * t628 - t341 * t446 + t354 * t604) * qJD(6)) * MDP(34) + (-t376 * t446 - t391 * t436 - t392 * t544 - t445 * t538) * MDP(19) + (t334 * t445 + t335 * t446 - t343 * t436 + t344 * t544 - t367 * t391 + t371 * t392 - t376 * t545 + t399 * t538) * MDP(25) + (t376 * t450 + t392 * t447 + t421 * t445 - t436 * t439 - t609) * MDP(23) + (-t339 * t445 + t342 * t436 - t376 * t398 - t387 * t392 + t609) * MDP(26) + (-t391 * t544 + t358) * MDP(18) + (-t391 * t628 + t358) * MDP(33) + (t445 * t369 - t349 * t446 + t391 * t415 + (t392 * t506 - t445 * t583) * t628) * MDP(32) + (-pkin(7) * t601 + t505 * t567) * MDP(9) - MDP(7) * t602 + (pkin(7) * t602 + t507 * t567) * MDP(10) + (t348 * t604 + t417 * t543) * MDP(29) + ((-t415 * t503 + t417 * t506) * t392 + (t347 - t349 * t503 + (-t415 * t506 - t417 * t503) * qJD(6)) * t445) * MDP(30) + MDP(6) * t601; t515 - t505 * MDP(4) * t600 + t510 * t632 + (t390 * t544 - t498 * t593 + t533) * MDP(27) + (t591 * t500 + (t463 * t585 - t500 * t570) * pkin(2) + t523) * MDP(17) + (t462 * t496 - t550 * t500 + (-t500 * t576 - t536) * qJD(3) + t537) * MDP(16) + (-t334 * t456 + t335 * t457 - t367 * t594 + t371 * t593 - t387 * t390) * MDP(28) + (t453 * t369 + t456 * t349 + t592 * t415 + (t503 * t560 - t506 * t597) * t628 + t546) * MDP(34) + (t456 * t348 + t560 * t638 + t592 * t417 + (-t453 * t538 + t597 * t628 - t565) * t503 + t566) * MDP(35) + (-t376 * t456 - t436 * t593 + t457 * t538 - t544 * t594 + t547) * MDP(25) + (t436 * t448 + t540 + t630) * MDP(23) + (-t390 * t436 + t627 - t630) * MDP(26) + (-t448 * t544 + t498 * t637 + t539) * MDP(24) + (MDP(9) * t505 * t510 + MDP(10) * t600) * pkin(1); (t537 + (-qJD(3) + t500) * t536) * MDP(16) + (t386 * t498 + (t463 * t544 - t498 * t571) * pkin(3) + t539) * MDP(24) + t515 + (-t376 * t492 - t436 * t595 + t493 * t538 + t544 * t552 + t547) * MDP(25) + (t393 * t544 - t498 * t595 + t533) * MDP(27) + (t491 * t369 + t492 * t349 + t589 * t415 + (t503 * t559 + t506 * t553) * t628 + t546) * MDP(34) + (t385 * t498 + (-t436 * t463 - t498 * t584) * pkin(3) + t540) * MDP(23) + (-t393 * t436 + t498 * t552 + t627) * MDP(26) + (-t334 * t492 + t335 * t493 + t367 * t552 + t371 * t595 - t387 * t393) * MDP(28) + (t492 * t348 + t559 * t638 + t589 * t417 + (-t491 * t538 - t553 * t628 - t565) * t503 + t566) * MDP(35) + (t500 * t548 + t523) * MDP(17); t624 * MDP(19) + (t540 + t606) * MDP(23) + (-t377 * t498 - t557) * MDP(24) + (-pkin(4) * t538 - qJ(5) * t376 + (-t371 - t378) * t544) * MDP(25) + (t627 - t606) * MDP(26) + (t401 * t544 + t498 * t580 - t334) * MDP(27) + (-pkin(4) * t335 - qJ(5) * t334 - t367 * t378 - t371 * t580 - t387 * t401) * MDP(28) + t644 * MDP(31) + t643 * MDP(32) + (qJ(5) * t349 + (-t599 + t631) * t506 + (-t357 * t506 + t503 * t558) * t628 + t581 * t415 + t598) * MDP(34) + (qJ(5) * t348 + t326 + t558 * t638 + t581 * t417 + (t357 * t628 - t565 + t599) * t503) * MDP(35) - (t447 * MDP(24) + (t367 - t580) * MDP(25) + t401 * MDP(26) - t387 * MDP(27) + t417 * MDP(31) - t415 * MDP(32) + t340 * MDP(34) - t341 * MDP(35) + MDP(19) * t436 - t629) * t436 + t642; t528 * MDP(25) + t544 * t436 * MDP(26) + (-t498 ^ 2 - t624) * MDP(27) + (t371 * t498 + t627) * MDP(28) + (-t415 * t498 + t369) * MDP(34) + (-t417 * t498 - t607) * MDP(35) + (-MDP(34) * t639 - MDP(35) * t638) * t628; t417 * t415 * MDP(29) + (-t415 ^ 2 + t417 ^ 2) * MDP(30) + (t596 + t640) * MDP(31) + (t417 * t628 + t370) * MDP(32) + t538 * MDP(33) + (-t330 * t503 + t341 * t628 - t354 * t417 + t328) * MDP(34) + (-t329 * t503 - t330 * t506 + t340 * t628 + t354 * t415) * MDP(35) + (-MDP(31) * t603 - MDP(32) * t417 - MDP(34) * t341 - MDP(35) * t340) * qJD(6);];
tauc  = t1;
