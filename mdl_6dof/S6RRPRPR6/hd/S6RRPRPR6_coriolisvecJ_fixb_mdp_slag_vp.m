% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:43:14
% EndTime: 2019-03-09 10:43:25
% DurationCPUTime: 8.07s
% Computational Cost: add. (7147->539), mult. (21556->700), div. (0->0), fcn. (17167->10), ass. (0->232)
t500 = sin(pkin(11));
t504 = sin(qJ(2));
t507 = cos(qJ(2));
t620 = cos(pkin(11));
t531 = t500 * t507 + t504 * t620;
t501 = sin(pkin(6));
t595 = qJD(1) * t501;
t472 = t531 * t595;
t506 = cos(qJ(4));
t571 = t620 * t507;
t554 = t501 * t571;
t486 = qJD(1) * t554;
t585 = qJD(1) * qJD(2);
t573 = t501 * t585;
t555 = t504 * t573;
t521 = qJD(2) * t486 - t500 * t555;
t621 = cos(pkin(6));
t570 = t621 * qJD(1);
t538 = t570 + qJD(2);
t528 = t506 * t538;
t503 = sin(qJ(4));
t593 = qJD(4) * t503;
t389 = -qJD(4) * t528 + t472 * t593 - t506 * t521;
t437 = t472 * t503 - t528;
t579 = t504 * t595;
t469 = -t500 * t579 + t486;
t465 = qJD(4) - t469;
t612 = t437 * t465;
t638 = t389 - t612;
t642 = t638 * MDP(20);
t505 = cos(qJ(6));
t502 = sin(qJ(6));
t609 = t465 * t502;
t391 = -t505 * t437 + t609;
t513 = -t506 * t472 - t503 * t538;
t636 = qJD(6) - t513;
t641 = t391 * t636;
t393 = t437 * t502 + t465 * t505;
t640 = t393 * t636;
t639 = t636 * t502;
t608 = t469 * t503;
t424 = t472 * t502 - t505 * t608;
t637 = t505 * t593 + t424;
t580 = pkin(1) * t621;
t557 = t507 * t580;
t491 = qJD(1) * t557;
t622 = pkin(8) + qJ(3);
t575 = t622 * t504;
t556 = t501 * t575;
t518 = pkin(2) * t621 - t556;
t441 = qJD(2) * pkin(2) + qJD(1) * t518 + t491;
t558 = t504 * t580;
t604 = t501 * t507;
t468 = t604 * t622 + t558;
t457 = t468 * qJD(1);
t572 = t620 * t457;
t386 = t500 * t441 + t572;
t376 = pkin(9) * t538 + t386;
t553 = (-pkin(2) * t507 - pkin(1)) * t501;
t534 = qJD(1) * t553;
t480 = qJD(3) + t534;
t403 = -pkin(3) * t469 - pkin(9) * t472 + t480;
t353 = t376 * t503 - t506 * t403;
t587 = qJD(5) + t353;
t635 = MDP(18) - MDP(21);
t634 = t465 ^ 2;
t497 = t501 ^ 2;
t633 = -0.2e1 * t497 * t585;
t594 = qJD(2) * t501;
t578 = t504 * t594;
t632 = pkin(2) * t578;
t631 = MDP(5) * (t504 ^ 2 - t507 ^ 2);
t630 = MDP(6) * t507;
t615 = t391 * t465;
t456 = -qJD(1) * t556 + t491;
t407 = t456 * t500 + t572;
t629 = qJD(5) * t503 + t407 + (-t593 + t608) * pkin(4);
t588 = -pkin(5) * t513 + t587;
t583 = -MDP(19) + MDP(22);
t354 = t506 * t376 + t503 * t403;
t347 = -qJ(5) * t465 - t354;
t624 = pkin(5) * t437;
t338 = -t347 - t624;
t626 = pkin(4) + pkin(10);
t628 = t626 * t389 + (t338 - t354 + t624) * t636;
t390 = -qJD(4) * t513 + t503 * t521;
t627 = t513 ^ 2;
t509 = qJD(1) ^ 2;
t475 = t531 * t501;
t470 = qJD(2) * t475;
t463 = qJD(1) * t470;
t625 = pkin(4) * t463;
t494 = pkin(2) * t500 + pkin(9);
t623 = pkin(5) + t494;
t619 = qJ(5) * t437;
t618 = qJ(5) * t506;
t590 = qJD(6) * t505;
t581 = t502 * t390 + t437 * t590 + t505 * t463;
t591 = qJD(6) * t502;
t349 = -t465 * t591 + t581;
t617 = t349 * t505;
t446 = t500 * t457;
t385 = t441 * t620 - t446;
t375 = -pkin(3) * t538 - t385;
t510 = qJ(5) * t513 + t375;
t352 = t437 * pkin(4) + t510;
t616 = t352 * t513;
t614 = t393 * t469;
t611 = t513 * t437;
t610 = t513 * t465;
t458 = t463 * qJ(5);
t607 = t494 * t463;
t606 = t497 * t509;
t605 = t501 * t504;
t603 = t502 * t389;
t382 = t505 * t389;
t592 = qJD(4) * t506;
t602 = t349 * t503 + t393 * t592;
t601 = qJ(5) * t592 - t469 * t618 + t629;
t455 = t557 + t518;
t417 = t500 * t455 + t620 * t468;
t406 = pkin(9) * t621 + t417;
t474 = t500 * t605 - t554;
t423 = pkin(3) * t474 - pkin(9) * t475 + t553;
t600 = t506 * t406 + t503 * t423;
t408 = t456 * t620 - t446;
t419 = pkin(2) * t579 + pkin(3) * t472 - pkin(9) * t469;
t599 = t506 * t408 + t503 * t419;
t598 = t506 * t463 + t465 * t608;
t355 = -qJ(5) * t472 - t599;
t597 = pkin(5) * t608 - t623 * t593 + t355;
t589 = qJD(6) * t506;
t582 = pkin(1) * t606;
t577 = t494 * t593;
t483 = t623 * t506;
t548 = qJD(2) * t570;
t536 = pkin(1) * t548;
t489 = t507 * t536;
t515 = (-qJD(2) * t575 + qJD(3) * t507) * t501;
t431 = qJD(1) * t515 + t489;
t443 = -t468 * qJD(2) - qJD(3) * t605;
t432 = t443 * qJD(1);
t369 = t431 * t620 + t500 * t432;
t488 = pkin(2) * t555;
t404 = t463 * pkin(3) - pkin(9) * t521 + t488;
t559 = t503 * t369 + t376 * t592 + t403 * t593 - t506 * t404;
t325 = -pkin(5) * t389 - t463 * t626 + t559;
t368 = t431 * t500 - t620 * t432;
t517 = qJ(5) * t389 + qJD(5) * t513 + t368;
t331 = t390 * t626 + t517;
t569 = t505 * t325 - t331 * t502;
t560 = t506 * t369 - t376 * t593 + t403 * t592 + t503 * t404;
t329 = -t465 * qJD(5) - t458 - t560;
t344 = -pkin(4) * t465 + t587;
t568 = -t344 * t469 - t329;
t330 = t559 - t625;
t567 = -t347 * t469 + t330;
t566 = -t505 * t390 + t463 * t502;
t565 = -t503 * t406 + t423 * t506;
t492 = qJD(2) * t557;
t442 = t492 + t515;
t378 = t442 * t500 - t620 * t443;
t564 = t465 * t506;
t561 = t497 * t504 * t507 * MDP(4);
t495 = -t620 * pkin(2) - pkin(3);
t552 = t501 * t509 * t621;
t551 = pkin(1) * t633;
t359 = -qJ(5) * t474 - t600;
t550 = t589 * t639 + t636 * t637;
t425 = t472 * t505 + t502 * t608;
t547 = t502 * t593 - t425;
t401 = t503 * t408;
t524 = -t503 * qJ(5) + t495;
t476 = -t506 * t626 + t524;
t546 = qJD(6) * t476 + t401 + (pkin(5) * t469 - t419) * t506 - t626 * t472 - qJD(4) * t483;
t482 = t623 * t503;
t545 = -qJD(6) * t482 + t629 - t465 * (pkin(10) * t503 - t618);
t416 = t455 * t620 - t500 * t468;
t543 = t325 * t502 + t331 * t505;
t336 = -t465 * t626 + t588;
t341 = t437 * t626 + t510;
t322 = t336 * t505 - t341 * t502;
t323 = t336 * t502 + t341 * t505;
t450 = t475 * t506 + t503 * t621;
t340 = pkin(5) * t450 - t474 * t626 - t565;
t449 = t475 * t503 - t506 * t621;
t405 = -pkin(3) * t621 - t416;
t512 = -t450 * qJ(5) + t405;
t351 = t449 * t626 + t512;
t542 = t340 * t505 - t351 * t502;
t541 = t340 * t502 + t351 * t505;
t539 = t449 * t505 - t474 * t502;
t422 = t449 * t502 + t474 * t505;
t379 = t442 * t620 + t500 * t443;
t471 = (-t500 * t504 + t571) * t594;
t420 = pkin(3) * t470 - pkin(9) * t471 + t632;
t535 = -t503 * t379 - t406 * t592 + t420 * t506 - t423 * t593;
t532 = -t636 * t639 - t382;
t530 = t354 * t465 - t559;
t529 = t506 * t379 - t406 * t593 + t503 * t420 + t423 * t592;
t326 = -pkin(5) * t390 - t329;
t527 = t326 + (t636 * t626 + t619) * t636;
t526 = t375 * t465 - t607;
t525 = -t352 * t465 + t607;
t523 = -t505 * t636 ^ 2 + t603;
t522 = -pkin(8) * t604 - t558;
t519 = -t505 * t589 + t547;
t415 = -t475 * t593 + (qJD(4) * t621 + t471) * t506;
t516 = -qJ(5) * t415 - qJD(5) * t450 + t378;
t332 = -qJ(5) * t470 - qJD(5) * t474 - t529;
t511 = t538 * t522;
t481 = -t506 * pkin(4) + t524;
t414 = qJD(4) * t450 + t471 * t503;
t384 = t389 * t503;
t377 = -pkin(4) * t513 + t619;
t365 = t389 * t450;
t362 = t449 * pkin(4) + t512;
t360 = -pkin(4) * t474 - t565;
t358 = qJD(6) * t539 + t414 * t502 + t470 * t505;
t357 = qJD(6) * t422 - t414 * t505 + t470 * t502;
t356 = -pkin(4) * t472 - t419 * t506 + t401;
t350 = t393 * qJD(6) + t566;
t345 = -pkin(5) * t449 - t359;
t337 = pkin(4) * t414 + t516;
t335 = pkin(4) * t390 + t517;
t334 = t414 * t626 + t516;
t333 = -pkin(4) * t470 - t535;
t328 = -pkin(5) * t414 - t332;
t327 = pkin(5) * t415 - t470 * t626 - t535;
t321 = -qJD(6) * t323 + t569;
t320 = qJD(6) * t322 + t543;
t1 = [(-MDP(7) * t578 + t594 * t630) * (0.2e1 * t570 + qJD(2)) + (t349 * t539 - t350 * t422 - t357 * t393 - t358 * t391) * MDP(25) + (-t415 * t513 - t365) * MDP(13) + (-t354 * t470 + t368 * t450 + t375 * t415 - t378 * t513 - t405 * t389 - t463 * t600 - t465 * t529 - t474 * t560) * MDP(19) + (-t389 * t474 + t415 * t465 + t450 * t463 - t470 * t513) * MDP(15) + (-t329 * t474 - t332 * t465 - t335 * t450 + t337 * t513 - t347 * t470 - t352 * t415 - t359 * t463 + t362 * t389) * MDP(22) + (t329 * t449 + t330 * t450 + t332 * t437 - t333 * t513 + t344 * t415 + t347 * t414 + t359 * t390 - t360 * t389) * MDP(20) + (t389 * t449 - t390 * t450 + t414 * t513 - t415 * t437) * MDP(14) + ((-qJD(6) * t541 + t327 * t505 - t334 * t502) * t636 - t542 * t389 + t321 * t450 + t322 * t415 + t328 * t391 + t345 * t350 - t326 * t539 + t338 * t357) * MDP(29) + (-t350 * t450 - t357 * t636 - t389 * t539 - t391 * t415) * MDP(27) + (t415 * t636 - t365) * MDP(28) + (-(qJD(6) * t542 + t327 * t502 + t334 * t505) * t636 + t541 * t389 - t320 * t450 - t323 * t415 + t328 * t393 + t345 * t349 + t326 * t422 + t338 * t358) * MDP(30) + (t349 * t450 + t358 * t636 - t389 * t422 + t393 * t415) * MDP(26) + (t368 * t475 - t369 * t474 + t378 * t472 + t379 * t469 - t385 * t471 - t386 * t470 - t416 * t521 - t417 * t463) * MDP(11) + (-t353 * t470 + t368 * t449 + t375 * t414 + t378 * t437 + t405 * t390 + t463 * t565 + t465 * t535 - t474 * t559) * MDP(18) + 0.2e1 * t561 * t585 + (-t368 * t416 + t369 * t417 - t385 * t378 + t386 * t379 + (t480 + t534) * t632) * MDP(12) + (t507 * t551 - (-pkin(8) * t578 + t492) * t538 - (-pkin(8) * t555 + t489) * t621) * MDP(10) + (qJD(2) * t511 + t504 * t551 + t522 * t548) * MDP(9) + (t463 * t474 + t465 * t470) * MDP(17) + (-t390 * t474 - t414 * t465 - t437 * t470 - t449 * t463) * MDP(16) + (t330 * t474 + t333 * t465 - t335 * t449 - t337 * t437 + t344 * t470 - t352 * t414 + t360 * t463 - t362 * t390) * MDP(21) + t631 * t633 + (t329 * t359 + t330 * t360 + t332 * t347 + t333 * t344 + t335 * t362 + t337 * t352) * MDP(23) + (t349 * t422 + t358 * t393) * MDP(24); (-t355 * t437 + t356 * t513 + (-t390 * t494 + (-t494 * t513 + t344) * qJD(4) + t568) * t506 + (-t389 * t494 + (t437 * t494 + t347) * qJD(4) + t567) * t503) * MDP(20) + (-t513 * t564 - t384) * MDP(13) + (t354 * t472 + t368 * t503 - t495 * t389 + t407 * t513 + (t577 + t599) * t465 + t526 * t506) * MDP(19) + (-t335 * t503 + t347 * t472 + t389 * t481 + (t355 - t577) * t465 - t601 * t513 + t525 * t506) * MDP(22) + (t463 * t503 + t465 * t564 + t472 * t513) * MDP(15) + (t564 * t636 - t384) * MDP(28) + ((t476 * t505 + t482 * t502) * t389 - t320 * t503 + t483 * t349 + (t502 * t546 + t505 * t545) * t636 + t597 * t393 + t547 * t338 + (-t323 * t465 - t326 * t502 - t338 * t590) * t506) * MDP(30) + (-(-t476 * t502 + t482 * t505) * t389 + t321 * t503 + t483 * t350 + (t502 * t545 - t505 * t546) * t636 + t597 * t391 - t637 * t338 + (t322 * t465 + t326 * t505 - t338 * t591) * t506) * MDP(29) + ((t603 - t614) * t506 + t519 * t636 + t602) * MDP(26) + (MDP(7) * t552 + (-t536 + t582) * MDP(9)) * t504 + (-pkin(8) * t507 * t573 - qJD(1) * t511) * MDP(9) + ((t386 - t407) * t472 + (-t408 + t385) * t469 + (-t500 * t463 - t521 * t620) * pkin(2)) * MDP(11) + (t437 * t472 - t465 * t593 + t598) * MDP(16) + (t335 * t506 - t344 * t472 - t390 * t481 + (t494 * t592 - t356) * t465 + t601 * t437 + t525 * t503) * MDP(21) + (t353 * t472 - t368 * t506 + t495 * t390 - t407 * t437 + (t401 + (-qJD(4) * t494 - t419) * t506) * t465 + t526 * t503) * MDP(18) - t509 * t561 + ((-t389 - t612) * t506 + (-t390 + t610) * t503) * MDP(14) - t465 * t472 * MDP(17) + (-t489 + t507 * t582 + (-pkin(8) * t579 + t491) * t570 + t491 * qJD(2)) * MDP(10) - t552 * t630 + (t385 * t407 - t386 * t408 + (-t368 * t620 + t369 * t500 - t480 * t579) * pkin(2)) * MDP(12) + (t391 * t425 + t393 * t424 + (-t391 * t502 + t393 * t505) * t593 + (-t617 + t350 * t502 + (t391 * t505 + t393 * t502) * qJD(6)) * t506) * MDP(25) + (-t350 * t503 + (t382 - t615) * t506 + t550) * MDP(27) + (t335 * t481 - t344 * t356 - t347 * t355 - t601 * t352 + (-t329 * t506 + t330 * t503 + (t344 * t506 + t347 * t503) * qJD(4)) * t494) * MDP(23) + (-t349 * t502 * t506 + t393 * t519) * MDP(24) + t606 * t631; -t469 ^ 2 * MDP(11) + (-t386 * t469 + t488) * MDP(12) + t550 * MDP(29) + (t425 * t636 + t602) * MDP(30) + (-MDP(11) * t472 + t385 * MDP(12) - t352 * MDP(23) - t437 * t635 - t513 * t583) * t472 + ((t469 * t513 - t390) * MDP(20) + t568 * MDP(23) + t350 * MDP(29) + t583 * t463 + (-MDP(20) * t513 + t344 * MDP(23) - MDP(30) * t639 - t465 * t635) * qJD(4)) * t503 + (t642 + (-qJD(4) * t347 - t567) * MDP(23) + (t382 + t615) * MDP(29) + (t590 * t636 - t603 - t614) * MDP(30) + t583 * t634) * t506 + t635 * t598; -MDP(13) * t611 + (-t437 ^ 2 + t627) * MDP(14) - t638 * MDP(15) + (-t390 - t610) * MDP(16) + t463 * MDP(17) + (t375 * t513 + t530) * MDP(18) + (-t353 * t465 + t375 * t437 - t560) * MDP(19) + (pkin(4) * t389 - qJ(5) * t390 - (-t347 - t354) * t513 + (t344 - t587) * t437) * MDP(20) + (t377 * t437 - t530 - t616 - 0.2e1 * t625) * MDP(21) + (-t352 * t437 - t377 * t513 + t465 * t587 - t329 + t458) * MDP(22) + (-pkin(4) * t330 - qJ(5) * t329 - t344 * t354 - t347 * t587 - t352 * t377) * MDP(23) + (-t393 * t639 + t617) * MDP(24) + ((-t350 - t640) * t505 + (-t349 + t641) * t502) * MDP(25) + (t393 * t437 + t532) * MDP(26) + (-t391 * t437 + t523) * MDP(27) + t636 * t437 * MDP(28) + (qJ(5) * t350 + t322 * t437 + t588 * t391 + t527 * t502 + t628 * t505) * MDP(29) + (qJ(5) * t349 - t323 * t437 + t588 * t393 - t628 * t502 + t527 * t505) * MDP(30); -t642 + (t463 + t611) * MDP(21) + (-t627 - t634) * MDP(22) + (t347 * t465 + t330 - t616) * MDP(23) + (t532 - t615) * MDP(29) + (-t393 * t465 + t523) * MDP(30); t393 * t391 * MDP(24) + (-t391 ^ 2 + t393 ^ 2) * MDP(25) + (t581 + t641) * MDP(26) + (-t566 + t640) * MDP(27) - t389 * MDP(28) + (t323 * t636 - t338 * t393 + t569) * MDP(29) + (t322 * t636 + t338 * t391 - t543) * MDP(30) + (-MDP(26) * t609 - MDP(27) * t393 - MDP(29) * t323 - MDP(30) * t322) * qJD(6);];
tauc  = t1;
