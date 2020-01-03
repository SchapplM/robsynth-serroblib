% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:30:08
% EndTime: 2019-12-31 22:30:25
% DurationCPUTime: 9.23s
% Computational Cost: add. (4693->465), mult. (11683->645), div. (0->0), fcn. (8425->8), ass. (0->205)
t504 = sin(qJ(2));
t508 = cos(qJ(2));
t468 = -pkin(2) * t508 - pkin(7) * t504 - pkin(1);
t448 = t468 * qJD(1);
t578 = qJD(1) * t508;
t495 = pkin(6) * t578;
t474 = qJD(2) * pkin(7) + t495;
t503 = sin(qJ(3));
t507 = cos(qJ(3));
t412 = t507 * t448 - t474 * t503;
t577 = qJD(2) * t503;
t579 = qJD(1) * t504;
t457 = t507 * t579 + t577;
t382 = -pkin(8) * t457 + t412;
t487 = -qJD(3) + t578;
t374 = -pkin(3) * t487 + t382;
t605 = t503 * t448;
t413 = t474 * t507 + t605;
t559 = t503 * t579;
t566 = t507 * qJD(2);
t455 = t559 - t566;
t383 = -pkin(8) * t455 + t413;
t506 = cos(qJ(4));
t379 = t506 * t383;
t502 = sin(qJ(4));
t336 = t374 * t502 + t379;
t404 = t506 * t455 + t457 * t502;
t635 = pkin(9) * t404;
t331 = t336 - t635;
t501 = sin(qJ(5));
t570 = qJD(5) * t501;
t329 = t331 * t570;
t505 = cos(qJ(5));
t526 = t455 * t502 - t506 * t457;
t614 = t526 * t501;
t361 = -t505 * t404 + t614;
t473 = -qJD(2) * pkin(2) + pkin(6) * t579;
t425 = pkin(3) * t455 + t473;
t372 = pkin(4) * t404 + t425;
t642 = -t372 * t361 + t329;
t528 = t404 * t501 + t505 * t526;
t549 = MDP(29) * t579;
t641 = qJD(2) * t549 + (-t361 ^ 2 + t528 ^ 2) * MDP(26) + t361 * t528 * MDP(25);
t565 = qJD(1) * qJD(2);
t552 = t508 * t565;
t574 = qJD(3) * t503;
t556 = t504 * t574;
t564 = qJD(2) * qJD(3);
t422 = -qJD(1) * t556 + (t552 + t564) * t507;
t573 = qJD(3) * t507;
t554 = t504 * t573;
t575 = qJD(2) * t508;
t558 = t503 * t575;
t516 = t554 + t558;
t423 = qJD(1) * t516 + t503 * t564;
t571 = qJD(4) * t506;
t572 = qJD(4) * t502;
t348 = t506 * t422 - t502 * t423 - t455 * t571 - t457 * t572;
t532 = pkin(2) * t504 - pkin(7) * t508;
t466 = t532 * qJD(2);
t449 = qJD(1) * t466;
t553 = t504 * t565;
t535 = pkin(6) * t553;
t587 = -t507 * t449 - t503 * t535;
t515 = -qJD(3) * t413 - t587;
t343 = pkin(3) * t553 - pkin(8) * t422 + t515;
t523 = t448 * t573 + t503 * t449 - t474 * t574;
t511 = -t507 * t535 + t523;
t353 = -pkin(8) * t423 + t511;
t545 = t506 * t343 - t502 * t353;
t514 = -t336 * qJD(4) + t545;
t318 = pkin(4) * t553 - pkin(9) * t348 + t514;
t512 = qJD(4) * t526 - t422 * t502 - t506 * t423;
t534 = -t502 * t343 - t506 * t353 - t374 * t571 + t383 * t572;
t319 = pkin(9) * t512 - t534;
t640 = -t501 * t318 - t505 * t319 + t642;
t569 = qJD(5) * t505;
t561 = t505 * t348 - t404 * t569 + t501 * t512;
t324 = t526 * t570 + t561;
t478 = -qJD(4) + t487;
t544 = t348 * t501 - t505 * t512;
t513 = qJD(5) * t528 - t544;
t550 = MDP(22) * t579;
t471 = -qJD(5) + t478;
t631 = t471 * t528;
t632 = t361 * t471;
t639 = qJD(2) * t550 + (-t404 ^ 2 + t526 ^ 2) * MDP(19) + (-t404 * t478 + t348) * MDP(20) + (t478 * t526 + t512) * MDP(21) - t404 * t526 * MDP(18) + (t513 + t631) * MDP(28) + (t324 + t632) * MDP(27) + t641;
t546 = t505 * t318 - t501 * t319;
t628 = t372 * t528 + t546;
t598 = t507 * t508;
t524 = pkin(3) * t504 - pkin(8) * t598;
t616 = pkin(7) + pkin(8);
t560 = qJD(3) * t616;
t463 = t532 * qJD(1);
t583 = pkin(6) * t559 + t507 * t463;
t636 = qJD(1) * t524 + t507 * t560 + t583;
t442 = t503 * t463;
t602 = t504 * t507;
t603 = t503 * t508;
t629 = -t442 - (-t602 * pkin(6) - pkin(8) * t603) * qJD(1) - t503 * t560;
t634 = pkin(9) * t526;
t458 = t502 * t503 - t506 * t507;
t519 = t458 * t508;
t617 = qJD(3) + qJD(4);
t589 = qJD(1) * t519 - t617 * t458;
t459 = t502 * t507 + t503 * t506;
t588 = (-t578 + t617) * t459;
t627 = t404 * t425 + t534;
t626 = t425 * t526 + t514;
t623 = -0.2e1 * t565;
t622 = MDP(4) * t504;
t499 = t504 ^ 2;
t621 = MDP(5) * (-t508 ^ 2 + t499);
t435 = t459 * t504;
t620 = t636 * t506;
t454 = t507 * t468;
t615 = pkin(6) * t503;
t411 = -pkin(8) * t602 + t454 + (-pkin(3) - t615) * t508;
t489 = pkin(6) * t598;
t581 = t503 * t468 + t489;
t604 = t503 * t504;
t418 = -pkin(8) * t604 + t581;
t590 = t502 * t411 + t506 * t418;
t533 = -t495 + (-t503 * t578 + t574) * pkin(3);
t475 = t616 * t503;
t476 = t616 * t507;
t584 = -t502 * t475 + t506 * t476;
t619 = -t475 * t571 - t476 * t572 - t636 * t502 + t629 * t506;
t618 = t508 * t566 - t556;
t613 = t422 * t503;
t612 = t455 * t487;
t611 = t457 * t487;
t610 = t473 * t503;
t609 = t473 * t507;
t608 = t487 * t507;
t607 = t501 * t502;
t377 = t502 * t383;
t606 = t502 * t505;
t509 = qJD(2) ^ 2;
t601 = t504 * t509;
t335 = t506 * t374 - t377;
t330 = t335 + t634;
t328 = -pkin(4) * t478 + t330;
t600 = t505 * t328;
t599 = t505 * t331;
t597 = t508 * t509;
t510 = qJD(1) ^ 2;
t596 = t508 * t510;
t408 = t505 * t458 + t459 * t501;
t595 = -qJD(5) * t408 - t501 * t588 + t505 * t589;
t409 = -t458 * t501 + t459 * t505;
t594 = qJD(5) * t409 + t501 * t589 + t505 * t588;
t593 = t506 * t382 - t377;
t592 = pkin(4) * t588 + t533;
t586 = t503 * t466 + t468 * t573;
t576 = qJD(2) * t504;
t585 = t507 * t466 + t576 * t615;
t467 = pkin(3) * t604 + t504 * pkin(6);
t562 = pkin(3) * qJD(4) * t471;
t426 = pkin(3) * t516 + pkin(6) * t575;
t493 = -pkin(3) * t507 - pkin(2);
t555 = t508 * t574;
t548 = MDP(15) * t576;
t401 = pkin(3) * t423 + pkin(6) * t552;
t547 = pkin(1) * t623;
t365 = t524 * qJD(2) + (-t489 + (pkin(8) * t504 - t468) * t503) * qJD(3) + t585;
t368 = -t516 * pkin(8) + (-t566 * t504 - t555) * pkin(6) + t586;
t543 = t506 * t365 - t368 * t502;
t542 = -t382 * t502 - t379;
t541 = t506 * t411 - t418 * t502;
t539 = -t506 * t475 - t476 * t502;
t538 = t455 + t566;
t537 = -t457 + t577;
t536 = qJD(5) * t328 + t319;
t389 = -pkin(9) * t458 + t584;
t531 = pkin(4) * t579 + t589 * pkin(9) + qJD(4) * t584 + qJD(5) * t389 + t502 * t629 + t620;
t388 = -pkin(9) * t459 + t539;
t530 = -t588 * pkin(9) + qJD(5) * t388 + t619;
t321 = t501 * t328 + t599;
t436 = t458 * t504;
t352 = -pkin(4) * t508 + pkin(9) * t436 + t541;
t354 = -pkin(9) * t435 + t590;
t529 = t352 * t501 + t354 * t505;
t384 = t505 * t435 - t436 * t501;
t385 = -t435 * t501 - t436 * t505;
t525 = qJD(1) * t499 - t487 * t508;
t492 = pkin(3) * t506 + pkin(4);
t522 = pkin(3) * t606 + t492 * t501;
t521 = -pkin(3) * t607 + t492 * t505;
t518 = t502 * t365 + t506 * t368 + t411 * t571 - t418 * t572;
t430 = pkin(4) * t458 + t493;
t414 = pkin(4) * t435 + t467;
t376 = pkin(3) * t457 - pkin(4) * t526;
t370 = -t572 * t604 + (t602 * t617 + t558) * t506 + t618 * t502;
t369 = -qJD(2) * t519 - t435 * t617;
t355 = pkin(4) * t370 + t426;
t334 = -pkin(4) * t512 + t401;
t333 = t593 + t634;
t332 = t542 + t635;
t327 = qJD(5) * t385 + t369 * t501 + t505 * t370;
t326 = -qJD(5) * t384 + t369 * t505 - t370 * t501;
t323 = -pkin(9) * t370 + t518;
t322 = pkin(4) * t576 - pkin(9) * t369 - qJD(4) * t590 + t543;
t320 = -t331 * t501 + t600;
t1 = [(-t348 * t508 - t369 * t478) * MDP(20) + (-t324 * t508 - t326 * t471) * MDP(27) + ((-t455 * t507 - t457 * t503) * t575 + (-t613 - t423 * t507 + (t455 * t503 - t457 * t507) * qJD(3)) * t504) * MDP(12) + ((-pkin(6) * t555 + t586) * t487 + t523 * t508 + (pkin(6) * t422 - t473 * t574) * t504 + ((pkin(6) * t457 + t609) * t508 + (-pkin(6) * t608 - qJD(1) * t581 - t413) * t504) * qJD(2)) * MDP(17) + (-(-t468 * t574 + t585) * t487 + (t473 * t573 + pkin(6) * t423 + (t454 * qJD(1) + t412) * qJD(2)) * t504 + ((pkin(6) * t455 + t610) * qJD(2) + (t605 + (pkin(6) * t487 + t474) * t507) * qJD(3) + t587) * t508) * MDP(16) - MDP(7) * t601 + (pkin(6) * t601 + t508 * t547) * MDP(10) + (-pkin(6) * t597 + t504 * t547) * MDP(9) + (-t487 - t578) * t548 + (-(t322 * t505 - t323 * t501) * t471 - t546 * t508 - t355 * t361 - t414 * t513 + t334 * t384 + t372 * t327 + (t321 * t508 + t471 * t529) * qJD(5)) * MDP(30) + (-t324 * t384 + t326 * t361 + t327 * t528 + t385 * t513) * MDP(26) + ((-qJD(1) * t590 - t336) * MDP(24) + ((t352 * t505 - t354 * t501) * qJD(1) + t320) * MDP(30) + (-qJD(1) * t436 - t526) * MDP(20) + (-qJD(1) * t435 - t404) * MDP(21) + (qJD(1) * t385 - t528) * MDP(27) + (-qJD(1) * t384 + t361) * MDP(28) + (qJD(1) * t541 + t335) * MDP(23) + (-qJD(1) * t529 - t321) * MDP(31) + (-t478 - t578) * MDP(22) + (-t471 - t578) * MDP(29)) * t576 + (t487 * t554 + t423 * t508 + (-t455 * t504 - t503 * t525) * qJD(2)) * MDP(14) + (t487 * t556 - t422 * t508 + (t457 * t504 + t507 * t525) * qJD(2)) * MDP(13) + 0.2e1 * t552 * t622 + (t422 * t602 + t457 * t618) * MDP(11) + (t327 * t471 - t508 * t513) * MDP(28) + (t370 * t478 - t508 * t512) * MDP(21) + (-t348 * t435 - t369 * t404 + t370 * t526 - t436 * t512) * MDP(19) + (-t543 * t478 - t545 * t508 + t426 * t404 - t467 * t512 + t401 * t435 + t425 * t370 + (t336 * t508 + t478 * t590) * qJD(4)) * MDP(23) + t621 * t623 + (t414 * t324 + t372 * t326 - t329 * t508 + t334 * t385 - t355 * t528 + ((-qJD(5) * t354 + t322) * t471 + t318 * t508) * t501 + ((qJD(5) * t352 + t323) * t471 + t536 * t508) * t505) * MDP(31) + (t324 * t385 - t326 * t528) * MDP(25) + (t467 * t348 + t425 * t369 - t401 * t436 - t426 * t526 + t518 * t478 - t534 * t508) * MDP(24) + (-t348 * t436 - t369 * t526) * MDP(18) + MDP(6) * t597; (-t457 * t608 + t613) * MDP(11) + (-pkin(2) * t422 - t442 * t487 + (-pkin(7) * t487 * t503 + t609) * qJD(3) + (-t473 * t598 + (-pkin(7) * t566 + t413) * t504 + (t487 * t602 + t508 * t537) * pkin(6)) * qJD(1)) * MDP(17) + (-pkin(2) * t423 + t583 * t487 + (pkin(7) * t608 + t610) * qJD(3) + ((-pkin(7) * t577 - t412) * t504 + (-pkin(6) * t538 - t610) * t508) * qJD(1)) * MDP(16) + ((t422 + t612) * t507 + (-t423 + t611) * t503) * MDP(12) + (t487 * t574 + (-t487 * t603 + t504 * t538) * qJD(1)) * MDP(14) + (t324 * t409 - t528 * t595) * MDP(25) + (t430 * t324 + t334 * t409 + (-t501 * t531 + t505 * t530) * t471 + t595 * t372 - t592 * t528 + (-(t388 * t501 + t389 * t505) * qJD(2) + t321) * t579) * MDP(31) + (-t324 * t408 + t361 * t595 + t409 * t513 + t528 * t594) * MDP(26) + (-t595 * t471 + (qJD(2) * t409 + t528) * t579) * MDP(27) + (-t487 * t573 + (t487 * t598 + t504 * t537) * qJD(1)) * MDP(13) + (t588 * t478 + (-qJD(2) * t458 + t404) * t579) * MDP(21) + (-t493 * t512 + t401 * t458 + (t476 * t571 + (-qJD(4) * t475 + t629) * t502 + t620) * t478 + t588 * t425 + t533 * t404 + (qJD(2) * t539 - t335) * t579) * MDP(23) + (t348 * t459 - t526 * t589) * MDP(18) + (-t348 * t458 - t404 * t589 + t459 * t512 + t526 * t588) * MDP(19) + (-t589 * t478 + (qJD(2) * t459 + t526) * t579) * MDP(20) + (t493 * t348 + t401 * t459 + t619 * t478 + t589 * t425 - t533 * t526 + (-qJD(2) * t584 + t336) * t579) * MDP(24) + (-t430 * t513 + t334 * t408 + (t501 * t530 + t505 * t531) * t471 + t594 * t372 - t592 * t361 + ((t388 * t505 - t389 * t501) * qJD(2) - t320) * t579) * MDP(30) + (t594 * t471 + (-qJD(2) * t408 - t361) * t579) * MDP(28) - t596 * t622 + t510 * t621 + t471 * t549 + t478 * t550 + t487 * MDP(15) * t579 + (MDP(9) * t504 * t510 + MDP(10) * t596) * pkin(1); (t422 - t612) * MDP(13) + (t521 * t553 + (t332 * t505 - t333 * t501) * t471 + t376 * t361 - (-t501 * t506 - t606) * t562 + (t471 * t522 - t321) * qJD(5) + t628) * MDP(30) + (-t455 ^ 2 + t457 ^ 2) * MDP(12) + (-t412 * t487 + t455 * t473 - t511) * MDP(17) + (-t423 - t611) * MDP(14) + qJD(1) * t548 + (-t522 * t553 - (t332 * t501 + t333 * t505) * t471 + t376 * t528 + (t505 * t506 - t607) * t562 + (t471 * t521 - t600) * qJD(5) + t640) * MDP(31) + (t542 * t478 + (-t404 * t457 + t478 * t572 + t506 * t553) * pkin(3) + t626) * MDP(23) + (-t593 * t478 + (t457 * t526 + t478 * t571 - t502 * t553) * pkin(3) + t627) * MDP(24) + t457 * t455 * MDP(11) + (-t413 * t487 - t457 * t473 + t515) * MDP(16) + t639; (-t336 * t478 + t626) * MDP(23) + (-t335 * t478 + t627) * MDP(24) + ((-t330 * t501 - t599) * t471 - t321 * qJD(5) + (-t361 * t526 + t471 * t570 + t505 * t553) * pkin(4) + t628) * MDP(30) + ((t331 * t471 - t318) * t501 + (-t330 * t471 - t536) * t505 + (t471 * t569 - t501 * t553 - t526 * t528) * pkin(4) + t642) * MDP(31) + t639; (t561 + t632) * MDP(27) + (-t544 + t631) * MDP(28) + (-t321 * t471 + t628) * MDP(30) + (-t320 * t471 + t640) * MDP(31) + (MDP(27) * t614 + MDP(28) * t528 - MDP(30) * t321 - MDP(31) * t600) * qJD(5) + t641;];
tauc = t1;
