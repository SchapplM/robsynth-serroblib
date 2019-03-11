% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:30
% EndTime: 2019-03-09 16:37:41
% DurationCPUTime: 5.95s
% Computational Cost: add. (10683->479), mult. (27369->606), div. (0->0), fcn. (20273->8), ass. (0->222)
t506 = cos(qJ(5));
t576 = qJD(5) * t506;
t504 = sin(qJ(3));
t505 = sin(qJ(2));
t507 = cos(qJ(2));
t625 = cos(qJ(3));
t534 = -t504 * t505 + t507 * t625;
t461 = t534 * qJD(1);
t476 = t504 * t507 + t625 * t505;
t463 = qJD(1) * t476;
t502 = sin(pkin(10));
t618 = cos(pkin(10));
t559 = t618 * t461 - t463 * t502;
t601 = t559 * t506;
t637 = t576 - t601;
t503 = sin(qJ(5));
t577 = qJD(5) * t503;
t636 = t559 * t503 - t577;
t628 = qJD(5) - t559;
t640 = t628 ^ 2;
t457 = t463 * qJ(4);
t626 = pkin(7) + pkin(8);
t484 = t626 * t507;
t479 = qJD(1) * t484;
t464 = t504 * t479;
t483 = t626 * t505;
t477 = qJD(1) * t483;
t619 = qJD(2) * pkin(2);
t471 = -t477 + t619;
t558 = t625 * t471 - t464;
t418 = -t457 + t558;
t572 = qJD(2) + qJD(3);
t412 = pkin(3) * t572 + t418;
t468 = t625 * t479;
t536 = -t504 * t471 - t468;
t598 = t461 * qJ(4);
t419 = -t536 + t598;
t565 = t618 * t419;
t378 = t502 * t412 + t565;
t376 = pkin(9) * t572 + t378;
t496 = -pkin(2) * t507 - pkin(1);
t482 = qJD(1) * t496;
t445 = -pkin(3) * t461 + qJD(4) + t482;
t530 = t502 * t461 + t463 * t618;
t385 = -pkin(4) * t559 - pkin(9) * t530 + t445;
t350 = t376 * t506 + t385 * t503;
t338 = qJ(6) * t628 + t350;
t639 = t338 * t628;
t423 = t503 * t572 + t506 * t530;
t557 = t423 * t628;
t516 = t463 * t572;
t425 = t477 * t504 - t468 - t598;
t581 = -t625 * t477 - t464;
t426 = -t457 + t581;
t391 = t502 * t425 + t426 * t618;
t596 = t502 * t504;
t620 = pkin(2) * qJD(3);
t454 = (t618 * t625 - t596) * t620;
t582 = t454 - t391;
t638 = t582 * t503;
t556 = t628 * t503;
t443 = t572 * t534;
t564 = t618 * t443;
t568 = qJD(3) * t625;
t578 = qJD(3) * t504;
t629 = qJD(1) * (t564 + t502 * (-qJD(2) * t476 - t505 * t568 - t507 * t578));
t635 = qJD(5) * t572 + t629;
t573 = qJD(1) * qJD(2);
t634 = -0.2e1 * t573;
t633 = MDP(5) * (t505 ^ 2 - t507 ^ 2);
t571 = qJD(2) * t626;
t553 = qJD(1) * t571;
t472 = t505 * t553;
t473 = t507 * t553;
t523 = t471 * t568 - t472 * t625 - t504 * t473 - t479 * t578;
t371 = -qJ(4) * t516 + t461 * qJD(4) + t523;
t515 = t443 * qJD(1);
t452 = t625 * t473;
t540 = -t479 * t568 - t452;
t594 = t504 * t472;
t631 = -qJ(4) * t515 - t463 * qJD(4) - t471 * t578 + t540 + t594;
t344 = t618 * t371 + t502 * t631;
t406 = t502 * t515 + t618 * t516;
t567 = t505 * t573;
t490 = pkin(2) * t567;
t513 = pkin(3) * t516 + t490;
t361 = t406 * pkin(4) - pkin(9) * t629 + t513;
t527 = t506 * t344 + t503 * t361 - t376 * t577 + t385 * t576;
t617 = qJ(6) * t406;
t324 = qJD(6) * t628 + t527 + t617;
t554 = t503 * t344 - t506 * t361 + t376 * t576 + t385 * t577;
t622 = pkin(5) * t406;
t326 = t554 - t622;
t632 = t324 * t506 + t326 * t503;
t441 = t476 * t502 - t534 * t618;
t442 = t476 * t618 + t502 * t534;
t543 = -pkin(3) * t534 + t496;
t399 = pkin(4) * t441 - pkin(9) * t442 + t543;
t432 = -t476 * qJ(4) - t483 * t625 - t504 * t484;
t535 = t504 * t483 - t484 * t625;
t433 = qJ(4) * t534 - t535;
t401 = t502 * t432 + t433 * t618;
t584 = t503 * t399 + t506 * t401;
t563 = t618 * t504;
t583 = t618 * t425 - t426 * t502 + (t502 * t625 + t563) * t620;
t630 = t636 * pkin(5) + qJ(6) * t637 + t503 * qJD(6);
t627 = t423 ^ 2;
t624 = pkin(3) * t463;
t623 = pkin(3) * t502;
t621 = pkin(5) * t530;
t343 = t371 * t502 - t618 * t631;
t368 = -t506 * t635 + t530 * t577;
t369 = t503 * t635 + t530 * t576;
t329 = pkin(5) * t369 + qJ(6) * t368 - qJD(6) * t423 + t343;
t616 = t329 * t503;
t615 = t350 * t628;
t413 = t502 * t419;
t377 = t412 * t618 - t413;
t375 = -pkin(4) * t572 - t377;
t421 = t503 * t530 - t506 * t572;
t353 = t421 * pkin(5) - t423 * qJ(6) + t375;
t614 = t353 * t559;
t613 = t368 * t503;
t612 = t369 * t506;
t611 = t375 * t559;
t495 = pkin(2) * t625 + pkin(3);
t456 = pkin(2) * t563 + t502 * t495;
t450 = pkin(9) + t456;
t610 = t406 * t450;
t492 = pkin(9) + t623;
t609 = t406 * t492;
t608 = t421 * t530;
t607 = t421 * t559;
t606 = t421 * t503;
t605 = t423 * t421;
t604 = t423 * t530;
t603 = t423 * t506;
t600 = t442 * t506;
t599 = t454 * t506;
t597 = t482 * t463;
t383 = t418 * t618 - t413;
t595 = t503 * t383;
t402 = t503 * t406;
t508 = qJD(2) ^ 2;
t592 = t505 * t508;
t403 = t506 * t406;
t591 = t507 * t508;
t509 = qJD(1) ^ 2;
t590 = t507 * t509;
t382 = t418 * t502 + t565;
t589 = t382 + t630;
t588 = t630 - t583;
t587 = -t503 * t369 - t421 * t576;
t395 = pkin(4) * t530 - pkin(9) * t559 + t624;
t586 = t506 * t383 + t503 * t395;
t579 = qJD(1) * t505;
t497 = pkin(2) * t579;
t392 = t395 + t497;
t585 = t506 * t391 + t503 * t392;
t349 = -t376 * t503 + t385 * t506;
t574 = qJD(6) - t349;
t499 = t505 * t619;
t569 = t618 * pkin(3);
t444 = t572 * t476;
t566 = pkin(3) * t444 + t499;
t562 = pkin(1) * t634;
t430 = t530 * qJ(6);
t347 = t430 + t585;
t561 = -t347 + t599;
t560 = t392 * t506 + t621 + t638;
t478 = t505 * t571;
t480 = t507 * t571;
t528 = -t625 * t478 - t504 * t480 - t483 * t568 - t484 * t578;
t388 = -qJ(4) * t444 + qJD(4) * t534 + t528;
t524 = qJD(3) * t535 + t504 * t478 - t625 * t480;
t389 = -t443 * qJ(4) - t476 * qJD(4) + t524;
t357 = t388 * t502 - t618 * t389;
t400 = -t618 * t432 + t433 * t502;
t493 = -t569 - pkin(4);
t551 = t343 * t503 + t350 * t530 + t375 * t576;
t550 = t506 * pkin(5) + t503 * qJ(6);
t549 = pkin(5) * t503 - qJ(6) * t506;
t455 = -pkin(2) * t596 + t495 * t618;
t337 = -pkin(5) * t628 + t574;
t548 = t337 * t506 - t338 * t503;
t547 = t343 * t442 - t401 * t406;
t546 = -t610 - t611;
t545 = -t609 - t611;
t544 = t377 * t559 + t378 * t530;
t542 = t628 * t637 + t402;
t449 = -pkin(4) - t455;
t541 = t628 * t636 + t403;
t539 = -t329 * t506 + t337 * t530 + t353 * t577;
t538 = -t338 * t530 + t353 * t601 - t616;
t537 = -t343 * t506 - t349 * t530 + t375 * t577;
t408 = -t502 * t444 + t564;
t533 = t408 * t503 + t442 * t576;
t532 = -t408 * t506 + t442 * t577;
t531 = -t450 * t577 + t599;
t529 = t353 * t423 + t554;
t358 = t388 * t618 + t502 * t389;
t407 = t443 * t502 + t444 * t618;
t365 = pkin(4) * t407 - pkin(9) * t408 + t566;
t526 = t506 * t358 + t503 * t365 + t399 * t576 - t401 * t577;
t525 = t337 * t637 + t338 * t636 + t632;
t522 = qJD(5) * t548 + t632;
t521 = -t613 - t612 + (t603 + t606) * qJD(5);
t520 = -t463 * t461 * MDP(11) - t628 * t530 * MDP(24) + ((-t368 + t607) * t506 - t423 * t556 + t587) * MDP(21) + (t541 + t608) * MDP(23) + (t542 - t604) * MDP(22) + (t506 * t557 - t613) * MDP(20) + (-t461 * t572 + t515) * MDP(13) + (-t461 ^ 2 + t463 ^ 2) * MDP(12);
t519 = -t482 * t461 - t523;
t470 = -t550 + t493;
t446 = t449 - t550;
t393 = pkin(5) * t423 + qJ(6) * t421;
t366 = t442 * t549 + t400;
t356 = -pkin(5) * t441 - t399 * t506 + t401 * t503;
t355 = qJ(6) * t441 + t584;
t352 = t421 * t628 - t368;
t346 = -t395 * t506 + t595 - t621;
t345 = t430 + t586;
t331 = t549 * t408 + (qJD(5) * t550 - qJD(6) * t506) * t442 + t357;
t328 = -pkin(5) * t407 + qJD(5) * t584 + t358 * t503 - t365 * t506;
t327 = qJ(6) * t407 + qJD(6) * t441 + t526;
t1 = [(0.2e1 * t482 * t444 - t461 * t499 - t490 * t534) * MDP(16) + (-t344 * t441 + t357 * t530 + t358 * t559 - t377 * t408 - t378 * t407 + t400 * t629 + t547) * MDP(18) + t633 * t634 + (t443 * t461 - t463 * t444) * MDP(12) + (-t327 * t421 + t328 * t423 - t355 * t369 - t356 * t368 + t548 * t408 + (-t324 * t503 + t326 * t506 + (-t337 * t503 - t338 * t506) * qJD(5)) * t442) * MDP(28) + (t463 * t443 + t476 * t515) * MDP(11) + (t406 * t441 + t407 * t628) * MDP(24) + MDP(6) * t591 + (t324 * t355 + t326 * t356 + t327 * t338 + t328 * t337 + t329 * t366 + t331 * t353) * MDP(30) + 0.2e1 * t507 * MDP(4) * t567 + (-t554 * t441 + t349 * t407 + t357 * t421 + t400 * t369 + ((-qJD(5) * t401 + t365) * t628 + t399 * t406 + t375 * qJD(5) * t442) * t506 + ((-qJD(5) * t399 - t358) * t628 + t375 * t408 + t547) * t503) * MDP(25) + (t343 * t400 + t344 * t401 - t377 * t357 + t378 * t358 + t445 * t566 + t513 * t543) * MDP(19) + (t482 * t443 + t463 * t499 + t476 * t490 + t496 * t515) * MDP(17) + (-pkin(7) * t591 + t505 * t562) * MDP(9) + (-t368 * t441 + t403 * t442 + t407 * t423 - t532 * t628) * MDP(22) - MDP(7) * t592 + (pkin(7) * t592 + t507 * t562) * MDP(10) + (-t369 * t441 - t402 * t442 - t407 * t421 - t533 * t628) * MDP(23) + (-t368 * t600 - t423 * t532) * MDP(20) + (t324 * t441 + t327 * t628 - t329 * t600 - t331 * t423 + t338 * t407 + t353 * t532 + t355 * t406 + t366 * t368) * MDP(29) + (t343 * t600 - t350 * t407 + t357 * t423 - t400 * t368 - t375 * t532 - t406 * t584 - t441 * t527 - t526 * t628) * MDP(26) + ((-t421 * t506 - t423 * t503) * t408 + (t613 - t612 + (-t603 + t606) * qJD(5)) * t442) * MDP(21) + (-t326 * t441 - t328 * t628 + t331 * t421 - t337 * t407 + t353 * t533 - t356 * t406 + t366 * t369 + t442 * t616) * MDP(27) + (t524 * MDP(16) + qJD(1) * (-t476 ^ 2 + t534 ^ 2) * MDP(12) + t443 * MDP(13) - t444 * MDP(14) - t528 * MDP(17)) * t572; t520 + (t344 * t456 - t343 * t455 - t445 * (t497 + t624) + t582 * t378 - t583 * t377) * MDP(19) + (t581 * t572 + (-t463 * t579 - t568 * t572) * pkin(2) + t519) * MDP(17) + (-t456 * t406 - t455 * t629 + t530 * t583 + t559 * t582 + t544) * MDP(18) + t509 * t633 + (-t449 * t368 + t546 * t506 + t583 * t423 + (-t531 + t585) * t628 + t551) * MDP(26) + (t449 * t369 + t546 * t503 + t583 * t421 + ((-qJD(5) * t450 - t392) * t506 - t638) * t628 + t537) * MDP(25) + (t461 * t497 - t597 + (-qJD(3) * t471 + t472) * t504 + t540 + (t468 + (-t477 - t620) * t504) * t572) * MDP(16) + (-t421 * t561 + t423 * t560 + t450 * t521 + t525) * MDP(28) + (t369 * t446 + (-t610 - t614) * t503 - t588 * t421 + (-t450 * t576 - t560) * t628 + t539) * MDP(27) + (t368 * t446 + (-qJD(5) * t353 + t610) * t506 + t588 * t423 + (-t347 + t531) * t628 + t538) * MDP(29) + (t329 * t446 + t337 * t560 + t338 * t561 - t353 * t588 + t450 * t522) * MDP(30) - t505 * MDP(4) * t590 + (MDP(9) * t505 * t509 + MDP(10) * t590) * pkin(1); t520 + (t493 * t369 - t382 * t421 + t545 * t503 + (t595 + (-qJD(5) * t492 - t395) * t506) * t628 + t537) * MDP(25) + (t558 * t572 + t519) * MDP(17) + (-t382 * t530 - t383 * t559 - t406 * t623 - t569 * t629 + t544) * MDP(18) + (-t493 * t368 - t382 * t423 + t545 * t506 + (t492 * t577 + t586) * t628 + t551) * MDP(26) + (t329 * t470 - t337 * t346 - t338 * t345 - t353 * t589 + t492 * t522) * MDP(30) + (t377 * t382 - t378 * t383 + (-t343 * t618 + t344 * t502 - t445 * t463) * pkin(3)) * MDP(19) + (-qJD(2) * t536 - t452 + t594 - t597) * MDP(16) + (t345 * t421 - t346 * t423 + t492 * t521 + t525) * MDP(28) + (t492 * t403 - t345 * t628 + t368 * t470 + t589 * t423 + (-t353 * t506 - t492 * t556) * qJD(5) + t538) * MDP(29) + (t369 * t470 + (-t609 - t614) * t503 + (-t492 * t576 + t346) * t628 - t589 * t421 + t539) * MDP(27); (-t530 ^ 2 - t559 ^ 2) * MDP(18) + (t377 * t530 - t378 * t559 + t513) * MDP(19) + (t541 - t608) * MDP(25) + (-t506 * t640 - t402 - t604) * MDP(26) + (-t556 * t628 + t403 - t608) * MDP(27) + ((t368 + t607) * t506 + t503 * t557 + t587) * MDP(28) + (t542 + t604) * MDP(29) + (-t353 * t530 + (-t326 + t639) * t506 + (t337 * t628 + t324) * t503) * MDP(30); MDP(20) * t605 + (-t421 ^ 2 + t627) * MDP(21) + t352 * MDP(22) + (-t369 + t557) * MDP(23) + t406 * MDP(24) + (-t375 * t423 - t554 + t615) * MDP(25) + (t349 * t628 + t375 * t421 - t527) * MDP(26) + (-t393 * t421 - t529 + t615 + 0.2e1 * t622) * MDP(27) + (pkin(5) * t368 - qJ(6) * t369 + (t338 - t350) * t423 + (t337 - t574) * t421) * MDP(28) + (0.2e1 * t617 - t353 * t421 + t393 * t423 + (0.2e1 * qJD(6) - t349) * t628 + t527) * MDP(29) + (-pkin(5) * t326 + qJ(6) * t324 - t337 * t350 + t338 * t574 - t353 * t393) * MDP(30); (-t406 + t605) * MDP(27) + t352 * MDP(28) + (-t627 - t640) * MDP(29) + (t529 - t622 - t639) * MDP(30);];
tauc  = t1;
