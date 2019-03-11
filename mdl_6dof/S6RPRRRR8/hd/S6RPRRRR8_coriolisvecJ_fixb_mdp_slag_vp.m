% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:26
% EndTime: 2019-03-09 07:21:36
% DurationCPUTime: 5.39s
% Computational Cost: add. (5455->455), mult. (11595->613), div. (0->0), fcn. (8368->8), ass. (0->202)
t507 = sin(qJ(4));
t511 = cos(qJ(3));
t614 = cos(qJ(4));
t558 = t614 * t511;
t508 = sin(qJ(3));
t578 = qJD(1) * t508;
t452 = -qJD(1) * t558 + t507 * t578;
t500 = qJD(3) + qJD(4);
t506 = sin(qJ(5));
t510 = cos(qJ(5));
t428 = -t452 * t506 - t510 * t500;
t528 = -t507 * t511 - t508 * t614;
t453 = t528 * qJD(1);
t627 = qJD(5) - t453;
t635 = t428 * t627;
t530 = t452 * t510 - t500 * t506;
t634 = t530 * t627;
t512 = -pkin(1) - pkin(7);
t478 = qJD(1) * t512 + qJD(2);
t448 = -pkin(8) * t578 + t478 * t508;
t439 = t507 * t448;
t577 = qJD(1) * t511;
t449 = -pkin(8) * t577 + t511 * t478;
t403 = t449 * t614 - t439;
t552 = t614 * qJD(4);
t629 = -pkin(3) * t552 + t403;
t416 = -pkin(4) * t452 - pkin(9) * t453;
t404 = pkin(3) * t577 + t416;
t633 = -t510 * t404 + t629 * t506;
t509 = cos(qJ(6));
t505 = sin(qJ(6));
t598 = t530 * t505;
t378 = t509 * t428 - t598;
t447 = qJD(6) + t627;
t632 = t378 * t447;
t531 = t428 * t505 + t509 * t530;
t631 = t447 * t531;
t593 = t505 * t510;
t463 = t506 * t509 + t593;
t616 = qJD(5) + qJD(6);
t421 = t616 * t463;
t586 = t463 * t453 - t421;
t460 = t505 * t506 - t509 * t510;
t630 = (t453 - t616) * t460;
t574 = qJD(4) * t507;
t576 = qJD(3) * t508;
t628 = -t507 * t576 - t508 * t574;
t440 = t614 * t448;
t441 = qJD(3) * pkin(3) + t449;
t395 = t507 * t441 + t440;
t387 = pkin(9) * t500 + t395;
t472 = pkin(3) * t578 + qJD(1) * qJ(2);
t397 = -pkin(4) * t453 + pkin(9) * t452 + t472;
t353 = t387 * t510 + t397 * t506;
t344 = -pkin(10) * t428 + t353;
t571 = qJD(6) * t505;
t342 = t344 * t571;
t394 = t441 * t614 - t439;
t386 = -t500 * pkin(4) - t394;
t363 = t428 * pkin(5) + t386;
t626 = t363 * t378 + t342;
t414 = t453 * t500;
t572 = qJD(5) * t510;
t573 = qJD(5) * t506;
t371 = t510 * t414 + t452 * t573 + t500 * t572;
t519 = t614 * qJD(3) + t552;
t518 = t519 * t511;
t582 = t628 * qJD(1);
t415 = qJD(1) * t518 + t582;
t549 = pkin(8) * qJD(1) - t478;
t442 = t549 * t576;
t575 = qJD(3) * t511;
t443 = t549 * t575;
t357 = t441 * t552 + t507 * t442 - t443 * t614 - t448 * t574;
t501 = qJD(1) * qJD(2);
t567 = qJD(1) * qJD(3);
t551 = t511 * t567;
t468 = pkin(3) * t551 + t501;
t362 = pkin(4) * t415 - pkin(9) * t414 + t468;
t360 = t510 * t362;
t517 = -qJD(5) * t353 - t357 * t506 + t360;
t322 = pkin(5) * t415 - pkin(10) * t371 + t517;
t372 = -qJD(5) * t530 + t414 * t506;
t523 = t510 * t357 + t506 * t362 - t387 * t573 + t397 * t572;
t323 = -pkin(10) * t372 + t523;
t548 = t509 * t322 - t505 * t323;
t625 = t363 * t531 + t548;
t624 = t415 * MDP(32) + (-t378 ^ 2 + t531 ^ 2) * MDP(29) - t378 * MDP(28) * t531;
t623 = MDP(8) * (t508 ^ 2 - t511 ^ 2);
t402 = t507 * t449 + t440;
t540 = pkin(3) * t574 - t402;
t622 = t506 * t404 + t510 * t629;
t621 = t505 * t573 + t506 * t571;
t612 = pkin(8) - t512;
t469 = t612 * t508;
t470 = t612 * t511;
t620 = t507 * t469 - t614 * t470;
t597 = t453 * t506;
t619 = (t573 - t597) * pkin(5);
t618 = -qJD(6) * t510 - t572;
t617 = MDP(12) * qJ(2);
t547 = t371 * t505 + t509 * t372;
t332 = -qJD(6) * t531 + t547;
t615 = -pkin(9) - pkin(10);
t613 = t510 * pkin(5);
t492 = pkin(3) * t507 + pkin(9);
t611 = -pkin(10) - t492;
t514 = qJD(1) ^ 2;
t610 = qJ(2) * t514;
t352 = -t387 * t506 + t510 * t397;
t343 = pkin(10) * t530 + t352;
t339 = pkin(5) * t627 + t343;
t609 = t339 * t509;
t608 = t344 * t509;
t607 = t371 * t506;
t605 = t415 * t460;
t604 = t415 * t463;
t422 = -t507 * t575 - t508 * t519 - t511 * t574;
t603 = t422 * t500;
t602 = t422 * t506;
t601 = t422 * t510;
t423 = t518 + t628;
t600 = t423 * t447;
t599 = t423 * t500;
t596 = t453 * t510;
t461 = t507 * t508 - t558;
t595 = t461 * t506;
t594 = t461 * t510;
t592 = t506 * t415;
t591 = t510 * t415;
t427 = -t469 * t614 - t507 * t470;
t418 = t510 * t427;
t590 = t511 * t514;
t513 = qJD(3) ^ 2;
t589 = t512 * t513;
t489 = t508 * pkin(3) + qJ(2);
t588 = t510 * t394 + t506 * t416;
t417 = -pkin(4) * t528 + pkin(9) * t461 + t489;
t584 = t506 * t417 + t418;
t583 = t619 + t540;
t570 = qJD(6) * t509;
t479 = pkin(3) * t575 + qJD(2);
t566 = 0.2e1 * qJD(1);
t565 = pkin(10) * t597;
t561 = t509 * t371 - t505 * t372 - t428 * t570;
t560 = qJD(5) * t615;
t385 = t386 * t572;
t550 = qJD(5) * t611;
t546 = -t394 * t506 + t510 * t416;
t545 = t510 * t627;
t544 = qJD(6) * t339 + t323;
t543 = -qJD(5) * t528 + qJD(1);
t358 = t441 * t574 - t614 * t442 - t507 * t443 + t448 * t552;
t493 = -pkin(3) * t614 - pkin(4);
t541 = -t452 * pkin(5) - pkin(10) * t596;
t539 = -t395 + t619;
t538 = -t353 * t452 + t358 * t506 + t385;
t499 = t510 * pkin(10);
t456 = t492 * t510 + t499;
t537 = qJD(6) * t456 - t510 * t550 + t541 - t633;
t474 = pkin(9) * t510 + t499;
t536 = qJD(6) * t474 - t510 * t560 + t541 + t546;
t455 = t611 * t506;
t535 = -qJD(6) * t455 - t506 * t550 - t565 + t622;
t473 = t615 * t506;
t534 = -qJD(6) * t473 - t506 * t560 - t565 + t588;
t327 = t339 * t505 + t608;
t532 = -t386 * t453 - t415 * t492;
t529 = t352 * t452 - t358 * t510 + t386 * t573;
t527 = t461 * t572 - t602;
t526 = t461 * t573 + t601;
t525 = t452 * t472 - t358;
t524 = t447 * t460;
t369 = pkin(4) * t423 - pkin(9) * t422 + t479;
t458 = t612 * t576;
t459 = qJD(3) * t470;
t374 = qJD(4) * t620 + t507 * t458 - t614 * t459;
t522 = t506 * t369 + t510 * t374 + t417 * t572 - t427 * t573;
t331 = t530 * t571 + t561;
t326 = -t344 * t505 + t609;
t335 = pkin(5) * t372 + t358;
t521 = t326 * t452 + t335 * t460 - t586 * t363;
t520 = -t327 * t452 + t335 * t463 + t630 * t363;
t516 = -t472 * t453 - t357;
t375 = qJD(4) * t427 - t458 * t614 - t507 * t459;
t515 = (-t331 * t460 - t332 * t463 - t378 * t630 - t531 * t586) * MDP(29) + (t331 * t463 - t531 * t630) * MDP(28) + ((t371 - t635) * t510 + (-t372 + t634) * t506) * MDP(22) + (t447 * t630 - t452 * t531 + t604) * MDP(30) + (-t378 * t452 + t447 * t586 - t605) * MDP(31) + (-t530 * t545 + t607) * MDP(21) + (-t506 * t627 ^ 2 - t428 * t452 + t591) * MDP(24) + (-t452 * t530 + t545 * t627 + t592) * MDP(23) + (-t452 * t500 - t519 * t577 - t582) * MDP(17) + (t452 ^ 2 - t453 ^ 2) * MDP(15) + (MDP(14) * t453 + MDP(25) * t627 + MDP(32) * t447) * t452;
t494 = -pkin(4) - t613;
t471 = t493 - t613;
t412 = t510 * t417;
t410 = t460 * t461;
t409 = t463 * t461;
t396 = -pkin(5) * t595 - t620;
t390 = t415 * t528;
t365 = t510 * t369;
t361 = pkin(10) * t595 + t584;
t354 = -pkin(5) * t528 + pkin(10) * t594 - t427 * t506 + t412;
t347 = -pkin(5) * t527 + t375;
t337 = t422 * t593 + (-t594 * t616 + t602) * t509 + t621 * t461;
t336 = t421 * t461 - t422 * t460;
t328 = pkin(10) * t527 + t522;
t324 = -pkin(10) * t601 + pkin(5) * t423 - t374 * t506 + t365 + (-t418 + (-pkin(10) * t461 - t417) * t506) * qJD(5);
t1 = [0.2e1 * (MDP(6) * qJ(2) + MDP(5)) * t501 + ((qJD(2) * t566 - t589) * MDP(12) - 0.2e1 * MDP(7) * t551 - t513 * MDP(9)) * t508 - MDP(17) * t599 + (-t511 * t589 + (-qJ(2) * t576 + qJD(2) * t511) * t566) * MDP(13) + ((-t427 * t572 + t365) * t627 + t412 * t415 - (-t387 * t572 + t360) * t528 + t352 * t423 + t375 * t428 - t620 * t372 - t461 * t385 + ((-qJD(5) * t417 - t374) * t627 - t427 * t415 - (-qJD(5) * t397 - t357) * t528 - t358 * t461 + t386 * t422) * t506) * MDP(26) + (-t353 * t423 - t358 * t594 - t371 * t620 - t375 * t530 + t386 * t526 - t415 * t584 - t522 * t627 + t523 * t528) * MDP(27) + (-t371 * t528 - t423 * t530 - t461 * t591 + t526 * t627) * MDP(23) + (t423 * t627 - t390) * MDP(25) + (t372 * t528 - t423 * t428 + t461 * t592 + t527 * t627) * MDP(24) + 0.2e1 * t567 * t623 + t575 * t566 * t617 + (-t390 + t600) * MDP(32) + (t331 * t410 - t336 * t531) * MDP(28) + (t331 * t409 - t332 * t410 - t336 * t378 + t337 * t531) * MDP(29) + ((-t428 * t510 + t506 * t530) * t422 + (t607 + t372 * t510 + (-t428 * t506 - t510 * t530) * qJD(5)) * t461) * MDP(22) + (-t371 * t594 - t526 * t530) * MDP(21) + (t332 * t528 - t337 * t447 - t378 * t423 + t409 * t415) * MDP(31) + ((t324 * t509 - t328 * t505) * t447 + (t354 * t509 - t361 * t505) * t415 - t548 * t528 + t326 * t423 + t347 * t378 + t396 * t332 - t335 * t409 + t363 * t337 + ((-t354 * t505 - t361 * t509) * t447 + t327 * t528) * qJD(6)) * MDP(33) + (-t375 * t500 + t415 * t489 + t423 * t472 - t453 * t479 - t468 * t528) * MDP(19) + (t414 * t528 + t415 * t461 + t422 * t453 + t423 * t452) * MDP(15) + (-t327 * t423 + t396 * t331 + t335 * t410 + t363 * t336 - t342 * t528 - t347 * t531 + (-(-qJD(6) * t361 + t324) * t447 - t354 * t415 + t322 * t528) * t505 + (-(qJD(6) * t354 + t328) * t447 - t361 * t415 + t544 * t528) * t509) * MDP(34) + (-t331 * t528 + t336 * t447 + t410 * t415 - t423 * t531) * MDP(30) + (-t374 * t500 + t414 * t489 + t422 * t472 - t452 * t479 - t461 * t468) * MDP(20) + (-t414 * t461 - t422 * t452) * MDP(14) - t513 * t511 * MDP(10) + MDP(16) * t603; -t514 * MDP(5) - MDP(6) * t610 + (qJD(1) * t453 + t603) * MDP(19) + (qJD(1) * t452 - t599) * MDP(20) + (t528 * t592 + t372 * t461 - t422 * t428 + (-t423 * t506 - t510 * t543) * t627) * MDP(26) + (t528 * t591 + t371 * t461 + t422 * t530 + (-t423 * t510 + t506 * t543) * t627) * MDP(27) + (t461 * t332 - t422 * t378 - t463 * t600 + qJD(1) * t524 - ((t509 * t618 + t621) * t447 - t604) * t528) * MDP(33) + (t461 * t331 + t422 * t531 + t423 * t524 + t463 * t447 * qJD(1) - (-(t505 * t618 - t506 * t570 - t509 * t573) * t447 + t605) * t528) * MDP(34) + (t508 * MDP(12) + t511 * MDP(13)) * (-t513 - t514); -t514 * t623 + (-(t455 * t505 + t456 * t509) * t415 + t471 * t331 + (t505 * t537 + t509 * t535) * t447 - t583 * t531 + t520) * MDP(34) + ((t455 * t509 - t456 * t505) * t415 + t471 * t332 + (t505 * t535 - t509 * t537) * t447 + t583 * t378 + t521) * MDP(33) + t515 - t590 * t617 + (t402 * t500 + (t453 * t577 - t500 * t574) * pkin(3) + t525) * MDP(19) + (t403 * t500 + (t452 * t577 - t500 * t552) * pkin(3) + t516) * MDP(20) + (t493 * t372 + t532 * t506 + t540 * t428 + (-t492 * t572 + t633) * t627 + t529) * MDP(26) + (t493 * t371 + t532 * t510 - t540 * t530 + (t492 * t573 + t622) * t627 + t538) * MDP(27) + (MDP(13) * t610 + MDP(7) * t590) * t508; (-(t473 * t505 + t474 * t509) * t415 + t494 * t331 + (t505 * t536 + t509 * t534) * t447 - t539 * t531 + t520) * MDP(34) + t515 + (t395 * t500 + t525) * MDP(19) + ((t473 * t509 - t474 * t505) * t415 + t494 * t332 + (t505 * t534 - t509 * t536) * t447 + t539 * t378 + t521) * MDP(33) + (t394 * t500 + t516) * MDP(20) + (-pkin(4) * t372 - t546 * t627 - t395 * t428 - t386 * t597 + (-t572 * t627 - t592) * pkin(9) + t529) * MDP(26) + (-pkin(4) * t371 + t588 * t627 + t395 * t530 - t386 * t596 + (t573 * t627 - t591) * pkin(9) + t538) * MDP(27); -t530 * t428 * MDP(21) + (-t428 ^ 2 + t530 ^ 2) * MDP(22) + (t371 + t635) * MDP(23) + (-t372 - t634) * MDP(24) + t415 * MDP(25) + (t353 * t627 + t386 * t530 + t517) * MDP(26) + (t352 * t627 + t386 * t428 - t523) * MDP(27) + (t331 + t632) * MDP(30) + (-t332 - t631) * MDP(31) + (-(-t343 * t505 - t608) * t447 - t327 * qJD(6) + (t378 * t530 + t509 * t415 - t447 * t571) * pkin(5) + t625) * MDP(33) + ((-t344 * t447 - t322) * t505 + (t343 * t447 - t544) * t509 + (-t505 * t415 - t447 * t570 - t530 * t531) * pkin(5) + t626) * MDP(34) + t624; (t561 + t632) * MDP(30) + (-t547 - t631) * MDP(31) + (t327 * t447 + t625) * MDP(33) + (-t505 * t322 - t509 * t323 + t326 * t447 + t626) * MDP(34) + (MDP(30) * t598 + MDP(31) * t531 - MDP(33) * t327 - MDP(34) * t609) * qJD(6) + t624;];
tauc  = t1;
