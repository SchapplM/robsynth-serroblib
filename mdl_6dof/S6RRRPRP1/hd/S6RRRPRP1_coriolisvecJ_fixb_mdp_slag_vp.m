% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP1
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:10
% EndTime: 2019-03-09 16:33:20
% DurationCPUTime: 4.83s
% Computational Cost: add. (8158->423), mult. (20941->552), div. (0->0), fcn. (15482->8), ass. (0->213)
t496 = cos(qJ(2));
t598 = pkin(7) + pkin(8);
t467 = t598 * t496;
t462 = qJD(1) * t467;
t495 = cos(qJ(3));
t450 = t495 * t462;
t493 = sin(qJ(2));
t466 = t598 * t493;
t460 = qJD(1) * t466;
t492 = sin(qJ(3));
t555 = qJD(1) * t496;
t542 = t495 * t555;
t556 = qJD(1) * t493;
t544 = t492 * t556;
t443 = -t542 + t544;
t591 = qJ(4) * t443;
t406 = t460 * t492 - t450 + t591;
t445 = -t492 * t555 - t495 * t556;
t440 = t445 * qJ(4);
t446 = t492 * t462;
t558 = -t495 * t460 - t446;
t407 = t440 + t558;
t489 = sin(pkin(10));
t490 = cos(pkin(10));
t577 = t489 * t492;
t593 = pkin(2) * qJD(3);
t559 = -t406 * t489 - t407 * t490 + (t490 * t495 - t577) * t593;
t549 = qJD(1) * qJD(2);
t541 = t496 * t549;
t548 = qJD(2) + qJD(3);
t419 = qJD(3) * t542 + t495 * t541 - t544 * t548;
t459 = t492 * t496 + t493 * t495;
t601 = qJD(1) * t459;
t501 = t548 * t601;
t499 = t490 * t419 - t489 * t501;
t612 = qJD(5) * t548 + t499;
t491 = sin(qJ(5));
t552 = qJD(5) * t491;
t537 = -t490 * t443 + t445 * t489;
t582 = t537 * t491;
t611 = t552 - t582;
t384 = t419 * t489 + t490 * t501;
t494 = cos(qJ(5));
t600 = qJD(5) - t537;
t531 = t600 * t494;
t610 = -t384 * t491 - t600 * t531;
t609 = -0.2e1 * t549;
t608 = MDP(5) * (t493 ^ 2 - t496 ^ 2);
t606 = t493 * MDP(4);
t592 = qJD(2) * pkin(2);
t452 = -t460 + t592;
t545 = qJD(2) * t598;
t528 = qJD(1) * t545;
t453 = t493 * t528;
t530 = qJD(3) * t452 - t453;
t605 = t530 * t495;
t519 = -t443 * t489 - t490 * t445;
t404 = t491 * t548 + t494 * t519;
t583 = t404 * t491;
t604 = t600 * t583;
t575 = t490 * t492;
t560 = t490 * t406 - t489 * t407 + (t489 * t495 + t575) * t593;
t536 = t495 * t452 - t446;
t398 = t440 + t536;
t603 = qJ(6) * t582 + t494 * qJD(6);
t597 = pkin(3) * t445;
t373 = pkin(4) * t519 - pkin(9) * t537 - t597;
t483 = pkin(2) * t556;
t370 = t373 + t483;
t602 = t491 * t370 - t559 * t494;
t599 = t404 ^ 2;
t596 = pkin(3) * t489;
t595 = pkin(3) * t490;
t594 = pkin(5) * t494;
t390 = pkin(3) * t548 + t398;
t518 = -t452 * t492 - t450;
t399 = -t518 - t591;
t576 = t490 * t399;
t355 = t489 * t390 + t576;
t353 = pkin(9) * t548 + t355;
t482 = -pkin(2) * t496 - pkin(1);
t465 = t482 * qJD(1);
t425 = pkin(3) * t443 + qJD(4) + t465;
t362 = -pkin(4) * t537 - pkin(9) * t519 + t425;
t332 = -t353 * t491 + t494 * t362;
t322 = -qJ(6) * t404 + t332;
t317 = pkin(5) * t600 + t322;
t590 = t317 * t494;
t346 = -t612 * t494 + t519 * t552;
t589 = t346 * t491;
t393 = t489 * t399;
t354 = t490 * t390 - t393;
t352 = -pkin(4) * t548 - t354;
t588 = t352 * t537;
t402 = t491 * t519 - t494 * t548;
t586 = t402 * t519;
t585 = t402 * t537;
t584 = t404 * t519;
t458 = t492 * t493 - t495 * t496;
t422 = -t458 * t489 + t459 * t490;
t581 = t422 * t491;
t580 = t422 * t494;
t578 = t465 * t445;
t497 = qJD(2) ^ 2;
t574 = t493 * t497;
t486 = t494 * qJ(6);
t412 = -qJ(4) * t459 - t466 * t495 - t467 * t492;
t517 = t466 * t492 - t467 * t495;
t413 = -qJ(4) * t458 - t517;
t380 = t412 * t489 + t413 * t490;
t377 = t494 * t380;
t381 = t494 * t384;
t573 = t496 * t497;
t498 = qJD(1) ^ 2;
t572 = t496 * t498;
t481 = pkin(2) * t495 + pkin(3);
t439 = pkin(2) * t575 + t489 * t481;
t433 = pkin(9) + t439;
t571 = -qJ(6) - t433;
t478 = pkin(9) + t596;
t570 = -qJ(6) - t478;
t569 = t317 - t322;
t551 = qJD(5) * t494;
t347 = t612 * t491 + t519 * t551;
t568 = -t491 * t347 - t402 * t551;
t360 = t398 * t490 - t393;
t567 = t494 * t360 + t491 * t373;
t421 = t490 * t458 + t459 * t489;
t515 = pkin(3) * t458 + t482;
t378 = pkin(4) * t421 - pkin(9) * t422 + t515;
t565 = t491 * t378 + t377;
t533 = qJD(5) * t571;
t564 = t491 * t533 - t602 + t603;
t365 = t494 * t370;
t526 = pkin(5) * t519 - t486 * t537;
t563 = t494 * t533 - t365 - t526 + (-qJD(6) - t559) * t491;
t532 = qJD(5) * t570;
t562 = t491 * t532 - t567 + t603;
t538 = t360 * t491 - t494 * t373;
t561 = -qJD(6) * t491 + t494 * t532 - t526 + t538;
t554 = qJD(3) * t492;
t553 = qJD(3) * t495;
t484 = t493 * t592;
t507 = t459 * qJD(3);
t424 = qJD(2) * t459 + t507;
t461 = t493 * t545;
t463 = t496 * t545;
t509 = -t495 * t461 - t492 * t463 - t466 * t553 - t467 * t554;
t366 = -qJ(4) * t424 - qJD(4) * t458 + t509;
t423 = t548 * t458;
t505 = qJD(3) * t517 + t492 * t461 - t495 * t463;
t367 = qJ(4) * t423 - qJD(4) * t459 + t505;
t337 = t366 * t490 + t367 * t489;
t385 = -t423 * t489 + t490 * t424;
t386 = -t423 * t490 - t424 * t489;
t540 = pkin(3) * t424 + t484;
t343 = pkin(4) * t385 - pkin(9) * t386 + t540;
t546 = t494 * t337 + t491 * t343 + t378 * t551;
t479 = -pkin(4) - t595;
t543 = t422 * t551;
t351 = t352 * t551;
t539 = pkin(1) * t609;
t454 = t496 * t528;
t534 = -t492 * t454 - t462 * t554;
t349 = -qJ(4) * t501 - t443 * qJD(4) + t534 + t605;
t435 = t495 * t454;
t535 = t492 * t453 - t435;
t502 = -qJ(4) * t419 + qJD(3) * t518 + qJD(4) * t445 + t535;
t327 = t349 * t489 - t490 * t502;
t336 = t366 * t489 - t490 * t367;
t359 = t398 * t489 + t576;
t379 = -t490 * t412 + t413 * t489;
t438 = -pkin(2) * t577 + t481 * t490;
t432 = -pkin(4) - t438;
t527 = t611 * pkin(5);
t333 = t353 * t494 + t362 * t491;
t525 = t327 * t491 + t333 * t519 + t351;
t323 = -qJ(6) * t402 + t333;
t524 = -t323 * t491 - t590;
t523 = t327 * t422 - t380 * t384;
t522 = -t384 * t433 - t588;
t521 = -t384 * t478 - t588;
t520 = t354 * t537 + t355 * t519;
t516 = -qJ(6) * t386 - qJD(6) * t422;
t514 = -t611 * t600 + t381;
t314 = pkin(5) * t347 + t327;
t513 = -t327 * t494 - t332 * t519 + t352 * t552;
t512 = t465 * t443 - t534;
t511 = t386 * t491 + t543;
t510 = t386 * t494 - t422 * t552;
t328 = t490 * t349 + t489 * t502;
t500 = pkin(3) * t501 + qJD(2) * t483;
t340 = t384 * pkin(4) - pkin(9) * t499 + t500;
t508 = t494 * t328 + t491 * t340 - t353 * t552 + t362 * t551;
t339 = t494 * t340;
t506 = -qJD(5) * t333 - t328 * t491 + t339;
t308 = pkin(5) * t384 + qJ(6) * t346 - qJD(6) * t404 + t506;
t310 = -qJ(6) * t347 - qJD(6) * t402 + t508;
t504 = qJD(5) * t524 - t308 * t491 + t310 * t494 + t323 * t582 + t537 * t590;
t503 = -t445 * t443 * MDP(11) - t600 * t519 * MDP(24) + ((-t346 + t585) * t494 - t604 + t568) * MDP(21) + (t514 + t586) * MDP(23) + (-t584 - t610) * MDP(22) + (t404 * t531 - t589) * MDP(20) + (t443 * t548 + t419) * MDP(13) + (-t445 * t548 - t501) * MDP(14) + (-t443 ^ 2 + t445 ^ 2) * MDP(12);
t456 = t478 * t494 + t486;
t455 = t570 * t491;
t427 = t433 * t494 + t486;
t426 = t571 * t491;
t401 = t402 ^ 2;
t376 = t494 * t378;
t344 = t402 * pkin(5) + qJD(6) + t352;
t342 = t494 * t343;
t334 = -qJ(6) * t581 + t565;
t331 = pkin(5) * t421 - t380 * t491 - t422 * t486 + t376;
t312 = -qJ(6) * t543 + (-qJD(5) * t380 + t516) * t491 + t546;
t311 = pkin(5) * t385 - t337 * t491 + t342 + t516 * t494 + (-t377 + (qJ(6) * t422 - t378) * t491) * qJD(5);
t1 = [0.2e1 * t541 * t606 + t608 * t609 + MDP(6) * t573 - MDP(7) * t574 + (-pkin(7) * t573 + t493 * t539) * MDP(9) + (pkin(7) * t574 + t496 * t539) * MDP(10) + (t419 * t459 + t423 * t445) * MDP(11) + (-t419 * t458 + t423 * t443 + t445 * t424 - t459 * t501) * MDP(12) + (t443 * t484 + t465 * t424 + (t482 * t507 + (t493 * pkin(2) * t458 + t459 * t482) * qJD(2)) * qJD(1)) * MDP(16) + (t482 * t419 - t465 * t423 + (-t445 + t601) * t484) * MDP(17) + (-t328 * t421 + t336 * t519 + t337 * t537 - t354 * t386 - t355 * t385 + t379 * t499 + t523) * MDP(18) + (t327 * t379 + t328 * t380 - t354 * t336 + t355 * t337 + t425 * t540 + t500 * t515) * MDP(19) + (-t346 * t580 + t404 * t510) * MDP(20) + ((-t402 * t494 - t583) * t386 + (t589 - t347 * t494 + (t402 * t491 - t404 * t494) * qJD(5)) * t422) * MDP(21) + (-t346 * t421 + t381 * t422 + t385 * t404 + t510 * t600) * MDP(22) + (-t347 * t421 - t384 * t581 - t385 * t402 - t511 * t600) * MDP(23) + (t384 * t421 + t385 * t600) * MDP(24) + ((-t380 * t551 + t342) * t600 + t376 * t384 + (-t353 * t551 + t339) * t421 + t332 * t385 + t336 * t402 + t379 * t347 + t422 * t351 + ((-qJD(5) * t378 - t337) * t600 + (-qJD(5) * t362 - t328) * t421 + t352 * t386 + t523) * t491) * MDP(25) + (-(-t380 * t552 + t546) * t600 - t565 * t384 - t508 * t421 - t333 * t385 + t336 * t404 - t379 * t346 + t327 * t580 + t510 * t352) * MDP(26) + (-t311 * t404 - t312 * t402 + t331 * t346 - t334 * t347 + t524 * t386 + (-t308 * t494 - t310 * t491 + (t317 * t491 - t323 * t494) * qJD(5)) * t422) * MDP(27) + (t310 * t334 + t323 * t312 + t308 * t331 + t317 * t311 + t314 * (pkin(5) * t581 + t379) + t344 * (pkin(5) * t511 + t336)) * MDP(28) + (-t423 * MDP(13) - t424 * MDP(14) + MDP(16) * t505 - MDP(17) * t509) * t548; (t432 * t347 + t522 * t491 + t560 * t402 + (-t433 * t551 - t491 * t559 - t365) * t600 + t513) * MDP(25) + (-t432 * t346 + t522 * t494 + t560 * t404 + (t433 * t552 + t602) * t600 + t525) * MDP(26) + (t310 * t427 + t308 * t426 + t314 * (t432 - t594) + (t527 + t560) * t344 + t564 * t323 + t563 * t317) * MDP(28) + (t346 * t426 - t347 * t427 - t402 * t564 - t404 * t563 + t504) * MDP(27) + (t445 * t483 + t558 * t548 + (-t548 * t593 - t530) * t495 + t512) * MDP(17) + (-t443 * t483 - t462 * t553 - t530 * t492 - t435 + t578 + (t450 + (-t460 - t593) * t492) * t548) * MDP(16) + (t328 * t439 - t327 * t438 - t425 * (t483 - t597) + t559 * t355 - t560 * t354) * MDP(19) + t503 + (-t439 * t384 - t438 * t499 + t519 * t560 + t537 * t559 + t520) * MDP(18) + t498 * t608 - t572 * t606 + (MDP(9) * t493 * t498 + MDP(10) * t572) * pkin(1); (t346 * t455 - t347 * t456 - t402 * t562 - t404 * t561 + t504) * MDP(27) + (t354 * t359 - t355 * t360 + (-t327 * t490 + t328 * t489 + t425 * t445) * pkin(3)) * MDP(19) + (t479 * t347 - t359 * t402 + t521 * t491 + (-t478 * t551 + t538) * t600 + t513) * MDP(25) + (-t479 * t346 - t359 * t404 + t521 * t494 + (t478 * t552 + t567) * t600 + t525) * MDP(26) + (t310 * t456 + t308 * t455 + t314 * (t479 - t594) + (-t359 + t527) * t344 + t562 * t323 + t561 * t317) * MDP(28) + (t536 * t548 + t512 - t605) * MDP(17) + (-qJD(2) * t518 + t535 + t578) * MDP(16) + t503 + (-t359 * t519 - t360 * t537 - t384 * t596 - t499 * t595 + t520) * MDP(18); (-t519 ^ 2 - t537 ^ 2) * MDP(18) + (t354 * t519 - t355 * t537 + t500) * MDP(19) + (t514 - t586) * MDP(25) + (-t584 + t610) * MDP(26) + ((t346 + t585) * t494 + t604 + t568) * MDP(27) + (-t344 * t519 + (t323 * t600 + t308) * t494 + (-t317 * t600 + t310) * t491) * MDP(28); t404 * t402 * MDP(20) + (-t401 + t599) * MDP(21) + (t402 * t600 - t346) * MDP(22) + (t404 * t600 - t347) * MDP(23) + t384 * MDP(24) + (t333 * t600 - t352 * t404 + t506) * MDP(25) + (t332 * t600 + t352 * t402 - t508) * MDP(26) + (pkin(5) * t346 - t402 * t569) * MDP(27) + (t569 * t323 + (-t344 * t404 + t308) * pkin(5)) * MDP(28); (-t401 - t599) * MDP(27) + (t317 * t404 + t323 * t402 + t314) * MDP(28);];
tauc  = t1;
