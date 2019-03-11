% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:22
% EndTime: 2019-03-09 05:07:32
% DurationCPUTime: 5.96s
% Computational Cost: add. (3026->473), mult. (6870->623), div. (0->0), fcn. (4334->8), ass. (0->212)
t452 = sin(qJ(4));
t455 = cos(qJ(4));
t524 = t455 * qJD(3);
t453 = sin(qJ(3));
t539 = qJD(1) * t453;
t405 = t452 * t539 - t524;
t537 = qJD(3) * t452;
t407 = t455 * t539 + t537;
t451 = sin(qJ(6));
t454 = cos(qJ(6));
t349 = t405 * t451 + t407 * t454;
t474 = -t454 * t405 + t407 * t451;
t525 = t453 * MDP(27);
t497 = qJD(3) * t525;
t601 = t474 * MDP(23) * t349 + (t349 ^ 2 - t474 ^ 2) * MDP(24) - qJD(1) * t497;
t456 = cos(qJ(3));
t538 = qJD(1) * t456;
t599 = qJD(4) - t538;
t521 = -qJD(6) + t599;
t600 = t349 * t521;
t598 = qJD(6) - qJD(4);
t437 = sin(pkin(10)) * pkin(1) + pkin(7);
t422 = t437 * qJD(1);
t589 = t456 * qJD(2) - t453 * t422;
t488 = qJD(3) * pkin(3) + t589;
t464 = qJ(5) * t407 + t488;
t582 = pkin(4) + pkin(5);
t322 = -t405 * t582 + t464;
t426 = t599 * qJD(5);
t520 = qJD(1) * qJD(3);
t501 = t453 * t520;
t432 = qJ(5) * t501;
t443 = t453 * qJD(2);
t384 = t456 * t422 + t443;
t374 = qJD(3) * pkin(8) + t384;
t375 = t589 * qJD(3);
t438 = -cos(pkin(10)) * pkin(1) - pkin(2);
t398 = -pkin(3) * t456 - pkin(8) * t453 + t438;
t377 = t398 * qJD(1);
t487 = pkin(3) * t453 - pkin(8) * t456;
t416 = t487 * qJD(3);
t397 = qJD(1) * t416;
t532 = qJD(4) * t455;
t509 = t455 * t375 + t377 * t532 + t452 * t397;
t533 = qJD(4) * t452;
t465 = t374 * t533 - t509;
t307 = t426 + t432 - t465;
t503 = t453 * t532;
t535 = qJD(3) * t456;
t461 = t452 * t535 + t503;
t519 = qJD(3) * qJD(4);
t366 = qJD(1) * t461 + t452 * t519;
t304 = pkin(9) * t366 + t307;
t500 = t456 * t520;
t504 = t453 * t533;
t365 = qJD(1) * t504 + (-t500 - t519) * t455;
t490 = t374 * t532 + t452 * t375 + t377 * t533 - t455 * t397;
t305 = pkin(9) * t365 - t501 * t582 + t490;
t495 = t451 * t304 - t454 * t305;
t596 = t322 * t349 + t495;
t594 = MDP(5) * t453;
t447 = t453 ^ 2;
t593 = MDP(6) * (-t456 ^ 2 + t447);
t592 = t521 * t474;
t491 = pkin(4) * t501;
t308 = t490 - t491;
t334 = t455 * t374 + t452 * t377;
t429 = t599 * qJ(5);
t328 = t429 + t334;
t591 = -t328 * t599 + t308;
t590 = qJD(5) * t452 + t384;
t505 = t456 * t524;
t588 = -t504 + t505;
t333 = -t452 * t374 + t455 * t377;
t522 = qJD(5) - t333;
t518 = MDP(17) + MDP(19);
t523 = pkin(9) * t407 - t522;
t316 = -t582 * t599 - t523;
t528 = qJD(6) * t454;
t511 = t454 * t304 + t451 * t305 + t316 * t528;
t586 = -t322 * t474 + t511;
t578 = qJ(5) * t452;
t585 = -t455 * t582 - t578;
t583 = t407 ^ 2;
t581 = pkin(8) - pkin(9);
t580 = pkin(9) * t453;
t579 = qJ(5) * t405;
t577 = qJ(5) * t455;
t376 = qJD(3) * t443 + t422 * t535;
t462 = -qJ(5) * t365 + qJD(5) * t407 - t376;
t314 = pkin(4) * t366 - t462;
t576 = t314 * t452;
t575 = t314 * t455;
t409 = t451 * t452 + t454 * t455;
t466 = t409 * t456;
t559 = t452 * t454;
t560 = t451 * t455;
t473 = -t559 + t560;
t326 = -t453 * t473 * t598 + qJD(3) * t466;
t574 = t326 * t521;
t332 = pkin(4) * t405 - t464;
t572 = t332 * t407;
t571 = t488 * t452;
t570 = t488 * t455;
t569 = t376 * t452;
t568 = t376 * t455;
t567 = t405 * t407;
t566 = t405 * t599;
t565 = t407 * t599;
t564 = t599 * t455;
t563 = t437 * t599;
t562 = t437 * t452;
t324 = pkin(9) * t405 + t334;
t318 = t324 + t429;
t561 = t451 * t318;
t558 = t452 * t456;
t557 = t453 * t455;
t458 = qJD(3) ^ 2;
t556 = t453 * t458;
t555 = t455 * t456;
t554 = t456 * t458;
t482 = pkin(4) * t452 - t577;
t553 = -t599 * t482 + t590;
t552 = -t366 * t557 - t405 * t505;
t551 = -qJD(1) * t466 - t409 * t598;
t507 = t452 * t538;
t529 = qJD(6) * t451;
t550 = t451 * t532 + t452 * t528 - t455 * t529 - t538 * t560 + (t507 - t533) * t454;
t468 = -t452 * t582 + t577;
t549 = t599 * t468 + t590;
t413 = t487 * qJD(1);
t548 = t452 * t413 + t455 * t589;
t547 = t398 * t532 + t452 * t416;
t502 = t447 * t520;
t427 = t455 * t502;
t546 = t505 * t599 + t427;
t412 = t437 * t555;
t545 = t452 * t398 + t412;
t543 = MDP(20) * t452;
t387 = t451 * t557 - t453 * t559;
t541 = qJD(1) * t387;
t388 = t409 * t453;
t540 = qJD(1) * t388;
t423 = qJD(1) * t438;
t536 = qJD(3) * t453;
t534 = qJD(4) * t405;
t530 = qJD(5) * t455;
t527 = t488 * qJD(4);
t526 = t453 * MDP(16);
t517 = MDP(18) - MDP(21);
t516 = pkin(8) * t599 * t452;
t515 = pkin(8) * t564;
t514 = pkin(9) * t555;
t513 = pkin(8) * t536;
t512 = pkin(8) * t524;
t425 = t581 * t455;
t510 = -t454 * t365 + t451 * t366 + t405 * t528;
t338 = qJ(5) * t539 + t548;
t508 = qJ(5) * t536 + t547;
t498 = qJD(3) * t526;
t496 = -pkin(4) - t562;
t494 = -t365 * t451 - t454 * t366;
t493 = t413 * t455 - t452 * t589;
t411 = t437 * t558;
t492 = t398 * t455 - t411;
t343 = -qJ(5) * t456 + t545;
t486 = -qJD(4) * t412 - t398 * t533 + t416 * t455;
t485 = (-t453 * t582 - t514) * qJD(1) - t493 + t598 * t425;
t424 = t581 * t452;
t484 = pkin(9) * t507 - qJD(6) * t424 + t581 * t533 + t338;
t483 = pkin(4) * t455 + t578;
t481 = qJ(5) * t454 - t451 * t582;
t480 = qJ(5) * t451 + t454 * t582;
t479 = t307 * t455 + t308 * t452;
t300 = t451 * t316 + t454 * t318;
t327 = -pkin(4) * t599 + t522;
t478 = t327 * t455 - t328 * t452;
t477 = t327 * t452 + t328 * t455;
t446 = t456 * pkin(4);
t335 = pkin(5) * t456 + t411 + t446 + (-t398 - t580) * t455;
t341 = t452 * t580 + t343;
t476 = t335 * t454 - t341 * t451;
t475 = t335 * t451 + t341 * t454;
t471 = 0.2e1 * qJD(3) * t423;
t469 = t437 + t482;
t467 = t334 * t599 - t490;
t310 = -t407 * t529 + t510;
t463 = -t437 + t468;
t460 = t333 * t599 + t465;
t311 = qJD(6) * t349 + t494;
t417 = -pkin(3) - t483;
t399 = pkin(3) - t585;
t390 = t407 * t536;
t367 = t469 * t453;
t354 = pkin(4) * t407 + t579;
t350 = t463 * t453;
t344 = t446 - t492;
t340 = -pkin(4) * t539 - t493;
t337 = -t407 * t582 - t579;
t336 = -t365 + t566;
t330 = (qJD(4) * t483 - t530) * t453 + t469 * t535;
t325 = qJD(6) * t388 + t451 * t588 - t454 * t461;
t321 = t496 * t536 - t486;
t320 = (qJD(4) * t585 + t530) * t453 + t463 * t535;
t319 = t325 * t521;
t317 = -qJD(5) * t456 + (-t453 * t524 - t533 * t456) * t437 + t508;
t313 = (pkin(9) * qJD(4) - qJD(3) * t437) * t557 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t437) * t452) * t456 + t508;
t312 = pkin(9) * t504 + (-t514 + (-pkin(5) + t496) * t453) * qJD(3) - t486;
t309 = t310 * t456;
t306 = -t366 * t582 + t462;
t299 = t316 * t454 - t561;
t1 = [-0.2e1 * t520 * t593 + (-t437 * t554 + t453 * t471) * MDP(10) + (t437 * t556 + t456 * t471) * MDP(11) + (-t365 * t557 + t407 * t588) * MDP(12) + (-t407 * t503 + (-t407 * t535 + (t365 + t534) * t453) * t452 + t552) * MDP(13) + (t365 * t456 - t504 * t599 + t390 + t546) * MDP(14) + (-t599 * t503 + t366 * t456 + (-t405 * t453 + (-qJD(1) * t447 - t456 * t599) * t452) * qJD(3)) * MDP(15) + (t599 - t538) * t498 + (t486 * t599 + ((t405 * t437 - t571) * qJD(3) + t490) * t456 + (-t455 * t527 + t437 * t366 + t569 + (qJD(1) * t492 + t562 * t599 + t333) * qJD(3)) * t453) * MDP(17) + (-t547 * t599 + ((-t374 + t563) * t533 + (t407 * t437 - t570) * qJD(3) + t509) * t456 + (t452 * t527 - t437 * t365 + t568 + (-qJD(1) * t545 + t455 * t563 - t334) * qJD(3)) * t453) * MDP(18) + (-t321 * t599 + t330 * t405 + t366 * t367 + (t332 * t537 + t308) * t456 + (t332 * t532 + t576 + (-qJD(1) * t344 - t327) * qJD(3)) * t453) * MDP(19) + (-t317 * t405 + t321 * t407 - t343 * t366 - t344 * t365 + t478 * t535 + (-qJD(4) * t477 - t307 * t452 + t308 * t455) * t453) * MDP(20) + (t317 * t599 - t330 * t407 + t365 * t367 + (-t332 * t524 - t307) * t456 + (t332 * t533 - t575 + (qJD(1) * t343 + t328) * qJD(3)) * t453) * MDP(21) + (t307 * t343 + t308 * t344 + t314 * t367 + t317 * t328 + t321 * t327 + t330 * t332) * MDP(22) + (t310 * t388 + t326 * t349) * MDP(23) + (-t310 * t387 - t311 * t388 - t325 * t349 - t326 * t474) * MDP(24) + (-t574 + t309 + (-t349 - t540) * t536) * MDP(25) + (-t311 * t456 + t319 + (t474 + t541) * t536) * MDP(26) + (t521 - t538) * t497 + (-(t312 * t454 - t313 * t451) * t521 - t495 * t456 + t320 * t474 + t350 * t311 + t306 * t387 + t322 * t325 + (-t300 * t456 + t475 * t521) * qJD(6) + (-qJD(1) * t476 - t299) * t536) * MDP(28) + ((qJD(6) * t476 + t312 * t451 + t313 * t454) * t521 - (-t318 * t529 + t511) * t456 + t320 * t349 + t350 * t310 + t306 * t388 + t322 * t326 + (qJD(1) * t475 + t300) * t536) * MDP(29) + MDP(7) * t554 - MDP(8) * t556 + 0.2e1 * t500 * t594; (t390 - t427) * MDP(18) + t552 * MDP(20) + t546 * MDP(21) + t319 * MDP(28) + (t574 + t309) * MDP(29) + (-t458 * MDP(11) - t314 * MDP(22) + t311 * MDP(28) - t518 * t366 + t517 * t365 + (t407 * t543 + t477 * MDP(22) - (t455 * MDP(18) + t452 * t518) * t599) * qJD(3)) * t456 + (-t458 * MDP(10) - t365 * t543 + t479 * MDP(22) + (-MDP(21) * t407 + t332 * MDP(22) + (-t474 + t541) * MDP(28) + (-t349 + t540) * MDP(29)) * qJD(3) + ((t405 * t452 + t407 * t455) * MDP(20) + t478 * MDP(22) - (-t452 * t517 + t455 * t518) * t599) * qJD(4)) * t453 + t518 * (t405 * t536 - t452 * t502); (qJD(3) * t384 - t423 * t539 - t376) * MDP(10) - t423 * t538 * MDP(11) + (-t365 * t452 + t407 * t564) * MDP(12) + ((-t365 - t566) * t455 + (-t366 - t565) * t452) * MDP(13) + t599 * t532 * MDP(14) - t599 * t533 * MDP(15) + (-pkin(3) * t366 - t568 - t493 * t599 - t384 * t405 + (-t515 - t571) * qJD(4)) * MDP(17) + (pkin(3) * t365 + t569 + t548 * t599 - t384 * t407 + (t516 - t570) * qJD(4)) * MDP(18) + (-t575 + t340 * t599 + t366 * t417 - t553 * t405 + (t332 * t452 - t515) * qJD(4)) * MDP(19) + (t338 * t405 - t340 * t407 + (t307 + t599 * t327 + (qJD(4) * t407 - t366) * pkin(8)) * t455 + ((-t365 + t534) * pkin(8) + t591) * t452) * MDP(20) + (-t576 - t338 * t599 + t365 * t417 + t553 * t407 + (-t332 * t455 - t516) * qJD(4)) * MDP(21) + (t314 * t417 - t327 * t340 - t328 * t338 - t553 * t332 + (qJD(4) * t478 + t479) * pkin(8)) * MDP(22) + (-t310 * t473 + t349 * t551) * MDP(23) + (-t310 * t409 + t311 * t473 - t349 * t550 - t474 * t551) * MDP(24) + (-t551 * t521 + (qJD(3) * t473 + t349) * t539) * MDP(25) + (t550 * t521 + (qJD(3) * t409 - t474) * t539) * MDP(26) + (t306 * t409 + t399 * t311 - (t451 * t484 - t454 * t485) * t521 + t549 * t474 + t550 * t322 + (-(t424 * t454 - t425 * t451) * qJD(3) + t299) * t539) * MDP(28) + (-t306 * t473 + t399 * t310 - (t451 * t485 + t454 * t484) * t521 + t549 * t349 + t551 * t322 + ((t424 * t451 + t425 * t454) * qJD(3) - t300) * t539) * MDP(29) + ((-t599 * t555 + (-t407 + t537) * t453) * MDP(14) + (t599 * t558 + (t405 + t524) * t453) * MDP(15) - t599 * t526 + (-t333 * t453 + (t456 * t488 - t513) * t452) * MDP(17) + (t488 * t555 + (t334 - t512) * t453) * MDP(18) + (t327 * t453 + (-t332 * t456 - t513) * t452) * MDP(19) + (t332 * t555 + (-t328 + t512) * t453) * MDP(21) - t521 * t525 + (-t456 * t594 + t593) * qJD(1)) * qJD(1); MDP(12) * t567 + (-t405 ^ 2 + t583) * MDP(13) + t336 * MDP(14) + (-t366 + t565) * MDP(15) + qJD(1) * t498 + (t407 * t488 + t467) * MDP(17) + (-t405 * t488 + t460) * MDP(18) + (-t354 * t405 + t467 + 0.2e1 * t491 - t572) * MDP(19) + (pkin(4) * t365 - qJ(5) * t366 + (t328 - t334) * t407 + (t327 - t522) * t405) * MDP(20) + (-t332 * t405 + t354 * t407 + 0.2e1 * t426 + 0.2e1 * t432 - t460) * MDP(21) + (-pkin(4) * t308 + qJ(5) * t307 - t327 * t334 + t328 * t522 - t332 * t354) * MDP(22) + (-t310 + t592) * MDP(25) + (t311 + t600) * MDP(26) + (t480 * t501 - t337 * t474 - (-t454 * t324 + t451 * t523) * t521 + (t481 * t521 + t300) * qJD(6) + t596) * MDP(28) + (t481 * t501 - t337 * t349 - (t451 * t324 + t454 * t523) * t521 + (-t480 * t521 - t561) * qJD(6) + t586) * MDP(29) - t601; (-t501 + t567) * MDP(19) + t336 * MDP(20) + (-t599 ^ 2 - t583) * MDP(21) + (t572 + t591) * MDP(22) + (-t407 * t474 - t454 * t501) * MDP(28) + (-t349 * t407 + t451 * t501) * MDP(29) - (MDP(28) * t451 + MDP(29) * t454) * t521 ^ 2; (t510 - t592) * MDP(25) + (-t494 - t600) * MDP(26) + (-t300 * t521 - t596) * MDP(28) + (-t299 * t521 - t586) * MDP(29) + ((-MDP(26) * t407 - MDP(28) * t318) * t454 + (-t407 * MDP(25) - t405 * MDP(26) - MDP(28) * t316 + MDP(29) * t318) * t451) * qJD(6) + t601;];
tauc  = t1;
