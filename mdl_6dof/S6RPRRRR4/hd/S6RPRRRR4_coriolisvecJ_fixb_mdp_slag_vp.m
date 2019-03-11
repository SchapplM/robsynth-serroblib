% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:06:04
% EndTime: 2019-03-09 07:06:15
% DurationCPUTime: 7.04s
% Computational Cost: add. (10261->413), mult. (27681->537), div. (0->0), fcn. (23052->10), ass. (0->192)
t458 = cos(qJ(6));
t520 = qJD(6) * t458;
t455 = sin(qJ(5));
t459 = cos(qJ(5));
t452 = sin(pkin(11));
t457 = sin(qJ(3));
t453 = cos(pkin(11));
t461 = cos(qJ(3));
t540 = t453 * t461;
t478 = t452 * t457 - t540;
t428 = t478 * qJD(1);
t541 = t452 * t461;
t434 = t453 * t457 + t541;
t429 = t434 * qJD(1);
t456 = sin(qJ(4));
t460 = cos(qJ(4));
t482 = -t428 * t456 + t460 * t429;
t483 = -t428 * t460 - t456 * t429;
t368 = t455 * t482 - t459 * t483;
t594 = t368 * t458;
t598 = t520 + t594;
t526 = qJD(3) * t461;
t441 = t453 * qJD(1) * t526;
t527 = qJD(1) * t457;
t513 = t452 * t527;
t422 = -qJD(3) * t513 + t441;
t431 = t434 * qJD(3);
t423 = qJD(1) * t431;
t365 = qJD(4) * t482 + t422 * t456 + t460 * t423;
t560 = pkin(7) + qJ(2);
t438 = t560 * t452;
t435 = qJD(1) * t438;
t439 = t560 * t453;
t436 = qJD(1) * t439;
t574 = -t461 * t435 - t436 * t457;
t379 = -pkin(8) * t423 - qJD(2) * t428 + t574 * qJD(3);
t395 = -pkin(8) * t429 + t574;
t394 = qJD(3) * pkin(3) + t395;
t473 = t434 * qJD(2);
t471 = qJD(1) * t473;
t480 = t435 * t457 - t436 * t461;
t380 = -pkin(8) * t422 + qJD(3) * t480 - t471;
t396 = -pkin(8) * t428 - t480;
t525 = qJD(4) * t456;
t504 = t456 * t380 - t396 * t525;
t567 = t460 * (qJD(4) * t394 + t379) + t504;
t315 = -pkin(9) * t365 + t567;
t391 = t456 * t396;
t503 = t460 * t394 - t391;
t579 = pkin(9) * t482;
t346 = t503 - t579;
t451 = qJD(3) + qJD(4);
t344 = pkin(4) * t451 + t346;
t524 = qJD(4) * t460;
t364 = t460 * t422 - t456 * t423 - t428 * t524 - t429 * t525;
t393 = t460 * t396;
t486 = -t394 * t456 - t393;
t505 = -t456 * t379 + t460 * t380;
t467 = qJD(4) * t486 + t505;
t316 = -pkin(9) * t364 + t467;
t578 = pkin(9) * t483;
t347 = -t486 + t578;
t523 = qJD(5) * t455;
t507 = t316 * t455 - t347 * t523;
t300 = (qJD(5) * t344 + t315) * t459 + t507;
t444 = -pkin(2) * t453 - pkin(1);
t437 = qJD(1) * t444 + qJD(2);
t413 = pkin(3) * t428 + t437;
t381 = -pkin(4) * t483 + t413;
t588 = t368 * t381;
t597 = -t300 + t588;
t454 = sin(qJ(6));
t521 = qJD(6) * t454;
t522 = qJD(5) * t459;
t332 = t459 * t364 - t455 * t365 - t482 * t523 + t483 * t522;
t448 = qJD(5) + t451;
t533 = t458 * t332 + t448 * t520;
t566 = t455 * t483 + t459 * t482;
t312 = -t521 * t566 + t533;
t311 = t312 * t458;
t357 = t448 * t454 + t458 * t566;
t557 = t332 * t454;
t313 = t357 * qJD(6) + t557;
t545 = t566 * t454;
t355 = -t458 * t448 + t545;
t596 = -t454 * t313 - t355 * t598 + t311;
t310 = t312 * t454;
t333 = qJD(5) * t566 + t364 * t455 + t459 * t365;
t329 = t454 * t333;
t592 = -qJD(6) - t368;
t534 = -t520 * t592 + t329;
t547 = t368 * t448;
t550 = t566 * t448;
t552 = t357 * t566;
t595 = (-t333 + t550) * MDP(25) - t368 ^ 2 * MDP(23) + (MDP(22) * t368 + MDP(23) * t566 + MDP(33) * t592) * t566 + (t332 + t547) * MDP(24) + (t357 * t598 + t310) * MDP(29) + (-t592 * t594 + t534 - t552) * MDP(31);
t555 = t347 * t455;
t321 = t344 * t459 - t555;
t319 = -pkin(5) * t448 - t321;
t589 = t319 * t368;
t554 = t347 * t459;
t322 = t344 * t455 + t554;
t508 = t315 * t455 - t459 * t316;
t301 = qJD(5) * t322 + t508;
t575 = t566 * t381;
t593 = -t301 - t575;
t343 = pkin(5) * t566 + pkin(10) * t368;
t331 = t458 * t333;
t553 = t355 * t566;
t586 = t592 * t454;
t590 = (-t586 * t592 + t331 + t553) * MDP(32) + (t357 * t586 + t596) * MDP(30) + t595;
t585 = (t452 ^ 2 + t453 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t320 = pkin(10) * t448 + t322;
t336 = pkin(5) * t368 - pkin(10) * t566 + t381;
t307 = t320 * t458 + t336 * t454;
t570 = t301 * t454 + t307 * t566 + t319 * t520;
t488 = t320 * t454 - t336 * t458;
t571 = -t301 * t458 + t319 * t521 + t488 * t566;
t543 = t483 * t451;
t544 = t482 * t451;
t583 = (t482 ^ 2 - t483 ^ 2) * MDP(16) - t483 * MDP(15) * t482 + (t364 - t543) * MDP(17) + (-t365 + t544) * MDP(18);
t561 = pkin(4) * t482;
t577 = t413 * t482;
t576 = t413 * t483;
t573 = -t521 * t592 - t331;
t412 = t434 * t460 - t456 * t478;
t430 = t478 * qJD(3);
t376 = qJD(4) * t412 - t430 * t456 + t460 * t431;
t469 = -t438 * t526 + qJD(2) * t540 + (-qJD(2) * t452 - qJD(3) * t439) * t457;
t385 = -pkin(8) * t431 + t469;
t479 = t438 * t457 - t439 * t461;
t463 = qJD(3) * t479 - t473;
t386 = pkin(8) * t430 + t463;
t402 = -pkin(8) * t434 - t438 * t461 - t439 * t457;
t403 = -pkin(8) * t478 - t479;
t474 = t460 * t385 + t456 * t386 + t402 * t524 - t403 * t525;
t327 = -pkin(9) * t376 + t474;
t481 = -t434 * t456 - t460 * t478;
t375 = qJD(4) * t481 - t430 * t460 - t431 * t456;
t485 = -t402 * t456 - t403 * t460;
t466 = qJD(4) * t485 - t385 * t456 + t460 * t386;
t328 = -pkin(9) * t375 + t466;
t350 = -pkin(9) * t412 + t402 * t460 - t403 * t456;
t351 = pkin(9) * t481 - t485;
t487 = t350 * t459 - t351 * t455;
t302 = qJD(5) * t487 + t327 * t459 + t328 * t455;
t335 = t350 * t455 + t351 * t459;
t377 = t412 * t455 - t459 * t481;
t337 = -qJD(5) * t377 + t375 * t459 - t376 * t455;
t378 = t412 * t459 + t455 * t481;
t417 = pkin(3) * t478 + t444;
t389 = -pkin(4) * t481 + t417;
t341 = pkin(5) * t377 - pkin(10) * t378 + t389;
t565 = t301 * t378 + t319 * t337 - t335 * t333 + (qJD(6) * t341 + t302) * t592 - (qJD(6) * t336 + t300) * t377;
t562 = pkin(3) * t429;
t558 = t319 * t378;
t556 = t341 * t333;
t551 = t357 * t454;
t539 = t455 * t456;
t538 = t456 * t459;
t502 = -t395 * t456 - t393;
t348 = t502 - t578;
t532 = t460 * t395 - t391;
t349 = t532 - t579;
t447 = pkin(3) * t460 + pkin(4);
t535 = t348 * t455 + t349 * t459 - t447 * t522 - (-t456 * t523 + (t459 * t460 - t539) * qJD(4)) * pkin(3);
t531 = t348 * t459 - t349 * t455 + t447 * t523 + (t456 * t522 + (t455 * t460 + t538) * qJD(4)) * pkin(3);
t515 = qJD(1) * qJD(2);
t511 = -pkin(3) * t451 - t394;
t510 = -pkin(4) * t448 - t344;
t387 = t562 + t561;
t425 = pkin(3) * t538 + t447 * t455 + pkin(10);
t494 = qJD(6) * t425 + t343 + t387;
t445 = pkin(4) * t455 + pkin(10);
t493 = qJD(6) * t445 + t343 + t561;
t323 = t346 * t455 + t554;
t492 = pkin(4) * t523 - t323;
t324 = t346 * t459 - t555;
t491 = -pkin(4) * t522 + t324;
t352 = pkin(3) * t423 + pkin(4) * t365;
t358 = pkin(3) * t431 + pkin(4) * t376;
t490 = -t333 * t445 + t589;
t489 = -t333 * t425 + t589;
t476 = t368 * t586 - t573;
t475 = t337 * t458 - t378 * t521;
t446 = -pkin(4) * t459 - pkin(5);
t424 = pkin(3) * t539 - t447 * t459 - pkin(5);
t338 = qJD(5) * t378 + t375 * t455 + t459 * t376;
t308 = pkin(5) * t338 - pkin(10) * t337 + t358;
t305 = pkin(5) * t333 - pkin(10) * t332 + t352;
t304 = t458 * t305;
t303 = qJD(5) * t335 + t327 * t455 - t328 * t459;
t1 = [(t417 * t364 + t413 * t375 + (t412 * t423 + t431 * t482) * pkin(3)) * MDP(21) + (t444 * t422 - t437 * t430) * MDP(14) + (t417 * t365 + t413 * t376 + (-t423 * t481 - t431 * t483) * pkin(3)) * MDP(20) + (t333 * t389 + t338 * t381 + t352 * t377 + t358 * t368) * MDP(27) + (t332 * t389 + t337 * t381 + t352 * t378 + t358 * t566) * MDP(28) + (t444 * t423 + t437 * t431) * MDP(13) + ((-t355 * t458 - t551) * t337 + (-t310 - t313 * t458 + (t355 * t454 - t357 * t458) * qJD(6)) * t378) * MDP(30) + (t303 * t357 - t307 * t338 - t487 * t312 + ((-qJD(6) * t335 + t308) * t592 - t556 - (-qJD(6) * t320 + t305) * t377 - qJD(6) * t558) * t454 + t565 * t458) * MDP(35) + (t303 * t355 + t304 * t377 - t488 * t338 - t487 * t313 + (-t308 * t592 + t556 + (-t320 * t377 + t335 * t592 + t558) * qJD(6)) * t458 + t565 * t454) * MDP(34) + (t311 * t378 + t357 * t475) * MDP(29) + (-t378 * t329 - t313 * t377 - t338 * t355 - (-t337 * t454 - t378 * t520) * t592) * MDP(32) + (t333 * t377 - t338 * t592) * MDP(33) + (t312 * t377 + t331 * t378 + t338 * t357 - t475 * t592) * MDP(31) + (t332 * t378 + t337 * t566) * MDP(22) + (-t332 * t377 - t333 * t378 - t337 * t368 - t338 * t566) * MDP(23) + (t364 * t412 + t375 * t482) * MDP(15) + (t364 * t481 - t365 * t412 + t375 * t483 - t376 * t482) * MDP(16) + (t422 * t434 - t429 * t430) * MDP(8) + (-t422 * t478 - t423 * t434 + t428 * t430 - t429 * t431) * MDP(9) + 0.2e1 * t515 * t585 + (t375 * MDP(17) - t376 * MDP(18) + MDP(20) * t466 - MDP(21) * t474) * t451 + (MDP(24) * t337 - MDP(25) * t338 - MDP(27) * t303 - MDP(28) * t302) * t448 + (-t430 * MDP(10) - t431 * MDP(11) + MDP(13) * t463 - MDP(14) * t469) * qJD(3); t441 * MDP(14) + (t365 + t544) * MDP(20) + (t364 + t543) * MDP(21) + (t333 + t550) * MDP(27) + (t332 - t547) * MDP(28) + (t476 - t553) * MDP(34) + (-t458 * t592 ^ 2 - t329 - t552) * MDP(35) + ((qJD(1) * t541 + t453 * t527 + t429) * MDP(13) + (-t428 - t513) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t585; (-t368 * t387 - t448 * t531 + t593) * MDP(27) + (-t387 * t566 + t448 * t535 + t597) * MDP(28) + (t441 + (t428 - t513) * qJD(3)) * MDP(10) + (t437 * t428 + t478 * t515) * MDP(14) + (t424 * t313 + t489 * t454 + t531 * t355 - (t454 * t535 - t458 * t494) * t592 + t571) * MDP(34) + (t424 * t312 + t489 * t458 + t531 * t357 - (t454 * t494 + t458 * t535) * t592 + t570) * MDP(35) + (-t482 * t562 - t576 + t532 * t451 + (qJD(4) * t511 - t379) * t460 - t504) * MDP(21) + (t483 * t562 - t577 - t502 * t451 + (t456 * t511 - t393) * qJD(4) + t505) * MDP(20) + t583 + (-t437 * t429 - t471) * MDP(13) + (t551 * t592 + t596) * MDP(30) + (t476 + t553) * MDP(32) + (-t428 ^ 2 + t429 ^ 2) * MDP(9) + t429 * t428 * MDP(8) + t595; (-t451 * t486 + t467 - t577) * MDP(20) + (t451 * t503 - t567 - t576) * MDP(21) + (-t368 * t561 + t323 * t448 - t575 + (t455 * t510 - t554) * qJD(5) - t508) * MDP(27) + (-t566 * t561 + t324 * t448 + t588 + (qJD(5) * t510 - t315) * t459 - t507) * MDP(28) + (t446 * t313 + t490 * t454 + t492 * t355 - (t454 * t491 - t458 * t493) * t592 + t571) * MDP(34) + (t446 * t312 + t490 * t458 + t492 * t357 - (t454 * t493 + t458 * t491) * t592 + t570) * MDP(35) + t583 + t590; (t322 * t448 + t593) * MDP(27) + (t321 * t448 + t597) * MDP(28) + (-pkin(5) * t313 + (-t321 * t454 + t343 * t458) * t592 - t322 * t355 + t454 * t589 - t534 * pkin(10) + t571) * MDP(34) + (-pkin(5) * t312 - (t321 * t458 + t343 * t454) * t592 - t322 * t357 + t319 * t594 + t573 * pkin(10) + t570) * MDP(35) + t590; t357 * t355 * MDP(29) + (-t355 ^ 2 + t357 ^ 2) * MDP(30) + (-t355 * t592 + t533) * MDP(31) + (-t357 * t592 - t557) * MDP(32) + t333 * MDP(33) + (-t300 * t454 - t307 * t592 - t319 * t357 + t304) * MDP(34) + (-t300 * t458 - t305 * t454 + t319 * t355 + t488 * t592) * MDP(35) + (-MDP(31) * t545 - MDP(32) * t357 - MDP(34) * t307 + MDP(35) * t488) * qJD(6);];
tauc  = t1;
