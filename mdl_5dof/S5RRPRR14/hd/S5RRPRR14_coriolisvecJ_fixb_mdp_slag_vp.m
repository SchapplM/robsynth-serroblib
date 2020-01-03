% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR14_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR14_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:49
% EndTime: 2019-12-31 20:39:01
% DurationCPUTime: 6.13s
% Computational Cost: add. (4834->436), mult. (13194->624), div. (0->0), fcn. (10630->10), ass. (0->186)
t442 = cos(pkin(5));
t513 = qJD(1) * t442;
t431 = qJD(2) + t513;
t439 = sin(pkin(10));
t441 = cos(pkin(10));
t445 = sin(qJ(2));
t440 = sin(pkin(5));
t514 = qJD(1) * t440;
t498 = t445 * t514;
t386 = t431 * t441 - t439 * t498;
t387 = t431 * t439 + t441 * t498;
t444 = sin(qJ(4));
t447 = cos(qJ(4));
t341 = -t447 * t386 + t387 * t444;
t338 = qJD(5) + t341;
t547 = t338 ^ 2;
t448 = cos(qJ(2));
t512 = qJD(1) * t448;
t497 = t440 * t512;
t424 = -qJD(4) + t497;
t546 = t341 * t424;
t522 = t447 * t441;
t407 = t439 * t444 - t522;
t525 = t440 * t448;
t455 = t407 * t525;
t545 = -qJD(1) * t455 + t407 * qJD(4);
t408 = t439 * t447 + t441 * t444;
t456 = t408 * t525;
t516 = -qJD(1) * t456 + t408 * qJD(4);
t436 = t440 ^ 2;
t504 = qJD(1) * qJD(2);
t544 = -0.2e1 * t436 * t504;
t543 = MDP(5) * (t445 ^ 2 - t448 ^ 2);
t393 = pkin(7) * t525 + (pkin(1) * t445 + qJ(3)) * t442;
t394 = (-pkin(2) * t448 - qJ(3) * t445 - pkin(1)) * t440;
t345 = -t393 * t439 + t441 * t394;
t526 = t440 * t445;
t402 = t439 * t442 + t441 * t526;
t316 = -pkin(3) * t525 - pkin(8) * t402 + t345;
t346 = t441 * t393 + t439 * t394;
t401 = t439 * t526 - t442 * t441;
t328 = -pkin(8) * t401 + t346;
t542 = t444 * t316 + t447 * t328;
t475 = pkin(2) * t445 - qJ(3) * t448;
t397 = t475 * t514;
t502 = pkin(1) * t513;
t398 = -pkin(7) * t498 + t448 * t502;
t349 = t441 * t397 - t398 * t439;
t524 = t441 * t448;
t459 = (pkin(3) * t445 - pkin(8) * t524) * t440;
t329 = qJD(1) * t459 + t349;
t350 = t439 * t397 + t441 * t398;
t482 = t439 * t497;
t335 = -pkin(8) * t482 + t350;
t536 = pkin(8) + qJ(3);
t419 = t536 * t439;
t420 = t536 * t441;
t465 = -t419 * t447 - t420 * t444;
t541 = -qJD(3) * t407 + qJD(4) * t465 - t444 * t329 - t447 * t335;
t373 = -t419 * t444 + t420 * t447;
t540 = qJD(3) * t408 + qJD(4) * t373 + t329 * t447 - t335 * t444;
t466 = t386 * t444 + t387 * t447;
t539 = qJD(4) * t466;
t399 = pkin(7) * t497 + t445 * t502;
t376 = qJ(3) * t431 + t399;
t381 = qJD(1) * t394;
t330 = -t376 * t439 + t441 * t381;
t302 = -pkin(3) * t497 - pkin(8) * t387 + t330;
t331 = t441 * t376 + t439 * t381;
t310 = pkin(8) * t386 + t331;
t286 = t302 * t444 + t310 * t447;
t379 = (qJD(2) * t475 - qJD(3) * t445) * t440;
t363 = qJD(1) * t379;
t493 = t440 * t504;
t480 = t445 * t493;
t501 = pkin(1) * qJD(2) * t442;
t484 = qJD(1) * t501;
t462 = -pkin(7) * t480 + t448 * t484;
t364 = qJD(3) * t431 + t462;
t325 = t441 * t363 - t364 * t439;
t454 = qJD(2) * t459;
t303 = qJD(1) * t454 + t325;
t326 = t439 * t363 + t441 * t364;
t479 = t448 * t493;
t464 = t439 * t479;
t311 = -pkin(8) * t464 + t326;
t450 = -qJD(4) * t286 + t447 * t303 - t444 * t311;
t276 = -pkin(4) * t480 - t450;
t538 = t338 * (pkin(4) * t466 + pkin(9) * t338) + t276;
t511 = qJD(2) * t445;
t496 = t440 * t511;
t463 = -pkin(7) * t496 + t448 * t501;
t385 = qJD(3) * t442 + t463;
t333 = t441 * t379 - t385 * t439;
t314 = t333 + t454;
t334 = t439 * t379 + t441 * t385;
t495 = qJD(2) * t525;
t481 = t439 * t495;
t324 = -pkin(8) * t481 + t334;
t537 = -qJD(4) * t542 + t314 * t447 - t324 * t444;
t509 = qJD(4) * t447;
t518 = t386 * t509 + t479 * t522;
t308 = (-qJD(4) * t387 - t464) * t444 + t518;
t443 = sin(qJ(5));
t446 = cos(qJ(5));
t507 = qJD(5) * t446;
t499 = t446 * t308 - t424 * t507 + t443 * t480;
t508 = qJD(5) * t443;
t288 = -t466 * t508 + t499;
t535 = t288 * t443;
t530 = t466 * t443;
t321 = t446 * t424 + t530;
t534 = t321 * t338;
t323 = -t424 * t443 + t446 * t466;
t533 = t323 * t338;
t532 = t330 * t445;
t531 = t331 * t445;
t529 = t408 * t446;
t449 = qJD(1) ^ 2;
t528 = t436 * t449;
t527 = t439 * t448;
t453 = qJD(2) * t456;
t309 = qJD(1) * t453 + t539;
t523 = t443 * t309;
t305 = t446 * t309;
t519 = pkin(4) * t498 + t540;
t392 = pkin(7) * t479 + t445 * t484;
t400 = pkin(7) * t495 + t445 * t501;
t510 = qJD(4) * t444;
t506 = qJD(2) - t431;
t503 = pkin(1) * t528;
t500 = t443 * t525;
t360 = pkin(3) * t464 + t392;
t367 = pkin(3) * t482 + t399;
t368 = pkin(3) * t481 + t400;
t435 = -pkin(3) * t441 - pkin(2);
t461 = t302 * t509 + t444 * t303 - t310 * t510 + t447 * t311;
t275 = pkin(9) * t480 + t461;
t282 = pkin(4) * t309 - pkin(9) * t308 + t360;
t491 = -t275 * t443 + t446 * t282;
t489 = t308 * t443 - t446 * t480;
t488 = t443 * t545 - t446 * t498;
t487 = t443 * t498 + t446 * t545;
t486 = t338 * t446;
t483 = t436 * t445 * t448 * MDP(4);
t478 = pkin(1) * t544;
t357 = pkin(4) * t407 - pkin(9) * t408 + t435;
t477 = pkin(9) * t498 - qJD(5) * t357 - t541;
t476 = -pkin(4) * t516 - pkin(9) * t545 + qJD(5) * t373 + t367;
t474 = t275 * t446 + t282 * t443;
t284 = -pkin(9) * t424 + t286;
t369 = -pkin(2) * t431 + qJD(3) - t398;
t339 = -pkin(3) * t386 + t369;
t290 = pkin(4) * t341 - pkin(9) * t466 + t339;
t278 = t284 * t446 + t290 * t443;
t473 = t284 * t443 - t290 * t446;
t292 = -pkin(9) * t525 + t542;
t353 = t447 * t401 + t402 * t444;
t354 = -t401 * t444 + t402 * t447;
t396 = pkin(7) * t526 + (-pkin(1) * t448 - pkin(2)) * t442;
t355 = pkin(3) * t401 + t396;
t297 = pkin(4) * t353 - pkin(9) * t354 + t355;
t472 = t292 * t446 + t297 * t443;
t471 = -t292 * t443 + t297 * t446;
t285 = t302 * t447 - t310 * t444;
t469 = t316 * t447 - t328 * t444;
t336 = t354 * t443 + t446 * t525;
t460 = t444 * t314 + t316 * t509 + t447 * t324 - t328 * t510;
t458 = t408 * t507 - t488;
t457 = -t408 * t508 - t487;
t283 = pkin(4) * t424 - t285;
t452 = -pkin(9) * t309 + (t283 + t285) * t338;
t451 = -qJ(3) * t511 + (-pkin(2) * qJD(2) + qJD(3) - t369) * t448;
t337 = t354 * t446 - t500;
t319 = qJD(4) * t354 + t453;
t318 = -qJD(2) * t455 - qJD(4) * t353;
t294 = -qJD(5) * t500 + t318 * t443 + t354 * t507 - t446 * t496;
t293 = -qJD(5) * t336 + t318 * t446 + t443 * t496;
t291 = pkin(4) * t525 - t469;
t289 = qJD(5) * t323 + t489;
t287 = pkin(4) * t319 - pkin(9) * t318 + t368;
t280 = -pkin(4) * t496 - t537;
t279 = pkin(9) * t496 + t460;
t274 = -qJD(5) * t278 + t491;
t273 = -qJD(5) * t473 + t474;
t1 = [t543 * t544 + (-t392 * t442 - t400 * t431 + t445 * t478) * MDP(9) + (-t431 * t463 - t442 * t462 + t448 * t478) * MDP(10) + (-t386 * t400 + t392 * t401 + ((-qJD(1) * t333 - t325) * t448 + (t369 * t527 + t532 + (t345 * t445 + t396 * t527) * qJD(1)) * qJD(2)) * t440) * MDP(11) + (t387 * t400 + t392 * t402 + ((qJD(1) * t334 + t326) * t448 + (t369 * t524 - t531 + (-t346 * t445 + t396 * t524) * qJD(1)) * qJD(2)) * t440) * MDP(12) + (-t325 * t402 - t326 * t401 - t333 * t387 + t334 * t386 + (-t330 * t441 - t331 * t439 + (-t345 * t441 - t346 * t439) * qJD(1)) * t495) * MDP(13) + (t325 * t345 + t326 * t346 + t330 * t333 + t331 * t334 + t369 * t400 + t392 * t396) * MDP(14) + (t308 * t354 + t318 * t466) * MDP(15) + (-t308 * t353 - t309 * t354 - t318 * t341 - t319 * t466) * MDP(16) + (-t318 * t424 + (-t308 * t448 + (qJD(1) * t354 + t466) * t511) * t440) * MDP(17) + (t319 * t424 + (t309 * t448 + (-qJD(1) * t353 - t341) * t511) * t440) * MDP(18) + (-t424 * t440 - t436 * t512) * MDP(19) * t511 + (-t537 * t424 + t368 * t341 + t355 * t309 + t360 * t353 + t339 * t319 + (-t450 * t448 + (qJD(1) * t469 + t285) * t511) * t440) * MDP(20) + (t460 * t424 + t368 * t466 + t355 * t308 + t360 * t354 + t339 * t318 + (t461 * t448 + (-qJD(1) * t542 - t286) * t511) * t440) * MDP(21) + (t288 * t337 + t293 * t323) * MDP(22) + (-t288 * t336 - t289 * t337 - t293 * t321 - t294 * t323) * MDP(23) + (t288 * t353 + t293 * t338 + t309 * t337 + t319 * t323) * MDP(24) + (-t289 * t353 - t294 * t338 - t309 * t336 - t319 * t321) * MDP(25) + (t309 * t353 + t319 * t338) * MDP(26) + ((-qJD(5) * t472 - t279 * t443 + t287 * t446) * t338 + t471 * t309 + t274 * t353 - t473 * t319 + t280 * t321 + t291 * t289 + t276 * t336 + t283 * t294) * MDP(27) + (-(qJD(5) * t471 + t279 * t446 + t287 * t443) * t338 - t472 * t309 - t273 * t353 - t278 * t319 + t280 * t323 + t291 * t288 + t276 * t337 + t283 * t293) * MDP(28) + 0.2e1 * t483 * t504 + (MDP(6) * t495 - MDP(7) * t496) * (t431 + t513); -t449 * t483 + t528 * t543 + t506 * MDP(6) * t497 + (t399 * t431 + t445 * t503 - t392) * MDP(9) + (t398 * t431 + t448 * t503 - t462) * MDP(10) + (t386 * t399 - t392 * t441 + (t349 * t448 + t439 * t451 - t532) * t514) * MDP(11) + (-t387 * t399 + t392 * t439 + (-t350 * t448 + t441 * t451 + t531) * t514) * MDP(12) + (t349 * t387 - t350 * t386 + (qJD(3) * t386 + t330 * t497 + t326) * t441 + (qJD(3) * t387 + t331 * t497 - t325) * t439) * MDP(13) + (-pkin(2) * t392 - t330 * t349 - t331 * t350 - t369 * t399 + (-t330 * t439 + t331 * t441) * qJD(3) + (-t325 * t439 + t326 * t441) * qJ(3)) * MDP(14) + (t308 * t408 - t466 * t545) * MDP(15) + (-t308 * t407 - t309 * t408 + t341 * t545 - t466 * t516) * MDP(16) + (t435 * t309 + t516 * t339 - t367 * t341 + t360 * t407) * MDP(20) + (t435 * t308 - t339 * t545 + t360 * t408 - t367 * t466) * MDP(21) + (t288 * t529 + t323 * t457) * MDP(22) + (t488 * t323 + t487 * t321 + (-t535 - t289 * t446 + (t321 * t443 - t323 * t446) * qJD(5)) * t408) * MDP(23) + (t288 * t407 + t305 * t408 + t323 * t516 + t338 * t457) * MDP(24) + (-t289 * t407 - t321 * t516 - t338 * t458 - t408 * t523) * MDP(25) + (t309 * t407 + t338 * t516) * MDP(26) + ((t357 * t446 - t373 * t443) * t309 + t274 * t407 - t465 * t289 + t276 * t443 * t408 + (t443 * t477 - t446 * t476) * t338 + t519 * t321 - t516 * t473 + t458 * t283) * MDP(27) + (-(t357 * t443 + t373 * t446) * t309 - t273 * t407 - t465 * t288 + t276 * t529 + (t443 * t476 + t446 * t477) * t338 + t519 * t323 - t516 * t278 + t457 * t283) * MDP(28) + (MDP(17) * t545 + t516 * MDP(18) + MDP(20) * t540 + MDP(21) * t541) * t424 + (t424 * MDP(19) - t506 * MDP(7) + (qJD(2) * t408 - t466) * MDP(17) + (-qJD(2) * t407 + t341) * MDP(18) + (qJD(2) * t465 - t285) * MDP(20) + (-qJD(2) * t373 + t286) * MDP(21)) * t498; (-t386 ^ 2 - t387 ^ 2) * MDP(13) + (t330 * t387 - t331 * t386 + t392) * MDP(14) + (t386 * t510 + t387 * t509 - t424 * t466) * MDP(20) + (-t387 * t510 + t518 + t546) * MDP(21) + (-t321 * t466 + t305) * MDP(27) + (-t323 * t466 - t523) * MDP(28) + (-t387 * MDP(11) - t386 * MDP(12) + ((MDP(20) * t444 + MDP(12)) * t441 + (MDP(20) * t447 - MDP(21) * t444 + MDP(11)) * t439) * qJD(2)) * t497 - (MDP(27) * t443 + MDP(28) * t446) * t547; -t341 ^ 2 * MDP(16) + (t308 - t546) * MDP(17) + (-t408 * t479 - t539) * MDP(18) + MDP(19) * t480 + (-t286 * t424 + t450) * MDP(20) + (-t285 * t424 + t339 * t341 - t461) * MDP(21) + (t323 * t486 + t535) * MDP(22) + ((t288 - t534) * t446 + (-t289 - t533) * t443) * MDP(23) + (t338 * t486 + t523) * MDP(24) + (-t443 * t547 + t305) * MDP(25) + (-pkin(4) * t289 - t286 * t321 + t452 * t443 - t446 * t538) * MDP(27) + (-pkin(4) * t288 - t286 * t323 + t443 * t538 + t452 * t446) * MDP(28) + (MDP(15) * t341 + MDP(16) * t466 - t424 * MDP(18) - t339 * MDP(20) - t323 * MDP(24) + t321 * MDP(25) - t338 * MDP(26) + MDP(27) * t473 + t278 * MDP(28)) * t466; t323 * t321 * MDP(22) + (-t321 ^ 2 + t323 ^ 2) * MDP(23) + (t499 + t534) * MDP(24) + (-t489 + t533) * MDP(25) + t309 * MDP(26) + (t278 * t338 - t283 * t323 + t491) * MDP(27) + (t283 * t321 - t338 * t473 - t474) * MDP(28) + (-MDP(24) * t530 - MDP(25) * t323 - MDP(27) * t278 + MDP(28) * t473) * qJD(5);];
tauc = t1;
