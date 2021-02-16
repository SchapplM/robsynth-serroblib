% Calculate Coriolis joint torque vector for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:37
% EndTime: 2021-01-16 03:28:53
% DurationCPUTime: 5.92s
% Computational Cost: add. (4488->396), mult. (11754->556), div. (0->0), fcn. (9508->12), ass. (0->185)
t452 = cos(qJ(6));
t502 = qJD(6) * t452;
t444 = sin(pkin(12));
t446 = cos(pkin(12));
t454 = cos(qJ(3));
t506 = qJD(2) * t454;
t490 = t446 * t506;
t450 = sin(qJ(3));
t508 = qJD(2) * t450;
t413 = -t444 * t508 + t490;
t453 = cos(qJ(5));
t401 = t453 * t413;
t421 = t444 * t454 + t446 * t450;
t415 = t421 * qJD(2);
t449 = sin(qJ(5));
t368 = -t415 * t449 + t401;
t561 = t368 * t452;
t565 = t502 - t561;
t540 = qJ(4) + pkin(8);
t486 = qJD(3) * t540;
t409 = qJD(4) * t454 - t450 * t486;
t410 = -qJD(4) * t450 - t454 * t486;
t455 = cos(qJ(2));
t445 = sin(pkin(6));
t510 = qJD(1) * t445;
t493 = t455 * t510;
t515 = -t409 * t444 + t410 * t446 + t421 * t493;
t420 = t444 * t450 - t446 * t454;
t514 = t409 * t446 + t410 * t444 + t420 * t493;
t468 = t413 * t449 + t415 * t453;
t448 = sin(qJ(6));
t503 = qJD(6) * t448;
t414 = t421 * qJD(3);
t404 = qJD(2) * t414;
t496 = qJD(2) * qJD(3);
t488 = t450 * t496;
t429 = t444 * t488;
t487 = t454 * t496;
t405 = t446 * t487 - t429;
t504 = qJD(5) * t449;
t328 = qJD(5) * t401 - t404 * t449 + t405 * t453 - t415 * t504;
t441 = qJD(3) + qJD(5);
t516 = t328 * t452 + t441 * t502;
t309 = -t468 * t503 + t516;
t308 = t309 * t452;
t356 = t441 * t448 + t452 * t468;
t536 = t328 * t448;
t310 = qJD(6) * t356 + t536;
t528 = t468 * t448;
t354 = -t441 * t452 + t528;
t564 = -t448 * t310 - t354 * t565 + t308;
t307 = t309 * t448;
t329 = qJD(5) * t468 + t404 * t453 + t405 * t449;
t497 = -qJD(6) + t368;
t322 = t448 * t329;
t517 = -t497 * t502 + t322;
t524 = t452 * t497;
t530 = t368 * t441;
t532 = t468 * t441;
t534 = t356 * t468;
t563 = (-t329 + t532) * MDP(19) - t368 ^ 2 * MDP(17) + (t368 * t524 + t517 - t534) * MDP(25) + (-MDP(16) * t368 + MDP(17) * t468 + MDP(27) * t497) * t468 + (t328 - t530) * MDP(18) + (t356 * t565 + t307) * MDP(23);
t562 = t368 * t448;
t417 = t420 * qJD(3);
t560 = -pkin(9) * t417 - t515;
t559 = -pkin(9) * t414 + t514;
t451 = sin(qJ(2));
t494 = t451 * t510;
t505 = qJD(3) * t450;
t558 = pkin(3) * t505 - t494;
t334 = pkin(5) * t468 - pkin(10) * t368;
t425 = qJD(2) * pkin(8) + t494;
t447 = cos(pkin(6));
t509 = qJD(1) * t447;
t433 = t454 * t509;
t473 = qJD(4) + t493;
t352 = (-t425 * t450 + t433) * qJD(3) + (-qJ(4) * t505 + t454 * t473) * qJD(2);
t492 = t450 * t509;
t353 = (-t425 * t454 - t492) * qJD(3) + (-qJ(4) * qJD(3) * t454 - t450 * t473) * qJD(2);
t315 = -t352 * t444 + t353 * t446;
t312 = -pkin(9) * t405 + t315;
t316 = t352 * t446 + t353 * t444;
t313 = -pkin(9) * t404 + t316;
t480 = qJ(4) * qJD(2) + t425;
t384 = t454 * t480 + t492;
t376 = t444 * t384;
t383 = -t450 * t480 + t433;
t380 = qJD(3) * pkin(3) + t383;
t337 = t380 * t446 - t376;
t541 = pkin(9) * t415;
t325 = qJD(3) * pkin(4) + t337 - t541;
t526 = t446 * t384;
t338 = t380 * t444 + t526;
t542 = pkin(9) * t413;
t330 = t338 + t542;
t287 = t453 * (qJD(5) * t325 + t313) + t312 * t449 - t330 * t504;
t438 = -pkin(3) * t454 - pkin(2);
t406 = qJD(2) * t438 + qJD(4) - t493;
t371 = -pkin(4) * t413 + t406;
t557 = -t368 * t371 - t287;
t553 = MDP(5) * t450;
t552 = MDP(6) * (t450 ^ 2 - t454 ^ 2);
t535 = t354 * t468;
t475 = pkin(4) * t414 + t558;
t324 = t452 * t329;
t550 = -t497 * t503 - t324;
t549 = MDP(10) * t450 + MDP(11) * t454;
t302 = t325 * t449 + t330 * t453;
t288 = qJD(5) * t302 - t312 * t453 + t313 * t449;
t301 = t325 * t453 - t330 * t449;
t299 = -pkin(5) * t441 - t301;
t300 = pkin(10) * t441 + t302;
t314 = -pkin(5) * t368 - pkin(10) * t468 + t371;
t471 = t300 * t448 - t314 * t452;
t548 = -t288 * t452 + t299 * t503 + t468 * t471;
t290 = t300 * t452 + t314 * t448;
t547 = t288 * t448 + t290 * t468 + t299 * t502;
t546 = -t371 * t468 - t288;
t374 = t420 * t453 + t421 * t449;
t375 = -t420 * t449 + t421 * t453;
t395 = pkin(4) * t420 + t438;
t318 = pkin(5) * t374 - pkin(10) * t375 + t395;
t427 = t540 * t450;
t428 = t540 * t454;
t387 = -t427 * t446 - t428 * t444;
t360 = -pkin(9) * t421 + t387;
t388 = -t427 * t444 + t428 * t446;
t361 = -pkin(9) * t420 + t388;
t321 = t360 * t449 + t361 * t453;
t335 = -qJD(5) * t374 - t414 * t449 - t417 * t453;
t470 = t360 * t453 - t361 * t449;
t521 = -qJD(5) * t470 + t449 * t560 - t453 * t559;
t545 = -(qJD(6) * t314 + t287) * t374 + t288 * t375 + t299 * t335 - (-qJD(6) * t318 + t521) * t497 - t321 * t329;
t544 = pkin(3) * t444;
t543 = pkin(3) * t450;
t539 = qJD(2) * pkin(2);
t538 = t299 * t375;
t537 = t318 * t329;
t533 = t356 * t448;
t527 = t445 * t451;
t456 = qJD(3) ^ 2;
t525 = t450 * t456;
t523 = t454 * t456;
t520 = qJD(5) * t321 + t449 * t559 + t453 * t560;
t340 = -t383 * t444 - t526;
t331 = t340 - t542;
t342 = t383 * t446 - t376;
t332 = t342 - t541;
t437 = pkin(3) * t446 + pkin(4);
t465 = t437 * t453 - t449 * t544;
t519 = -qJD(5) * t465 + t331 * t449 + t332 * t453;
t466 = t437 * t449 + t453 * t544;
t513 = qJD(5) * t466 + t331 * t453 - t332 * t449;
t411 = pkin(3) * t488 + qJD(2) * t494;
t507 = qJD(2) * t451;
t501 = qJD(6) * t455;
t439 = pkin(3) * t508;
t491 = qJD(2) * t445 * t455;
t389 = pkin(4) * t415 + t439;
t481 = t497 * t448;
t408 = pkin(10) + t466;
t477 = qJD(6) * t408 + t334 + t389;
t370 = pkin(4) * t404 + t411;
t336 = qJD(5) * t375 + t414 * t453 - t417 * t449;
t476 = pkin(5) * t336 - pkin(10) * t335 + t475;
t472 = -t299 * t368 - t329 * t408;
t418 = t447 * t454 - t450 * t527;
t419 = t447 * t450 + t454 * t527;
t363 = t418 * t446 - t419 * t444;
t364 = t418 * t444 + t419 * t446;
t469 = t363 * t453 - t364 * t449;
t327 = t363 * t449 + t364 * t453;
t467 = -t497 * t562 - t550;
t464 = t335 * t452 - t375 * t503;
t462 = -0.2e1 * qJD(3) * t539;
t457 = qJD(2) ^ 2;
t407 = -pkin(5) - t465;
t382 = -qJD(3) * t419 - t450 * t491;
t381 = qJD(3) * t418 + t454 * t491;
t341 = t381 * t446 + t382 * t444;
t339 = -t381 * t444 + t382 * t446;
t296 = pkin(5) * t329 - pkin(10) * t328 + t370;
t295 = t452 * t296;
t292 = qJD(5) * t327 - t339 * t453 + t341 * t449;
t291 = qJD(5) * t469 + t339 * t449 + t341 * t453;
t1 = [(-t339 * t415 + t341 * t413 - t363 * t405 - t364 * t404) * MDP(14) + (t315 * t363 + t316 * t364 + t337 * t339 + t338 * t341) * MDP(15) + (-(-t291 * t448 - t327 * t502) * t497 - t327 * t322 + t292 * t354 - t469 * t310) * MDP(28) + ((t291 * t452 - t327 * t503) * t497 - t327 * t324 + t292 * t356 - t469 * t309) * MDP(29) + (-MDP(21) * t292 - MDP(22) * t291) * t441 + (MDP(10) * t382 - MDP(11) * t381 + MDP(12) * t339 - MDP(13) * t341) * qJD(3) + ((-t404 * t455 - t413 * t507) * MDP(12) + (-t405 * t455 + t415 * t507) * MDP(13) + (t406 * t507 - t411 * t455) * MDP(15) + (-t329 * t455 - t368 * t507) * MDP(21) + (-t328 * t455 + t468 * t507) * MDP(22) + (-(t448 * t501 + t452 * t507) * t497 - t455 * t324) * MDP(28) + ((t448 * t507 - t452 * t501) * t497 + t455 * t322) * MDP(29) - t549 * t455 * t496 + (-t455 * MDP(4) + (-MDP(10) * t454 + MDP(11) * t450 - MDP(3)) * t451) * t457) * t445; 0.2e1 * t487 * t553 - 0.2e1 * t496 * t552 + MDP(7) * t523 - MDP(8) * t525 + (-pkin(8) * t523 + t450 * t462) * MDP(10) + (pkin(8) * t525 + t454 * t462) * MDP(11) + (t413 * t494 + t404 * t438 + t406 * t414 + t411 * t420 + (-t413 * t543 + t515) * qJD(3)) * MDP(12) + (-t415 * t494 + t405 * t438 - t406 * t417 + t411 * t421 + (t415 * t543 - t514) * qJD(3)) * MDP(13) + (-t315 * t421 - t316 * t420 + t337 * t417 - t338 * t414 - t387 * t405 - t388 * t404 + t413 * t514 - t415 * t515) * MDP(14) + (t315 * t387 + t316 * t388 + t515 * t337 + t514 * t338 + t406 * t558 + t411 * t438) * MDP(15) + (t328 * t375 + t335 * t468) * MDP(16) + (-t328 * t374 - t329 * t375 + t335 * t368 - t336 * t468) * MDP(17) + (t329 * t395 + t336 * t371 - t368 * t475 + t370 * t374) * MDP(21) + (t328 * t395 + t335 * t371 + t370 * t375 + t468 * t475) * MDP(22) + (t308 * t375 + t356 * t464) * MDP(23) + ((-t354 * t452 - t533) * t335 + (-t307 - t310 * t452 + (t354 * t448 - t356 * t452) * qJD(6)) * t375) * MDP(24) + (t309 * t374 + t324 * t375 + t336 * t356 - t464 * t497) * MDP(25) + (-t375 * t322 - t310 * t374 - t336 * t354 - (-t335 * t448 - t375 * t502) * t497) * MDP(26) + (t329 * t374 - t336 * t497) * MDP(27) + (-t471 * t336 + t295 * t374 - t470 * t310 + t520 * t354 + (t537 - t476 * t497 + (-t300 * t374 + t321 * t497 + t538) * qJD(6)) * t452 + t545 * t448) * MDP(28) + (-t290 * t336 - t470 * t309 + t520 * t356 + (-t537 - (-qJD(6) * t300 + t296) * t374 - qJD(6) * t538 - (qJD(6) * t321 - t476) * t497) * t448 + t545 * t452) * MDP(29) + (MDP(18) * t335 - MDP(19) * t336 - MDP(21) * t520 + MDP(22) * t521) * t441; t549 * t539 * qJD(2) + ((t338 + t340) * t415 + (t337 - t342) * t413 + (-t404 * t444 - t405 * t446) * pkin(3)) * MDP(14) + (t368 * t389 - t441 * t513 + t546) * MDP(21) + (-t389 * t468 + t441 * t519 + t557) * MDP(22) + (t407 * t310 + t472 * t448 + t513 * t354 - (t448 * t519 - t452 * t477) * t497 + t548) * MDP(28) + (t407 * t309 + t472 * t452 + t513 * t356 - (t448 * t477 + t452 * t519) * t497 + t547) * MDP(29) + (-t454 * t553 + t552) * t457 + (t467 + t535) * MDP(26) + (t497 * t533 + t564) * MDP(24) + (-t337 * t340 - t338 * t342 + (t315 * t446 + t316 * t444 - t406 * t508) * pkin(3)) * MDP(15) + (-qJD(3) * t340 - t406 * t415 + t413 * t439 + t315) * MDP(12) + (qJD(3) * t342 - t406 * t413 - t415 * t439 - t316) * MDP(13) + t563; -t429 * MDP(13) + (-t413 ^ 2 - t415 ^ 2) * MDP(14) + (t337 * t415 - t338 * t413 + t411) * MDP(15) + (t329 + t532) * MDP(21) + (t328 + t530) * MDP(22) + (t467 - t535) * MDP(28) + (-t497 * t524 - t322 - t534) * MDP(29) + ((t444 * t506 + t446 * t508 + t415) * MDP(12) + (t413 + t490) * MDP(13)) * qJD(3); (t302 * t441 + t546) * MDP(21) + (t301 * t441 + t557) * MDP(22) + (t356 * t481 + t564) * MDP(24) + (-t481 * t497 + t324 + t535) * MDP(26) + (-pkin(5) * t310 + (-t301 * t448 + t334 * t452) * t497 - t302 * t354 - t299 * t562 - t517 * pkin(10) + t548) * MDP(28) + (-pkin(5) * t309 - (t301 * t452 + t334 * t448) * t497 - t302 * t356 - t299 * t561 + t550 * pkin(10) + t547) * MDP(29) + t563; t356 * t354 * MDP(23) + (-t354 ^ 2 + t356 ^ 2) * MDP(24) + (-t354 * t497 + t516) * MDP(25) + (-t356 * t497 - t536) * MDP(26) + t329 * MDP(27) + (-t287 * t448 - t290 * t497 - t299 * t356 + t295) * MDP(28) + (-t287 * t452 - t296 * t448 + t299 * t354 + t471 * t497) * MDP(29) + (-MDP(25) * t528 - MDP(26) * t356 - MDP(28) * t290 + MDP(29) * t471) * qJD(6);];
tauc = t1;
