% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:52
% EndTime: 2019-03-09 08:40:00
% DurationCPUTime: 5.09s
% Computational Cost: add. (3978->490), mult. (9617->631), div. (0->0), fcn. (6286->6), ass. (0->198)
t439 = cos(pkin(9));
t442 = sin(qJ(2));
t508 = qJD(1) * t442;
t481 = t439 * t508;
t438 = sin(pkin(9));
t504 = qJD(2) * t438;
t391 = t481 + t504;
t443 = cos(qJ(5));
t483 = t438 * t508;
t494 = t439 * qJD(2);
t390 = -t483 + t494;
t441 = sin(qJ(5));
t528 = t390 * t441;
t337 = -t391 * t443 + t528;
t545 = t337 ^ 2;
t393 = t438 * t441 + t439 * t443;
t377 = t393 * qJD(5);
t444 = cos(qJ(2));
t455 = t393 * t444;
t513 = -qJD(1) * t455 - t377;
t507 = qJD(1) * t444;
t480 = t439 * t507;
t482 = t438 * t507;
t496 = qJD(5) * t443;
t497 = qJD(5) * t441;
t512 = t438 * t496 - t439 * t497 - t441 * t480 + t443 * t482;
t492 = qJD(1) * qJD(2);
t544 = -0.2e1 * t492;
t543 = MDP(4) * t442;
t437 = t444 ^ 2;
t542 = MDP(5) * (t442 ^ 2 - t437);
t425 = qJD(5) + t507;
t460 = t443 * t390 + t391 * t441;
t541 = t425 * t460;
t485 = -pkin(7) * t438 - pkin(3);
t522 = t439 * t444;
t489 = pkin(8) * t522;
t450 = -t489 + (-pkin(4) + t485) * t442;
t465 = pkin(2) * t442 - qJ(3) * t444;
t397 = t465 * qJD(1);
t527 = t397 * t439;
t319 = qJD(1) * t450 - t527;
t379 = t438 * t397;
t426 = qJ(4) * t508;
t523 = t439 * t442;
t524 = t438 * t444;
t457 = -pkin(7) * t523 + pkin(8) * t524;
t332 = qJD(1) * t457 + t379 + t426;
t533 = -pkin(8) + qJ(3);
t408 = t533 * t438;
t409 = t533 * t439;
t459 = t408 * t443 - t409 * t441;
t540 = -qJD(3) * t393 - qJD(5) * t459 + t441 * t319 + t443 * t332;
t354 = t408 * t441 + t409 * t443;
t525 = t438 * t443;
t394 = -t439 * t441 + t525;
t539 = qJD(3) * t394 - qJD(5) * t354 - t319 * t443 + t332 * t441;
t532 = qJ(3) * t442;
t404 = -pkin(2) * t444 - pkin(1) - t532;
t421 = pkin(7) * t524;
t435 = t444 * pkin(3);
t334 = pkin(4) * t444 + t421 + t435 + (-pkin(8) * t442 - t404) * t439;
t422 = pkin(7) * t522;
t363 = t438 * t404 + t422;
t356 = -qJ(4) * t444 + t363;
t526 = t438 * t442;
t343 = pkin(8) * t526 + t356;
t538 = t441 * t334 + t443 * t343;
t499 = qJD(4) * t438;
t429 = pkin(7) * t507;
t511 = qJ(4) * t480 - t429;
t534 = -pkin(3) - pkin(4);
t473 = -t482 * t534 + t499 - t511;
t498 = qJD(4) * t444;
t503 = qJD(2) * t442;
t536 = qJ(4) * t503 - t498;
t491 = MDP(24) + MDP(26);
t490 = MDP(25) - MDP(28);
t375 = qJD(2) * t465 - qJD(3) * t442;
t529 = t375 * t439;
t312 = qJD(2) * t450 - t529;
t367 = t438 * t375;
t313 = qJD(2) * t457 + t367 + t536;
t535 = -qJD(5) * t538 + t312 * t443 - t313 * t441;
t388 = t391 ^ 2;
t440 = qJD(2) * pkin(2);
t383 = t404 * qJD(1);
t411 = qJD(2) * qJ(3) + t429;
t344 = t383 * t439 - t438 * t411;
t331 = pkin(3) * t507 + qJD(4) - t344;
t531 = t331 * t442;
t345 = t438 * t383 + t439 * t411;
t333 = -qJ(4) * t507 + t345;
t530 = t333 * t442;
t445 = qJD(2) ^ 2;
t521 = t442 * t445;
t520 = t444 * t445;
t446 = qJD(1) ^ 2;
t519 = t444 * t446;
t518 = -qJ(6) * t508 + t540;
t517 = pkin(5) * t508 + t539;
t516 = -pkin(5) * t512 + qJ(6) * t513 + qJD(6) * t394 - t473;
t366 = t375 * qJD(1);
t428 = pkin(7) * t508;
t401 = (qJD(3) - t428) * qJD(2);
t329 = t438 * t366 + t439 * t401;
t479 = t444 * t494;
t510 = -qJ(4) * t479 - qJD(4) * t523;
t506 = qJD(2) * t459;
t505 = qJD(2) * t354;
t502 = qJD(2) * t444;
t501 = qJD(3) * t391;
t500 = qJD(3) * t439;
t495 = qJD(6) * t425;
t301 = pkin(4) * t507 - pkin(8) * t391 + t331;
t308 = -pkin(8) * t390 + t333;
t283 = t301 * t443 - t308 * t441;
t493 = qJD(6) - t283;
t402 = -t439 * pkin(3) - t438 * qJ(4) - pkin(2);
t488 = pkin(7) * t503;
t478 = t442 * t492;
t487 = qJ(4) * t478 + t329;
t477 = t444 * t492;
t413 = t438 * t477;
t469 = t439 * t477;
t486 = -t390 * t496 + t441 * t413 + t443 * t469;
t484 = pkin(3) * t438 + pkin(7);
t403 = qJD(3) + t428 - t440;
t476 = MDP(23) * t508;
t475 = -t403 - t440;
t474 = pkin(1) * t544;
t328 = t366 * t439 - t438 * t401;
t362 = t404 * t439 - t421;
t380 = t439 * pkin(4) - t402;
t361 = pkin(3) * t482 - t511;
t472 = t361 + t499;
t299 = (t442 * t534 - t489) * t492 - t328;
t300 = (pkin(8) * t504 - qJD(4)) * t507 + t487;
t471 = -t443 * t299 + t441 * t300 + t301 * t497 + t308 * t496;
t424 = pkin(7) * t477;
t320 = pkin(3) * t413 - qJ(4) * t469 - t391 * qJD(4) + t424;
t470 = qJ(6) * t478;
t468 = t438 * t534 - pkin(7);
t467 = t485 * t442;
t423 = pkin(5) * t478;
t275 = t423 + t471;
t325 = -t390 * pkin(3) - t391 * qJ(4) + t403;
t284 = t301 * t441 + t308 * t443;
t462 = t334 * t443 - t343 * t441;
t359 = -pkin(7) * t481 + t379;
t351 = -t439 * t488 + t367;
t456 = t284 * t425 - t471;
t454 = -t441 * t299 - t443 * t300 - t301 * t496 + t308 * t497;
t453 = t441 * t312 + t443 * t313 + t334 * t496 - t343 * t497;
t302 = t391 * t497 - t486;
t452 = -t443 * t413 + t441 * t469;
t419 = qJ(4) * t523;
t355 = t442 * t468 + t419;
t306 = pkin(4) * t390 - t325;
t451 = t441 * t490 - t443 * t491;
t449 = t283 * t425 + t454;
t330 = t468 * t502 - t510;
t309 = -pkin(4) * t413 - t320;
t303 = -qJD(5) * t337 + t452;
t276 = pkin(5) * t303 + qJ(6) * t302 + qJD(6) * t337 + t309;
t280 = -pkin(5) * t425 + t493;
t281 = qJ(6) * t425 + t284;
t447 = (-t337 * t441 - t443 * t460) * MDP(27) + (t280 * t441 + t281 * t443) * MDP(29) + (-t441 * t491 - t443 * t490) * t425;
t412 = qJD(3) * t482;
t373 = t390 * t507;
t372 = t390 * t500;
t371 = t393 * t442;
t370 = t441 * t523 - t442 * t525;
t368 = t442 * t484 - t419;
t358 = pkin(7) * t483 + t527;
t357 = -t362 + t435;
t350 = t438 * t488 + t529;
t349 = qJD(1) * t467 - t527;
t348 = t359 + t426;
t347 = t484 * t502 + t510;
t336 = qJD(2) * t467 - t529;
t327 = t351 + t536;
t324 = qJD(5) * t394 * t442 + qJD(2) * t455;
t323 = t377 * t442 + t441 * t479 - t502 * t525;
t318 = pkin(5) * t393 - qJ(6) * t394 + t380;
t315 = -pkin(3) * t478 - t328;
t307 = -qJD(1) * t498 + t487;
t294 = pkin(5) * t370 - qJ(6) * t371 + t355;
t291 = -pkin(5) * t337 + qJ(6) * t460;
t289 = -pkin(5) * t444 - t462;
t288 = qJ(6) * t444 + t538;
t285 = -t302 + t541;
t282 = pkin(5) * t460 + qJ(6) * t337 + t306;
t279 = pkin(5) * t323 - qJ(6) * t324 - qJD(6) * t371 + t330;
t278 = pkin(5) * t503 - t535;
t277 = -qJ(6) * t503 + qJD(6) * t444 + t453;
t274 = -t454 - t470 + t495;
t1 = [(-t303 * t444 - t323 * t425 + (qJD(1) * t370 + t460) * t503) * MDP(22) + (-t275 * t444 + t276 * t370 - t278 * t425 + t279 * t460 + t282 * t323 + t294 * t303 + (qJD(1) * t289 + t280) * t503) * MDP(26) + 0.2e1 * t477 * t543 + (t274 * t444 - t276 * t371 + t277 * t425 + t279 * t337 - t282 * t324 + t294 * t302 + (-qJD(1) * t288 - t281) * t503) * MDP(28) + (-t302 * t371 - t324 * t337) * MDP(19) + (-t453 * t425 + t454 * t444 - t330 * t337 - t355 * t302 + t309 * t371 + t306 * t324 + (qJD(1) * t538 + t284) * t503) * MDP(25) + (-t320 * t523 - t347 * t391 + (-qJD(1) * t327 - t307) * t444 + (-t325 * t522 + t530 + (t356 * t442 - t368 * t522) * qJD(1)) * qJD(2)) * MDP(17) + (t320 * t526 - t347 * t390 + (qJD(1) * t336 + t315) * t444 + (t325 * t524 - t531 + (-t357 * t442 + t368 * t524) * qJD(1)) * qJD(2)) * MDP(15) + (-t274 * t370 + t275 * t371 - t277 * t460 - t278 * t337 + t280 * t324 - t281 * t323 - t288 * t303 - t289 * t302) * MDP(27) + (t302 * t370 - t303 * t371 + t323 * t337 - t324 * t460) * MDP(20) + (-t302 * t444 + t324 * t425 + (-qJD(1) * t371 + t337) * t503) * MDP(21) + (t535 * t425 - t471 * t444 + t330 * t460 + t355 * t303 + t309 * t370 + t306 * t323 + (-qJD(1) * t462 - t283) * t503) * MDP(24) + (-pkin(7) * t520 + t442 * t474) * MDP(9) - MDP(7) * t521 + (pkin(7) * t521 + t444 * t474) * MDP(10) + (-t350 * t391 + t351 * t390 + (-t328 * t439 - t329 * t438) * t442 + (-t344 * t439 - t345 * t438 + (-t362 * t439 - t363 * t438) * qJD(1)) * t502) * MDP(13) + (t327 * t390 + t336 * t391 + (-t307 * t438 + t315 * t439) * t442 + (t331 * t439 - t333 * t438 + (-t356 * t438 + t357 * t439) * qJD(1)) * t502) * MDP(16) + (t328 * t362 + t329 * t363 + t344 * t350 + t345 * t351 + (t403 + t428) * pkin(7) * t502) * MDP(14) + ((qJD(1) * t351 + t329) * t444 + ((pkin(7) * t391 + t403 * t439) * t444 + (-t345 + (-t363 + 0.2e1 * t422) * qJD(1)) * t442) * qJD(2)) * MDP(12) + ((-qJD(1) * t350 - t328) * t444 + ((-pkin(7) * t390 + t403 * t438) * t444 + (t344 + (t362 + 0.2e1 * t421) * qJD(1)) * t442) * qJD(2)) * MDP(11) + t542 * t544 + (t274 * t288 + t275 * t289 + t276 * t294 + t277 * t281 + t278 * t280 + t279 * t282) * MDP(29) + (-t425 - t507) * MDP(23) * t503 + (t307 * t356 + t315 * t357 + t320 * t368 + t325 * t347 + t327 * t333 + t331 * t336) * MDP(18) + MDP(6) * t520; -t519 * t543 + t446 * t542 + (t412 + ((-qJ(3) * t504 - t344) * t442 + (t358 + t475 * t438 + (t390 - t494) * pkin(7)) * t444) * qJD(1)) * MDP(11) + ((-qJ(3) * t494 + t345) * t442 + (-t359 + (-t391 + t504) * pkin(7) + (qJD(3) + t475) * t439) * t444) * qJD(1) * MDP(12) + (t358 * t391 - t359 * t390 + t372 + (t344 * t507 + t329) * t439 + (t345 * t507 - t328 + t501) * t438) * MDP(13) + (-t344 * t358 - t345 * t359 + (-t344 * t438 + t345 * t439) * qJD(3) + (-t328 * t438 + t329 * t439) * qJ(3) + t475 * t429) * MDP(14) + (-t320 * t439 + t412 + t472 * t390 + (t531 - t349 * t444 + (-t325 * t444 + (t402 * t444 - t532) * qJD(2)) * t438) * qJD(1)) * MDP(15) + (-t348 * t390 - t349 * t391 + t372 + (-t331 * t507 + t307) * t439 + (t333 * t507 + t315 + t501) * t438) * MDP(16) + (-t320 * t438 + t472 * t391 + (-t530 + t348 * t444 + (qJ(3) * t503 + (-qJD(2) * t402 - qJD(3) + t325) * t444) * t439) * qJD(1)) * MDP(17) + (qJ(3) * t307 * t439 + t320 * t402 - t325 * t361 - t331 * t349 + (-t348 + t500) * t333 + (qJ(3) * t315 + qJD(3) * t331 - qJD(4) * t325) * t438) * MDP(18) + (-t302 * t394 - t337 * t513) * MDP(19) + (t302 * t393 - t303 * t394 + t337 * t512 - t460 * t513) * MDP(20) + (t513 * t425 + (-qJD(2) * t394 - t337) * t508) * MDP(21) + (-t512 * t425 + (qJD(2) * t393 - t460) * t508) * MDP(22) + t425 * t476 + (t380 * t303 + t309 * t393 + t539 * t425 + t473 * t460 + t512 * t306 + (t283 - t506) * t508) * MDP(24) + (-t380 * t302 + t309 * t394 + t540 * t425 - t473 * t337 + t513 * t306 + (-t284 + t505) * t508) * MDP(25) + (t276 * t393 + t303 * t318 + t517 * t425 - t516 * t460 + t512 * t282 + (-t280 - t506) * t508) * MDP(26) + (-t274 * t393 + t275 * t394 + t280 * t513 - t281 * t512 + t302 * t459 - t303 * t354 + t337 * t517 + t460 * t518) * MDP(27) + (-t276 * t394 + t302 * t318 - t518 * t425 - t516 * t337 - t513 * t282 + (t281 - t505) * t508) * MDP(28) + (t274 * t354 - t275 * t459 + t276 * t318 - t280 * t517 - t281 * t518 - t282 * t516) * MDP(29) + (MDP(9) * t442 * t446 + MDP(10) * t519) * pkin(1); (t344 * t391 - t345 * t390 + t424) * MDP(14) + (-t331 * t391 - t333 * t390 + t320) * MDP(18) + (t460 ^ 2 + t545) * MDP(27) + (-t280 * t337 - t281 * t460 - t276) * MDP(29) + (MDP(12) - MDP(17)) * (-t373 + t469) + (t391 * t451 + t491 * t528) * qJD(5) + (MDP(11) + MDP(15)) * (-t391 * t507 + t413) + t490 * (-t486 + t541) + t491 * (t337 * t425 - t452) + (MDP(13) + MDP(16)) * (-t390 ^ 2 - t388); t373 * MDP(16) + (-t437 * t446 - t388) * MDP(17) - t328 * MDP(18) + (t302 * t443 - t303 * t441) * MDP(27) + (t274 * t441 - t275 * t443) * MDP(29) + (-t390 * MDP(15) + t325 * MDP(18) - t282 * MDP(29) + t337 * t490 - t460 * t491) * t391 + t447 * qJD(5) + ((-MDP(18) * pkin(3) - MDP(15) + t451) * t503 + (MDP(16) * t494 + t333 * MDP(18) + t447) * t444) * qJD(1); t285 * MDP(21) + (t390 * t497 - t391 * t496 - t452) * MDP(22) - qJD(2) * t476 + t456 * MDP(24) + t449 * MDP(25) + (-0.2e1 * t423 + t456) * MDP(26) + (pkin(5) * t302 - qJ(6) * t303) * MDP(27) + (-t449 - 0.2e1 * t470 + 0.2e1 * t495) * MDP(28) + (-pkin(5) * t275 + qJ(6) * t274 - t280 * t284 + t281 * t493 - t282 * t291) * MDP(29) - (t425 * MDP(22) - t306 * MDP(24) - t282 * MDP(26) + (t281 - t284) * MDP(27) + t291 * MDP(28) - MDP(20) * t337) * t337 + (-t337 * MDP(19) + t306 * MDP(25) - t291 * MDP(26) + (t280 - t493) * MDP(27) - t282 * MDP(28) - MDP(20) * t460) * t460; (-t337 * t460 + t478) * MDP(26) + t285 * MDP(27) + (-t425 ^ 2 - t545) * MDP(28) + (-t281 * t425 - t282 * t337 + t275) * MDP(29);];
tauc  = t1;
