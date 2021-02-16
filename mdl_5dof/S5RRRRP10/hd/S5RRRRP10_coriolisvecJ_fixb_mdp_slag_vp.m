% Calculate Coriolis joint torque vector for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:36:13
% EndTime: 2021-01-16 00:36:32
% DurationCPUTime: 8.24s
% Computational Cost: add. (5849->505), mult. (15524->685), div. (0->0), fcn. (11664->8), ass. (0->194)
t437 = sin(pkin(5));
t443 = cos(qJ(2));
t509 = qJD(1) * t443;
t428 = t437 * t509;
t464 = t428 - qJD(3);
t440 = sin(qJ(2));
t546 = cos(pkin(5));
t486 = t546 * qJD(1);
t473 = pkin(1) * t486;
t426 = t440 * t473;
t399 = pkin(7) * t428 + t426;
t439 = sin(qJ(3));
t442 = cos(qJ(3));
t559 = t399 + t464 * (pkin(3) * t439 - pkin(9) * t442);
t438 = sin(qJ(4));
t441 = cos(qJ(4));
t474 = t442 * t428;
t510 = qJD(1) * t437;
t495 = t440 * t510;
t369 = t438 * t474 - t441 * t495;
t505 = qJD(3) * t442;
t466 = t438 * t505 - t369;
t460 = t486 + qJD(2);
t381 = t439 * t495 - t442 * t460;
t499 = qJD(1) * qJD(2);
t489 = t437 * t499;
t471 = t443 * t489;
t461 = t442 * t471;
t445 = -qJD(3) * t381 + t461;
t558 = -qJD(4) * t464 + t445;
t434 = t437 ^ 2;
t557 = -0.2e1 * t434 * t499;
t556 = MDP(5) * (t440 ^ 2 - t443 ^ 2);
t383 = t439 * t460 + t442 * t495;
t472 = t440 * t489;
t504 = qJD(4) * t438;
t308 = t383 * t504 - t438 * t472 - t441 * t558;
t347 = t383 * t438 + t441 * t464;
t376 = qJD(4) + t381;
t538 = t347 * t376;
t555 = -t308 - t538;
t503 = qJD(4) * t441;
t309 = t383 * t503 + t438 * t558 - t441 * t472;
t349 = t441 * t383 - t438 * t464;
t537 = t349 * t376;
t554 = t309 + t537;
t496 = pkin(1) * t546;
t531 = t437 * t440;
t552 = -pkin(7) * t531 + t443 * t496;
t393 = -pkin(2) * t546 - t552;
t404 = t439 * t531 - t442 * t546;
t405 = t439 * t546 + t442 * t531;
t334 = t404 * pkin(3) - t405 * pkin(9) + t393;
t530 = t437 * t443;
t450 = pkin(7) * t530 + t440 * t496;
t394 = pkin(8) * t546 + t450;
t395 = (-pkin(2) * t443 - pkin(8) * t440 - pkin(1)) * t437;
t516 = t442 * t394 + t439 * t395;
t336 = -pkin(9) * t530 + t516;
t518 = t438 * t334 + t441 * t336;
t507 = qJD(3) * t439;
t498 = pkin(8) * t507;
t553 = t438 * t498 - t441 * t559;
t396 = -pkin(7) * t495 + t443 * t473;
t457 = (pkin(2) * t440 - pkin(8) * t443) * t437;
t397 = qJD(1) * t457;
t515 = t442 * t396 + t439 * t397;
t338 = pkin(9) * t495 + t515;
t421 = -pkin(3) * t442 - pkin(9) * t439 - pkin(2);
t551 = t441 * t338 - t421 * t503 + t438 * t559;
t550 = t349 ^ 2;
t444 = qJD(1) ^ 2;
t462 = t439 * t471;
t351 = qJD(3) * t383 + t462;
t549 = pkin(4) * t351;
t548 = t347 * pkin(4);
t547 = qJ(5) + pkin(9);
t545 = qJ(5) * t439;
t365 = pkin(8) * t460 + t399;
t375 = qJD(1) * t395;
t398 = qJD(2) * t457;
t390 = qJD(1) * t398;
t400 = t552 * qJD(2);
t391 = qJD(1) * t400;
t477 = -t365 * t505 - t375 * t507 + t442 * t390 - t439 * t391;
t300 = -pkin(3) * t472 - t477;
t286 = pkin(4) * t309 + t300;
t544 = t286 * t438;
t543 = t286 * t441;
t364 = -pkin(2) * t460 - t396;
t322 = t381 * pkin(3) - t383 * pkin(9) + t364;
t331 = t442 * t365 + t439 * t375;
t325 = -pkin(9) * t464 + t331;
t293 = t322 * t438 + t325 * t441;
t288 = -qJ(5) * t347 + t293;
t542 = t288 * t376;
t541 = t300 * t438;
t540 = t300 * t441;
t539 = t308 * t438;
t536 = t349 * t438;
t535 = t351 * t441;
t534 = t381 * t438;
t533 = t381 * t441;
t532 = t434 * t444;
t529 = t438 * t351;
t528 = t439 * t441;
t527 = t441 * t442;
t526 = t441 * t443;
t292 = t441 * t322 - t325 * t438;
t287 = -qJ(5) * t349 + t292;
t285 = pkin(4) * t376 + t287;
t525 = t285 - t287;
t330 = -t439 * t365 + t442 * t375;
t341 = pkin(3) * t383 + pkin(9) * t381;
t340 = t441 * t341;
t488 = qJD(4) * t547;
t524 = pkin(4) * t383 + qJ(5) * t533 + t441 * t488 + t340 + (qJD(5) - t330) * t438;
t502 = qJD(5) * t441;
t519 = t441 * t330 + t438 * t341;
t523 = qJ(5) * t534 + t438 * t488 - t502 + t519;
t370 = (t438 * t440 + t442 * t526) * t510;
t431 = pkin(8) * t527;
t475 = t439 * t428;
t522 = pkin(4) * t475 - qJ(5) * t370 - t338 * t438 + t439 * t502 - (pkin(4) * t439 - qJ(5) * t527) * qJD(3) - (-t431 + (-t421 + t545) * t438) * qJD(4) - t553;
t521 = -qJ(5) * t369 - (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t528 - (-qJD(5) * t439 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t442) * t438 + t551;
t384 = t439 * t396;
t337 = -pkin(3) * t495 - t397 * t442 + t384;
t520 = pkin(8) * t505 - t337 + (t439 * t503 + t466) * pkin(4);
t512 = t438 * t421 + t431;
t508 = qJD(2) * t440;
t506 = qJD(3) * t441;
t501 = t364 * qJD(3);
t500 = t437 * qJD(2);
t497 = t438 * t530;
t493 = t376 * t504;
t492 = t440 * t500;
t491 = t443 * t500;
t487 = qJD(5) + t548;
t485 = t441 * t334 - t336 * t438;
t484 = -t439 * t394 + t395 * t442;
t483 = t443 * t464;
t482 = t376 * t441;
t481 = qJD(3) * t464;
t479 = t434 * t440 * t443 * MDP(4);
t453 = t365 * t507 - t375 * t505 - t439 * t390 - t442 * t391;
t299 = pkin(9) * t472 - t453;
t463 = pkin(7) * t471;
t307 = t351 * pkin(3) - pkin(9) * t445 + qJD(2) * t426 + t463;
t478 = -t441 * t299 - t438 * t307 - t322 * t503 + t325 * t504;
t470 = MDP(15) * t495;
t468 = pkin(1) * t557;
t335 = pkin(3) * t530 - t484;
t465 = t441 * t505 - t370;
t458 = -t394 * t505 - t395 * t507 + t398 * t442 - t439 * t400;
t358 = t405 * t438 + t437 * t526;
t456 = -t376 * t503 - t529;
t324 = pkin(3) * t464 - t330;
t455 = -pkin(9) * t351 + t324 * t376;
t454 = qJ(5) * t309 + t478;
t452 = -t394 * t507 + t395 * t505 + t439 * t398 + t442 * t400;
t303 = pkin(9) * t492 + t452;
t356 = qJD(3) * t405 + t439 * t491;
t357 = -qJD(3) * t404 + t442 * t491;
t401 = t450 * qJD(2);
t314 = t356 * pkin(3) - t357 * pkin(9) + t401;
t451 = t441 * t303 + t438 * t314 + t334 * t503 - t336 * t504;
t449 = pkin(1) * (-qJD(2) * t486 + t532);
t304 = -pkin(3) * t492 - t458;
t284 = -qJD(4) * t293 - t299 * t438 + t441 * t307;
t448 = -t518 * qJD(4) - t303 * t438 + t441 * t314;
t446 = qJ(5) * t308 + t284;
t433 = -pkin(4) * t441 - pkin(3);
t423 = t547 * t441;
t422 = t547 * t438;
t414 = (pkin(4) * t438 + pkin(8)) * t439;
t411 = t441 * t421;
t392 = qJD(1) * t401;
t362 = -t438 * t545 + t512;
t359 = t405 * t441 - t497;
t352 = -qJ(5) * t528 + t411 + (-pkin(8) * t438 - pkin(4)) * t442;
t346 = t347 ^ 2;
t320 = -qJD(4) * t358 + t357 * t441 + t438 * t492;
t319 = -qJD(4) * t497 + t357 * t438 + t405 * t503 - t441 * t492;
t311 = -pkin(4) * t534 + t331;
t310 = pkin(4) * t358 + t335;
t301 = t324 + t487;
t294 = -qJ(5) * t358 + t518;
t290 = pkin(4) * t404 - qJ(5) * t359 + t485;
t289 = pkin(4) * t319 + t304;
t282 = -qJ(5) * t319 - qJD(5) * t358 + t451;
t281 = pkin(4) * t356 - qJ(5) * t320 - qJD(5) * t359 + t448;
t280 = -qJD(5) * t347 - t454;
t279 = -qJD(5) * t349 + t446 + t549;
t1 = [t556 * t557 + (-t392 * t546 - t401 * t460 + t440 * t468) * MDP(9) + (-t391 * t546 - t400 * t460 + t443 * t468) * MDP(10) + (t383 * t357 + t405 * t445) * MDP(11) + (-t405 * t351 - t383 * t356 - t357 * t381 - t404 * t445) * MDP(12) + (-t357 * t464 + t383 * t492 + t405 * t472 - t445 * t530) * MDP(13) + (t356 * t464 + (t351 * t443 + (-qJD(1) * t404 - t381) * t508) * t437) * MDP(14) + (-t434 * t509 - t437 * t464) * MDP(15) * t508 + (-t458 * t464 + t401 * t381 + t393 * t351 + t392 * t404 + t364 * t356 + (-t477 * t443 + (qJD(1) * t484 + t330) * t508) * t437) * MDP(16) + (-t331 * t492 + t364 * t357 + t401 * t383 + t392 * t405 + t393 * t445 + t452 * t464 - t453 * t530 - t472 * t516) * MDP(17) + (-t308 * t359 + t320 * t349) * MDP(18) + (t308 * t358 - t309 * t359 - t319 * t349 - t320 * t347) * MDP(19) + (-t308 * t404 + t320 * t376 + t349 * t356 + t351 * t359) * MDP(20) + (-t309 * t404 - t319 * t376 - t347 * t356 - t351 * t358) * MDP(21) + (t351 * t404 + t356 * t376) * MDP(22) + (t284 * t404 + t292 * t356 + t300 * t358 + t304 * t347 + t335 * t309 + t324 * t319 + t351 * t485 + t376 * t448) * MDP(23) + (-t293 * t356 + t300 * t359 + t304 * t349 - t335 * t308 + t324 * t320 - t351 * t518 - t376 * t451 + t404 * t478) * MDP(24) + (t279 * t404 + t281 * t376 + t285 * t356 + t286 * t358 + t289 * t347 + t290 * t351 + t301 * t319 + t309 * t310) * MDP(25) + (-t280 * t404 - t282 * t376 + t286 * t359 - t288 * t356 + t289 * t349 - t294 * t351 + t301 * t320 - t308 * t310) * MDP(26) + (-t279 * t359 - t280 * t358 - t281 * t349 - t282 * t347 - t285 * t320 - t288 * t319 + t290 * t308 - t294 * t309) * MDP(27) + (t279 * t290 + t280 * t294 + t281 * t285 + t282 * t288 + t286 * t310 + t289 * t301) * MDP(28) + 0.2e1 * t479 * t499 + (MDP(6) * t491 - MDP(7) * t492) * (0.2e1 * t486 + qJD(2)); t532 * t556 + (t399 * t460 + t440 * t449 - t463) * MDP(9) + (pkin(7) * t472 + t396 * t460 + t443 * t449) * MDP(10) + (-qJD(3) * t439 ^ 2 * t495 + ((qJD(3) * t460 + t471) * t439 - t464 * t383) * t442) * MDP(11) + (-t439 * t351 + t442 * t445 + (t475 - t507) * t383 + (t474 - t505) * t381) * MDP(12) + (-t442 * t481 + (t442 * t483 + (qJD(2) * t439 - t383) * t440) * t510) * MDP(13) + (t439 * t481 + (-t439 * t483 + (qJD(2) * t442 + t381) * t440) * t510) * MDP(14) + (-pkin(2) * t351 + t439 * t501 - t384 * t464 - t399 * t381 + (pkin(8) * t481 + t397 * t464 - t392) * t442 + (-t330 * t440 + (-pkin(8) * t508 - t364 * t443) * t439) * t510) * MDP(16) + (-pkin(2) * t445 + t331 * t495 - t364 * t474 - t399 * t383 + t392 * t439 + (-t498 - t515) * t464 + (-pkin(8) * t472 + t501) * t442) * MDP(17) + (-t308 * t528 + (-t439 * t504 + t465) * t349) * MDP(18) + (t347 * t370 + t349 * t369 + (-t347 * t441 - t536) * t505 + (t539 - t309 * t441 + (t347 * t438 - t349 * t441) * qJD(4)) * t439) * MDP(19) + (t308 * t442 + t465 * t376 + (-t349 * t464 - t493 + t535) * t439) * MDP(20) + (t309 * t442 - t466 * t376 + (t347 * t464 + t456) * t439) * MDP(21) + (-t376 * t439 * t464 - t351 * t442) * MDP(22) + (-t324 * t369 - t337 * t347 + t411 * t351 + ((-qJD(4) * t421 + t338) * t438 + t553) * t376 + (t324 * t438 * qJD(3) - t284 + (qJD(3) * t347 + t456) * pkin(8)) * t442 + (pkin(8) * t309 - t292 * t464 + t324 * t503 + t541) * t439) * MDP(23) + (-t512 * t351 - t337 * t349 - t324 * t370 + t551 * t376 + (t324 * t506 - t478 + (qJD(3) * t349 + t493) * pkin(8)) * t442 + (-t324 * t504 + t540 + t464 * t293 + (t376 * t506 - t308) * pkin(8)) * t439) * MDP(24) + (-t279 * t442 + t309 * t414 + t351 * t352 - t522 * t376 + t520 * t347 + t466 * t301 + (-t285 * t464 + t301 * t503 + t544) * t439) * MDP(25) + (t280 * t442 - t308 * t414 - t351 * t362 + t521 * t376 + t520 * t349 + t465 * t301 + (t288 * t464 - t301 * t504 + t543) * t439) * MDP(26) + (t285 * t370 + t288 * t369 + t308 * t352 - t309 * t362 + t522 * t349 + t521 * t347 + (-t285 * t441 - t288 * t438) * t505 + (-t279 * t441 - t280 * t438 + (t285 * t438 - t288 * t441) * qJD(4)) * t439) * MDP(27) + (t279 * t352 + t280 * t362 - t285 * t522 + t286 * t414 - t288 * t521 + t301 * t520) * MDP(28) + t464 * t470 + (-t479 + (-MDP(6) * t443 + MDP(7) * t440) * t437 * t546) * t444; -t381 ^ 2 * MDP(12) + (-t381 * t428 + t461) * MDP(13) - t462 * MDP(14) + qJD(2) * t470 + (-t331 * t464 + t477) * MDP(16) + (-t330 * t464 + t364 * t381 + t453) * MDP(17) + (t349 * t482 - t539) * MDP(18) + (-t554 * t438 + t555 * t441) * MDP(19) + (t376 * t482 + t529) * MDP(20) + (-t376 ^ 2 * t438 + t535) * MDP(21) + (-pkin(3) * t309 - t540 - t331 * t347 + (-pkin(9) * t503 - t340) * t376 + (t330 * t376 + t455) * t438) * MDP(23) + (pkin(3) * t308 + t541 - t331 * t349 + (pkin(9) * t504 + t519) * t376 + t455 * t441) * MDP(24) + (-t543 + t309 * t433 - t311 * t347 - t351 * t422 - t524 * t376 + (t301 * t381 + (t301 + t548) * qJD(4)) * t438) * MDP(25) + (t301 * t533 + t544 - t308 * t433 - t311 * t349 - t351 * t423 + t523 * t376 + (pkin(4) * t536 + t301 * t441) * qJD(4)) * MDP(26) + (-t308 * t422 - t309 * t423 + t524 * t349 + t523 * t347 + (-t285 * t376 + t280) * t441 + (-t279 - t542) * t438) * MDP(27) + (-t279 * t422 + t280 * t423 + t286 * t433 + (pkin(4) * t504 - t311) * t301 - t523 * t288 - t524 * t285) * MDP(28) + (MDP(11) * t381 + MDP(12) * t383 - MDP(14) * t428 - t364 * MDP(16) - t349 * MDP(20) + t347 * MDP(21) - t376 * MDP(22) - t292 * MDP(23) + t293 * MDP(24) - t285 * MDP(25) + t288 * MDP(26)) * t383; t349 * t347 * MDP(18) + (-t346 + t550) * MDP(19) + (-t308 + t538) * MDP(20) + (-t309 + t537) * MDP(21) + t351 * MDP(22) + (t293 * t376 - t324 * t349 + t284) * MDP(23) + (t292 * t376 + t324 * t347 + t478) * MDP(24) + (0.2e1 * t549 + t542 + (-t301 - t487) * t349 + t446) * MDP(25) + (-pkin(4) * t550 + t287 * t376 + (qJD(5) + t301) * t347 + t454) * MDP(26) + (pkin(4) * t308 - t347 * t525) * MDP(27) + (t525 * t288 + (-t301 * t349 + t279) * pkin(4)) * MDP(28); t554 * MDP(25) + t555 * MDP(26) + (-t346 - t550) * MDP(27) + (t285 * t349 + t288 * t347 + t286) * MDP(28);];
tauc = t1;
