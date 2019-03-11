% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:06
% EndTime: 2019-03-09 08:44:14
% DurationCPUTime: 4.67s
% Computational Cost: add. (4243->438), mult. (9598->564), div. (0->0), fcn. (6040->6), ass. (0->191)
t443 = cos(pkin(9));
t447 = cos(qJ(2));
t512 = qJD(1) * t447;
t492 = t443 * t512;
t442 = sin(pkin(9));
t504 = t442 * qJD(2);
t390 = -t492 - t504;
t493 = t442 * t512;
t509 = qJD(2) * t443;
t391 = -t493 + t509;
t445 = sin(qJ(5));
t538 = cos(qJ(5));
t342 = -t538 * t390 + t391 * t445;
t558 = t342 ^ 2;
t463 = -t445 * t390 - t391 * t538;
t540 = t463 ^ 2;
t539 = pkin(3) + pkin(7);
t446 = sin(qJ(2));
t513 = qJD(1) * t446;
t426 = qJD(5) + t513;
t557 = t342 * t426;
t556 = t426 * t463;
t460 = -t445 * t442 + t538 * t443;
t371 = t460 * t447;
t497 = t446 * t538;
t481 = qJD(1) * t497;
t495 = t445 * t513;
t363 = -t442 * t481 - t443 * t495;
t491 = qJD(5) * t538;
t505 = qJD(5) * t445;
t382 = -t442 * t491 - t443 * t505;
t516 = t382 + t363;
t515 = (t481 + t491) * t443 + (-t495 - t505) * t442;
t429 = pkin(7) * t513;
t555 = qJD(3) + t429;
t500 = qJD(1) * qJD(2);
t490 = t446 * t500;
t425 = pkin(2) * t490;
t475 = -qJ(3) * t447 + qJ(4) * t446;
t506 = qJD(3) * t446;
t452 = qJD(2) * t475 - qJD(4) * t447 - t506;
t340 = qJD(1) * t452 + t425;
t489 = t447 * t500;
t424 = pkin(7) * t489;
t431 = pkin(3) * t512;
t375 = t424 + (-qJD(4) + t431) * qJD(2);
t311 = -t340 * t442 + t443 * t375;
t529 = t442 * t446;
t469 = pkin(4) * t447 - pkin(8) * t529;
t454 = t469 * qJD(2);
t298 = qJD(1) * t454 + t311;
t312 = t443 * t340 + t442 * t375;
t528 = t443 * t446;
t499 = pkin(8) * t528;
t484 = qJD(2) * t499;
t304 = qJD(1) * t484 + t312;
t444 = -pkin(2) - qJ(4);
t487 = -qJ(3) * t446 - pkin(1);
t389 = t444 * t447 + t487;
t364 = t389 * qJD(1);
t502 = pkin(3) * t513 + t555;
t369 = qJD(2) * t444 + t502;
t324 = -t364 * t442 + t443 * t369;
t305 = pkin(4) * t513 - pkin(8) * t391 + t324;
t325 = t443 * t364 + t442 * t369;
t309 = pkin(8) * t390 + t325;
t457 = -t445 * t298 - t538 * t304 - t305 * t491 + t309 * t505;
t479 = qJ(6) * t489;
t280 = qJD(6) * t426 - t457 + t479;
t482 = -t538 * t298 + t445 * t304 + t305 * t505 + t309 * t491;
t483 = pkin(5) * t489;
t281 = t482 - t483;
t287 = t305 * t538 - t445 * t309;
t501 = qJD(6) - t287;
t285 = -t426 * pkin(5) + t501;
t288 = t445 * t305 + t309 * t538;
t286 = t426 * qJ(6) + t288;
t461 = -t442 * t538 - t445 * t443;
t554 = -t280 * t461 - t281 * t460 - t285 * t516 + t286 * t515;
t553 = -0.2e1 * t500;
t440 = t446 ^ 2;
t552 = MDP(5) * (-t447 ^ 2 + t440);
t549 = t446 * MDP(4);
t433 = pkin(2) * t513;
t376 = qJD(1) * t475 + t433;
t430 = pkin(7) * t512;
t400 = t430 + t431;
t337 = -t376 * t442 + t443 * t400;
t320 = qJD(1) * t469 + t337;
t338 = t443 * t376 + t442 * t400;
t329 = qJD(1) * t499 + t338;
t537 = -pkin(8) + t444;
t402 = t537 * t442;
t403 = t537 * t443;
t462 = -t445 * t402 + t403 * t538;
t548 = -qJD(4) * t461 - qJD(5) * t462 + t445 * t320 + t538 * t329;
t352 = t402 * t538 + t445 * t403;
t547 = -qJD(4) * t460 - qJD(5) * t352 - t320 * t538 + t445 * t329;
t415 = t539 * t446;
t396 = t443 * t415;
t331 = pkin(4) * t446 + t396 + (pkin(8) * t447 - t389) * t442;
t348 = t443 * t389 + t442 * t415;
t527 = t443 * t447;
t335 = -pkin(8) * t527 + t348;
t546 = t445 * t331 + t538 * t335;
t498 = -pkin(4) * t443 - pkin(3);
t503 = -t498 * t513 + t555;
t439 = qJD(2) * qJ(3);
t545 = qJD(4) + t439;
t544 = -MDP(25) + MDP(28);
t480 = qJD(2) * t497;
t470 = qJD(1) * t480;
t478 = t445 * t490;
t313 = -t390 * t491 + t391 * t505 - t442 * t470 - t443 * t478;
t543 = -t313 * t460 - t463 * t516;
t380 = t400 + t545;
t507 = qJD(2) * t447;
t542 = t446 * (-t380 + t545) - t444 * t507;
t314 = -qJD(5) * t463 + t442 * t478 - t443 * t470;
t508 = qJD(2) * t446;
t432 = pkin(2) * t508;
t354 = t432 + t452;
t401 = t539 * t507;
t327 = -t354 * t442 + t443 * t401;
t310 = t327 + t454;
t328 = t443 * t354 + t442 * t401;
t317 = t484 + t328;
t541 = -qJD(5) * t546 + t310 * t538 - t445 * t317;
t536 = qJD(2) * pkin(2);
t349 = -pkin(4) * t390 + t380;
t292 = pkin(5) * t342 + qJ(6) * t463 + t349;
t535 = t292 * t463;
t533 = t324 * t447;
t532 = t325 * t447;
t531 = t342 * t463;
t399 = t539 * t508;
t438 = qJD(2) * qJD(3);
t374 = -qJD(1) * t399 + t438;
t530 = t374 * t442;
t448 = qJD(2) ^ 2;
t525 = t446 * t448;
t524 = t447 * t448;
t449 = qJD(1) ^ 2;
t523 = t447 * t449;
t427 = t442 * pkin(4) + qJ(3);
t522 = qJ(6) * t512 + t548;
t521 = -pkin(5) * t512 + t547;
t520 = -t515 * pkin(5) + t516 * qJ(6) + qJD(6) * t460 - t503;
t517 = t382 * t426 + t460 * t489;
t416 = t539 * t447;
t409 = -pkin(2) * t447 + t487;
t388 = qJD(1) * t409;
t511 = qJD(2) * t462;
t510 = qJD(2) * t352;
t381 = pkin(4) * t527 + t416;
t488 = MDP(23) * t512;
t486 = pkin(1) * t553;
t485 = qJD(3) - t536;
t473 = t311 * t443 + t312 * t442;
t472 = -t324 * t442 + t325 * t443;
t471 = -0.2e1 * qJD(2) * t388;
t459 = -qJ(3) * t507 - t506;
t365 = qJD(1) * t459 + t425;
t378 = t432 + t459;
t468 = pkin(7) * t448 + qJD(1) * t378 + t365;
t466 = t331 * t538 - t445 * t335;
t367 = (-pkin(7) + t498) * t508;
t458 = t288 * t426 - t482;
t456 = t445 * t310 + t538 * t317 + t331 * t491 - t335 * t505;
t372 = t461 * t447;
t353 = qJD(1) * t367 + t438;
t407 = pkin(7) * t490 - t438;
t408 = t429 + t485;
t414 = -t430 - t439;
t451 = -t407 * t447 + (t408 * t447 + (t414 + t430) * t446) * qJD(2);
t284 = pkin(5) * t314 + qJ(6) * t313 + qJD(6) * t463 + t353;
t397 = -qJ(3) * t512 + t433;
t368 = t388 * t513;
t347 = -t389 * t442 + t396;
t339 = -pkin(5) * t461 - qJ(6) * t460 + t427;
t334 = -t445 * t446 * t504 - qJD(5) * t372 + t443 * t480;
t333 = -qJD(5) * t371 - t461 * t508;
t316 = pkin(5) * t371 - qJ(6) * t372 + t381;
t303 = -pkin(5) * t463 + qJ(6) * t342;
t295 = -t446 * pkin(5) - t466;
t294 = qJ(6) * t446 + t546;
t293 = -t313 + t557;
t289 = -pkin(5) * t334 - qJ(6) * t333 - qJD(6) * t372 + t367;
t283 = -pkin(5) * t507 - t541;
t282 = qJ(6) * t507 + qJD(6) * t446 + t456;
t1 = [(-t313 * t372 - t333 * t463) * MDP(19) + (-t313 * t446 + t333 * t426 + (qJD(1) * t372 - t463) * t507) * MDP(21) + (t280 * t446 + t282 * t426 - t284 * t372 + t289 * t463 - t292 * t333 + t313 * t316 + (qJD(1) * t294 + t286) * t507) * MDP(28) + (-t456 * t426 + t457 * t446 - t367 * t463 - t381 * t313 + t353 * t372 + t349 * t333 + (-qJD(1) * t546 - t288) * t507) * MDP(25) + (t287 * t507 + t381 * t314 - t349 * t334 + t367 * t342 + t353 * t371 + t426 * t541 - t482 * t446 + t466 * t489) * MDP(24) + (-t280 * t371 + t281 * t372 - t282 * t342 - t283 * t463 + t285 * t333 + t286 * t334 - t294 * t314 - t295 * t313) * MDP(27) + (t313 * t371 - t314 * t372 - t333 * t342 - t334 * t463) * MDP(20) + (-t314 * t446 + t334 * t426 + (-qJD(1) * t371 - t342) * t507) * MDP(22) + (-t281 * t446 - t283 * t426 + t284 * t371 + t289 * t342 - t292 * t334 + t314 * t316 + (-qJD(1) * t295 - t285) * t507) * MDP(26) + t451 * MDP(11) + (pkin(7) * t451 + t365 * t409 + t378 * t388) * MDP(14) + (t280 * t294 + t281 * t295 + t282 * t286 + t283 * t285 + t284 * t316 + t289 * t292) * MDP(29) + 0.2e1 * t489 * t549 + (t374 * t527 + t390 * t399 + (qJD(1) * t327 + t311) * t446 + (-t380 * t528 + t533 + (t347 * t447 - t416 * t528) * qJD(1)) * qJD(2)) * MDP(15) + (t311 * t347 + t312 * t348 + t324 * t327 + t325 * t328 + t374 * t416 - t380 * t399) * MDP(18) + (-t447 * t530 - t391 * t399 + (-qJD(1) * t328 - t312) * t446 + (t380 * t529 - t532 + (-t348 * t447 + t416 * t529) * qJD(1)) * qJD(2)) * MDP(16) + (-pkin(7) * t524 + t446 * t486) * MDP(9) - MDP(7) * t525 + (pkin(7) * t525 + t447 * t486) * MDP(10) + (-t327 * t391 + t328 * t390 + (t311 * t442 - t312 * t443) * t447 + ((-t347 * t442 + t348 * t443) * qJD(1) + t472) * t508) * MDP(17) + t552 * t553 + (t426 + t513) * MDP(23) * t507 + (-t446 * t468 + t447 * t471) * MDP(13) + (t446 * t471 + t447 * t468) * MDP(12) + MDP(6) * t524; -t523 * t549 + t449 * t552 + ((-t414 - t439) * t446 + (-t408 + t485) * t447) * qJD(1) * MDP(11) + (-t397 * t512 + t368) * MDP(12) + (0.2e1 * t438 + (t388 * t447 + t397 * t446) * qJD(1)) * MDP(13) + (-qJ(3) * t407 - qJD(3) * t414 - t388 * t397 + (-t414 * t446 + (-t408 - t536) * t447) * qJD(1) * pkin(7)) * MDP(14) + (t530 - t502 * t390 + (-t337 * t446 - t443 * t542 - t533) * qJD(1)) * MDP(15) + (t374 * t443 + t502 * t391 + (t338 * t446 + t442 * t542 + t532) * qJD(1)) * MDP(16) + (t337 * t391 - t338 * t390 + (qJD(4) * t391 - t325 * t513 - t311) * t443 + (-qJD(4) * t390 + t324 * t513 - t312) * t442) * MDP(17) + (qJ(3) * t374 - t324 * t337 - t325 * t338 + t473 * t444 + t502 * t380 + (-t324 * t443 - t325 * t442) * qJD(4)) * MDP(18) + t543 * MDP(19) + (-t313 * t461 - t314 * t460 - t342 * t516 + t463 * t515) * MDP(20) + (t363 * t426 + t463 * t512 + t517) * MDP(21) + (-t515 * t426 + (qJD(2) * t461 + t342) * t512) * MDP(22) - t426 * t488 + (t427 * t314 - t353 * t461 + t547 * t426 + t515 * t349 + t503 * t342 + (-t287 + t511) * t512) * MDP(24) + (-t427 * t313 + t353 * t460 + t548 * t426 + t516 * t349 - t503 * t463 + (t288 - t510) * t512) * MDP(25) + (-t284 * t461 + t314 * t339 + t521 * t426 - t520 * t342 + t515 * t292 + (t285 + t511) * t512) * MDP(26) + (t313 * t462 - t314 * t352 + t342 * t522 + t463 * t521 - t554) * MDP(27) + (-t284 * t460 + t313 * t339 - t522 * t426 - t520 * t463 - t516 * t292 + (-t286 + t510) * t512) * MDP(28) + (t280 * t352 - t281 * t462 + t284 * t339 - t285 * t521 - t286 * t522 - t292 * t520) * MDP(29) + (MDP(9) * t446 * t449 + MDP(10) * t523) * pkin(1); -t448 * MDP(13) + (t368 + t424) * MDP(14) + t473 * MDP(18) + t517 * MDP(24) + (t314 * t461 - t342 * t515 - t543) * MDP(27) + t554 * MDP(29) + (-MDP(15) * t442 - MDP(16) * t443 - MDP(13)) * t449 * t440 + (MDP(12) * t523 + ((t390 * t443 + t391 * t442) * MDP(17) + t472 * MDP(18)) * qJD(1)) * t446 + (t363 * MDP(24) + MDP(26) * t516 + t515 * t544) * t426 + (t414 * MDP(14) + (t390 + t492) * MDP(15) + (-t391 - t493) * MDP(16) - t380 * MDP(18) - t342 * MDP(24) + (t460 * t512 - t342) * MDP(26) - t292 * MDP(29) + t544 * (-t461 * t512 - t463)) * qJD(2); (-t390 ^ 2 - t391 ^ 2) * MDP(17) + (t324 * t391 - t325 * t390 + t374) * MDP(18) + (-t540 - t558) * MDP(27) + (t285 * t463 + t286 * t342 + t284) * MDP(29) + ((t391 - t509) * MDP(15) + (t390 + t504) * MDP(16)) * t513 + t544 * (t313 + t557) + (MDP(24) + MDP(26)) * (t314 - t556); -MDP(19) * t531 + (t540 - t558) * MDP(20) + t293 * MDP(21) + (-t314 - t556) * MDP(22) + qJD(2) * t488 + (t349 * t463 + t458) * MDP(24) + (t287 * t426 + t342 * t349 + t457) * MDP(25) + (-t303 * t342 + t458 + 0.2e1 * t483 + t535) * MDP(26) + (pkin(5) * t313 - qJ(6) * t314 - (t286 - t288) * t463 + (t285 - t501) * t342) * MDP(27) + (0.2e1 * t479 - t292 * t342 - t303 * t463 + (0.2e1 * qJD(6) - t287) * t426 - t457) * MDP(28) + (-pkin(5) * t281 + qJ(6) * t280 - t285 * t288 + t286 * t501 - t292 * t303) * MDP(29); (-t489 - t531) * MDP(26) + t293 * MDP(27) + (-t426 ^ 2 - t540) * MDP(28) + (-t286 * t426 + t281 - t535) * MDP(29);];
tauc  = t1;
