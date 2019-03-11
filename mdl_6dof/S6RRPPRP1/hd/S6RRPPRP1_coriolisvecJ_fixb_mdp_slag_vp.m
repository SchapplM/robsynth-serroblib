% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:58
% EndTime: 2019-03-09 08:28:09
% DurationCPUTime: 7.23s
% Computational Cost: add. (6956->441), mult. (17736->571), div. (0->0), fcn. (13212->8), ass. (0->182)
t460 = sin(pkin(9));
t463 = sin(qJ(2));
t465 = cos(qJ(2));
t536 = cos(pkin(9));
t440 = t460 * t465 + t463 * t536;
t427 = t440 * qJD(1);
t459 = sin(pkin(10));
t461 = cos(pkin(10));
t400 = qJD(2) * t459 + t427 * t461;
t462 = sin(qJ(5));
t464 = cos(qJ(5));
t521 = t459 * t427;
t556 = qJD(2) * t461 - t521;
t477 = t464 * t556;
t352 = -t400 * t462 + t477;
t564 = t352 ^ 2;
t492 = t536 * t465;
t450 = qJD(1) * t492;
t505 = qJD(1) * t463;
t424 = t460 * t505 - t450;
t420 = qJD(5) + t424;
t563 = t352 * t420;
t558 = t464 * t400 + t462 * t556;
t542 = t558 ^ 2;
t562 = MDP(22) + MDP(24);
t561 = -MDP(23) + MDP(26);
t426 = t440 * qJD(2);
t416 = qJD(1) * t426;
t517 = t464 * t461;
t439 = t459 * t462 - t517;
t441 = t459 * t464 + t461 * t462;
t431 = t441 * qJD(5);
t507 = t441 * t424 + t431;
t487 = -t439 * t416 - t420 * t507;
t502 = qJD(5) * t464;
t503 = qJD(5) * t462;
t545 = -t459 * t503 + t461 * t502;
t508 = -t439 * t424 + t545;
t486 = t441 * t416 + t420 * t508;
t500 = qJD(1) * qJD(2);
t551 = -0.2e1 * t500;
t550 = MDP(5) * (t463 ^ 2 - t465 ^ 2);
t495 = t463 * t500;
t417 = qJD(2) * t450 - t460 * t495;
t549 = t441 * t417;
t475 = -t460 * t463 + t492;
t498 = -pkin(2) * t465 - pkin(1);
t385 = -pkin(3) * t475 - qJ(4) * t440 + t498;
t538 = -qJ(3) - pkin(7);
t447 = t538 * t463;
t448 = t538 * t465;
t398 = t460 * t447 - t448 * t536;
t338 = t461 * t385 - t398 * t459;
t540 = pkin(8) * t461;
t326 = -pkin(4) * t475 - t440 * t540 + t338;
t339 = t459 * t385 + t461 * t398;
t522 = t440 * t459;
t331 = -pkin(8) * t522 + t339;
t548 = t462 * t326 + t464 * t331;
t372 = pkin(2) * t505 + pkin(3) * t427 + qJ(4) * t424;
t444 = qJD(1) * t448;
t432 = t460 * t444;
t443 = qJD(1) * t447;
t393 = t443 * t536 + t432;
t334 = t461 * t372 - t393 * t459;
t317 = pkin(4) * t427 + t424 * t540 + t334;
t335 = t459 * t372 + t461 * t393;
t524 = t424 * t459;
t327 = pkin(8) * t524 + t335;
t453 = pkin(2) * t460 + qJ(4);
t539 = pkin(8) + t453;
t435 = t539 * t459;
t436 = t539 * t461;
t478 = -t435 * t464 - t436 * t462;
t547 = qJD(4) * t439 - qJD(5) * t478 + t462 * t317 + t464 * t327;
t383 = -t435 * t462 + t436 * t464;
t546 = -qJD(4) * t441 - qJD(5) * t383 - t317 * t464 + t327 * t462;
t525 = t417 * t459;
t319 = t462 * (qJD(5) * t400 + t525) - qJD(5) * t477 - t417 * t517;
t544 = t319 * t439 - t507 * t558;
t429 = t475 * qJD(2);
t537 = qJD(2) * pkin(2);
t499 = t463 * t537;
t354 = pkin(3) * t426 - qJ(4) * t429 - qJD(4) * t440 + t499;
t494 = qJD(2) * t538;
t421 = qJD(3) * t465 + t463 * t494;
t422 = -qJD(3) * t463 + t465 * t494;
t374 = t421 * t536 + t460 * t422;
t322 = t461 * t354 - t374 * t459;
t304 = pkin(4) * t426 - t429 * t540 + t322;
t323 = t459 * t354 + t461 * t374;
t523 = t429 * t459;
t312 = -pkin(8) * t523 + t323;
t543 = -qJD(5) * t548 + t304 * t464 - t312 * t462;
t423 = t424 ^ 2;
t541 = pkin(5) * t416;
t535 = qJ(6) * t416;
t485 = t498 * qJD(1);
t446 = qJD(3) + t485;
t364 = pkin(3) * t424 - qJ(4) * t427 + t446;
t437 = t443 + t537;
t493 = t536 * t444;
t387 = t460 * t437 - t493;
t381 = qJD(2) * qJ(4) + t387;
t328 = t461 * t364 - t381 * t459;
t307 = pkin(4) * t424 - pkin(8) * t400 + t328;
t329 = t459 * t364 + t461 * t381;
t318 = pkin(8) * t556 + t329;
t289 = t307 * t462 + t318 * t464;
t534 = t289 * t420;
t532 = t352 * t427;
t531 = t558 * t352;
t530 = t558 * t427;
t409 = t421 * qJD(1);
t410 = t422 * qJD(1);
t360 = t460 * t409 - t536 * t410;
t397 = -t536 * t447 - t448 * t460;
t528 = t360 * t397;
t527 = t478 * t416;
t526 = t383 * t416;
t519 = t461 * t417;
t466 = qJD(2) ^ 2;
t518 = t463 * t466;
t516 = t465 * t466;
t467 = qJD(1) ^ 2;
t515 = t465 * t467;
t514 = qJ(6) * t427 + t547;
t513 = -pkin(5) * t427 + t546;
t392 = t443 * t460 - t493;
t358 = -pkin(4) * t524 + t392;
t512 = -t507 * pkin(5) + t508 * qJ(6) + qJD(6) * t441 + t358;
t452 = pkin(2) * t495;
t342 = pkin(3) * t416 - qJ(4) * t417 - qJD(4) * t427 + t452;
t361 = t536 * t409 + t460 * t410;
t357 = qJD(2) * qJD(4) + t361;
t310 = t459 * t342 + t461 * t357;
t288 = t307 * t464 - t318 * t462;
t501 = qJD(6) - t288;
t337 = pkin(4) * t525 + t360;
t491 = pkin(1) * t551;
t309 = t461 * t342 - t357 * t459;
t373 = t421 * t460 - t536 * t422;
t298 = pkin(4) * t416 - pkin(8) * t519 + t309;
t301 = -pkin(8) * t525 + t310;
t489 = -t464 * t298 + t462 * t301 + t307 * t503 + t318 * t502;
t456 = -pkin(2) * t536 - pkin(3);
t320 = qJD(5) * t558 + t549;
t488 = -t441 * t320 + t352 * t508;
t344 = pkin(4) * t523 + t373;
t370 = pkin(4) * t522 + t397;
t386 = t437 * t536 + t432;
t482 = t326 * t464 - t331 * t462;
t480 = -t328 * t459 + t329 * t461;
t479 = t360 * t440 + t397 * t417;
t376 = -qJD(2) * pkin(3) + qJD(4) - t386;
t343 = -pkin(4) * t556 + t376;
t295 = -pkin(5) * t352 - qJ(6) * t558 + t343;
t474 = t295 * t558 + t489;
t445 = -t461 * pkin(4) + t456;
t473 = t462 * t298 + t464 * t301 + t307 * t502 - t318 * t503;
t472 = t462 * t304 + t464 * t312 + t326 * t502 - t331 * t503;
t471 = t376 * t429 + t479;
t469 = -t453 * t416 + t456 * t417 + (-qJD(4) + t376) * t424;
t284 = pkin(5) * t320 + qJ(6) * t319 - qJD(6) * t558 + t337;
t378 = t439 * t440;
t377 = t441 * t440;
t375 = t439 * pkin(5) - t441 * qJ(6) + t445;
t333 = t441 * t429 + t440 * t545;
t332 = t429 * t439 + t440 * t431;
t316 = pkin(5) * t558 - qJ(6) * t352;
t311 = pkin(5) * t377 + qJ(6) * t378 + t370;
t294 = -t319 - t563;
t293 = pkin(5) * t475 - t482;
t292 = -qJ(6) * t475 + t548;
t287 = pkin(5) * t333 + qJ(6) * t332 + qJD(6) * t378 + t344;
t286 = qJ(6) * t420 + t289;
t285 = -pkin(5) * t420 + t501;
t283 = -pkin(5) * t426 - t543;
t282 = qJ(6) * t426 - qJD(6) * t475 + t472;
t281 = t489 - t541;
t280 = qJD(6) * t420 + t473 + t535;
t1 = [0.2e1 * t465 * MDP(4) * t495 + t550 * t551 + MDP(6) * t516 - MDP(7) * t518 + (-pkin(7) * t516 + t463 * t491) * MDP(9) + (pkin(7) * t518 + t465 * t491) * MDP(10) + (t361 * t475 + t373 * t427 - t374 * t424 - t386 * t429 - t387 * t426 - t398 * t416 + t479) * MDP(11) + (t528 + t361 * t398 - t386 * t373 + t387 * t374 + (t446 + t485) * t499) * MDP(12) + (-t309 * t475 + t322 * t424 + t328 * t426 + t338 * t416 - t373 * t556 + t459 * t471) * MDP(13) + (t310 * t475 - t323 * t424 - t329 * t426 - t339 * t416 + t373 * t400 + t461 * t471) * MDP(14) + (-t322 * t400 - t323 * t521 + (qJD(2) * t323 - t309 * t440 - t328 * t429 - t338 * t417) * t461 + (-t310 * t440 - t329 * t429 - t339 * t417) * t459) * MDP(15) + (t309 * t338 + t310 * t339 + t322 * t328 + t323 * t329 + t373 * t376 + t528) * MDP(16) + (t319 * t378 - t332 * t558) * MDP(17) + (t319 * t377 + t320 * t378 - t332 * t352 - t333 * t558) * MDP(18) + (t319 * t475 - t332 * t420 - t378 * t416 + t426 * t558) * MDP(19) + (t320 * t475 - t333 * t420 + t352 * t426 - t377 * t416) * MDP(20) + (-t416 * t475 + t420 * t426) * MDP(21) + (t288 * t426 + t370 * t320 + t343 * t333 + t337 * t377 - t344 * t352 + t482 * t416 + t420 * t543 + t475 * t489) * MDP(22) + (-t289 * t426 - t370 * t319 - t343 * t332 - t337 * t378 + t344 * t558 - t416 * t548 - t420 * t472 + t473 * t475) * MDP(23) + (t281 * t475 - t283 * t420 + t284 * t377 - t285 * t426 - t287 * t352 - t293 * t416 + t295 * t333 + t311 * t320) * MDP(24) + (-t280 * t377 - t281 * t378 + t282 * t352 + t283 * t558 - t285 * t332 - t286 * t333 - t292 * t320 - t293 * t319) * MDP(25) + (-t280 * t475 + t282 * t420 + t284 * t378 + t286 * t426 - t287 * t558 + t292 * t416 + t295 * t332 + t311 * t319) * MDP(26) + (t280 * t292 + t281 * t293 + t282 * t286 + t283 * t285 + t284 * t311 + t287 * t295) * MDP(27); -t463 * MDP(4) * t515 + t467 * t550 + ((t387 - t392) * t427 + (-t386 + t393) * t424 + (-t416 * t460 - t417 * t536) * pkin(2)) * MDP(11) + (t386 * t392 - t387 * t393 + (-t360 * t536 + t361 * t460 - t446 * t505) * pkin(2)) * MDP(12) + (-t328 * t427 - t334 * t424 - t360 * t461 + t392 * t556 + t459 * t469) * MDP(13) + (t329 * t427 + t335 * t424 + t360 * t459 - t392 * t400 + t461 * t469) * MDP(14) + (t334 * t400 + t335 * t521 + (-qJD(4) * t521 - t328 * t424 + t310 + (qJD(4) * t461 - t335) * qJD(2)) * t461 + (qJD(4) * t400 - t329 * t424 - t309) * t459) * MDP(15) + (-t328 * t334 - t329 * t335 + t360 * t456 - t376 * t392 + (-t309 * t459 + t310 * t461) * t453 + t480 * qJD(4)) * MDP(16) + (-t319 * t441 + t508 * t558) * MDP(17) + (t488 + t544) * MDP(18) + (t486 - t530) * MDP(19) + (t487 - t532) * MDP(20) - t420 * t427 * MDP(21) + (-t288 * t427 + t445 * t320 + t337 * t439 + t507 * t343 + t352 * t358 + t420 * t546 + t527) * MDP(22) + (t289 * t427 - t445 * t319 + t337 * t441 + t508 * t343 - t358 * t558 + t420 * t547 - t526) * MDP(23) + (t284 * t439 + t285 * t427 + t295 * t507 + t320 * t375 + t352 * t512 + t420 * t513 + t527) * MDP(24) + (-t280 * t439 + t281 * t441 + t285 * t508 - t286 * t507 + t319 * t478 - t320 * t383 - t352 * t514 - t513 * t558) * MDP(25) + (-t284 * t441 - t286 * t427 - t295 * t508 + t319 * t375 - t420 * t514 + t512 * t558 + t526) * MDP(26) + (t280 * t383 - t281 * t478 + t284 * t375 - t285 * t513 - t286 * t514 - t295 * t512) * MDP(27) + (MDP(9) * t463 * t467 + MDP(10) * t515) * pkin(1); (-t427 ^ 2 - t423) * MDP(11) + (t386 * t427 + t387 * t424 + t452) * MDP(12) + (t461 * t416 - t423 * t459 + t427 * t556) * MDP(13) + (-t400 * t427 - t416 * t459 - t423 * t461) * MDP(14) + ((t459 * t400 + t461 * t556) * t424 + (-t459 ^ 2 - t461 ^ 2) * t417) * MDP(15) + (t309 * t461 + t310 * t459 - t376 * t427 + t424 * t480) * MDP(16) + (t488 - t544) * MDP(25) + (t280 * t441 + t281 * t439 + t285 * t507 + t286 * t508 - t295 * t427) * MDP(27) + t561 * (t486 + t530) + t562 * (t487 + t532); (t400 * t424 + t525) * MDP(13) + (t424 * t556 + t519) * MDP(14) + (-t400 ^ 2 - t556 ^ 2) * MDP(15) + (t328 * t400 - t329 * t556 + t360) * MDP(16) + (-t542 - t564) * MDP(25) + (-t285 * t558 - t286 * t352 + t284) * MDP(27) + t562 * (t420 * t558 + t320) + t561 * (t319 - t563); -MDP(17) * t531 + (t542 - t564) * MDP(18) + t294 * MDP(19) + (-t549 + (-qJD(5) + t420) * t558) * MDP(20) + t416 * MDP(21) + (-t343 * t558 - t489 + t534) * MDP(22) + (t288 * t420 - t343 * t352 - t473) * MDP(23) + (t316 * t352 - t474 + t534 + 0.2e1 * t541) * MDP(24) + (pkin(5) * t319 - qJ(6) * t320 + (t286 - t289) * t558 - (t285 - t501) * t352) * MDP(25) + (0.2e1 * t535 + t295 * t352 + t316 * t558 + (0.2e1 * qJD(6) - t288) * t420 + t473) * MDP(26) + (-pkin(5) * t281 + qJ(6) * t280 - t285 * t289 + t286 * t501 - t295 * t316) * MDP(27); (-qJD(2) * t427 - t531) * MDP(24) + t294 * MDP(25) + (-t420 ^ 2 - t542) * MDP(26) + (-t286 * t420 + t474 - t541) * MDP(27);];
tauc  = t1;
