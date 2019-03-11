% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:37
% EndTime: 2019-03-09 09:52:46
% DurationCPUTime: 4.48s
% Computational Cost: add. (5342->471), mult. (13299->573), div. (0->0), fcn. (9289->6), ass. (0->192)
t424 = sin(pkin(9));
t426 = sin(qJ(2));
t428 = cos(qJ(2));
t518 = cos(pkin(9));
t407 = t424 * t428 + t426 * t518;
t394 = t407 * qJD(1);
t425 = sin(qJ(4));
t427 = cos(qJ(4));
t362 = qJD(2) * t425 + t394 * t427;
t393 = t407 * qJD(2);
t382 = qJD(1) * t393;
t377 = t382 * pkin(4);
t520 = -qJ(3) - pkin(7);
t468 = qJD(2) * t520;
t389 = qJD(3) * t428 + t426 * t468;
t374 = t389 * qJD(1);
t390 = -t426 * qJD(3) + t428 * t468;
t435 = t390 * qJD(1);
t320 = t374 * t518 + t424 * t435;
t479 = qJD(1) * qJD(2);
t469 = t426 * t479;
t416 = pkin(2) * t469;
t466 = t518 * t428;
t415 = qJD(1) * t466;
t440 = qJD(2) * t415 - t424 * t469;
t325 = t382 * pkin(3) - pkin(8) * t440 + t416;
t484 = t426 * qJD(1);
t392 = -t424 * t484 + t415;
t474 = -pkin(2) * t428 - pkin(1);
t455 = t474 * qJD(1);
t411 = qJD(3) + t455;
t326 = -pkin(3) * t392 - pkin(8) * t394 + t411;
t412 = t520 * t426;
t409 = qJD(1) * t412;
t519 = qJD(2) * pkin(2);
t401 = t409 + t519;
t413 = t520 * t428;
t410 = qJD(1) * t413;
t467 = t518 * t410;
t348 = t424 * t401 - t467;
t342 = qJD(2) * pkin(8) + t348;
t487 = qJD(4) * t427;
t488 = qJD(4) * t425;
t457 = t425 * t320 - t427 * t325 + t326 * t488 + t342 * t487;
t278 = -t377 + t457;
t483 = t427 * qJD(2);
t450 = qJD(4) * t483 - t394 * t488 + t427 * t440;
t439 = -qJ(6) * t450 + t278;
t434 = -pkin(5) * t382 + t439;
t272 = -qJD(6) * t362 + t434;
t300 = t425 * t326 + t427 * t342;
t360 = t394 * t425 - t483;
t290 = qJ(6) * t360 + t300;
t388 = qJD(4) - t392;
t380 = t388 * qJ(5);
t285 = t290 + t380;
t514 = t285 * t388;
t539 = -t272 + t514;
t445 = -t427 * t320 - t425 * t325 - t326 * t487 + t342 * t488;
t537 = t382 * qJ(5) + t388 * qJD(5);
t276 = -t445 + t537;
t322 = qJD(4) * t362 + t425 * t440;
t530 = t322 * qJ(6) + t360 * qJD(6);
t273 = t276 + t530;
t299 = t427 * t326 - t425 * t342;
t289 = qJ(6) * t362 + t299;
t482 = qJD(5) - t289;
t522 = pkin(4) + pkin(5);
t283 = -t388 * t522 + t482;
t538 = t283 * t388 + t273;
t536 = MDP(22) + MDP(25);
t535 = -0.2e1 * t479;
t534 = MDP(5) * (t426 ^ 2 - t428 ^ 2);
t532 = t427 * t522;
t397 = t424 * t410;
t347 = t518 * t401 + t397;
t456 = qJD(2) * pkin(3) + t347;
t441 = qJ(5) * t362 + t456;
t288 = -t360 * t522 + qJD(6) + t441;
t531 = (qJD(6) + t288) * t362;
t349 = t409 * t424 - t467;
t529 = -qJD(5) * t425 - t349;
t528 = 0.2e1 * t537;
t446 = -t424 * t426 + t466;
t526 = t393 * qJ(5) - qJD(5) * t446;
t481 = qJD(5) - t299;
t478 = MDP(20) + MDP(24);
t525 = MDP(21) - MDP(26);
t301 = pkin(4) * t360 - t441;
t418 = pkin(2) * t424 + pkin(8);
t504 = t418 * t382;
t524 = t301 * t388 - t504;
t523 = t360 ^ 2;
t359 = t362 ^ 2;
t521 = t427 * pkin(4);
t517 = qJ(5) * t322;
t516 = qJ(5) * t360;
t515 = qJ(5) * t427;
t513 = t300 * t388;
t512 = t450 * t425;
t511 = t360 * t388;
t510 = t360 * t392;
t509 = t362 * t360;
t508 = t362 * t388;
t458 = t388 * t427;
t507 = t392 * t425;
t506 = t407 * t425;
t505 = t407 * t427;
t503 = t425 * qJ(5);
t370 = t425 * t382;
t430 = qJD(2) ^ 2;
t502 = t426 * t430;
t371 = t427 * t382;
t501 = t428 * t430;
t431 = qJD(1) ^ 2;
t500 = t428 * t431;
t499 = qJ(6) - t418;
t336 = pkin(2) * t484 + pkin(3) * t394 - pkin(8) * t392;
t350 = t409 * t518 + t397;
t344 = t425 * t350;
t405 = t499 * t427;
t498 = t344 + (-qJ(6) * t392 - t336) * t427 - t522 * t394 + qJD(4) * t405 + qJD(6) * t425;
t492 = t425 * t336 + t427 * t350;
t295 = t394 * qJ(5) + t492;
t497 = -qJ(6) * t507 - qJD(6) * t427 + t488 * t499 - t295;
t449 = -t425 * t522 + t515;
t496 = -t388 * t449 + t529;
t454 = pkin(4) * t425 - t515;
t495 = -t388 * t454 - t529;
t494 = -t425 * t322 - t360 * t487;
t493 = t388 * t507 + t371;
t346 = -pkin(3) * t446 - pkin(8) * t407 + t474;
t356 = t424 * t412 - t413 * t518;
t491 = t425 * t346 + t427 * t356;
t490 = t388 * t487 + t370;
t485 = qJD(5) * t427;
t477 = t426 * t519;
t333 = t389 * t518 + t424 * t390;
t476 = t425 * t333 + t346 * t488 + t356 * t487;
t396 = t446 * qJD(2);
t337 = pkin(3) * t393 - pkin(8) * t396 + t477;
t475 = t427 * t333 + t425 * t337 + t346 * t487;
t302 = -qJ(5) * t446 + t491;
t473 = t407 * t488;
t472 = t407 * t487;
t471 = t418 * t488;
t470 = t388 * t488;
t465 = MDP(19) - t536;
t464 = pkin(1) * t535;
t292 = -pkin(4) * t388 + t481;
t463 = -t292 * t392 + t276;
t293 = t380 + t300;
t462 = t293 * t392 + t278;
t461 = -t450 + t510;
t352 = t425 * t356;
t460 = t346 * t427 - t352;
t319 = t374 * t424 - t518 * t435;
t332 = t389 * t424 - t518 * t390;
t355 = -t518 * t412 - t413 * t424;
t459 = t388 * t425;
t419 = -pkin(2) * t518 - pkin(3);
t453 = t292 * t427 - t293 * t425;
t452 = -qJ(6) * t396 - qJD(6) * t407;
t451 = t337 * t427 - t476;
t448 = t396 * t425 + t472;
t447 = -t396 * t427 + t473;
t281 = t322 * pkin(4) - qJ(5) * t450 - t362 * qJD(5) + t319;
t444 = -t356 * t488 + t475;
t443 = -t388 * t456 - t504;
t442 = t419 - t503;
t438 = t301 * t362 + t278;
t298 = t450 + t511;
t277 = -pkin(5) * t322 - t281;
t436 = t299 * t388 + t445;
t404 = t499 * t425;
t403 = t442 - t521;
t381 = -t442 + t532;
t308 = pkin(4) * t362 + t516;
t307 = t407 * t454 + t355;
t305 = t407 * t449 - t355;
t304 = -t362 * t522 - t516;
t303 = pkin(4) * t446 - t460;
t296 = -pkin(4) * t394 - t336 * t427 + t344;
t294 = qJ(6) * t506 + t302;
t291 = t352 + (-qJ(6) * t407 - t346) * t427 + t522 * t446;
t284 = t454 * t396 + (-t485 + (t503 + t521) * qJD(4)) * t407 + t332;
t282 = t449 * t396 + (t485 + (-t503 - t532) * qJD(4)) * t407 - t332;
t280 = -pkin(4) * t393 - t451;
t279 = t444 + t526;
t275 = qJ(6) * t472 + (-qJD(4) * t356 - t452) * t425 + t475 + t526;
t274 = qJ(6) * t473 - t522 * t393 + (-t337 + t452) * t427 + t476;
t1 = [0.2e1 * t428 * MDP(4) * t469 + t534 * t535 + MDP(6) * t501 - MDP(7) * t502 + (-pkin(7) * t501 + t426 * t464) * MDP(9) + (pkin(7) * t502 + t428 * t464) * MDP(10) + (t319 * t407 + t320 * t446 + t332 * t394 + t333 * t392 - t347 * t396 - t348 * t393 + t355 * t440 - t356 * t382) * MDP(11) + (t319 * t355 + t320 * t356 - t347 * t332 + t348 * t333 + (t411 + t455) * t477) * MDP(12) + (-t447 * t362 + t450 * t505) * MDP(13) + ((-t360 * t427 - t362 * t425) * t396 + (-t512 - t322 * t427 + (t360 * t425 - t362 * t427) * qJD(4)) * t407) * MDP(14) + (t362 * t393 + t371 * t407 - t388 * t447 - t446 * t450) * MDP(15) + (t322 * t446 - t360 * t393 - t370 * t407 - t388 * t448) * MDP(16) + (-t382 * t446 + t388 * t393) * MDP(17) + (t299 * t393 + t319 * t506 + t355 * t322 + t332 * t360 + t382 * t460 + t388 * t451 + t446 * t457 - t448 * t456) * MDP(18) + (-t300 * t393 + t319 * t505 + t332 * t362 + t355 * t450 - t382 * t491 - t388 * t444 - t445 * t446 + t447 * t456) * MDP(19) + (t278 * t446 - t280 * t388 + t281 * t506 + t284 * t360 - t292 * t393 + t301 * t448 - t303 * t382 + t307 * t322) * MDP(20) + (-t279 * t360 + t280 * t362 - t302 * t322 + t303 * t450 + t453 * t396 + (-t276 * t425 + t278 * t427 + (-t292 * t425 - t293 * t427) * qJD(4)) * t407) * MDP(21) + (-t276 * t446 + t279 * t388 - t281 * t505 - t284 * t362 + t293 * t393 + t301 * t447 + t302 * t382 - t307 * t450) * MDP(22) + (t276 * t302 + t278 * t303 + t279 * t293 + t280 * t292 + t281 * t307 + t284 * t301) * MDP(23) + (t272 * t446 - t274 * t388 - t277 * t506 - t282 * t360 - t283 * t393 - t288 * t448 - t291 * t382 - t305 * t322) * MDP(24) + (-t273 * t446 + t275 * t388 + t277 * t505 + t282 * t362 + t285 * t393 - t288 * t447 + t294 * t382 + t305 * t450) * MDP(25) + (-t274 * t362 + t275 * t360 - t291 * t450 + t294 * t322 + (-t283 * t427 + t285 * t425) * t396 + (-t272 * t427 + t273 * t425 + (t283 * t425 + t285 * t427) * qJD(4)) * t407) * MDP(26) + (t272 * t291 + t273 * t294 + t274 * t283 + t275 * t285 + t277 * t305 + t282 * t288) * MDP(27); -t426 * MDP(4) * t500 + t431 * t534 + ((t348 - t349) * t394 + (-t350 + t347) * t392 + (-t424 * t382 - t440 * t518) * pkin(2)) * MDP(11) + (t347 * t349 - t348 * t350 + (-t319 * t518 + t320 * t424 - t411 * t484) * pkin(2)) * MDP(12) + (t362 * t458 + t512) * MDP(13) + ((t450 + t510) * t427 - t362 * t459 + t494) * MDP(14) + (-t362 * t394 - t392 * t458 + t490) * MDP(15) + (t360 * t394 - t470 + t493) * MDP(16) - t388 * t394 * MDP(17) + (-t299 * t394 - t319 * t427 + t419 * t322 - t349 * t360 + (t344 + (-qJD(4) * t418 - t336) * t427) * t388 + t443 * t425) * MDP(18) + (t300 * t394 + t319 * t425 + t419 * t450 - t349 * t362 + (t471 + t492) * t388 + t443 * t427) * MDP(19) + (-t281 * t427 + t292 * t394 + t322 * t403 + (-t418 * t487 + t296) * t388 - t495 * t360 + t524 * t425) * MDP(20) + (t295 * t360 - t296 * t362 + (-t322 * t418 + (t362 * t418 + t292) * qJD(4) + t463) * t427 + (t450 * t418 + (t360 * t418 - t293) * qJD(4) + t462) * t425) * MDP(21) + (-t281 * t425 - t293 * t394 - t450 * t403 + (-t295 - t471) * t388 + t495 * t362 - t524 * t427) * MDP(22) + (t281 * t403 - t292 * t296 - t293 * t295 - t495 * t301 + (qJD(4) * t453 + t276 * t427 + t278 * t425) * t418) * MDP(23) + (t277 * t427 + t283 * t394 - t288 * t459 - t322 * t381 + t360 * t496 + t382 * t404 + t388 * t498) * MDP(24) + (t277 * t425 - t285 * t394 + t288 * t458 - t362 * t496 + t381 * t450 - t382 * t405 + t388 * t497) * MDP(25) + (-t322 * t405 + t497 * t360 + t498 * t362 + t450 * t404 + t425 * t539 - t538 * t427) * MDP(26) + (-t272 * t404 - t273 * t405 + t277 * t381 - t283 * t498 + t285 * t497 - t288 * t496) * MDP(27) + (MDP(9) * t426 * t431 + MDP(10) * t500) * pkin(1); -t392 ^ 2 * MDP(11) + (-t348 * t392 + t416) * MDP(12) + t493 * MDP(18) + t494 * MDP(21) - t478 * t470 + (-MDP(11) * t394 + t347 * MDP(12) - t301 * MDP(23) + t288 * MDP(27) - t465 * t362 + (-MDP(18) - t478) * t360) * t394 + (-t382 * MDP(19) + (qJD(4) * t292 + t463) * MDP(23) + t322 * MDP(26) + t538 * MDP(27) + (-qJD(4) * MDP(18) + t362 * t525 + t478 * t392) * t388) * t425 + (t461 * MDP(21) + (qJD(4) * t293 - t462) * MDP(23) + (qJD(4) * t360 - t461) * MDP(26) + t539 * MDP(27) + t478 * t382 + (-qJD(4) * MDP(19) + t392 * t465) * t388) * t427 + t536 * t490; MDP(13) * t509 + (t359 - t523) * MDP(14) + t298 * MDP(15) + (-t322 + t508) * MDP(16) + t382 * MDP(17) + (t362 * t456 - t457 + t513) * MDP(18) + (-t360 * t456 + t436) * MDP(19) + (-t308 * t360 + t377 - t438 + t513) * MDP(20) + (-pkin(4) * t450 - t517 + (t293 - t300) * t362 + (t292 - t481) * t360) * MDP(21) + (-t301 * t360 + t308 * t362 - t436 + t528) * MDP(22) + (-pkin(4) * t278 + qJ(5) * t276 - t292 * t300 + t293 * t481 - t301 * t308) * MDP(23) + (t290 * t388 + t304 * t360 + (pkin(5) + t522) * t382 + t531 - t439) * MDP(24) + (t288 * t360 - t289 * t388 - t304 * t362 - t445 + t528 + t530) * MDP(25) + (t517 + t450 * t522 + (-t285 + t290) * t362 + (-t283 + t482) * t360) * MDP(26) + (qJ(5) * t273 - t272 * t522 - t283 * t290 + t285 * t482 - t288 * t304) * MDP(27); (-t293 * t388 + t438) * MDP(23) + (t434 - t514 - t531) * MDP(27) + t478 * (-qJD(2) * t394 + t509) + t525 * t298 + t536 * (-t388 ^ 2 - t359); (-t322 - t508) * MDP(24) + (t450 - t511) * MDP(25) + (-t359 - t523) * MDP(26) + (t283 * t362 - t285 * t360 + t277) * MDP(27);];
tauc  = t1;
