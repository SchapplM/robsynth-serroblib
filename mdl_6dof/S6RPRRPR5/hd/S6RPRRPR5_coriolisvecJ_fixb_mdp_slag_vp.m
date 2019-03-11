% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR5
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
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:41
% EndTime: 2019-03-09 05:13:51
% DurationCPUTime: 5.50s
% Computational Cost: add. (5201->388), mult. (13723->489), div. (0->0), fcn. (10873->8), ass. (0->164)
t387 = qJD(3) + qJD(4);
t391 = sin(qJ(6));
t449 = qJD(6) * t391;
t389 = sin(pkin(10));
t390 = cos(pkin(10));
t480 = sin(qJ(3));
t482 = cos(qJ(3));
t414 = t389 * t480 - t390 * t482;
t367 = t414 * qJD(1);
t413 = -t389 * t482 - t390 * t480;
t368 = t413 * qJD(1);
t392 = sin(qJ(4));
t401 = t367 * qJD(3);
t481 = cos(qJ(4));
t436 = qJD(4) * t481;
t450 = qJD(4) * t392;
t437 = qJD(3) * t480;
t438 = qJD(3) * t482;
t456 = (t389 * t438 + t390 * t437) * qJD(1);
t310 = -t367 * t450 - t368 * t436 - t392 * t401 + t481 * t456;
t346 = t481 * t367 - t368 * t392;
t393 = cos(qJ(6));
t448 = qJD(6) * t393;
t457 = t391 * t310 + t346 * t448;
t287 = -t387 * t449 + t457;
t286 = t287 * t393;
t308 = t393 * t310;
t337 = t346 * t391 + t387 * t393;
t288 = t337 * qJD(6) - t308;
t402 = t481 * t414;
t400 = qJD(3) * t402;
t309 = qJD(1) * t400 + t367 * t436 - t368 * t450 + t392 * t456;
t465 = t346 * t387;
t405 = -t309 + t465;
t418 = -t392 * t367 - t368 * t481;
t467 = t418 * t387;
t488 = qJD(6) + t418;
t501 = t488 * t393;
t502 = t488 * t391;
t462 = t387 * t391;
t335 = -t393 * t346 + t462;
t503 = t335 * t488;
t506 = t405 * MDP(17) + (-t310 + t467) * MDP(18) + (-t337 * t502 + t286) * MDP(26) + (-t337 * t501 + (-t287 + t503) * t391 - t393 * t288) * MDP(27);
t307 = t393 * t309;
t500 = -t488 * t502 - t307;
t461 = t391 * t309;
t411 = -t488 * t501 + t461;
t479 = pkin(5) * t346;
t476 = qJ(5) * t346;
t478 = pkin(7) + qJ(2);
t374 = t478 * t389;
t371 = qJD(1) * t374;
t375 = t478 * t390;
t372 = qJD(1) * t375;
t334 = -t367 * pkin(8) - t371 * t480 + t372 * t482;
t331 = t481 * t334;
t489 = -t482 * t371 - t480 * t372;
t333 = t368 * pkin(8) + t489;
t332 = qJD(3) * pkin(3) + t333;
t298 = t392 * t332 + t331;
t296 = -qJ(5) * t387 - t298;
t279 = -t296 - t479;
t504 = t279 * t488;
t330 = t392 * t334;
t297 = -t481 * t332 + t330;
t446 = qJD(5) + t297;
t499 = (t389 ^ 2 + t390 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t496 = t418 ^ 2;
t495 = pkin(4) * t418;
t494 = pkin(5) * t418;
t493 = t279 * t418;
t483 = pkin(4) + pkin(9);
t491 = t418 * t483;
t300 = t392 * t333 + t331;
t422 = pkin(3) * t450 - t300;
t301 = t333 * t481 - t330;
t458 = -pkin(3) * t436 - qJD(5) + t301;
t439 = qJD(2) * t480;
t440 = qJD(2) * t482;
t490 = -t389 * t440 - t390 * t439;
t447 = t494 + t446;
t318 = -t456 * pkin(8) - qJD(2) * t367 + qJD(3) * t489;
t319 = pkin(8) * t401 + qJD(1) * t490 + t371 * t437 - t372 * t438;
t272 = t392 * t318 - t481 * t319 + t332 * t450 + t334 * t436;
t382 = -t390 * pkin(2) - pkin(1);
t373 = qJD(1) * t382 + qJD(2);
t353 = pkin(3) * t367 + t373;
t399 = -qJ(5) * t418 + t353;
t299 = pkin(4) * t346 + t399;
t416 = t299 * t418 + t272;
t487 = -t418 * t353 - t272;
t277 = -t387 * t483 + t447;
t282 = t346 * t483 + t399;
t266 = t277 * t393 - t282 * t391;
t267 = t277 * t391 + t282 * t393;
t486 = MDP(15) * t418 + t353 * MDP(21) - t299 * MDP(24) + MDP(30) * t488 + t266 * MDP(31) - t267 * MDP(32);
t398 = -t374 * t438 - t375 * t437 - t389 * t439 + t390 * t440;
t408 = qJD(3) * t413;
t322 = pkin(8) * t408 + t398;
t409 = qJD(3) * t414;
t323 = pkin(8) * t409 + t374 * t437 - t375 * t438 + t490;
t342 = pkin(8) * t413 - t374 * t482 - t375 * t480;
t415 = t374 * t480 - t375 * t482;
t343 = -pkin(8) * t414 - t415;
t274 = -t481 * t322 - t392 * t323 - t342 * t436 + t343 * t450;
t475 = t274 * t387;
t419 = t392 * t342 + t343 * t481;
t275 = qJD(4) * t419 + t392 * t322 - t323 * t481;
t474 = t275 * t387;
t351 = -t392 * t413 + t402;
t410 = t392 * t414;
t352 = -t413 * t481 - t410;
t357 = pkin(3) * t414 + t382;
t397 = -t352 * qJ(5) + t357;
t290 = t351 * t483 + t397;
t473 = t290 * t309;
t472 = t298 * t387;
t471 = t335 * t346;
t470 = t337 * t346;
t466 = t346 ^ 2;
t464 = t351 * t391;
t383 = -pkin(3) * t481 - pkin(4);
t380 = -pkin(9) + t383;
t463 = t380 * t309;
t460 = t483 * t309;
t427 = t481 * t318 + t392 * t319 + t332 * t436 - t334 * t450;
t271 = -t387 * qJD(5) - t427;
t262 = -pkin(5) * t310 - t271;
t459 = t262 * t391 + t279 * t448;
t455 = -t458 + t494;
t451 = qJD(3) * t368;
t443 = qJD(1) * qJD(2);
t435 = t456 * pkin(3);
t421 = -pkin(3) * t368 + t476;
t429 = -qJD(6) * t380 + t421 + t491;
t428 = qJD(6) * t483 + t476 + t491;
t423 = t422 + t479;
t304 = -t342 * t481 + t392 * t343;
t317 = -qJD(4) * t410 - t392 * t409 - t408 * t481 - t413 * t436;
t417 = t391 * t317 + t351 * t448;
t291 = t352 * pkin(5) + t304;
t412 = t262 * t351 + t279 * t317 + t291 * t309;
t407 = t413 * qJD(2);
t406 = t309 * qJ(5) - qJD(5) * t418 + t435;
t403 = pkin(3) * t408;
t273 = t310 * pkin(4) + t406;
t316 = qJD(4) * t402 - t392 * t408 - t413 * t450 + t400;
t276 = t317 * pkin(4) + t316 * qJ(5) - t352 * qJD(5) - t403;
t381 = pkin(3) * t392 + qJ(5);
t311 = t476 + t495;
t303 = t351 * pkin(4) + t397;
t302 = t421 + t495;
t295 = -pkin(4) * t387 + t446;
t294 = t309 * t352;
t292 = -t351 * pkin(5) + t419;
t281 = t298 - t479;
t270 = t317 * pkin(9) + t276;
t269 = -t316 * pkin(5) + t275;
t268 = -pkin(5) * t317 - t274;
t265 = t310 * t483 + t406;
t264 = -pkin(5) * t309 + t272;
t263 = t393 * t264;
t261 = t262 * t393;
t1 = [(t263 * t352 - t266 * t316 + t268 * t335 + t292 * t288 + (-t265 * t352 - t270 * t488 + t473) * t391 + (t269 * t488 - t412) * t393 + ((-t290 * t393 - t291 * t391) * t488 - t267 * t352 + t279 * t464) * qJD(6)) * MDP(31) + (t267 * t316 + t268 * t337 + t292 * t287 + (-(qJD(6) * t291 + t270) * t488 + t473 - (qJD(6) * t277 + t265) * t352 + t279 * qJD(6) * t351) * t393 + (-(-qJD(6) * t290 + t269) * t488 - (-qJD(6) * t282 + t264) * t352 + t412) * t391) * MDP(32) + ((-t335 * t391 + t337 * t393) * t317 + (t286 - t288 * t391 + (-t335 * t393 - t337 * t391) * qJD(6)) * t351) * MDP(27) + (-t273 * t351 - t276 * t346 - t299 * t317 - t303 * t310 + t474) * MDP(23) + (t287 * t464 + t337 * t417) * MDP(26) + (-t351 * t307 - t288 * t352 + t316 * t335 + (t393 * t317 - t351 * t449) * t488) * MDP(29) + (t287 * t352 - t316 * t337 - t351 * t461 + t417 * t488) * MDP(28) + (-t316 * t488 - t294) * MDP(30) + (t367 * t409 - t368 * t408 + t413 * t456) * MDP(9) + (-t373 * t409 - t382 * t401) * MDP(14) + (t357 * t310 + t353 * t317 - t346 * t403 + t351 * t435 - t474) * MDP(20) + (-t273 * t352 - t276 * t418 + t299 * t316 + t303 * t309 - t475) * MDP(24) + (-t316 * t418 - t294) * MDP(15) + (t309 * t351 - t310 * t352 + t316 * t346 - t317 * t418) * MDP(16) + (-t271 * t419 + t272 * t304 + t273 * t303 + t274 * t296 + t275 * t295 + t276 * t299) * MDP(25) + (t271 * t351 + t272 * t352 + t274 * t346 + t275 * t418 - t295 * t316 + t296 * t317 - t304 * t309 - t310 * t419) * MDP(22) + (-t357 * t309 - t353 * t316 + t352 * t435 - t403 * t418 + t475) * MDP(21) + t382 * t456 * MDP(13) + 0.2e1 * t443 * t499 + (-MDP(17) * t316 - MDP(18) * t317) * t387 + ((t367 * t413 + t368 * t414) * MDP(8) + qJD(1) * t414 ^ 2 * MDP(9) - t398 * MDP(14) + (-t373 * t413 + t407) * MDP(13) + (-t414 * MDP(10) + t413 * MDP(11) + t415 * MDP(13)) * qJD(3)) * qJD(3); (-t451 + t456) * MDP(13) - 0.2e1 * MDP(14) * t401 + (-t466 - t496) * MDP(22) + (-t295 * t418 - t296 * t346 + t273) * MDP(25) + (t411 + t471) * MDP(31) + (t470 - t500) * MDP(32) - qJD(1) ^ 2 * t499 + (-MDP(21) + MDP(24)) * (t309 + t465) + (MDP(20) - MDP(23)) * (t310 + t467); (t470 + t500) * MDP(28) + (-t271 * t381 + t272 * t383 + t295 * t422 + t296 * t458 - t299 * t302) * MDP(25) + (t381 * t287 + t261 + t429 * t501 + t455 * t337 + (-t423 * t488 + t463 - t504) * t391) * MDP(32) + (t381 * t288 + (-t463 + t493) * t393 + t455 * t335 + (t391 * t429 + t393 * t423) * t488 + t459) * MDP(31) - (-t295 * MDP(22) - t486) * t346 + (-t466 + t496) * MDP(16) + (t373 * t367 + t414 * t443) * MDP(14) + (qJD(1) * t407 + t373 * t368) * MDP(13) + (t301 * t387 + (t368 * t418 - t387 * t436) * pkin(3) - t427) * MDP(21) + (t302 * t418 - t387 * t458 - t271) * MDP(24) + (-t309 * t383 - t310 * t381 + t346 * t458 + (-t296 + t422) * t418) * MDP(22) + (-t451 - t456) * MDP(11) + (t411 - t471) * MDP(29) + (t300 * t387 + (t346 * t368 - t387 * t450) * pkin(3) + t487) * MDP(20) + (t302 * t346 + t387 * t422 + t416) * MDP(23) - t368 * t367 * MDP(8) + (-t367 ^ 2 + t368 ^ 2) * MDP(9) + t506; t496 * MDP(16) + (t472 + t487) * MDP(20) + (-t297 * t387 - t427) * MDP(21) + (pkin(4) * t309 - qJ(5) * t310 + (-t296 - t298) * t418) * MDP(22) + (t416 - t472) * MDP(23) + (t311 * t418 + t387 * t446 - t271) * MDP(24) + (-pkin(4) * t272 - qJ(5) * t271 - t295 * t298 - t296 * t446 - t299 * t311) * MDP(25) + t500 * MDP(28) + t411 * MDP(29) + (qJ(5) * t288 + (t460 + t493) * t393 + (-t393 * t281 + t391 * t428) * t488 + t447 * t335 + t459) * MDP(31) + (qJ(5) * t287 + t261 + t428 * t501 + t447 * t337 + (t281 * t488 - t460 - t504) * t391) * MDP(32) + ((t295 - t446) * MDP(22) + t311 * MDP(23) + t337 * MDP(28) - t335 * MDP(29) - MDP(16) * t346 + t486) * t346 + t506; t405 * MDP(22) - t418 * t346 * MDP(23) + (-t387 ^ 2 - t496) * MDP(24) + (t296 * t387 + t416) * MDP(25) + (-t335 * t387 - t307) * MDP(31) + (-t337 * t387 + t461) * MDP(32) + (-MDP(31) * t502 - MDP(32) * t501) * t488; t337 * t335 * MDP(26) + (-t335 ^ 2 + t337 ^ 2) * MDP(27) + (t457 + t503) * MDP(28) + (t337 * t488 + t308) * MDP(29) - t309 * MDP(30) + (-t265 * t391 + t267 * t488 - t279 * t337 + t263) * MDP(31) + (-t264 * t391 - t265 * t393 + t266 * t488 + t279 * t335) * MDP(32) + (-MDP(28) * t462 - MDP(29) * t337 - MDP(31) * t267 - MDP(32) * t266) * qJD(6);];
tauc  = t1;
