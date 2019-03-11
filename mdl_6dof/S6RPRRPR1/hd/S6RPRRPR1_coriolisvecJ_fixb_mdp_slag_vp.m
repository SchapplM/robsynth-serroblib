% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:52
% EndTime: 2019-03-09 04:58:59
% DurationCPUTime: 3.34s
% Computational Cost: add. (4923->330), mult. (11844->453), div. (0->0), fcn. (8622->10), ass. (0->169)
t368 = sin(pkin(10)) * pkin(1) + pkin(7);
t474 = pkin(8) + t368;
t376 = qJD(3) + qJD(4);
t385 = cos(qJ(6));
t382 = sin(qJ(6));
t386 = cos(qJ(4));
t387 = cos(qJ(3));
t447 = qJD(1) * t387;
t437 = t386 * t447;
t383 = sin(qJ(4));
t384 = sin(qJ(3));
t448 = qJD(1) * t384;
t438 = t383 * t448;
t346 = -t437 + t438;
t348 = -t383 * t447 - t386 * t448;
t379 = sin(pkin(11));
t471 = cos(pkin(11));
t402 = -t379 * t346 - t348 * t471;
t462 = t402 * t382;
t304 = -t385 * t376 + t462;
t427 = -t471 * t346 + t348 * t379;
t477 = qJD(6) - t427;
t485 = t304 * t477;
t306 = t376 * t382 + t385 * t402;
t484 = t306 * t477;
t426 = t477 * t385;
t441 = qJD(1) * qJD(3);
t436 = t387 * t441;
t313 = qJD(4) * t437 - t376 * t438 + t386 * t436;
t353 = t383 * t387 + t384 * t386;
t478 = qJD(1) * t353;
t392 = t376 * t478;
t289 = t313 * t379 + t392 * t471;
t457 = t382 * t289;
t483 = -t477 * t426 - t457;
t482 = MDP(5) * t384;
t481 = MDP(6) * (t384 ^ 2 - t387 ^ 2);
t431 = t474 * qJD(1);
t328 = t387 * qJD(2) - t431 * t384;
t322 = t328 * qJD(3);
t472 = qJD(3) * pkin(3);
t327 = t328 + t472;
t479 = t386 * (qJD(4) * t327 + t322);
t342 = t348 * qJ(5);
t329 = qJD(2) * t384 + t387 * t431;
t324 = t383 * t329;
t430 = t386 * t327 - t324;
t294 = t342 + t430;
t323 = t329 * qJD(3);
t445 = qJD(4) * t383;
t429 = -t383 * t323 - t329 * t445;
t260 = -qJ(5) * t392 - t346 * qJD(5) + t429 + t479;
t326 = t386 * t329;
t411 = -t327 * t383 - t326;
t412 = -t383 * t322 - t386 * t323;
t395 = qJD(4) * t411 + t412;
t391 = -qJ(5) * t313 + qJD(5) * t348 + t395;
t249 = t260 * t471 + t379 * t391;
t400 = t353 * qJD(4);
t320 = qJD(3) * t353 + t400;
t352 = t383 * t384 - t386 * t387;
t433 = qJD(3) * t474;
t343 = t384 * t433;
t344 = t387 * t433;
t351 = t474 * t387;
t350 = t474 * t384;
t460 = t350 * t386;
t401 = -qJD(4) * t460 - t386 * t343 - t383 * t344 - t351 * t445;
t276 = -qJ(5) * t320 - qJD(5) * t352 + t401;
t319 = t376 * t352;
t410 = t350 * t383 - t351 * t386;
t396 = qJD(4) * t410 + t343 * t383 - t386 * t344;
t393 = qJ(5) * t319 - qJD(5) * t353 + t396;
t256 = t276 * t471 + t379 * t393;
t288 = pkin(4) * t376 + t294;
t470 = qJ(5) * t346;
t295 = -t411 - t470;
t459 = t379 * t295;
t267 = t288 * t471 - t459;
t264 = -t376 * pkin(5) - t267;
t370 = -cos(pkin(10)) * pkin(1) - pkin(2);
t354 = -pkin(3) * t387 + t370;
t349 = t354 * qJD(1);
t315 = pkin(4) * t346 + qJD(5) + t349;
t278 = -pkin(5) * t427 - pkin(9) * t402 + t315;
t316 = t352 * t471 + t353 * t379;
t317 = -t379 * t352 + t353 * t471;
t399 = pkin(4) * t352 + t354;
t283 = pkin(5) * t316 - pkin(9) * t317 + t399;
t298 = -t319 * t471 - t379 * t320;
t248 = t260 * t379 - t391 * t471;
t302 = -qJ(5) * t352 - t410;
t398 = -qJ(5) * t353 - t351 * t383 - t460;
t280 = t302 * t471 + t379 * t398;
t418 = t248 * t317 - t280 * t289;
t476 = t264 * t298 - (qJD(6) * t283 + t256) * t477 - (qJD(6) * t278 + t249) * t316 + t418;
t475 = pkin(4) * t348;
t473 = pkin(3) * qJD(4);
t469 = t264 * t427;
t468 = t264 * t317;
t444 = qJD(6) * t382;
t290 = t313 * t471 - t379 * t392;
t443 = qJD(6) * t385;
t453 = t385 * t290 + t376 * t443;
t273 = -t402 * t444 + t453;
t467 = t273 * t382;
t466 = t283 * t289;
t465 = t290 * t382;
t464 = t304 * t402;
t463 = t306 * t402;
t461 = t349 * t348;
t458 = t379 * t383;
t388 = qJD(3) ^ 2;
t456 = t384 * t388;
t285 = t385 * t289;
t455 = t387 * t388;
t297 = -t319 * t379 + t320 * t471;
t454 = t273 * t316 + t306 * t297;
t291 = t471 * t295;
t268 = t379 * t288 + t291;
t452 = t386 * t328 - t324;
t296 = t342 + t452;
t428 = -t328 * t383 - t326;
t407 = t428 + t470;
t432 = t471 * t383;
t451 = -t296 * t379 + t407 * t471 + (t379 * t386 + t432) * t473;
t450 = -t296 * t471 - t379 * t407 + (t386 * t471 - t458) * t473;
t372 = pkin(3) * t386 + pkin(4);
t341 = pkin(3) * t432 + t379 * t372;
t356 = qJD(1) * t370;
t374 = t384 * t472;
t373 = pkin(3) * t448;
t440 = t317 * t457;
t439 = t317 * t285;
t435 = -pkin(3) * t376 - t327;
t434 = pkin(4) * t320 + t374;
t282 = pkin(5) * t402 - pkin(9) * t427 - t475;
t337 = pkin(9) + t341;
t422 = qJD(6) * t337 + t282 + t373;
t367 = pkin(4) * t379 + pkin(9);
t421 = qJD(6) * t367 + t282;
t265 = pkin(9) * t376 + t268;
t253 = t265 * t385 + t278 * t382;
t420 = t248 * t382 + t253 * t402 + t264 * t443;
t417 = -t289 * t337 - t469;
t416 = -t289 * t367 - t469;
t415 = t265 * t382 - t278 * t385;
t414 = t267 * t427 + t268 * t402;
t274 = qJD(6) * t306 + t465;
t413 = -t274 * t316 - t297 * t304;
t409 = 0.2e1 * qJD(3) * t356;
t408 = t285 + (t382 * t427 - t444) * t477;
t406 = -t248 * t385 + t264 * t444 + t402 * t415;
t405 = t349 * t346 - t429;
t404 = -t298 * t382 - t317 * t443;
t403 = -t298 * t385 + t317 * t444;
t340 = -pkin(3) * t458 + t372 * t471;
t394 = -t348 * t346 * MDP(12) - t477 * t402 * MDP(25) + ((t273 - t485) * t385 + (-t274 - t484) * t382) * MDP(22) + (t408 + t464) * MDP(24) + (-t463 - t483) * MDP(23) + (t306 * t426 + t467) * MDP(21) + (t346 * t376 + t313) * MDP(14) + (-t348 * t376 - t392) * MDP(15) + (-t346 ^ 2 + t348 ^ 2) * MDP(13);
t390 = pkin(4) * t392 + qJD(3) * t373;
t369 = -pkin(4) * t471 - pkin(5);
t336 = -pkin(5) - t340;
t279 = t302 * t379 - t398 * t471;
t270 = t294 * t471 - t459;
t269 = t294 * t379 + t291;
t261 = pkin(5) * t297 - pkin(9) * t298 + t434;
t258 = t289 * pkin(5) - t290 * pkin(9) + t390;
t257 = t385 * t258;
t255 = t276 * t379 - t393 * t471;
t1 = [0.2e1 * t436 * t482 - 0.2e1 * t441 * t481 + MDP(7) * t455 - MDP(8) * t456 + (-t368 * t455 + t384 * t409) * MDP(10) + (t368 * t456 + t387 * t409) * MDP(11) + (t313 * t353 + t319 * t348) * MDP(12) + (-t313 * t352 + t319 * t346 + t348 * t320 - t353 * t392) * MDP(13) + (t346 * t374 + t349 * t320 + (t354 * t400 + (t384 * pkin(3) * t352 + t353 * t354) * qJD(3)) * qJD(1)) * MDP(17) + (t354 * t313 - t349 * t319 + (-t348 + t478) * t374) * MDP(18) + (-t249 * t316 + t255 * t402 + t256 * t427 - t267 * t298 - t268 * t297 + t279 * t290 + t418) * MDP(19) + (t248 * t279 + t249 * t280 - t267 * t255 + t268 * t256 + t315 * t434 + t390 * t399) * MDP(20) + (t273 * t317 * t385 - t306 * t403) * MDP(21) + ((-t304 * t385 - t306 * t382) * t298 + (-t467 - t274 * t385 + (t304 * t382 - t306 * t385) * qJD(6)) * t317) * MDP(22) + (-t403 * t477 + t439 + t454) * MDP(23) + (t404 * t477 + t413 - t440) * MDP(24) + (t289 * t316 + t297 * t477) * MDP(25) + (-t415 * t297 + t255 * t304 + t257 * t316 + t279 * t274 + (t261 * t477 + t466 + (-t265 * t316 - t280 * t477 + t468) * qJD(6)) * t385 + t476 * t382) * MDP(26) + (-t253 * t297 + t255 * t306 + t279 * t273 + (-(-qJD(6) * t280 + t261) * t477 - t466 - (-qJD(6) * t265 + t258) * t316 - qJD(6) * t468) * t382 + t476 * t385) * MDP(27) + (-t319 * MDP(14) - t320 * MDP(15) + MDP(17) * t396 - MDP(18) * t401) * t376; (-t289 * t317 + t290 * t316 + t297 * t402 + t298 * t427) * MDP(19) + (t248 * t316 + t249 * t317 - t267 * t297 + t268 * t298) * MDP(20) + (-t413 - t440) * MDP(26) + (-t439 + t454) * MDP(27) + (-MDP(10) * t384 - MDP(11) * t387) * t388 + (-t320 * MDP(17) + t319 * MDP(18)) * t376 + (MDP(26) * t404 + MDP(27) * t403) * t477; (t249 * t341 - t248 * t340 - t315 * (t373 - t475) + t450 * t268 - t451 * t267) * MDP(20) + (t336 * t273 + t417 * t385 + t451 * t306 + (t382 * t422 - t385 * t450) * t477 + t420) * MDP(27) + (t348 * t373 + t452 * t376 + (qJD(4) * t435 - t322) * t386 + t405) * MDP(18) + (-t346 * t373 + t461 - t428 * t376 + (t383 * t435 - t326) * qJD(4) + t412) * MDP(17) + (t336 * t274 + t417 * t382 + t451 * t304 + (-t382 * t450 - t385 * t422) * t477 + t406) * MDP(26) + t394 + (-t289 * t341 - t290 * t340 + t402 * t451 + t427 * t450 + t414) * MDP(19) + (-t387 * t482 + t481) * qJD(1) ^ 2 + (-MDP(10) * t448 - MDP(11) * t447) * t356; (-t269 * t402 - t270 * t427 + (-t289 * t379 - t290 * t471) * pkin(4) + t414) * MDP(19) + (-t376 * t411 + t395 + t461) * MDP(17) + (t376 * t430 + t405 - t479) * MDP(18) + (t267 * t269 - t268 * t270 + (-t248 * t471 + t249 * t379 + t315 * t348) * pkin(4)) * MDP(20) + (-t269 * t306 + t369 * t273 + t416 * t385 + (t385 * t270 + t382 * t421) * t477 + t420) * MDP(27) + (-t269 * t304 + t369 * t274 + t416 * t382 + (t382 * t270 - t385 * t421) * t477 + t406) * MDP(26) + t394; (-t402 ^ 2 - t427 ^ 2) * MDP(19) + (t267 * t402 - t268 * t427 + t390) * MDP(20) + (t408 - t464) * MDP(26) + (-t463 + t483) * MDP(27); t306 * t304 * MDP(21) + (-t304 ^ 2 + t306 ^ 2) * MDP(22) + (t453 + t485) * MDP(23) + (-t465 + t484) * MDP(24) + t289 * MDP(25) + (-t249 * t382 + t253 * t477 - t264 * t306 + t257) * MDP(26) + (-t249 * t385 - t258 * t382 + t264 * t304 - t415 * t477) * MDP(27) + (-MDP(23) * t462 - MDP(24) * t306 - MDP(26) * t253 + MDP(27) * t415) * qJD(6);];
tauc  = t1;
