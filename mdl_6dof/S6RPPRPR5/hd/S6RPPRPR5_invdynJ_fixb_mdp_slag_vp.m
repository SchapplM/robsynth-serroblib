% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:30
% EndTime: 2019-03-09 01:49:36
% DurationCPUTime: 4.31s
% Computational Cost: add. (2175->445), mult. (3957->570), div. (0->0), fcn. (2488->10), ass. (0->184)
t376 = sin(qJ(4));
t452 = qJD(1) * t376;
t346 = qJD(6) + t452;
t371 = sin(pkin(9));
t379 = cos(qJ(4));
t451 = qJD(1) * t379;
t427 = t371 * t451;
t372 = cos(pkin(9));
t443 = t372 * qJD(4);
t319 = t427 - t443;
t426 = t372 * t451;
t448 = qJD(4) * t371;
t321 = t426 + t448;
t375 = sin(qJ(6));
t378 = cos(qJ(6));
t397 = t319 * t375 - t321 * t378;
t486 = t346 * t397;
t408 = pkin(4) * t379 + qJ(5) * t376;
t306 = t408 * qJD(4) - qJD(5) * t379 + qJD(3);
t474 = qJ(5) * t379;
t407 = pkin(4) * t376 - t474;
t364 = qJDD(1) * qJ(3);
t370 = qJDD(1) * pkin(1);
t432 = qJDD(2) - t370;
t418 = -t364 + t432;
t264 = qJD(1) * t306 + qJDD(1) * t407 - t418;
t365 = qJD(1) * qJD(2);
t366 = qJ(2) * qJDD(1);
t417 = qJDD(3) + t365 + t366;
t327 = -pkin(7) * qJDD(1) + t417;
t349 = qJ(2) * qJD(1) + qJD(3);
t340 = -pkin(7) * qJD(1) + t349;
t466 = t340 * t379;
t273 = qJDD(4) * qJ(5) + t327 * t376 + (qJD(5) + t466) * qJD(4);
t254 = t372 * t264 - t273 * t371;
t255 = t371 * t264 + t372 * t273;
t377 = sin(qJ(1));
t380 = cos(qJ(1));
t482 = g(1) * t377 - g(2) * t380;
t485 = -t254 * t372 - t255 * t371 - t482;
t373 = -pkin(7) + qJ(2);
t446 = qJD(4) * t379;
t450 = qJD(2) * t376;
t484 = t373 * t446 + t450;
t460 = t378 * t372;
t324 = t371 * t375 - t460;
t391 = t346 * t324;
t410 = g(1) * t380 + g(2) * t377;
t483 = t410 * t372;
t444 = qJD(6) * t378;
t445 = qJD(6) * t375;
t481 = -t371 * t445 + t372 * t444;
t374 = pkin(1) + qJ(3);
t480 = qJD(1) * t374;
t325 = t371 * t378 + t372 * t375;
t389 = t325 * qJD(6);
t477 = g(3) * t376;
t479 = t410 * t379 - t477;
t341 = -qJD(2) + t480;
t478 = qJD(4) * (qJD(2) + t341 + t480) + qJDD(4) * t373;
t358 = 0.2e1 * t365;
t476 = g(3) * t379;
t475 = pkin(8) + qJ(5);
t447 = qJD(4) * t376;
t406 = -qJDD(4) * pkin(4) + t340 * t447 + qJDD(5);
t275 = -t327 * t379 + t406;
t473 = t275 * t371;
t472 = t275 * t372;
t471 = t275 * t379;
t469 = t321 * t375;
t276 = t378 * t319 + t469;
t470 = t276 * t346;
t441 = qJD(1) * qJD(4);
t420 = t379 * t441;
t436 = qJDD(1) * t376;
t390 = t420 + t436;
t323 = qJDD(6) + t390;
t468 = t323 * t324;
t467 = t323 * t325;
t465 = t371 * t379;
t464 = t372 * t379;
t463 = t373 * t376;
t330 = t376 * t340;
t462 = t376 * t377;
t461 = t376 * t380;
t331 = t407 + t374;
t298 = qJD(1) * t331 - qJD(2);
t314 = qJD(4) * qJ(5) + t330;
t266 = t371 * t298 + t372 * t314;
t428 = t371 * t452;
t459 = -t375 * t428 + t452 * t460 + t481;
t388 = t325 * qJD(1);
t458 = t376 * t388 + t389;
t329 = t408 * qJD(1);
t282 = t371 * t329 + t340 * t464;
t287 = t371 * t331 + t372 * t463;
t435 = qJDD(1) * t379;
t457 = t371 * qJDD(4) + t372 * t435;
t456 = t380 * pkin(1) + t377 * qJ(2);
t368 = t376 ^ 2;
t369 = t379 ^ 2;
t454 = t368 - t369;
t381 = qJD(4) ^ 2;
t382 = qJD(1) ^ 2;
t453 = -t381 - t382;
t449 = qJD(2) * t379;
t308 = -qJD(4) * pkin(4) + qJD(5) - t466;
t442 = -qJD(5) + t308;
t440 = qJD(3) * qJD(1);
t439 = qJDD(1) * t371;
t438 = qJDD(1) * t372;
t437 = qJDD(1) * t374;
t433 = qJDD(4) * t376;
t431 = pkin(8) * t372 * t376;
t412 = -qJDD(4) * t372 + t371 * t435;
t421 = t376 * t441;
t290 = t371 * t421 - t412;
t291 = t372 * t421 - t457;
t430 = t375 * t290 - t378 * t291 - t319 * t444;
t270 = t371 * t306 + t484 * t372;
t429 = t380 * qJ(3) + t456;
t425 = t371 * t447;
t419 = qJDD(2) - t482;
t416 = -t371 * t373 + pkin(5);
t415 = pkin(5) * t371 - t373;
t250 = pkin(5) * t390 + pkin(8) * t291 + t254;
t251 = pkin(8) * t290 + t255;
t414 = t378 * t250 - t251 * t375;
t413 = -t378 * t290 - t375 * t291;
t265 = t372 * t298 - t314 * t371;
t411 = -t370 + t419;
t281 = t372 * t329 - t340 * t465;
t404 = t250 * t375 + t251 * t378;
t402 = -t254 * t371 + t255 * t372;
t256 = pkin(5) * t452 - pkin(8) * t321 + t265;
t258 = -pkin(8) * t319 + t266;
t247 = t256 * t378 - t258 * t375;
t248 = t256 * t375 + t258 * t378;
t401 = t265 * t372 + t266 * t371;
t400 = t265 * t371 - t266 * t372;
t316 = t372 * t331;
t274 = -pkin(8) * t464 + t376 * t416 + t316;
t279 = -pkin(8) * t465 + t287;
t399 = t274 * t378 - t279 * t375;
t398 = t274 * t375 + t279 * t378;
t396 = -t364 + t411;
t395 = t410 * t371;
t338 = t475 * t372;
t394 = qJD(5) * t371 + qJD(6) * t338 + (pkin(5) * t379 + t431) * qJD(1) + t281;
t337 = t475 * t371;
t393 = pkin(8) * t428 - qJD(5) * t372 + qJD(6) * t337 + t282;
t392 = t358 + 0.2e1 * t366 - t410;
t252 = -t321 * t445 + t430;
t387 = qJD(1) * t341 + t410;
t386 = -t327 + t387;
t384 = -t410 * t376 - t476;
t253 = -t397 * qJD(6) + t413;
t328 = -t418 + t440;
t383 = -t373 * t381 + t328 + t437 + t440 + t482;
t363 = pkin(9) + qJ(6);
t357 = t380 * qJ(2);
t354 = qJDD(4) * t379;
t353 = cos(t363);
t352 = sin(t363);
t348 = -pkin(5) * t372 - pkin(4);
t318 = t415 * t379;
t305 = t324 * t379;
t304 = t325 * t379;
t302 = -t352 * t377 + t353 * t461;
t301 = -t352 * t461 - t353 * t377;
t300 = -t352 * t380 - t353 * t462;
t299 = t352 * t462 - t353 * t380;
t297 = -pkin(5) * t428 + t330;
t294 = -t415 * t447 - t449;
t293 = t372 * t306;
t286 = -t371 * t463 + t316;
t280 = pkin(5) * t319 + t308;
t269 = -t484 * t371 + t293;
t268 = -t375 * t376 * t443 - t378 * t425 + t481 * t379;
t267 = t324 * t447 - t379 * t389;
t260 = pkin(8) * t425 + t270;
t259 = -t371 * t450 + t293 + (t379 * t416 + t431) * qJD(4);
t257 = -pkin(5) * t290 + t275;
t1 = [(-t269 * t321 - t270 * t319 + t286 * t291 + t287 * t290 + t485 * t379 + t401 * t447) * MDP(19) + (-t252 * t304 + t253 * t305 - t267 * t276 + t268 * t397) * MDP(22) + (-t252 * t305 - t267 * t397) * MDP(21) + (-(t259 * t375 + t260 * t378) * t346 - t398 * t323 - t404 * t376 - t248 * t446 - t294 * t397 + t318 * t252 - t257 * t305 + t280 * t267 - g(1) * t299 - g(2) * t301 + (-t247 * t376 - t346 * t399) * qJD(6)) * MDP(27) + (t252 * t376 + t267 * t346 - t305 * t323 - t397 * t446) * MDP(23) + (t395 + (-qJD(2) * t319 + t473 + t373 * t290 + (qJD(1) * t286 + t265) * qJD(4)) * t379 + (t269 * qJD(1) + t286 * qJDD(1) + t254 + t482 * t372 + (-t308 * t371 + t319 * t373) * qJD(4)) * t376) * MDP(17) + (t483 + (-qJD(2) * t321 + t472 + t373 * t291 + (-qJD(1) * t287 - t266) * qJD(4)) * t379 + (-t270 * qJD(1) - t287 * qJDD(1) - t255 - t482 * t371 + (-t308 * t372 + t321 * t373) * qJD(4)) * t376) * MDP(18) + t482 * MDP(2) + (-t376 * t381 + t354) * MDP(12) + (-t432 * pkin(1) - g(1) * (-pkin(1) * t377 + t357) - g(2) * t456 + (t358 + t366) * qJ(2)) * MDP(6) + (t328 * t374 + t341 * qJD(3) + t417 * qJ(2) + t349 * qJD(2) - g(1) * (-t374 * t377 + t357) - g(2) * t429) * MDP(9) + (t255 * t287 + t266 * t270 + t254 * t286 + t265 * t269 - t373 * t471 - g(1) * (-pkin(7) * t380 + t357) - g(2) * (pkin(4) * t461 - t380 * t474 + t429) + (t373 * t447 - t449) * t308 + (g(2) * pkin(7) + g(1) * t331) * t377) * MDP(20) + 0.2e1 * (-t376 * t435 + t441 * t454) * MDP(11) + (-t253 * t376 - t268 * t346 - t276 * t446 - t304 * t323) * MDP(24) + ((t259 * t378 - t260 * t375) * t346 + t399 * t323 + t414 * t376 + t247 * t446 + t294 * t276 + t318 * t253 + t257 * t304 + t280 * t268 - g(1) * t300 - g(2) * t302 + (-t248 * t376 - t346 * t398) * qJD(6)) * MDP(26) + (t323 * t376 + t346 * t446) * MDP(25) + (-t396 + t437 + 0.2e1 * t440) * MDP(8) + (-t379 * t381 - t433) * MDP(13) + (qJDD(1) * t369 - 0.2e1 * t376 * t420) * MDP(10) + (-0.2e1 * t370 + t419) * MDP(4) + t410 * MDP(3) + (t383 * t376 + t478 * t379) * MDP(15) + (-t478 * t376 + t383 * t379) * MDP(16) + qJDD(1) * MDP(1) + (qJDD(3) + t392) * MDP(7) + t392 * MDP(5); t411 * MDP(6) + t396 * MDP(9) + (-t290 * t371 - t291 * t372) * MDP(19) + t485 * MDP(20) + (MDP(26) * t324 + MDP(27) * t325) * t323 + (MDP(26) * t458 + MDP(27) * t459) * t346 + (-qJ(2) * MDP(6) - MDP(5) - MDP(7) + (t371 * MDP(17) + t372 * MDP(18)) * t368) * t382 + (-t379 * MDP(16) + MDP(4) - MDP(8) + (-MDP(17) * t372 + MDP(18) * t371 - MDP(15)) * t376) * qJDD(1) + ((-qJD(3) - t349) * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + (t319 * t372 - t321 * t371) * MDP(19) + t400 * MDP(20)) * t376 + (-0.2e1 * qJD(4) * MDP(15) + (t319 - t443) * MDP(17) + (t321 + t448) * MDP(18) + t308 * MDP(20) + t276 * MDP(26) - t397 * MDP(27)) * t379) * qJD(1); qJDD(1) * MDP(7) - t382 * MDP(8) + (-t387 + t417) * MDP(9) + (t376 * t453 + t354) * MDP(15) + (t379 * t453 - t433) * MDP(16) + (-t368 * t439 + t290 * t379 + (-t372 * t382 + (t319 - 0.2e1 * t427) * qJD(4)) * t376) * MDP(17) + (-t368 * t438 + t291 * t379 + (t371 * t382 + (t321 - 0.2e1 * t426) * qJD(4)) * t376) * MDP(18) + ((qJD(1) * t321 + t290 * t376 - t319 * t446) * t372 + (qJD(1) * t319 - t291 * t376 + t321 * t446) * t371) * MDP(19) + (-t471 + t402 * t376 - t401 * qJD(1) + (t308 * t376 - t379 * t400) * qJD(4) - t410) * MDP(20) + (qJD(1) * t391 + (-qJD(4) * t325 * t346 - t253) * t379 + (qJD(4) * t276 + qJD(6) * t391 - t467) * t376) * MDP(26) + (t346 * t388 + (qJD(4) * t391 - t252) * t379 + (-qJD(4) * t397 + t346 * t389 + t468) * t376) * MDP(27); MDP(12) * t435 - MDP(13) * t436 + qJDD(4) * MDP(14) + (-t379 * t386 + t477) * MDP(15) + (t376 * t386 + t476) * MDP(16) + (pkin(4) * t290 - t472 + (-t483 + (-qJ(5) * t448 - t265) * qJD(1)) * t379 + (-qJ(5) * t439 + g(3) * t372 - t319 * t340 + (t371 * t442 - t281) * qJD(1)) * t376) * MDP(17) + (pkin(4) * t291 + t473 + (t395 + (-qJ(5) * t443 + t266) * qJD(1)) * t379 + (-qJ(5) * t438 - g(3) * t371 - t321 * t340 + (t372 * t442 + t282) * qJD(1)) * t376) * MDP(18) + (t281 * t321 + t282 * t319 + (qJ(5) * t290 - qJD(5) * t319 - t265 * t452 + t255) * t372 + (-qJ(5) * t291 + qJD(5) * t321 - t266 * t452 - t254) * t371 + t384) * MDP(19) + (-t308 * t330 - t265 * t281 - t266 * t282 - t400 * qJD(5) + (-t275 - t479) * pkin(4) + (t384 + t402) * qJ(5)) * MDP(20) + (t252 * t325 - t397 * t459) * MDP(21) + (-t252 * t324 - t253 * t325 - t276 * t459 + t397 * t458) * MDP(22) + (t346 * t459 + t397 * t451 + t467) * MDP(23) + (t276 * t451 - t346 * t458 - t468) * MDP(24) - t346 * MDP(25) * t451 + ((-t337 * t378 - t338 * t375) * t323 + t348 * t253 + t257 * t324 - t247 * t451 - t297 * t276 + (t375 * t393 - t378 * t394) * t346 + t458 * t280 - t479 * t353) * MDP(26) + (-(-t337 * t375 + t338 * t378) * t323 + t348 * t252 + t257 * t325 + t248 * t451 + t297 * t397 + (t375 * t394 + t378 * t393) * t346 + t459 * t280 + t479 * t352) * MDP(27) + (t379 * t376 * MDP(10) - t454 * MDP(11)) * t382; ((t321 - t448) * t452 + t412) * MDP(17) + ((-t319 - t443) * t452 + t457) * MDP(18) + (-t319 ^ 2 - t321 ^ 2) * MDP(19) + (-t477 + t265 * t321 + t266 * t319 + (-t327 + t410) * t379 + t406) * MDP(20) + (t253 - t486) * MDP(26) + (t252 - t470) * MDP(27); -t397 * t276 * MDP(21) + (-t276 ^ 2 + t397 ^ 2) * MDP(22) + (t430 + t470) * MDP(23) + (-t413 - t486) * MDP(24) + t323 * MDP(25) + (-g(1) * t301 + g(2) * t299 + t248 * t346 + t280 * t397 + t352 * t476 + t414) * MDP(26) + (g(1) * t302 - g(2) * t300 + t247 * t346 + t276 * t280 + t353 * t476 - t404) * MDP(27) + (-MDP(23) * t469 + MDP(24) * t397 - MDP(26) * t248 - MDP(27) * t247) * qJD(6);];
tau  = t1;
