% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:26
% EndTime: 2019-03-09 03:35:34
% DurationCPUTime: 3.66s
% Computational Cost: add. (4133->318), mult. (10037->436), div. (0->0), fcn. (7632->10), ass. (0->159)
t397 = cos(qJ(6));
t442 = qJD(6) * t397;
t390 = sin(pkin(11));
t392 = cos(pkin(11));
t396 = sin(qJ(3));
t399 = cos(qJ(3));
t373 = -t390 * t396 + t392 * t399;
t366 = t373 * qJD(1);
t398 = cos(qJ(5));
t354 = t398 * t366;
t374 = t390 * t399 + t392 * t396;
t368 = t374 * qJD(1);
t395 = sin(qJ(5));
t318 = -t368 * t395 + t354;
t489 = t318 * t397;
t493 = t442 - t489;
t416 = t366 * t395 + t398 * t368;
t394 = sin(qJ(6));
t443 = qJD(6) * t394;
t367 = t374 * qJD(3);
t357 = qJD(1) * t367;
t439 = qJD(1) * qJD(3);
t432 = t399 * t439;
t433 = t396 * t439;
t358 = -t390 * t433 + t392 * t432;
t444 = qJD(5) * t395;
t285 = qJD(5) * t354 - t395 * t357 + t398 * t358 - t368 * t444;
t387 = qJD(3) + qJD(5);
t451 = t397 * t285 + t387 * t442;
t271 = -t416 * t443 + t451;
t270 = t271 * t397;
t312 = t387 * t394 + t397 * t416;
t469 = t285 * t394;
t272 = t312 * qJD(6) + t469;
t461 = t416 * t394;
t310 = -t397 * t387 + t461;
t492 = -t394 * t272 - t493 * t310 + t270;
t269 = t271 * t394;
t286 = t416 * qJD(5) + t398 * t357 + t358 * t395;
t440 = -qJD(6) + t318;
t282 = t394 * t286;
t452 = -t440 * t442 + t282;
t463 = t318 * t387;
t465 = t416 * t387;
t467 = t312 * t416;
t491 = (-t286 + t465) * MDP(17) - t318 ^ 2 * MDP(15) + (-t318 * MDP(14) + MDP(15) * t416 + MDP(25) * t440) * t416 + (t285 - t463) * MDP(16) + (t493 * t312 + t269) * MDP(21) + (t440 * t489 + t452 - t467) * MDP(23);
t490 = t318 * t394;
t381 = sin(pkin(10)) * pkin(1) + pkin(7);
t457 = qJ(4) + t381;
t288 = pkin(5) * t416 - pkin(9) * t318;
t425 = t457 * qJD(1);
t342 = t399 * qJD(2) - t425 * t396;
t438 = qJD(1) * qJD(4);
t325 = t342 * qJD(3) + t399 * t438;
t343 = qJD(2) * t396 + t425 * t399;
t326 = -t343 * qJD(3) - t396 * t438;
t292 = -t325 * t390 + t392 * t326;
t279 = -pkin(8) * t358 + t292;
t293 = t392 * t325 + t390 * t326;
t280 = -pkin(8) * t357 + t293;
t333 = t390 * t343;
t472 = qJD(3) * pkin(3);
t337 = t342 + t472;
t300 = t392 * t337 - t333;
t473 = pkin(8) * t368;
t290 = qJD(3) * pkin(4) + t300 - t473;
t460 = t392 * t343;
t301 = t390 * t337 + t460;
t474 = pkin(8) * t366;
t291 = t301 + t474;
t250 = (qJD(5) * t290 + t280) * t398 + t279 * t395 - t291 * t444;
t383 = -cos(pkin(10)) * pkin(1) - pkin(2);
t413 = -pkin(3) * t399 + t383;
t406 = t413 * qJD(1);
t364 = qJD(4) + t406;
t327 = -pkin(4) * t366 + t364;
t487 = -t318 * t327 - t250;
t483 = MDP(5) * t396;
t468 = t310 * t416;
t481 = (t396 ^ 2 - t399 ^ 2) * MDP(6);
t284 = t397 * t286;
t480 = -t440 * t443 - t284;
t265 = t290 * t395 + t291 * t398;
t251 = t265 * qJD(5) - t398 * t279 + t280 * t395;
t264 = t290 * t398 - t291 * t395;
t262 = -pkin(5) * t387 - t264;
t263 = pkin(9) * t387 + t265;
t274 = -pkin(5) * t318 - pkin(9) * t416 + t327;
t419 = t263 * t394 - t274 * t397;
t479 = -t251 * t397 + t262 * t443 + t416 * t419;
t253 = t263 * t397 + t274 * t394;
t478 = t251 * t394 + t253 * t416 + t262 * t442;
t477 = -t327 * t416 - t251;
t427 = qJD(3) * t457;
t346 = qJD(4) * t399 - t396 * t427;
t347 = -qJD(4) * t396 - t399 * t427;
t306 = -t346 * t390 + t392 * t347;
t370 = t373 * qJD(3);
t298 = -pkin(8) * t370 + t306;
t307 = t392 * t346 + t390 * t347;
t299 = -pkin(8) * t367 + t307;
t371 = t457 * t396;
t372 = t457 * t399;
t323 = -t392 * t371 - t372 * t390;
t308 = -pkin(8) * t374 + t323;
t324 = -t390 * t371 + t392 * t372;
t309 = pkin(8) * t373 + t324;
t417 = t308 * t398 - t309 * t395;
t254 = t417 * qJD(5) + t298 * t395 + t299 * t398;
t276 = t308 * t395 + t309 * t398;
t329 = t373 * t395 + t374 * t398;
t341 = -pkin(4) * t373 + t413;
t415 = t398 * t373 - t374 * t395;
t281 = -pkin(5) * t415 - pkin(9) * t329 + t341;
t294 = t415 * qJD(5) - t367 * t395 + t370 * t398;
t476 = t251 * t329 + t262 * t294 - t276 * t286 + (qJD(6) * t281 + t254) * t440 + (qJD(6) * t274 + t250) * t415;
t475 = pkin(3) * t390;
t471 = t262 * t329;
t470 = t281 * t286;
t466 = t312 * t394;
t400 = qJD(3) ^ 2;
t459 = t396 * t400;
t458 = t399 * t400;
t295 = t329 * qJD(5) + t398 * t367 + t370 * t395;
t455 = -t271 * t415 + t312 * t295;
t302 = -t342 * t390 - t460;
t296 = t302 - t474;
t303 = t392 * t342 - t333;
t297 = t303 - t473;
t382 = pkin(3) * t392 + pkin(4);
t409 = t382 * t398 - t395 * t475;
t454 = -t409 * qJD(5) + t296 * t395 + t297 * t398;
t410 = t382 * t395 + t398 * t475;
t450 = t410 * qJD(5) + t296 * t398 - t297 * t395;
t448 = MDP(11) * t399;
t376 = qJD(1) * t383;
t446 = qJD(1) * t396;
t385 = t396 * t472;
t437 = t329 * t282;
t436 = t329 * t284;
t379 = pkin(3) * t433;
t332 = pkin(4) * t357 + t379;
t345 = pkin(4) * t367 + t385;
t344 = pkin(3) * t446 + pkin(4) * t368;
t426 = t440 * t394;
t363 = pkin(9) + t410;
t421 = qJD(6) * t363 + t288 + t344;
t420 = -t262 * t318 - t286 * t363;
t418 = t272 * t415 - t295 * t310;
t412 = 0.2e1 * qJD(3) * t376;
t411 = -t440 * t490 - t480;
t408 = -t294 * t394 - t329 * t442;
t407 = -t294 * t397 + t329 * t443;
t362 = -pkin(5) - t409;
t258 = pkin(5) * t295 - pkin(9) * t294 + t345;
t257 = pkin(5) * t286 - pkin(9) * t285 + t332;
t256 = t397 * t257;
t255 = t276 * qJD(5) - t298 * t398 + t299 * t395;
t1 = [0.2e1 * t432 * t483 - 0.2e1 * t439 * t481 + MDP(7) * t458 - MDP(8) * t459 + (-t381 * t458 + t396 * t412) * MDP(10) + (t381 * t459 + t399 * t412) * MDP(11) + (-t292 * t374 + t293 * t373 - t300 * t370 - t301 * t367 - t306 * t368 + t307 * t366 - t323 * t358 - t324 * t357) * MDP(12) + (t292 * t323 + t293 * t324 + t300 * t306 + t301 * t307 + (t364 + t406) * t385) * MDP(13) + (t285 * t329 + t294 * t416) * MDP(14) + (t285 * t415 - t286 * t329 + t294 * t318 - t295 * t416) * MDP(15) + (t286 * t341 + t295 * t327 - t318 * t345 - t332 * t415) * MDP(19) + (t285 * t341 + t294 * t327 + t329 * t332 + t345 * t416) * MDP(20) + (t329 * t270 - t407 * t312) * MDP(21) + ((-t310 * t397 - t466) * t294 + (-t269 - t272 * t397 + (t310 * t394 - t312 * t397) * qJD(6)) * t329) * MDP(22) + (t407 * t440 + t436 + t455) * MDP(23) + (-t408 * t440 + t418 - t437) * MDP(24) + (-t286 * t415 - t295 * t440) * MDP(25) + (-t419 * t295 + t255 * t310 - t256 * t415 - t417 * t272 + (-t258 * t440 + t470 + (t263 * t415 + t276 * t440 + t471) * qJD(6)) * t397 + t476 * t394) * MDP(26) + (-t253 * t295 + t255 * t312 - t417 * t271 + ((-qJD(6) * t276 + t258) * t440 - t470 + (-qJD(6) * t263 + t257) * t415 - qJD(6) * t471) * t394 + t476 * t397) * MDP(27) + (t294 * MDP(16) - t295 * MDP(17) - t255 * MDP(19) - t254 * MDP(20)) * t387; (-t357 * t374 - t358 * t373 + t366 * t370 + t367 * t368) * MDP(12) + (t292 * t373 + t293 * t374 - t300 * t367 + t301 * t370) * MDP(13) + (-t418 - t437) * MDP(26) + (-t436 + t455) * MDP(27) + (-MDP(10) * t396 - t448) * t400 + (-MDP(19) * t295 - MDP(20) * t294) * t387 - (t408 * MDP(26) + t407 * MDP(27)) * t440; ((t301 + t302) * t368 + (t300 - t303) * t366 + (-t357 * t390 - t358 * t392) * pkin(3)) * MDP(12) + (-t300 * t302 - t301 * t303 + (t292 * t392 + t293 * t390 - t364 * t446) * pkin(3)) * MDP(13) + (t318 * t344 - t450 * t387 + t477) * MDP(19) + (-t344 * t416 + t454 * t387 + t487) * MDP(20) + (t440 * t466 + t492) * MDP(22) + (t411 + t468) * MDP(24) + (t362 * t272 + t420 * t394 + t450 * t310 - (t454 * t394 - t421 * t397) * t440 + t479) * MDP(26) + (t362 * t271 + t420 * t397 + t450 * t312 - (t421 * t394 + t454 * t397) * t440 + t478) * MDP(27) + (-t399 * t483 + t481) * qJD(1) ^ 2 + (-MDP(10) * t446 - qJD(1) * t448) * t376 + t491; (-t366 ^ 2 - t368 ^ 2) * MDP(12) + (t300 * t368 - t301 * t366 + t379) * MDP(13) + (t286 + t465) * MDP(19) + (t285 + t463) * MDP(20) + (t411 - t468) * MDP(26) + (-t397 * t440 ^ 2 - t282 - t467) * MDP(27); (t265 * t387 + t477) * MDP(19) + (t264 * t387 + t487) * MDP(20) + (t312 * t426 + t492) * MDP(22) + (-t426 * t440 + t284 + t468) * MDP(24) + (-pkin(5) * t272 + (-t264 * t394 + t288 * t397) * t440 - t265 * t310 - t262 * t490 - t452 * pkin(9) + t479) * MDP(26) + (-pkin(5) * t271 - (t264 * t397 + t288 * t394) * t440 - t265 * t312 - t262 * t489 + t480 * pkin(9) + t478) * MDP(27) + t491; t312 * t310 * MDP(21) + (-t310 ^ 2 + t312 ^ 2) * MDP(22) + (-t310 * t440 + t451) * MDP(23) + (-t312 * t440 - t469) * MDP(24) + t286 * MDP(25) + (-t250 * t394 - t253 * t440 - t262 * t312 + t256) * MDP(26) + (-t250 * t397 - t257 * t394 + t262 * t310 + t419 * t440) * MDP(27) + (-MDP(23) * t461 - t312 * MDP(24) - t253 * MDP(26) + t419 * MDP(27)) * qJD(6);];
tauc  = t1;
