% Calculate Coriolis joint torque vector for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:59
% EndTime: 2021-01-15 22:03:15
% DurationCPUTime: 4.64s
% Computational Cost: add. (4225->424), mult. (13014->591), div. (0->0), fcn. (10383->10), ass. (0->185)
t368 = sin(pkin(10));
t370 = cos(pkin(5));
t373 = sin(qJ(2));
t469 = pkin(1) * t373;
t428 = t370 * t469;
t369 = sin(pkin(5));
t376 = cos(qJ(2));
t451 = t369 * t376;
t467 = pkin(7) + qJ(3);
t339 = t451 * t467 + t428;
t437 = qJD(3) * t373;
t378 = -t339 * qJD(2) - t369 * t437;
t478 = t368 * t378;
t466 = cos(pkin(10));
t417 = t466 * t376;
t401 = t369 * t417;
t355 = qJD(1) * t401;
t441 = qJD(1) * t369;
t425 = t373 * t441;
t341 = t368 * t425 - t355;
t338 = qJD(4) + t341;
t389 = t368 * t376 + t373 * t466;
t344 = t389 * t441;
t375 = cos(qJ(4));
t440 = qJD(1) * t370;
t406 = qJD(2) + t440;
t356 = t375 * t406;
t372 = sin(qJ(4));
t312 = t344 * t372 - t356;
t311 = qJD(5) + t312;
t477 = MDP(5) * (t373 ^ 2 - t376 ^ 2);
t374 = cos(qJ(5));
t371 = sin(qJ(5));
t382 = -t344 * t375 - t372 * t406;
t456 = t382 * t371;
t278 = -t338 * t374 - t456;
t476 = t278 * t338;
t475 = t373 * MDP(4);
t454 = t341 * t375;
t305 = t344 * t371 - t374 * t454;
t433 = qJD(5) * t371;
t434 = qJD(4) * t375;
t383 = -t372 * t433 + t374 * t434 - t305;
t474 = t383 * t311;
t421 = t467 * t373;
t404 = t369 * t421;
t468 = pkin(1) * t376;
t328 = (pkin(2) + t468) * t370 - t404;
t297 = t328 * t368 + t339 * t466;
t288 = pkin(8) * t370 + t297;
t452 = t368 * t373;
t347 = t369 * t452 - t401;
t348 = t389 * t369;
t351 = (-pkin(2) * t376 - pkin(1)) * t369;
t303 = pkin(3) * t347 - pkin(8) * t348 + t351;
t473 = t288 * t375 + t372 * t303;
t365 = t369 ^ 2;
t388 = -pkin(7) * t451 - t428;
t472 = -t365 * t469 + t388 * t370;
t343 = qJD(2) * t348;
t336 = qJD(1) * t343;
t430 = t370 * t468;
t360 = qJD(1) * t430;
t315 = qJD(2) * pkin(2) + t360 + (t370 * pkin(2) - t404) * qJD(1);
t330 = t339 * qJD(1);
t418 = t466 * t330;
t274 = t315 * t368 + t418;
t267 = pkin(8) * t406 + t274;
t442 = qJD(1) * t351;
t349 = qJD(3) + t442;
t285 = pkin(3) * t341 - pkin(8) * t344 + t349;
t253 = t267 * t375 + t285 * t372;
t358 = t360 * qJD(2);
t380 = (-qJD(2) * t421 + qJD(3) * t376) * t369;
t308 = qJD(1) * t380 + t358;
t419 = t466 * t308;
t264 = qJD(1) * t478 + t419;
t439 = qJD(2) * t373;
t424 = t369 * t439;
t402 = qJD(1) * t424;
t337 = qJD(2) * t355 - t368 * t402;
t357 = pkin(2) * t402;
t286 = pkin(3) * t336 - pkin(8) * t337 + t357;
t413 = t264 * t372 - t286 * t375;
t470 = -qJD(4) * t253 - t413;
t241 = -pkin(4) * t336 - t470;
t471 = t311 * (-pkin(4) * t382 + pkin(9) * t311) + t241;
t277 = -qJD(4) * t382 + t372 * t337;
t377 = qJD(1) ^ 2;
t435 = qJD(4) * t372;
t276 = qJD(4) * t356 + t337 * t375 - t344 * t435;
t432 = qJD(5) * t374;
t426 = t276 * t374 + t336 * t371 + t338 * t432;
t249 = t382 * t433 + t426;
t465 = t249 * t371;
t280 = t338 * t371 - t374 * t382;
t412 = t276 * t371 - t336 * t374;
t250 = qJD(5) * t280 + t412;
t464 = t250 * t375;
t463 = t278 * t311;
t462 = t280 * t311;
t461 = t280 * t341;
t410 = t311 * t374;
t460 = t312 * t338;
t459 = t312 * t344;
t458 = t382 * t338;
t457 = t382 * t344;
t455 = t338 * t372;
t453 = t365 * t377;
t318 = t368 * t330;
t450 = t371 * t277;
t449 = t372 * t336;
t447 = t374 * t277;
t329 = -qJD(1) * t404 + t360;
t289 = t329 * t368 + t418;
t446 = t289 - t338 * (pkin(4) * t372 - pkin(9) * t375);
t290 = t329 * t466 - t318;
t408 = pkin(2) * t425;
t299 = pkin(3) * t344 + pkin(8) * t341 + t408;
t444 = t290 * t375 + t299 * t372;
t438 = qJD(2) * t376;
t252 = -t267 * t372 + t285 * t375;
t247 = -pkin(4) * t338 - t252;
t436 = qJD(4) * t247;
t431 = qJD(1) * qJD(2);
t427 = t376 * t453;
t420 = t365 * t431;
t387 = t264 * t375 - t267 * t435 + t285 * t434 + t286 * t372;
t240 = pkin(9) * t336 + t387;
t310 = t466 * t378;
t263 = -qJD(1) * t310 + t368 * t308;
t245 = pkin(4) * t277 - pkin(9) * t276 + t263;
t415 = -t240 * t371 + t245 * t374;
t414 = -t249 * t375 + t280 * t435;
t361 = qJD(2) * t430;
t316 = t361 + t380;
t268 = t316 * t368 - t310;
t411 = t338 * t375;
t364 = -pkin(2) * t466 - pkin(3);
t350 = -pkin(4) * t375 - pkin(9) * t372 + t364;
t409 = pkin(9) * t344 - qJD(5) * t350 + t444;
t407 = pkin(2) * t424;
t403 = t376 * t420;
t273 = t315 * t466 - t318;
t296 = t328 * t466 - t339 * t368;
t399 = t240 * t374 + t245 * t371;
t248 = pkin(9) * t338 + t253;
t266 = -pkin(3) * t406 - t273;
t251 = pkin(4) * t312 + pkin(9) * t382 + t266;
t239 = t248 * t374 + t251 * t371;
t398 = t248 * t371 - t251 * t374;
t259 = pkin(9) * t347 + t473;
t287 = -pkin(3) * t370 - t296;
t321 = t348 * t372 - t370 * t375;
t322 = t348 * t375 + t370 * t372;
t260 = pkin(4) * t321 - pkin(9) * t322 + t287;
t397 = t259 * t374 + t260 * t371;
t396 = -t259 * t371 + t260 * t374;
t269 = t316 * t466 + t478;
t346 = (t417 - t452) * t369 * qJD(2);
t300 = pkin(3) * t343 - pkin(8) * t346 + t407;
t395 = -t269 * t372 + t300 * t375;
t393 = -t288 * t372 + t303 * t375;
t302 = t322 * t374 + t347 * t371;
t301 = t322 * t371 - t347 * t374;
t392 = t336 * t375 - t338 * t435 - t341 * t455;
t304 = -t344 * t374 - t371 * t454;
t391 = (-t371 * t434 + t304) * t311;
t390 = -t311 * t432 - t450;
t386 = t269 * t375 - t288 * t435 + t300 * t372 + t303 * t434;
t363 = pkin(2) * t368 + pkin(8);
t385 = t266 * t338 - t336 * t363;
t381 = -pkin(9) * t277 + (t247 + t252) * t311;
t295 = qJD(4) * t322 + t346 * t372;
t294 = -qJD(4) * t321 + t346 * t375;
t258 = -pkin(4) * t347 - t393;
t257 = qJD(5) * t302 + t294 * t371 - t343 * t374;
t256 = -qJD(5) * t301 + t294 * t374 + t343 * t371;
t254 = -pkin(4) * t344 + t290 * t372 - t299 * t375;
t246 = pkin(4) * t295 - pkin(9) * t294 + t268;
t243 = -pkin(4) * t343 + qJD(4) * t473 - t395;
t242 = pkin(9) * t343 + t386;
t237 = -qJD(5) * t239 + t415;
t236 = -qJD(5) * t398 + t399;
t1 = [-0.2e1 * t420 * t477 + (t388 * qJD(2) ^ 2 + 0.2e1 * t431 * t472) * MDP(9) + (-0.2e1 * pkin(1) * t403 - (-pkin(7) * t424 + t361) * t406 - (-pkin(7) * t402 + t358) * t370) * MDP(10) + (-t268 * t406 - t263 * t370 + t351 * t336 + t349 * t343 + (qJD(1) * t347 + t341) * t407) * MDP(11) + (-t269 * t406 - t264 * t370 + t351 * t337 + t349 * t346 + (qJD(1) * t348 + t344) * t407) * MDP(12) + (t263 * t348 - t264 * t347 + t268 * t344 - t269 * t341 - t273 * t346 - t274 * t343 - t296 * t337 - t297 * t336) * MDP(13) + (-t263 * t296 + t264 * t297 - t268 * t273 + t269 * t274 + (t349 + t442) * t407) * MDP(14) + (t276 * t322 - t294 * t382) * MDP(15) + (-t276 * t321 - t277 * t322 - t294 * t312 + t295 * t382) * MDP(16) + (t276 * t347 + t294 * t338 + t322 * t336 - t343 * t382) * MDP(17) + (-t277 * t347 - t295 * t338 - t312 * t343 - t321 * t336) * MDP(18) + (t336 * t347 + t338 * t343) * MDP(19) + (t395 * t338 + t393 * t336 - t413 * t347 + t252 * t343 + t268 * t312 + t287 * t277 + t263 * t321 + t266 * t295 + (-t253 * t347 - t338 * t473) * qJD(4)) * MDP(20) + (-t253 * t343 + t263 * t322 + t266 * t294 - t268 * t382 + t287 * t276 - t336 * t473 - t338 * t386 - t347 * t387) * MDP(21) + (t249 * t302 + t256 * t280) * MDP(22) + (-t249 * t301 - t250 * t302 - t256 * t278 - t257 * t280) * MDP(23) + (t249 * t321 + t256 * t311 + t277 * t302 + t280 * t295) * MDP(24) + (-t250 * t321 - t257 * t311 - t277 * t301 - t278 * t295) * MDP(25) + (t277 * t321 + t295 * t311) * MDP(26) + ((-qJD(5) * t397 - t242 * t371 + t246 * t374) * t311 + t396 * t277 + t237 * t321 - t398 * t295 + t243 * t278 + t258 * t250 + t241 * t301 + t247 * t257) * MDP(27) + (-(qJD(5) * t396 + t242 * t374 + t246 * t371) * t311 - t397 * t277 - t236 * t321 - t239 * t295 + t243 * t280 + t258 * t249 + t241 * t302 + t247 * t256) * MDP(28) + 0.2e1 * t403 * t475 + (MDP(6) * t369 * t438 - MDP(7) * t424) * (qJD(2) + 0.2e1 * t440); -t338 * t344 * MDP(19) + t453 * t477 + (pkin(1) * t427 + (-pkin(7) * t425 + t360) * t440) * MDP(10) + (t289 * t406 - t341 * t408 - t349 * t344 - t263) * MDP(11) + (-t419 + t290 * qJD(2) + t349 * t341 + ((pkin(1) * t368 * t439 + t290) * t370 + (-t368 * (-t438 * t467 - t437) - t373 * pkin(2) * t344) * t369) * qJD(1)) * MDP(12) + ((t274 - t289) * t344 + (-t273 + t290) * t341 + (-t336 * t368 - t337 * t466) * pkin(2)) * MDP(13) + (t273 * t289 - t274 * t290 + (-t263 * t466 + t264 * t368 - t349 * t425) * pkin(2)) * MDP(14) + (t276 * t372 - t382 * t411) * MDP(15) + ((t276 - t460) * t375 + (-t277 + t458) * t372) * MDP(16) + (t338 * t411 + t449 + t457) * MDP(17) + (t392 + t459) * MDP(18) + (-t252 * t344 + t364 * t277 - t289 * t312 + (-t263 + (-qJD(4) * t363 - t299) * t338) * t375 + (t290 * t338 + t385) * t372) * MDP(20) + (t253 * t344 + t263 * t372 + t364 * t276 + t289 * t382 + (t363 * t435 + t444) * t338 + t385 * t375) * MDP(21) + (t249 * t372 * t374 + t280 * t383) * MDP(22) + (t278 * t305 + t280 * t304 + (-t278 * t374 - t280 * t371) * t434 + (-t465 - t250 * t374 + (t278 * t371 - t280 * t374) * qJD(5)) * t372) * MDP(23) + ((t447 + t461) * t372 + t474 + t414) * MDP(24) + (t464 + t391 + (t390 - t476) * t372) * MDP(25) + (-t277 * t375 + t311 * t455) * MDP(26) + (t350 * t447 - t247 * t304 - t254 * t278 + (t371 * t409 - t374 * t446) * t311 + (t371 * t436 - t237 + (qJD(4) * t278 + t390) * t363) * t375 + (t247 * t432 - t398 * t341 + t241 * t371 + t363 * t250 + (t311 * t363 * t371 - t398) * qJD(4)) * t372) * MDP(27) + (-t350 * t450 - t247 * t305 - t254 * t280 + (t371 * t446 + t374 * t409) * t311 + (t374 * t436 + t236 + (qJD(4) * t280 + t311 * t433 - t447) * t363) * t375 + (-t247 * t433 - t239 * t341 + t241 * t374 + t363 * t249 + (t363 * t410 - t239) * qJD(4)) * t372) * MDP(28) - t427 * t475 + ((-MDP(6) * t376 + MDP(7) * t373) * t369 * t370 - t472 * MDP(9)) * t377; (t344 * qJD(2) + (t344 * t370 + t343) * qJD(1)) * MDP(11) + (-t341 * t406 + t337) * MDP(12) + (-t341 ^ 2 - t344 ^ 2) * MDP(13) + (t273 * t344 + t274 * t341 + t357) * MDP(14) + (t392 - t459) * MDP(20) + (-t338 ^ 2 * t375 - t449 + t457) * MDP(21) + (-t464 + t391 + (t390 + t476) * t372) * MDP(27) + ((-t447 + t461) * t372 - t474 + t414) * MDP(28); -t312 ^ 2 * MDP(16) + (t276 + t460) * MDP(17) + (-t277 - t458) * MDP(18) + t336 * MDP(19) + (t253 * t338 + t470) * MDP(20) + (t252 * t338 + t266 * t312 - t387) * MDP(21) + (t280 * t410 + t465) * MDP(22) + ((t249 - t463) * t374 + (-t250 - t462) * t371) * MDP(23) + (t311 * t410 + t450) * MDP(24) + (-t311 ^ 2 * t371 + t447) * MDP(25) + (-pkin(4) * t250 - t253 * t278 + t381 * t371 - t374 * t471) * MDP(27) + (-pkin(4) * t249 - t253 * t280 + t371 * t471 + t381 * t374) * MDP(28) - (MDP(15) * t312 - MDP(16) * t382 - t266 * MDP(20) - MDP(24) * t280 + t278 * MDP(25) - t311 * MDP(26) + MDP(27) * t398 + MDP(28) * t239) * t382; t280 * t278 * MDP(22) + (-t278 ^ 2 + t280 ^ 2) * MDP(23) + (t426 + t463) * MDP(24) + (-t412 + t462) * MDP(25) + t277 * MDP(26) + (t239 * t311 - t247 * t280 + t415) * MDP(27) + (t247 * t278 - t311 * t398 - t399) * MDP(28) + (MDP(24) * t456 - MDP(25) * t280 - MDP(27) * t239 + MDP(28) * t398) * qJD(5);];
tauc = t1;
