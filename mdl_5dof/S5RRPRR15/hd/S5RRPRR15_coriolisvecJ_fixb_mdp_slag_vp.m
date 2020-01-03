% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR15_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:15
% EndTime: 2019-12-31 20:43:22
% DurationCPUTime: 3.46s
% Computational Cost: add. (1784->354), mult. (4150->497), div. (0->0), fcn. (2596->6), ass. (0->176)
t438 = pkin(3) + pkin(6);
t350 = cos(qJ(4));
t347 = sin(qJ(4));
t408 = qJD(2) * t347;
t351 = cos(qJ(2));
t409 = qJD(1) * t351;
t299 = t350 * t409 + t408;
t349 = cos(qJ(5));
t387 = t347 * t409;
t406 = qJD(2) * t350;
t301 = -t387 + t406;
t346 = sin(qJ(5));
t428 = t301 * t346;
t248 = t349 * t299 + t428;
t348 = sin(qJ(2));
t410 = qJD(1) * t348;
t333 = qJD(4) + t410;
t325 = qJD(5) + t333;
t448 = t248 * t325;
t369 = t299 * t346 - t349 * t301;
t447 = t325 * t369;
t440 = qJD(4) + qJD(5);
t352 = -pkin(2) - pkin(7);
t383 = -qJ(3) * t348 - pkin(1);
t295 = t352 * t351 + t383;
t269 = t295 * qJD(1);
t335 = pkin(6) * t410;
t441 = qJD(3) + t335;
t396 = pkin(3) * t410 + t441;
t274 = t352 * qJD(2) + t396;
t241 = t269 * t350 + t274 * t347;
t236 = -pkin(8) * t299 + t241;
t400 = qJD(5) * t346;
t234 = t236 * t400;
t336 = pkin(6) * t409;
t309 = pkin(3) * t409 + t336;
t343 = qJD(2) * qJ(3);
t288 = t343 + t309;
t255 = pkin(4) * t299 + t288;
t446 = t248 * t255 + t234;
t395 = qJD(1) * qJD(2);
t386 = t348 * t395;
t262 = -qJD(4) * t299 + t347 * t386;
t332 = pkin(2) * t386;
t371 = pkin(7) * t348 - qJ(3) * t351;
t404 = qJD(3) * t348;
t357 = t371 * qJD(2) - t404;
t257 = qJD(1) * t357 + t332;
t385 = t351 * t395;
t331 = pkin(6) * t385;
t294 = pkin(3) * t385 + t331;
t378 = -t257 * t347 + t350 * t294;
t356 = -t241 * qJD(4) + t378;
t227 = pkin(4) * t385 - pkin(8) * t262 + t356;
t402 = qJD(4) * t350;
t263 = qJD(2) * t402 - qJD(4) * t387 - t350 * t386;
t393 = -t350 * t257 - t274 * t402 - t347 * t294;
t403 = qJD(4) * t347;
t359 = -t269 * t403 - t393;
t228 = -pkin(8) * t263 + t359;
t379 = t349 * t227 - t346 * t228;
t445 = t255 * t369 + t379;
t384 = MDP(26) * t409;
t444 = qJD(2) * t384 + (-t248 ^ 2 + t369 ^ 2) * MDP(23) - t248 * t369 * MDP(22);
t443 = -0.2e1 * t395;
t344 = t348 ^ 2;
t345 = t351 ^ 2;
t442 = (t344 - t345) * MDP(5);
t318 = t438 * t348;
t304 = t347 * t318;
t413 = t350 * t295 + t304;
t377 = t262 * t346 + t349 * t263;
t232 = -t369 * qJD(5) + t377;
t439 = t440 * t351;
t437 = pkin(8) - t352;
t436 = qJD(2) * pkin(2);
t240 = -t269 * t347 + t350 * t274;
t235 = -pkin(8) * t301 + t240;
t233 = pkin(4) * t333 + t235;
t435 = t233 * t349;
t434 = t236 * t349;
t433 = t262 * t350;
t407 = qJD(2) * t348;
t308 = t438 * t407;
t342 = qJD(2) * qJD(3);
t276 = -qJD(1) * t308 + t342;
t432 = t276 * t347;
t431 = t276 * t350;
t430 = t299 * t333;
t429 = t301 * t333;
t427 = t301 * t351;
t426 = t333 * t350;
t425 = t333 * t352;
t424 = t346 * t347;
t423 = t347 * t348;
t353 = qJD(2) ^ 2;
t422 = t348 * t353;
t421 = t349 * t350;
t420 = t350 * t351;
t419 = t351 * t353;
t354 = qJD(1) ^ 2;
t418 = t351 * t354;
t302 = t346 * t350 + t347 * t349;
t253 = t440 * t302;
t368 = -t421 + t424;
t417 = -t253 * t325 - t368 * t385;
t360 = t302 * t348;
t271 = qJD(1) * t360;
t416 = -t253 - t271;
t390 = t350 * t410;
t415 = -t346 * t403 - t347 * t400 + t349 * t390 - t410 * t424 + t440 * t421;
t339 = pkin(2) * t410;
t279 = t371 * qJD(1) + t339;
t414 = t350 * t279 + t347 * t309;
t391 = -pkin(4) * t350 - pkin(3);
t412 = pkin(4) * t402 - t391 * t410 + t441;
t319 = t438 * t351;
t315 = -pkin(2) * t351 + t383;
t289 = qJD(1) * t315;
t405 = qJD(2) * t351;
t401 = qJD(4) * t351;
t399 = qJD(5) * t349;
t397 = t288 * qJD(4);
t394 = t348 * t418;
t392 = t349 * t262 - t346 * t263 - t299 * t399;
t389 = t347 * t401;
t388 = t350 * t401;
t314 = t437 * t350;
t382 = pkin(8) * t351 - t295;
t381 = pkin(1) * t443;
t380 = qJD(3) - t436;
t338 = pkin(2) * t407;
t265 = t338 + t357;
t310 = t438 * t405;
t376 = -t265 * t347 + t350 * t310;
t375 = -t279 * t347 + t350 * t309;
t374 = qJD(5) * t233 + t228;
t321 = t348 * t385;
t313 = t437 * t347;
t364 = pkin(4) * t351 - pkin(8) * t423;
t373 = qJD(1) * t364 - qJD(5) * t313 - t437 * t403 + t375;
t372 = pkin(8) * t390 + t440 * t314 + t414;
t225 = t233 * t346 + t434;
t305 = t350 * t318;
t244 = pkin(4) * t348 + t382 * t347 + t305;
t246 = -pkin(8) * t420 + t413;
t370 = t244 * t346 + t246 * t349;
t367 = -qJD(1) * t345 + t333 * t348;
t366 = -0.2e1 * qJD(2) * t289;
t365 = t333 * t347;
t361 = -qJ(3) * t405 - t404;
t267 = qJD(1) * t361 + t332;
t284 = t338 + t361;
t363 = pkin(6) * t353 + qJD(1) * t284 + t267;
t362 = t288 * t348 + t352 * t405;
t358 = t350 * t265 - t295 * t403 + t347 * t310 + t318 * t402;
t231 = -t301 * t400 + t392;
t311 = pkin(6) * t386 - t342;
t312 = t335 + t380;
t317 = -t336 - t343;
t355 = -t311 * t351 + (t312 * t351 + (t317 + t336) * t348) * qJD(2);
t334 = pkin(4) * t347 + qJ(3);
t322 = t350 * t385;
t306 = -qJ(3) * t409 + t339;
t287 = pkin(4) * t420 + t319;
t281 = t302 * t351;
t280 = t368 * t351;
t273 = t289 * t410;
t256 = -pkin(4) * t389 + (-pkin(6) + t391) * t407;
t243 = pkin(4) * t263 + t276;
t238 = t302 * t439 - t368 * t407;
t237 = qJD(2) * t360 + t368 * t439;
t230 = (t348 * t406 + t389) * pkin(8) + t358;
t229 = t364 * qJD(2) + (t382 * t350 - t304) * qJD(4) + t376;
t224 = -t236 * t346 + t435;
t1 = [MDP(6) * t419 - MDP(7) * t422 + (t348 * t366 + t351 * t363) * MDP(12) + (-t348 * t363 + t351 * t366) * MDP(13) + (pkin(6) * t355 + t267 * t315 + t284 * t289) * MDP(14) + (-t262 * t347 * t351 + (t347 * t407 - t388) * t301) * MDP(15) + ((-t299 * t347 + t301 * t350) * t407 + (-t433 + t263 * t347 + (t299 * t350 + t301 * t347) * qJD(4)) * t351) * MDP(16) + (-t333 * t388 + t262 * t348 + (t347 * t367 + t427) * qJD(2)) * MDP(17) + (t333 * t389 - t263 * t348 + (-t299 * t351 + t350 * t367) * qJD(2)) * MDP(18) + (t333 * t405 + t321) * MDP(19) + (t376 * t333 - t308 * t299 + t319 * t263 + (-t288 * t406 + t378) * t348 + (-t241 * t348 - t413 * t333) * qJD(4) + (-t347 * t397 + t431 + ((-t295 * t347 + t305) * qJD(1) + t240) * qJD(2)) * t351) * MDP(20) + (-t358 * t333 - t308 * t301 + t319 * t262 + ((qJD(2) * t288 + qJD(4) * t269) * t347 + t393) * t348 + (-t350 * t397 - t432 + (-t413 * qJD(1) - t241) * qJD(2)) * t351) * MDP(21) + (-t231 * t281 - t237 * t369) * MDP(22) + (t231 * t280 + t232 * t281 - t237 * t248 - t238 * t369) * MDP(23) + (t231 * t348 + t237 * t325 + (-qJD(1) * t281 - t369) * t405) * MDP(24) + (-t232 * t348 + t238 * t325 + (qJD(1) * t280 - t248) * t405) * MDP(25) + (t325 * t405 + t321) * MDP(26) + ((t229 * t349 - t230 * t346) * t325 + t379 * t348 + t256 * t248 + t287 * t232 - t243 * t280 - t255 * t238 + (-t225 * t348 - t370 * t325) * qJD(5) + ((t244 * t349 - t246 * t346) * qJD(1) + t224) * t405) * MDP(27) + (t287 * t231 + t234 * t348 + t255 * t237 - t243 * t281 - t256 * t369 + (-(-qJD(5) * t246 + t229) * t325 - t227 * t348) * t346 + (-(qJD(5) * t244 + t230) * t325 - t374 * t348) * t349 + (-t370 * qJD(1) - t225) * t405) * MDP(28) + t442 * t443 + (-pkin(6) * t419 + t348 * t381) * MDP(9) + (pkin(6) * t422 + t351 * t381) * MDP(10) + t355 * MDP(11) + 0.2e1 * MDP(4) * t321; -MDP(4) * t394 + t354 * t442 + ((-t317 - t343) * t348 + (-t312 + t380) * t351) * qJD(1) * MDP(11) + (-t306 * t409 + t273) * MDP(12) + (0.2e1 * t342 + (t289 * t351 + t306 * t348) * qJD(1)) * MDP(13) + (-qJ(3) * t311 - qJD(3) * t317 - t289 * t306 + (-t317 * t348 + (-t312 - t436) * t351) * qJD(1) * pkin(6)) * MDP(14) + (-t301 * t365 + t433) * MDP(15) + ((-t263 - t429) * t350 + (-t262 + t430) * t347) * MDP(16) + (-t333 * t403 + t322 + (-t333 * t423 - t427) * qJD(1)) * MDP(17) + (-t333 * t402 + (-t348 * t426 + (t299 - t408) * t351) * qJD(1)) * MDP(18) - t333 * MDP(19) * t409 + (qJ(3) * t263 + t432 - t375 * t333 + t396 * t299 + (t288 * t350 - t347 * t425) * qJD(4) + (-t240 * t351 + t350 * t362) * qJD(1)) * MDP(20) + (qJ(3) * t262 + t431 + t414 * t333 + t396 * t301 + (-t288 * t347 - t350 * t425) * qJD(4) + (t241 * t351 - t347 * t362) * qJD(1)) * MDP(21) + (-t231 * t368 - t369 * t416) * MDP(22) + (-t231 * t302 + t232 * t368 - t416 * t248 + t369 * t415) * MDP(23) + (-t271 * t325 + t369 * t409 + t417) * MDP(24) + (-t415 * t325 + (-qJD(2) * t302 + t248) * t409) * MDP(25) - t325 * t384 + (t334 * t232 + t243 * t302 + (t372 * t346 - t373 * t349) * t325 + t415 * t255 + t412 * t248 + ((t313 * t346 - t314 * t349) * qJD(2) - t224) * t409) * MDP(27) + (t334 * t231 - t243 * t368 + (t373 * t346 + t372 * t349) * t325 + t416 * t255 - t412 * t369 + (-(-t313 * t349 - t314 * t346) * qJD(2) + t225) * t409) * MDP(28) + (t354 * t348 * MDP(9) + MDP(10) * t418) * pkin(1); MDP(12) * t394 + (-t344 * t354 - t353) * MDP(13) + (t273 + t331) * MDP(14) + t322 * MDP(20) + t417 * MDP(27) + (-t271 * MDP(27) - t415 * MDP(28)) * t325 + (-MDP(20) * t365 - MDP(21) * t426) * t333 + (t317 * MDP(14) - t299 * MDP(20) + (-t301 - t387) * MDP(21) - t248 * MDP(27) + (-t302 * t409 + t369) * MDP(28)) * qJD(2); t301 * t299 * MDP(15) + (-t299 ^ 2 + t301 ^ 2) * MDP(16) + (t262 + t430) * MDP(17) + (-t263 + t429) * MDP(18) + MDP(19) * t385 + (t241 * t333 - t288 * t301 + t356) * MDP(20) + (t240 * t333 + t288 * t299 - t359) * MDP(21) + (t231 + t448) * MDP(24) + (-t232 - t447) * MDP(25) + (-(-t235 * t346 - t434) * t325 - t225 * qJD(5) + (-t248 * t301 - t325 * t400 + t349 * t385) * pkin(4) + t445) * MDP(27) + ((-t236 * t325 - t227) * t346 + (t235 * t325 - t374) * t349 + (t301 * t369 - t325 * t399 - t346 * t385) * pkin(4) + t446) * MDP(28) + t444; (t392 + t448) * MDP(24) + (-t377 - t447) * MDP(25) + (t225 * t325 + t445) * MDP(27) + (t224 * t325 - t346 * t227 - t349 * t228 + t446) * MDP(28) + (-MDP(24) * t428 + t369 * MDP(25) - t225 * MDP(27) - MDP(28) * t435) * qJD(5) + t444;];
tauc = t1;
