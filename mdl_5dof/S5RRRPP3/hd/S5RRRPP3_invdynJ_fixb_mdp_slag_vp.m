% Calculate vector of inverse dynamics joint torques for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:49
% EndTime: 2019-12-31 20:53:52
% DurationCPUTime: 2.81s
% Computational Cost: add. (1713->362), mult. (2394->393), div. (0->0), fcn. (1163->8), ass. (0->175)
t338 = qJDD(1) + qJDD(2);
t346 = sin(qJ(3));
t329 = t346 * qJ(4);
t349 = cos(qJ(3));
t416 = t349 * pkin(3) + t329;
t454 = -pkin(2) - t416;
t456 = t338 * t454;
t339 = qJD(1) + qJD(2);
t455 = t339 * t454;
t453 = pkin(2) + t329;
t447 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t437 = pkin(3) + qJ(5);
t452 = t437 * qJD(3);
t451 = t437 * qJDD(3);
t424 = t339 * t349;
t347 = sin(qJ(2));
t404 = qJDD(1) * t347;
t350 = cos(qJ(2));
t411 = qJD(2) * t350;
t264 = pkin(7) * t338 + (qJD(1) * t411 + t404) * pkin(1);
t258 = t349 * t264;
t436 = pkin(1) * qJD(1);
t402 = t347 * t436;
t284 = pkin(7) * t339 + t402;
t408 = qJD(3) * t346;
t236 = t284 * t408 - t258 - t447;
t257 = t346 * t264;
t407 = qJD(3) * t349;
t268 = t284 * t407;
t394 = qJDD(4) + t257 + t268;
t434 = qJDD(3) * pkin(3);
t237 = t394 - t434;
t450 = -t236 * t349 + t237 * t346;
t395 = t339 * t407;
t427 = t338 * t346;
t449 = (-t395 - t427) * pkin(4);
t344 = qJ(1) + qJ(2);
t327 = sin(t344);
t328 = cos(t344);
t418 = g(1) * t328 + g(2) * t327;
t448 = 0.2e1 * t447;
t273 = t349 * t284;
t255 = pkin(4) * t424 + t273;
t410 = qJD(3) * qJ(4);
t388 = -qJD(5) - t410;
t246 = -t388 + t255;
t272 = t346 * t284;
t446 = -qJD(4) - t272;
t429 = t328 * t346;
t431 = t327 * t346;
t445 = -g(1) * t429 - g(2) * t431 + g(3) * t349;
t364 = -t437 * t349 - t453;
t413 = qJD(1) * t350;
t401 = pkin(1) * t413;
t241 = t364 * t339 - t401;
t444 = t241 * t424 - t349 * t418;
t443 = pkin(1) * t350;
t442 = pkin(2) * t338;
t441 = pkin(2) * t339;
t352 = qJD(3) ^ 2;
t440 = pkin(7) * t352;
t317 = g(1) * t327;
t348 = sin(qJ(1));
t439 = g(1) * t348;
t438 = g(3) * t346;
t435 = pkin(7) * qJDD(3);
t320 = pkin(1) * t347 + pkin(7);
t432 = t320 * t352;
t430 = t327 * t349;
t428 = t328 * t349;
t426 = t338 * t349;
t425 = t339 * t346;
t423 = t346 * t349;
t422 = -g(1) * t431 + g(2) * t429;
t421 = g(1) * t430 - g(2) * t428;
t420 = qJ(4) * t426 + qJD(4) * t424;
t313 = t328 * pkin(7);
t419 = t328 * pkin(4) + t313;
t412 = qJD(2) * t347;
t326 = pkin(1) * t412;
t417 = -qJD(1) * t326 + qJDD(1) * t443;
t342 = t346 ^ 2;
t343 = t349 ^ 2;
t415 = t342 - t343;
t414 = t342 + t343;
t409 = qJD(3) * t339;
t406 = qJD(4) * t346;
t254 = -pkin(4) * t425 - t272;
t405 = qJD(4) - t254;
t403 = qJDD(3) * t320;
t400 = pkin(1) * t411;
t324 = pkin(3) * t408;
t337 = t339 ^ 2;
t399 = t337 * t423;
t389 = t339 * t402;
t282 = t346 * t389;
t398 = t282 - t422;
t397 = t339 * t324 - t417;
t321 = -pkin(2) - t443;
t396 = t339 * t412;
t393 = pkin(4) * t426 + qJDD(5) + t258;
t392 = -pkin(4) * t339 - t284;
t391 = t414 * t338;
t248 = -t401 + t455;
t390 = -t248 - t455;
t387 = t257 + t445;
t263 = -t417 - t442;
t285 = -t401 - t441;
t386 = t263 * t346 + t285 * t407 + t422;
t385 = pkin(3) * t428 + t327 * pkin(7) + t453 * t328;
t384 = -t349 * t389 - t401 * t408 - t421;
t383 = t440 - t442;
t381 = qJ(5) * t349 + t416;
t379 = t392 * qJD(3);
t378 = 0.2e1 * (t338 * t423 - t415 * t409) * MDP(8) + (t338 * t342 + 0.2e1 * t346 * t395) * MDP(7) + (qJDD(3) * t349 - t346 * t352) * MDP(10) + (qJDD(3) * t346 + t349 * t352) * MDP(9) + t338 * MDP(4);
t377 = -qJ(4) * t349 + qJ(5) * t346;
t261 = -qJD(3) * pkin(3) - t446;
t265 = -t273 - t410;
t376 = t261 * t349 + t265 * t346;
t375 = t261 * t346 - t265 * t349;
t374 = t393 + t447;
t373 = -qJDD(4) - t387;
t271 = -pkin(2) - t381;
t372 = -g(2) * t328 + t317 + t417;
t371 = t414 * t401;
t369 = t327 * pkin(4) + qJ(5) * t428 + t385;
t228 = (t377 * qJD(3) - qJD(5) * t349 - t406) * t339 + t364 * t338 + t397;
t247 = qJ(5) * t408 + t388 * t349 + t324 - t406;
t243 = t326 + t247;
t262 = t271 - t443;
t368 = -t243 * t339 - t262 * t338 - t228;
t367 = -t247 * t339 - t271 * t338 - t228;
t366 = -t268 + t373;
t365 = -qJ(4) * t407 - t406;
t266 = t324 + t365;
t363 = t266 * t339 + t440 + t456;
t362 = t261 * t407 + t265 * t408 - t418 + t450;
t256 = t266 + t326;
t274 = t321 - t416;
t361 = t256 * t339 + t274 * t338 + t432;
t358 = t394 - t449 - t451;
t231 = -qJD(3) * qJD(5) + t358;
t233 = t346 * t379 + t374;
t244 = t405 - t452;
t360 = t231 * t346 + t233 * t349 + t244 * t407 - t246 * t408 - t418;
t359 = t364 * t317;
t357 = pkin(1) * t396 + t321 * t338 + t432;
t356 = -g(1) * t313 - t317 * t454;
t355 = -t403 + (t321 * t339 - t400) * qJD(3);
t354 = t403 + (-t274 * t339 - t248 + t400) * qJD(3);
t353 = t376 * qJD(3) + t450;
t351 = cos(qJ(1));
t335 = t351 * pkin(1);
t333 = t349 * pkin(4);
t332 = t346 * pkin(4);
t325 = pkin(4) * t407;
t309 = pkin(3) * t425;
t296 = qJ(4) * t428;
t293 = qJ(4) * t430;
t290 = pkin(7) * t349 + t333;
t289 = pkin(7) * t346 + t332;
t280 = pkin(7) * t407 + t325;
t279 = (-pkin(4) - pkin(7)) * t408;
t276 = t320 * t349 + t333;
t275 = t320 * t346 + t332;
t269 = t285 * t408;
t267 = -qJ(4) * t424 + t309;
t251 = t377 * t339 + t309;
t250 = t320 * t407 + t346 * t400 + t325;
t249 = t349 * t400 + (-pkin(4) - t320) * t408;
t242 = t248 * t425;
t239 = t241 * t408;
t234 = t365 * t339 + t397 + t456;
t232 = t234 * t349;
t1 = [(t234 * t274 + t248 * t256 - g(2) * (t335 + t385) + (t375 * t411 + t439) * pkin(1) + t353 * t320 + t356) * MDP(17) + ((t338 * t350 - t396) * pkin(1) + t372) * MDP(5) + (qJDD(3) * t276 + t368 * t346 + (t249 + (-t262 * t339 - t241) * t349) * qJD(3) - t422) * MDP(19) + (-qJDD(3) * t275 + t239 + (t262 * t425 - t250) * qJD(3) + t368 * t349 + t421) * MDP(20) + (t228 * t262 + t241 * t243 + t231 * t275 + t244 * t250 + t233 * t276 + t246 * t249 - g(1) * (-pkin(1) * t348 + t419) - g(2) * (t335 + t369) - t359) * MDP(21) + (t354 * t349 + (-t234 - t361) * t346 - t422) * MDP(16) + (t354 * t346 + t361 * t349 + t232 - t421) * MDP(15) + (t414 * t339 * t400 + t320 * t391 + t362) * MDP(14) + ((t275 * t346 + t276 * t349) * t338 + (t249 * t349 + t250 * t346 + (t275 * t349 - t276 * t346) * qJD(3)) * t339 + t360) * MDP(18) + (t269 + t355 * t346 + (-t263 - t357) * t349 + t421) * MDP(12) + (t357 * t346 + t355 * t349 + t386) * MDP(13) + (((-qJDD(1) - t338) * t347 + (-qJD(1) - t339) * t411) * pkin(1) + t418) * MDP(6) + t378 + qJDD(1) * MDP(1) + (-g(2) * t351 + t439) * MDP(2) + (g(1) * t351 + g(2) * t348) * MDP(3); (t372 + t389) * MDP(5) + ((-t404 + (-qJD(2) + t339) * t413) * pkin(1) + t418) * MDP(6) + (t269 + (-pkin(2) * t409 - t435) * t346 + (-t263 - t383) * t349 - t384) * MDP(12) + (-t282 + t383 * t346 + (-t435 + (t401 - t441) * qJD(3)) * t349 + t386) * MDP(13) + (pkin(7) * t391 - t339 * t371 + t362) * MDP(14) + (t232 + t363 * t349 + (t390 * qJD(3) + t435) * t346 + t384) * MDP(15) + ((t435 + (t390 - t401) * qJD(3)) * t349 + (-t234 - t363) * t346 + t398) * MDP(16) + (t234 * t454 + t248 * t266 - g(2) * t385 + (-t248 * t347 - t375 * t350) * t436 + t353 * pkin(7) + t356) * MDP(17) + ((t289 * t346 + t290 * t349) * t338 + (t279 * t349 + t280 * t346 + (t289 * t349 - t290 * t346) * qJD(3) - t371) * t339 + t360) * MDP(18) + (qJDD(3) * t290 + t367 * t346 + (t279 + (-t271 * t339 - t241 - t401) * t349) * qJD(3) + t398) * MDP(19) + (-qJDD(3) * t289 + t239 + (t271 * t425 - t280) * qJD(3) + t367 * t349 - t384) * MDP(20) + (t228 * t271 + t241 * t247 + t231 * t289 + t244 * t280 + t233 * t290 + t246 * t279 - g(1) * t419 - g(2) * t369 - t359 + (-t241 * t347 + (-t244 * t346 - t246 * t349) * t350) * t436) * MDP(21) + t378; -MDP(7) * t399 + t415 * MDP(8) * t337 + MDP(9) * t427 + MDP(10) * t426 + qJDD(3) * MDP(11) + (-t285 * t425 - t387) * MDP(12) + (t438 - t258 + (-t285 * t339 + t418) * t349) * MDP(13) + (-pkin(3) * t427 + (-qJD(3) * t416 - t376) * t339 + t420) * MDP(14) + (-t267 * t424 + t242 - t373 - 0.2e1 * t434) * MDP(15) + (t258 + (t267 * t339 - g(3)) * t346 + (t248 * t339 - t418) * t349 + t448) * MDP(16) + (-t236 * qJ(4) - t237 * pkin(3) - t248 * t267 - t261 * t273 - g(1) * (-pkin(3) * t429 + t296) - g(2) * (-pkin(3) * t431 + t293) - g(3) * t416 + t446 * t265) * MDP(17) + (-t437 * t427 + (-t244 - t254 - t452) * t424 + t420) * MDP(18) + (-qJD(3) * t254 + (t251 * t339 - g(3) + t379) * t346 + t393 + t444 + t448) * MDP(19) + ((-t241 * t346 + t251 * t349) * t339 + (0.2e1 * qJD(5) + t255) * qJD(3) + 0.2e1 * t451 + t366 + t449) * MDP(20) + (t233 * qJ(4) - t241 * t251 - g(1) * t296 - g(2) * t293 - g(3) * t381 + t405 * t246 + (-qJD(5) - t255) * t244 + (t418 * t346 - t231) * t437) * MDP(21); (qJD(3) * t265 + t242 - t366 - t434) * MDP(17) + (t241 * t425 + (-qJD(5) - t246) * qJD(3) + t358 + t445) * MDP(21) + (MDP(14) + MDP(18)) * t427 + (MDP(16) + MDP(19)) * (-t337 * t342 - t352) + (MDP(15) - MDP(20)) * (qJDD(3) + t399); qJDD(3) * MDP(19) + (-t337 * t343 - t352) * MDP(20) + (-t337 * t346 * MDP(19) + t338 * MDP(18)) * t349 + (t374 - t438 + (t392 * t346 + t244) * qJD(3) + t444) * MDP(21);];
tau = t1;
