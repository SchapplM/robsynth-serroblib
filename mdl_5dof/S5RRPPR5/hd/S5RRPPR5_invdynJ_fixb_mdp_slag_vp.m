% Calculate vector of inverse dynamics joint torques for
% S5RRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:37
% EndTime: 2021-01-15 19:36:47
% DurationCPUTime: 4.14s
% Computational Cost: add. (1889->352), mult. (4278->430), div. (0->0), fcn. (3007->10), ass. (0->151)
t359 = sin(pkin(8));
t360 = cos(pkin(8));
t363 = sin(qJ(2));
t366 = cos(qJ(2));
t321 = t359 * t366 + t360 * t363;
t310 = t321 * qJD(2);
t411 = qJDD(1) * t366;
t338 = t360 * t411;
t412 = qJDD(1) * t363;
t281 = qJD(1) * t310 + t359 * t412 - t338;
t413 = qJD(1) * qJD(2);
t404 = t363 * t413;
t376 = qJDD(1) * t321 - t359 * t404;
t403 = t366 * t413;
t282 = t360 * t403 + t376;
t362 = sin(qJ(5));
t365 = cos(qJ(5));
t425 = t360 * t366;
t406 = qJD(1) * t425;
t418 = qJD(1) * t363;
t308 = t359 * t418 - t406;
t311 = t321 * qJD(1);
t458 = t365 * t308 - t311 * t362;
t234 = qJD(5) * t458 + t362 * t281 + t365 * t282;
t351 = qJDD(2) - qJDD(5);
t352 = qJD(2) - qJD(5);
t431 = t458 * t352;
t449 = t308 * t362 + t365 * t311;
t461 = -MDP(19) * t458 * t449 + (t449 ^ 2 - t458 ^ 2) * MDP(20) + (t234 + t431) * MDP(21) - t351 * MDP(23);
t432 = t449 * t352;
t442 = pkin(2) * t366;
t347 = pkin(1) + t442;
t328 = -qJD(1) * t347 + qJD(3);
t459 = -qJ(4) * t311 + t328;
t364 = sin(qJ(1));
t367 = cos(qJ(1));
t450 = g(1) * t364 - g(2) * t367;
t398 = g(1) * t367 + g(2) * t364;
t361 = -qJ(3) - pkin(6);
t400 = qJD(2) * t361;
t377 = -qJD(3) * t363 + t366 * t400;
t405 = t361 * t363;
t277 = qJDD(2) * pkin(2) + qJD(1) * t377 + qJDD(1) * t405;
t304 = qJD(3) * t366 + t363 * t400;
t330 = t361 * t366;
t285 = qJD(1) * t304 - qJDD(1) * t330;
t421 = -t360 * t277 + t359 * t285;
t402 = -qJDD(4) - t421;
t444 = -pkin(3) - pkin(4);
t232 = -pkin(7) * t282 + qJDD(2) * t444 - t402;
t240 = t359 * t277 + t360 * t285;
t354 = qJDD(2) * qJ(4);
t355 = qJD(2) * qJD(4);
t237 = t354 + t355 + t240;
t233 = pkin(7) * t281 + t237;
t244 = t308 * t444 - t459;
t353 = qJ(2) + pkin(8);
t348 = sin(t353);
t349 = cos(t353);
t302 = -t348 * t362 - t349 * t365;
t303 = t348 * t365 - t349 * t362;
t457 = -g(3) * t303 + t362 * t232 + t365 * t233 + t244 * t458 + t398 * t302;
t395 = t347 * qJDD(1);
t455 = g(3) * t302 - t365 * t232 + t362 * t233 + t244 * t449 + t398 * t303;
t325 = qJD(1) * t405;
t326 = qJD(1) * t330;
t428 = t359 * t326;
t287 = t360 * t325 + t428;
t414 = qJD(4) - t287;
t291 = -t360 * t330 + t359 * t405;
t447 = -qJDD(2) * t291 - t348 * t450;
t446 = g(3) * t348 + qJD(2) * t287 + t398 * t349 - t240;
t399 = -t365 * t281 + t282 * t362;
t235 = qJD(5) * t449 + t399;
t306 = t311 ^ 2;
t443 = pkin(2) * t363;
t441 = pkin(7) * t308;
t440 = pkin(7) * t311;
t436 = g(3) * t349;
t435 = g(3) * t366;
t257 = pkin(3) * t308 + t459;
t433 = t257 * t311;
t427 = t359 * t363;
t426 = t360 * t326;
t424 = t361 * t364;
t260 = t360 * t304 + t359 * t377;
t319 = qJD(2) * pkin(2) + t325;
t280 = t359 * t319 - t426;
t357 = t363 ^ 2;
t420 = -t366 ^ 2 + t357;
t417 = qJD(2) * t363;
t415 = -t440 + t414;
t408 = pkin(2) * t404 + qJDD(3);
t407 = pkin(2) * t417;
t271 = qJD(2) * qJ(4) + t280;
t345 = -pkin(2) * t360 - pkin(3);
t259 = t304 * t359 - t360 * t377;
t279 = t319 * t360 + t428;
t286 = t325 * t359 - t426;
t290 = -t330 * t359 - t360 * t405;
t396 = qJD(4) - t279;
t246 = qJD(2) * t444 + t396 - t440;
t252 = t271 + t441;
t394 = t246 * t365 - t252 * t362;
t393 = -t246 * t362 - t252 * t365;
t261 = -pkin(7) * t321 + t290;
t320 = -t425 + t427;
t262 = pkin(7) * t320 + t291;
t392 = t261 * t365 - t262 * t362;
t391 = t261 * t362 + t262 * t365;
t388 = t365 * t320 - t321 * t362;
t284 = t320 * t362 + t321 * t365;
t341 = -pkin(4) + t345;
t343 = pkin(2) * t359 + qJ(4);
t387 = t341 * t365 - t343 * t362;
t386 = t341 * t362 + t343 * t365;
t385 = qJ(4) * t321 + t347;
t238 = -qJDD(2) * pkin(3) - t402;
t384 = -pkin(2) * t418 - qJ(4) * t308;
t383 = -0.2e1 * pkin(1) * t413 - pkin(6) * qJDD(2);
t381 = -qJDD(2) * t290 + t349 * t450;
t378 = qJ(4) * t282 + qJD(4) * t311 - t408;
t313 = qJD(2) * t425 - t359 * t417;
t375 = qJ(4) * t313 + qJD(4) * t321 - t407;
t374 = qJD(2) * t286 + t348 * t398 - t421 - t436;
t368 = qJD(2) ^ 2;
t373 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t368 + t450;
t369 = qJD(1) ^ 2;
t372 = pkin(1) * t369 - pkin(6) * qJDD(1) + t398;
t371 = pkin(3) * t281 - t378;
t370 = t259 * t311 - t260 * t308 - t291 * t281 + t282 * t290 - t398;
t346 = t361 * t367;
t329 = -pkin(3) * t359 + qJ(4) * t360;
t327 = pkin(3) * t360 + qJ(4) * t359 + pkin(2);
t301 = -t395 + t408;
t288 = t327 * t366 + t329 * t363 + pkin(1);
t278 = pkin(3) * t320 - t385;
t263 = -qJD(2) * pkin(3) + t396;
t258 = pkin(3) * t311 - t384;
t255 = t286 + t441;
t254 = t320 * t444 + t385;
t251 = pkin(3) * t310 - t375;
t249 = pkin(7) * t310 + t260;
t248 = -pkin(7) * t313 + t259;
t247 = t311 * t444 + t384;
t243 = t310 * t444 + t375;
t242 = qJD(5) * t284 - t365 * t310 + t313 * t362;
t241 = qJD(5) * t388 + t310 * t362 + t313 * t365;
t236 = -t395 + t371;
t231 = t281 * t444 + t378 + t395;
t1 = [qJDD(1) * MDP(1) + t450 * MDP(2) + t398 * MDP(3) + (qJDD(1) * t357 + 0.2e1 * t363 * t403) * MDP(4) + 0.2e1 * (t363 * t411 - t413 * t420) * MDP(5) + (qJDD(2) * t363 + t366 * t368) * MDP(6) + (qJDD(2) * t366 - t363 * t368) * MDP(7) + (t363 * t383 + t366 * t373) * MDP(9) + (-t363 * t373 + t366 * t383) * MDP(10) + (-t281 * t347 + t301 * t320 + t310 * t328 + (t308 * t443 - t259) * qJD(2) + t381) * MDP(11) + (-t282 * t347 + t301 * t321 + t313 * t328 + (t311 * t443 - t260) * qJD(2) + t447) * MDP(12) + (-t240 * t320 - t279 * t313 - t280 * t310 + t321 * t421 + t370) * MDP(13) + (t240 * t291 + t280 * t260 + t421 * t290 - t279 * t259 - t301 * t347 + t328 * t407 - g(1) * (-t347 * t364 - t346) - g(2) * (t347 * t367 - t424)) * MDP(14) + (-qJD(2) * t259 + t236 * t320 + t251 * t308 + t257 * t310 + t278 * t281 + t381) * MDP(15) + (-t237 * t320 + t238 * t321 + t263 * t313 - t271 * t310 + t370) * MDP(16) + (qJD(2) * t260 - t236 * t321 - t251 * t311 - t257 * t313 - t278 * t282 - t447) * MDP(17) + (t237 * t291 + t271 * t260 + t236 * t278 + t257 * t251 + t238 * t290 + t263 * t259 - g(1) * (-t288 * t364 - t346) - g(2) * (t288 * t367 - t424)) * MDP(18) + (t234 * t284 + t241 * t449) * MDP(19) + (t234 * t388 - t235 * t284 + t241 * t458 - t242 * t449) * MDP(20) + (-t241 * t352 - t284 * t351) * MDP(21) + (t242 * t352 - t351 * t388) * MDP(22) + (-t243 * t458 + t254 * t235 - t231 * t388 + t244 * t242 - (-qJD(5) * t391 + t248 * t365 - t249 * t362) * t352 - t392 * t351 - t450 * t302) * MDP(24) + (t243 * t449 + t254 * t234 + t231 * t284 + t244 * t241 + (qJD(5) * t392 + t248 * t362 + t249 * t365) * t352 + t391 * t351 + t450 * t303) * MDP(25); MDP(6) * t412 + MDP(7) * t411 + qJDD(2) * MDP(8) + (t363 * t372 - t435) * MDP(9) + (g(3) * t363 + t366 * t372) * MDP(10) + (-t311 * t328 + (qJDD(2) * t360 - t308 * t418) * pkin(2) + t374) * MDP(11) + (t308 * t328 + (-qJDD(2) * t359 - t311 * t418) * pkin(2) + t446) * MDP(12) + ((t280 - t286) * t311 + (-t279 + t287) * t308 + (-t281 * t359 - t282 * t360) * pkin(2)) * MDP(13) + (t279 * t286 - t280 * t287 + (-t435 - t421 * t360 + t240 * t359 + (-qJD(1) * t328 + t398) * t363) * pkin(2)) * MDP(14) + (-t433 - t258 * t308 - qJDD(4) + (pkin(3) - t345) * qJDD(2) + t374) * MDP(15) + (-t281 * t343 + t282 * t345 + (t271 - t286) * t311 + (t263 - t414) * t308) * MDP(16) + (qJDD(2) * t343 - t257 * t308 + t258 * t311 + t354 + 0.2e1 * t355 - t446) * MDP(17) + (t237 * t343 + t238 * t345 - t257 * t258 - t263 * t286 - g(3) * (pkin(3) * t349 + qJ(4) * t348 + t442) - t398 * (-t327 * t363 + t329 * t366) + t414 * t271) * MDP(18) + (t235 + t432) * MDP(22) + (-t387 * t351 + t247 * t458 + (t255 * t365 + t362 * t415) * t352 + (t352 * t386 - t393) * qJD(5) + t455) * MDP(24) + (t386 * t351 - t247 * t449 + (-t255 * t362 + t365 * t415) * t352 + (t352 * t387 + t394) * qJD(5) + t457) * MDP(25) + (-MDP(4) * t363 * t366 + MDP(5) * t420) * t369 - t461; (t279 * t311 + t280 * t308 + t408 - t450) * MDP(14) + (-t263 * t311 + t271 * t308 + t371 - t450) * MDP(18) + (-t235 + t432) * MDP(24) + (-t234 + t431) * MDP(25) + (MDP(12) - MDP(17)) * ((-t308 + t406) * qJD(2) + t376) + (-MDP(14) - MDP(18)) * t395 + (MDP(13) + MDP(16)) * (-t308 ^ 2 - t306) + (0.2e1 * qJD(2) * t311 + qJDD(1) * t427 - t338) * (MDP(11) + MDP(15)); (t308 * t311 - qJDD(2)) * MDP(15) + ((t308 + t406) * qJD(2) + t376) * MDP(16) + (-t306 - t368) * MDP(17) + (-qJD(2) * t271 - t321 * t398 + t238 + t433 + t436) * MDP(18) + (t311 * t458 - t365 * t351) * MDP(24) + (-t311 * t449 + t362 * t351) * MDP(25) + (-MDP(24) * t362 - MDP(25) * t365) * t352 ^ 2; (-t399 - t432) * MDP(22) + (t352 * t393 - t455) * MDP(24) + (-t352 * t394 - t457) * MDP(25) + (-MDP(22) * t449 + MDP(24) * t393 - MDP(25) * t394) * qJD(5) + t461;];
tau = t1;
