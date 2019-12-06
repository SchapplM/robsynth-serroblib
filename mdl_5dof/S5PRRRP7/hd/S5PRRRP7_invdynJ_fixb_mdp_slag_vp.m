% Calculate vector of inverse dynamics joint torques for
% S5PRRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:36
% EndTime: 2019-12-05 16:56:42
% DurationCPUTime: 3.96s
% Computational Cost: add. (1705->343), mult. (3999->473), div. (0->0), fcn. (3004->10), ass. (0->159)
t319 = cos(pkin(5));
t326 = cos(qJ(2));
t420 = cos(pkin(9));
t356 = t420 * t326;
t317 = sin(pkin(9));
t323 = sin(qJ(2));
t411 = t317 * t323;
t267 = -t319 * t356 + t411;
t357 = t420 * t323;
t410 = t317 * t326;
t269 = t319 * t410 + t357;
t353 = g(1) * t269 + g(2) * t267;
t318 = sin(pkin(5));
t407 = t318 * t326;
t444 = -g(3) * t407 + t353;
t268 = t319 * t357 + t410;
t270 = -t319 * t411 + t356;
t352 = g(1) * t270 + g(2) * t268;
t409 = t318 * t323;
t336 = -g(3) * t409 - t352;
t325 = cos(qJ(3));
t375 = qJD(2) * qJD(3);
t363 = t325 * t375;
t322 = sin(qJ(3));
t373 = qJDD(2) * t322;
t442 = -t363 - t373;
t380 = qJD(4) * t322;
t441 = qJD(2) * t380 - qJDD(3);
t321 = sin(qJ(4));
t324 = cos(qJ(4));
t387 = qJD(2) * t325;
t241 = t321 * (qJD(3) * (qJD(4) + t387) + t373) + t441 * t324;
t440 = pkin(4) * t322;
t354 = pkin(3) * t322 - pkin(8) * t325;
t290 = t354 * qJD(3);
t384 = qJD(3) * t322;
t391 = qJD(1) * t318;
t403 = t325 * t326;
t428 = pkin(7) * t321;
t438 = (-t321 * t403 + t323 * t324) * t391 - t324 * t290 - t384 * t428;
t293 = -pkin(3) * t325 - pkin(8) * t322 - pkin(2);
t379 = qJD(4) * t324;
t437 = -(t321 * t323 + t324 * t403) * t391 + t321 * t290 + t293 * t379;
t291 = qJD(2) * pkin(7) + t323 * t391;
t406 = t319 * t325;
t436 = qJD(1) * t406 - t322 * t291;
t358 = t318 * t420;
t244 = t268 * t325 - t322 * t358;
t246 = t317 * t318 * t322 + t270 * t325;
t408 = t318 * t325;
t273 = t319 * t322 + t323 * t408;
t249 = -t273 * t321 - t324 * t407;
t435 = -g(1) * (-t246 * t321 + t269 * t324) - g(2) * (-t244 * t321 + t267 * t324) - g(3) * t249;
t434 = pkin(4) * t241 + qJDD(5);
t327 = qJD(3) ^ 2;
t376 = qJD(1) * qJD(2);
t364 = t323 * t376;
t351 = -qJDD(1) * t407 + t318 * t364;
t433 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t327 + t318 * (-g(3) * t326 + t364) - t351 + t353;
t390 = qJD(1) * t322;
t304 = t319 * t390;
t261 = t325 * t291 + t304;
t255 = qJD(3) * pkin(8) + t261;
t306 = -qJD(4) + t387;
t431 = qJD(4) * (pkin(7) * t306 + t255) + t444;
t385 = qJD(3) * t321;
t388 = qJD(2) * t322;
t287 = t324 * t388 + t385;
t430 = t287 ^ 2;
t422 = qJ(5) + pkin(8);
t421 = qJD(2) * pkin(2);
t419 = qJ(5) * t322;
t417 = qJDD(3) * pkin(3);
t377 = t324 * qJD(3);
t240 = -qJD(4) * t377 + t441 * t321 + t442 * t324;
t416 = t240 * t321;
t285 = t321 * t388 - t377;
t415 = t285 * t306;
t414 = t287 * t306;
t413 = t287 * t324;
t412 = t306 * t321;
t405 = t322 * t324;
t404 = t324 * t325;
t402 = qJDD(1) - g(3);
t307 = pkin(7) * t404;
t350 = -qJ(5) * t404 + t440;
t378 = qJD(5) * t324;
t401 = -t322 * t378 + t350 * qJD(3) + (-t307 + (-t293 + t419) * t321) * qJD(4) - t438;
t400 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t405 + (-qJD(5) * t322 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t325) * t321 + t437;
t369 = t326 * t391;
t262 = qJD(2) * t293 - t369;
t233 = -t255 * t321 + t324 * t262;
t227 = -qJ(5) * t287 + t233;
t222 = -pkin(4) * t306 + t227;
t399 = -t227 + t222;
t289 = t354 * qJD(2);
t398 = t321 * t289 + t324 * t436;
t359 = qJD(4) * t422;
t397 = t378 - t398 + (qJ(5) * t387 - t359) * t321;
t276 = t324 * t289;
t396 = -qJD(2) * t350 - t324 * t359 - t276 + (-qJD(5) + t436) * t321;
t393 = t321 * t293 + t307;
t315 = t322 ^ 2;
t392 = -t325 ^ 2 + t315;
t389 = qJD(2) * t318;
t386 = qJD(3) * t285;
t383 = qJD(3) * t325;
t382 = qJD(4) * t306;
t381 = qJD(4) * t321;
t374 = qJDD(1) * t319;
t372 = t325 * qJDD(2);
t264 = qJDD(2) * pkin(7) + (qJDD(1) * t323 + t326 * t376) * t318;
t361 = t322 * t374;
t230 = qJDD(3) * pkin(8) + qJD(3) * t436 + t264 * t325 + t361;
t238 = qJD(2) * t290 + qJDD(2) * t293 + t351;
t371 = t324 * t230 + t321 * t238 + t262 * t379;
t367 = t323 * t389;
t366 = t326 * t389;
t365 = t306 * t377;
t234 = t255 * t324 + t262 * t321;
t346 = -t273 * t324 + t321 * t407;
t272 = t322 * t409 - t406;
t342 = -t322 * t375 + t372;
t282 = qJDD(4) - t342;
t344 = t282 * t321 - t306 * t379;
t343 = t282 * t324 + t306 * t381;
t341 = t255 * t381 - t371;
t243 = t268 * t322 + t325 * t358;
t245 = t270 * t322 - t317 * t408;
t340 = g(1) * t245 + g(2) * t243 + g(3) * t272;
t339 = g(1) * t246 + g(2) * t244 + g(3) * t273;
t338 = qJD(3) * t304 + t322 * t264 + t291 * t383 - t325 * t374;
t254 = -qJD(3) * pkin(3) - t436;
t335 = -pkin(8) * t282 - t254 * t306;
t231 = t338 - t417;
t237 = t324 * t238;
t333 = -qJD(4) * t234 - t321 * t230 + t237;
t332 = pkin(8) * t382 - t231 + t340;
t292 = -t369 - t421;
t331 = -pkin(7) * qJDD(3) + (t292 + t369 - t421) * qJD(3);
t330 = -t338 + t340;
t328 = qJD(2) ^ 2;
t310 = pkin(4) * t324 + pkin(3);
t295 = t422 * t324;
t294 = t422 * t321;
t284 = t324 * t293;
t281 = t285 ^ 2;
t251 = -t321 * t419 + t393;
t248 = qJD(3) * t273 + t322 * t366;
t247 = -qJD(3) * t272 + t325 * t366;
t242 = -qJ(5) * t405 + t284 + (-pkin(4) - t428) * t325;
t239 = pkin(4) * t285 + qJD(5) + t254;
t228 = -qJ(5) * t285 + t234;
t225 = qJD(4) * t249 + t247 * t324 + t321 * t367;
t224 = qJD(4) * t346 - t247 * t321 + t324 * t367;
t221 = t231 + t434;
t220 = -qJ(5) * t241 - qJD(5) * t285 - t341;
t219 = pkin(4) * t282 + qJ(5) * t240 - qJD(5) * t287 + t333;
t1 = [t402 * MDP(1) + (-qJD(3) * t248 - qJDD(3) * t272) * MDP(10) + (-qJD(3) * t247 - qJDD(3) * t273) * MDP(11) + (-t224 * t306 + t241 * t272 + t248 * t285 + t249 * t282) * MDP(17) + (t225 * t306 - t240 * t272 + t248 * t287 + t282 * t346) * MDP(18) + (-t224 * t287 - t225 * t285 + t240 * t249 + t241 * t346) * MDP(19) + (t219 * t249 - t220 * t346 + t221 * t272 + t222 * t224 + t225 * t228 + t239 * t248 - g(3)) * MDP(20) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t325 + MDP(11) * t322 - MDP(3)) * t328) * t323 + (t342 * MDP(10) + t442 * MDP(11) + qJDD(2) * MDP(3) - t328 * MDP(4)) * t326) * t318; qJDD(2) * MDP(2) + (t402 * t407 + t353) * MDP(3) + (-t402 * t409 + t352) * MDP(4) + (qJDD(2) * t315 + 0.2e1 * t322 * t363) * MDP(5) + 0.2e1 * (t322 * t372 - t375 * t392) * MDP(6) + (qJDD(3) * t322 + t325 * t327) * MDP(7) + (qJDD(3) * t325 - t322 * t327) * MDP(8) + (t331 * t322 + t325 * t433) * MDP(10) + (-t322 * t433 + t331 * t325) * MDP(11) + (-t240 * t405 + (-t321 * t380 + t325 * t377) * t287) * MDP(12) + ((-t285 * t324 - t287 * t321) * t383 + (t416 - t241 * t324 + (t285 * t321 - t413) * qJD(4)) * t322) * MDP(13) + ((t240 - t365) * t325 + (qJD(3) * t287 + t343) * t322) * MDP(14) + ((t306 * t385 + t241) * t325 + (-t344 - t386) * t322) * MDP(15) + (-t282 * t325 - t306 * t384) * MDP(16) + (t284 * t282 + t438 * t306 + (t293 * t382 + t336) * t321 + (pkin(7) * t386 - t237 + (-pkin(7) * t282 + qJD(3) * t254 + qJD(4) * t262 + t230) * t321 + t431 * t324) * t325 + (pkin(7) * t241 + qJD(3) * t233 + t231 * t321 + t254 * t379 - t285 * t369) * t322) * MDP(17) + (-t393 * t282 + t437 * t306 + t336 * t324 + ((pkin(7) * t287 + t254 * t324) * qJD(3) - t431 * t321 + t371) * t325 + (-t287 * t369 - t254 * t381 - t234 * qJD(3) + t231 * t324 + (-t240 - t365) * pkin(7)) * t322) * MDP(18) + (t240 * t242 - t241 * t251 - t401 * t287 - t400 * t285 + (-t222 * t324 - t228 * t321) * t383 + (-t219 * t324 - t220 * t321 + (t222 * t321 - t228 * t324) * qJD(4) + t444) * t322) * MDP(19) + (t219 * t242 + t220 * t251 + t239 * t379 * t440 + t400 * t228 + t401 * t222 - t239 * t390 * t407 + (t221 * t322 + t239 * t383 + t336) * (pkin(4) * t321 + pkin(7)) + t444 * (t310 * t325 + t322 * t422 + pkin(2))) * MDP(20); MDP(7) * t373 + MDP(8) * t372 + qJDD(3) * MDP(9) + (qJD(3) * t261 - t292 * t388 + t330) * MDP(10) + (-t361 + (-qJD(2) * t292 - t264) * t325 + t339) * MDP(11) + (-t306 * t413 - t416) * MDP(12) + ((-t240 + t415) * t324 + (-t241 + t414) * t321) * MDP(13) + ((-t287 * t322 + t306 * t404) * qJD(2) + t344) * MDP(14) + ((t285 * t322 - t325 * t412) * qJD(2) + t343) * MDP(15) + t306 * MDP(16) * t388 + (-t233 * t388 - pkin(3) * t241 - t261 * t285 + t276 * t306 + (-t306 * t436 + t335) * t321 + t332 * t324) * MDP(17) + (pkin(3) * t240 + t234 * t388 - t261 * t287 - t306 * t398 - t321 * t332 + t324 * t335) * MDP(18) + (-t240 * t294 - t241 * t295 - t396 * t287 - t397 * t285 + (t222 * t306 + t220) * t324 + (t228 * t306 - t219) * t321 - t339) * MDP(19) + (t220 * t295 - t219 * t294 - t221 * t310 - g(1) * (-t245 * t310 + t246 * t422) - g(2) * (-t243 * t310 + t244 * t422) - g(3) * (-t272 * t310 + t273 * t422) + (-pkin(4) * t412 - t261) * t239 + t397 * t228 + t396 * t222) * MDP(20) + (-MDP(5) * t322 * t325 + MDP(6) * t392) * t328; t287 * t285 * MDP(12) + (-t281 + t430) * MDP(13) + (-t240 - t415) * MDP(14) + (-t241 - t414) * MDP(15) + t282 * MDP(16) + (-t234 * t306 - t254 * t287 + t333 + t435) * MDP(17) + (-t233 * t306 + t254 * t285 - g(1) * (-t246 * t324 - t269 * t321) - g(2) * (-t244 * t324 - t267 * t321) - g(3) * t346 + t341) * MDP(18) + (pkin(4) * t240 - t285 * t399) * MDP(19) + (t399 * t228 + (-t239 * t287 + t219 + t435) * pkin(4)) * MDP(20); (-t281 - t430) * MDP(19) + (t222 * t287 + t228 * t285 - t330 - t417 + t434) * MDP(20);];
tau = t1;
