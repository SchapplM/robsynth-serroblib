% Calculate vector of inverse dynamics joint torques for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:04
% EndTime: 2019-12-31 20:50:09
% DurationCPUTime: 2.63s
% Computational Cost: add. (2616->320), mult. (3909->370), div. (0->0), fcn. (2359->12), ass. (0->158)
t448 = qJ(4) + pkin(7);
t359 = cos(qJ(3));
t340 = t359 * qJD(4);
t356 = sin(qJ(3));
t394 = qJD(3) * t448;
t297 = -t356 * t394 + t340;
t354 = cos(pkin(8));
t424 = t354 * t359;
t353 = sin(pkin(8));
t427 = t353 * t356;
t303 = -t424 + t427;
t375 = -qJD(4) * t356 - t359 * t394;
t360 = cos(qJ(2));
t415 = qJD(1) * t360;
t407 = pkin(1) * t415;
t420 = t354 * t297 + t303 * t407 + t353 * t375;
t352 = qJ(1) + qJ(2);
t342 = cos(t352);
t328 = g(2) * t342;
t357 = sin(qJ(2));
t414 = qJD(2) * t357;
t337 = pkin(1) * t414;
t438 = pkin(1) * t360;
t417 = -qJD(1) * t337 + qJDD(1) * t438;
t346 = qJDD(1) + qJDD(2);
t437 = pkin(2) * t346;
t447 = -t417 - t437 + t328;
t304 = t353 * t359 + t354 * t356;
t347 = qJD(1) + qJD(2);
t292 = t304 * t347;
t341 = sin(t352);
t445 = g(1) * t342 + g(2) * t341;
t329 = g(1) * t341;
t444 = t329 - t328;
t434 = pkin(3) * t359;
t332 = pkin(2) + t434;
t289 = -t332 * t347 + qJD(4) - t407;
t405 = t347 * t424;
t290 = t347 * t427 - t405;
t251 = pkin(4) * t290 - qJ(5) * t292 + t289;
t348 = qJ(3) + pkin(8);
t338 = sin(t348);
t339 = cos(t348);
t409 = qJDD(1) * t357;
t413 = qJD(2) * t360;
t295 = pkin(7) * t346 + (qJD(1) * t413 + t409) * pkin(1);
t379 = qJ(4) * t346 + qJD(4) * t347 + t295;
t440 = pkin(1) * t357;
t408 = qJD(1) * t440;
t392 = t448 * t347 + t408;
t380 = qJD(3) * t392;
t247 = qJDD(3) * pkin(3) - t356 * t379 - t359 * t380;
t250 = -t356 * t380 + t359 * t379;
t232 = t354 * t247 - t353 * t250;
t396 = -qJDD(5) + t232;
t443 = -g(3) * t339 - t251 * t292 + t445 * t338 + t396;
t233 = t353 * t247 + t354 * t250;
t404 = qJDD(3) * qJ(5) + t233;
t230 = qJD(3) * qJD(5) + t404;
t432 = qJDD(3) * pkin(4);
t231 = -t396 - t432;
t280 = t392 * t356;
t279 = qJD(3) * pkin(3) - t280;
t281 = t392 * t359;
t431 = t281 * t353;
t256 = t279 * t354 - t431;
t254 = -qJD(3) * pkin(4) + qJD(5) - t256;
t275 = t354 * t281;
t257 = t353 * t279 + t275;
t255 = qJD(3) * qJ(5) + t257;
t298 = t304 * qJD(3);
t412 = qJD(3) * t356;
t398 = t353 * t412;
t411 = qJD(3) * t359;
t399 = t354 * t411;
t299 = -t398 + t399;
t442 = -t230 * t303 + t231 * t304 + t254 * t299 - t255 * t298;
t441 = -t232 * t304 - t233 * t303 - t256 * t299 - t257 * t298;
t286 = t292 ^ 2;
t358 = sin(qJ(1));
t439 = pkin(1) * t358;
t436 = pkin(2) * t347;
t435 = pkin(3) * t356;
t433 = g(3) * t359;
t430 = t338 * t342;
t429 = t339 * t342;
t428 = t342 * t448;
t423 = t356 * t359;
t331 = pkin(7) + t440;
t422 = -qJ(4) - t331;
t421 = t297 * t353 - t304 * t407 - t354 * t375;
t309 = -t407 - t436;
t419 = t309 * t412 + t359 * t329;
t350 = t356 ^ 2;
t416 = -t359 ^ 2 + t350;
t260 = -t280 * t354 - t431;
t410 = qJD(5) - t260;
t406 = pkin(1) * t413;
t336 = pkin(3) * t412;
t403 = t309 * t411 + t447 * t356;
t315 = t342 * t332;
t402 = pkin(4) * t429 + qJ(5) * t430 + t315;
t401 = t347 * t414;
t400 = t347 * t412;
t397 = t448 * t356;
t393 = t422 * t356;
t391 = t341 * t448 + t315;
t390 = qJD(3) * t422;
t389 = t347 * t408;
t388 = t417 + t444;
t387 = pkin(3) * t400 + qJDD(4) - t417;
t258 = pkin(4) * t298 - qJ(5) * t299 - qJD(5) * t304 + t336;
t386 = -t258 + t408;
t385 = t332 * t346;
t383 = -pkin(4) * t339 - qJ(5) * t338;
t362 = qJD(3) ^ 2;
t382 = 0.2e1 * (-qJD(3) * t347 * t416 + t346 * t423) * MDP(8) + (t346 * t350 + 0.2e1 * t359 * t400) * MDP(7) + (qJDD(3) * t359 - t356 * t362) * MDP(10) + (qJDD(3) * t356 + t359 * t362) * MDP(9) + t346 * MDP(4);
t381 = -t332 * t341 + t428;
t316 = t346 * t424;
t262 = t298 * t347 + t346 * t427 - t316;
t311 = t347 * t398;
t372 = t304 * t346 - t311;
t263 = t347 * t399 + t372;
t365 = pkin(4) * t262 - qJ(5) * t263 - qJD(5) * t292 + t387;
t229 = -t385 + t365;
t378 = -g(2) * t429 + t229 * t303 + t251 * t298 + t339 * t329;
t377 = -g(2) * t430 - t229 * t304 - t251 * t299 + t338 * t329;
t270 = pkin(4) * t303 - qJ(5) * t304 - t332;
t374 = -t309 * t347 - t295 + t445;
t371 = pkin(7) * t362 - t389 - t437;
t273 = t356 * t390 + t359 * t406 + t340;
t364 = (-qJD(4) - t406) * t356 + t359 * t390;
t248 = t273 * t353 - t354 * t364;
t249 = t354 * t273 + t353 * t364;
t343 = t359 * qJ(4);
t302 = t331 * t359 + t343;
t268 = t302 * t353 - t354 * t393;
t269 = t354 * t302 + t353 * t393;
t370 = t248 * t292 - t249 * t290 - t269 * t262 + t263 * t268 - t445;
t333 = -pkin(2) - t438;
t369 = t401 * pkin(1) + t331 * t362 + t333 * t346;
t368 = -pkin(7) * qJDD(3) + (t407 - t436) * qJD(3);
t367 = -qJDD(3) * t331 + (t333 * t347 - t406) * qJD(3);
t366 = (-g(1) * (-t332 + t383) - g(2) * t448) * t341;
t318 = pkin(7) * t359 + t343;
t277 = t318 * t353 + t354 * t397;
t278 = t354 * t318 - t353 * t397;
t363 = -t278 * t262 + t263 * t277 - t420 * t290 + t292 * t421 - t445;
t361 = cos(qJ(1));
t344 = t361 * pkin(1);
t326 = -pkin(3) * t354 - pkin(4);
t324 = pkin(3) * t353 + qJ(5);
t267 = -t385 + t387;
t266 = t270 - t438;
t261 = pkin(4) * t292 + qJ(5) * t290 + t347 * t435;
t259 = -t280 * t353 + t275;
t252 = t258 + t337;
t1 = [qJDD(1) * MDP(1) + (g(1) * t358 - g(2) * t361) * MDP(2) + (g(1) * t361 + g(2) * t358) * MDP(3) + ((t346 * t360 - t401) * pkin(1) + t388) * MDP(5) + (((-qJDD(1) - t346) * t357 + (-qJD(1) - t347) * t413) * pkin(1) + t445) * MDP(6) + (t367 * t356 + (-t369 - t447) * t359 + t419) * MDP(12) + (t367 * t359 + (t369 - t329) * t356 + t403) * MDP(13) + (t370 + t441) * MDP(14) + (t233 * t269 + t257 * t249 - t232 * t268 - t256 * t248 + t267 * (-t332 - t438) + t289 * (t337 + t336) - g(1) * (t381 - t439) - g(2) * (t344 + t391)) * MDP(15) + (-qJD(3) * t248 - qJDD(3) * t268 + t252 * t290 + t262 * t266 + t378) * MDP(16) + (t370 + t442) * MDP(17) + (qJD(3) * t249 + qJDD(3) * t269 - t252 * t292 - t263 * t266 + t377) * MDP(18) + (t230 * t269 + t255 * t249 + t229 * t266 + t251 * t252 + t231 * t268 + t254 * t248 - g(1) * (t428 - t439) - g(2) * (t344 + t402) + t366) * MDP(19) + t382; (t388 + t389) * MDP(5) + ((-t409 + (-qJD(2) + t347) * t415) * pkin(1) + t445) * MDP(6) + (t368 * t356 + (-t371 - t447) * t359 + t419) * MDP(12) + (t368 * t359 + (t371 - t329) * t356 + t403) * MDP(13) + (t363 + t441) * MDP(14) + (t233 * t278 - t232 * t277 - t267 * t332 - g(1) * t381 - g(2) * t391 + (t336 - t408) * t289 + t420 * t257 - t421 * t256) * MDP(15) + (-qJD(3) * t421 - qJDD(3) * t277 + t262 * t270 - t290 * t386 + t378) * MDP(16) + (t363 + t442) * MDP(17) + (qJD(3) * t420 + qJDD(3) * t278 - t263 * t270 + t292 * t386 + t377) * MDP(18) + (-g(1) * t428 - g(2) * t402 + t229 * t270 + t230 * t278 + t231 * t277 - t251 * t386 + t254 * t421 + t255 * t420 + t366) * MDP(19) + t382; qJDD(3) * MDP(11) + (t356 * t374 - t433) * MDP(12) + (g(3) * t356 + t359 * t374) * MDP(13) + ((t257 - t259) * t292 + (-t256 + t260) * t290 + (-t262 * t353 - t263 * t354) * pkin(3)) * MDP(14) + (t256 * t259 - t257 * t260 + (-t433 + t232 * t354 + t233 * t353 + (-t289 * t347 + t445) * t356) * pkin(3)) * MDP(15) + (qJD(3) * t259 - t261 * t290 + (pkin(4) - t326) * qJDD(3) + t443) * MDP(16) + (-t262 * t324 + t263 * t326 + (t255 - t259) * t292 + (t254 - t410) * t290) * MDP(17) + (-g(3) * t338 + qJDD(3) * t324 - t251 * t290 + t261 * t292 - t445 * t339 + (0.2e1 * qJD(5) - t260) * qJD(3) + t404) * MDP(18) + (t230 * t324 + t231 * t326 - t251 * t261 - t254 * t259 - g(3) * (-t383 + t434) + t410 * t255 + t445 * (pkin(4) * t338 - qJ(5) * t339 + t435)) * MDP(19) + (t359 * MDP(10) + t356 * MDP(9)) * t346 + (-MDP(7) * t423 + t416 * MDP(8)) * t347 ^ 2; (t256 * t292 + t257 * t290 + t387 - t444) * MDP(15) - t316 * MDP(16) + t311 * MDP(18) + (-t254 * t292 + t255 * t290 + t365 - t444) * MDP(19) + (MDP(16) * t427 - MDP(18) * t304 + (-MDP(15) - MDP(19)) * t332) * t346 + (0.2e1 * t292 * MDP(16) + (t290 - t405) * MDP(18)) * qJD(3) + (MDP(14) + MDP(17)) * (-t290 ^ 2 - t286); (t290 * t292 - qJDD(3)) * MDP(16) + ((t290 + t405) * qJD(3) + t372) * MDP(17) + (-t286 - t362) * MDP(18) + (-qJD(3) * t255 - t432 - t443) * MDP(19);];
tau = t1;
