% Calculate vector of inverse dynamics joint torques for
% S5RRPPR7
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
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:09
% EndTime: 2021-01-15 19:59:20
% DurationCPUTime: 4.47s
% Computational Cost: add. (2050->392), mult. (4713->483), div. (0->0), fcn. (3265->10), ass. (0->171)
t325 = sin(pkin(8));
t329 = sin(qJ(2));
t386 = qJD(1) * t329;
t326 = cos(pkin(8));
t332 = cos(qJ(2));
t394 = t326 * t332;
t283 = -qJD(1) * t394 + t325 * t386;
t328 = sin(qJ(5));
t331 = cos(qJ(5));
t266 = qJD(2) * t328 - t331 * t283;
t395 = t326 * t329;
t397 = t325 * t332;
t294 = t395 + t397;
t286 = t294 * qJD(1);
t424 = qJD(5) + t286;
t425 = t266 * t424;
t268 = qJD(2) * t331 + t283 * t328;
t358 = t424 * t268;
t423 = t328 * t424;
t330 = sin(qJ(1));
t333 = cos(qJ(1));
t421 = g(1) * t330 - g(2) * t333;
t376 = qJD(1) * qJD(2);
t365 = t329 * t376;
t303 = t325 * t365;
t364 = t332 * t376;
t258 = qJDD(1) * t294 + t326 * t364 - t303;
t420 = -qJ(4) * t258 - qJD(4) * t286;
t327 = -qJ(3) - pkin(6);
t302 = t327 * t332;
t297 = qJD(1) * t302;
t289 = t325 * t297;
t301 = t327 * t329;
t296 = qJD(1) * t301;
t261 = t296 * t326 + t289;
t379 = -qJD(4) + t261;
t292 = qJD(2) * pkin(2) + t296;
t396 = t326 * t297;
t256 = t325 * t292 - t396;
t247 = -qJD(2) * qJ(4) - t256;
t410 = pkin(4) * t283;
t233 = -t247 - t410;
t253 = qJDD(5) + t258;
t260 = t296 * t325 - t396;
t312 = -pkin(2) * t326 - pkin(3);
t307 = -pkin(7) + t312;
t419 = t307 * t253 + (t233 - t260 + t410) * t424;
t265 = t301 * t325 - t302 * t326;
t321 = qJ(2) + pkin(8);
t317 = sin(t321);
t418 = -qJDD(2) * t265 - t317 * t421;
t355 = g(1) * t333 + g(2) * t330;
t417 = -qJD(2) * t260 - t355 * t317;
t264 = -t326 * t301 - t302 * t325;
t318 = cos(t321);
t416 = -qJDD(2) * t264 + t318 * t421;
t342 = -g(3) * t317 - t355 * t318;
t282 = t286 ^ 2;
t414 = pkin(3) + pkin(7);
t413 = pkin(2) * t329;
t412 = pkin(2) * t332;
t285 = t294 * qJD(2);
t373 = qJDD(1) * t332;
t304 = t326 * t373;
t374 = qJDD(1) * t329;
t257 = qJD(1) * t285 + t325 * t374 - t304;
t411 = pkin(3) * t257;
t409 = pkin(4) * t286;
t310 = g(3) * t318;
t406 = g(3) * t332;
t404 = qJDD(2) * pkin(3);
t380 = qJD(5) * t331;
t366 = t331 * qJDD(2) + t328 * t257 + t283 * t380;
t375 = qJD(2) * qJD(5);
t225 = -t328 * t375 + t366;
t403 = t225 * t331;
t293 = t325 * t329 - t394;
t314 = pkin(1) + t412;
t350 = -qJ(4) * t294 - t314;
t234 = t414 * t293 + t350;
t402 = t234 * t253;
t401 = t253 * t328;
t400 = t266 * t283;
t399 = t268 * t283;
t398 = t293 * t328;
t393 = t327 * t330;
t392 = t328 * t330;
t391 = t328 * t333;
t390 = t330 * t331;
t246 = t331 * t253;
t389 = t331 * t333;
t361 = qJD(2) * t327;
t281 = -qJD(3) * t329 + t332 * t361;
t252 = qJDD(2) * pkin(2) + t281 * qJD(1) + qJDD(1) * t301;
t280 = qJD(3) * t332 + t329 * t361;
t259 = t280 * qJD(1) - qJDD(1) * t302;
t388 = -t326 * t252 + t325 * t259;
t223 = t325 * t252 + t326 * t259;
t323 = t329 ^ 2;
t387 = -t332 ^ 2 + t323;
t384 = qJD(2) * t329;
t299 = -t314 * qJD(1) + qJD(3);
t341 = -qJ(4) * t286 + t299;
t227 = t414 * t283 + t341;
t382 = qJD(5) * t227;
t381 = qJD(5) * t293;
t378 = t409 - t379;
t370 = MDP(11) - MDP(16);
t369 = MDP(12) - MDP(17);
t368 = pkin(2) * t365 + qJDD(3);
t316 = pkin(2) * t384;
t367 = qJDD(2) * qJ(4) + t223;
t363 = qJDD(4) + t388;
t360 = pkin(2) * t386 + qJ(4) * t283;
t239 = t280 * t325 - t326 * t281;
t255 = t292 * t326 + t289;
t359 = t331 * t424;
t273 = -t314 * qJDD(1) + t368;
t338 = t273 + t420;
t213 = t257 * t414 + t338;
t353 = qJD(4) - t255;
t228 = -t414 * qJD(2) + t353 + t409;
t357 = qJD(5) * t228 + t213;
t356 = MDP(24) * t424;
t352 = -t382 + t310;
t215 = t227 * t331 + t228 * t328;
t240 = t280 * t326 + t281 * t325;
t349 = -t421 + t368;
t220 = -qJD(2) * qJD(4) - t367;
t348 = -0.2e1 * pkin(1) * t376 - pkin(6) * qJDD(2);
t347 = t285 * t328 + t293 * t380;
t288 = qJD(2) * t394 - t325 * t384;
t346 = -qJ(4) * t288 - qJD(4) * t294 + t316;
t237 = pkin(3) * t283 + t341;
t345 = t237 * t286 + t310 + t363;
t218 = -pkin(4) * t257 - t220;
t241 = pkin(4) * t294 + t264;
t344 = t218 * t293 + t233 * t285 - t241 * t253;
t343 = t370 * t329 + t369 * t332;
t334 = qJD(2) ^ 2;
t340 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t334 + t421;
t335 = qJD(1) ^ 2;
t339 = pkin(1) * t335 - pkin(6) * qJDD(1) + t355;
t337 = t239 * t286 - t240 * t283 - t257 * t265 + t258 * t264 - t355;
t336 = t218 + (-qJD(5) * t307 + t414 * t286 + t360) * t424 + t342;
t313 = t327 * t333;
t309 = pkin(2) * t325 + qJ(4);
t300 = -pkin(3) * t325 + qJ(4) * t326;
t298 = pkin(3) * t326 + qJ(4) * t325 + pkin(2);
t279 = -t317 * t392 + t389;
t278 = t317 * t390 + t391;
t277 = t317 * t391 + t390;
t276 = t317 * t389 - t392;
t274 = qJD(2) * t283;
t262 = t298 * t332 + t300 * t329 + pkin(1);
t254 = pkin(3) * t293 + t350;
t249 = t331 * t257;
t243 = -qJD(2) * pkin(3) + t353;
t242 = -pkin(4) * t293 + t265;
t238 = pkin(3) * t286 + t360;
t232 = pkin(3) * t285 + t346;
t231 = -pkin(4) * t285 + t240;
t230 = pkin(4) * t288 + t239;
t226 = qJD(5) * t268 + qJDD(2) * t328 - t249;
t224 = t414 * t285 + t346;
t221 = t363 - t404;
t219 = t338 + t411;
t217 = pkin(4) * t258 - qJDD(2) * t414 + t363;
t216 = t331 * t217;
t214 = -t227 * t328 + t228 * t331;
t1 = [(qJDD(1) * t323 + 0.2e1 * t329 * t364) * MDP(4) + 0.2e1 * (t329 * t373 - t376 * t387) * MDP(5) + (qJDD(2) * t329 + t332 * t334) * MDP(6) + (qJDD(2) * t332 - t329 * t334) * MDP(7) + (t329 * t348 + t332 * t340) * MDP(9) + (-t329 * t340 + t332 * t348) * MDP(10) + (-t257 * t314 + t273 * t293 + t285 * t299 + (t283 * t413 - t239) * qJD(2) + t416) * MDP(11) + (-t258 * t314 + t273 * t294 + t288 * t299 + (t286 * t413 - t240) * qJD(2) + t418) * MDP(12) + (-t223 * t293 - t255 * t288 - t256 * t285 + t294 * t388 + t337) * MDP(13) + (t223 * t265 + t256 * t240 + t388 * t264 - t255 * t239 - t273 * t314 + t299 * t316 - g(1) * (-t314 * t330 - t313) - g(2) * (t314 * t333 - t393)) * MDP(14) + (t220 * t293 + t221 * t294 + t243 * t288 + t247 * t285 + t337) * MDP(15) + (qJD(2) * t239 - t219 * t293 - t232 * t283 - t237 * t285 - t254 * t257 - t416) * MDP(16) + (qJD(2) * t240 - t219 * t294 - t232 * t286 - t237 * t288 - t254 * t258 - t418) * MDP(17) + (t219 * t254 + t237 * t232 - t220 * t265 - t247 * t240 + t221 * t264 + t243 * t239 - g(1) * (-t262 * t330 - t313) - g(2) * (t262 * t333 - t393)) * MDP(18) + (t225 * t398 + t268 * t347) * MDP(19) + ((-t266 * t328 + t268 * t331) * t285 + (t403 - t226 * t328 + (-t266 * t331 - t268 * t328) * qJD(5)) * t293) * MDP(20) + (t225 * t294 + t253 * t398 + t268 * t288 + t347 * t424) * MDP(21) + (t293 * t246 - t226 * t294 - t266 * t288 + (t285 * t331 - t328 * t381) * t424) * MDP(22) + (t253 * t294 + t288 * t424) * MDP(23) + (-g(1) * t279 - g(2) * t277 + t214 * t288 + t216 * t294 + t242 * t226 + t231 * t266 + (-t213 * t294 - t224 * t424 - t402) * t328 + (t230 * t424 - t344) * t331 + ((-t234 * t331 - t241 * t328) * t424 - t215 * t294 + t233 * t398) * qJD(5)) * MDP(24) + (g(1) * t278 - g(2) * t276 - t215 * t288 + t242 * t225 + t231 * t268 + (-(qJD(5) * t241 + t224) * t424 - t402 - t357 * t294 + t233 * t381) * t331 + (-(-qJD(5) * t234 + t230) * t424 - (t217 - t382) * t294 + t344) * t328) * MDP(25) + qJDD(1) * MDP(1) + t421 * MDP(2) + t355 * MDP(3); MDP(6) * t374 + MDP(7) * t373 + qJDD(2) * MDP(8) + (t329 * t339 - t406) * MDP(9) + (g(3) * t329 + t332 * t339) * MDP(10) + (-t286 * t299 - t310 + (qJDD(2) * t326 - t283 * t386) * pkin(2) - t388 - t417) * MDP(11) + (qJD(2) * t261 + t283 * t299 + (-qJDD(2) * t325 - t286 * t386) * pkin(2) - t223 - t342) * MDP(12) + ((t256 - t260) * t286 + (-t255 + t261) * t283 + (-t257 * t325 - t258 * t326) * pkin(2)) * MDP(13) + (t255 * t260 - t256 * t261 + (-t406 - t388 * t326 + t223 * t325 + (-qJD(1) * t299 + t355) * t329) * pkin(2)) * MDP(14) + (-t257 * t309 + t258 * t312 + (-t247 - t260) * t286 + (t243 + t379) * t283) * MDP(15) + (t238 * t283 + (-pkin(3) + t312) * qJDD(2) + t345 + t417) * MDP(16) + (qJDD(2) * t309 - t237 * t283 + t238 * t286 + (0.2e1 * qJD(4) - t261) * qJD(2) + t342 + t367) * MDP(17) + (-t220 * t309 + t221 * t312 - t237 * t238 - t243 * t260 - g(3) * (pkin(3) * t318 + qJ(4) * t317 + t412) - t355 * (-t298 * t329 + t300 * t332) + t379 * t247) * MDP(18) + (-t328 * t358 + t403) * MDP(19) + ((-t226 - t358) * t331 + (-t225 + t425) * t328) * MDP(20) + (-t423 * t424 + t246 + t399) * MDP(21) + (-t359 * t424 - t400 - t401) * MDP(22) + t424 * t283 * MDP(23) + (t214 * t283 + t309 * t226 + t378 * t266 + t336 * t328 + t419 * t331) * MDP(24) + (-t215 * t283 + t309 * t225 + t378 * t268 - t419 * t328 + t336 * t331) * MDP(25) + (-MDP(4) * t329 * t332 + MDP(5) * t387) * t335; -t303 * MDP(12) + (t255 * t286 + t256 * t283 + t349) * MDP(14) + (t274 + t303) * MDP(17) + (-t243 * t286 - t247 * t283 + t349 + t411 + t420) * MDP(18) + (t400 - t401) * MDP(24) + (-t246 + t399) * MDP(25) - t370 * t304 + (MDP(25) * t423 - t331 * t356) * t424 + (t343 * t325 + t369 * t395 + (-MDP(14) - MDP(18)) * t314) * qJDD(1) + (-t283 * MDP(12) + t370 * t286 + (t326 * t343 + t370 * t397) * qJD(1)) * qJD(2) + (MDP(13) + MDP(15)) * (-t283 ^ 2 - t282); (t274 + t258) * MDP(15) + (-t283 * t286 + qJDD(2)) * MDP(16) + (-t282 - t334) * MDP(17) + (qJD(2) * t247 - t294 * t355 + t345 - t404) * MDP(18) + (-qJD(2) * t266 + t246) * MDP(24) + (-qJD(2) * t268 - t401) * MDP(25) + (-MDP(25) * t359 - t328 * t356) * t424; t268 * t266 * MDP(19) + (-t266 ^ 2 + t268 ^ 2) * MDP(20) + (t366 + t425) * MDP(21) + (t249 + t358) * MDP(22) + t253 * MDP(23) + (-g(1) * t276 - g(2) * t278 + t215 * t424 - t233 * t268 + t216) * MDP(24) + (g(1) * t277 - g(2) * t279 + t214 * t424 + t233 * t266) * MDP(25) + (-MDP(22) * t375 + MDP(24) * t352 - MDP(25) * t357) * t331 + (-MDP(21) * t375 + (-qJD(5) * t283 - qJDD(2)) * MDP(22) - t357 * MDP(24) + (-t217 - t352) * MDP(25)) * t328;];
tau = t1;
