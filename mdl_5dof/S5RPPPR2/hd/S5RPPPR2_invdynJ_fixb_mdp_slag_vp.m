% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:43
% EndTime: 2019-12-05 17:31:49
% DurationCPUTime: 3.57s
% Computational Cost: add. (1461->341), mult. (3651->500), div. (0->0), fcn. (2888->10), ass. (0->159)
t322 = sin(pkin(9));
t325 = cos(pkin(9));
t327 = cos(pkin(7));
t324 = sin(pkin(7));
t326 = cos(pkin(8));
t402 = t324 * t326;
t278 = t322 * t402 + t325 * t327;
t270 = t278 * qJD(1);
t265 = qJD(5) + t270;
t323 = sin(pkin(8));
t329 = sin(qJ(1));
t331 = cos(qJ(1));
t397 = t327 * t331;
t285 = t323 * t329 + t326 * t397;
t400 = t324 * t331;
t252 = t285 * t325 + t322 * t400;
t394 = t329 * t326;
t284 = t323 * t397 - t394;
t328 = sin(qJ(5));
t330 = cos(qJ(5));
t416 = t252 * t328 - t284 * t330;
t415 = t252 * t330 + t284 * t328;
t380 = qJ(2) * qJDD(1);
t414 = g(2) * t329;
t413 = g(2) * t331;
t412 = g(3) * t331;
t411 = qJDD(1) * pkin(1);
t389 = qJD(1) * t324;
t366 = t325 * t389;
t387 = qJD(1) * t327;
t367 = t322 * t387;
t273 = t326 * t366 - t367;
t410 = t273 * t328;
t406 = t322 * t323;
t405 = t323 * t324;
t404 = t323 * t328;
t403 = t323 * t330;
t401 = t324 * t329;
t399 = t326 * t327;
t398 = t327 * t329;
t332 = qJD(1) ^ 2;
t396 = t327 * t332;
t268 = t278 * qJDD(1);
t264 = qJDD(5) + t268;
t395 = t328 * t264;
t393 = t330 * t264;
t344 = pkin(2) * t327 + qJ(3) * t324 + pkin(1);
t385 = qJD(3) * t324;
t259 = -qJD(1) * t385 - qJDD(1) * t344 + qJDD(2);
t374 = qJDD(1) * t327;
t360 = t326 * t374;
t379 = qJD(1) * qJD(2);
t365 = t327 * t379;
t232 = qJ(2) * t360 + t323 * t259 + t326 * t365;
t378 = qJD(1) * qJD(4);
t225 = (-qJ(4) * qJDD(1) - t378) * t327 + t232;
t376 = qJDD(1) * t324;
t286 = qJ(2) * t376 + t324 * t379 + qJDD(3);
t351 = pkin(3) * t323 - qJ(4) * t326;
t236 = (qJDD(1) * t351 - t326 * t378) * t324 + t286;
t211 = t325 * t225 + t322 * t236;
t277 = -qJD(1) * t344 + qJD(2);
t371 = qJ(2) * t387;
t247 = t323 * t277 + t326 * t371;
t239 = -qJ(4) * t387 + t247;
t300 = qJ(2) * t389 + qJD(3);
t257 = t351 * t389 + t300;
t218 = t325 * t239 + t322 * t257;
t262 = qJ(2) * t399 - t323 * t344;
t254 = -qJ(4) * t327 + t262;
t266 = (qJ(2) + t351) * t324;
t227 = t325 * t254 + t322 * t266;
t392 = g(1) * t327 + g(2) * t401;
t319 = t324 ^ 2;
t391 = t327 ^ 2 + t319;
t390 = qJD(1) * t323;
t388 = qJD(1) * t326;
t386 = qJD(2) * t327;
t384 = qJD(5) * t265;
t383 = qJD(5) * t328;
t382 = t264 * MDP(20);
t377 = qJDD(1) * t323;
t375 = qJDD(1) * t326;
t301 = t322 * t374;
t361 = t325 * t375;
t269 = t324 * t361 - t301;
t368 = t330 * t390;
t354 = t324 * t368;
t363 = t323 * t376;
t372 = qJD(5) * t354 + t330 * t269 + t328 * t363;
t370 = t323 * t389;
t369 = t328 * t390;
t318 = t323 ^ 2;
t364 = t318 * t376;
t362 = t323 * t374;
t359 = t391 * t332;
t209 = pkin(6) * t363 + t211;
t231 = -qJ(2) * t362 + t259 * t326 - t323 * t365;
t228 = pkin(3) * t374 + qJDD(4) - t231;
t213 = pkin(4) * t268 - pkin(6) * t269 + t228;
t358 = -t328 * t209 + t330 * t213;
t357 = t269 * t328 - t330 * t363;
t246 = t277 * t326 - t323 * t371;
t261 = -t323 * t327 * qJ(2) - t326 * t344;
t356 = (-t326 ^ 2 - t318) * MDP(10);
t355 = 0.2e1 * t391;
t282 = t323 * t398 + t326 * t331;
t353 = -g(2) * t284 - g(3) * t282;
t352 = g(3) * t329 + t413;
t256 = t327 * pkin(3) - t261;
t281 = -t323 * t385 + t326 * t386;
t350 = t330 * t209 + t328 * t213;
t238 = pkin(3) * t387 + qJD(4) - t246;
t214 = pkin(4) * t270 - pkin(6) * t273 + t238;
t216 = pkin(6) * t370 + t218;
t349 = t214 * t330 - t216 * t328;
t348 = -t214 * t328 - t216 * t330;
t279 = -t322 * t327 + t325 * t402;
t221 = pkin(4) * t278 - pkin(6) * t279 + t256;
t223 = pkin(6) * t405 + t227;
t347 = t221 * t330 - t223 * t328;
t346 = t221 * t328 + t223 * t330;
t210 = -t225 * t322 + t236 * t325;
t217 = -t239 * t322 + t257 * t325;
t226 = -t254 * t322 + t266 * t325;
t345 = (-qJD(4) * t326 + qJD(2)) * t324;
t343 = -t279 * t328 + t324 * t403;
t249 = t279 * t330 + t324 * t404;
t342 = t325 * t403 - t326 * t328;
t341 = t325 * t404 + t326 * t330;
t280 = t323 * t386 + t326 * t385;
t340 = qJD(1) * t280 - qJDD(1) * t261 - t231;
t339 = qJD(1) * t281 + qJDD(1) * t262 + t232;
t338 = -t352 - t411;
t242 = t273 * t330 + t324 * t369;
t337 = t342 * t265;
t314 = qJDD(2) - t411;
t336 = -t314 - t338;
t335 = t355 * t379 + t414;
t334 = t286 * t324 + (t379 + t380) * t319;
t316 = t331 * qJ(2);
t283 = t323 * t331 - t327 * t394;
t275 = (t322 * t324 + t325 * t399) * qJD(1);
t272 = t326 * t367 - t366;
t263 = -qJD(4) * t327 + t281;
t251 = t283 * t325 - t322 * t401;
t244 = t249 * qJD(5);
t243 = t343 * qJD(5);
t240 = -t354 + t410;
t235 = t325 * t263 + t322 * t345;
t234 = t263 * t322 - t325 * t345;
t230 = t251 * t330 - t282 * t328;
t229 = -t251 * t328 - t282 * t330;
t222 = -pkin(4) * t405 - t226;
t220 = qJD(5) * t242 + t357;
t219 = -t273 * t383 + t372;
t215 = -pkin(4) * t370 - t217;
t208 = -pkin(4) * t363 - t210;
t1 = [qJDD(1) * MDP(1) + t352 * MDP(2) + (t412 - t414) * MDP(3) + t336 * t327 * MDP(4) + (t355 * t380 + t335 - t412) * MDP(6) + (-g(3) * t316 + (-t314 + t352) * pkin(1) + (t380 * t391 + t335) * qJ(2)) * MDP(7) + (g(2) * t285 - g(3) * t283 + t323 * t334 + t327 * t340) * MDP(8) + (t326 * t334 + t327 * t339 + t353) * MDP(9) + (t232 * t262 + t247 * t281 + t231 * t261 - t246 * t280 - g(2) * (-pkin(1) * t331 - pkin(2) * t397 - qJ(2) * t329) - g(3) * (-pkin(1) * t329 - pkin(2) * t398 + t316)) * MDP(11) + (g(2) * t252 - g(3) * t251 + t228 * t278 + t256 * t268 + t270 * t280 + (-qJD(1) * t234 + qJDD(1) * t226 + t210) * t405) * MDP(12) + (t280 * t273 + t256 * t269 + t228 * t279 - g(2) * (t285 * t322 - t325 * t400) - g(3) * (-t283 * t322 - t325 * t401) + (-qJD(1) * t235 - qJDD(1) * t227 - t211) * t405) * MDP(13) + (-t210 * t279 - t211 * t278 - t226 * t269 - t227 * t268 + t234 * t273 - t235 * t270 - t353) * MDP(14) + (t211 * t227 + t218 * t235 + t210 * t226 - t217 * t234 + t228 * t256 + t238 * t280 - g(2) * (-pkin(3) * t285 - qJ(4) * t284) - g(3) * (pkin(3) * t283 - qJ(4) * t282 + t316) + t344 * t413 + (g(2) * qJ(2) + g(3) * t344) * t329) * MDP(15) + (t219 * t249 + t242 * t243) * MDP(16) + (t219 * t343 - t220 * t249 - t240 * t243 - t242 * t244) * MDP(17) + (t219 * t278 + t243 * t265 + t249 * t264) * MDP(18) + (-t220 * t278 - t244 * t265 + t264 * t343) * MDP(19) + t278 * t382 + ((-t235 * t328 + t280 * t330) * t265 + t347 * t264 + t358 * t278 + t234 * t240 + t222 * t220 - t208 * t343 + t215 * t244 + g(2) * t415 - g(3) * t230 + (-t265 * t346 + t278 * t348) * qJD(5)) * MDP(21) + (-(t235 * t330 + t280 * t328) * t265 - t346 * t264 - t350 * t278 + t234 * t242 + t222 * t219 + t208 * t249 + t215 * t243 - g(2) * t416 - g(3) * t229 + (-t265 * t347 - t278 * t349) * qJD(5)) * MDP(22) + (-t336 * MDP(5) + (-t323 * t339 + t326 * t340 + t352) * MDP(10) + (t286 * qJ(2) + qJ(3) * t352 + t300 * qJD(2)) * MDP(11)) * t324; -MDP(4) * t374 - MDP(6) * t359 + (-qJ(2) * t359 + qJDD(2) + t338) * MDP(7) + (-t323 * t359 - t360) * MDP(8) + (-t326 * t359 + t362) * MDP(9) + (t231 * t326 + t232 * t323 + (-t300 * t324 + (t246 * t323 - t247 * t326) * t327) * qJD(1) - t352) * MDP(11) + (-t322 * t364 - t268 * t326 + (-t270 * t327 + t272 * t324) * t390) * MDP(12) + (-t325 * t364 - t269 * t326 + (-t273 * t327 + t275 * t324) * t390) * MDP(13) + (t270 * t275 - t272 * t273 + (-t268 * t325 + t269 * t322) * t323) * MDP(14) + (t217 * t272 - t218 * t275 - t228 * t326 + (-t210 * t322 + t211 * t325 - t238 * t387) * t323 - t352) * MDP(15) + (-t341 * t264 + t220 * t406 - (-t275 * t328 + t327 * t368) * t265 - t272 * t240 - qJD(5) * t337) * MDP(21) + (-t342 * t264 + t219 * t406 + (t275 * t330 + t327 * t369) * t265 - t272 * t242 + t341 * t384) * MDP(22) + (MDP(5) + t356) * t376; (t286 + t392) * MDP(11) + (-t268 * t322 - t269 * t325 + (-t270 * t325 + t273 * t322) * t370) * MDP(14) + (t210 * t325 + t211 * t322 + t392) * MDP(15) + (-t325 * t220 + (-t330 * t384 - t395) * t322 + (t240 * t406 - t265 * t341) * t389) * MDP(21) + (-t325 * t219 + (t265 * t383 - t393) * t322 + (t242 * t406 - t337) * t389) * MDP(22) + ((-t326 * t396 + t377) * MDP(8) + (t323 * t396 + t375) * MDP(9) + (-t412 + (t246 * t326 + t247 * t323) * qJD(1)) * MDP(11) + (-t270 * t388 + t325 * t377) * MDP(12) + (-t273 * t388 - t322 * t377) * MDP(13) + (-t412 + (-t238 * t326 + (-t217 * t322 + t218 * t325) * t323) * qJD(1)) * MDP(15)) * t324 + (t356 + (-MDP(12) * t322 - MDP(13) * t325) * t318) * t319 * t332; t325 * MDP(12) * t374 - t301 * MDP(13) + (-t270 ^ 2 - t273 ^ 2) * MDP(14) + (g(2) * t282 - g(3) * t284 + t217 * t273 + t218 * t270 + t228) * MDP(15) + (-t240 * t273 + t393) * MDP(21) + (-t242 * t273 - t395) * MDP(22) + ((t273 * t390 + t322 * t375) * MDP(12) + (-t270 * t390 + t361) * MDP(13) - g(1) * t323 * MDP(15)) * t324 - (MDP(21) * t328 + MDP(22) * t330) * t265 ^ 2; t242 * t240 * MDP(16) + (-t240 ^ 2 + t242 ^ 2) * MDP(17) + (t240 * t265 + t372) * MDP(18) + (t242 * t265 - t357) * MDP(19) + t382 + (-g(1) * t343 - g(2) * t229 + g(3) * t416 - t215 * t242 - t348 * t265 + t358) * MDP(21) + (g(1) * t249 + g(2) * t230 + g(3) * t415 + t215 * t240 + t349 * t265 - t350) * MDP(22) + (-MDP(18) * t410 - MDP(19) * t242 + MDP(21) * t348 - MDP(22) * t349) * qJD(5);];
tau = t1;
