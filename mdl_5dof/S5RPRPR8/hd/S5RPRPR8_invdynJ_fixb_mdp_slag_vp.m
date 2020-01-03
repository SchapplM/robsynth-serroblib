% Calculate vector of inverse dynamics joint torques for
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:20
% EndTime: 2019-12-31 18:22:25
% DurationCPUTime: 4.04s
% Computational Cost: add. (1876->354), mult. (3944->491), div. (0->0), fcn. (2662->14), ass. (0->159)
t339 = cos(qJ(3));
t332 = sin(pkin(8));
t315 = pkin(1) * t332 + pkin(6);
t311 = t315 * qJD(1);
t336 = sin(qJ(3));
t393 = qJD(3) * t339;
t309 = t315 * qJDD(1);
t417 = -qJD(2) * qJD(3) - t309;
t383 = -t311 * t393 + t336 * t417;
t352 = -qJDD(3) * pkin(3) + qJDD(4) - t383;
t234 = -qJDD(2) * t339 + t352;
t328 = qJ(1) + pkin(8);
t322 = sin(t328);
t324 = cos(t328);
t373 = g(1) * t324 + g(2) * t322;
t359 = t373 * t336;
t348 = -g(3) * t339 + t359;
t419 = t234 - t348;
t397 = qJD(1) * t339;
t313 = -qJD(5) + t397;
t331 = sin(pkin(9));
t333 = cos(pkin(9));
t390 = t333 * qJD(3);
t398 = qJD(1) * t336;
t292 = t331 * t398 - t390;
t395 = qJD(3) * t331;
t294 = t333 * t398 + t395;
t335 = sin(qJ(5));
t338 = cos(qJ(5));
t363 = t292 * t335 - t294 * t338;
t418 = t313 * t363;
t403 = qJDD(2) - g(3);
t416 = t403 * t339;
t279 = qJD(2) * t339 - t336 * t311;
t334 = cos(pkin(8));
t415 = pkin(1) * t334;
t414 = g(3) * t336;
t412 = pkin(7) + qJ(4);
t410 = t294 * t335;
t247 = t338 * t292 + t410;
t411 = t247 * t313;
t409 = t322 * t339;
t408 = t324 * t339;
t407 = t331 * t335;
t406 = t331 * t339;
t405 = t333 * t336;
t404 = t333 * t339;
t231 = qJDD(3) * qJ(4) + qJDD(2) * t336 + t309 * t339 + (qJD(4) + t279) * qJD(3);
t370 = pkin(3) * t336 - qJ(4) * t339;
t283 = qJD(3) * t370 - qJD(4) * t336;
t361 = pkin(3) * t339 + qJ(4) * t336 + pkin(2);
t289 = -t361 - t415;
t237 = qJD(1) * t283 + qJDD(1) * t289;
t218 = t333 * t231 + t331 * t237;
t298 = t331 * t338 + t333 * t335;
t356 = t298 * t339;
t391 = qJD(5) * t338;
t392 = qJD(5) * t336;
t239 = qJD(3) * t356 + t391 * t405 - t392 * t407;
t275 = t298 * t336;
t385 = t339 * qJDD(1);
t389 = qJD(1) * qJD(3);
t353 = t336 * t389 - t385;
t296 = qJDD(5) + t353;
t402 = t239 * t313 - t275 * t296;
t280 = t336 * qJD(2) + t339 * t311;
t268 = qJD(3) * qJ(4) + t280;
t271 = t289 * qJD(1);
t226 = t333 * t268 + t331 * t271;
t303 = t370 * qJD(1);
t243 = t333 * t279 + t331 * t303;
t297 = -t338 * t333 + t407;
t355 = t297 * t339;
t401 = qJD(1) * t355 - t297 * qJD(5);
t400 = -qJD(1) * t356 + t298 * qJD(5);
t394 = qJD(3) * t336;
t381 = t315 * t394;
t245 = t333 * t283 + t331 * t381;
t252 = t331 * t289 + t315 * t404;
t329 = t336 ^ 2;
t399 = -t339 ^ 2 + t329;
t316 = -pkin(2) - t415;
t312 = qJD(1) * t316;
t387 = qJDD(1) * t336;
t386 = qJDD(3) * t331;
t319 = t333 * qJDD(3);
t379 = t339 * t389;
t354 = t379 + t387;
t264 = t331 * t354 - t319;
t265 = t333 * t354 + t386;
t384 = -t335 * t264 + t338 * t265 - t292 * t391;
t382 = t331 * t397;
t380 = qJ(4) * t385;
t378 = t331 * t387;
t377 = t333 * t387;
t376 = pkin(4) * t331 + t315;
t217 = -t231 * t331 + t333 * t237;
t215 = pkin(4) * t353 - pkin(7) * t265 + t217;
t216 = -pkin(7) * t264 + t218;
t375 = t338 * t215 - t216 * t335;
t374 = t338 * t264 + t335 * t265;
t225 = -t268 * t331 + t333 * t271;
t242 = -t279 * t331 + t333 * t303;
t372 = g(1) * t322 - g(2) * t324;
t337 = sin(qJ(1));
t340 = cos(qJ(1));
t371 = g(1) * t337 - g(2) * t340;
t369 = t215 * t335 + t216 * t338;
t368 = -t217 * t331 + t218 * t333;
t221 = -pkin(4) * t397 - pkin(7) * t294 + t225;
t223 = -pkin(7) * t292 + t226;
t212 = t221 * t338 - t223 * t335;
t213 = t221 * t335 + t223 * t338;
t367 = -t225 * t331 + t226 * t333;
t278 = t333 * t289;
t236 = -pkin(7) * t405 + t278 + (-t315 * t331 - pkin(4)) * t339;
t241 = -pkin(7) * t331 * t336 + t252;
t366 = t236 * t338 - t241 * t335;
t365 = t236 * t335 + t241 * t338;
t238 = -qJD(3) * t355 - t298 * t392;
t276 = t297 * t336;
t364 = t238 * t313 + t276 * t296;
t362 = pkin(4) * t336 - pkin(7) * t404;
t360 = t371 * pkin(1);
t308 = t412 * t333;
t358 = qJD(1) * t362 + qJD(4) * t331 + qJD(5) * t308 + t242;
t307 = t412 * t331;
t357 = pkin(7) * t382 + qJD(4) * t333 - qJD(5) * t307 - t243;
t219 = -qJD(5) * t410 + t384;
t351 = -qJD(1) * t312 + t373;
t266 = -qJD(3) * pkin(3) + qJD(4) - t279;
t350 = -qJ(4) * t394 + (qJD(4) - t266) * t339;
t349 = 0.2e1 * qJD(3) * t312 - qJDD(3) * t315;
t346 = -t339 * t373 - t414;
t220 = -qJD(5) * t363 + t374;
t341 = qJD(3) ^ 2;
t344 = -0.2e1 * qJDD(1) * t316 - t315 * t341 + t372;
t327 = pkin(9) + qJ(5);
t323 = cos(t327);
t321 = sin(t327);
t317 = -pkin(4) * t333 - pkin(3);
t305 = qJDD(3) * t339 - t336 * t341;
t304 = qJDD(3) * t336 + t339 * t341;
t282 = t376 * t336;
t274 = t376 * t393;
t272 = t331 * t283;
t263 = t321 * t322 + t323 * t408;
t262 = -t321 * t408 + t322 * t323;
t261 = t321 * t324 - t323 * t409;
t260 = t321 * t409 + t323 * t324;
t253 = pkin(4) * t382 + t280;
t251 = -t315 * t406 + t278;
t246 = -t333 * t381 + t272;
t244 = pkin(4) * t292 + t266;
t240 = t363 * t394;
t235 = t272 + (-pkin(7) * t406 - t315 * t405) * qJD(3);
t229 = qJD(3) * t362 + t245;
t222 = pkin(4) * t264 + t234;
t1 = [qJDD(1) * MDP(1) + t371 * MDP(2) + (g(1) * t340 + g(2) * t337) * MDP(3) + ((t332 ^ 2 + t334 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t360) * MDP(4) + (qJDD(1) * t329 + 0.2e1 * t336 * t379) * MDP(5) + 0.2e1 * (t336 * t385 - t389 * t399) * MDP(6) + t304 * MDP(7) + t305 * MDP(8) + (t336 * t349 + t339 * t344) * MDP(10) + (-t336 * t344 + t339 * t349) * MDP(11) + (-t373 * t331 + (t234 * t331 + t315 * t264 + (qJD(1) * t251 + t225) * qJD(3)) * t336 + (-t245 * qJD(1) - t251 * qJDD(1) - t217 + t372 * t333 + (t266 * t331 + t292 * t315) * qJD(3)) * t339) * MDP(12) + (-t373 * t333 + (t234 * t333 + t315 * t265 + (-qJD(1) * t252 - t226) * qJD(3)) * t336 + (t246 * qJD(1) + t252 * qJDD(1) + t218 - t372 * t331 + (t266 * t333 + t294 * t315) * qJD(3)) * t339) * MDP(13) + (-t245 * t294 - t246 * t292 - t251 * t265 - t252 * t264 + (-t225 * t333 - t226 * t331) * t393 + (-t217 * t333 - t218 * t331 + t372) * t336) * MDP(14) + (t217 * t251 + t218 * t252 + t225 * t245 + t226 * t246 + t360 + (-g(1) * pkin(6) - g(2) * t361) * t324 + (-g(2) * pkin(6) + g(1) * t361) * t322 + (t234 * t336 + t266 * t393) * t315) * MDP(15) + (-t219 * t276 - t238 * t363) * MDP(16) + (-t219 * t275 + t220 * t276 - t238 * t247 + t239 * t363) * MDP(17) + (-t219 * t339 - t240 - t364) * MDP(18) + (t220 * t339 - t247 * t394 + t402) * MDP(19) + (-t296 * t339 - t313 * t394) * MDP(20) + (-(t229 * t338 - t235 * t335) * t313 + t366 * t296 - t375 * t339 + t212 * t394 + t274 * t247 + t282 * t220 + t222 * t275 + t244 * t239 - g(1) * t261 - g(2) * t263 + (t213 * t339 + t313 * t365) * qJD(5)) * MDP(21) + ((t229 * t335 + t235 * t338) * t313 - t365 * t296 + t369 * t339 - t213 * t394 - t274 * t363 + t282 * t219 - t222 * t276 + t244 * t238 - g(1) * t260 - g(2) * t262 + (t212 * t339 + t313 * t366) * qJD(5)) * MDP(22); t403 * MDP(4) + t305 * MDP(10) - t304 * MDP(11) - g(3) * MDP(15) + t402 * MDP(21) + (-t240 + t364) * MDP(22) + ((-t264 * t333 + t265 * t331) * MDP(14) + t368 * MDP(15)) * t336 + ((-t264 + t378) * MDP(12) + (-t265 + t377) * MDP(13) - t234 * MDP(15) - t220 * MDP(21) - t219 * MDP(22)) * t339 + ((t292 * MDP(12) + t294 * MDP(13) + t266 * MDP(15) + t247 * MDP(21)) * t336 + ((-t292 * t333 + t294 * t331) * MDP(14) + t367 * MDP(15)) * t339 + (-MDP(12) * t331 - MDP(13) * t333) * qJD(1) * t399) * qJD(3); MDP(7) * t387 + MDP(8) * t385 + qJDD(3) * MDP(9) + (qJD(3) * t280 + t351 * t336 + t383 + t416) * MDP(10) + (qJD(3) * t279 + (qJD(3) * t311 - t403) * t336 + (t351 + t417) * t339) * MDP(11) + (t331 * t380 - pkin(3) * t264 - t280 * t292 - t419 * t333 + (-t225 * t336 + t242 * t339 + t331 * t350) * qJD(1)) * MDP(12) + (t333 * t380 - pkin(3) * t265 - t280 * t294 + t419 * t331 + (t226 * t336 - t243 * t339 + t333 * t350) * qJD(1)) * MDP(13) + (t242 * t294 + t243 * t292 + (-qJ(4) * t264 - qJD(4) * t292 + t225 * t397 + t218) * t333 + (qJ(4) * t265 + qJD(4) * t294 + t226 * t397 - t217) * t331 + t346) * MDP(14) + (-t225 * t242 - t226 * t243 - t266 * t280 + t367 * qJD(4) - t419 * pkin(3) + (t346 + t368) * qJ(4)) * MDP(15) + (t219 * t298 - t363 * t401) * MDP(16) + (-t219 * t297 - t220 * t298 - t401 * t247 + t363 * t400) * MDP(17) + (t296 * t298 - t401 * t313 + t363 * t398) * MDP(18) + (t247 * t398 - t296 * t297 + t400 * t313) * MDP(19) + t313 * MDP(20) * t398 + ((-t307 * t338 - t308 * t335) * t296 + t317 * t220 + t222 * t297 - t212 * t398 - t253 * t247 + (t335 * t357 + t338 * t358) * t313 + t400 * t244 + t348 * t323) * MDP(21) + (-(-t307 * t335 + t308 * t338) * t296 + t317 * t219 + t222 * t298 + t213 * t398 + t253 * t363 + (-t335 * t358 + t338 * t357) * t313 + t401 * t244 - t348 * t321) * MDP(22) + (-t336 * t339 * MDP(5) + t399 * MDP(6)) * qJD(1) ^ 2; (t378 - t319 + (-t294 + t395) * t397) * MDP(12) + (t377 + t386 + (t292 + t390) * t397) * MDP(13) + (-t292 ^ 2 - t294 ^ 2) * MDP(14) + (t225 * t294 + t226 * t292 + t352 - t359 - t416) * MDP(15) + (t220 + t418) * MDP(21) + (t219 + t411) * MDP(22); -t363 * t247 * MDP(16) + (-t247 ^ 2 + t363 ^ 2) * MDP(17) + (t384 - t411) * MDP(18) + (-t374 + t418) * MDP(19) + t296 * MDP(20) + (-g(1) * t262 + g(2) * t260 - t213 * t313 + t244 * t363 + t321 * t414 + t375) * MDP(21) + (g(1) * t263 - g(2) * t261 - t212 * t313 + t244 * t247 + t323 * t414 - t369) * MDP(22) + (-MDP(18) * t410 + MDP(19) * t363 - t213 * MDP(21) - t212 * MDP(22)) * qJD(5);];
tau = t1;
