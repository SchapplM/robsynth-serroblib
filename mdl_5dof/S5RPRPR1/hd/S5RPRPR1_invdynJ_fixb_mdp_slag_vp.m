% Calculate vector of inverse dynamics joint torques for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:34:00
% EndTime: 2021-01-15 11:34:08
% DurationCPUTime: 2.81s
% Computational Cost: add. (1649->296), mult. (3201->375), div. (0->0), fcn. (2180->12), ass. (0->137)
t332 = sin(qJ(1));
t335 = cos(qJ(1));
t398 = g(1) * t332 - g(2) * t335;
t406 = qJDD(2) - t398;
t328 = sin(pkin(8));
t329 = cos(pkin(8));
t331 = sin(qJ(3));
t334 = cos(qJ(3));
t284 = t328 * t334 + t329 * t331;
t274 = t284 * qJD(1);
t333 = cos(qJ(5));
t379 = t333 * t274;
t372 = qJD(1) * t331;
t360 = t328 * t372;
t371 = qJD(1) * t334;
t278 = t329 * t371 - t360;
t330 = sin(qJ(5));
t383 = t278 * t330;
t235 = t379 + t383;
t321 = qJD(3) + qJD(5);
t385 = t235 * t321;
t367 = qJD(1) * qJD(3);
t358 = t334 * t367;
t365 = qJDD(1) * t331;
t405 = t358 + t365;
t323 = qJDD(1) * qJ(2);
t324 = qJD(1) * qJD(2);
t355 = g(1) * t335 + g(2) * t332;
t343 = -t355 + 0.2e1 * t324;
t404 = 0.2e1 * t323 + t343;
t388 = qJDD(1) * pkin(1);
t403 = t388 - t406;
t336 = -pkin(1) - pkin(6);
t293 = t336 * qJDD(1) + qJDD(2);
t285 = t334 * t293;
t294 = t336 * qJD(1) + qJD(2);
t359 = t331 * t367;
t364 = qJDD(1) * t334;
t366 = qJD(1) * qJD(4);
t370 = qJD(3) * t331;
t224 = -t334 * t366 - t294 * t370 + qJDD(3) * pkin(3) + t285 + (t359 - t364) * qJ(4);
t369 = qJD(3) * t334;
t229 = (-qJ(4) * qJD(1) + t294) * t369 + (-qJ(4) * qJDD(1) + t293 - t366) * t331;
t209 = t329 * t224 - t229 * t328;
t353 = -t328 * t365 + t329 * t364;
t243 = t284 * t367 - t353;
t205 = qJDD(3) * pkin(4) + pkin(7) * t243 + t209;
t210 = t328 * t224 + t329 * t229;
t361 = t328 * t364 + t329 * t405;
t242 = t328 * t359 - t361;
t206 = pkin(7) * t242 + t210;
t271 = -qJ(4) * t371 + t334 * t294;
t261 = qJD(3) * pkin(3) + t271;
t270 = -qJ(4) * t372 + t294 * t331;
t381 = t329 * t270;
t226 = t328 * t261 + t381;
t392 = pkin(7) * t274;
t215 = t226 - t392;
t287 = pkin(3) * t372 + qJD(1) * qJ(2) + qJD(4);
t248 = pkin(4) * t274 + t287;
t322 = qJ(3) + pkin(8);
t316 = qJ(5) + t322;
t304 = sin(t316);
t305 = cos(t316);
t368 = qJD(5) * t330;
t402 = g(3) * t305 - t330 * t205 - t333 * t206 + t215 * t368 + t248 * t235 + t304 * t398;
t319 = qJDD(3) + qJDD(5);
t348 = -t274 * t330 + t278 * t333;
t401 = t319 * MDP(22) + t235 * MDP(18) * t348 + (-t235 ^ 2 + t348 ^ 2) * MDP(19);
t386 = t348 * t321;
t399 = qJ(4) - t336;
t283 = -t328 * t331 + t329 * t334;
t245 = t283 * t333 - t284 * t330;
t244 = t283 * t330 + t284 * t333;
t276 = t328 * t370 - t329 * t369;
t277 = t284 * qJD(3);
t211 = -qJD(5) * t244 + t276 * t330 - t277 * t333;
t397 = t211 * t321 + t245 * t319;
t396 = g(3) * t304 + t333 * t205 - t330 * t206 - t248 * t348 - t305 * t398;
t357 = t242 * t333 + t243 * t330;
t208 = qJD(5) * t348 - t357;
t393 = pkin(3) * t328;
t391 = pkin(7) * t278;
t390 = g(3) * t331;
t338 = qJD(1) ^ 2;
t389 = qJ(2) * t338;
t257 = t328 * t270;
t225 = t261 * t329 - t257;
t214 = qJD(3) * pkin(4) + t225 - t391;
t380 = t333 * t214;
t268 = -qJD(4) * t334 + t370 * t399;
t289 = t399 * t334;
t269 = -qJD(3) * t289 - qJD(4) * t331;
t228 = t268 * t328 + t269 * t329;
t233 = t271 * t329 - t257;
t288 = t399 * t331;
t247 = -t288 * t329 - t289 * t328;
t327 = t334 ^ 2;
t376 = t331 ^ 2 - t327;
t337 = qJD(3) ^ 2;
t375 = -t337 - t338;
t373 = qJD(1) * t287;
t298 = pkin(3) * t369 + qJD(2);
t363 = qJDD(3) * t331;
t362 = -qJD(5) * t379 + t242 * t330 - t243 * t333;
t306 = pkin(3) * t331 + qJ(2);
t227 = t268 * t329 - t269 * t328;
t232 = -t271 * t328 - t381;
t246 = t288 * t328 - t289 * t329;
t212 = qJD(5) * t245 - t333 * t276 - t277 * t330;
t352 = -t212 * t321 - t244 * t319;
t351 = -t330 * t214 - t333 * t215;
t230 = -pkin(7) * t283 + t246;
t231 = -pkin(7) * t284 + t247;
t350 = t230 * t333 - t231 * t330;
t349 = t230 * t330 + t231 * t333;
t256 = pkin(3) * t405 + qJDD(4) + t323 + t324;
t307 = pkin(3) * t329 + pkin(4);
t347 = t307 * t330 + t333 * t393;
t346 = t307 * t333 - t330 * t393;
t345 = 0.2e1 * qJ(2) * t367 + qJDD(3) * t336;
t207 = -t278 * t368 + t362;
t344 = -t398 - t389;
t340 = t209 * t283 + t210 * t284 - t225 * t277 - t226 * t276 - t398;
t339 = -t336 * t337 + t404;
t315 = qJDD(3) * t334;
t314 = cos(t322);
t313 = sin(t322);
t262 = pkin(4) * t284 + t306;
t252 = pkin(3) * t371 + pkin(4) * t278;
t249 = -pkin(4) * t276 + t298;
t220 = -pkin(4) * t242 + t256;
t219 = t233 - t391;
t218 = t232 + t392;
t217 = pkin(7) * t276 + t228;
t216 = pkin(7) * t277 + t227;
t1 = [qJDD(1) * MDP(1) + t398 * MDP(2) + t355 * MDP(3) + (-0.2e1 * t388 + t406) * MDP(4) + t404 * MDP(5) + (t403 * pkin(1) + (t343 + t323) * qJ(2)) * MDP(6) + (qJDD(1) * t327 - 0.2e1 * t331 * t358) * MDP(7) + 0.2e1 * (-t331 * t364 + t367 * t376) * MDP(8) + (-t331 * t337 + t315) * MDP(9) + (-t334 * t337 - t363) * MDP(10) + (t331 * t339 + t334 * t345) * MDP(12) + (-t331 * t345 + t334 * t339) * MDP(13) + (qJD(3) * t227 + qJDD(3) * t246 - t242 * t306 + t256 * t284 + t274 * t298 - t276 * t287 - t313 * t355) * MDP(14) + (-qJD(3) * t228 - qJDD(3) * t247 - t243 * t306 + t256 * t283 - t277 * t287 + t278 * t298 - t314 * t355) * MDP(15) + (-t227 * t278 - t228 * t274 + t242 * t247 + t243 * t246 - t340) * MDP(16) + (t210 * t247 + t226 * t228 + t209 * t246 + t225 * t227 + t256 * t306 + t287 * t298 - g(1) * (t306 * t335 - t332 * t399) - g(2) * (t306 * t332 + t335 * t399)) * MDP(17) + (t207 * t245 + t211 * t348) * MDP(18) + (-t207 * t244 - t208 * t245 - t211 * t235 - t212 * t348) * MDP(19) + t397 * MDP(20) + t352 * MDP(21) + (t249 * t235 + t262 * t208 + t220 * t244 + t248 * t212 + (-qJD(5) * t349 + t216 * t333 - t217 * t330) * t321 + t350 * t319 - t355 * t304) * MDP(23) + (t249 * t348 + t262 * t207 + t220 * t245 + t248 * t211 - (qJD(5) * t350 + t216 * t330 + t217 * t333) * t321 - t349 * t319 - t355 * t305) * MDP(24); qJDD(1) * MDP(4) - t338 * MDP(5) + (-t389 - t403) * MDP(6) + (t331 * t375 + t315) * MDP(12) + (t334 * t375 - t363) * MDP(13) + (-qJD(1) * t274 - qJD(3) * t277 + qJDD(3) * t283) * MDP(14) + (-qJD(1) * t278 + qJD(3) * t276 - qJDD(3) * t284) * MDP(15) + (t242 * t284 + t243 * t283 + t274 * t276 + t277 * t278) * MDP(16) + (t340 - t373) * MDP(17) + (-qJD(1) * t235 + t397) * MDP(23) + (-qJD(1) * t348 + t352) * MDP(24); MDP(9) * t364 - MDP(10) * t365 + qJDD(3) * MDP(11) + (t334 * t344 + t285 + t390) * MDP(12) + (g(3) * t334 + (-t293 - t344) * t331) * MDP(13) + (g(3) * t313 - qJD(3) * t232 - t278 * t287 - t398 * t314 + (qJDD(3) * t329 - t274 * t371) * pkin(3) + t209) * MDP(14) + (g(3) * t314 + qJD(3) * t233 + t274 * t287 + t398 * t313 + (-qJDD(3) * t328 - t278 * t371) * pkin(3) - t210) * MDP(15) + ((t226 + t232) * t278 + (-t225 + t233) * t274 + (t242 * t328 + t243 * t329) * pkin(3)) * MDP(16) + (-t225 * t232 - t226 * t233 + (t390 + t209 * t329 + t210 * t328 + (-t398 - t373) * t334) * pkin(3)) * MDP(17) + (t207 + t385) * MDP(20) + (-t208 + t386) * MDP(21) + (t346 * t319 - t252 * t235 - (t218 * t333 - t219 * t330) * t321 + (-t321 * t347 + t351) * qJD(5) + t396) * MDP(23) + (-t347 * t319 - t252 * t348 + (t218 * t330 + t219 * t333) * t321 + (-t321 * t346 - t380) * qJD(5) + t402) * MDP(24) + (MDP(7) * t331 * t334 - MDP(8) * t376) * t338 + t401; t361 * MDP(14) + t353 * MDP(15) + (-t274 ^ 2 - t278 ^ 2) * MDP(16) + (t225 * t278 + t226 * t274 + t256 - t355) * MDP(17) + (t208 + t386) * MDP(23) + (t207 - t385) * MDP(24) + ((t278 - t360) * MDP(14) + (-t328 * t371 - t329 * t372 - t274) * MDP(15)) * qJD(3); (t362 + t385) * MDP(20) + (t357 + t386) * MDP(21) + (-t321 * t351 + t396) * MDP(23) + ((-t215 * t330 + t380) * t321 + t402) * MDP(24) + (-MDP(20) * t383 - t348 * MDP(21) + t351 * MDP(23) - MDP(24) * t380) * qJD(5) + t401;];
tau = t1;
