% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:21
% EndTime: 2021-01-15 11:44:28
% DurationCPUTime: 2.52s
% Computational Cost: add. (1688->287), mult. (3589->380), div. (0->0), fcn. (2526->16), ass. (0->140)
t306 = sin(pkin(9));
t308 = cos(pkin(9));
t315 = cos(qJ(3));
t350 = qJD(1) * t315;
t342 = t308 * t350;
t312 = sin(qJ(3));
t351 = qJD(1) * t312;
t261 = t306 * t351 - t342;
t314 = cos(qJ(5));
t253 = t314 * t261;
t270 = t306 * t315 + t308 * t312;
t264 = t270 * qJD(1);
t311 = sin(qJ(5));
t359 = t264 * t311;
t222 = -t253 - t359;
t301 = qJD(3) + qJD(5);
t360 = t222 * t301;
t329 = t261 * t311 - t314 * t264;
t361 = t329 * t301;
t307 = sin(pkin(8));
t288 = pkin(1) * t307 + pkin(6);
t355 = qJ(4) + t288;
t303 = qJ(1) + pkin(8);
t294 = sin(t303);
t296 = cos(t303);
t335 = g(2) * t294 - g(3) * t296;
t297 = t315 * qJDD(2);
t276 = t288 * qJDD(1);
t323 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t276;
t337 = t355 * qJD(1);
t328 = t337 * qJD(3);
t208 = qJDD(3) * pkin(3) - t312 * t323 - t315 * t328 + t297;
t212 = (qJDD(2) - t328) * t312 + t323 * t315;
t191 = t308 * t208 - t212 * t306;
t348 = qJD(1) * qJD(3);
t341 = t312 * t348;
t280 = t306 * t341;
t340 = t315 * t348;
t230 = t270 * qJDD(1) + t308 * t340 - t280;
t189 = qJDD(3) * pkin(4) - pkin(7) * t230 + t191;
t192 = t306 * t208 + t308 * t212;
t263 = t270 * qJD(3);
t346 = qJDD(1) * t315;
t281 = t308 * t346;
t347 = qJDD(1) * t312;
t229 = qJD(1) * t263 + t306 * t347 - t281;
t190 = -pkin(7) * t229 + t192;
t246 = t315 * qJD(2) - t312 * t337;
t362 = qJD(3) * pkin(3);
t240 = t246 + t362;
t247 = qJD(2) * t312 + t315 * t337;
t357 = t308 * t247;
t211 = t306 * t240 + t357;
t367 = pkin(7) * t261;
t198 = t211 - t367;
t292 = pkin(3) * t315 + pkin(2);
t309 = cos(pkin(8));
t370 = pkin(1) * t309;
t273 = -t292 - t370;
t259 = qJD(1) * t273 + qJD(4);
t231 = pkin(4) * t261 + t259;
t302 = qJ(3) + pkin(9);
t298 = qJ(5) + t302;
t286 = sin(t298);
t287 = cos(t298);
t349 = qJD(5) * t311;
t373 = g(1) * t286 - t311 * t189 - t314 * t190 + t198 * t349 - t231 * t222 + t335 * t287;
t372 = -g(1) * t287 + t314 * t189 - t311 * t190 + t231 * t329 + t335 * t286;
t300 = qJDD(3) + qJDD(5);
t371 = t300 * MDP(20) + t222 * MDP(16) * t329 + (-t222 ^ 2 + t329 ^ 2) * MDP(17);
t339 = t314 * t229 + t230 * t311;
t195 = -qJD(5) * t329 + t339;
t369 = pkin(3) * t306;
t368 = pkin(3) * t312;
t366 = pkin(7) * t264;
t365 = g(1) * t315;
t236 = t306 * t247;
t358 = t306 * t312;
t210 = t308 * t240 - t236;
t197 = qJD(3) * pkin(4) + t210 - t366;
t356 = t314 * t197;
t354 = qJDD(2) - g(1);
t214 = t308 * t246 - t236;
t338 = qJD(3) * t355;
t250 = qJD(4) * t315 - t312 * t338;
t251 = -qJD(4) * t312 - t315 * t338;
t216 = t308 * t250 + t306 * t251;
t267 = t355 * t312;
t268 = t355 * t315;
t228 = -t306 * t267 + t308 * t268;
t304 = t312 ^ 2;
t353 = -t315 ^ 2 + t304;
t290 = -pkin(2) - t370;
t279 = qJD(1) * t290;
t345 = pkin(3) * t341 + qJDD(4);
t344 = t312 * t362;
t343 = -qJD(5) * t253 - t311 * t229 + t314 * t230;
t213 = -t246 * t306 - t357;
t215 = -t250 * t306 + t308 * t251;
t227 = -t308 * t267 - t268 * t306;
t336 = g(2) * t296 + g(3) * t294;
t313 = sin(qJ(1));
t316 = cos(qJ(1));
t334 = -g(2) * t316 - g(3) * t313;
t333 = -t311 * t197 - t314 * t198;
t269 = -t308 * t315 + t358;
t232 = t314 * t269 + t270 * t311;
t266 = t269 * qJD(3);
t199 = -qJD(5) * t232 - t263 * t311 - t266 * t314;
t233 = -t269 * t311 + t270 * t314;
t332 = t199 * t301 + t233 * t300;
t217 = -pkin(7) * t270 + t227;
t218 = -pkin(7) * t269 + t228;
t331 = t217 * t314 - t218 * t311;
t330 = t217 * t311 + t218 * t314;
t289 = pkin(3) * t308 + pkin(4);
t327 = t289 * t311 + t314 * t369;
t326 = t289 * t314 - t311 * t369;
t194 = -t264 * t349 + t343;
t325 = -qJD(1) * t279 - t276 + t335;
t324 = 0.2e1 * qJD(3) * t279 - qJDD(3) * t288;
t245 = qJDD(1) * t273 + t345;
t317 = qJD(3) ^ 2;
t321 = 0.2e1 * qJDD(1) * t290 + t288 * t317 + t336;
t310 = -qJ(4) - pkin(6);
t295 = cos(t302);
t293 = sin(t302);
t275 = qJDD(3) * t315 - t312 * t317;
t274 = qJDD(3) * t312 + t315 * t317;
t249 = pkin(4) * t263 + t344;
t248 = pkin(3) * t351 + pkin(4) * t264;
t244 = pkin(4) * t269 + t273;
t209 = pkin(4) * t229 + t245;
t204 = -pkin(7) * t263 + t216;
t203 = pkin(7) * t266 + t215;
t202 = t214 - t366;
t201 = t213 + t367;
t200 = qJD(5) * t233 + t314 * t263 - t266 * t311;
t193 = -t200 * t301 - t232 * t300;
t1 = [qJDD(1) * MDP(1) + t334 * MDP(2) + (g(2) * t313 - g(3) * t316) * MDP(3) + (t334 + (t307 ^ 2 + t309 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t304 + 0.2e1 * t312 * t340) * MDP(5) + 0.2e1 * (t312 * t346 - t348 * t353) * MDP(6) + t274 * MDP(7) + t275 * MDP(8) + (t312 * t324 - t315 * t321) * MDP(10) + (t312 * t321 + t315 * t324) * MDP(11) + (qJDD(3) * t227 + t229 * t273 + t245 * t269 + t259 * t263 - t336 * t295 + (t261 * t368 + t215) * qJD(3)) * MDP(12) + (-qJDD(3) * t228 + t230 * t273 + t245 * t270 - t259 * t266 + t336 * t293 + (t264 * t368 - t216) * qJD(3)) * MDP(13) + (-t191 * t270 - t192 * t269 + t210 * t266 - t211 * t263 - t215 * t264 - t216 * t261 - t227 * t230 - t228 * t229 - t335) * MDP(14) + (t192 * t228 + t211 * t216 + t191 * t227 + t210 * t215 + t245 * t273 + t259 * t344 - g(2) * (pkin(1) * t316 + t292 * t296 - t294 * t310) - g(3) * (pkin(1) * t313 + t292 * t294 + t296 * t310)) * MDP(15) + (t194 * t233 - t199 * t329) * MDP(16) + (-t194 * t232 - t195 * t233 + t199 * t222 + t200 * t329) * MDP(17) + t332 * MDP(18) + t193 * MDP(19) + (-t249 * t222 + t244 * t195 + t209 * t232 + t231 * t200 + (-qJD(5) * t330 + t203 * t314 - t204 * t311) * t301 + t331 * t300 - t336 * t287) * MDP(21) + (-t249 * t329 + t244 * t194 + t209 * t233 + t231 * t199 - (qJD(5) * t331 + t203 * t311 + t204 * t314) * t301 - t330 * t300 + t336 * t286) * MDP(22); t354 * MDP(4) + t275 * MDP(10) - t274 * MDP(11) + (-qJD(3) * t263 - qJDD(3) * t269) * MDP(12) + (qJD(3) * t266 - qJDD(3) * t270) * MDP(13) + (-t229 * t270 + t230 * t269 + t261 * t266 + t263 * t264) * MDP(14) + (-t191 * t269 + t192 * t270 - t210 * t263 - t211 * t266 - g(1)) * MDP(15) + t193 * MDP(21) - t332 * MDP(22); MDP(7) * t347 + MDP(8) * t346 + qJDD(3) * MDP(9) + (t312 * t325 + t297 - t365) * MDP(10) + (-t312 * t354 + t315 * t325) * MDP(11) + (-g(1) * t295 - qJD(3) * t213 - t259 * t264 + t335 * t293 + (qJDD(3) * t308 - t261 * t351) * pkin(3) + t191) * MDP(12) + (g(1) * t293 + qJD(3) * t214 + t259 * t261 + t335 * t295 + (-qJDD(3) * t306 - t264 * t351) * pkin(3) - t192) * MDP(13) + ((t211 + t213) * t264 + (-t210 + t214) * t261 + (-t229 * t306 - t230 * t308) * pkin(3)) * MDP(14) + (-t210 * t213 - t211 * t214 + (-t365 + t191 * t308 + t192 * t306 + (-qJD(1) * t259 + t335) * t312) * pkin(3)) * MDP(15) + (t194 - t360) * MDP(18) + (-t195 - t361) * MDP(19) + (t326 * t300 + t248 * t222 - (t201 * t314 - t202 * t311) * t301 + (-t301 * t327 + t333) * qJD(5) + t372) * MDP(21) + (-t327 * t300 + t248 * t329 + (t201 * t311 + t202 * t314) * t301 + (-t301 * t326 - t356) * qJD(5) + t373) * MDP(22) + (-MDP(5) * t312 * t315 + MDP(6) * t353) * qJD(1) ^ 2 + t371; -t281 * MDP(12) - t280 * MDP(13) + (-t261 ^ 2 - t264 ^ 2) * MDP(14) + (t210 * t264 + t211 * t261 + t336 + t345) * MDP(15) + (t195 - t361) * MDP(21) + (t194 + t360) * MDP(22) + (MDP(12) * t358 + t270 * MDP(13) + t273 * MDP(15)) * qJDD(1) + ((t306 * t350 + t308 * t351 + t264) * MDP(12) + (-t261 + t342) * MDP(13)) * qJD(3); (t343 - t360) * MDP(18) + (-t339 - t361) * MDP(19) + (-t301 * t333 + t372) * MDP(21) + ((-t198 * t311 + t356) * t301 + t373) * MDP(22) + (-MDP(18) * t359 + t329 * MDP(19) + t333 * MDP(21) - MDP(22) * t356) * qJD(5) + t371;];
tau = t1;
