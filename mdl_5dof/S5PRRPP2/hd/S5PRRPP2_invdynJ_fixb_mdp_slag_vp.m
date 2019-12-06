% Calculate vector of inverse dynamics joint torques for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:17
% EndTime: 2019-12-05 16:10:21
% DurationCPUTime: 2.34s
% Computational Cost: add. (1385->287), mult. (3019->364), div. (0->0), fcn. (2065->10), ass. (0->133)
t347 = qJ(4) + pkin(6);
t274 = sin(qJ(3));
t276 = cos(qJ(3));
t309 = qJD(3) * t347;
t228 = qJD(4) * t276 - t274 * t309;
t269 = sin(pkin(8));
t271 = cos(pkin(8));
t289 = -qJD(4) * t274 - t276 * t309;
t340 = t271 * t276;
t342 = t269 * t274;
t240 = -t340 + t342;
t277 = cos(qJ(2));
t293 = t240 * t277;
t333 = qJD(1) * t293 + t271 * t228 + t269 * t289;
t270 = sin(pkin(7));
t272 = cos(pkin(7));
t305 = g(1) * t272 + g(2) * t270;
t295 = t305 * t277;
t275 = sin(qJ(2));
t348 = g(3) * t275;
t286 = t295 + t348;
t357 = MDP(12) + MDP(15);
t296 = t305 * t275;
t264 = g(3) * t277;
t285 = -t264 + t296;
t327 = qJD(3) * t274;
t355 = -qJD(3) * t340 + t269 * t327;
t328 = qJD(2) * t276;
t329 = qJD(2) * t274;
t354 = -t269 * t328 - t271 * t329;
t241 = t269 * t276 + t271 * t274;
t233 = t241 * qJD(3);
t353 = -MDP(13) - MDP(17);
t278 = qJD(3) ^ 2;
t323 = qJDD(1) * t277;
t325 = qJD(1) * qJD(2);
t306 = t275 * t325 - t323;
t352 = (2 * qJDD(2) * pkin(2)) - pkin(6) * t278 + t275 * (t305 + t325) - t264 - t306;
t234 = t241 * qJD(2);
t230 = t234 ^ 2;
t351 = pkin(3) * t276;
t346 = qJD(2) * pkin(2);
t344 = qJDD(3) * pkin(4);
t331 = qJD(1) * t275;
t308 = t347 * qJD(2) + t331;
t227 = t308 * t276;
t343 = t227 * t269;
t341 = t270 * t277;
t214 = t271 * t227;
t339 = t272 * t276;
t338 = t272 * t277;
t337 = t274 * t277;
t279 = qJD(2) ^ 2;
t336 = t276 * t279;
t335 = qJDD(1) - g(3);
t239 = qJDD(2) * pkin(6) + qJDD(1) * t275 + t277 * t325;
t292 = (qJ(4) * qJDD(2)) + qJD(2) * qJD(4) + t239;
t300 = qJD(3) * t308;
t192 = qJDD(3) * pkin(3) - t274 * t292 - t276 * t300;
t195 = -t274 * t300 + t276 * t292;
t182 = t271 * t192 - t269 * t195;
t183 = t269 * t192 + t271 * t195;
t330 = qJD(1) * t277;
t334 = t228 * t269 - t241 * t330 - t271 * t289;
t226 = t308 * t274;
t216 = qJD(3) * pkin(3) - t226;
t199 = t269 * t216 + t214;
t267 = t274 ^ 2;
t332 = -t276 ^ 2 + t267;
t201 = -t226 * t271 - t343;
t326 = qJD(5) - t201;
t324 = qJD(2) * qJD(3);
t322 = qJDD(2) * t274;
t321 = qJDD(2) * t276;
t320 = pkin(3) * t327;
t319 = qJDD(3) * qJ(5) + t183;
t261 = pkin(2) + t351;
t314 = t271 * t328;
t313 = t347 * t274;
t312 = t274 * t324;
t311 = t276 * t324;
t310 = -qJDD(5) + t182;
t189 = pkin(4) * t233 + qJ(5) * t355 - qJD(5) * t241 + t320;
t307 = -t189 + t331;
t304 = g(1) * t270 - g(2) * t272;
t303 = t261 * qJDD(2);
t265 = qJ(3) + pkin(8);
t262 = sin(t265);
t263 = cos(t265);
t302 = pkin(4) * t263 + qJ(5) * t262;
t198 = t216 * t271 - t343;
t299 = qJDD(3) * t276 - t274 * t278;
t298 = qJDD(3) * t274 + t276 * t278;
t294 = pkin(3) * t312 + qJDD(4) + t306;
t247 = t269 * t312;
t287 = t241 * qJDD(2) - t247;
t237 = -qJD(2) * t261 + qJD(4) - t330;
t251 = -t330 - t346;
t284 = -pkin(6) * qJDD(3) + (t251 + t330 - t346) * qJD(3);
t283 = -t251 * qJD(2) - t239 + t286;
t252 = t271 * t321;
t208 = qJD(2) * t233 + t269 * t322 - t252;
t209 = t271 * t311 + t287;
t282 = pkin(4) * t208 - qJ(5) * t209 - qJD(5) * t234 + t294;
t231 = t269 * t329 - t314;
t196 = pkin(4) * t231 - qJ(5) * t234 + t237;
t220 = t262 * t341 + t263 * t272;
t222 = t262 * t338 - t270 * t263;
t281 = g(1) * t222 + g(2) * t220 - t196 * t234 + t262 * t348 + t310;
t246 = t347 * t276;
t211 = t246 * t269 + t271 * t313;
t212 = t271 * t246 - t269 * t313;
t280 = -t212 * t208 + t209 * t211 - t333 * t231 + t234 * t334 - t286;
t259 = -pkin(3) * t271 - pkin(4);
t257 = pkin(3) * t269 + qJ(5);
t255 = t270 * t351;
t248 = t277 * t261;
t225 = t240 * t275;
t224 = t241 * t275;
t223 = t262 * t270 + t263 * t338;
t221 = -t262 * t272 + t263 * t341;
t210 = -t303 + t294;
t207 = pkin(4) * t240 - qJ(5) * t241 - t261;
t204 = pkin(3) * t329 + pkin(4) * t234 + qJ(5) * t231;
t203 = -qJD(2) * t293 - t233 * t275;
t202 = t275 * t355 + t277 * t354;
t200 = -t226 * t269 + t214;
t197 = qJD(3) * qJ(5) + t199;
t194 = -qJD(3) * pkin(4) + qJD(5) - t198;
t181 = -t310 - t344;
t180 = -t303 + t282;
t179 = qJD(3) * qJD(5) + t319;
t1 = [t335 * MDP(1) + (-t182 * t224 - t183 * t225 + t198 * t202 + t199 * t203 - g(3)) * MDP(13) + (qJD(3) * t202 - qJDD(3) * t224) * MDP(14) + (qJD(3) * t203 - qJDD(3) * t225) * MDP(16) + (-t179 * t225 + t181 * t224 - t194 * t202 + t197 * t203 - g(3)) * MDP(17) + ((qJDD(2) * MDP(3)) - t279 * MDP(4) + (-0.2e1 * t312 + t321) * MDP(10) + (-0.2e1 * t311 - t322) * MDP(11) - t210 * MDP(13) - t208 * MDP(14) + t209 * MDP(16) - t180 * MDP(17)) * t277 + (-t279 * MDP(3) - qJDD(2) * MDP(4) + (-t298 - t336) * MDP(10) + (t274 * t279 - t299) * MDP(11) + (t237 * MDP(13) + t231 * MDP(14) - t234 * MDP(16) + t196 * MDP(17)) * qJD(2)) * t275 + t357 * (-t202 * t234 - t203 * t231 + t225 * t208 + t209 * t224); (qJDD(2) * MDP(2)) + (t323 + t285) * MDP(3) + (-t275 * t335 + t295) * MDP(4) + (qJDD(2) * t267 + 0.2e1 * t274 * t311) * MDP(5) + 0.2e1 * (t274 * t321 - t324 * t332) * MDP(6) + t298 * MDP(7) + t299 * MDP(8) + (t284 * t274 + t276 * t352) * MDP(10) + (-t274 * t352 + t284 * t276) * MDP(11) + (-t182 * t241 - t183 * t240 + t198 * t355 - t199 * t233 + t280) * MDP(12) + (t183 * t212 - t182 * t211 - t210 * t261 - g(3) * (t275 * t347 + t248) + (t320 - t331) * t237 + t333 * t199 - t334 * t198 + t305 * (t261 * t275 - t277 * t347)) * MDP(13) + (-qJD(3) * t334 - qJDD(3) * t211 + t180 * t240 + t196 * t233 + t207 * t208 - t231 * t307 + t263 * t285) * MDP(14) + (-t179 * t240 + t181 * t241 - t194 * t355 - t197 * t233 + t280) * MDP(15) + (qJD(3) * t333 + qJDD(3) * t212 - t180 * t241 + t196 * t355 - t207 * t209 + t234 * t307 + t262 * t285) * MDP(16) + (-g(3) * t248 + t179 * t212 + t180 * t207 + t181 * t211 + t196 * t189 + t333 * t197 + t334 * t194 + (-g(3) * t302 - t305 * t347) * t277 + (-g(3) * t347 - t196 * qJD(1) + t305 * (t261 + t302)) * t275) * MDP(17); -t274 * MDP(5) * t336 + t332 * MDP(6) * t279 + MDP(7) * t322 + MDP(8) * t321 + qJDD(3) * MDP(9) + (t274 * t283 - t276 * t304) * MDP(10) + (t274 * t304 + t276 * t283) * MDP(11) + ((t199 - t200) * t234 + (-t198 + t201) * t231 + (-t208 * t269 - t209 * t271) * pkin(3)) * MDP(12) + (-g(1) * t255 + t198 * t200 - t199 * t201 + (g(2) * t339 + t182 * t271 + t183 * t269 + (-t237 * qJD(2) + t286) * t274) * pkin(3)) * MDP(13) + (qJD(3) * t200 - t204 * t231 + (pkin(4) - t259) * qJDD(3) + t281) * MDP(14) + (-t208 * t257 + t209 * t259 + (t197 - t200) * t234 + (t194 - t326) * t231) * MDP(15) + (-t263 * t348 - g(1) * t223 - g(2) * t221 + qJDD(3) * t257 - t196 * t231 + t204 * t234 + (0.2e1 * qJD(5) - t201) * qJD(3) + t319) * MDP(16) + (t179 * t257 + t181 * t259 - t196 * t204 - t194 * t200 - g(1) * (-pkin(3) * t272 * t337 - pkin(4) * t222 + qJ(5) * t223 + t255) - g(2) * (-pkin(4) * t220 + qJ(5) * t221 + (-t270 * t337 - t339) * pkin(3)) + t326 * t197 - (-pkin(3) * t274 - pkin(4) * t262 + qJ(5) * t263) * t348) * MDP(17); (t198 * t234 + t199 * t231 + t264 + t294) * MDP(13) - t252 * MDP(14) + t247 * MDP(16) + (-t194 * t234 + t197 * t231 + t264 + t282) * MDP(17) + (MDP(14) * t342 - MDP(16) * t241 + t261 * t353) * qJDD(2) + ((t234 - t354) * MDP(14) + (t231 - t314) * MDP(16)) * qJD(3) + t353 * t296 + t357 * (-t231 ^ 2 - t230); (t231 * t234 - qJDD(3)) * MDP(14) + ((t231 + t314) * qJD(3) + t287) * MDP(15) + (-t230 - t278) * MDP(16) + (-qJD(3) * t197 - t281 - t344) * MDP(17);];
tau = t1;
