% Calculate vector of inverse dynamics joint torques for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:20
% EndTime: 2019-12-31 17:26:24
% DurationCPUTime: 3.03s
% Computational Cost: add. (1666->275), mult. (3782->378), div. (0->0), fcn. (2676->10), ass. (0->132)
t280 = qJD(2) + qJD(3);
t285 = sin(qJ(3));
t289 = cos(qJ(2));
t360 = cos(qJ(3));
t328 = qJD(1) * t360;
t286 = sin(qJ(2));
t340 = qJD(1) * t286;
t367 = -t285 * t340 + t289 * t328;
t369 = t280 * t367;
t283 = qJ(2) + qJ(3);
t277 = sin(t283);
t287 = sin(qJ(1));
t290 = cos(qJ(1));
t314 = g(1) * t290 + g(2) * t287;
t368 = t314 * t277;
t279 = qJDD(2) + qJDD(3);
t361 = pkin(5) + pkin(6);
t251 = t361 * t286;
t244 = qJD(1) * t251;
t356 = qJD(2) * pkin(2);
t241 = -t244 + t356;
t252 = t361 * t289;
t246 = qJD(1) * t252;
t330 = t360 * t246;
t218 = t285 * t241 + t330;
t336 = qJD(1) * qJD(2);
t325 = t289 * t336;
t335 = qJDD(1) * t286;
t224 = qJDD(2) * pkin(2) + t361 * (-t325 - t335);
t326 = t286 * t336;
t334 = qJDD(1) * t289;
t225 = t361 * (-t326 + t334);
t363 = t218 * qJD(3) - t360 * t224 + t285 * t225;
t194 = -t279 * pkin(3) + t363;
t278 = cos(t283);
t357 = g(3) * t278;
t366 = t194 + t357;
t345 = t285 * t289;
t240 = -qJD(1) * t345 - t286 * t328;
t215 = -pkin(3) * t240 - pkin(7) * t367;
t235 = qJD(4) - t367;
t364 = (pkin(7) * qJD(4) + t215) * t235;
t331 = qJD(2) * t361;
t245 = t286 * t331;
t247 = t289 * t331;
t305 = -t360 * t251 - t285 * t252;
t202 = t305 * qJD(3) - t360 * t245 - t285 * t247;
t243 = t360 * t286 + t345;
t223 = t280 * t243;
t323 = qJDD(1) * t360;
t312 = t285 * t335 - t289 * t323;
t207 = t223 * qJD(1) + t312;
t205 = qJDD(4) + t207;
t346 = t285 * t246;
t217 = t360 * t241 - t346;
t213 = -t280 * pkin(3) - t217;
t276 = -pkin(2) * t289 - pkin(1);
t304 = -t285 * t286 + t360 * t289;
t216 = -pkin(3) * t304 - pkin(7) * t243 + t276;
t222 = t280 * t304;
t230 = -t285 * t251 + t360 * t252;
t327 = qJD(3) * t360;
t339 = qJD(3) * t285;
t296 = t285 * t224 + t360 * t225 + t241 * t327 - t246 * t339;
t193 = t279 * pkin(7) + t296;
t250 = t276 * qJD(1);
t211 = -pkin(3) * t367 + pkin(7) * t240 + t250;
t318 = qJD(4) * t211 + t193;
t362 = t194 * t243 - t230 * t205 + t213 * t222 - (qJD(4) * t216 + t202) * t235 + t318 * t304;
t271 = g(3) * t277;
t206 = t285 * t334 + t286 * t323 + t369;
t284 = sin(qJ(4));
t288 = cos(qJ(4));
t337 = qJD(4) * t288;
t332 = t288 * t206 + t284 * t279 + t280 * t337;
t338 = qJD(4) * t284;
t195 = t240 * t338 + t332;
t355 = t195 * t284;
t354 = t213 * t367;
t353 = t213 * t243;
t352 = t216 * t205;
t226 = -t240 * t284 - t288 * t280;
t351 = t226 * t235;
t228 = -t240 * t288 + t280 * t284;
t350 = t228 * t235;
t349 = t284 * t205;
t348 = t284 * t287;
t347 = t284 * t290;
t344 = t287 * t288;
t343 = t288 * t205;
t342 = t288 * t290;
t281 = t286 ^ 2;
t341 = -t289 ^ 2 + t281;
t333 = t286 * t356;
t320 = t235 * t288;
t236 = pkin(2) * t326 + t276 * qJDD(1);
t191 = pkin(3) * t207 - pkin(7) * t206 + t236;
t214 = t280 * pkin(7) + t218;
t319 = qJD(4) * t214 - t191;
t274 = pkin(2) * t285 + pkin(7);
t316 = pkin(2) * t340 + qJD(4) * t274 + t215;
t220 = -t285 * t244 + t330;
t315 = pkin(2) * t339 - t220;
t313 = g(1) * t287 - g(2) * t290;
t311 = -t205 * t274 - t354;
t221 = -t360 * t244 - t346;
t310 = -pkin(2) * t327 + t221;
t198 = t211 * t284 + t214 * t288;
t309 = -t198 * t240 + t213 * t337 + t366 * t284;
t197 = t211 * t288 - t214 * t284;
t308 = t197 * t240 + t213 * t338 + t288 * t368;
t306 = -0.2e1 * pkin(1) * t336 - pkin(5) * qJDD(2);
t303 = t222 * t288 - t243 * t338;
t300 = -pkin(7) * t205 + t217 * t235 - t354;
t291 = qJD(2) ^ 2;
t298 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t291 + t313;
t292 = qJD(1) ^ 2;
t297 = pkin(1) * t292 - pkin(5) * qJDD(1) + t314;
t268 = t288 * t279;
t196 = t228 * qJD(4) + t206 * t284 - t268;
t295 = ((t195 - t351) * t288 + (-t196 - t350) * t284) * MDP(19) + (t228 * t320 + t355) * MDP(18) + (-t235 ^ 2 * t284 - t226 * t240 + t343) * MDP(21) + (t228 * t240 + t235 * t320 + t349) * MDP(20) + (t206 - t369) * MDP(13) + (-t312 + (-qJD(1) * t243 - t240) * t280) * MDP(14) + (t240 ^ 2 - t367 ^ 2) * MDP(12) + t279 * MDP(15) + (MDP(11) * t367 + t235 * MDP(22)) * t240;
t294 = -t250 * t367 + t314 * t278 + t271 - t296;
t293 = t250 * t240 - t357 - t363 + t368;
t275 = -t360 * pkin(2) - pkin(3);
t234 = t278 * t342 + t348;
t233 = -t278 * t347 + t344;
t232 = -t278 * t344 + t347;
t231 = t278 * t348 + t342;
t203 = t230 * qJD(3) - t285 * t245 + t360 * t247;
t201 = pkin(3) * t223 - pkin(7) * t222 + t333;
t190 = t288 * t191;
t1 = [qJDD(1) * MDP(1) + t313 * MDP(2) + t314 * MDP(3) + (qJDD(1) * t281 + 0.2e1 * t286 * t325) * MDP(4) + 0.2e1 * (t286 * t334 - t341 * t336) * MDP(5) + (qJDD(2) * t286 + t289 * t291) * MDP(6) + (qJDD(2) * t289 - t286 * t291) * MDP(7) + (t306 * t286 + t298 * t289) * MDP(9) + (-t298 * t286 + t306 * t289) * MDP(10) + (t206 * t243 - t222 * t240) * MDP(11) + (t206 * t304 - t207 * t243 + t222 * t367 + t223 * t240) * MDP(12) + (t222 * t280 + t243 * t279) * MDP(13) + (-t223 * t280 + t279 * t304) * MDP(14) + (-t203 * t280 + t207 * t276 + t223 * t250 - t236 * t304 + t313 * t278 + t279 * t305 - t333 * t367) * MDP(16) + (-t202 * t280 + t206 * t276 + t222 * t250 - t230 * t279 + t236 * t243 - t240 * t333 - t313 * t277) * MDP(17) + (t195 * t243 * t288 + t303 * t228) * MDP(18) + ((-t226 * t288 - t228 * t284) * t222 + (-t355 - t196 * t288 + (t226 * t284 - t228 * t288) * qJD(4)) * t243) * MDP(19) + (-t195 * t304 + t223 * t228 + t303 * t235 + t243 * t343) * MDP(20) + (-t243 * t349 + t196 * t304 - t223 * t226 + (-t222 * t284 - t243 * t337) * t235) * MDP(21) + (-t205 * t304 + t223 * t235) * MDP(22) + (-g(1) * t232 - g(2) * t234 - t190 * t304 - t305 * t196 + t197 * t223 + t203 * t226 + (t201 * t235 + t352 + (t214 * t304 - t230 * t235 + t353) * qJD(4)) * t288 + t362 * t284) * MDP(23) + (-g(1) * t231 - g(2) * t233 - t305 * t195 - t198 * t223 + t203 * t228 + (-(-qJD(4) * t230 + t201) * t235 - t352 - t319 * t304 - qJD(4) * t353) * t284 + t362 * t288) * MDP(24); MDP(7) * t334 + MDP(6) * t335 + (-g(3) * t289 + t297 * t286) * MDP(9) + (t221 * t280 + (t240 * t340 - t279 * t285 - t280 * t327) * pkin(2) + t294) * MDP(17) + (t220 * t280 + (t360 * t279 - t280 * t339 + t340 * t367) * pkin(2) + t293) * MDP(16) + qJDD(2) * MDP(8) + (g(3) * t286 + t297 * t289) * MDP(10) + (t275 * t196 - t366 * t288 + t311 * t284 + t315 * t226 + (t310 * t284 - t316 * t288) * t235 + t308) * MDP(23) + (t275 * t195 + t311 * t288 - t284 * t368 + t315 * t228 + (t316 * t284 + t310 * t288) * t235 + t309) * MDP(24) + t295 + (-t286 * t289 * MDP(4) + t341 * MDP(5)) * t292; (t217 * t280 + t294) * MDP(17) + (t218 * t280 + t293) * MDP(16) + (-pkin(3) * t196 - t218 * t226 + t300 * t284 + (-t366 - t364) * t288 + t308) * MDP(23) + (-pkin(3) * t195 - t218 * t228 + t300 * t288 + (-t368 + t364) * t284 + t309) * MDP(24) + t295; t228 * t226 * MDP(18) + (-t226 ^ 2 + t228 ^ 2) * MDP(19) + (t332 + t351) * MDP(20) + (t268 + t350) * MDP(21) + t205 * MDP(22) + (-g(1) * t233 + g(2) * t231 + t198 * t235 - t213 * t228 + t190) * MDP(23) + (g(1) * t234 - g(2) * t232 + t197 * t235 + t213 * t226) * MDP(24) + ((-t193 + t271) * MDP(24) + (MDP(21) * t240 - MDP(23) * t214 - MDP(24) * t211) * qJD(4)) * t288 + (qJD(4) * t240 * MDP(20) + (-qJD(4) * t280 - t206) * MDP(21) + (-t318 + t271) * MDP(23) + t319 * MDP(24)) * t284;];
tau = t1;
