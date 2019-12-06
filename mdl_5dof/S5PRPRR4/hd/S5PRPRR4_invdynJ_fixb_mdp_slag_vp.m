% Calculate vector of inverse dynamics joint torques for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:42
% EndTime: 2019-12-05 15:51:47
% DurationCPUTime: 2.27s
% Computational Cost: add. (1078->274), mult. (2595->397), div. (0->0), fcn. (2207->12), ass. (0->137)
t274 = sin(qJ(2));
t277 = cos(qJ(2));
t268 = sin(pkin(5));
t332 = qJD(2) * t268;
t312 = qJD(1) * t332;
t318 = qJDD(1) * t268;
t357 = t274 * t318 + t277 * t312;
t276 = cos(qJ(4));
t320 = qJD(2) * qJD(4);
t311 = t276 * t320;
t273 = sin(qJ(4));
t317 = qJDD(2) * t273;
t356 = qJD(4) * qJD(5) + t311 + t317;
t249 = t277 * t318;
t223 = qJDD(2) * pkin(2) - t274 * t312 + t249;
t266 = sin(pkin(10));
t269 = cos(pkin(10));
t197 = t266 * t223 + t357 * t269;
t195 = qJDD(2) * pkin(7) + t197;
t333 = qJD(1) * t268;
t314 = t277 * t333;
t240 = qJD(2) * pkin(2) + t314;
t315 = t274 * t333;
t244 = t269 * t315;
t216 = t266 * t240 + t244;
t214 = qJD(2) * pkin(7) + t216;
t271 = cos(pkin(5));
t252 = qJD(1) * t271 + qJD(3);
t201 = t214 * t276 + t252 * t273;
t250 = t271 * qJDD(1) + qJDD(3);
t346 = t250 * t276;
t181 = -qJDD(4) * pkin(4) + qJD(4) * t201 + t195 * t273 - t346;
t330 = qJD(2) * t276;
t253 = -qJD(5) + t330;
t296 = t266 * t277 + t269 * t274;
t225 = t296 * t268;
t211 = t225 * t273 - t271 * t276;
t267 = sin(pkin(9));
t270 = cos(pkin(9));
t226 = t296 * t271;
t337 = t277 * t269;
t231 = t266 * t274 - t337;
t297 = -t226 * t267 - t231 * t270;
t298 = t226 * t270 - t231 * t267;
t343 = t268 * t276;
t291 = g(1) * (t267 * t343 - t273 * t297) + g(2) * (-t270 * t343 - t273 * t298) - g(3) * t211;
t302 = pkin(4) * t273 - pkin(8) * t276;
t355 = (pkin(8) * qJD(5) + t302 * qJD(2)) * t253 - t181 - t291;
t219 = t266 * t314 + t244;
t295 = -pkin(4) * t276 - pkin(8) * t273 - pkin(3);
t352 = pkin(2) * t269;
t229 = t295 - t352;
t316 = t276 * qJDD(2);
t230 = t273 * t320 + qJDD(5) - t316;
t239 = t302 * qJD(4);
t354 = (t219 - t239) * t253 + t229 * t230;
t199 = qJD(4) * pkin(8) + t201;
t256 = pkin(2) * t266 + pkin(7);
t344 = t268 * t274;
t224 = t266 * t344 - t268 * t337;
t292 = t231 * t271;
t300 = g(1) * (t267 * t292 - t270 * t296) + g(2) * (-t267 * t296 - t270 * t292);
t290 = -g(3) * t224 + t300;
t353 = (t253 * t256 + t199) * qJD(5) - t290;
t272 = sin(qJ(5));
t275 = cos(qJ(5));
t304 = t272 * qJDD(4) + t356 * t275;
t323 = qJD(5) * t273;
t310 = qJD(2) * t323;
t209 = -t272 * t310 + t304;
t351 = t209 * t272;
t321 = t275 * qJD(4);
t331 = qJD(2) * t273;
t233 = t272 * t331 - t321;
t349 = t233 * t253;
t327 = qJD(4) * t272;
t235 = t275 * t331 + t327;
t348 = t235 * t253;
t347 = t250 * t273;
t345 = t268 * t273;
t342 = t271 * t274;
t341 = t271 * t277;
t340 = t272 * t253;
t339 = t275 * t230;
t338 = t275 * t253;
t336 = qJDD(1) - g(3);
t264 = t273 ^ 2;
t334 = -t276 ^ 2 + t264;
t329 = qJD(4) * t233;
t328 = qJD(4) * t256;
t326 = qJD(4) * t273;
t325 = qJD(4) * t276;
t324 = qJD(5) * t272;
t322 = qJD(5) * t275;
t313 = t253 * t327;
t243 = t266 * t315;
t215 = t240 * t269 - t243;
t213 = -qJD(2) * pkin(3) - t215;
t306 = -qJD(2) * t213 - t195;
t202 = t295 * qJD(2) - t215;
t299 = t214 * t273 - t252 * t276;
t305 = -qJDD(4) * pkin(8) + t299 * qJD(4) - qJD(5) * t202 - t195 * t276 - t347;
t301 = g(1) * (t267 * t345 + t276 * t297) + g(2) * (-t270 * t345 + t276 * t298);
t196 = t223 * t269 - t357 * t266;
t212 = t225 * t276 + t271 * t273;
t189 = t212 * t275 + t224 * t272;
t188 = -t212 * t272 + t224 * t275;
t294 = -t272 * t230 + t253 * t322;
t293 = -t253 * t324 - t339;
t289 = -qJD(5) * t199 + t300;
t288 = -t272 * t323 + t276 * t321;
t286 = t301 + t305;
t222 = t269 * t314 - t243;
t257 = -pkin(3) - t352;
t285 = -qJDD(4) * t256 + (qJD(2) * t257 + t213 + t222) * qJD(4);
t284 = qJD(5) * t229 * t253 - g(1) * t297 - g(2) * t298 - g(3) * t225;
t198 = -qJD(4) * pkin(4) + t299;
t283 = -pkin(8) * t230 + (-t198 + t299) * t253;
t282 = qJD(4) * t198 - t222 * t253 - t230 * t256 - t305;
t278 = qJD(4) ^ 2;
t281 = -qJD(2) * t219 + t256 * t278 - t196 + t290 + (-pkin(3) + t257) * qJDD(2);
t280 = -g(1) * (-t267 * t341 - t270 * t274) - g(2) * (-t267 * t274 + t270 * t341) - g(3) * t268 * t277;
t279 = qJD(2) ^ 2;
t261 = t275 * qJDD(4);
t246 = qJDD(4) * t276 - t273 * t278;
t245 = qJDD(4) * t273 + t276 * t278;
t227 = t235 * t326;
t221 = t231 * t332;
t220 = qJD(2) * t225;
t210 = t275 * t310 - t261 + (t317 + (qJD(5) + t330) * qJD(4)) * t272;
t187 = -t211 * qJD(4) - t221 * t276;
t186 = t212 * qJD(4) - t221 * t273;
t185 = qJD(2) * t239 + t295 * qJDD(2) - t196;
t184 = t275 * t185;
t183 = t199 * t275 + t202 * t272;
t182 = -t199 * t272 + t202 * t275;
t1 = [t336 * MDP(1) + (-t196 * t224 + t197 * t225 - t215 * t220 - t216 * t221 + t250 * t271 - g(3)) * MDP(5) + (-qJD(4) * t186 - qJDD(4) * t211 - t224 * t316) * MDP(11) + (-qJD(4) * t187 - qJDD(4) * t212 + t224 * t317) * MDP(12) + (-(-t189 * qJD(5) - t187 * t272 + t220 * t275) * t253 + t188 * t230 + t186 * t233 + t211 * t210) * MDP(18) + ((t188 * qJD(5) + t187 * t275 + t220 * t272) * t253 - t189 * t230 + t186 * t235 + t211 * t209) * MDP(19) + ((-t220 * t276 + t224 * t326) * MDP(11) + (t220 * t273 + t224 * t325) * MDP(12)) * qJD(2) + ((qJDD(2) * t277 - t274 * t279) * MDP(3) + (-qJDD(2) * t274 - t277 * t279) * MDP(4)) * t268; qJDD(2) * MDP(2) + (t249 + t280) * MDP(3) + (-g(1) * (t267 * t342 - t270 * t277) - g(2) * (-t267 * t277 - t270 * t342) - t336 * t344) * MDP(4) + (t215 * t219 - t216 * t222 + (t196 * t269 + t197 * t266 + t280) * pkin(2)) * MDP(5) + (qJDD(2) * t264 + 0.2e1 * t273 * t311) * MDP(6) + 0.2e1 * (t273 * t316 - t334 * t320) * MDP(7) + t245 * MDP(8) + t246 * MDP(9) + (t285 * t273 - t281 * t276) * MDP(11) + (t281 * t273 + t285 * t276) * MDP(12) + (t209 * t273 * t275 + t288 * t235) * MDP(13) + ((-t233 * t275 - t235 * t272) * t325 + (-t351 - t210 * t275 + (t233 * t272 - t235 * t275) * qJD(5)) * t273) * MDP(14) + (-t209 * t276 - t288 * t253 + t273 * t339 + t227) * MDP(15) + ((t210 + t313) * t276 + (t294 - t329) * t273) * MDP(16) + (-t230 * t276 - t253 * t326) * MDP(17) + (t354 * t275 + t284 * t272 + (t233 * t328 + t282 * t272 + t275 * t353 - t184) * t276 + (t198 * t322 + t181 * t272 + t256 * t210 - t222 * t233 + (-t256 * t340 + t182) * qJD(4)) * t273) * MDP(18) + (-t354 * t272 + t284 * t275 + (t235 * t328 + t282 * t275 + (t185 - t353) * t272) * t276 + (-t198 * t324 + t181 * t275 + t256 * t209 - t222 * t235 + (-t256 * t338 - t183) * qJD(4)) * t273) * MDP(19); (-g(3) * t271 + (-g(1) * t267 + g(2) * t270) * t268 + t250) * MDP(5) + t246 * MDP(11) - t245 * MDP(12) + t227 * MDP(19) + ((-t210 + t313) * MDP(18) + (t253 * t321 - t209) * MDP(19)) * t276 + ((t294 + t329) * MDP(18) + t293 * MDP(19)) * t273; MDP(8) * t317 + MDP(9) * t316 + qJDD(4) * MDP(10) + (t306 * t273 - t291 + t346) * MDP(11) + (g(3) * t212 + t306 * t276 + t301 - t347) * MDP(12) + (-t235 * t338 + t351) * MDP(13) + ((t209 + t349) * t275 + (-t210 + t348) * t272) * MDP(14) + ((-t235 * t273 + t276 * t338) * qJD(2) - t294) * MDP(15) + ((t233 * t273 - t276 * t340) * qJD(2) - t293) * MDP(16) + t253 * MDP(17) * t331 + (-pkin(4) * t210 - t182 * t331 - t201 * t233 + t283 * t272 + t275 * t355) * MDP(18) + (-pkin(4) * t209 + t183 * t331 - t201 * t235 - t272 * t355 + t283 * t275) * MDP(19) + (-t273 * t276 * MDP(6) + t334 * MDP(7)) * t279; t235 * t233 * MDP(13) + (-t233 ^ 2 + t235 ^ 2) * MDP(14) + (t304 - t349) * MDP(15) + (t261 - t348) * MDP(16) + t230 * MDP(17) + (-g(3) * t188 - t183 * t253 - t198 * t235 + t184) * MDP(18) + (g(3) * t189 - t182 * t253 + t198 * t233) * MDP(19) + (-MDP(16) * t310 + t289 * MDP(18) + t286 * MDP(19)) * t275 + (-MDP(15) * t310 - t356 * MDP(16) + t286 * MDP(18) + (-t185 - t289) * MDP(19)) * t272;];
tau = t1;
