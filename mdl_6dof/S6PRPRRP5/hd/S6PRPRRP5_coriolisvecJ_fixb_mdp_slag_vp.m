% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:47
% EndTime: 2019-03-08 20:16:52
% DurationCPUTime: 1.96s
% Computational Cost: add. (1687->305), mult. (3756->433), div. (0->0), fcn. (2466->8), ass. (0->147)
t264 = cos(qJ(4));
t257 = t264 ^ 2;
t261 = sin(qJ(4));
t359 = (t261 ^ 2 - t257) * MDP(9);
t281 = pkin(4) * t264 + pkin(9) * t261;
t235 = t281 * qJD(4) + qJD(3);
t263 = cos(qJ(5));
t265 = cos(qJ(2));
t258 = sin(pkin(6));
t319 = qJD(1) * t258;
t260 = sin(qJ(5));
t262 = sin(qJ(2));
t335 = t260 * t262;
t358 = -(-t261 * t335 + t263 * t265) * t319 + t263 * t235;
t244 = pkin(4) * t261 - pkin(9) * t264 + qJ(3);
t266 = -pkin(2) - pkin(8);
t305 = t263 * qJD(4);
t309 = qJD(5) * t263;
t331 = t262 * t263;
t357 = (t260 * t265 + t261 * t331) * t319 - t264 * t266 * t305 - t260 * t235 - t244 * t309;
t332 = t261 * t266;
t322 = t260 * t244 + t263 * t332;
t298 = t265 * t319;
t279 = qJD(3) - t298;
t228 = t266 * qJD(2) + t279;
t259 = cos(pkin(6));
t318 = qJD(1) * t261;
t211 = t228 * t264 - t259 * t318;
t314 = qJD(2) * t264;
t294 = t263 * t314;
t313 = qJD(4) * t260;
t240 = t294 + t313;
t356 = t240 ^ 2;
t355 = -qJ(6) - pkin(9);
t354 = qJD(2) * pkin(2);
t316 = qJD(2) * t258;
t297 = t262 * t316;
t282 = t264 * t297;
t337 = t259 * t264;
t252 = qJD(1) * t337;
t312 = qJD(4) * t261;
t325 = -qJD(4) * t252 - t228 * t312;
t198 = -qJD(1) * t282 - t325;
t353 = t198 * t260;
t352 = t198 * t263;
t203 = -qJD(4) * pkin(4) - t211;
t351 = t203 * t260;
t350 = t203 * t263;
t302 = qJD(4) * qJD(5);
t254 = t263 * t302;
t308 = qJD(5) * t264;
t292 = t260 * t308;
t293 = t261 * t305;
t273 = t292 + t293;
t215 = t273 * qJD(2) - t254;
t349 = t215 * t260;
t315 = qJD(2) * t261;
t290 = t260 * t315;
t249 = qJD(4) * t290;
t216 = qJD(5) * t240 - t249;
t348 = t216 * t263;
t295 = t260 * t314;
t238 = t295 - t305;
t346 = t238 * t260;
t345 = t238 * t263;
t253 = qJD(5) + t315;
t344 = t240 * t253;
t343 = t240 * t260;
t342 = t240 * t263;
t341 = t253 * t260;
t340 = t253 * t261;
t339 = t253 * t263;
t338 = t258 * t265;
t299 = t262 * t319;
t219 = t244 * qJD(2) + t299;
t336 = t260 * t219;
t334 = t260 * t266;
t333 = t261 * t263;
t330 = t263 * t264;
t212 = t261 * t228 + t252;
t204 = qJD(4) * pkin(9) + t212;
t192 = -t204 * t260 + t263 * t219;
t186 = -qJ(6) * t240 + t192;
t183 = pkin(5) * t253 + t186;
t329 = t183 - t186;
t288 = pkin(5) - t334;
t307 = qJD(6) * t263;
t310 = qJD(5) * t260;
t328 = qJ(6) * t293 - t322 * qJD(5) + (qJ(6) * t310 + t288 * qJD(4) - t307) * t264 + t358;
t291 = t263 * t308;
t327 = -qJ(6) * t291 + (-qJD(6) * t264 + (qJ(6) * qJD(4) - qJD(5) * t266) * t261) * t260 - t357;
t242 = t281 * qJD(2);
t326 = t263 * t211 + t260 * t242;
t286 = qJD(5) * t355;
t324 = -qJ(6) * t290 + t260 * t286 + t307 - t326;
t283 = -t211 * t260 + t263 * t242;
t323 = -qJD(6) * t260 + t263 * t286 - (pkin(5) * t264 + qJ(6) * t333) * qJD(2) - t283;
t267 = qJD(4) ^ 2;
t268 = qJD(2) ^ 2;
t320 = -t267 - t268;
t317 = qJD(2) * qJ(3);
t311 = qJD(4) * t264;
t306 = t203 * qJD(5);
t304 = qJD(2) * qJD(4);
t303 = qJD(4) * MDP(19);
t197 = t228 * t311 + (-qJD(4) * t259 + t297) * t318;
t213 = (t235 + t298) * qJD(2);
t301 = -t263 * t197 - t260 * t213 - t219 * t309;
t296 = t265 * t316;
t289 = t264 * t304;
t287 = pkin(5) * t260 - t266;
t285 = -t260 * t197 + t263 * t213;
t284 = t253 * t266 + t204;
t243 = t299 + t317;
t280 = -t243 + t299;
t193 = t204 * t263 + t336;
t187 = -qJ(6) * t238 + t193;
t278 = t183 * t263 + t187 * t260;
t277 = t183 * t260 - t187 * t263;
t276 = qJD(2) * t257 - t340;
t275 = -pkin(9) * t311 + t203 * t261;
t227 = -t261 * t338 + t337;
t208 = -t227 * t260 + t258 * t331;
t209 = t227 * t263 + t258 * t335;
t226 = t259 * t261 + t264 * t338;
t274 = t204 * t310 + t301;
t272 = t280 - t317;
t190 = pkin(5) * t216 + t198;
t236 = (qJD(3) + t298) * qJD(2);
t271 = t279 * qJD(2) - t266 * t267 + t236;
t270 = -t193 * qJD(5) + t285;
t269 = (t342 + t346) * MDP(22) - t278 * MDP(23) + (-t263 * MDP(20) + t260 * MDP(21)) * t253;
t247 = t355 * t263;
t246 = t355 * t260;
t237 = t279 - t354;
t234 = t238 ^ 2;
t233 = t263 * t244;
t207 = t227 * qJD(4) - t282;
t206 = -t226 * qJD(4) + t261 * t297;
t205 = -qJ(6) * t260 * t264 + t322;
t199 = -qJ(6) * t330 + t288 * t261 + t233;
t195 = pkin(5) * t238 + qJD(6) + t203;
t189 = t208 * qJD(5) + t206 * t263 + t260 * t296;
t188 = -t209 * qJD(5) - t206 * t260 + t263 * t296;
t182 = -qJ(6) * t216 - qJD(6) * t238 - t274;
t181 = pkin(5) * t289 + qJ(6) * t215 - qJD(6) * t240 + t270;
t1 = [(t188 * t253 + t207 * t238 + t216 * t226) * MDP(20) + (-t189 * t253 + t207 * t240 - t215 * t226) * MDP(21) + (-t188 * t240 - t189 * t238 + t208 * t215 - t209 * t216) * MDP(22) + (t181 * t208 + t182 * t209 + t183 * t188 + t187 * t189 + t190 * t226 + t195 * t207) * MDP(23) + (-t207 * MDP(13) - t206 * MDP(14) + (MDP(20) * t208 - MDP(21) * t209) * t314) * qJD(4) + ((t243 * qJD(2) * MDP(7) + (t261 * MDP(13) + MDP(14) * t264 - MDP(4) + MDP(6)) * t268) * t265 + (t236 * MDP(7) + (-MDP(3) + MDP(5)) * t268 + ((MDP(13) * t264 - MDP(14) * t261) * qJD(4) + (t237 - t298) * MDP(7)) * qJD(2)) * t262) * t258; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t236 + qJD(3) * t243 + (-t243 * t265 + (-t237 - t354) * t262) * t319) * MDP(7) + 0.2e1 * t304 * t359 - t267 * t264 * MDP(11) - t272 * t311 * MDP(13) + (t271 * t264 + t272 * t312) * MDP(14) + (-t215 * t330 - t273 * t240) * MDP(15) + ((t343 + t345) * t312 + (t349 - t348 + (-t342 + t346) * qJD(5)) * t264) * MDP(16) + (-t253 * t292 + (t240 * t264 + t276 * t263) * qJD(4)) * MDP(17) + (-t253 * t291 + (-t238 * t264 - t276 * t260) * qJD(4)) * MDP(18) + (t253 + t315) * t264 * t303 + ((-t244 * t310 + t358) * t253 + (t238 * t299 + t263 * t306 + t353 - t266 * t216 + (-t253 * t334 + (-t260 * t332 + t233) * qJD(2) + t192) * qJD(4)) * t264) * MDP(20) + (t357 * t253 + (t240 * t299 - t260 * t306 + t352 + t266 * t215 + (-t322 * qJD(2) - t193) * qJD(4)) * t264) * MDP(21) + (t199 * t215 - t205 * t216 - t328 * t240 - t327 * t238 + t278 * t312 + (t277 * qJD(5) - t181 * t263 - t182 * t260) * t264) * MDP(22) + (t181 * t199 + t182 * t205 + t327 * t187 + t328 * t183 - t195 * t287 * t312 + (t190 * t287 + (pkin(5) * t309 + t299) * t195) * t264) * MDP(23) + (-0.2e1 * MDP(8) * t289 - t267 * MDP(10) + t271 * MDP(13) - t215 * MDP(17) - t216 * MDP(18) + ((t238 * t266 - t351) * qJD(4) + (-t284 * t263 - t336) * qJD(5) + t285) * MDP(20) + (t284 * t310 + (t240 * t266 - t350) * qJD(4) + t301) * MDP(21)) * t261; -t268 * MDP(6) + (t280 * MDP(7) + t269) * qJD(2) + (t320 * MDP(14) - t216 * MDP(20) + t215 * MDP(21) - t190 * MDP(23) + ((t343 - t345) * MDP(22) - t277 * MDP(23) + (-t260 * MDP(20) - t263 * MDP(21)) * t253) * qJD(4)) * t264 + (t320 * MDP(13) + (-t348 - t349) * MDP(22) + (-t181 * t260 + t182 * t263) * MDP(23) + t269 * qJD(5) + ((t238 - t295) * MDP(20) + (t240 - t294) * MDP(21) + t195 * MDP(23)) * qJD(4)) * t261; (qJD(4) * t212 + t280 * t314 + t325) * MDP(13) - t280 * t315 * MDP(14) + (t240 * t339 - t349) * MDP(15) + ((-t238 * t253 - t215) * t263 + (-t216 - t344) * t260) * MDP(16) + (t253 * t309 + (t253 * t333 + (-t240 + t313) * t264) * qJD(2)) * MDP(17) + (-t253 * t310 + (-t260 * t340 + (t238 + t305) * t264) * qJD(2)) * MDP(18) - t253 * MDP(19) * t314 + (-pkin(4) * t216 - t352 - t283 * t253 - t212 * t238 + (-pkin(9) * t339 + t351) * qJD(5) + (-t192 * t264 + t275 * t260) * qJD(2)) * MDP(20) + (pkin(4) * t215 + t353 + t326 * t253 - t212 * t240 + (pkin(9) * t341 + t350) * qJD(5) + (t193 * t264 + t275 * t263) * qJD(2)) * MDP(21) + (t215 * t246 + t216 * t247 - t323 * t240 - t324 * t238 + (-t253 * t183 + t182) * t263 + (-t253 * t187 - t181) * t260) * MDP(22) + (-t182 * t247 + t181 * t246 + t190 * (-pkin(5) * t263 - pkin(4)) + (pkin(5) * t341 - t212) * t195 + t324 * t187 + t323 * t183) * MDP(23) + (t264 * t261 * MDP(8) - t359) * t268; (-t234 + t356) * MDP(16) + t254 * MDP(17) + (-t260 * t302 + t249 + t344) * MDP(18) + (t193 * t253 - t203 * t240 + t270) * MDP(20) + (t192 * t253 + t274) * MDP(21) + t329 * MDP(23) * t187 + (t215 * MDP(22) + (-t195 * t240 + t181) * MDP(23)) * pkin(5) + (t240 * MDP(15) + t253 * MDP(17) + t203 * MDP(21) - t329 * MDP(22)) * t238 + (-MDP(17) * t293 + (t303 + (-t260 * MDP(17) - t263 * MDP(18)) * qJD(5)) * t264) * qJD(2); (-t234 - t356) * MDP(22) + (t183 * t240 + t187 * t238 + t190) * MDP(23);];
tauc  = t1;
