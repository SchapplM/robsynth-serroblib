% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:47
% EndTime: 2019-03-08 19:37:51
% DurationCPUTime: 2.19s
% Computational Cost: add. (1242->286), mult. (3151->399), div. (0->0), fcn. (2295->10), ass. (0->143)
t279 = cos(qJ(2));
t271 = sin(pkin(6));
t337 = qJD(1) * t271;
t313 = t279 * t337;
t250 = qJD(2) * pkin(2) + t313;
t272 = cos(pkin(11));
t276 = sin(qJ(2));
t314 = t276 * t337;
t252 = t272 * t314;
t270 = sin(pkin(11));
t219 = t270 * t250 + t252;
t214 = qJD(2) * pkin(8) + t219;
t273 = cos(pkin(6));
t260 = qJD(1) * t273 + qJD(3);
t275 = sin(qJ(4));
t278 = cos(qJ(4));
t340 = -t275 * t214 + t278 * t260;
t360 = -qJD(5) + t340;
t199 = -qJD(4) * pkin(4) - t360;
t335 = qJD(2) * t278;
t268 = t275 ^ 2;
t269 = t278 ^ 2;
t359 = (t268 - t269) * MDP(7);
t230 = (t270 * t276 - t272 * t279) * t271;
t226 = t270 * t313 + t252;
t329 = qJD(5) * t275;
t332 = qJD(4) * t275;
t358 = -pkin(4) * t332 + t226 + t329;
t205 = t278 * t214 + t275 * t260;
t200 = -qJD(4) * qJ(5) - t205;
t336 = qJD(2) * t275;
t319 = pkin(5) * t336;
t323 = t319 - t360;
t322 = qJD(2) * qJD(4);
t308 = t275 * t322;
t261 = pkin(4) * t308;
t231 = (t270 * t279 + t272 * t276) * t271;
t284 = qJD(1) * t231 - t329;
t330 = qJD(4) * t278;
t315 = qJ(5) * t330;
t203 = t261 + (t284 - t315) * qJD(2);
t342 = t315 + t358;
t263 = pkin(2) * t270 + pkin(8);
t281 = qJD(4) ^ 2;
t345 = t263 * t281;
t357 = qJD(2) * t342 - t203 - t345;
t280 = -pkin(4) - pkin(9);
t356 = pkin(2) * t272;
t355 = pkin(5) + t263;
t228 = qJD(2) * t230;
t221 = qJD(1) * t228;
t316 = t214 * t332 + t278 * t221 - t260 * t330;
t187 = (qJD(5) - t319) * qJD(4) - t316;
t274 = sin(qJ(6));
t354 = t187 * t274;
t277 = cos(qJ(6));
t353 = t187 * t277;
t324 = t277 * qJD(2);
t333 = qJD(4) * t274;
t244 = t278 * t324 + t333;
t255 = t274 * t308;
t222 = -qJD(6) * t244 + t255;
t352 = t222 * t277;
t327 = qJD(6) * t277;
t309 = t274 * t335;
t339 = qJD(6) * t309 + t277 * t308;
t223 = qJD(4) * t327 - t339;
t351 = t223 * t275;
t262 = qJD(6) + t336;
t350 = t244 * t262;
t349 = t244 * t278;
t331 = qJD(4) * t277;
t246 = -t309 + t331;
t348 = t246 * t262;
t347 = t262 * t275;
t346 = t262 * t280;
t344 = t222 * t275 + t246 * t330;
t299 = pkin(9) * t275 - qJ(5) * t278;
t286 = t299 * qJD(4);
t343 = -t286 + t358;
t328 = qJD(6) * t274;
t311 = t262 * t328;
t318 = t277 * t347;
t341 = qJD(4) * t318 + t278 * t311;
t326 = t274 * MDP(23);
t325 = t277 * MDP(23);
t321 = MDP(11) - MDP(14);
t320 = MDP(12) - MDP(15);
t282 = qJD(2) ^ 2;
t317 = t275 * t278 * t282;
t191 = t214 * t330 - t275 * t221 + t260 * t332;
t312 = t269 * t324;
t310 = t278 * t327;
t241 = t355 * t278;
t307 = t278 * t322;
t306 = MDP(21) * t335;
t305 = -qJ(5) * t275 - pkin(3);
t189 = pkin(5) * t307 + t191;
t196 = t261 + (t286 + t284) * qJD(2);
t303 = t277 * t189 - t196 * t274;
t251 = t270 * t314;
t218 = t250 * t272 - t251;
t300 = t262 * t310;
t194 = qJD(4) * t280 + t323;
t285 = t278 * t280 + t305;
t201 = qJD(2) * t285 - t218;
t185 = t194 * t277 - t201 * t274;
t186 = t194 * t274 + t201 * t277;
t298 = t199 * t275 - t200 * t278;
t212 = t231 * t278 + t273 * t275;
t232 = t285 - t356;
t240 = t355 * t275;
t297 = t232 * t277 + t240 * t274;
t295 = qJD(2) * t269 - t347;
t294 = t262 * t274;
t292 = -pkin(4) * t278 + t305;
t291 = qJD(4) * t205 - t191;
t265 = pkin(5) * t335;
t195 = -t200 + t265;
t290 = t195 * t275 + t280 * t330;
t227 = qJD(2) * t231;
t220 = qJD(1) * t227;
t289 = qJD(2) * t226 - t220 - t345;
t206 = qJD(2) * t292 - t218;
t229 = t272 * t313 - t251;
t239 = t292 - t356;
t288 = qJD(4) * (-qJD(2) * t239 - t206 - t229);
t213 = -qJD(2) * pkin(3) - t218;
t287 = qJD(4) * (qJD(2) * (-pkin(3) - t356) + t213 + t229);
t190 = -qJD(4) * qJD(5) + t316;
t283 = -t190 * t278 + t191 * t275 + (t199 * t278 + t200 * t275) * qJD(4);
t267 = pkin(4) * t336;
t257 = t277 * t307;
t247 = -qJ(5) * t335 + t267;
t236 = qJD(4) * t241;
t235 = t355 * t332;
t233 = qJD(2) * t299 + t267;
t211 = t231 * t275 - t273 * t278;
t202 = t206 * t336;
t198 = t265 + t205;
t193 = -t231 * t332 + (qJD(4) * t273 - t228) * t278;
t192 = qJD(4) * t212 - t228 * t275;
t1 = [(-t218 * t227 - t219 * t228 + t220 * t230 - t221 * t231) * MDP(5) + (-t190 * t212 + t191 * t211 + t192 * t199 - t193 * t200 + t203 * t230 + t206 * t227) * MDP(16) + ((t192 * t277 - t227 * t274 + (-t211 * t274 - t230 * t277) * qJD(6)) * t262 + t193 * t244 + t212 * t223) * MDP(22) + (-(t192 * t274 + t227 * t277 + (t211 * t277 - t230 * t274) * qJD(6)) * t262 + t193 * t246 + t212 * t222) * MDP(23) + (-MDP(3) * t276 - MDP(4) * t279) * t282 * t271 + (-t321 * t192 - t320 * t193) * qJD(4) + ((t192 * t275 + t193 * t278) * MDP(13) + (t275 * t320 - t278 * t321) * t227 + ((-t212 * MDP(13) + t230 * t321) * t275 + ((t277 * MDP(22) + MDP(13) - t326) * t211 + (-t274 * MDP(22) + t320 - t325) * t230) * t278) * qJD(4)) * qJD(2); (t218 * t226 - t219 * t229 + (-t220 * t272 - t221 * t270) * pkin(2)) * MDP(5) + 0.2e1 * t275 * MDP(6) * t307 - 0.2e1 * t322 * t359 + (t275 * t287 + t278 * t289) * MDP(11) + (-t275 * t289 + t278 * t287) * MDP(12) + ((-t268 - t269) * t229 * qJD(2) + t283) * MDP(13) + (t275 * t288 - t357 * t278) * MDP(14) + (t357 * t275 + t278 * t288) * MDP(15) + (t203 * t239 - t206 * t342 - t229 * t298 + t263 * t283) * MDP(16) + (-t222 * t274 * t278 + (t274 * t332 - t310) * t246) * MDP(17) + ((-t244 * t274 + t246 * t277) * t332 + (-t352 + t223 * t274 + (t244 * t277 + t246 * t274) * qJD(6)) * t278) * MDP(18) + (-t295 * t333 - t300 + t344) * MDP(19) + (-t351 + (-t312 - t349) * qJD(4) + t341) * MDP(20) + (t262 + t336) * MDP(21) * t330 + (t241 * t223 - t235 * t244 + (-t195 * t331 + t303) * t275 + ((-t229 * t275 + t236) * t277 + t343 * t274) * t262 + (-t186 * t275 - t262 * t297) * qJD(6) + (-t195 * t328 + t353 - t229 * t244 + ((-t232 * t274 + t240 * t277) * qJD(2) + t185) * qJD(4)) * t278) * MDP(22) + (t241 * t222 - t235 * t246 + (-(qJD(6) * t194 + t196) * t275 + (-qJD(6) * t240 + t343) * t262) * t277 + (-(-qJD(6) * t232 + t236) * t262 + (qJD(4) * t195 + qJD(6) * t201 + t229 * t262 - t189) * t275) * t274 + (-t195 * t327 - t354 - t229 * t246 + (-qJD(2) * t297 - t186) * qJD(4)) * t278) * MDP(23) + (MDP(8) * t278 - MDP(9) * t275) * t281; (-t190 * t275 - t191 * t278) * MDP(16) + (t341 + t351) * MDP(22) + (t300 + t344) * MDP(23) + (t298 * MDP(16) + (-t312 + t349) * MDP(22) + t295 * t326) * qJD(4) + (-t321 * t275 - t278 * t320) * t281; -MDP(6) * t317 + t282 * t359 + (-t213 * t336 + t291) * MDP(11) + (qJD(4) * t340 - t213 * t335 + t316) * MDP(12) + (-t247 * t335 + t202 - t291) * MDP(14) + ((0.2e1 * qJD(5) - t340) * qJD(4) + (t206 * t278 + t247 * t275) * qJD(2) - t316) * MDP(15) + (-pkin(4) * t191 - qJ(5) * t190 - t199 * t205 + t200 * t360 - t206 * t247) * MDP(16) + (-t246 * t294 + t352) * MDP(17) + ((-t223 - t348) * t277 + (-t222 + t350) * t274) * MDP(18) + (-t311 + t257 + (-t246 * t278 - t274 * t347) * qJD(2)) * MDP(19) + (-t262 * t327 + (-t318 + (t244 - t333) * t278) * qJD(2)) * MDP(20) - t262 * t306 + (qJ(5) * t223 + t354 - (t198 * t277 - t233 * t274) * t262 + t323 * t244 + (t195 * t277 - t274 * t346) * qJD(6) + (-t185 * t278 + t277 * t290) * qJD(2)) * MDP(22) + (qJ(5) * t222 + t353 + (t198 * t274 + t233 * t277) * t262 + t323 * t246 + (-t195 * t274 - t277 * t346) * qJD(6) + (t186 * t278 - t274 * t290) * qJD(2)) * MDP(23); MDP(14) * t317 + (-t268 * t282 - t281) * MDP(15) + (t202 + t191) * MDP(16) + t257 * MDP(22) + (t200 * MDP(16) - t244 * MDP(22) + (-t246 - t309) * MDP(23)) * qJD(4) + (-MDP(22) * t294 - t262 * t325) * t262; t246 * t244 * MDP(17) + (-t244 ^ 2 + t246 ^ 2) * MDP(18) + (t255 + t350) * MDP(19) + (t339 + t348) * MDP(20) + qJD(4) * t306 + (t186 * t262 - t195 * t246 + t303) * MDP(22) + (t185 * t262 - t189 * t274 + t195 * t244 - t196 * t277) * MDP(23) + (-MDP(19) * t244 - MDP(20) * t331 - MDP(22) * t186 - MDP(23) * t185) * qJD(6);];
tauc  = t1;
