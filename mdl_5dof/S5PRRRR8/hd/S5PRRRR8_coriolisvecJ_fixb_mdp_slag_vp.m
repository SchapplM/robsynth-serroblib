% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:49
% EndTime: 2019-12-05 17:16:59
% DurationCPUTime: 3.00s
% Computational Cost: add. (1789->260), mult. (4538->390), div. (0->0), fcn. (3430->10), ass. (0->128)
t275 = sin(qJ(4));
t279 = cos(qJ(3));
t356 = cos(qJ(4));
t313 = qJD(2) * t356;
t276 = sin(qJ(3));
t331 = qJD(2) * t276;
t362 = -t275 * t331 + t279 * t313;
t357 = pkin(7) + pkin(8);
t269 = qJD(3) + qJD(4);
t361 = MDP(5) * t279;
t360 = (t276 ^ 2 - t279 ^ 2) * MDP(6);
t277 = sin(qJ(2));
t272 = sin(pkin(5));
t334 = qJD(1) * t272;
t318 = t277 * t334;
t309 = t357 * qJD(2) + t318;
t273 = cos(pkin(5));
t333 = qJD(1) * t273;
t230 = -t309 * t276 + t279 * t333;
t354 = qJD(3) * pkin(3);
t287 = t276 * t354 - t318;
t359 = t276 * MDP(10) + t279 * MDP(11);
t280 = cos(qJ(2));
t332 = qJD(2) * t272;
t311 = qJD(1) * t332;
t304 = t280 * t311;
t214 = qJD(3) * t230 + t279 * t304;
t231 = t276 * t333 + t309 * t279;
t215 = -qJD(3) * t231 - t276 * t304;
t225 = t230 + t354;
t312 = qJD(4) * t356;
t329 = qJD(4) * t275;
t185 = t356 * t214 + t275 * t215 + t225 * t312 - t231 * t329;
t319 = t356 * t231;
t204 = t275 * t225 + t319;
t308 = t275 * t214 - t356 * t215;
t186 = t204 * qJD(4) + t308;
t343 = t275 * t231;
t203 = t356 * t225 - t343;
t198 = -t269 * pkin(4) - t203;
t268 = -pkin(3) * t279 - pkin(2);
t317 = t280 * t334;
t239 = t268 * qJD(2) - t317;
t342 = t275 * t279;
t246 = -qJD(2) * t342 - t276 * t313;
t211 = -pkin(4) * t362 + pkin(9) * t246 + t239;
t248 = t356 * t276 + t342;
t227 = t269 * t248;
t221 = t227 * qJD(2);
t290 = -t275 * t276 + t356 * t279;
t223 = -pkin(4) * t290 - pkin(9) * t248 + t268;
t226 = t269 * t290;
t253 = t357 * t276;
t254 = t357 * t279;
t238 = -t275 * t253 + t356 * t254;
t240 = qJD(5) - t362;
t320 = qJD(3) * t357;
t249 = t276 * t320;
t250 = t279 * t320;
t291 = -t356 * t253 - t275 * t254;
t338 = -t291 * qJD(4) + t356 * t249 + t275 * t250 + t290 * t317;
t358 = (qJD(5) * t211 + t185) * t290 + t186 * t248 + t198 * t226 + (-qJD(5) * t223 + t338) * t240 - t238 * t221;
t355 = qJD(2) * pkin(2);
t353 = t198 * t362;
t352 = t198 * t248;
t274 = sin(qJ(5));
t328 = qJD(5) * t274;
t220 = t362 * t269;
t278 = cos(qJ(5));
t327 = qJD(5) * t278;
t336 = t278 * t220 + t269 * t327;
t201 = t246 * t328 + t336;
t351 = t201 * t274;
t350 = t220 * t274;
t349 = t223 * t221;
t346 = t246 * t274;
t234 = -t278 * t269 - t346;
t348 = t234 * t240;
t296 = t246 * t278 - t269 * t274;
t347 = t296 * t240;
t345 = t272 * t277;
t344 = t274 * t221;
t281 = qJD(3) ^ 2;
t341 = t276 * t281;
t340 = t278 * t221;
t339 = t279 * t281;
t337 = t238 * qJD(4) - t248 * t317 - t275 * t249 + t356 * t250;
t323 = qJD(2) * qJD(3);
t310 = t276 * t323;
t241 = pkin(3) * t310 + t277 * t311;
t330 = qJD(2) * t277;
t326 = qJD(5) * t280;
t322 = pkin(3) * t331;
t315 = t280 * t332;
t307 = t278 * t240;
t222 = -pkin(4) * t246 - pkin(9) * t362;
t266 = pkin(3) * t275 + pkin(9);
t305 = qJD(5) * t266 + t222 + t322;
t205 = t275 * t230 + t319;
t303 = pkin(3) * t329 - t205;
t199 = t269 * pkin(9) + t204;
t191 = t199 * t278 + t211 * t274;
t302 = t186 * t274 - t191 * t246 + t198 * t327;
t301 = pkin(4) * t227 - pkin(9) * t226 + t287;
t298 = -t221 * t266 - t353;
t297 = t199 * t274 - t211 * t278;
t206 = t356 * t230 - t343;
t295 = -pkin(3) * t312 + t206;
t243 = t273 * t276 + t279 * t345;
t242 = t273 * t279 - t276 * t345;
t294 = -t186 * t278 + t198 * t328 - t246 * t297;
t293 = t239 * t246 - t308;
t292 = t356 * t242 - t275 * t243;
t218 = t275 * t242 + t356 * t243;
t289 = t226 * t278 - t248 * t328;
t286 = -0.2e1 * qJD(3) * t355;
t202 = -t296 * qJD(5) + t350;
t285 = ((t201 - t348) * t278 + (-t202 + t347) * t274) * MDP(20) + (-t296 * t307 + t351) * MDP(19) + (-t240 ^ 2 * t274 - t234 * t246 + t340) * MDP(22) + (t240 * t307 - t246 * t296 + t344) * MDP(21) + t220 * MDP(14) + (t246 ^ 2 - t362 ^ 2) * MDP(13) + (MDP(12) * t362 + t240 * MDP(23)) * t246 + (-t362 * MDP(14) + (-qJD(2) * t248 - t246) * MDP(15)) * t269;
t283 = -t239 * t362 - t185;
t282 = qJD(2) ^ 2;
t267 = -t356 * pkin(3) - pkin(4);
t229 = -t243 * qJD(3) - t276 * t315;
t228 = t242 * qJD(3) + t279 * t315;
t195 = pkin(4) * t221 - pkin(9) * t220 + t241;
t194 = t278 * t195;
t193 = t218 * qJD(4) + t275 * t228 - t356 * t229;
t192 = t292 * qJD(4) + t356 * t228 + t275 * t229;
t1 = [((-t192 * t274 - t218 * t327) * t240 - t218 * t344 + t193 * t234 - t292 * t202) * MDP(24) + (-(t192 * t278 - t218 * t328) * t240 - t218 * t340 - t193 * t296 - t292 * t201) * MDP(25) + (-MDP(17) * t193 - MDP(18) * t192) * t269 + (MDP(10) * t229 - MDP(11) * t228) * qJD(3) + ((-t221 * t280 - t330 * t362) * MDP(17) + (-t220 * t280 - t246 * t330) * MDP(18) + ((t274 * t326 + t278 * t330) * t240 - t280 * t340) * MDP(24) + (-(t274 * t330 - t278 * t326) * t240 + t280 * t344) * MDP(25) - t359 * t280 * t323 + (-MDP(4) * t280 + (-MDP(10) * t279 + MDP(11) * t276 - MDP(3)) * t277) * t282) * t272; 0.2e1 * t310 * t361 - 0.2e1 * t323 * t360 + MDP(7) * t339 - MDP(8) * t341 + (-pkin(7) * t339 + t276 * t286) * MDP(10) + (pkin(7) * t341 + t279 * t286) * MDP(11) + (t220 * t248 - t226 * t246) * MDP(12) + (t220 * t290 - t221 * t248 + t226 * t362 + t227 * t246) * MDP(13) + (t221 * t268 + t227 * t239 - t241 * t290 - t287 * t362) * MDP(17) + (t220 * t268 + t226 * t239 + t241 * t248 - t287 * t246) * MDP(18) + (t201 * t248 * t278 - t289 * t296) * MDP(19) + ((-t234 * t278 + t274 * t296) * t226 + (-t351 - t202 * t278 + (t234 * t274 + t278 * t296) * qJD(5)) * t248) * MDP(20) + (-t201 * t290 - t227 * t296 + t289 * t240 + t248 * t340) * MDP(21) + (-t248 * t344 + t202 * t290 - t227 * t234 + (-t226 * t274 - t248 * t327) * t240) * MDP(22) + (-t221 * t290 + t227 * t240) * MDP(23) + (-t297 * t227 - t194 * t290 - t291 * t202 + t337 * t234 + (t349 + t301 * t240 + (t199 * t290 - t238 * t240 + t352) * qJD(5)) * t278 + t358 * t274) * MDP(24) + (-t191 * t227 - t291 * t201 - t337 * t296 + (-t349 + (-qJD(5) * t199 + t195) * t290 - qJD(5) * t352 + (qJD(5) * t238 - t301) * t240) * t274 + t358 * t278) * MDP(25) + (t226 * MDP(14) - t227 * MDP(15) - t337 * MDP(17) + t338 * MDP(18)) * t269; t285 + (t267 * t202 + t298 * t274 + t303 * t234 + (t295 * t274 - t305 * t278) * t240 + t294) * MDP(24) + (t267 * t201 + t298 * t278 - t303 * t296 + (t305 * t274 + t295 * t278) * t240 + t302) * MDP(25) + (t362 * t322 + t205 * t269 + (-t319 + (-pkin(3) * t269 - t225) * t275) * qJD(4) + t293) * MDP(17) + (t206 * t269 + (t246 * t331 - t269 * t312) * pkin(3) + t283) * MDP(18) + t359 * t355 * qJD(2) + (-t276 * t361 + t360) * t282; (t293 + (-qJD(4) + t269) * t204) * MDP(17) + (t203 * t269 + t283) * MDP(18) + (-pkin(4) * t202 - (-t203 * t274 + t222 * t278) * t240 - t204 * t234 - t274 * t353 + (-t240 * t327 - t344) * pkin(9) + t294) * MDP(24) + (-pkin(4) * t201 + (t203 * t278 + t222 * t274) * t240 + t204 * t296 - t278 * t353 + (t240 * t328 - t340) * pkin(9) + t302) * MDP(25) + t285; -t296 * t234 * MDP(19) + (-t234 ^ 2 + t296 ^ 2) * MDP(20) + (t336 + t348) * MDP(21) + (-t347 - t350) * MDP(22) + t221 * MDP(23) + (-t185 * t274 + t191 * t240 + t198 * t296 + t194) * MDP(24) + (-t185 * t278 - t195 * t274 + t198 * t234 - t240 * t297) * MDP(25) + (MDP(21) * t346 + t296 * MDP(22) - t191 * MDP(24) + t297 * MDP(25)) * qJD(5);];
tauc = t1;
