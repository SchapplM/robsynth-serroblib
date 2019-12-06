% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:34:08
% EndTime: 2019-12-05 18:34:11
% DurationCPUTime: 1.57s
% Computational Cost: add. (1744->209), mult. (3035->280), div. (0->0), fcn. (2238->8), ass. (0->115)
t284 = cos(qJ(2));
t340 = pkin(1) * qJD(1);
t318 = t284 * t340;
t301 = qJD(3) - t318;
t276 = qJD(1) + qJD(2);
t278 = cos(pkin(9));
t283 = cos(qJ(4));
t331 = t283 * t278;
t314 = t276 * t331;
t277 = sin(pkin(9));
t280 = sin(qJ(4));
t334 = t277 * t280;
t315 = t276 * t334;
t235 = -t314 + t315;
t282 = cos(qJ(5));
t227 = t282 * t235;
t259 = qJD(4) * t314;
t231 = -qJD(4) * t315 + t259;
t256 = t277 * t283 + t278 * t280;
t248 = t256 * qJD(4);
t232 = t276 * t248;
t237 = t256 * t276;
t279 = sin(qJ(5));
t321 = qJD(5) * t279;
t186 = -qJD(5) * t227 + t282 * t231 - t279 * t232 - t237 * t321;
t206 = t237 * t279 + t227;
t296 = t235 * t279 - t282 * t237;
t287 = t296 * qJD(5) - t231 * t279 - t282 * t232;
t275 = qJD(4) + qJD(5);
t336 = t206 * t275;
t337 = t296 * t275;
t356 = (-t206 ^ 2 + t296 ^ 2) * MDP(19) - t206 * MDP(18) * t296 + (t186 + t336) * MDP(20) + (t287 - t337) * MDP(21);
t339 = pkin(1) * qJD(2);
t316 = qJD(1) * t339;
t254 = t276 * qJD(3) + t284 * t316;
t324 = t277 ^ 2 + t278 ^ 2;
t355 = t324 * t254;
t291 = t256 * t254;
t281 = sin(qJ(2));
t319 = t281 * t340;
t258 = qJ(3) * t276 + t319;
t312 = pkin(7) * t276 + t258;
t224 = t312 * t277;
t225 = t312 * t278;
t297 = t224 * t280 - t225 * t283;
t183 = -pkin(8) * t231 + t297 * qJD(4) - t291;
t195 = -pkin(8) * t235 - t297;
t270 = -t278 * pkin(3) - pkin(2);
t234 = t270 * t276 + t301;
t210 = t235 * pkin(4) + t234;
t354 = t210 * t206 + t195 * t321 + (-t195 * t275 - t183) * t279;
t350 = -t283 * t224 - t225 * t280;
t255 = -t331 + t334;
t351 = t255 * t254;
t182 = -pkin(8) * t232 + qJD(4) * t350 - t351;
t353 = -t279 * t182 + t282 * t183 + t210 * t296;
t261 = (-pkin(7) - qJ(3)) * t277;
t271 = t278 * pkin(7);
t262 = qJ(3) * t278 + t271;
t294 = -t261 * t280 - t262 * t283;
t349 = t294 * qJD(4) - t301 * t256;
t348 = qJD(5) - t275;
t322 = qJD(4) * t283;
t347 = -t255 * t318 + (qJD(3) * t277 + qJD(4) * t262) * t280 - qJD(3) * t331 - t261 * t322;
t346 = t278 * t281 * MDP(7) + t284 * MDP(6);
t345 = pkin(1) * t284;
t344 = pkin(4) * t237;
t343 = pkin(4) * t248;
t247 = t255 * qJD(4);
t342 = pkin(8) * t247;
t341 = pkin(8) * t256;
t338 = t195 * t282;
t215 = t255 * t282 + t256 * t279;
t197 = -t215 * qJD(5) - t247 * t282 - t248 * t279;
t216 = -t255 * t279 + t256 * t282;
t267 = t281 * t316;
t217 = pkin(4) * t232 + t267;
t330 = t210 * t197 + t217 * t216;
t198 = t216 * qJD(5) - t247 * t279 + t282 * t248;
t329 = t210 * t198 + t217 * t215;
t328 = -t234 * t247 + t256 * t267;
t327 = t234 * t248 + t255 * t267;
t317 = t281 * t339;
t194 = -pkin(8) * t237 + t350;
t193 = qJD(4) * pkin(4) + t194;
t313 = -pkin(4) * t275 - t193;
t265 = t284 * t339 + qJD(3);
t311 = t265 * t324;
t310 = t324 * t284;
t306 = t324 * qJD(3);
t305 = t276 * t319;
t304 = t276 * t317;
t244 = t248 * pkin(8);
t300 = -qJD(5) * (t261 * t283 - t262 * t280 - t341) + t244 + t347;
t252 = t255 * pkin(8);
t299 = qJD(5) * (-t252 - t294) - t342 - t349;
t269 = pkin(1) * t281 + qJ(3);
t249 = (-pkin(7) - t269) * t277;
t250 = t269 * t278 + t271;
t295 = -t249 * t280 - t250 * t283;
t233 = t255 * pkin(4) + t270;
t292 = -t319 + t343;
t289 = (-t186 * t215 - t197 * t206 + t198 * t296 + t216 * t287) * MDP(19) + (t186 * t216 - t197 * t296) * MDP(18) + (-t231 * t255 - t232 * t256 + t235 * t247 - t237 * t248) * MDP(12) + (t231 * t256 - t237 * t247) * MDP(11) + (MDP(20) * t197 - MDP(21) * t198) * t275 + (-t247 * MDP(13) - t248 * MDP(14)) * qJD(4);
t288 = t249 * t322 + t265 * t331 + (-qJD(4) * t250 - t265 * t277) * t280;
t286 = t295 * qJD(4) - t256 * t265;
t263 = t277 * t267;
t260 = t270 - t345;
t257 = -t276 * pkin(2) + t301;
t226 = t317 + t343;
t222 = t233 - t345;
t204 = -t252 - t295;
t203 = t249 * t283 - t250 * t280 - t341;
t189 = t286 + t342;
t188 = -t244 + t288;
t1 = [(t276 * t311 + t355) * MDP(9) + (t286 * qJD(4) + t260 * t232 + t235 * t317 + t327) * MDP(16) + (-t226 * t296 + t222 * t186 - (t188 * t282 + t189 * t279 + (t203 * t282 - t204 * t279) * qJD(5)) * t275 + t330) * MDP(24) + (-t288 * qJD(4) + t260 * t231 + t237 * t317 + t328) * MDP(17) + (t258 * t311 + t269 * t355 + (t257 + (-pkin(2) - t345) * qJD(1)) * t317) * MDP(10) + t289 + (t226 * t206 - t222 * t287 + (-t188 * t279 + t189 * t282 + (-t203 * t279 - t204 * t282) * qJD(5)) * t275 + t329) * MDP(23) + (-t267 - t304) * MDP(5) + (t277 * t304 + t263) * MDP(8) + t346 * (-qJD(1) - t276) * t339; (qJD(4) * t349 + t270 * t232 - t235 * t319 + t327) * MDP(16) + (t233 * t186 + (t299 * t279 + t300 * t282) * t275 - t292 * t296 + t330) * MDP(24) + (-t233 * t287 + (t300 * t279 - t299 * t282) * t275 + t292 * t206 + t329) * MDP(23) + t289 + ((-t310 * t340 + t306) * t276 + t355) * MDP(9) + (-t267 + t305) * MDP(5) + (-t277 * t305 + t263) * MDP(8) + (qJD(4) * t347 + t270 * t231 - t237 * t319 + t328) * MDP(17) + (t258 * t306 + qJ(3) * t355 + ((-pkin(2) * qJD(2) - t257) * t281 - t258 * t310) * t340) * MDP(10) + t346 * (-qJD(2) + t276) * t340; (-t324 * t276 * t258 + t267) * MDP(10) + t259 * MDP(17) + (-t287 - t337) * MDP(23) + (t186 - t336) * MDP(24) - t324 * MDP(9) * t276 ^ 2 + (0.2e1 * t237 * MDP(16) + (-t235 - t315) * MDP(17)) * qJD(4); t237 * t235 * MDP(11) + (-t235 ^ 2 + t237 ^ 2) * MDP(12) + (t259 + (t235 - t315) * qJD(4)) * MDP(13) + (-t234 * t237 - t291) * MDP(16) + (t234 * t235 + t351) * MDP(17) + (-t206 * t344 - (-t194 * t279 - t338) * t275 + (t313 * t279 - t338) * qJD(5) + t353) * MDP(23) + (t296 * t344 + (t313 * qJD(5) + t194 * t275 - t182) * t282 + t354) * MDP(24) + t356; (t348 * (-t193 * t279 - t338) + t353) * MDP(23) + ((-t193 * t348 - t182) * t282 + t354) * MDP(24) + t356;];
tauc = t1;
