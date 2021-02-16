% Calculate Coriolis joint torque vector for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:30:06
% EndTime: 2021-01-15 20:30:16
% DurationCPUTime: 2.90s
% Computational Cost: add. (2848->342), mult. (7254->445), div. (0->0), fcn. (5009->6), ass. (0->145)
t298 = sin(qJ(4));
t300 = cos(qJ(4));
t297 = sin(pkin(8));
t299 = sin(qJ(2));
t301 = cos(qJ(2));
t365 = cos(pkin(8));
t278 = t297 * t301 + t299 * t365;
t339 = qJD(1) * t278;
t251 = qJD(2) * t298 + t300 * t339;
t357 = t251 * t298;
t333 = qJD(1) * qJD(2);
t375 = -0.2e1 * t333;
t374 = MDP(5) * (t299 ^ 2 - t301 ^ 2);
t325 = t365 * t301;
t287 = qJD(1) * t325;
t337 = qJD(1) * t299;
t265 = t297 * t337 - t287;
t263 = qJD(4) + t265;
t373 = t263 * t357;
t293 = -pkin(2) * t301 - pkin(1);
t338 = qJD(1) * t293;
t282 = qJD(3) + t338;
t220 = pkin(3) * t265 - pkin(7) * t339 + t282;
t367 = -qJ(3) - pkin(6);
t284 = t367 * t299;
t280 = qJD(1) * t284;
t366 = qJD(2) * pkin(2);
t274 = t280 + t366;
t285 = t367 * t301;
t281 = qJD(1) * t285;
t326 = t365 * t281;
t241 = t297 * t274 - t326;
t236 = qJD(2) * pkin(7) + t241;
t204 = t220 * t298 + t236 * t300;
t267 = t278 * qJD(2);
t261 = qJD(1) * t267;
t328 = t299 * t333;
t288 = pkin(2) * t328;
t286 = t297 * t328;
t308 = t287 * qJD(2) - t286;
t219 = t261 * pkin(3) - pkin(7) * t308 + t288;
t212 = t300 * t219;
t327 = qJD(2) * t367;
t264 = qJD(3) * t301 + t299 * t327;
t258 = t264 * qJD(1);
t307 = -qJD(3) * t299 + t301 * t327;
t259 = t307 * qJD(1);
t215 = t258 * t365 + t297 * t259;
t305 = -qJD(4) * t204 - t215 * t298 + t212;
t334 = t300 * qJD(2);
t336 = qJD(4) * t298;
t314 = qJD(4) * t334 + t300 * t308 - t336 * t339;
t304 = -qJ(5) * t314 + t305;
t369 = pkin(4) * t261;
t190 = -qJD(5) * t251 + t304 + t369;
t249 = t298 * t339 - t334;
t199 = -qJ(5) * t249 + t204;
t364 = t199 * t263;
t372 = t190 + t364;
t217 = t251 * qJD(4) + t298 * t308;
t371 = t251 ^ 2;
t370 = pkin(2) * t299;
t368 = t249 * pkin(4);
t363 = t314 * t298;
t362 = t249 * t263;
t361 = t249 * t265;
t360 = t249 * t339;
t359 = t251 * t263;
t358 = t251 * t339;
t356 = t261 * t298;
t354 = t265 * t298;
t353 = t265 * t300;
t352 = t278 * t298;
t351 = t278 * t300;
t271 = t297 * t281;
t302 = qJD(2) ^ 2;
t350 = t299 * t302;
t247 = t297 * t284 - t285 * t365;
t244 = t300 * t247;
t256 = t300 * t261;
t349 = t301 * t302;
t303 = qJD(1) ^ 2;
t348 = t301 * t303;
t290 = pkin(2) * t297 + pkin(7);
t347 = qJ(5) + t290;
t203 = t300 * t220 - t236 * t298;
t198 = -qJ(5) * t251 + t203;
t194 = pkin(4) * t263 + t198;
t346 = t194 - t198;
t332 = pkin(2) * t337;
t230 = pkin(3) * t339 + pkin(7) * t265 + t332;
t224 = t300 * t230;
t243 = t280 * t365 + t271;
t322 = qJD(4) * t347;
t345 = pkin(4) * t339 + qJ(5) * t353 + t300 * t322 + t224 + (qJD(5) - t243) * t298;
t342 = t298 * t230 + t300 * t243;
t344 = qJ(5) * t354 - qJD(5) * t300 + t298 * t322 + t342;
t335 = qJD(4) * t300;
t343 = -t298 * t217 - t249 * t335;
t311 = -t297 * t299 + t325;
t239 = -pkin(3) * t311 - pkin(7) * t278 + t293;
t341 = t298 * t239 + t244;
t331 = t299 * t366;
t229 = t264 * t365 + t297 * t307;
t270 = t311 * qJD(2);
t231 = pkin(3) * t267 - pkin(7) * t270 + t331;
t330 = t300 * t229 + t298 * t231 + t239 * t335;
t329 = t278 * t335;
t324 = qJD(5) + t368;
t323 = pkin(1) * t375;
t214 = t258 * t297 - t365 * t259;
t228 = t264 * t297 - t365 * t307;
t242 = t280 * t297 - t326;
t246 = -t365 * t284 - t285 * t297;
t321 = t263 * t300;
t320 = 0.2e1 * t339;
t319 = t300 * t215 + t298 * t219 + t220 * t335 - t236 * t336;
t291 = -pkin(2) * t365 - pkin(3);
t240 = t274 * t365 + t271;
t318 = t214 * t278 - t247 * t261;
t310 = qJ(5) * t217 - t319;
t191 = -qJD(5) * t249 - t310;
t317 = -t194 * t263 + t191;
t316 = -qJ(5) * t270 - qJD(5) * t278;
t315 = t256 + (-t336 - t354) * t263;
t202 = pkin(4) * t217 + t214;
t313 = t270 * t298 + t329;
t312 = t270 * t300 - t278 * t336;
t235 = -qJD(2) * pkin(3) - t240;
t309 = t235 * t263 - t290 * t261;
t283 = -t300 * pkin(4) + t291;
t276 = t347 * t300;
t275 = t347 * t298;
t248 = t249 ^ 2;
t234 = t300 * t239;
t226 = pkin(4) * t352 + t246;
t225 = t300 * t231;
t210 = -pkin(4) * t354 + t242;
t207 = t235 + t324;
t206 = pkin(4) * t313 + t228;
t205 = -qJ(5) * t352 + t341;
t200 = -pkin(4) * t311 - qJ(5) * t351 - t247 * t298 + t234;
t193 = -qJ(5) * t329 + (-qJD(4) * t247 + t316) * t298 + t330;
t192 = pkin(4) * t267 - t229 * t298 + t225 + t316 * t300 + (-t244 + (qJ(5) * t278 - t239) * t298) * qJD(4);
t1 = [0.2e1 * t301 * MDP(4) * t328 + t374 * t375 + MDP(6) * t349 - MDP(7) * t350 + (-pkin(6) * t349 + t299 * t323) * MDP(9) + (pkin(6) * t350 + t301 * t323) * MDP(10) + (t261 * t293 + t267 * t282 + (-t228 + (-qJD(1) * t311 + t265) * t370) * qJD(2)) * MDP(11) + (t282 * t270 - t293 * t286 + (t287 * t293 + t320 * t370 - t229) * qJD(2)) * MDP(12) + (t215 * t311 + t228 * t339 - t229 * t265 - t240 * t270 - t241 * t267 + t246 * t308 + t318) * MDP(13) + (t214 * t246 + t215 * t247 - t228 * t240 + t229 * t241 + (t282 + t338) * t331) * MDP(14) + (t251 * t312 + t314 * t351) * MDP(15) + ((-t249 * t300 - t357) * t270 + (-t363 - t217 * t300 + (t249 * t298 - t251 * t300) * qJD(4)) * t278) * MDP(16) + (t251 * t267 + t256 * t278 + t263 * t312 - t311 * t314) * MDP(17) + (t217 * t311 - t249 * t267 - t261 * t352 - t263 * t313) * MDP(18) + (-t261 * t311 + t263 * t267) * MDP(19) + ((-t247 * t335 + t225) * t263 + t234 * t261 - (-t236 * t335 + t212) * t311 + t203 * t267 + t228 * t249 + t246 * t217 + t235 * t329 + ((-qJD(4) * t239 - t229) * t263 - (-qJD(4) * t220 - t215) * t311 + t235 * t270 + t318) * t298) * MDP(20) + (-(-t247 * t336 + t330) * t263 - t341 * t261 + t319 * t311 - t204 * t267 + t228 * t251 + t246 * t314 + t214 * t351 + t312 * t235) * MDP(21) + (-t190 * t311 + t192 * t263 + t194 * t267 + t200 * t261 + t202 * t352 + t206 * t249 + t207 * t313 + t217 * t226) * MDP(22) + (t191 * t311 - t193 * t263 - t199 * t267 + t202 * t351 - t205 * t261 + t206 * t251 + t207 * t312 + t226 * t314) * MDP(23) + (-t192 * t251 - t193 * t249 - t200 * t314 - t205 * t217 + (-t194 * t300 - t199 * t298) * t270 + (-t190 * t300 - t191 * t298 + (t194 * t298 - t199 * t300) * qJD(4)) * t278) * MDP(24) + (t190 * t200 + t191 * t205 + t192 * t194 + t193 * t199 + t202 * t226 + t206 * t207) * MDP(25); -t299 * MDP(4) * t348 + t303 * t374 + (qJD(2) * t242 - t265 * t332 - t282 * t339 - t214) * MDP(11) + (t243 * qJD(2) + t282 * t265 - t332 * t339 - t215) * MDP(12) + ((t241 - t242) * t339 + (t243 - t240) * t265 + (-t297 * t261 - t308 * t365) * pkin(2)) * MDP(13) + (t240 * t242 - t241 * t243 + (-t214 * t365 + t215 * t297 - t282 * t337) * pkin(2)) * MDP(14) + (t251 * t321 + t363) * MDP(15) + ((t314 - t361) * t300 - t373 + t343) * MDP(16) + (t263 * t321 + t356 - t358) * MDP(17) + (t315 + t360) * MDP(18) - t263 * t339 * MDP(19) + (-t203 * t339 - t214 * t300 + t291 * t217 - t242 * t249 + (-t290 * t335 - t224) * t263 + (t243 * t263 + t309) * t298) * MDP(20) + (t204 * t339 + t214 * t298 + t291 * t314 - t242 * t251 + (t290 * t336 + t342) * t263 + t309 * t300) * MDP(21) + (-t194 * t339 - t202 * t300 - t210 * t249 + t217 * t283 - t261 * t275 - t345 * t263 + (t207 * t265 + (t207 + t368) * qJD(4)) * t298) * MDP(22) + (t207 * t353 + t199 * t339 + t202 * t298 - t210 * t251 + t314 * t283 - t261 * t276 + t344 * t263 + (pkin(4) * t357 + t207 * t300) * qJD(4)) * MDP(23) + (-t217 * t276 + t344 * t249 + t345 * t251 + t275 * t314 - t372 * t298 + t317 * t300) * MDP(24) + (-t190 * t275 + t191 * t276 + t202 * t283 + (pkin(4) * t336 - t210) * t207 - t344 * t199 - t345 * t194) * MDP(25) + (MDP(9) * t299 * t303 + MDP(10) * t348) * pkin(1); t320 * qJD(2) * MDP(11) + (-t286 + (t287 - t265) * qJD(2)) * MDP(12) + (-t265 ^ 2 - t339 ^ 2) * MDP(13) + (t240 * t339 + t241 * t265 + t288) * MDP(14) + ((-t314 - t361) * t300 + t373 + t343) * MDP(24) + (-t207 * t339 + t317 * t298 + t372 * t300) * MDP(25) + (MDP(21) + MDP(23)) * (-t263 ^ 2 * t300 - t356 - t358) + (MDP(20) + MDP(22)) * (t315 - t360); t251 * t249 * MDP(15) + (-t248 + t371) * MDP(16) + (t314 + t362) * MDP(17) + (-t217 + t359) * MDP(18) + t261 * MDP(19) + (t204 * t263 - t235 * t251 + t305) * MDP(20) + (t203 * t263 + t235 * t249 - t319) * MDP(21) + (0.2e1 * t369 + t364 + (-t207 - t324) * t251 + t304) * MDP(22) + (-pkin(4) * t371 + t198 * t263 + (qJD(5) + t207) * t249 + t310) * MDP(23) + (-pkin(4) * t314 - t249 * t346) * MDP(24) + (t346 * t199 + (-t207 * t251 + t190) * pkin(4)) * MDP(25); (t217 + t359) * MDP(22) + (t314 - t362) * MDP(23) + (-t248 - t371) * MDP(24) + (t194 * t251 + t199 * t249 + t202) * MDP(25);];
tauc = t1;
