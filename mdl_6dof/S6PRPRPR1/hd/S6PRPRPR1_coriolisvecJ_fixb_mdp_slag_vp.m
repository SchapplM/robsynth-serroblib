% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:28:08
% EndTime: 2019-03-08 19:28:11
% DurationCPUTime: 2.14s
% Computational Cost: add. (2064->268), mult. (5309->389), div. (0->0), fcn. (4219->12), ass. (0->136)
t285 = cos(pkin(11));
t289 = sin(qJ(2));
t283 = sin(pkin(6));
t338 = qJD(1) * t283;
t326 = t289 * t338;
t266 = t285 * t326;
t282 = sin(pkin(11));
t292 = cos(qJ(2));
t325 = t292 * t338;
t243 = t282 * t325 + t266;
t235 = qJD(2) * t243;
t247 = (t282 * t289 - t285 * t292) * t283;
t288 = sin(qJ(4));
t291 = cos(qJ(4));
t363 = (t288 ^ 2 - t291 ^ 2) * MDP(7);
t336 = qJD(4) * t288;
t362 = pkin(4) * t336 - t243;
t286 = cos(pkin(6));
t270 = qJD(1) * t286 + qJD(3);
t264 = qJD(2) * pkin(2) + t325;
t233 = t282 * t264 + t266;
t319 = t233 + (pkin(8) + qJ(5)) * qJD(2);
t217 = t291 * t270 - t319 * t288;
t218 = t288 * t270 + t319 * t291;
t245 = qJD(2) * t247;
t236 = qJD(1) * t245;
t315 = qJD(2) * qJD(5) - t236;
t361 = -t218 * qJD(4) - t315 * t288;
t197 = t217 * qJD(4) + t315 * t291;
t281 = sin(pkin(12));
t284 = cos(pkin(12));
t186 = t197 * t281 - t284 * t361;
t337 = qJD(2) * t288;
t347 = t284 * t291;
t253 = qJD(2) * t347 - t281 * t337;
t252 = qJD(6) - t253;
t261 = t281 * t291 + t284 * t288;
t255 = t261 * qJD(2);
t273 = pkin(4) * t281 + pkin(9);
t360 = (pkin(4) * t337 + pkin(5) * t255 - pkin(9) * t253 + qJD(6) * t273) * t252 + t186;
t187 = t284 * t197 + t361 * t281;
t212 = qJD(4) * pkin(4) + t217;
t356 = t218 * t281;
t192 = t212 * t284 - t356;
t190 = -qJD(4) * pkin(5) - t192;
t265 = t282 * t326;
t232 = t264 * t285 - t265;
t327 = -pkin(4) * t291 - pkin(3);
t226 = t327 * qJD(2) + qJD(5) - t232;
t200 = -pkin(5) * t253 - pkin(9) * t255 + t226;
t260 = t281 * t288 - t347;
t358 = pkin(2) * t285;
t305 = t327 - t358;
t222 = pkin(5) * t260 - pkin(9) * t261 + t305;
t257 = t260 * qJD(4);
t274 = pkin(2) * t282 + pkin(8);
t345 = qJ(5) + t274;
t259 = t345 * t291;
t322 = t345 * t288;
t224 = t284 * t259 - t281 * t322;
t254 = t261 * qJD(4);
t249 = qJD(2) * t254;
t313 = t186 * t261 - t224 * t249;
t321 = qJD(4) * t345;
t240 = qJD(5) * t291 - t288 * t321;
t246 = t285 * t325 - t265;
t298 = -qJD(5) * t288 - t291 * t321;
t342 = t284 * t240 + t260 * t246 + t281 * t298;
t359 = -(qJD(6) * t200 + t187) * t260 - t190 * t257 + (-qJD(6) * t222 - t342) * t252 + t313;
t287 = sin(qJ(6));
t334 = qJD(6) * t287;
t331 = qJD(2) * qJD(4);
t323 = t291 * t331;
t324 = t288 * t331;
t250 = -t281 * t324 + t284 * t323;
t290 = cos(qJ(6));
t332 = t290 * qJD(4);
t340 = qJD(6) * t332 + t290 * t250;
t215 = -t255 * t334 + t340;
t357 = t215 * t287;
t355 = t222 * t249;
t348 = t255 * t287;
t237 = -t332 + t348;
t354 = t237 * t252;
t353 = t237 * t255;
t239 = qJD(4) * t287 + t255 * t290;
t352 = t239 * t252;
t351 = t239 * t255;
t350 = t249 * t287;
t349 = t250 * t287;
t210 = t284 * t218;
t241 = t290 * t249;
t344 = t215 * t260 + t239 * t254;
t193 = t281 * t212 + t210;
t343 = t240 * t281 - t261 * t246 - t284 * t298;
t341 = pkin(5) * t254 + pkin(9) * t257 + t362;
t335 = qJD(6) * t261;
t333 = t288 * MDP(11);
t329 = t261 * t350;
t328 = t261 * t241;
t228 = pkin(4) * t324 + t235;
t320 = t252 * t290;
t230 = -qJD(2) * pkin(3) - t232;
t318 = -qJD(2) * t230 + t236;
t191 = qJD(4) * pkin(9) + t193;
t185 = t191 * t290 + t200 * t287;
t312 = t191 * t287 - t200 * t290;
t248 = (t282 * t292 + t285 * t289) * t283;
t229 = t248 * t291 + t286 * t288;
t308 = -t248 * t288 + t286 * t291;
t202 = t284 * t229 + t281 * t308;
t311 = t202 * t290 + t247 * t287;
t310 = -t202 * t287 + t247 * t290;
t216 = t239 * qJD(6) + t349;
t309 = -t216 * t260 - t237 * t254;
t303 = t241 + (t253 * t287 - t334) * t252;
t302 = t257 * t287 - t290 * t335;
t301 = t257 * t290 + t261 * t334;
t293 = qJD(4) ^ 2;
t300 = t274 * t293;
t299 = qJD(4) * (qJD(2) * (-pkin(3) - t358) + t230 + t246);
t195 = t217 * t284 - t356;
t297 = -t273 * t249 + (t190 + t195) * t252;
t294 = qJD(2) ^ 2;
t275 = -pkin(4) * t284 - pkin(5);
t244 = qJD(2) * t248;
t223 = t259 * t281 + t284 * t322;
t206 = t308 * qJD(4) - t245 * t291;
t205 = -t229 * qJD(4) + t245 * t288;
t201 = t229 * t281 - t284 * t308;
t199 = pkin(5) * t249 - pkin(9) * t250 + t228;
t198 = t290 * t199;
t194 = t217 * t281 + t210;
t189 = t205 * t281 + t206 * t284;
t188 = -t284 * t205 + t206 * t281;
t1 = [(-t232 * t244 - t233 * t245 + t235 * t247 - t236 * t248) * MDP(5) + (t188 * t255 + t189 * t253 + t201 * t250 - t202 * t249) * MDP(13) + (t186 * t201 + t187 * t202 - t188 * t192 + t189 * t193 + t226 * t244 + t228 * t247) * MDP(14) + ((-t311 * qJD(6) - t189 * t287 + t244 * t290) * t252 + t310 * t249 + t188 * t237 + t201 * t216) * MDP(20) + (-(t310 * qJD(6) + t189 * t290 + t244 * t287) * t252 - t311 * t249 + t188 * t239 + t201 * t215) * MDP(21) + (-MDP(3) * t289 - MDP(4) * t292) * t294 * t283 + (MDP(11) * t205 - MDP(12) * t206) * qJD(4) + ((-t244 * t291 + t247 * t336) * MDP(11) + (qJD(4) * t247 * t291 + t244 * t288) * MDP(12)) * qJD(2); (t232 * t243 - t233 * t246 + (-t235 * t285 - t236 * t282) * pkin(2)) * MDP(5) - 0.2e1 * t331 * t363 + (-t187 * t260 + t192 * t257 - t193 * t254 + t223 * t250 + t342 * t253 + t343 * t255 + t313) * MDP(13) + (t186 * t223 + t187 * t224 - t343 * t192 + t342 * t193 + t362 * t226 + t228 * t305) * MDP(14) + (t215 * t261 * t290 - t301 * t239) * MDP(15) + (-(-t237 * t290 - t239 * t287) * t257 + (-t357 - t216 * t290 + (t237 * t287 - t239 * t290) * qJD(6)) * t261) * MDP(16) + (-t301 * t252 + t328 + t344) * MDP(17) + (t302 * t252 + t309 - t329) * MDP(18) + (t249 * t260 + t252 * t254) * MDP(19) + (-t312 * t254 + t198 * t260 + t223 * t216 + t343 * t237 + (t355 + t341 * t252 + (t190 * t261 - t191 * t260 - t224 * t252) * qJD(6)) * t290 + t359 * t287) * MDP(20) + (-t185 * t254 + t223 * t215 + t343 * t239 + (-t355 - (-qJD(6) * t191 + t199) * t260 - t190 * t335 + (qJD(6) * t224 - t341) * t252) * t287 + t359 * t290) * MDP(21) + (-t300 * MDP(11) + t299 * MDP(12) + t293 * MDP(8)) * t291 + (t299 * MDP(11) + t300 * MDP(12) + 0.2e1 * MDP(6) * t323 - t293 * MDP(9)) * t288; (-t249 * t261 + t250 * t260 - t253 * t257 + t254 * t255) * MDP(13) + (t186 * t260 + t187 * t261 - t192 * t254 - t193 * t257) * MDP(14) + (-t309 - t329) * MDP(20) + (-t328 + t344) * MDP(21) + (-t291 * MDP(12) - t333) * t293 + (t302 * MDP(20) + t301 * MDP(21)) * t252; t294 * t363 + t318 * t333 + ((t193 - t194) * t255 + (t192 - t195) * t253 + (-t249 * t281 - t250 * t284) * pkin(4)) * MDP(13) + (t192 * t194 - t193 * t195 + (-t186 * t284 + t187 * t281 - t226 * t337) * pkin(4)) * MDP(14) + (t239 * t320 + t357) * MDP(15) + ((t215 - t354) * t290 + (-t216 - t352) * t287) * MDP(16) + (t252 * t320 + t350 - t351) * MDP(17) + (t303 + t353) * MDP(18) - t252 * t255 * MDP(19) + (-t194 * t237 + t275 * t216 + t255 * t312 + t297 * t287 - t360 * t290) * MDP(20) + (t185 * t255 - t194 * t239 + t275 * t215 + t360 * t287 + t297 * t290) * MDP(21) + (-t288 * t294 * MDP(6) + t318 * MDP(12)) * t291; (-t253 ^ 2 - t255 ^ 2) * MDP(13) + (t192 * t255 - t193 * t253 + t228) * MDP(14) + (t303 - t353) * MDP(20) + (-t252 ^ 2 * t290 - t350 - t351) * MDP(21); t239 * t237 * MDP(15) + (-t237 ^ 2 + t239 ^ 2) * MDP(16) + (t340 + t354) * MDP(17) + (-t349 + t352) * MDP(18) + t249 * MDP(19) + (t185 * t252 - t187 * t287 - t190 * t239 + t198) * MDP(20) + (-t187 * t290 + t190 * t237 - t199 * t287 - t252 * t312) * MDP(21) + (-MDP(17) * t348 - t239 * MDP(18) - t185 * MDP(20) + t312 * MDP(21)) * qJD(6);];
tauc  = t1;
