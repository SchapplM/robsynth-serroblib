% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:41
% EndTime: 2019-03-09 01:42:45
% DurationCPUTime: 2.35s
% Computational Cost: add. (1734->277), mult. (4161->350), div. (0->0), fcn. (3000->8), ass. (0->125)
t277 = cos(pkin(10));
t345 = cos(qJ(4));
t314 = t345 * t277;
t303 = qJD(1) * t314;
t275 = sin(pkin(10));
t280 = sin(qJ(4));
t325 = qJD(1) * t280;
t313 = t275 * t325;
t248 = -t303 + t313;
t281 = cos(qJ(6));
t279 = sin(qJ(6));
t322 = qJD(4) * t279;
t231 = -t248 * t281 + t322;
t258 = t275 * t345 + t277 * t280;
t352 = t258 * qJD(1);
t354 = qJD(6) + t352;
t357 = t231 * t354;
t233 = qJD(4) * t281 + t248 * t279;
t305 = t354 * t233;
t308 = qJD(1) * (t275 ^ 2 + t277 ^ 2);
t356 = MDP(7) * t308;
t265 = qJD(4) * t303;
t321 = qJD(4) * t280;
t312 = t275 * t321;
t242 = qJD(1) * t312 - t265;
t235 = t281 * t242;
t307 = t279 * t354;
t355 = -t307 * t354 - t235;
t269 = sin(pkin(9)) * pkin(1) + qJ(3);
t260 = t269 * qJD(1);
t272 = t277 * qJD(2);
t341 = pkin(7) * qJD(1);
t229 = t272 + (-t260 - t341) * t275;
t240 = qJD(2) * t275 + t260 * t277;
t230 = t277 * t341 + t240;
t329 = t229 * t345 - t230 * t280;
t351 = qJD(5) - t329;
t202 = t229 * t280 + t230 * t345;
t199 = -qJD(4) * qJ(5) - t202;
t343 = pkin(5) * t248;
t190 = -t199 - t343;
t346 = pkin(4) + pkin(8);
t350 = t346 * t242 + (t190 - t202 + t343) * t354;
t253 = t258 * qJD(4);
t257 = t275 * t280 - t314;
t319 = qJD(6) * t281;
t292 = t253 * t279 + t257 * t319;
t331 = t257 * t279;
t349 = t242 * t331 - t292 * t354;
t342 = pkin(7) + t269;
t254 = t342 * t275;
t255 = t342 * t277;
t310 = qJD(4) * t345;
t311 = qJD(3) * t345;
t203 = (qJD(3) * t275 + qJD(4) * t255) * t280 + t254 * t310 - t277 * t311;
t348 = t248 ^ 2;
t347 = t352 ^ 2;
t243 = qJD(1) * t253;
t344 = pkin(4) * t243;
t340 = qJ(5) * t242;
t339 = qJ(5) * t248;
t259 = -cos(pkin(9)) * pkin(1) - pkin(3) * t277 - pkin(2);
t286 = -qJ(5) * t258 + t259;
t205 = t257 * t346 + t286;
t338 = t205 * t242;
t247 = qJD(1) * t259 + qJD(3);
t284 = -qJ(5) * t352 + t247;
t209 = pkin(4) * t248 + t284;
t337 = t209 * t352;
t320 = qJD(6) * t279;
t328 = t243 * t279 + t248 * t319;
t211 = -qJD(4) * t320 + t328;
t336 = t211 * t281;
t335 = t231 * t248;
t334 = t233 * t248;
t333 = t242 * t279;
t332 = t248 * t352;
t252 = -t277 * t310 + t312;
t330 = t211 * t258 - t233 * t252;
t324 = qJD(4) * t203;
t293 = -t254 * t280 + t255 * t345;
t204 = qJD(3) * t258 + qJD(4) * t293;
t323 = qJD(4) * t204;
t317 = pkin(5) * t352 + t351;
t309 = qJD(3) * t325;
t306 = t281 * t354;
t302 = qJD(1) * t311;
t304 = -t229 * t310 + t230 * t321 + t275 * t309 - t277 * t302;
t189 = t229 * t321 + t230 * t310 + t275 * t302 + t277 * t309;
t215 = t254 * t345 + t255 * t280;
t186 = -qJD(4) * t346 + t317;
t195 = t248 * t346 + t284;
t181 = t186 * t281 - t195 * t279;
t182 = t186 * t279 + t195 * t281;
t236 = t281 * t243;
t212 = qJD(6) * t233 - t236;
t300 = -t212 * t258 + t231 * t252;
t299 = (-t260 * t275 + t272) * t275 - t240 * t277;
t298 = -t242 * t257 + t253 * t352;
t297 = -t243 * t258 + t248 * t252;
t296 = -qJD(5) * t352 + t340;
t295 = qJ(5) * t252 - qJD(5) * t258;
t291 = -t257 * t235 + (t253 * t281 - t257 * t320) * t354;
t290 = qJD(4) * t202 - t189;
t187 = -qJD(4) * qJD(5) + t304;
t183 = -pkin(5) * t243 - t187;
t289 = t183 + (t346 * t354 + t339) * t354;
t207 = pkin(5) * t258 + t215;
t288 = t183 * t257 + t190 * t253 + t207 * t242;
t287 = -t306 * t354 + t333;
t245 = qJD(4) * t248;
t220 = t242 * t258;
t218 = pkin(4) * t352 + t339;
t214 = pkin(4) * t257 + t286;
t210 = pkin(4) * t253 + t295;
t208 = -pkin(5) * t257 + t293;
t200 = t296 + t344;
t198 = -qJD(4) * pkin(4) + t351;
t196 = t253 * t346 + t295;
t194 = -pkin(5) * t252 + t204;
t193 = -pkin(5) * t253 - t203;
t188 = t243 * t346 + t296;
t185 = -pkin(5) * t242 + t189;
t184 = t281 * t185;
t1 = [(-t252 * t352 - t220) * MDP(9) + (t297 - t298) * MDP(10) + (t243 * t259 + t247 * t253 - t323) * MDP(14) + (-t242 * t259 - t247 * t252 + t324) * MDP(15) + (t187 * t257 + t189 * t258 - t198 * t252 + t199 * t253 + t203 * t248 + t204 * t352 - t215 * t242 - t243 * t293) * MDP(16) + (-t200 * t257 - t209 * t253 - t210 * t248 - t214 * t243 + t323) * MDP(17) + (-t200 * t258 + t209 * t252 - t210 * t352 + t214 * t242 - t324) * MDP(18) + (-t187 * t293 + t189 * t215 + t198 * t204 + t199 * t203 + t200 * t214 + t209 * t210) * MDP(19) + (t211 * t331 + t233 * t292) * MDP(20) + ((-t231 * t279 + t233 * t281) * t253 + (t336 - t212 * t279 + (-t231 * t281 - t233 * t279) * qJD(6)) * t257) * MDP(21) + (t330 - t349) * MDP(22) + (t291 + t300) * MDP(23) + (-t252 * t354 - t220) * MDP(24) + (-t181 * t252 + t184 * t258 + t193 * t231 + t208 * t212 + (-t188 * t258 - t196 * t354 + t338) * t279 + (t194 * t354 - t288) * t281 + ((-t205 * t281 - t207 * t279) * t354 - t182 * t258 + t190 * t331) * qJD(6)) * MDP(25) + (t182 * t252 + t193 * t233 + t208 * t211 + (-(qJD(6) * t207 + t196) * t354 + t338 - (qJD(6) * t186 + t188) * t258 + t190 * qJD(6) * t257) * t281 + (-(-qJD(6) * t205 + t194) * t354 - (-qJD(6) * t195 + t185) * t258 + t288) * t279) * MDP(26) + (-MDP(11) * t252 - MDP(12) * t253) * qJD(4) + (0.2e1 * t356 + (t269 * t308 - t299) * MDP(8)) * qJD(3); (t297 + t298) * MDP(16) + (-t187 * t258 + t189 * t257 + t198 * t253 + t199 * t252) * MDP(19) + (t291 - t300) * MDP(25) + (t330 + t349) * MDP(26) + ((-MDP(14) + MDP(17)) * t253 + (MDP(15) - MDP(18)) * t252) * qJD(4); t265 * MDP(15) + (-t347 - t348) * MDP(16) + (t242 + t245) * MDP(18) + (t344 + t340 - t199 * t248 + (-qJD(5) - t198) * t352) * MDP(19) + (t287 + t335) * MDP(25) + (t334 - t355) * MDP(26) + ((-t248 - t313) * MDP(15) + (0.2e1 * MDP(14) - 0.2e1 * MDP(17)) * t352) * qJD(4) + (t299 * MDP(8) - t356) * qJD(1); MDP(9) * t332 + (t347 - t348) * MDP(10) + (t265 + (t248 - t313) * qJD(4)) * MDP(11) + (-t247 * t352 + t290) * MDP(14) + (qJD(4) * t329 + t247 * t248 + t304) * MDP(15) + (pkin(4) * t242 - qJ(5) * t243 + (-t199 - t202) * t352 + (t198 - t351) * t248) * MDP(16) + (t218 * t248 - t290 + t337) * MDP(17) + (-t209 * t248 + t218 * t352 + (0.2e1 * qJD(5) - t329) * qJD(4) - t304) * MDP(18) + (-pkin(4) * t189 - qJ(5) * t187 - t198 * t202 - t199 * t351 - t209 * t218) * MDP(19) + (-t279 * t305 + t336) * MDP(20) + ((-t212 - t305) * t281 + (-t211 + t357) * t279) * MDP(21) + (t334 + t355) * MDP(22) + (t287 - t335) * MDP(23) + t354 * t248 * MDP(24) + (qJ(5) * t212 + t181 * t248 + t317 * t231 + t289 * t279 + t281 * t350) * MDP(25) + (qJ(5) * t211 - t182 * t248 + t317 * t233 - t279 * t350 + t289 * t281) * MDP(26); (-t242 + t245) * MDP(16) - MDP(17) * t332 + (-qJD(4) ^ 2 - t347) * MDP(18) + (qJD(4) * t199 + t189 + t337) * MDP(19) + (-qJD(4) * t231 - t235) * MDP(25) + (-qJD(4) * t233 + t333) * MDP(26) + (-MDP(25) * t307 - MDP(26) * t306) * t354; t233 * t231 * MDP(20) + (-t231 ^ 2 + t233 ^ 2) * MDP(21) + (t328 + t357) * MDP(22) + (t236 + t305) * MDP(23) - t242 * MDP(24) + (t182 * t354 - t188 * t279 - t190 * t233 + t184) * MDP(25) + (t181 * t354 - t185 * t279 - t188 * t281 + t190 * t231) * MDP(26) + (-MDP(22) * t322 - MDP(23) * t233 - MDP(25) * t182 - MDP(26) * t181) * qJD(6);];
tauc  = t1;
