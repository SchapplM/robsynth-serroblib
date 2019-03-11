% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:43
% EndTime: 2019-03-09 02:42:49
% DurationCPUTime: 2.59s
% Computational Cost: add. (2004->295), mult. (4762->385), div. (0->0), fcn. (3224->8), ass. (0->139)
t280 = sin(pkin(10));
t285 = sin(qJ(3));
t324 = t285 * qJD(1);
t282 = cos(pkin(10));
t287 = cos(qJ(3));
t337 = t282 * t287;
t253 = -qJD(1) * t337 + t280 * t324;
t286 = cos(qJ(6));
t284 = sin(qJ(6));
t329 = qJD(3) * t284;
t234 = -t286 * t253 + t329;
t262 = t280 * t287 + t282 * t285;
t256 = t262 * qJD(1);
t359 = qJD(6) + t256;
t360 = t234 * t359;
t236 = qJD(3) * t286 + t253 * t284;
t310 = t359 * t236;
t272 = sin(pkin(9)) * pkin(1) + pkin(7);
t334 = qJ(4) + t272;
t358 = MDP(12) + MDP(14);
t357 = MDP(5) * t285;
t356 = MDP(6) * (t285 ^ 2 - t287 ^ 2);
t355 = t284 * t359;
t309 = t334 * qJD(1);
t237 = t287 * qJD(2) - t309 * t285;
t238 = qJD(2) * t285 + t287 * t309;
t229 = t280 * t238;
t207 = t237 * t282 - t229;
t322 = -qJD(5) + t207;
t232 = qJD(3) * pkin(3) + t237;
t338 = t282 * t238;
t204 = t280 * t232 + t338;
t201 = -qJD(3) * qJ(5) - t204;
t348 = pkin(5) * t253;
t189 = -t201 - t348;
t206 = t237 * t280 + t338;
t319 = qJD(1) * qJD(3);
t316 = t285 * t319;
t265 = t280 * t316;
t315 = t287 * t319;
t247 = t282 * t315 - t265;
t273 = -pkin(3) * t282 - pkin(4);
t268 = -pkin(8) + t273;
t353 = t268 * t247 + (t189 - t206 + t348) * t359;
t255 = t262 * qJD(3);
t261 = t280 * t285 - t337;
t326 = qJD(6) * t286;
t300 = t255 * t284 + t261 * t326;
t339 = t261 * t284;
t352 = -t247 * t339 - t300 * t359;
t251 = t256 ^ 2;
t350 = pkin(4) + pkin(8);
t246 = qJD(1) * t255;
t349 = pkin(4) * t246;
t347 = pkin(5) * t256;
t318 = qJD(1) * qJD(4);
t226 = qJD(3) * t237 + t287 * t318;
t227 = -qJD(3) * t238 - t285 * t318;
t190 = t226 * t280 - t282 * t227;
t259 = t334 * t285;
t260 = t334 * t287;
t222 = t282 * t259 + t260 * t280;
t346 = t190 * t222;
t345 = t190 * t261;
t274 = -cos(pkin(9)) * pkin(1) - pkin(2);
t303 = -pkin(3) * t287 + t274;
t293 = -qJ(5) * t262 + t303;
t208 = t261 * t350 + t293;
t344 = t208 * t247;
t327 = qJD(6) * t284;
t332 = t284 * t246 + t253 * t326;
t214 = -qJD(3) * t327 + t332;
t343 = t214 * t286;
t342 = t234 * t253;
t341 = t236 * t253;
t340 = t247 * t284;
t288 = qJD(3) ^ 2;
t336 = t285 * t288;
t243 = t286 * t247;
t335 = t287 * t288;
t328 = qJD(3) * t285;
t258 = qJD(3) * t337 - t280 * t328;
t333 = t214 * t262 + t236 * t258;
t191 = t282 * t226 + t280 * t227;
t264 = qJD(1) * t274;
t325 = t256 * MDP(15);
t323 = t287 * MDP(11);
t321 = t347 - t322;
t276 = pkin(3) * t328;
t269 = pkin(3) * t316;
t314 = -qJ(5) * t247 + t269;
t313 = pkin(3) * t324 + qJ(5) * t253;
t203 = t232 * t282 - t229;
t312 = qJD(3) * t334;
t239 = qJD(4) * t287 - t285 * t312;
t240 = -qJD(4) * t285 - t287 * t312;
t209 = t239 * t280 - t282 * t240;
t311 = t286 * t359;
t308 = MDP(23) * t359;
t307 = qJD(5) - t203;
t188 = -qJD(3) * qJD(5) - t191;
t187 = -qJD(3) * t350 + t307 + t347;
t295 = t303 * qJD(1);
t252 = qJD(4) + t295;
t291 = -qJ(5) * t256 + t252;
t197 = t253 * t350 + t291;
t181 = t187 * t286 - t197 * t284;
t182 = t187 * t284 + t197 * t286;
t242 = t286 * t246;
t215 = qJD(6) * t236 - t242;
t305 = -t215 * t262 - t234 * t258;
t210 = t239 * t282 + t240 * t280;
t223 = -t259 * t280 + t260 * t282;
t302 = 0.2e1 * qJD(3) * t264;
t211 = pkin(4) * t253 + t291;
t301 = t211 * t256 + t190;
t299 = -qJD(5) * t256 + t314;
t298 = -qJ(5) * t258 - qJD(5) * t262 + t276;
t297 = t261 * t243 + (t255 * t286 - t261 * t327) * t359;
t183 = -pkin(5) * t246 - t188;
t296 = t183 + (-qJD(6) * t268 + t256 * t350 + t313) * t359;
t212 = pkin(5) * t262 + t222;
t294 = t183 * t261 + t189 * t255 - t212 * t247;
t290 = t190 * t262 + t209 * t256 - t210 * t253 + t222 * t247 - t223 * t246;
t270 = pkin(3) * t280 + qJ(5);
t249 = qJD(3) * t253;
t218 = pkin(4) * t261 + t293;
t216 = pkin(4) * t256 + t313;
t213 = -pkin(5) * t261 + t223;
t205 = pkin(4) * t255 + t298;
t199 = -qJD(3) * pkin(4) + t307;
t198 = t299 + t349;
t196 = -pkin(5) * t255 + t210;
t195 = pkin(5) * t258 + t209;
t192 = t255 * t350 + t298;
t186 = t246 * t350 + t299;
t185 = pkin(5) * t247 + t190;
t184 = t286 * t185;
t1 = [0.2e1 * t315 * t357 - 0.2e1 * t319 * t356 + MDP(7) * t335 - MDP(8) * t336 + (-t272 * t335 + t285 * t302) * MDP(10) + (t272 * t336 + t287 * t302) * MDP(11) + (-t191 * t261 - t203 * t258 - t204 * t255 + t290) * MDP(12) + (t346 + t191 * t223 - t203 * t209 + t204 * t210 + (t252 + t295) * t276) * MDP(13) + (t188 * t261 + t199 * t258 + t201 * t255 + t290) * MDP(14) + (qJD(3) * t209 - t198 * t261 - t205 * t253 - t211 * t255 - t218 * t246) * MDP(15) + (qJD(3) * t210 - t198 * t262 - t205 * t256 - t211 * t258 - t218 * t247) * MDP(16) + (-t188 * t223 + t198 * t218 + t199 * t209 - t201 * t210 + t205 * t211 + t346) * MDP(17) + (t214 * t339 + t236 * t300) * MDP(18) + ((-t234 * t284 + t236 * t286) * t255 + (t343 - t215 * t284 + (-t234 * t286 - t236 * t284) * qJD(6)) * t261) * MDP(19) + (t333 - t352) * MDP(20) + (t297 + t305) * MDP(21) + (t247 * t262 + t258 * t359) * MDP(22) + (t181 * t258 + t184 * t262 + t196 * t234 + t213 * t215 + (-t186 * t262 - t192 * t359 - t344) * t284 + (t195 * t359 - t294) * t286 + ((-t208 * t286 - t212 * t284) * t359 - t182 * t262 + t189 * t339) * qJD(6)) * MDP(23) + (-t182 * t258 + t196 * t236 + t213 * t214 + (-(qJD(6) * t212 + t192) * t359 - t344 - (qJD(6) * t187 + t186) * t262 + t189 * qJD(6) * t261) * t286 + (-(-qJD(6) * t208 + t195) * t359 - (-qJD(6) * t197 + t185) * t262 + t294) * t284) * MDP(24); (t191 * t262 - t203 * t255 + t204 * t258 + t345) * MDP(13) + (-t188 * t262 + t199 * t255 - t201 * t258 + t345) * MDP(17) + (t297 - t305) * MDP(23) + (t333 + t352) * MDP(24) + (-t285 * MDP(10) - t323) * t288 + (t255 * MDP(15) + t258 * MDP(16)) * qJD(3) + t358 * (-t246 * t262 + t247 * t261 - t253 * t258 + t255 * t256); ((t204 - t206) * t256 + (-t203 + t207) * t253 + (-t246 * t280 - t247 * t282) * pkin(3)) * MDP(12) + (t203 * t206 - t204 * t207 + (-t190 * t282 + t191 * t280 - t252 * t324) * pkin(3)) * MDP(13) + (-t246 * t270 + t247 * t273 + (-t201 - t206) * t256 + (t199 + t322) * t253) * MDP(14) + (-qJD(3) * t206 + t216 * t253 + t301) * MDP(15) + (-t211 * t253 + t216 * t256 + (0.2e1 * qJD(5) - t207) * qJD(3) + t191) * MDP(16) + (-t188 * t270 + t190 * t273 - t199 * t206 + t201 * t322 - t211 * t216) * MDP(17) + (-t284 * t310 + t343) * MDP(18) + ((-t215 - t310) * t286 + (-t214 + t360) * t284) * MDP(19) + (-t355 * t359 + t243 + t341) * MDP(20) + (-t311 * t359 - t340 - t342) * MDP(21) + t359 * t253 * MDP(22) + (t181 * t253 + t270 * t215 + t321 * t234 + t296 * t284 + t286 * t353) * MDP(23) + (-t182 * t253 + t270 * t214 + t321 * t236 - t284 * t353 + t296 * t286) * MDP(24) + (-t287 * t357 + t356) * qJD(1) ^ 2 + (-MDP(10) * t324 - qJD(1) * t323) * t264; (t203 * t256 + t204 * t253 + t269) * MDP(13) + (t249 + t265) * MDP(16) + (t349 - t201 * t253 + (-qJD(5) - t199) * t256 + t314) * MDP(17) + (-t340 + t342) * MDP(23) + (-t243 + t341) * MDP(24) + (MDP(24) * t355 - t286 * t308) * t359 + (-t325 + (-MDP(15) * t262 - MDP(16) * t337) * qJD(1)) * qJD(3) + t358 * (-t253 ^ 2 - t251); (t247 + t249) * MDP(14) - t253 * t325 + (-t251 - t288) * MDP(16) + (qJD(3) * t201 + t301) * MDP(17) + (-qJD(3) * t234 + t243) * MDP(23) + (-qJD(3) * t236 - t340) * MDP(24) + (-MDP(24) * t311 - t284 * t308) * t359; t236 * t234 * MDP(18) + (-t234 ^ 2 + t236 ^ 2) * MDP(19) + (t332 + t360) * MDP(20) + (t242 + t310) * MDP(21) + t247 * MDP(22) + (t182 * t359 - t186 * t284 - t189 * t236 + t184) * MDP(23) + (t181 * t359 - t185 * t284 - t186 * t286 + t189 * t234) * MDP(24) + (-MDP(20) * t329 - MDP(21) * t236 - MDP(23) * t182 - MDP(24) * t181) * qJD(6);];
tauc  = t1;
