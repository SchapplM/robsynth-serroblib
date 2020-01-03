% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:18
% EndTime: 2019-12-31 17:28:22
% DurationCPUTime: 2.50s
% Computational Cost: add. (1265->287), mult. (3251->419), div. (0->0), fcn. (2168->6), ass. (0->137)
t267 = sin(qJ(3));
t268 = sin(qJ(2));
t319 = qJD(1) * t268;
t304 = t267 * t319;
t270 = cos(qJ(3));
t309 = t270 * qJD(2);
t233 = t304 - t309;
t269 = cos(qJ(4));
t317 = qJD(2) * t267;
t235 = t270 * t319 + t317;
t266 = sin(qJ(4));
t339 = t235 * t266;
t194 = t269 * t233 + t339;
t271 = cos(qJ(2));
t318 = qJD(1) * t271;
t256 = -qJD(3) + t318;
t252 = -qJD(4) + t256;
t357 = t194 * t252;
t281 = t233 * t266 - t269 * t235;
t356 = t252 * t281;
t262 = pkin(5) * t318;
t248 = qJD(2) * pkin(6) + t262;
t243 = -pkin(2) * t271 - pkin(6) * t268 - pkin(1);
t227 = t243 * qJD(1);
t335 = t267 * t227;
t201 = t248 * t270 + t335;
t190 = -pkin(7) * t233 + t201;
t312 = qJD(4) * t266;
t188 = t190 * t312;
t247 = -qJD(2) * pkin(2) + pkin(5) * t319;
t210 = pkin(3) * t233 + t247;
t355 = t194 * t210 + t188;
t308 = qJD(1) * qJD(2);
t297 = t271 * t308;
t314 = qJD(3) * t267;
t301 = t268 * t314;
t307 = qJD(2) * qJD(3);
t208 = -qJD(1) * t301 + (t297 + t307) * t270;
t285 = pkin(2) * t268 - pkin(6) * t271;
t241 = t285 * qJD(2);
t228 = qJD(1) * t241;
t298 = t268 * t308;
t287 = pkin(5) * t298;
t325 = -t270 * t228 - t267 * t287;
t275 = -t201 * qJD(3) - t325;
t179 = pkin(3) * t298 - pkin(7) * t208 + t275;
t313 = qJD(3) * t270;
t299 = t268 * t313;
t315 = qJD(2) * t271;
t303 = t267 * t315;
t276 = t299 + t303;
t209 = t276 * qJD(1) + t267 * t307;
t278 = t227 * t313 + t267 * t228 - t248 * t314;
t274 = -t270 * t287 + t278;
t182 = -pkin(7) * t209 + t274;
t292 = t269 * t179 - t266 * t182;
t354 = t210 * t281 + t292;
t295 = MDP(22) * t319;
t353 = qJD(2) * t295 + (-t194 ^ 2 + t281 ^ 2) * MDP(19) - t194 * t281 * MDP(18);
t352 = -0.2e1 * t308;
t237 = t266 * t270 + t267 * t269;
t218 = t237 * t268;
t351 = t268 * MDP(4);
t264 = t268 ^ 2;
t350 = (-t271 ^ 2 + t264) * MDP(5);
t349 = t271 * t309 - t301;
t348 = qJD(3) + qJD(4);
t291 = t208 * t266 + t269 * t209;
t181 = -t281 * qJD(4) + t291;
t347 = pkin(6) + pkin(7);
t346 = pkin(3) * t267;
t345 = pkin(5) * t267;
t200 = t270 * t227 - t248 * t267;
t189 = -pkin(7) * t235 + t200;
t187 = -pkin(3) * t256 + t189;
t344 = t187 * t269;
t343 = t190 * t269;
t342 = t208 * t267;
t341 = t233 * t256;
t340 = t235 * t256;
t338 = t247 * t267;
t337 = t247 * t270;
t336 = t256 * t270;
t334 = t267 * t268;
t333 = t267 * t271;
t332 = t268 * t270;
t272 = qJD(2) ^ 2;
t331 = t268 * t272;
t330 = t270 * t271;
t329 = t271 * t272;
t273 = qJD(1) ^ 2;
t328 = t271 * t273;
t236 = t266 * t267 - t269 * t270;
t277 = t236 * t271;
t327 = qJD(1) * t277 - t348 * t236;
t326 = (-t318 + t348) * t237;
t324 = t267 * t241 + t243 * t313;
t316 = qJD(2) * t268;
t323 = t270 * t241 + t316 * t345;
t238 = t285 * qJD(1);
t322 = pkin(5) * t304 + t270 * t238;
t257 = pkin(5) * t330;
t321 = t267 * t243 + t257;
t311 = qJD(4) * t269;
t306 = t269 * t208 - t266 * t209 - t233 * t311;
t305 = qJD(3) * t347;
t300 = t271 * t314;
t294 = MDP(15) * t316;
t293 = pkin(1) * t352;
t290 = t233 + t309;
t289 = -t235 + t317;
t288 = qJD(4) * t187 + t182;
t286 = pkin(3) * t314 - t318 * t346 - t262;
t250 = t347 * t270;
t279 = pkin(3) * t268 - pkin(7) * t330;
t284 = t279 * qJD(1) + qJD(4) * t250 + t270 * t305 + t322;
t223 = t267 * t238;
t249 = t347 * t267;
t283 = -qJD(4) * t249 - t223 - (-pkin(5) * t332 - pkin(7) * t333) * qJD(1) - t267 * t305;
t177 = t187 * t266 + t343;
t232 = t270 * t243;
t199 = -pkin(7) * t332 + t232 + (-pkin(3) - t345) * t271;
t204 = -pkin(7) * t334 + t321;
t282 = t199 * t266 + t204 * t269;
t280 = qJD(1) * t264 - t256 * t271;
t180 = -t235 * t312 + t306;
t260 = -pkin(3) * t270 - pkin(2);
t242 = (pkin(5) + t346) * t268;
t219 = t236 * t268;
t211 = t276 * pkin(3) + pkin(5) * t315;
t192 = pkin(3) * t209 + pkin(5) * t297;
t186 = -t312 * t334 + (t348 * t332 + t303) * t269 + t349 * t266;
t185 = -qJD(2) * t277 - t348 * t218;
t184 = -t276 * pkin(7) + (-t268 * t309 - t300) * pkin(5) + t324;
t183 = t279 * qJD(2) + (-t257 + (pkin(7) * t268 - t243) * t267) * qJD(3) + t323;
t176 = -t190 * t266 + t344;
t1 = [0.2e1 * t297 * t351 + t350 * t352 + MDP(6) * t329 - MDP(7) * t331 + (-pkin(5) * t329 + t268 * t293) * MDP(9) + (pkin(5) * t331 + t271 * t293) * MDP(10) + (t208 * t332 + t349 * t235) * MDP(11) + ((-t233 * t270 - t235 * t267) * t315 + (-t342 - t209 * t270 + (t233 * t267 - t235 * t270) * qJD(3)) * t268) * MDP(12) + (t256 * t301 - t208 * t271 + (t235 * t268 + t280 * t270) * qJD(2)) * MDP(13) + (t256 * t299 + t209 * t271 + (-t233 * t268 - t280 * t267) * qJD(2)) * MDP(14) + (-t256 - t318) * t294 + (-(-t243 * t314 + t323) * t256 + (t247 * t313 + pkin(5) * t209 + (t232 * qJD(1) + t200) * qJD(2)) * t268 + ((pkin(5) * t233 + t338) * qJD(2) + (t335 + (pkin(5) * t256 + t248) * t270) * qJD(3) + t325) * t271) * MDP(16) + ((-pkin(5) * t300 + t324) * t256 + t278 * t271 + (pkin(5) * t208 - t247 * t314) * t268 + ((pkin(5) * t235 + t337) * t271 + (-pkin(5) * t336 - t321 * qJD(1) - t201) * t268) * qJD(2)) * MDP(17) + (-t180 * t219 - t185 * t281) * MDP(18) + (-t180 * t218 + t181 * t219 - t185 * t194 + t186 * t281) * MDP(19) + (-t180 * t271 - t185 * t252 + (-qJD(1) * t219 - t281) * t316) * MDP(20) + (t181 * t271 + t186 * t252 + (-qJD(1) * t218 - t194) * t316) * MDP(21) + (-t252 - t318) * MDP(22) * t316 + (-(t183 * t269 - t184 * t266) * t252 - t292 * t271 + t211 * t194 + t242 * t181 + t192 * t218 + t210 * t186 + (t177 * t271 + t282 * t252) * qJD(4) + ((t199 * t269 - t204 * t266) * qJD(1) + t176) * t316) * MDP(23) + (t242 * t180 + t210 * t185 - t188 * t271 - t192 * t219 - t211 * t281 + ((-qJD(4) * t204 + t183) * t252 + t179 * t271) * t266 + ((qJD(4) * t199 + t184) * t252 + t288 * t271) * t269 + (-t282 * qJD(1) - t177) * t316) * MDP(24); -t328 * t351 + t273 * t350 + (-t235 * t336 + t342) * MDP(11) + ((t208 + t341) * t270 + (-t209 + t340) * t267) * MDP(12) + (-t256 * t313 + (t256 * t330 + t289 * t268) * qJD(1)) * MDP(13) + (t256 * t314 + (-t256 * t333 + t290 * t268) * qJD(1)) * MDP(14) + t256 * MDP(15) * t319 + (-pkin(2) * t209 + t322 * t256 + (pkin(6) * t336 + t338) * qJD(3) + ((-pkin(6) * t317 - t200) * t268 + (-t290 * pkin(5) - t338) * t271) * qJD(1)) * MDP(16) + (-pkin(2) * t208 - t223 * t256 + (-pkin(6) * t256 * t267 + t337) * qJD(3) + (-t247 * t330 + (-pkin(6) * t309 + t201) * t268 + (t256 * t332 + t289 * t271) * pkin(5)) * qJD(1)) * MDP(17) + (t180 * t237 - t281 * t327) * MDP(18) + (-t180 * t236 - t181 * t237 - t327 * t194 + t281 * t326) * MDP(19) + (-t327 * t252 + (qJD(2) * t237 + t281) * t319) * MDP(20) + (t326 * t252 + (-qJD(2) * t236 + t194) * t319) * MDP(21) + t252 * t295 + (t260 * t181 + t192 * t236 + (t283 * t266 + t284 * t269) * t252 + t326 * t210 + t286 * t194 + ((-t249 * t269 - t250 * t266) * qJD(2) - t176) * t319) * MDP(23) + (t260 * t180 + t192 * t237 + (-t284 * t266 + t283 * t269) * t252 + t327 * t210 - t286 * t281 + (-(-t249 * t266 + t250 * t269) * qJD(2) + t177) * t319) * MDP(24) + (t273 * t268 * MDP(9) + MDP(10) * t328) * pkin(1); t235 * t233 * MDP(11) + (-t233 ^ 2 + t235 ^ 2) * MDP(12) + (t208 - t341) * MDP(13) + (-t209 - t340) * MDP(14) + qJD(1) * t294 + (-t201 * t256 - t235 * t247 + t275) * MDP(16) + (-t200 * t256 + t233 * t247 - t274) * MDP(17) + (t180 - t357) * MDP(20) + (-t181 + t356) * MDP(21) + ((-t189 * t266 - t343) * t252 - t177 * qJD(4) + (-t194 * t235 + t252 * t312 + t269 * t298) * pkin(3) + t354) * MDP(23) + ((t190 * t252 - t179) * t266 + (-t189 * t252 - t288) * t269 + (t235 * t281 + t252 * t311 - t266 * t298) * pkin(3) + t355) * MDP(24) + t353; (t306 - t357) * MDP(20) + (-t291 + t356) * MDP(21) + (-t177 * t252 + t354) * MDP(23) + (-t176 * t252 - t266 * t179 - t269 * t182 + t355) * MDP(24) + (-MDP(20) * t339 + t281 * MDP(21) - t177 * MDP(23) - MDP(24) * t344) * qJD(4) + t353;];
tauc = t1;
