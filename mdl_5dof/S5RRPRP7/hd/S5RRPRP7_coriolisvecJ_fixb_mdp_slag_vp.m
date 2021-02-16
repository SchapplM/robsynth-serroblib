% Calculate Coriolis joint torque vector for
% S5RRPRP7
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
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:31
% EndTime: 2021-01-15 20:41:42
% DurationCPUTime: 3.00s
% Computational Cost: add. (2851->346), mult. (7215->451), div. (0->0), fcn. (4931->6), ass. (0->135)
t285 = cos(qJ(2));
t344 = cos(pkin(8));
t309 = t344 * t285;
t271 = qJD(1) * t309;
t281 = sin(pkin(8));
t283 = sin(qJ(2));
t321 = qJD(1) * t283;
t251 = t281 * t321 - t271;
t248 = qJD(4) + t251;
t316 = qJD(1) * qJD(2);
t353 = -0.2e1 * t316;
t352 = MDP(5) * (t283 ^ 2 - t285 ^ 2);
t263 = t281 * t285 + t344 * t283;
t277 = -pkin(2) * t285 - pkin(1);
t295 = -t281 * t283 + t309;
t225 = -pkin(3) * t295 - pkin(7) * t263 + t277;
t346 = -qJ(3) - pkin(6);
t268 = t346 * t283;
t269 = t346 * t285;
t233 = t281 * t268 - t344 * t269;
t282 = sin(qJ(4));
t284 = cos(qJ(4));
t325 = t282 * t225 + t284 * t233;
t266 = qJD(1) * t269;
t257 = t281 * t266;
t265 = qJD(1) * t268;
t345 = qJD(2) * pkin(2);
t260 = t265 + t345;
t226 = t344 * t260 + t257;
t221 = -qJD(2) * pkin(3) - t226;
t318 = t284 * qJD(2);
t323 = qJD(1) * t263;
t234 = t282 * t323 - t318;
t236 = qJD(2) * t282 + t284 * t323;
t192 = t234 * pkin(4) - t236 * qJ(5) + t221;
t253 = t263 * qJD(2);
t246 = qJD(1) * t253;
t274 = pkin(2) * t281 + pkin(7);
t333 = t274 * t246;
t351 = t192 * t248 - t333;
t312 = t283 * t316;
t270 = t281 * t312;
t290 = qJD(2) * t271 - t270;
t208 = t236 * qJD(4) + t282 * t290;
t350 = t236 ^ 2;
t349 = t248 ^ 2;
t348 = pkin(2) * t283;
t347 = pkin(4) * t246;
t343 = qJ(5) * t246;
t311 = qJD(2) * t346;
t249 = qJD(3) * t285 + t283 * t311;
t243 = t249 * qJD(1);
t289 = -qJD(3) * t283 + t285 * t311;
t244 = t289 * qJD(1);
t205 = t243 * t281 - t344 * t244;
t320 = qJD(4) * t282;
t207 = -qJD(4) * t318 - t284 * t290 + t320 * t323;
t183 = pkin(4) * t208 + qJ(5) * t207 - qJD(5) * t236 + t205;
t342 = t183 * t282;
t322 = qJD(1) * t277;
t267 = qJD(3) + t322;
t210 = pkin(3) * t251 - pkin(7) * t323 + t267;
t310 = t344 * t266;
t227 = t281 * t260 - t310;
t222 = qJD(2) * pkin(7) + t227;
t191 = t210 * t282 + t222 * t284;
t186 = qJ(5) * t248 + t191;
t341 = t186 * t248;
t340 = t191 * t248;
t339 = t207 * t282;
t338 = t234 * t251;
t337 = t234 * t323;
t336 = t236 * t234;
t307 = t236 * t248;
t335 = t236 * t323;
t334 = t263 * t284;
t240 = t282 * t246;
t332 = t282 * t248;
t286 = qJD(2) ^ 2;
t331 = t283 * t286;
t241 = t284 * t246;
t330 = t285 * t286;
t287 = qJD(1) ^ 2;
t329 = t285 * t287;
t228 = t265 * t281 - t310;
t302 = pkin(4) * t282 - qJ(5) * t284;
t328 = qJD(5) * t282 - t248 * t302 + t228;
t319 = qJD(4) * t284;
t327 = -t282 * t208 - t234 * t319;
t315 = pkin(2) * t321;
t217 = pkin(3) * t323 + pkin(7) * t251 + t315;
t229 = t344 * t265 + t257;
t326 = t282 * t217 + t284 * t229;
t190 = t210 * t284 - t222 * t282;
t317 = qJD(5) - t190;
t314 = t283 * t345;
t313 = t274 * t320;
t308 = pkin(1) * t353;
t215 = t249 * t281 - t344 * t289;
t232 = -t344 * t268 - t269 * t281;
t305 = 0.2e1 * t323;
t206 = t344 * t243 + t281 * t244;
t272 = pkin(2) * t312;
t209 = t246 * pkin(3) - pkin(7) * t290 + t272;
t304 = t282 * t206 - t284 * t209 + t210 * t320 + t222 * t319;
t275 = -t344 * pkin(2) - pkin(3);
t303 = t284 * pkin(4) + t282 * qJ(5);
t185 = -pkin(4) * t248 + t317;
t301 = t185 * t284 - t186 * t282;
t300 = t205 * t263 - t233 * t246;
t299 = t240 + (t251 * t284 + t319) * t248;
t298 = -t248 * t320 - t251 * t332 + t241;
t256 = t295 * qJD(2);
t297 = t256 * t282 + t263 * t319;
t296 = -t256 * t284 + t263 * t320;
t294 = t192 * t236 + t304;
t293 = t284 * t206 + t282 * t209 + t210 * t319 - t222 * t320;
t216 = t344 * t249 + t281 * t289;
t218 = pkin(3) * t253 - pkin(7) * t256 + t314;
t292 = t284 * t216 + t282 * t218 + t225 * t319 - t233 * t320;
t291 = t248 * t221 - t333;
t261 = -t303 + t275;
t197 = pkin(4) * t236 + qJ(5) * t234;
t196 = t302 * t263 + t232;
t194 = pkin(4) * t295 - t225 * t284 + t233 * t282;
t193 = -qJ(5) * t295 + t325;
t189 = t234 * t248 - t207;
t188 = -pkin(4) * t323 - t217 * t284 + t229 * t282;
t187 = qJ(5) * t323 + t326;
t184 = t302 * t256 + (t303 * qJD(4) - qJD(5) * t284) * t263 + t215;
t182 = -pkin(4) * t253 + t325 * qJD(4) + t216 * t282 - t218 * t284;
t181 = qJ(5) * t253 - qJD(5) * t295 + t292;
t180 = t304 - t347;
t179 = qJD(5) * t248 + t293 + t343;
t1 = [0.2e1 * t285 * MDP(4) * t312 + t352 * t353 + MDP(6) * t330 - MDP(7) * t331 + (-pkin(6) * t330 + t283 * t308) * MDP(9) + (pkin(6) * t331 + t285 * t308) * MDP(10) + (t246 * t277 + t253 * t267 + (-t215 + (-qJD(1) * t295 + t251) * t348) * qJD(2)) * MDP(11) + (t267 * t256 - t277 * t270 + (t277 * t271 + t305 * t348 - t216) * qJD(2)) * MDP(12) + (t206 * t295 + t215 * t323 - t216 * t251 - t226 * t256 - t227 * t253 + t232 * t290 + t300) * MDP(13) + (t205 * t232 + t206 * t233 - t215 * t226 + t216 * t227 + (t267 + t322) * t314) * MDP(14) + (-t207 * t334 - t236 * t296) * MDP(15) + ((-t234 * t284 - t236 * t282) * t256 + (t339 - t208 * t284 + (t234 * t282 - t236 * t284) * qJD(4)) * t263) * MDP(16) + (t207 * t295 + t236 * t253 + t263 * t241 - t248 * t296) * MDP(17) + (t208 * t295 - t234 * t253 - t263 * t240 - t248 * t297) * MDP(18) + (-t246 * t295 + t248 * t253) * MDP(19) + (t304 * t295 + t190 * t253 + t215 * t234 + t232 * t208 + ((-qJD(4) * t233 + t218) * t248 + t225 * t246 + t221 * qJD(4) * t263) * t284 + ((-qJD(4) * t225 - t216) * t248 + t221 * t256 + t300) * t282) * MDP(20) + (-t191 * t253 + t205 * t334 - t232 * t207 + t215 * t236 - t296 * t221 - t325 * t246 - t292 * t248 + t293 * t295) * MDP(21) + (t180 * t295 - t182 * t248 + t184 * t234 - t185 * t253 + t192 * t297 - t194 * t246 + t196 * t208 + t263 * t342) * MDP(22) + (-t181 * t234 + t182 * t236 - t193 * t208 - t194 * t207 + t301 * t256 + (-t179 * t282 + t180 * t284 + (-t185 * t282 - t186 * t284) * qJD(4)) * t263) * MDP(23) + (-t179 * t295 + t181 * t248 - t183 * t334 - t184 * t236 + t186 * t253 + t192 * t296 + t193 * t246 + t196 * t207) * MDP(24) + (t179 * t193 + t180 * t194 + t181 * t186 + t182 * t185 + t183 * t196 + t184 * t192) * MDP(25); -t283 * MDP(4) * t329 + t287 * t352 + (qJD(2) * t228 - t251 * t315 - t267 * t323 - t205) * MDP(11) + (t229 * qJD(2) + t267 * t251 - t315 * t323 - t206) * MDP(12) + ((t227 - t228) * t323 + (t229 - t226) * t251 + (-t281 * t246 - t290 * t344) * pkin(2)) * MDP(13) + (t226 * t228 - t227 * t229 + (-t205 * t344 + t206 * t281 - t267 * t321) * pkin(2)) * MDP(14) + (t284 * t307 - t339) * MDP(15) + ((-t207 - t338) * t284 - t236 * t332 + t327) * MDP(16) + (t299 - t335) * MDP(17) + (t298 + t337) * MDP(18) - t248 * t323 * MDP(19) + (-t190 * t323 + t275 * t208 - t228 * t234 + (-t205 + (-qJD(4) * t274 - t217) * t248) * t284 + (t229 * t248 + t291) * t282) * MDP(20) + (t191 * t323 + t205 * t282 - t275 * t207 - t228 * t236 + (t313 + t326) * t248 + t291 * t284) * MDP(21) + (-t183 * t284 + t185 * t323 + t208 * t261 + (-t274 * t319 + t188) * t248 - t328 * t234 + t351 * t282) * MDP(22) + (t187 * t234 - t188 * t236 + (t185 * t251 - t208 * t274 + t179 + (t236 * t274 + t185) * qJD(4)) * t284 + (-t186 * t251 - t207 * t274 + t180 + (t234 * t274 - t186) * qJD(4)) * t282) * MDP(23) + (-t342 - t186 * t323 + t207 * t261 + (-t187 - t313) * t248 + t328 * t236 - t351 * t284) * MDP(24) + (t183 * t261 - t185 * t188 - t186 * t187 - t328 * t192 + (qJD(4) * t301 + t179 * t284 + t180 * t282) * t274) * MDP(25) + (MDP(9) * t283 * t287 + MDP(10) * t329) * pkin(1); t305 * qJD(2) * MDP(11) + (-t270 + (t271 - t251) * qJD(2)) * MDP(12) + (-t251 ^ 2 - t323 ^ 2) * MDP(13) + (t226 * t323 + t227 * t251 + t272) * MDP(14) + (t298 - t337) * MDP(20) + (-t284 * t349 - t240 - t335) * MDP(21) + (-t248 * t332 + t241 - t337) * MDP(22) + ((t207 - t338) * t284 + t282 * t307 + t327) * MDP(23) + (t299 + t335) * MDP(24) + (-t192 * t323 + (-t180 + t341) * t284 + (t185 * t248 + t179) * t282) * MDP(25); MDP(15) * t336 + (-t234 ^ 2 + t350) * MDP(16) + t189 * MDP(17) + (-t208 + t307) * MDP(18) + t246 * MDP(19) + (-t221 * t236 - t304 + t340) * MDP(20) + (t190 * t248 + t221 * t234 - t293) * MDP(21) + (-t197 * t234 - t294 + t340 + 0.2e1 * t347) * MDP(22) + (pkin(4) * t207 - qJ(5) * t208 + (t186 - t191) * t236 + (t185 - t317) * t234) * MDP(23) + (0.2e1 * t343 - t192 * t234 + t197 * t236 + (0.2e1 * qJD(5) - t190) * t248 + t293) * MDP(24) + (-pkin(4) * t180 + qJ(5) * t179 - t185 * t191 + t186 * t317 - t192 * t197) * MDP(25); (-qJD(2) * t323 + t336) * MDP(22) + t189 * MDP(23) + (-t349 - t350) * MDP(24) + (t294 - t341 - t347) * MDP(25);];
tauc = t1;
