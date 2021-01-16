% Calculate Coriolis joint torque vector for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:54
% EndTime: 2021-01-15 21:03:01
% DurationCPUTime: 2.53s
% Computational Cost: add. (1653->318), mult. (3725->427), div. (0->0), fcn. (2022->4), ass. (0->148)
t352 = pkin(3) + pkin(6);
t277 = cos(qJ(2));
t278 = -pkin(2) - pkin(7);
t275 = sin(qJ(2));
t303 = -qJ(3) * t275 - pkin(1);
t231 = t278 * t277 + t303;
t211 = t231 * qJD(1);
t327 = qJD(1) * t275;
t263 = pkin(6) * t327;
t354 = qJD(3) + t263;
t315 = pkin(3) * t327 + t354;
t214 = t278 * qJD(2) + t315;
t274 = sin(qJ(4));
t276 = cos(qJ(4));
t193 = t211 * t276 + t214 * t274;
t326 = qJD(1) * t277;
t307 = t274 * t326;
t323 = qJD(2) * t276;
t236 = -t307 + t323;
t314 = qJD(1) * qJD(2);
t306 = t275 * t314;
t253 = t274 * t306;
t325 = qJD(2) * t274;
t234 = t276 * t326 + t325;
t320 = qJD(4) * t234;
t204 = -t253 + t320;
t260 = pkin(2) * t306;
t294 = pkin(7) * t275 - qJ(3) * t277;
t321 = qJD(3) * t275;
t282 = t294 * qJD(2) - t321;
t201 = t282 * qJD(1) + t260;
t305 = t277 * t314;
t259 = pkin(6) * t305;
t230 = pkin(3) * t305 + t259;
t298 = -t201 * t274 + t276 * t230;
t289 = qJ(5) * t204 + t298;
t296 = pkin(4) * t305;
t182 = -t193 * qJD(4) - qJD(5) * t236 + t289 + t296;
t188 = -qJ(5) * t234 + t193;
t261 = qJD(4) + t327;
t350 = t188 * t261;
t362 = t182 + t350;
t318 = qJD(4) * t276;
t330 = qJD(4) * t307 + t276 * t306;
t205 = qJD(2) * t318 - t330;
t319 = qJD(4) * t274;
t295 = -t276 * t201 + t211 * t319 - t214 * t318 - t274 * t230;
t284 = -qJ(5) * t205 - t295;
t183 = -qJD(5) * t234 + t284;
t192 = -t211 * t274 + t276 * t214;
t187 = -qJ(5) * t236 + t192;
t186 = pkin(4) * t261 + t187;
t361 = -t261 * t186 + t183;
t360 = -0.2e1 * t314;
t345 = t236 * t261;
t272 = t275 ^ 2;
t273 = t277 ^ 2;
t358 = (t272 - t273) * MDP(5);
t357 = t186 - t187;
t356 = t234 * t261 + t204;
t355 = t205 + t345;
t251 = t352 * t275;
t331 = t276 * t231 + t274 * t251;
t353 = MDP(20) + MDP(22);
t351 = qJD(2) * pkin(2);
t349 = t204 * t276;
t324 = qJD(2) * t275;
t242 = t352 * t324;
t270 = qJD(2) * qJD(3);
t217 = -qJD(1) * t242 + t270;
t348 = t217 * t274;
t347 = t217 * t276;
t344 = t236 * t277;
t343 = t261 * t278;
t342 = t274 * t275;
t341 = t275 * t276;
t279 = qJD(2) ^ 2;
t340 = t275 * t279;
t339 = t276 * t277;
t338 = t277 * t279;
t280 = qJD(1) ^ 2;
t337 = t277 * t280;
t336 = qJ(5) - t278;
t267 = pkin(2) * t327;
t220 = t294 * qJD(1) + t267;
t264 = pkin(6) * t326;
t243 = pkin(3) * t326 + t264;
t297 = -t220 * t274 + t276 * t243;
t335 = (pkin(4) * t277 - qJ(5) * t342) * qJD(1) + t297 + qJD(5) * t276 - t336 * t319;
t247 = t336 * t276;
t333 = t276 * t220 + t274 * t243;
t334 = qJ(5) * t276 * t327 + qJD(4) * t247 + qJD(5) * t274 + t333;
t311 = -pkin(4) * t276 - pkin(3);
t332 = pkin(4) * t318 - t311 * t327 + t354;
t252 = t352 * t277;
t271 = qJD(2) * qJ(3);
t224 = t271 + t243;
t302 = -pkin(4) * t234 - qJD(5);
t199 = -t302 + t224;
t328 = MDP(25) * t199;
t249 = -pkin(2) * t277 + t303;
t225 = qJD(1) * t249;
t322 = qJD(2) * t277;
t317 = qJD(4) * t277;
t316 = qJD(5) * t277;
t313 = t261 * t341;
t312 = t275 * t337;
t310 = t274 * t317;
t309 = t261 * t318;
t308 = t276 * t317;
t304 = MDP(19) * t322;
t301 = pkin(1) * t360;
t300 = qJD(3) - t351;
t299 = qJ(5) * t277 - t231;
t293 = -qJD(1) * t273 + t261 * t275;
t292 = -0.2e1 * qJD(2) * t225;
t291 = t261 * t274;
t285 = -qJ(3) * t322 - t321;
t209 = t285 * qJD(1) + t260;
t266 = pkin(2) * t324;
t222 = t266 + t285;
t290 = pkin(6) * t279 + qJD(1) * t222 + t209;
t288 = t224 * t275 + t278 * t322;
t195 = pkin(4) * t205 + t217;
t287 = -t195 * t274 - t199 * t318;
t286 = t195 * t276 - t199 * t319;
t207 = t266 + t282;
t244 = t352 * t322;
t283 = t276 * t207 - t231 * t319 + t274 * t244 + t251 * t318;
t245 = pkin(6) * t306 - t270;
t248 = t263 + t300;
t250 = -t264 - t271;
t281 = -t245 * t277 + (t248 * t277 + (t250 + t264) * t275) * qJD(2);
t262 = pkin(4) * t274 + qJ(3);
t254 = t276 * t305;
t246 = t336 * t274;
t240 = -qJ(3) * t326 + t267;
t239 = t276 * t251;
t233 = t234 ^ 2;
t229 = t276 * t244;
t223 = pkin(4) * t339 + t252;
t213 = t225 * t327;
t200 = -pkin(4) * t310 + (-pkin(6) + t311) * t324;
t197 = -qJ(5) * t339 + t331;
t196 = pkin(4) * t275 + t299 * t274 + t239;
t185 = -t276 * t316 + (t275 * t323 + t310) * qJ(5) + t283;
t184 = pkin(4) * t322 + t229 + t299 * t318 + (-qJ(5) * t324 - qJD(4) * t251 - t207 + t316) * t274;
t1 = [0.2e1 * t275 * MDP(4) * t305 + t358 * t360 + MDP(6) * t338 - MDP(7) * t340 + (-pkin(6) * t338 + t275 * t301) * MDP(9) + (pkin(6) * t340 + t277 * t301) * MDP(10) + t281 * MDP(11) + (t275 * t292 + t290 * t277) * MDP(12) + (-t290 * t275 + t277 * t292) * MDP(13) + (t281 * pkin(6) + t209 * t249 + t222 * t225) * MDP(14) + (t204 * t274 * t277 + (t274 * t324 - t308) * t236) * MDP(15) + ((-t234 * t274 + t236 * t276) * t324 + (t349 + t205 * t274 + (t234 * t276 + t236 * t274) * qJD(4)) * t277) * MDP(16) + (-t261 * t308 - t204 * t275 + (t293 * t274 + t344) * qJD(2)) * MDP(17) + (t261 * t310 - t205 * t275 + (-t234 * t277 + t293 * t276) * qJD(2)) * MDP(18) + (t261 + t327) * t304 + ((-t207 * t274 + t229) * t261 - t242 * t234 + t252 * t205 + (-t224 * t323 + t298) * t275 + (-t193 * t275 - t261 * t331) * qJD(4) + (-t224 * t319 + t347 + ((-t231 * t274 + t239) * qJD(1) + t192) * qJD(2)) * t277) * MDP(20) + (-t283 * t261 - t242 * t236 - t252 * t204 + (t224 * t325 + t295) * t275 + (-t224 * t318 - t348 + (-t331 * qJD(1) - t193) * qJD(2)) * t277) * MDP(21) + (t184 * t261 + t200 * t234 + t205 * t223 + (-t199 * t323 + t182) * t275 + ((qJD(1) * t196 + t186) * qJD(2) + t286) * t277) * MDP(22) + (-t185 * t261 + t200 * t236 - t204 * t223 + (t199 * t325 - t183) * t275 + ((-qJD(1) * t197 - t188) * qJD(2) + t287) * t277) * MDP(23) + (-t184 * t236 - t185 * t234 + t196 * t204 - t197 * t205 + (-t186 * t274 + t188 * t276) * t324 + (t182 * t274 - t183 * t276 + (t186 * t276 + t188 * t274) * qJD(4)) * t277) * MDP(24) + (t182 * t196 + t183 * t197 + t184 * t186 + t185 * t188 + t195 * t223 + t199 * t200) * MDP(25); -MDP(4) * t312 + t280 * t358 + ((-t250 - t271) * t275 + (-t248 + t300) * t277) * qJD(1) * MDP(11) + (-t240 * t326 + t213) * MDP(12) + (0.2e1 * t270 + (t225 * t277 + t240 * t275) * qJD(1)) * MDP(13) + (-qJ(3) * t245 - qJD(3) * t250 - t225 * t240 + (-t250 * t275 + (-t248 - t351) * t277) * qJD(1) * pkin(6)) * MDP(14) + (-t236 * t291 - t349) * MDP(15) + (t274 * t356 - t276 * t355) * MDP(16) + (-t261 * t319 + t254 + (-t261 * t342 - t344) * qJD(1)) * MDP(17) + (-t309 + (-t313 + (t234 - t325) * t277) * qJD(1)) * MDP(18) - t261 * MDP(19) * t326 + (qJ(3) * t205 + t348 - t297 * t261 + t315 * t234 + (t224 * t276 - t274 * t343) * qJD(4) + (-t192 * t277 + t288 * t276) * qJD(1)) * MDP(20) + (-qJ(3) * t204 + t347 + t333 * t261 + t315 * t236 + (-t224 * t274 - t276 * t343) * qJD(4) + (t193 * t277 - t288 * t274) * qJD(1)) * MDP(21) + (t205 * t262 - t335 * t261 + t332 * t234 + (t199 * t341 + (-qJD(2) * t247 - t186) * t277) * qJD(1) - t287) * MDP(22) + (-t204 * t262 + t334 * t261 + t332 * t236 + (-t199 * t342 + (qJD(2) * t246 + t188) * t277) * qJD(1) + t286) * MDP(23) + (-t204 * t247 + t205 * t246 + t334 * t234 + t335 * t236 - t361 * t274 - t362 * t276) * MDP(24) + (-t182 * t247 - t183 * t246 - t335 * t186 - t334 * t188 + t195 * t262 + t332 * t199) * MDP(25) + (t280 * t275 * MDP(9) + MDP(10) * t337) * pkin(1); MDP(12) * t312 + (-t272 * t280 - t279) * MDP(13) + (qJD(2) * t250 + t213 + t259) * MDP(14) - qJD(2) * t328 + (MDP(21) + MDP(23)) * (-t309 - qJD(2) * t236 + (-t274 * t322 - t313) * qJD(1)) + t353 * (-qJD(2) * t234 - t261 * t291 + t254) + ((-t234 * t327 + t204 - t320) * MDP(24) + t362 * MDP(25)) * t276 + ((-t205 + t345) * MDP(24) + t361 * MDP(25)) * t274; -t233 * MDP(16) + t253 * MDP(17) + t330 * MDP(18) + qJD(1) * t304 + (t193 * t261 + t298) * MDP(20) + (t192 * t261 + t295) * MDP(21) + (t289 + 0.2e1 * t296 + t350) * MDP(22) + (t187 * t261 - t284) * MDP(23) + t204 * pkin(4) * MDP(24) + (pkin(4) * t182 + t357 * t188) * MDP(25) + (t261 * MDP(17) + t224 * MDP(21) + (qJD(5) + t199) * MDP(23) - t357 * MDP(24)) * t234 + (-t234 * MDP(17) - MDP(18) * t323 - t193 * t353) * qJD(4) + (t234 * MDP(15) + t261 * MDP(18) - t224 * MDP(20) + (-t199 + t302) * MDP(22) - pkin(4) * t328 + (-MDP(23) * pkin(4) + MDP(16)) * t236) * t236; t355 * MDP(22) - t356 * MDP(23) + (-t236 ^ 2 - t233) * MDP(24) + (t186 * t236 + t188 * t234 + t195) * MDP(25);];
tauc = t1;
