% Calculate Coriolis joint torque vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:53
% EndTime: 2021-01-15 22:14:58
% DurationCPUTime: 1.52s
% Computational Cost: add. (2066->240), mult. (3588->304), div. (0->0), fcn. (2144->6), ass. (0->132)
t359 = -qJ(4) - pkin(7);
t293 = cos(qJ(3));
t283 = t293 * qJD(4);
t291 = sin(qJ(3));
t317 = qJD(3) * t359;
t252 = t291 * t317 + t283;
t290 = cos(pkin(8));
t350 = t290 * t293;
t289 = sin(pkin(8));
t351 = t289 * t291;
t259 = -t350 + t351;
t300 = -qJD(4) * t291 + t293 * t317;
t294 = cos(qJ(2));
t358 = pkin(1) * qJD(1);
t326 = t294 * t358;
t364 = t290 * t252 + t259 * t326 + t289 * t300;
t368 = MDP(16) + MDP(19);
t367 = MDP(7) * t291;
t366 = MDP(8) * (t291 ^ 2 - t293 ^ 2);
t260 = t289 * t293 + t290 * t291;
t286 = qJD(1) + qJD(2);
t249 = t260 * t286;
t332 = t293 * MDP(12);
t365 = t286 * (-t291 * MDP(13) + MDP(5) + t332);
t344 = t252 * t289 - t260 * t326 - t290 * t300;
t331 = t293 * MDP(13);
t362 = t291 * MDP(12) + t331;
t328 = MDP(14) + MDP(18);
t327 = MDP(15) - MDP(20);
t244 = t249 ^ 2;
t361 = pkin(1) * t294;
t360 = pkin(2) * t286;
t357 = pkin(1) * qJD(2);
t324 = qJD(1) * t357;
t313 = t294 * t324;
t301 = qJD(4) * t286 + t313;
t292 = sin(qJ(2));
t315 = -t359 * t286 + t292 * t358;
t309 = qJD(3) * t315;
t217 = -t291 * t309 + t293 * t301;
t296 = -t291 * t301 - t293 * t309;
t188 = t217 * t289 - t290 * t296;
t278 = pkin(1) * t292 + pkin(7);
t284 = t293 * qJ(4);
t258 = t278 * t293 + t284;
t347 = -qJ(4) - t278;
t316 = t347 * t291;
t223 = t258 * t289 - t290 * t316;
t356 = t188 * t223;
t270 = pkin(7) * t293 + t284;
t318 = t359 * t291;
t233 = t270 * t289 - t290 * t318;
t355 = t188 * t233;
t354 = t188 * t260;
t237 = t315 * t293;
t353 = t237 * t289;
t352 = t286 * t291;
t231 = t290 * t237;
t295 = qJD(3) ^ 2;
t349 = t291 * t295;
t348 = t293 * t295;
t254 = t260 * qJD(3);
t240 = t286 * t254;
t336 = qJD(3) * t291;
t319 = t289 * t336;
t266 = t286 * t319;
t335 = qJD(3) * t293;
t320 = t290 * t335;
t241 = t286 * t320 - t266;
t273 = t292 * t324;
t322 = t286 * t336;
t253 = pkin(3) * t322 + t273;
t305 = pkin(4) * t240 - qJ(5) * t241 + t253;
t191 = -qJD(5) * t249 + t305;
t280 = -pkin(3) * t293 - pkin(2);
t246 = t280 * t286 + qJD(4) - t326;
t323 = t286 * t350;
t247 = t286 * t351 - t323;
t200 = pkin(4) * t247 - qJ(5) * t249 + t246;
t346 = t191 * t259 + t200 * t254;
t255 = -t319 + t320;
t345 = -t191 * t260 - t200 * t255;
t189 = t290 * t217 + t289 * t296;
t342 = t246 * t254 + t253 * t259;
t341 = t246 * t255 + t253 * t260;
t236 = t315 * t291;
t235 = qJD(3) * pkin(3) - t236;
t207 = t289 * t235 + t231;
t265 = -t326 - t360;
t340 = t265 * t335 + t291 * t273;
t314 = qJD(3) * t347;
t325 = t294 * t357;
t229 = t291 * t314 + t293 * t325 + t283;
t297 = (-qJD(4) - t325) * t291 + t293 * t314;
t198 = t229 * t289 - t290 * t297;
t338 = qJD(3) * t198;
t199 = t290 * t229 + t289 * t297;
t337 = qJD(3) * t199;
t334 = t246 * MDP(17);
t210 = -t236 * t290 - t353;
t330 = qJD(5) - t210;
t329 = qJD(3) * qJD(5);
t321 = t286 * t335;
t187 = t329 + t189;
t206 = t235 * t290 - t353;
t203 = -qJD(3) * pkin(4) + qJD(5) - t206;
t204 = qJD(3) * qJ(5) + t207;
t312 = -t187 * t259 + t203 * t255 - t204 * t254 + t354;
t311 = -t189 * t259 - t206 * t255 - t207 * t254 + t354;
t310 = qJD(3) * t210 - t189;
t308 = t326 - t360;
t224 = t290 * t258 + t289 * t316;
t303 = t198 * t249 - t199 * t247 + t223 * t241 - t224 * t240;
t234 = t290 * t270 + t289 * t318;
t302 = t233 * t241 - t234 * t240 - t364 * t247;
t225 = pkin(4) * t259 - qJ(5) * t260 + t280;
t299 = t247 * MDP(14) + t249 * MDP(15) + t334;
t281 = pkin(3) * t336;
t208 = pkin(4) * t254 - qJ(5) * t255 - qJD(5) * t260 + t281;
t298 = -0.2e1 * t286 * qJD(3) * t366 - MDP(10) * t349 - t273 * MDP(5) + MDP(9) * t348 + 0.2e1 * t321 * t367;
t282 = t292 * t357;
t279 = -pkin(2) - t361;
t276 = -pkin(3) * t290 - pkin(4);
t274 = pkin(3) * t289 + qJ(5);
t267 = t280 - t361;
t263 = t282 + t281;
t256 = t265 * t336;
t220 = t225 - t361;
t211 = pkin(3) * t352 + pkin(4) * t249 + qJ(5) * t247;
t209 = -t236 * t289 + t231;
t201 = t208 + t282;
t1 = [(-t278 * t348 + t279 * t322 + t256) * MDP(12) + (t278 * t349 + t279 * t321 + t340) * MDP(13) + (t240 * t267 + t247 * t263 - t338 + t342) * MDP(14) + (t241 * t267 + t249 * t263 - t337 + t341) * MDP(15) + (t303 + t311) * MDP(16) + (t189 * t224 - t198 * t206 + t199 * t207 + t246 * t263 + t253 * t267 + t356) * MDP(17) + (t201 * t247 + t220 * t240 - t338 + t346) * MDP(18) + (t303 + t312) * MDP(19) + (-t201 * t249 - t220 * t241 + t337 + t345) * MDP(20) + (t187 * t224 + t191 * t220 + t198 * t203 + t199 * t204 + t200 * t201 + t356) * MDP(21) + (((-qJD(1) - t286) * MDP(6) - t362 * qJD(3)) * t294 + (-qJD(1) * t332 - t365) * t292) * t357 + t298; (-pkin(7) * t348 + t256) * MDP(12) + (pkin(7) * t349 + t340) * MDP(13) + (t240 * t280 + t342) * MDP(14) + (t241 * t280 + t341) * MDP(15) + (t302 + t311) * MDP(16) + (t189 * t234 - t344 * t206 + t364 * t207 + t253 * t280 + t355) * MDP(17) + (t208 * t247 + t225 * t240 + t346) * MDP(18) + (t302 + t312) * MDP(19) + (-t225 * t241 + t345) * MDP(20) + (t187 * t234 + t191 * t225 + t200 * t208 + t344 * t203 + t364 * t204 + t355) * MDP(21) + (-t208 * MDP(20) + t368 * t344) * t249 + ((-qJD(2) + t286) * MDP(6) * t294 + (-t200 * MDP(21) - qJD(2) * t332 - t247 * t328 - t249 * t327 - t334 + t365) * t292) * t358 + (t308 * t331 + (MDP(12) * t308 + pkin(3) * t299) * t291 - t328 * t344 - t327 * t364) * qJD(3) + t298; t310 * MDP(15) + (t206 * t209 - t207 * t210) * MDP(17) + (-t240 * t274 + t241 * t276) * MDP(19) + (-t310 + 0.2e1 * t329) * MDP(20) + (t187 * t274 + t188 * t276 - t200 * t211 - t203 * t209 + t204 * t330) * MDP(21) + (-t293 * t367 + t366) * t286 ^ 2 + (-t246 * MDP(14) + (t207 - t209) * MDP(16) - t200 * MDP(18) + (t204 - t209) * MDP(19) + t211 * MDP(20)) * t249 + (t246 * MDP(15) + (-t206 + t210) * MDP(16) - t211 * MDP(18) + (t203 - t330) * MDP(19) - t200 * MDP(20)) * t247 + ((-t240 * t289 - t241 * t290) * MDP(16) + (-t188 * t290 + t189 * t289) * MDP(17) - t299 * t352) * pkin(3) + t328 * (qJD(3) * t209 - t188) + t362 * (-t265 * t286 - t313); (t206 * t249 + t207 * t247 + t253) * MDP(17) + (t204 * t247 + (-qJD(5) - t203) * t249 + t305) * MDP(21) + t327 * (-t266 + (-t247 + t323) * qJD(3)) + 0.2e1 * t328 * t249 * qJD(3) + t368 * (-t247 ^ 2 - t244); t249 * t247 * MDP(18) + (-t266 + (t247 + t323) * qJD(3)) * MDP(19) + (-t244 - t295) * MDP(20) + (-qJD(3) * t204 + t200 * t249 + t188) * MDP(21);];
tauc = t1;
