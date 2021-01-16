% Calculate Coriolis joint torque vector for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:31
% EndTime: 2021-01-15 16:45:40
% DurationCPUTime: 2.29s
% Computational Cost: add. (1626->315), mult. (4139->441), div. (0->0), fcn. (2849->8), ass. (0->139)
t271 = sin(qJ(3));
t364 = MDP(5) * t271;
t266 = t271 ^ 2;
t274 = cos(qJ(3));
t363 = (-t274 ^ 2 + t266) * MDP(6);
t272 = sin(qJ(2));
t268 = sin(pkin(5));
t320 = qJD(1) * t268;
t306 = t272 * t320;
t252 = qJD(2) * pkin(7) + t306;
t269 = cos(pkin(5));
t338 = t269 * t271;
t260 = qJD(1) * t338;
t226 = t274 * t252 + t260;
t218 = qJD(3) * pkin(8) + t226;
t254 = -pkin(3) * t274 - pkin(8) * t271 - pkin(2);
t275 = cos(qJ(2));
t305 = t275 * t320;
t228 = t254 * qJD(2) - t305;
t270 = sin(qJ(4));
t273 = cos(qJ(4));
t198 = -t218 * t270 + t273 * t228;
t315 = qJD(3) * t270;
t317 = qJD(2) * t271;
t246 = t273 * t317 + t315;
t194 = -qJ(5) * t246 + t198;
t316 = qJD(2) * t274;
t261 = -qJD(4) + t316;
t189 = -pkin(4) * t261 + t194;
t362 = t189 - t194;
t312 = qJD(4) * t270;
t298 = t271 * t312;
t308 = qJD(2) * qJD(3);
t295 = t274 * t308;
t307 = qJD(3) * qJD(4);
t323 = (t295 + t307) * t273;
t220 = qJD(2) * t298 - t323;
t301 = t270 * t317;
t309 = t273 * qJD(3);
t244 = t301 - t309;
t361 = t244 * t261 - t220;
t311 = qJD(4) * t273;
t297 = t271 * t311;
t313 = qJD(3) * t274;
t279 = t270 * t313 + t297;
t221 = t279 * qJD(2) + t270 * t307;
t360 = -t246 * t261 + t221;
t286 = pkin(3) * t271 - pkin(8) * t274;
t249 = t286 * qJD(3);
t332 = t274 * t275;
t359 = -(t270 * t272 + t273 * t332) * t320 + t270 * t249 + t254 * t311;
t319 = qJD(1) * t274;
t225 = -t271 * t252 + t269 * t319;
t314 = qJD(3) * t271;
t355 = pkin(7) * t270;
t358 = t273 * t249 + t314 * t355 - (-t270 * t332 + t272 * t273) * t320;
t357 = MDP(17) + MDP(19);
t356 = pkin(4) * t244;
t354 = -qJ(5) - pkin(8);
t353 = qJD(2) * pkin(2);
t352 = qJ(5) * t271;
t337 = t270 * t228;
t199 = t218 * t273 + t337;
t195 = -qJ(5) * t244 + t199;
t351 = t195 * t261;
t318 = qJD(2) * t268;
t302 = t275 * t318;
t289 = t271 * t302;
t204 = qJD(1) * t289 + qJD(3) * t260 + t252 * t313;
t196 = pkin(4) * t221 + t204;
t350 = t196 * t270;
t349 = t196 * t273;
t348 = t204 * t270;
t347 = t204 * t273;
t217 = -qJD(3) * pkin(3) - t225;
t346 = t217 * t270;
t345 = t220 * t270;
t342 = t246 * t270;
t341 = t261 * t273;
t340 = t268 * t272;
t339 = t268 * t275;
t336 = t270 * t274;
t335 = t271 * t273;
t276 = qJD(3) ^ 2;
t334 = t271 * t276;
t333 = t273 * t274;
t331 = t274 * t276;
t262 = pkin(7) * t333;
t283 = pkin(4) * t271 - qJ(5) * t333;
t310 = qJD(5) * t273;
t330 = -t271 * t310 + t283 * qJD(3) + (-t262 + (-t254 + t352) * t270) * qJD(4) + t358;
t329 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t335 + (-qJD(5) * t271 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t274) * t270 + t359;
t248 = t286 * qJD(2);
t291 = -t225 * t270 + t273 * t248;
t294 = qJD(4) * t354;
t328 = t283 * qJD(2) + qJD(5) * t270 - t273 * t294 + t291;
t300 = t270 * t316;
t326 = t273 * t225 + t270 * t248;
t327 = -qJ(5) * t300 - t270 * t294 - t310 + t326;
t322 = t270 * t254 + t262;
t303 = t272 * t318;
t299 = t261 * t312;
t296 = t271 * t308;
t293 = -qJD(5) - t356;
t203 = -t252 * t314 + (qJD(3) * t269 + t302) * t319;
t224 = (t249 + t306) * qJD(2);
t292 = t270 * t203 - t273 * t224;
t290 = -t273 * t203 + t218 * t312 - t270 * t224 - t228 * t311;
t288 = t244 * t305;
t287 = t246 * t305;
t253 = -t305 - t353;
t285 = -t253 - t305;
t284 = qJD(2) * t266 - t261 * t274;
t233 = t274 * t340 + t338;
t210 = -t233 * t270 - t273 * t339;
t282 = -t233 * t273 + t270 * t339;
t232 = -t269 * t274 + t271 * t340;
t281 = qJ(5) * t220 - t292;
t280 = -qJ(5) * t221 - t290;
t278 = qJD(3) * (-t285 - t353);
t277 = qJD(2) ^ 2;
t264 = -pkin(4) * t273 - pkin(3);
t256 = t354 * t273;
t255 = t354 * t270;
t250 = (pkin(4) * t270 + pkin(7)) * t271;
t243 = t273 * t254;
t241 = t244 ^ 2;
t227 = t279 * pkin(4) + pkin(7) * t313;
t212 = -t270 * t352 + t322;
t209 = t233 * qJD(3) + t289;
t208 = -t232 * qJD(3) + t274 * t302;
t206 = pkin(4) * t300 + t226;
t205 = -qJ(5) * t335 + t243 + (-pkin(4) - t355) * t274;
t201 = t217 - t293;
t192 = t210 * qJD(4) + t208 * t273 + t270 * t303;
t191 = t282 * qJD(4) - t208 * t270 + t273 * t303;
t188 = -qJD(5) * t244 + t280;
t187 = pkin(4) * t296 - t199 * qJD(4) - qJD(5) * t246 + t281;
t1 = [(-t191 * t246 - t192 * t244 + t210 * t220 + t221 * t282) * MDP(21) + (t187 * t210 - t188 * t282 + t189 * t191 + t192 * t195 + t196 * t232 + t201 * t209) * MDP(22) + (MDP(18) + MDP(20)) * (t192 * t261 + t209 * t246 - t220 * t232 + t282 * t296) + t357 * (-t191 * t261 + t209 * t244 + t210 * t296 + t221 * t232) + (-MDP(10) * t209 - MDP(11) * t208) * qJD(3) + ((-MDP(10) * t271 - MDP(11) * t274) * t275 * t308 + (-t275 * MDP(4) + (-MDP(10) * t274 + MDP(11) * t271 - MDP(3)) * t272) * t277) * t268; 0.2e1 * t295 * t364 - 0.2e1 * t308 * t363 + MDP(7) * t331 - MDP(8) * t334 + (-pkin(7) * t331 + t271 * t278) * MDP(10) + (pkin(7) * t334 + t274 * t278) * MDP(11) + (-t220 * t335 + (t274 * t309 - t298) * t246) * MDP(12) + ((-t244 * t273 - t342) * t313 + (t345 - t221 * t273 + (t244 * t270 - t246 * t273) * qJD(4)) * t271) * MDP(13) + (t261 * t298 + t220 * t274 + (t246 * t271 + t284 * t273) * qJD(3)) * MDP(14) + (t261 * t297 + t221 * t274 + (-t244 * t271 - t284 * t270) * qJD(3)) * MDP(15) + (-t261 - t316) * MDP(16) * t314 + ((t254 * t312 - t358) * t261 + ((pkin(7) * t244 + t346) * qJD(3) + (t337 + (pkin(7) * t261 + t218) * t273) * qJD(4) + t292) * t274 + (-t288 + t217 * t311 + pkin(7) * t221 + t348 + ((-pkin(7) * t336 + t243) * qJD(2) + t198) * qJD(3)) * t271) * MDP(17) + (t359 * t261 + (t217 * t309 + (qJD(3) * t246 - t299) * pkin(7) - t290) * t274 + (-t287 - t217 * t312 - pkin(7) * t220 + t347 + (-pkin(7) * t341 - t322 * qJD(2) - t199) * qJD(3)) * t271) * MDP(18) + (t221 * t250 + t227 * t244 + (t201 * t315 - t187) * t274 - t330 * t261 + (-t288 + t201 * t311 + t350 + (qJD(2) * t205 + t189) * qJD(3)) * t271) * MDP(19) + (-t220 * t250 + t227 * t246 + (t201 * t309 + t188) * t274 + t329 * t261 + (-t287 - t201 * t312 + t349 + (-qJD(2) * t212 - t195) * qJD(3)) * t271) * MDP(20) + (t205 * t220 - t212 * t221 - t330 * t246 - t329 * t244 + (-t189 * t273 - t195 * t270) * t313 + (-t187 * t273 - t188 * t270 + (t189 * t270 - t195 * t273) * qJD(4)) * t271) * MDP(21) + (t187 * t205 + t188 * t212 + t196 * t250 + (-t271 * t305 + t227) * t201 + t329 * t195 + t330 * t189) * MDP(22); (qJD(3) * t226 - t253 * t317 - t204) * MDP(10) + t285 * t316 * MDP(11) + (-t246 * t341 - t345) * MDP(12) + (-t360 * t270 + t361 * t273) * MDP(13) + (-t261 * t311 + (t261 * t333 + (-t246 + t315) * t271) * qJD(2)) * MDP(14) + (t299 + (-t261 * t336 + (t244 + t309) * t271) * qJD(2)) * MDP(15) + t261 * MDP(16) * t317 + (-pkin(3) * t221 - t347 + t291 * t261 - t226 * t244 + (pkin(8) * t341 + t346) * qJD(4) + (-t198 * t271 + (-pkin(8) * t314 - t217 * t274) * t270) * qJD(2)) * MDP(17) + (pkin(3) * t220 + t348 - t326 * t261 - t226 * t246 + (-pkin(8) * t261 * t270 + t217 * t273) * qJD(4) + (-t217 * t333 + (-pkin(8) * t309 + t199) * t271) * qJD(2)) * MDP(18) + (-t349 - t206 * t244 + t221 * t264 + t328 * t261 + (t201 + t356) * t312 + (-t201 * t336 + (qJD(3) * t255 - t189) * t271) * qJD(2)) * MDP(19) + (t350 - t206 * t246 - t220 * t264 - t327 * t261 + (pkin(4) * t342 + t201 * t273) * qJD(4) + (-t201 * t333 + (qJD(3) * t256 + t195) * t271) * qJD(2)) * MDP(20) + (t220 * t255 + t221 * t256 + t328 * t246 + t327 * t244 + (t261 * t189 + t188) * t273 + (-t187 + t351) * t270) * MDP(21) + (t187 * t255 - t188 * t256 + t196 * t264 + (pkin(4) * t312 - t206) * t201 - t327 * t195 - t328 * t189) * MDP(22) + (-t274 * t364 + t363) * t277; -t241 * MDP(13) + t323 * MDP(14) + (-t199 * t261 - t292) * MDP(17) + (-t198 * t261 + t290) * MDP(18) + (t281 - t351) * MDP(19) + (-t194 * t261 - t280) * MDP(20) + t220 * pkin(4) * MDP(21) + (pkin(4) * t187 + t362 * t195) * MDP(22) + (-t261 * MDP(14) + t217 * MDP(18) + (qJD(5) + t201) * MDP(20) - t362 * MDP(21)) * t244 + (-MDP(15) * t336 + (0.2e1 * MDP(19) * pkin(4) + MDP(16)) * t271) * t308 + (-MDP(14) * t301 - t246 * MDP(15) - t357 * t199) * qJD(4) + (t244 * MDP(12) - t261 * MDP(15) - t217 * MDP(17) + (-t201 + t293) * MDP(19) - pkin(4) * t201 * MDP(22) + (-MDP(20) * pkin(4) + MDP(13)) * t246) * t246; t360 * MDP(19) + t361 * MDP(20) + (-t246 ^ 2 - t241) * MDP(21) + (t189 * t246 + t195 * t244 + t196) * MDP(22);];
tauc = t1;
