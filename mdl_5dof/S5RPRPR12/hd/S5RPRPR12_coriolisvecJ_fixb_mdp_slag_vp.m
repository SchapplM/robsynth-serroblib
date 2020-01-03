% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:21
% EndTime: 2019-12-31 18:30:27
% DurationCPUTime: 3.35s
% Computational Cost: add. (2256->294), mult. (6051->416), div. (0->0), fcn. (4642->8), ass. (0->122)
t302 = cos(pkin(8));
t352 = cos(qJ(3));
t322 = qJD(1) * t352;
t291 = t302 * t322;
t300 = sin(pkin(8));
t304 = sin(qJ(3));
t337 = t304 * t300;
t323 = qJD(1) * t337;
t269 = -t291 + t323;
t264 = qJD(5) + t269;
t350 = pkin(6) + qJ(2);
t285 = t350 * t300;
t281 = qJD(1) * t285;
t287 = t350 * t302;
t282 = qJD(1) * t287;
t247 = -t304 * t281 + t352 * t282;
t358 = qJD(3) * t247;
t336 = t304 * t302;
t280 = t300 * t352 + t336;
t271 = t280 * qJD(1);
t299 = sin(pkin(9));
t301 = cos(pkin(9));
t254 = -t301 * qJD(3) + t271 * t299;
t305 = cos(qJ(5));
t357 = t305 * t254;
t303 = sin(qJ(5));
t279 = t299 * t305 + t301 * t303;
t332 = t264 * t279;
t256 = qJD(3) * t299 + t271 * t301;
t356 = t254 * t303 - t256 * t305;
t355 = -t352 * t285 - t304 * t287;
t313 = t352 * t281 + t304 * t282;
t354 = (t300 ^ 2 + t302 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t275 = t280 * qJD(3);
t261 = qJD(1) * t275;
t335 = t305 * t301;
t338 = t299 * t303;
t277 = -t335 + t338;
t333 = t264 * t277;
t353 = -t261 * t279 + t264 * t333;
t265 = t269 ^ 2;
t351 = pkin(7) * t301;
t349 = pkin(7) + qJ(4);
t215 = t256 * t303 + t357;
t347 = t215 * t271;
t346 = t356 * t271;
t290 = qJD(3) * t291;
t260 = -qJD(3) * t323 + t290;
t345 = t260 * t299;
t344 = t260 * t301;
t342 = t269 * t299;
t312 = t302 * t352 - t337;
t274 = t312 * qJD(3);
t341 = t274 * t299;
t340 = t280 * t299;
t339 = t280 * t301;
t214 = pkin(3) * t261 - qJ(4) * t260 - qJD(4) * t271;
t308 = t312 * qJD(2);
t218 = qJD(1) * t308 + (qJD(4) - t313) * qJD(3);
t190 = t299 * t214 + t301 * t218;
t223 = pkin(3) * t275 - qJ(4) * t274 - qJD(4) * t280;
t228 = t355 * qJD(3) + t308;
t196 = t299 * t223 + t301 * t228;
t295 = -pkin(2) * t302 - pkin(1);
t283 = qJD(1) * t295 + qJD(2);
t227 = pkin(3) * t269 - qJ(4) * t271 + t283;
t242 = qJD(3) * qJ(4) + t247;
t200 = t299 * t227 + t301 * t242;
t243 = pkin(3) * t271 + qJ(4) * t269;
t206 = t299 * t243 - t301 * t313;
t244 = -pkin(3) * t312 - qJ(4) * t280 + t295;
t253 = -t304 * t285 + t287 * t352;
t208 = t299 * t244 + t301 * t253;
t328 = qJD(5) * t305;
t334 = -t254 * t328 + t260 * t335;
t191 = -pkin(7) * t254 + t200;
t330 = qJD(5) * t191;
t329 = qJD(5) * t280;
t326 = qJD(1) * qJD(2);
t189 = t301 * t214 - t218 * t299;
t195 = t301 * t223 - t228 * t299;
t199 = t301 * t227 - t242 * t299;
t205 = t301 * t243 + t299 * t313;
t207 = t301 * t244 - t253 * t299;
t185 = -pkin(7) * t345 + t190;
t187 = pkin(4) * t269 - pkin(7) * t256 + t199;
t320 = -qJD(5) * t187 - t185;
t221 = t300 * qJD(2) * t322 + t326 * t336 + t358;
t319 = -t277 * t261 - t264 * t332;
t181 = t187 * t305 - t191 * t303;
t182 = t187 * t303 + t191 * t305;
t197 = -pkin(4) * t312 - pkin(7) * t339 + t207;
t201 = -pkin(7) * t340 + t208;
t318 = t197 * t305 - t201 * t303;
t317 = t197 * t303 + t201 * t305;
t316 = -t199 * t299 + t200 * t301;
t315 = -qJD(5) * t256 - t345;
t286 = t349 * t301;
t311 = pkin(4) * t271 + qJD(4) * t299 + qJD(5) * t286 + t269 * t351 + t205;
t284 = t349 * t299;
t310 = pkin(7) * t342 - qJD(4) * t301 + qJD(5) * t284 + t206;
t240 = -qJD(3) * pkin(3) + qJD(4) + t313;
t309 = t221 * t280 + t240 * t274 - t260 * t355;
t307 = -pkin(3) * t260 - qJ(4) * t261 + (-qJD(4) + t240) * t269;
t194 = -qJD(5) * t356 + t260 * t279;
t229 = qJD(2) * t280 + qJD(3) * t253;
t294 = -pkin(4) * t301 - pkin(3);
t237 = t277 * t280;
t236 = t279 * t280;
t230 = pkin(4) * t340 - t355;
t222 = -pkin(4) * t342 + t247;
t210 = t254 * pkin(4) + t240;
t209 = pkin(4) * t341 + t229;
t204 = pkin(4) * t345 + t221;
t203 = t274 * t279 + t328 * t339 - t329 * t338;
t202 = -t274 * t277 - t279 * t329;
t193 = t303 * t315 + t334;
t188 = -pkin(7) * t341 + t196;
t186 = pkin(4) * t275 - t274 * t351 + t195;
t184 = pkin(4) * t261 - pkin(7) * t344 + t189;
t183 = t305 * t184;
t1 = [(t260 * t280 + t271 * t274) * MDP(8) + (t260 * t312 - t261 * t280 - t269 * t274 - t271 * t275) * MDP(9) + (t261 * t295 + t275 * t283) * MDP(13) + (t260 * t295 + t274 * t283) * MDP(14) + (-t189 * t312 + t195 * t269 + t199 * t275 + t207 * t261 + t229 * t254 + t299 * t309) * MDP(15) + (t190 * t312 - t196 * t269 - t200 * t275 - t208 * t261 + t229 * t256 + t301 * t309) * MDP(16) + (-t195 * t256 - t196 * t254 + (-t189 * t280 - t199 * t274 - t207 * t260) * t301 + (-t190 * t280 - t200 * t274 - t208 * t260) * t299) * MDP(17) + (t189 * t207 + t190 * t208 + t195 * t199 + t196 * t200 - t221 * t355 + t229 * t240) * MDP(18) + (-t193 * t237 - t202 * t356) * MDP(19) + (-t193 * t236 + t194 * t237 - t202 * t215 + t203 * t356) * MDP(20) + (-t193 * t312 + t202 * t264 - t237 * t261 - t275 * t356) * MDP(21) + (t194 * t312 - t203 * t264 - t215 * t275 - t236 * t261) * MDP(22) + (-t261 * t312 + t264 * t275) * MDP(23) + ((t186 * t305 - t188 * t303) * t264 + t318 * t261 - (-t185 * t303 + t183) * t312 + t181 * t275 + t209 * t215 + t230 * t194 + t204 * t236 + t210 * t203 + (t182 * t312 - t264 * t317) * qJD(5)) * MDP(24) + (-(t186 * t303 + t188 * t305) * t264 - t317 * t261 + (t184 * t303 + t185 * t305) * t312 - t182 * t275 - t209 * t356 + t230 * t193 - t204 * t237 + t210 * t202 + (t181 * t312 - t264 * t318) * qJD(5)) * MDP(25) + 0.2e1 * t326 * t354 + (t274 * MDP(10) - t275 * MDP(11) - t229 * MDP(13) - t228 * MDP(14)) * qJD(3); 0.2e1 * t271 * qJD(3) * MDP(13) + (t290 + (-t269 - t323) * qJD(3)) * MDP(14) + (-t254 * t271 + t261 * t301 - t265 * t299) * MDP(15) + (-t256 * t271 - t261 * t299 - t265 * t301) * MDP(16) + ((-t254 * t301 + t256 * t299) * t269 + (-t299 ^ 2 - t301 ^ 2) * t260) * MDP(17) + (t189 * t301 + t190 * t299 - t240 * t271 + t269 * t316) * MDP(18) + (t319 - t347) * MDP(24) + (t346 + t353) * MDP(25) - qJD(1) ^ 2 * t354; -t265 * MDP(9) + (t290 + (t269 - t323) * qJD(3)) * MDP(10) + (-t221 + t358) * MDP(13) + (t283 * t269 - t312 * t326) * MDP(14) + (-t205 * t269 - t221 * t301 - t247 * t254 + t299 * t307) * MDP(15) + (t206 * t269 + t221 * t299 - t247 * t256 + t301 * t307) * MDP(16) + (t205 * t256 + t206 * t254 + (-qJD(4) * t254 - t199 * t269 + t190) * t301 + (qJD(4) * t256 - t200 * t269 - t189) * t299) * MDP(17) + (-pkin(3) * t221 - t199 * t205 - t200 * t206 - t240 * t247 + t316 * qJD(4) + (-t189 * t299 + t190 * t301) * qJ(4)) * MDP(18) + (t193 * t279 + t333 * t356) * MDP(19) + (-t193 * t277 - t194 * t279 + t215 * t333 + t332 * t356) * MDP(20) + (t346 - t353) * MDP(21) + (t319 + t347) * MDP(22) + ((-t284 * t305 - t286 * t303) * t261 + t294 * t194 + t204 * t277 - t222 * t215 + (t303 * t310 - t305 * t311) * t264 + t332 * t210) * MDP(24) + (-(-t284 * t303 + t286 * t305) * t261 + t294 * t193 + t204 * t279 + t222 * t356 + (t303 * t311 + t305 * t310) * t264 - t333 * t210) * MDP(25) + (-MDP(13) * t283 - MDP(15) * t199 + MDP(16) * t200 - MDP(23) * t264 - MDP(24) * t181 + MDP(25) * t182 + MDP(8) * t269 + MDP(9) * t271) * t271; (t256 * t269 + t345) * MDP(15) + (-t254 * t269 + t344) * MDP(16) + (-t254 ^ 2 - t256 ^ 2) * MDP(17) + (t199 * t256 + t200 * t254 + t221) * MDP(18) + (-t264 * t356 + t194) * MDP(24) + (-t264 * t357 + (-t256 * t264 + t315) * t303 + t334) * MDP(25); -t215 ^ 2 * MDP(20) + (t215 * t264 + t334) * MDP(21) + t261 * MDP(23) + (t182 * t264 + t183) * MDP(24) + (t181 * t264 + t210 * t215) * MDP(25) - (MDP(19) * t215 - MDP(20) * t356 + MDP(22) * t264 - MDP(24) * t210) * t356 + (MDP(22) * t315 - MDP(24) * t330 + MDP(25) * t320) * t305 + (t315 * MDP(21) + (qJD(5) * t254 - t344) * MDP(22) + t320 * MDP(24) + (-t184 + t330) * MDP(25)) * t303;];
tauc = t1;
