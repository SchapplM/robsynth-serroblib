% Calculate vector of inverse dynamics joint torques for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:33
% EndTime: 2021-01-15 12:06:41
% DurationCPUTime: 2.83s
% Computational Cost: add. (1837->312), mult. (3936->414), div. (0->0), fcn. (2712->14), ass. (0->143)
t357 = 2 * qJD(3);
t278 = sin(pkin(8));
t259 = pkin(1) * t278 + pkin(6);
t333 = qJ(4) + t259;
t285 = cos(qJ(3));
t346 = cos(pkin(9));
t317 = t346 * t285;
t255 = qJD(1) * t317;
t277 = sin(pkin(9));
t282 = sin(qJ(3));
t245 = t277 * t285 + t346 * t282;
t326 = qJD(1) * qJD(3);
t318 = t282 * t326;
t291 = t245 * qJDD(1) - t277 * t318;
t217 = qJD(3) * t255 + t291;
t356 = (qJD(3) * qJD(5)) + t217;
t270 = t285 * qJDD(2);
t249 = t259 * qJDD(1);
t292 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + (qJD(2) * qJD(3)) + t249;
t313 = t333 * qJD(1);
t303 = t313 * qJD(3);
t197 = qJDD(3) * pkin(3) - t292 * t282 - t285 * t303 + t270;
t200 = (qJDD(2) - t303) * t282 + t292 * t285;
t183 = t346 * t197 - t277 * t200;
t181 = -qJDD(3) * pkin(4) - t183;
t329 = qJD(1) * t282;
t237 = t277 * t329 - t255;
t235 = qJD(5) + t237;
t240 = t245 * qJD(1);
t258 = pkin(3) * t277 + pkin(7);
t273 = qJ(3) + pkin(9);
t264 = sin(t273);
t266 = cos(t273);
t274 = qJ(1) + pkin(8);
t265 = sin(t274);
t267 = cos(t274);
t308 = g(1) * t267 + g(2) * t265;
t294 = -g(3) * t266 + t308 * t264;
t355 = t294 - (pkin(3) * t329 + pkin(4) * t240 + pkin(7) * t237 + qJD(5) * t258) * t235 - t181;
t227 = t285 * qJD(2) - t313 * t282;
t347 = qJD(3) * pkin(3);
t221 = t227 + t347;
t228 = qJD(2) * t282 + t313 * t285;
t335 = t277 * t228;
t198 = t346 * t221 - t335;
t194 = -qJD(3) * pkin(4) - t198;
t315 = qJD(3) * t333;
t229 = qJD(4) * t285 - t282 * t315;
t297 = -qJD(4) * t282 - t285 * t315;
t204 = t346 * t229 + t277 * t297;
t263 = pkin(3) * t285 + pkin(2);
t279 = cos(pkin(8));
t352 = pkin(1) * t279;
t246 = -t263 - t352;
t299 = -t277 * t282 + t317;
t209 = -pkin(4) * t299 - pkin(7) * t245 + t246;
t243 = t333 * t285;
t316 = t333 * t282;
t213 = t346 * t243 - t277 * t316;
t239 = t245 * qJD(3);
t324 = qJDD(1) * t282;
t305 = -qJDD(1) * t317 + t277 * t324;
t216 = qJD(1) * t239 + t305;
t215 = qJDD(5) + t216;
t242 = t299 * qJD(3);
t184 = t277 * t197 + t346 * t200;
t182 = qJDD(3) * pkin(7) + t184;
t236 = t246 * qJD(1) + qJD(4);
t205 = pkin(4) * t237 - pkin(7) * t240 + t236;
t312 = qJD(5) * t205 + t182;
t353 = t181 * t245 + t194 * t242 - t213 * t215 - (qJD(5) * t209 + t204) * t235 + t312 * t299;
t351 = pkin(3) * t282;
t350 = g(3) * t264;
t348 = g(3) * t285;
t281 = sin(qJ(5));
t284 = cos(qJ(5));
t319 = t281 * qJDD(3) + t284 * t356;
t327 = qJD(5) * t281;
t190 = -t240 * t327 + t319;
t345 = t190 * t281;
t344 = t209 * t215;
t223 = -t284 * qJD(3) + t240 * t281;
t343 = t223 * t235;
t342 = t223 * t240;
t225 = qJD(3) * t281 + t240 * t284;
t341 = t225 * t235;
t340 = t225 * t240;
t339 = t265 * t281;
t338 = t265 * t284;
t337 = t267 * t281;
t336 = t267 * t284;
t334 = t281 * t215;
t210 = t284 * t215;
t332 = qJDD(2) - g(3);
t331 = -t190 * t299 + t225 * t239;
t219 = t346 * t228;
t199 = t277 * t221 + t219;
t275 = t282 ^ 2;
t330 = -t285 ^ 2 + t275;
t261 = -pkin(2) - t352;
t252 = qJD(1) * t261;
t328 = qJD(5) * t245;
t323 = qJDD(1) * t285;
t322 = t282 * t347;
t321 = t245 * t334;
t320 = t245 * t210;
t314 = t284 * t235;
t226 = pkin(3) * t318 + t246 * qJDD(1) + qJDD(4);
t188 = pkin(4) * t216 - pkin(7) * t217 + t226;
t195 = qJD(3) * pkin(7) + t199;
t311 = qJD(5) * t195 - t188;
t307 = g(1) * t265 - g(2) * t267;
t283 = sin(qJ(1));
t286 = cos(qJ(1));
t306 = g(1) * t283 - g(2) * t286;
t269 = t284 * qJDD(3);
t191 = t225 * qJD(5) + t217 * t281 - t269;
t304 = t191 * t299 - t223 * t239;
t302 = t210 + (-t237 * t281 - t327) * t235;
t301 = -t242 * t281 - t284 * t328;
t300 = -t242 * t284 + t245 * t327;
t296 = -qJD(1) * t252 - t249 + t308;
t295 = -qJDD(3) * t259 + t252 * t357;
t202 = t346 * t227 - t335;
t293 = -t258 * t215 + (t194 + t202) * t235;
t287 = qJD(3) ^ 2;
t290 = -0.2e1 * qJDD(1) * t261 - t259 * t287 + t307;
t280 = -qJ(4) - pkin(6);
t260 = -t346 * pkin(3) - pkin(4);
t248 = qJDD(3) * t285 - t282 * t287;
t247 = qJDD(3) * t282 + t285 * t287;
t233 = t266 * t336 + t339;
t232 = -t266 * t337 + t338;
t231 = -t266 * t338 + t337;
t230 = t266 * t339 + t336;
t212 = t243 * t277 + t346 * t316;
t207 = pkin(4) * t239 - pkin(7) * t242 + t322;
t203 = t229 * t277 - t346 * t297;
t201 = t227 * t277 + t219;
t187 = t284 * t188;
t186 = t195 * t284 + t205 * t281;
t185 = -t195 * t281 + t205 * t284;
t1 = [qJDD(1) * MDP(1) + t306 * MDP(2) + (g(1) * t286 + g(2) * t283) * MDP(3) + (t306 + (t278 ^ 2 + t279 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t275 + 0.2e1 * t285 * t318) * MDP(5) + 0.2e1 * (t282 * t323 - t330 * t326) * MDP(6) + t247 * MDP(7) + t248 * MDP(8) + (t295 * t282 + t290 * t285) * MDP(10) + (-t290 * t282 + t295 * t285) * MDP(11) + (-qJDD(3) * t212 + t216 * t246 - t226 * t299 + t236 * t239 + t307 * t266 + (t237 * t351 - t203) * qJD(3)) * MDP(12) + (-qJDD(3) * t213 + t217 * t246 + t226 * t245 + t236 * t242 - t307 * t264 + (t240 * t351 - t204) * qJD(3)) * MDP(13) + (-t183 * t245 + t184 * t299 - t198 * t242 - t199 * t239 + t203 * t240 - t204 * t237 + t212 * t217 - t213 * t216 - t308) * MDP(14) + (t184 * t213 + t199 * t204 - t183 * t212 - t198 * t203 + t226 * t246 + t236 * t322 - g(1) * (-pkin(1) * t283 - t263 * t265 - t267 * t280) - g(2) * (pkin(1) * t286 + t263 * t267 - t265 * t280)) * MDP(15) + (t190 * t245 * t284 - t300 * t225) * MDP(16) + ((-t223 * t284 - t225 * t281) * t242 + (-t345 - t191 * t284 + (t223 * t281 - t225 * t284) * qJD(5)) * t245) * MDP(17) + (-t300 * t235 + t320 + t331) * MDP(18) + (t301 * t235 + t304 - t321) * MDP(19) + (-t215 * t299 + t235 * t239) * MDP(20) + (-g(1) * t231 - g(2) * t233 + t185 * t239 - t187 * t299 + t212 * t191 + t203 * t223 + (t207 * t235 + t344 + (t194 * t245 + t195 * t299 - t213 * t235) * qJD(5)) * t284 + t353 * t281) * MDP(21) + (-g(1) * t230 - g(2) * t232 - t186 * t239 + t212 * t190 + t203 * t225 + (-(-qJD(5) * t213 + t207) * t235 - t344 - t311 * t299 - t194 * t328) * t281 + t353 * t284) * MDP(22); t332 * MDP(4) + t248 * MDP(10) - t247 * MDP(11) + (-qJD(3) * t239 + qJDD(3) * t299) * MDP(12) + (-qJD(3) * t242 - qJDD(3) * t245) * MDP(13) + (-t216 * t245 - t217 * t299 - t237 * t242 + t239 * t240) * MDP(14) + (t183 * t299 + t184 * t245 - t198 * t239 + t199 * t242 - g(3)) * MDP(15) + (-t304 - t321) * MDP(21) + (-t320 + t331) * MDP(22) + (t301 * MDP(21) + t300 * MDP(22)) * t235; MDP(7) * t324 + MDP(8) * t323 + qJDD(3) * MDP(9) + (t296 * t282 + t270 - t348) * MDP(10) + (-t332 * t282 + t296 * t285) * MDP(11) + (t201 * qJD(3) - t236 * t240 + (t346 * qJDD(3) - t237 * t329) * pkin(3) + t294 + t183) * MDP(12) + (t350 + qJD(3) * t202 + t236 * t237 + t308 * t266 + (-qJDD(3) * t277 - t240 * t329) * pkin(3) - t184) * MDP(13) + ((t199 - t201) * t240 + (-t198 + t202) * t237 + (-t216 * t277 - t346 * t217) * pkin(3)) * MDP(14) + (t198 * t201 - t199 * t202 + (t346 * t183 - t348 + t184 * t277 + (-qJD(1) * t236 + t308) * t282) * pkin(3)) * MDP(15) + (t225 * t314 + t345) * MDP(16) + ((t190 - t343) * t284 + (-t191 - t341) * t281) * MDP(17) + (t235 * t314 + t334 - t340) * MDP(18) + (t302 + t342) * MDP(19) - t235 * t240 * MDP(20) + (-t185 * t240 + t260 * t191 - t201 * t223 + t293 * t281 + t284 * t355) * MDP(21) + (t186 * t240 + t260 * t190 - t201 * t225 - t281 * t355 + t293 * t284) * MDP(22) + (-t282 * t285 * MDP(5) + t330 * MDP(6)) * qJD(1) ^ 2; (t240 * t357 + t305) * MDP(12) + ((t255 - t237) * qJD(3) + t291) * MDP(13) + (-t237 ^ 2 - t240 ^ 2) * MDP(14) + (t198 * t240 + t199 * t237 + t226 - t307) * MDP(15) + (t302 - t342) * MDP(21) + (-t235 ^ 2 * t284 - t334 - t340) * MDP(22); t225 * t223 * MDP(16) + (-t223 ^ 2 + t225 ^ 2) * MDP(17) + (t319 + t343) * MDP(18) + (t269 + t341) * MDP(19) + t215 * MDP(20) + (-g(1) * t232 + g(2) * t230 + t186 * t235 - t194 * t225 + t187) * MDP(21) + (g(1) * t233 - g(2) * t231 + t185 * t235 + t194 * t223) * MDP(22) + ((-t182 + t350) * MDP(22) + (-t240 * MDP(19) - MDP(21) * t195 - MDP(22) * t205) * qJD(5)) * t284 + (-qJD(5) * t240 * MDP(18) - t356 * MDP(19) + (-t312 + t350) * MDP(21) + t311 * MDP(22)) * t281;];
tau = t1;
