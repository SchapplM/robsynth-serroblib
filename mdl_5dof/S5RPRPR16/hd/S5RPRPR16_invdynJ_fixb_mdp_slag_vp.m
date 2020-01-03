% Calculate vector of inverse dynamics joint torques for
% S5RPRPR16
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR16_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR16_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR16_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:39
% EndTime: 2019-12-31 18:39:42
% DurationCPUTime: 2.44s
% Computational Cost: add. (975->318), mult. (1768->389), div. (0->0), fcn. (956->6), ass. (0->153)
t261 = cos(qJ(3));
t316 = qJD(1) * qJD(3);
t301 = t261 * t316;
t258 = sin(qJ(3));
t314 = qJDD(1) * t258;
t369 = t301 + t314;
t264 = -pkin(1) - pkin(6);
t230 = t264 * qJD(1) + qJD(2);
t201 = (pkin(4) * qJD(1) - t230) * t261;
t318 = qJD(4) + t201;
t320 = qJD(5) * t258;
t368 = qJD(1) * t320 + qJDD(3);
t262 = cos(qJ(1));
t250 = g(2) * t262;
t259 = sin(qJ(1));
t251 = g(1) * t259;
t367 = t251 - t250;
t302 = t258 * t316;
t309 = t261 * qJDD(1);
t271 = -t302 + t309;
t215 = -qJDD(5) - t271;
t260 = cos(qJ(5));
t202 = t260 * t215;
t330 = qJD(1) * t261;
t236 = qJD(5) + t330;
t257 = sin(qJ(5));
t346 = t236 * t257;
t273 = -qJD(5) * t346 - t202;
t248 = t258 * pkin(3);
t353 = qJ(4) * t261;
t366 = -t353 + t248;
t252 = qJDD(1) * qJ(2);
t253 = qJD(1) * qJD(2);
t307 = 0.2e1 * t253;
t365 = 0.2e1 * t252 + t307;
t364 = t264 * qJDD(1);
t222 = t258 * t230;
t331 = qJD(1) * t258;
t200 = -pkin(4) * t331 + t222;
t329 = qJD(3) * qJ(4);
t196 = t200 + t329;
t263 = -pkin(3) - pkin(7);
t363 = -t263 * t215 + (t196 - t200) * t236;
t229 = qJDD(2) + t364;
t221 = t258 * t229;
t313 = qJDD(3) * qJ(4);
t288 = -t221 - t313;
t347 = t230 * t261;
t191 = (-qJD(4) - t347) * qJD(3) + t288;
t324 = qJD(3) * t258;
t308 = t230 * t324 + qJDD(4);
t281 = -t229 * t261 + t308;
t352 = qJDD(3) * pkin(3);
t192 = t281 - t352;
t290 = -qJD(4) + t347;
t356 = qJD(3) * pkin(3);
t203 = -t290 - t356;
t205 = -t222 - t329;
t362 = -t191 * t258 - t192 * t261 + (t203 * t258 - t205 * t261) * qJD(3);
t338 = pkin(3) * t331 + qJD(1) * qJ(2);
t204 = -qJ(4) * t330 + t338;
t223 = qJ(2) + t366;
t310 = qJDD(3) * t264;
t361 = (qJD(1) * t223 + t204) * qJD(3) + t310;
t359 = g(3) * t258;
t358 = g(3) * t261;
t357 = pkin(4) - t264;
t355 = pkin(1) * qJDD(1);
t266 = qJD(1) ^ 2;
t354 = qJ(2) * t266;
t187 = -pkin(4) * t314 + (qJD(4) - t201) * qJD(3) - t288;
t351 = t187 * t258;
t287 = t257 * t369 + t368 * t260;
t315 = qJD(3) * qJD(5);
t189 = -t257 * t315 + t287;
t350 = t189 * t260;
t325 = qJD(3) * t257;
t216 = -t260 * t331 + t325;
t349 = t216 * t236;
t323 = qJD(3) * t260;
t218 = t257 * t331 + t323;
t348 = t218 * t236;
t345 = t257 * t215;
t344 = t257 * t258;
t343 = t258 * t266;
t342 = t259 * t261;
t341 = t260 * t261;
t340 = t261 * t262;
t339 = t369 * t260;
t220 = pkin(3) * t330 + qJ(4) * t331;
t337 = t262 * pkin(1) + t259 * qJ(2);
t255 = t258 ^ 2;
t256 = t261 ^ 2;
t335 = -t255 - t256;
t334 = t255 - t256;
t265 = qJD(3) ^ 2;
t333 = t265 + t266;
t332 = MDP(24) * t260;
t328 = qJD(3) * t205;
t327 = qJD(3) * t216;
t326 = qJD(3) * t218;
t322 = qJD(3) * t261;
t283 = pkin(7) * t258 - t353;
t195 = t283 * qJD(1) + t338;
t321 = qJD(5) * t195;
t319 = qJD(5) * t260;
t317 = qJ(4) * qJDD(1);
t312 = qJDD(3) * t258;
t311 = qJDD(3) * t261;
t306 = t236 * t325;
t305 = t236 * t323;
t303 = pkin(3) * t322 + qJ(4) * t324 + qJD(2);
t299 = qJDD(2) - t367;
t298 = -t229 - t250;
t237 = g(1) * t342;
t297 = -t237 + t359;
t296 = qJD(3) * t357;
t294 = qJD(1) * t220 - g(3);
t293 = qJD(3) * pkin(7) - qJD(4);
t280 = pkin(3) * t369 + qJ(4) * t302 + t252 + t253;
t184 = pkin(7) * t314 + (t293 * qJD(1) - t317) * t261 + t280;
t193 = t263 * qJD(3) + t318;
t289 = qJD(5) * t193 + t184;
t286 = g(1) * t262 + g(2) * t259;
t284 = pkin(3) * t261 + qJ(4) * t258;
t282 = -t354 + t250;
t279 = -t321 - t359;
t183 = t193 * t257 + t195 * t260;
t274 = t236 * t319 - t345;
t272 = 0.2e1 * qJ(2) * t316 + t310;
t270 = -t264 * t265 - t286;
t269 = t270 + t365;
t188 = (-qJD(1) * qJD(4) - t317) * t261 + t280;
t198 = -qJD(4) * t261 + t303;
t268 = -qJD(1) * t198 - qJDD(1) * t223 - t188 - t270;
t267 = -t358 + t187 - t367 * t258 + (pkin(7) * t330 - qJD(5) * t263 + t220) * t236;
t247 = t262 * qJ(2);
t225 = t357 * t261;
t224 = t357 * t258;
t214 = t261 * t296;
t213 = t258 * t296;
t212 = qJ(2) + t248 + t283;
t209 = -t257 * t342 + t260 * t262;
t208 = -t257 * t262 - t259 * t341;
t207 = -t257 * t340 - t259 * t260;
t206 = t257 * t259 - t260 * t340;
t197 = t204 * t330;
t194 = t293 * t261 + t303;
t190 = t218 * qJD(5) + qJDD(3) * t257 - t339;
t186 = t271 * pkin(4) + t263 * qJDD(3) + t281;
t185 = t260 * t186;
t182 = t193 * t260 - t195 * t257;
t1 = [qJDD(1) * MDP(1) + t367 * MDP(2) + t286 * MDP(3) + (t299 - 0.2e1 * t355) * MDP(4) + (-t286 + t365) * MDP(5) + (-(qJDD(2) - t355) * pkin(1) - g(1) * (-pkin(1) * t259 + t247) - g(2) * t337 + (t307 + t252) * qJ(2)) * MDP(6) + (qJDD(1) * t256 - 0.2e1 * t258 * t301) * MDP(7) + 0.2e1 * (-t258 * t309 + t334 * t316) * MDP(8) + (-t258 * t265 + t311) * MDP(9) + (-t261 * t265 - t312) * MDP(10) + (t269 * t258 + t272 * t261) * MDP(12) + (-t272 * t258 + t269 * t261) * MDP(13) + (t335 * t364 - t362 + t367) * MDP(14) + (t268 * t258 - t361 * t261) * MDP(15) + (t361 * t258 + t268 * t261) * MDP(16) + (t188 * t223 + t204 * t198 - g(1) * (-qJ(4) * t340 + t262 * t248 + t247) - g(2) * (pkin(6) * t262 + t337) + (-g(1) * t264 - g(2) * t366) * t259 + t362 * t264) * MDP(17) + (t189 * t344 + (t257 * t322 + t258 * t319) * t218) * MDP(18) + ((-t216 * t257 + t218 * t260) * t322 + (t350 - t190 * t257 + (-t216 * t260 - t218 * t257) * qJD(5)) * t258) * MDP(19) + ((t189 + t306) * t261 + (t274 - t326) * t258) * MDP(20) + ((-t190 + t305) * t261 + (t273 + t327) * t258) * MDP(21) + (-t215 * t261 - t236 * t324) * MDP(22) + ((-t194 * t257 - t213 * t260) * t236 - (-t212 * t257 + t225 * t260) * t215 + (-t184 * t257 + t185) * t261 - t214 * t216 - t224 * t190 - t260 * t351 - g(1) * t207 - g(2) * t209 + (-t182 * t258 - t196 * t341) * qJD(3) + ((-t212 * t260 - t225 * t257) * t236 - t183 * t261 + t196 * t344) * qJD(5)) * MDP(23) + (t183 * t324 - g(1) * t206 - g(2) * t208 - t224 * t189 - t214 * t218 + (-(qJD(5) * t225 + t194) * t236 + t212 * t215 - t289 * t261 + t196 * t320) * t260 + (-(-qJD(5) * t212 - t213) * t236 + t225 * t215 + t351 + (qJD(3) * t196 - t186 + t321) * t261) * t257) * MDP(24); -t266 * MDP(5) + (t299 - t354) * MDP(6) - t367 * MDP(17) + (MDP(12) - MDP(15)) * (-t333 * t258 + t311) + (-MDP(13) + MDP(16)) * (t333 * t261 + t312) + (-t204 * MDP(17) + (MDP(23) * t257 + t332) * t236) * qJD(1) + ((qJD(3) * t203 - t191) * MDP(17) + (t190 + t305) * MDP(23) + (t189 - t306) * MDP(24)) * t258 + (t335 * MDP(14) - pkin(1) * MDP(6) + MDP(4)) * qJDD(1) + ((-t192 - t328) * MDP(17) + (-t273 + t327) * MDP(23) + (t274 + t326) * MDP(24)) * t261; t261 * MDP(7) * t343 - t334 * MDP(8) * t266 + MDP(9) * t309 - MDP(10) * t314 + qJDD(3) * MDP(11) + ((t229 + t282) * t261 + t297) * MDP(12) + (t358 - t221 + (-t282 + t251) * t258) * MDP(13) + (-t284 * qJDD(1) + ((-t205 - t329) * t261 + (-qJD(4) + t203 + t356) * t258) * qJD(1)) * MDP(14) + (t294 * t258 + t298 * t261 + qJDD(4) + t197 + t237 - 0.2e1 * t352) * MDP(15) + (0.2e1 * t313 + 0.2e1 * qJD(3) * qJD(4) + t221 + t294 * t261 + (-qJD(1) * t204 - t367) * t258) * MDP(16) + (-t192 * pkin(3) + g(3) * t366 - t191 * qJ(4) - t203 * t222 - t204 * t220 + t290 * t205 - t367 * t284) * MDP(17) + (-t218 * t346 + t350) * MDP(18) + ((-t190 - t348) * t260 + (-t189 + t349) * t257) * MDP(19) + ((t218 * t258 - t261 * t346) * qJD(1) + t273) * MDP(20) + ((-t216 * t258 - t236 * t341) * qJD(1) - t274) * MDP(21) + t236 * MDP(22) * t331 + (qJ(4) * t190 + t182 * t331 + t318 * t216 + t267 * t257 + t363 * t260) * MDP(23) + (qJ(4) * t189 - t183 * t331 + t318 * t218 - t363 * t257 + t267 * t260) * MDP(24); qJDD(3) * MDP(15) + (-t256 * t266 - t265) * MDP(16) + (t197 - t297 + t308 + t328 - t352) * MDP(17) + (-t202 - t327) * MDP(23) + (-t326 + t345) * MDP(24) + (qJDD(1) * MDP(14) - MDP(15) * t343 + t298 * MDP(17)) * t261 + (-MDP(23) * t346 - t236 * t332) * t236; t218 * t216 * MDP(18) + (-t216 ^ 2 + t218 ^ 2) * MDP(19) + (t287 + t349) * MDP(20) + (t339 + t348) * MDP(21) - t215 * MDP(22) + (-g(1) * t208 + g(2) * t206 + t183 * t236 - t196 * t218 + t185) * MDP(23) + (g(1) * t209 - g(2) * t207 + t182 * t236 + t196 * t216) * MDP(24) + (-MDP(21) * t315 + t279 * MDP(23) - t289 * MDP(24)) * t260 + (-MDP(20) * t315 - t368 * MDP(21) - t289 * MDP(23) + (-t186 - t279) * MDP(24)) * t257;];
tau = t1;
