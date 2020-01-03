% Calculate vector of inverse dynamics joint torques for
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:27
% EndTime: 2019-12-31 18:57:30
% DurationCPUTime: 2.48s
% Computational Cost: add. (1478->314), mult. (2833->412), div. (0->0), fcn. (1631->6), ass. (0->145)
t260 = sin(qJ(3));
t263 = cos(qJ(3));
t236 = pkin(3) * t260 - pkin(7) * t263 + qJ(2);
t259 = sin(qJ(4));
t265 = -pkin(1) - pkin(6);
t262 = cos(qJ(4));
t336 = t260 * t262;
t327 = t259 * t236 + t265 * t336;
t261 = sin(qJ(1));
t264 = cos(qJ(1));
t367 = -g(1) * t261 + g(2) * t264;
t332 = t262 * t264;
t339 = t259 * t261;
t214 = -t260 * t339 + t332;
t334 = t261 * t262;
t335 = t260 * t264;
t216 = t259 * t335 + t334;
t366 = -g(1) * t214 - g(2) * t216;
t365 = -pkin(4) * t259 + t265;
t243 = qJD(1) * t265 + qJD(2);
t318 = qJD(3) * t260;
t293 = -qJDD(3) * pkin(3) + t243 * t318;
t241 = qJDD(1) * t265 + qJDD(2);
t342 = t241 * t263;
t205 = t293 - t342;
t323 = qJD(1) * t260;
t245 = qJD(4) + t323;
t357 = g(3) * t260;
t364 = qJD(4) * pkin(7) * t245 - t367 * t263 + t205 - t357;
t311 = t262 * qJD(3);
t300 = t260 * t311;
t315 = qJD(4) * t259;
t273 = t263 * t315 + t300;
t307 = qJDD(1) * t263;
t196 = qJD(1) * t273 - qJD(4) * t311 - t259 * qJDD(3) - t262 * t307;
t319 = qJD(3) * t259;
t322 = qJD(1) * t263;
t231 = t262 * t322 + t319;
t301 = t259 * t323;
t197 = -qJD(3) * t301 + qJD(4) * t231 - t262 * qJDD(3) + t259 * t307;
t363 = t231 ^ 2;
t359 = g(1) * t264;
t356 = g(3) * t263;
t355 = qJ(5) + pkin(7);
t354 = pkin(1) * qJDD(1);
t267 = qJD(1) ^ 2;
t353 = qJ(2) * t267;
t352 = t196 * t259;
t351 = t197 * t262;
t310 = qJD(1) * qJD(3);
t298 = t263 * t310;
t308 = qJDD(1) * t260;
t227 = qJDD(4) + t298 + t308;
t350 = t227 * t259;
t349 = t227 * t262;
t229 = t259 * t322 - t311;
t348 = t229 * t245;
t347 = t229 * t259;
t346 = t229 * t262;
t345 = t231 * t245;
t344 = t231 * t259;
t343 = t231 * t262;
t341 = t243 * t263;
t340 = t245 * t259;
t338 = t259 * t263;
t337 = t259 * t264;
t235 = t260 * t243;
t333 = t262 * t263;
t213 = t236 * qJD(1);
t220 = qJD(3) * pkin(7) + t235;
t194 = t262 * t213 - t220 * t259;
t190 = -qJ(5) * t231 + t194;
t189 = pkin(4) * t245 + t190;
t331 = -t190 + t189;
t294 = qJD(4) * t355;
t313 = qJD(5) * t262;
t290 = pkin(3) * t263 + pkin(7) * t260;
t234 = t290 * qJD(1);
t328 = t259 * t234 + t243 * t333;
t330 = -qJ(5) * t301 - t259 * t294 + t313 - t328;
t219 = t262 * t234;
t329 = -qJD(5) * t259 - t262 * t294 + t243 * t338 - t219 - (pkin(4) * t263 + qJ(5) * t336) * qJD(1);
t326 = t264 * pkin(1) + t261 * qJ(2);
t257 = t263 ^ 2;
t325 = t260 ^ 2 - t257;
t266 = qJD(3) ^ 2;
t324 = -t266 - t267;
t321 = qJD(3) * t229;
t320 = qJD(3) * t231;
t317 = qJD(3) * t263;
t316 = qJD(3) * t265;
t314 = qJD(4) * t262;
t221 = -qJD(3) * pkin(3) - t341;
t203 = pkin(4) * t229 + qJD(5) + t221;
t312 = t203 * qJD(3);
t309 = qJDD(1) * qJ(2);
t306 = 0.2e1 * qJD(1) * qJD(2);
t305 = qJ(5) * t333;
t228 = qJD(3) * t290 + qJD(2);
t201 = qJD(1) * t228 + qJDD(1) * t236;
t206 = qJDD(3) * pkin(7) + t241 * t260 + t243 * t317;
t304 = -t259 * t201 - t262 * t206 - t213 * t314;
t299 = t263 * t316;
t303 = t259 * t228 + t236 * t314 + t262 * t299;
t297 = -t259 * t265 + pkin(4);
t292 = t245 * t265 + t220;
t291 = -qJD(4) * t213 - t206;
t289 = g(2) * t261 + t359;
t287 = qJDD(2) + t367;
t195 = t213 * t259 + t220 * t262;
t199 = t262 * t201;
t184 = pkin(4) * t227 + qJ(5) * t196 - qJD(4) * t195 - qJD(5) * t231 - t206 * t259 + t199;
t275 = -t220 * t315 - t304;
t185 = -qJ(5) * t197 - qJD(5) * t229 + t275;
t286 = -t184 * t259 + t185 * t262;
t191 = -qJ(5) * t229 + t195;
t285 = t189 * t262 + t191 * t259;
t284 = t189 * t259 - t191 * t262;
t247 = pkin(4) * t262 + pkin(3);
t282 = t247 * t260 - t263 * t355;
t281 = -t241 - t367;
t279 = t245 * t314 + t350;
t278 = -t245 * t315 + t349;
t277 = 0.2e1 * qJ(2) * t310 + qJDD(3) * t265;
t274 = pkin(4) * t197 + qJDD(5) + t293;
t272 = t281 + t353;
t271 = -pkin(7) * t227 + t221 * t245;
t270 = -t289 + t306 + 0.2e1 * t309;
t269 = -t265 * t266 + t270;
t268 = (t343 + t347) * MDP(21) - t285 * MDP(22) + (-t262 * MDP(19) + t259 * MDP(20)) * t245;
t253 = t264 * qJ(2);
t250 = qJDD(3) * t263;
t239 = t355 * t262;
t238 = t355 * t259;
t226 = t229 ^ 2;
t225 = t262 * t236;
t217 = t260 * t332 - t339;
t215 = t260 * t334 + t337;
t212 = t262 * t228;
t204 = -qJ(5) * t338 + t327;
t200 = t260 * t297 + t225 - t305;
t188 = t274 - t342;
t187 = -qJD(4) * t305 + (-qJD(5) * t263 + (qJ(5) * qJD(3) - qJD(4) * t265) * t260) * t259 + t303;
t186 = qJ(5) * t300 + t212 - t327 * qJD(4) + (qJ(5) * t315 + qJD(3) * t297 - t313) * t263;
t1 = [qJDD(1) * MDP(1) - t367 * MDP(2) + t289 * MDP(3) + (t287 - 0.2e1 * t354) * MDP(4) + t270 * MDP(5) + (-(qJDD(2) - t354) * pkin(1) - g(1) * (-pkin(1) * t261 + t253) - g(2) * t326 + (t306 + t309) * qJ(2)) * MDP(6) + (qJDD(1) * t257 - 0.2e1 * t260 * t298) * MDP(7) + 0.2e1 * (-t260 * t307 + t310 * t325) * MDP(8) + (-t260 * t266 + t250) * MDP(9) + (-qJDD(3) * t260 - t263 * t266) * MDP(10) + (t260 * t269 + t263 * t277) * MDP(12) + (-t260 * t277 + t263 * t269) * MDP(13) + (-t196 * t333 - t231 * t273) * MDP(14) + ((t344 + t346) * t318 + (t352 - t351 + (-t343 + t347) * qJD(4)) * t263) * MDP(15) + ((-t245 * t311 - t196) * t260 + (t278 + t320) * t263) * MDP(16) + ((t245 * t319 - t197) * t260 + (-t279 - t321) * t263) * MDP(17) + (t227 * t260 + t245 * t317) * MDP(18) + (-g(1) * t217 - g(2) * t215 + t212 * t245 + t225 * t227 + (t229 * t316 - t292 * t314 + t199) * t260 + (qJD(3) * t194 - t197 * t265 + t221 * t314) * t263 + ((-qJD(4) * t236 - t299) * t245 + t205 * t263 + (-qJD(3) * t221 - t227 * t265 + t291) * t260) * t259) * MDP(19) + (-t303 * t245 - t327 * t227 + g(1) * t216 - g(2) * t214 + (t292 * t315 + (-t221 * t262 + t231 * t265) * qJD(3) + t304) * t260 + (-qJD(3) * t195 + t196 * t265 + t205 * t262 - t221 * t315) * t263) * MDP(20) + (-t186 * t231 - t187 * t229 + t196 * t200 - t197 * t204 + t285 * t318 + (qJD(4) * t284 - t184 * t262 - t185 * t259 + t289) * t263) * MDP(21) + (t185 * t204 + t191 * t187 + t184 * t200 + t189 * t186 - g(1) * (t247 * t335 + t253) - g(2) * (pkin(4) * t337 + pkin(6) * t264 + t326) + t365 * t260 * t312 + (pkin(4) * t203 * t314 - t188 * t365 + t355 * t359) * t263 + (-g(1) * t365 - g(2) * t282) * t261) * MDP(22); qJDD(1) * MDP(4) - t267 * MDP(5) + (t287 - t353 - t354) * MDP(6) + t250 * MDP(12) + t367 * MDP(22) + t268 * qJD(1) + (t324 * MDP(13) - t197 * MDP(19) + t196 * MDP(20) - t188 * MDP(22) + ((t344 - t346) * MDP(21) - t284 * MDP(22) + (-t259 * MDP(19) - t262 * MDP(20)) * t245) * qJD(3)) * t263 + (t324 * MDP(12) - qJDD(3) * MDP(13) + (t321 - t350) * MDP(19) + (t320 - t349) * MDP(20) + (-t351 - t352) * MDP(21) + (t286 + t312) * MDP(22) + t268 * qJD(4)) * t260; MDP(9) * t307 - MDP(10) * t308 + qJDD(3) * MDP(11) + (-t263 * t272 + t357) * MDP(12) + (t260 * t272 + t356) * MDP(13) + (t245 * t343 - t352) * MDP(14) + ((-t196 - t348) * t262 + (-t197 - t345) * t259) * MDP(15) + ((-t231 * t263 + t245 * t336) * qJD(1) + t279) * MDP(16) + ((t229 * t263 - t260 * t340) * qJD(1) + t278) * MDP(17) - t245 * MDP(18) * t322 + (-t194 * t322 - t229 * t235 - pkin(3) * t197 - t219 * t245 + (t245 * t341 + t271) * t259 - t364 * t262) * MDP(19) + (pkin(3) * t196 + t195 * t322 - t231 * t235 + t328 * t245 + t364 * t259 + t271 * t262) * MDP(20) + (-t356 - t196 * t238 - t197 * t239 - t329 * t231 - t330 * t229 - t285 * qJD(4) + (-qJD(1) * t285 + t367) * t260 + t286) * MDP(21) + (t185 * t239 - t184 * t238 - t188 * t247 + g(3) * t282 + (pkin(4) * t340 - t235) * t203 + t330 * t191 + t329 * t189 + t367 * (t247 * t263 + t260 * t355)) * MDP(22) + (t263 * t260 * MDP(7) - t325 * MDP(8)) * t267; t231 * t229 * MDP(14) + (-t226 + t363) * MDP(15) + (-t196 + t348) * MDP(16) + (-t197 + t345) * MDP(17) + t227 * MDP(18) + (-t220 * t314 + t195 * t245 - t221 * t231 + t199 + (t291 + t356) * t259 + t366) * MDP(19) + (g(1) * t215 - g(2) * t217 + g(3) * t333 + t194 * t245 + t221 * t229 - t275) * MDP(20) + (pkin(4) * t196 - t229 * t331) * MDP(21) + (t331 * t191 + (g(3) * t338 - t203 * t231 + t184 + t366) * pkin(4)) * MDP(22); (-t226 - t363) * MDP(21) + (t189 * t231 + t191 * t229 + t263 * t281 + t274 - t357) * MDP(22);];
tau = t1;
