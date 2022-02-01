% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:20:04
% EndTime: 2022-01-20 10:20:07
% DurationCPUTime: 1.71s
% Computational Cost: add. (1428->275), mult. (2253->325), div. (0->0), fcn. (1262->12), ass. (0->156)
t373 = qJDD(3) - g(3);
t290 = cos(qJ(2));
t358 = pkin(1) * qJD(2);
t327 = qJD(1) * t358;
t287 = sin(qJ(2));
t332 = qJDD(1) * t287;
t372 = pkin(1) * t332 + t290 * t327;
t289 = cos(qJ(4));
t279 = qJD(1) + qJD(2);
t359 = pkin(1) * qJD(1);
t329 = t290 * t359;
t236 = pkin(2) * t279 + t329;
t284 = cos(pkin(8));
t330 = t287 * t359;
t253 = t284 * t330;
t283 = sin(pkin(8));
t214 = t283 * t236 + t253;
t209 = pkin(7) * t279 + t214;
t319 = qJ(5) * t279 + t209;
t305 = t319 * t289;
t282 = qJ(1) + qJ(2);
t272 = sin(t282);
t261 = g(1) * t272;
t273 = cos(t282);
t371 = -g(2) * t273 + t261;
t370 = pkin(1) * t290;
t369 = pkin(2) * t284;
t286 = sin(qJ(4));
t280 = t286 ^ 2;
t368 = pkin(4) * t280;
t367 = pkin(4) * t289;
t269 = pkin(8) + t282;
t256 = sin(t269);
t366 = g(1) * t256;
t257 = cos(t269);
t365 = g(1) * t257;
t288 = sin(qJ(1));
t364 = g(1) * t288;
t363 = g(2) * t256;
t362 = g(2) * t257;
t360 = g(3) * t289;
t357 = qJD(4) * pkin(4);
t356 = qJDD(4) * pkin(4);
t252 = t283 * t330;
t213 = t236 * t284 - t252;
t208 = -pkin(3) * t279 - t213;
t355 = t208 * t279;
t265 = pkin(2) + t370;
t346 = t283 * t287;
t317 = -pkin(1) * t346 + t265 * t284;
t224 = -pkin(3) - t317;
t219 = t224 - t367;
t354 = t219 * t279;
t264 = pkin(3) + t367;
t237 = -t264 - t369;
t278 = qJDD(1) + qJDD(2);
t353 = t237 * t278;
t352 = t256 * t289;
t351 = t257 * t286;
t350 = t278 * t286;
t349 = t278 * t289;
t348 = t279 * t286;
t347 = t279 * t289;
t345 = t284 * t287;
t344 = t286 * t289;
t340 = pkin(1) * t345 + t283 * t265;
t225 = pkin(7) + t340;
t343 = -qJ(5) - t225;
t258 = pkin(2) * t283 + pkin(7);
t342 = -qJ(5) - t258;
t271 = t289 * qJD(3);
t197 = -t286 * t319 + t271;
t194 = t197 + t357;
t341 = t194 - t197;
t339 = g(1) * t273 + g(2) * t272;
t281 = t289 ^ 2;
t338 = -t280 - t281;
t337 = t280 - t281;
t336 = qJD(4) * t279;
t335 = qJD(4) * t286;
t334 = qJD(4) * t289;
t331 = qJDD(4) * t225;
t328 = pkin(4) * t335;
t266 = qJDD(1) * t370;
t223 = pkin(2) * t278 - t287 * t327 + t266;
t201 = t223 * t284 - t283 * t372;
t322 = t279 * t335;
t249 = pkin(4) * t322;
t190 = -t264 * t278 + qJDD(5) - t201 + t249;
t203 = -t264 * t279 + qJD(5) - t213;
t242 = g(2) * t351;
t325 = t190 * t286 + t203 * t334 + t242;
t195 = -pkin(3) * t278 - t201;
t324 = t195 * t286 + t208 * t334 + t242;
t226 = t283 * t329 + t253;
t228 = t284 * t329 - t252;
t243 = g(1) * t352;
t323 = t226 * t347 + t228 * t335 + t243;
t202 = t283 * t223 + t284 * t372;
t321 = -t190 - t362;
t320 = -t195 - t362;
t229 = (t284 * t290 - t346) * t358;
t318 = t224 * t279 - t229;
t316 = qJD(4) * t343;
t315 = qJD(4) * t342;
t314 = 0.2e1 * t279 * t334;
t313 = qJD(1) * (-qJD(2) + t279);
t196 = pkin(7) * t278 + t202;
t312 = -qJD(4) * qJD(3) - t196;
t310 = t266 + t371;
t309 = -t363 - t365;
t263 = pkin(2) * t273;
t285 = -qJ(5) - pkin(7);
t308 = -t256 * t285 + t257 * t264 + t263;
t307 = -t226 * t279 - t366;
t292 = qJD(4) ^ 2;
t238 = qJDD(4) * t286 + t289 * t292;
t239 = qJDD(4) * t289 - t286 * t292;
t306 = 0.2e1 * (t278 * t344 - t336 * t337) * MDP(9) + (t278 * t280 + t286 * t314) * MDP(8) + t238 * MDP(10) + t239 * MDP(11) + t278 * MDP(4);
t198 = qJD(3) * t286 + t305;
t304 = t194 * t286 - t198 * t289;
t227 = (t283 * t290 + t345) * t358;
t215 = t227 + t328;
t303 = t215 * t279 + t219 * t278;
t259 = -pkin(3) - t369;
t302 = -t258 * t292 - t259 * t278;
t268 = t289 * qJDD(3);
t301 = g(1) * t351 + t286 * t363 + t268 - t360;
t206 = t209 * t335;
t299 = -qJ(5) * t278 + t312;
t294 = qJD(5) * t279 - t299;
t188 = -t206 + (-qJ(5) * t336 + qJDD(3)) * t286 + t294 * t289;
t300 = t188 * t289 + t309;
t298 = -qJDD(4) * t258 + t259 * t336;
t297 = g(2) * t352 - t286 * t373 + t289 * t365 + t206;
t296 = t224 * t278 + t225 * t292 + t227 * t279;
t295 = g(1) * (-pkin(2) * t272 - t256 * t264 - t257 * t285);
t293 = (-qJD(5) - t203) * t279 + t299;
t291 = cos(qJ(1));
t277 = t279 ^ 2;
t275 = t291 * pkin(1);
t274 = t289 * qJ(5);
t270 = t289 * qJD(5);
t232 = t258 * t289 + t274;
t231 = t342 * t286;
t221 = -qJD(5) * t286 + t289 * t315;
t220 = t286 * t315 + t270;
t217 = t228 * t334;
t211 = t225 * t289 + t274;
t210 = t343 * t286;
t204 = t208 * t335;
t199 = t203 * t335;
t192 = (-qJD(5) - t229) * t286 + t289 * t316;
t191 = t229 * t289 + t286 * t316 + t270;
t187 = -qJD(4) * t305 - t286 * t294 + t268 + t356;
t1 = [(t278 * t290 * MDP(5) + (-qJDD(1) - t278) * MDP(6) * t287 + (MDP(18) + MDP(7)) * t364 + (MDP(5) * t287 + MDP(6) * t290) * qJD(2) * (-qJD(1) - t279)) * pkin(1) + (qJD(4) * t192 + qJDD(4) * t210 + t199 + t243) * MDP(15) + t310 * MDP(5) + t339 * MDP(6) + (t202 * t340 + t214 * t229 + t201 * t317 - t213 * t227 + pkin(2) * t261 - g(2) * (t263 + t275)) * MDP(7) + (t204 + t243) * MDP(13) + (-qJD(4) * t191 - qJDD(4) * t211 + t325) * MDP(16) + (-g(2) * t291 + t364) * MDP(2) + (g(1) * t291 + g(2) * t288) * MDP(3) + t306 + ((-t296 + t320) * MDP(13) - MDP(14) * t331 + (-t303 + t321) * MDP(15) + (t191 * t279 + t211 * t278) * MDP(17) + (t318 * MDP(14) + MDP(16) * t354 + (-t210 * t279 - t194) * MDP(17)) * qJD(4)) * t289 + (-MDP(13) * t331 + (t296 - t366) * MDP(14) + (t303 - t366) * MDP(16) + (-t192 * t279 - t210 * t278 - t187) * MDP(17) + (t318 * MDP(13) + MDP(15) * t354 + (-t211 * t279 - t198) * MDP(17)) * qJD(4)) * t286 + t300 * MDP(17) + (t188 * t211 + t198 * t191 + t187 * t210 + t194 * t192 + t190 * t219 + t203 * t215 - t295 - g(2) * (t275 + t308)) * MDP(18) + t324 * MDP(14) + qJDD(1) * MDP(1); (pkin(1) * t287 * t313 + t310) * MDP(5) + ((t290 * t313 - t332) * pkin(1) + t339) * MDP(6) + (t213 * t226 - t214 * t228 + (t201 * t284 + t202 * t283 + t371) * pkin(2)) * MDP(7) + (t204 + t298 * t286 + (t302 + t320) * t289 + t323) * MDP(13) + (t217 + t298 * t289 + (-t302 + t307) * t286 + t324) * MDP(14) + (qJDD(4) * t231 + t199 + (t237 * t348 + t221) * qJD(4) + (-t249 + t321 - t353) * t289 + t323) * MDP(15) + (-qJDD(4) * t232 + t217 + (t307 + t353) * t286 + (-t220 + (t237 * t289 + t368) * t279) * qJD(4) + t325) * MDP(16) + ((-qJD(4) * t194 + t232 * t278) * t289 + (-qJD(4) * t198 - t231 * t278 - t187) * t286 + (t220 * t289 - t221 * t286 + t338 * t228 + (-t231 * t289 - t232 * t286) * qJD(4)) * t279 + t300) * MDP(17) + (t188 * t232 + t187 * t231 + t190 * t237 - t295 - g(2) * t308 + (-t226 + t328) * t203 + (-t228 * t289 + t220) * t198 + (t228 * t286 + t221) * t194) * MDP(18) + t306; t373 * MDP(7) + (-qJD(4) * t304 + t187 * t289 + t188 * t286 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t239 + (-MDP(14) - MDP(16)) * t238; -t277 * MDP(8) * t344 + t337 * t277 * MDP(9) + MDP(10) * t350 + MDP(11) * t349 + qJDD(4) * MDP(12) + ((-t196 - t355) * t286 + t301) * MDP(13) + ((-t209 * t286 + t271) * qJD(4) + (t312 - t355) * t289 + t297) * MDP(14) + (0.2e1 * t356 + (t198 - t305) * qJD(4) + (t277 * t367 + t293) * t286 + t301) * MDP(15) + (-t277 * t368 + (qJ(5) * t348 + t197) * qJD(4) + t293 * t289 + t297) * MDP(16) + (-pkin(4) * t350 + (t341 - t357) * t347) * MDP(17) + (t341 * t198 + (-t360 + t187 + (-t203 * t279 - t309) * t286) * pkin(4)) * MDP(18); (0.2e1 * t322 - t349) * MDP(15) + (t314 + t350) * MDP(16) + (t279 * t304 - t321 - t366) * MDP(18) + t338 * MDP(17) * t277;];
tau = t1;
