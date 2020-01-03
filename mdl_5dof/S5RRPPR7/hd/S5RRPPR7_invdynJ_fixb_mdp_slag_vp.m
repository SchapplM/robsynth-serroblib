% Calculate vector of inverse dynamics joint torques for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:30
% EndTime: 2019-12-31 19:36:35
% DurationCPUTime: 3.94s
% Computational Cost: add. (1922->360), mult. (4429->452), div. (0->0), fcn. (3077->10), ass. (0->160)
t302 = sin(pkin(8));
t306 = sin(qJ(2));
t360 = qJD(1) * t306;
t303 = cos(pkin(8));
t309 = cos(qJ(2));
t366 = t303 * t309;
t262 = -qJD(1) * t366 + t302 * t360;
t305 = sin(qJ(5));
t308 = cos(qJ(5));
t246 = qJD(2) * t305 - t308 * t262;
t273 = t302 * t309 + t303 * t306;
t265 = t273 * qJD(1);
t392 = qJD(5) + t265;
t393 = t246 * t392;
t248 = qJD(2) * t308 + t262 * t305;
t338 = t392 * t248;
t391 = t305 * t392;
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t389 = g(1) * t307 - g(2) * t310;
t351 = qJD(1) * qJD(2);
t344 = t306 * t351;
t280 = t302 * t344;
t343 = t309 * t351;
t240 = qJDD(1) * t273 + t303 * t343 - t280;
t388 = -qJ(4) * t240 - qJD(4) * t265;
t335 = g(1) * t310 + g(2) * t307;
t304 = -qJ(3) - pkin(6);
t279 = t304 * t309;
t276 = qJD(1) * t279;
t268 = t302 * t276;
t278 = t304 * t306;
t275 = qJD(1) * t278;
t243 = t275 * t303 + t268;
t354 = -qJD(4) + t243;
t271 = qJD(2) * pkin(2) + t275;
t367 = t303 * t276;
t238 = t302 * t271 - t367;
t229 = -qJD(2) * qJ(4) - t238;
t382 = pkin(4) * t262;
t216 = -t229 - t382;
t235 = qJDD(5) + t240;
t242 = t275 * t302 - t367;
t290 = -pkin(2) * t303 - pkin(3);
t285 = -pkin(7) + t290;
t387 = t285 * t235 + (t216 - t242 + t382) * t392;
t261 = t265 ^ 2;
t385 = pkin(3) + pkin(7);
t384 = pkin(2) * t309;
t264 = t273 * qJD(2);
t348 = qJDD(1) * t309;
t282 = t303 * t348;
t349 = qJDD(1) * t306;
t239 = qJD(1) * t264 + t302 * t349 - t282;
t383 = pkin(3) * t239;
t381 = pkin(4) * t265;
t298 = qJ(2) + pkin(8);
t295 = cos(t298);
t288 = g(3) * t295;
t377 = g(3) * t309;
t375 = qJDD(2) * pkin(3);
t355 = qJD(5) * t308;
t345 = t308 * qJDD(2) + t305 * t239 + t262 * t355;
t350 = qJD(2) * qJD(5);
t208 = -t305 * t350 + t345;
t374 = t208 * t308;
t368 = t302 * t306;
t272 = -t366 + t368;
t291 = pkin(1) + t384;
t328 = -qJ(4) * t273 - t291;
t217 = t272 * t385 + t328;
t373 = t217 * t235;
t372 = t235 * t305;
t371 = t246 * t262;
t370 = t248 * t262;
t369 = t272 * t305;
t365 = t305 * t307;
t364 = t305 * t310;
t363 = t307 * t308;
t228 = t308 * t235;
t362 = t308 * t310;
t341 = qJD(2) * t304;
t260 = -qJD(3) * t306 + t309 * t341;
t234 = qJDD(2) * pkin(2) + qJD(1) * t260 + qJDD(1) * t278;
t259 = qJD(3) * t309 + t306 * t341;
t241 = qJD(1) * t259 - qJDD(1) * t279;
t206 = t302 * t234 + t303 * t241;
t300 = t306 ^ 2;
t361 = -t309 ^ 2 + t300;
t359 = qJD(2) * t306;
t277 = -qJD(1) * t291 + qJD(3);
t319 = -qJ(4) * t265 + t277;
t210 = t262 * t385 + t319;
t357 = qJD(5) * t210;
t356 = qJD(5) * t272;
t353 = t381 - t354;
t347 = pkin(2) * t344 + qJDD(3);
t293 = pkin(2) * t359;
t346 = qJDD(2) * qJ(4) + t206;
t340 = pkin(2) * t360 + qJ(4) * t262;
t205 = t234 * t303 - t302 * t241;
t222 = t259 * t302 - t303 * t260;
t237 = t271 * t303 + t268;
t244 = -t303 * t278 - t279 * t302;
t339 = t308 * t392;
t322 = -qJDD(1) * t291 + t347;
t316 = t322 + t388;
t196 = t239 * t385 + t316;
t333 = qJD(4) - t237;
t211 = -qJD(2) * t385 + t333 + t381;
t337 = qJD(5) * t211 + t196;
t336 = MDP(22) * t392;
t294 = sin(t298);
t332 = pkin(3) * t295 + qJ(4) * t294;
t331 = qJDD(4) - t205;
t330 = -t357 + t288;
t198 = t210 * t308 + t211 * t305;
t223 = t259 * t303 + t260 * t302;
t245 = t278 * t302 - t279 * t303;
t327 = -t389 + t347;
t203 = -qJD(2) * qJD(4) - t346;
t326 = -0.2e1 * pkin(1) * t351 - pkin(6) * qJDD(2);
t325 = t264 * t305 + t272 * t355;
t267 = qJD(2) * t366 - t302 * t359;
t324 = -qJ(4) * t267 - qJD(4) * t273 + t293;
t201 = -pkin(4) * t239 - t203;
t224 = pkin(4) * t273 + t244;
t321 = t201 * t272 + t216 * t264 - t224 * t235;
t320 = -g(3) * t294 - t295 * t335;
t311 = qJD(2) ^ 2;
t318 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t311 + t389;
t312 = qJD(1) ^ 2;
t317 = pkin(1) * t312 - pkin(6) * qJDD(1) + t335;
t315 = t222 * t265 - t223 * t262 - t239 * t245 + t240 * t244 - t335;
t220 = pkin(3) * t262 + t319;
t314 = t220 * t265 - t294 * t335 + t288 + t331;
t313 = t201 + (-qJD(5) * t285 + t265 * t385 + t340) * t392 + t320;
t287 = pkin(2) * t302 + qJ(4);
t281 = t310 * t291;
t258 = -t294 * t365 + t362;
t257 = t294 * t363 + t364;
t256 = t294 * t364 + t363;
t255 = t294 * t362 - t365;
t253 = qJD(2) * t262;
t236 = pkin(3) * t272 + t328;
t231 = t308 * t239;
t226 = -qJD(2) * pkin(3) + t333;
t225 = -pkin(4) * t272 + t245;
t221 = pkin(3) * t265 + t340;
t215 = pkin(3) * t264 + t324;
t214 = -pkin(4) * t264 + t223;
t213 = pkin(4) * t267 + t222;
t209 = t248 * qJD(5) + qJDD(2) * t305 - t231;
t207 = t264 * t385 + t324;
t204 = t331 - t375;
t202 = t316 + t383;
t200 = pkin(4) * t240 - qJDD(2) * t385 + t331;
t199 = t308 * t200;
t197 = -t210 * t305 + t211 * t308;
t1 = [qJDD(1) * MDP(1) + t389 * MDP(2) + t335 * MDP(3) + (qJDD(1) * t300 + 0.2e1 * t306 * t343) * MDP(4) + 0.2e1 * (t306 * t348 - t351 * t361) * MDP(5) + (qJDD(2) * t306 + t309 * t311) * MDP(6) + (qJDD(2) * t309 - t306 * t311) * MDP(7) + (t306 * t326 + t309 * t318) * MDP(9) + (-t306 * t318 + t309 * t326) * MDP(10) + (-t205 * t273 - t206 * t272 - t237 * t267 - t238 * t264 + t315) * MDP(11) + (t206 * t245 + t238 * t223 - t205 * t244 - t237 * t222 - t322 * t291 + t277 * t293 - g(1) * (-t291 * t307 - t304 * t310) - g(2) * (-t304 * t307 + t281)) * MDP(12) + (t203 * t272 + t204 * t273 + t226 * t267 + t229 * t264 + t315) * MDP(13) + (qJD(2) * t222 + qJDD(2) * t244 - t202 * t272 - t215 * t262 - t220 * t264 - t236 * t239 - t295 * t389) * MDP(14) + (qJD(2) * t223 + qJDD(2) * t245 - t202 * t273 - t215 * t265 - t220 * t267 - t236 * t240 + t294 * t389) * MDP(15) + (-g(2) * t281 + t202 * t236 - t203 * t245 + t204 * t244 + t220 * t215 + t226 * t222 - t229 * t223 + (g(1) * t304 - g(2) * t332) * t310 + (-g(1) * (-t291 - t332) + g(2) * t304) * t307) * MDP(16) + (t208 * t369 + t248 * t325) * MDP(17) + ((-t246 * t305 + t248 * t308) * t264 + (t374 - t209 * t305 + (-t246 * t308 - t248 * t305) * qJD(5)) * t272) * MDP(18) + (t208 * t273 + t235 * t369 + t248 * t267 + t325 * t392) * MDP(19) + (t272 * t228 - t209 * t273 - t246 * t267 + (t264 * t308 - t305 * t356) * t392) * MDP(20) + (t235 * t273 + t267 * t392) * MDP(21) + (-g(1) * t258 - g(2) * t256 + t197 * t267 + t199 * t273 + t225 * t209 + t214 * t246 + (-t196 * t273 - t207 * t392 - t373) * t305 + (t213 * t392 - t321) * t308 + ((-t217 * t308 - t224 * t305) * t392 - t198 * t273 + t216 * t369) * qJD(5)) * MDP(22) + (g(1) * t257 - g(2) * t255 - t198 * t267 + t225 * t208 + t214 * t248 + (-(qJD(5) * t224 + t207) * t392 - t373 - t337 * t273 + t216 * t356) * t308 + (-(-qJD(5) * t217 + t213) * t392 - (t200 - t357) * t273 + t321) * t305) * MDP(23); MDP(6) * t349 + MDP(7) * t348 + qJDD(2) * MDP(8) + (t306 * t317 - t377) * MDP(9) + (g(3) * t306 + t309 * t317) * MDP(10) + ((t238 - t242) * t265 + (-t237 + t243) * t262 + (-t239 * t302 - t240 * t303) * pkin(2)) * MDP(11) + (t237 * t242 - t238 * t243 + (-t377 + t205 * t303 + t206 * t302 + (-qJD(1) * t277 + t335) * t306) * pkin(2)) * MDP(12) + (-t239 * t287 + t240 * t290 + (-t229 - t242) * t265 + (t226 + t354) * t262) * MDP(13) + (-qJD(2) * t242 + t221 * t262 + (-pkin(3) + t290) * qJDD(2) + t314) * MDP(14) + (qJDD(2) * t287 - t220 * t262 + t221 * t265 + (0.2e1 * qJD(4) - t243) * qJD(2) + t320 + t346) * MDP(15) + (-t203 * t287 + t204 * t290 - t220 * t221 - t226 * t242 - g(3) * (t332 + t384) + t354 * t229 + t335 * (pkin(2) * t306 + pkin(3) * t294 - qJ(4) * t295)) * MDP(16) + (-t305 * t338 + t374) * MDP(17) + ((-t209 - t338) * t308 + (-t208 + t393) * t305) * MDP(18) + (-t391 * t392 + t228 + t370) * MDP(19) + (-t339 * t392 - t371 - t372) * MDP(20) + t392 * t262 * MDP(21) + (t197 * t262 + t287 * t209 + t353 * t246 + t313 * t305 + t308 * t387) * MDP(22) + (-t198 * t262 + t287 * t208 + t353 * t248 - t305 * t387 + t313 * t308) * MDP(23) + (-MDP(4) * t306 * t309 + MDP(5) * t361) * t312; (t237 * t265 + t238 * t262 + t327) * MDP(12) + t282 * MDP(14) + (t253 + t280) * MDP(15) + (-t226 * t265 - t229 * t262 + t327 + t383 + t388) * MDP(16) + (t371 - t372) * MDP(22) + (-t228 + t370) * MDP(23) + (MDP(23) * t391 - t308 * t336) * t392 + (-t265 * MDP(14) + (-MDP(14) * t273 - MDP(15) * t366) * qJD(1)) * qJD(2) + (-MDP(14) * t368 - t273 * MDP(15) + (-MDP(12) - MDP(16)) * t291) * qJDD(1) + (MDP(11) + MDP(13)) * (-t262 ^ 2 - t261); (t253 + t240) * MDP(13) + (-t262 * t265 + qJDD(2)) * MDP(14) + (-t261 - t311) * MDP(15) + (qJD(2) * t229 + t314 - t375) * MDP(16) + (-qJD(2) * t246 + t228) * MDP(22) + (-qJD(2) * t248 - t372) * MDP(23) + (-MDP(23) * t339 - t305 * t336) * t392; t248 * t246 * MDP(17) + (-t246 ^ 2 + t248 ^ 2) * MDP(18) + (t345 + t393) * MDP(19) + (t231 + t338) * MDP(20) + t235 * MDP(21) + (-g(1) * t255 - g(2) * t257 + t198 * t392 - t216 * t248 + t199) * MDP(22) + (g(1) * t256 - g(2) * t258 + t197 * t392 + t216 * t246) * MDP(23) + (-MDP(20) * t350 + MDP(22) * t330 - MDP(23) * t337) * t308 + (-MDP(19) * t350 + (-qJD(5) * t262 - qJDD(2)) * MDP(20) - t337 * MDP(22) + (-t200 - t330) * MDP(23)) * t305;];
tau = t1;
