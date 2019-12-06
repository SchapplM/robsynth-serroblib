% Calculate vector of inverse dynamics joint torques for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:53
% EndTime: 2019-12-05 15:19:57
% DurationCPUTime: 2.40s
% Computational Cost: add. (1377->294), mult. (3589->444), div. (0->0), fcn. (3395->14), ass. (0->151)
t280 = cos(pkin(5));
t263 = qJD(1) * t280 + qJD(2);
t276 = sin(pkin(5));
t274 = sin(pkin(11));
t283 = sin(qJ(3));
t286 = cos(qJ(3));
t277 = cos(pkin(11));
t279 = cos(pkin(6));
t369 = t277 * t279;
t313 = t274 * t286 + t283 * t369;
t301 = t313 * t276;
t275 = sin(pkin(6));
t372 = t275 * t283;
t225 = qJD(1) * t301 + t263 * t372;
t261 = qJDD(1) * t280 + qJDD(2);
t365 = qJD(1) * t276;
t341 = t277 * t365;
t323 = t279 * t341;
t349 = qJDD(1) * t276;
t333 = t277 * t349;
t363 = qJD(3) * t283;
t340 = t275 * t363;
t342 = t274 * t365;
t373 = t274 * t283;
t346 = t276 * t373;
t361 = qJD(3) * t286;
t290 = -(t261 * t275 + t279 * t333) * t286 + qJDD(1) * t346 + t263 * t340 + t323 * t363 + t342 * t361;
t371 = t275 * t286;
t343 = t280 * t371;
t345 = t286 * t369;
t230 = -t276 * t345 - t343 + t346;
t383 = sin(pkin(10));
t329 = t383 * t277;
t278 = cos(pkin(10));
t368 = t278 * t280;
t241 = t274 * t368 + t329;
t330 = t383 * t274;
t242 = t277 * t278 - t280 * t330;
t303 = t274 * t278 + t280 * t329;
t331 = t276 * t383;
t388 = -t275 * t331 + t279 * t303;
t302 = -t277 * t368 + t330;
t370 = t276 * t278;
t389 = t275 * t370 + t279 * t302;
t319 = g(1) * (t242 * t283 + t286 * t388) + g(2) * (t241 * t283 + t286 * t389);
t308 = g(3) * t230 + t319;
t391 = qJD(3) * t225 - t290 + t308;
t285 = cos(qJ(4));
t351 = qJD(3) * qJD(4);
t336 = t285 * t351;
t282 = sin(qJ(4));
t348 = qJDD(3) * t282;
t390 = qJD(4) * qJD(5) + t336 + t348;
t224 = -t283 * t342 + (t263 * t275 + t323) * t286;
t314 = t345 - t373;
t209 = qJDD(3) * pkin(8) + (t261 * t283 + t263 * t361) * t275 + (qJD(1) * qJD(3) * t314 + qJDD(1) * t313) * t276;
t223 = qJD(3) * pkin(8) + t225;
t239 = t263 * t279 - t275 * t341;
t214 = t223 * t285 + t239 * t282;
t238 = t261 * t279 - t275 * t333;
t379 = t238 * t285;
t196 = -qJDD(4) * pkin(4) + qJD(4) * t214 + t209 * t282 - t379;
t362 = qJD(3) * t285;
t264 = -qJD(5) + t362;
t216 = t241 * t286 - t283 * t389;
t218 = t242 * t286 - t283 * t388;
t231 = t280 * t372 + t301;
t312 = -t275 * t276 * t277 + t279 * t280;
t219 = t231 * t282 - t285 * t312;
t232 = t275 * t302 - t279 * t370;
t233 = t275 * t303 + t279 * t331;
t309 = g(1) * (-t218 * t282 + t233 * t285) + g(2) * (-t216 * t282 + t232 * t285) - g(3) * t219;
t321 = pkin(4) * t282 - pkin(9) * t285;
t387 = (pkin(9) * qJD(5) + qJD(3) * t321) * t264 - t196 - t309;
t347 = t285 * qJDD(3);
t247 = t282 * t351 + qJDD(5) - t347;
t255 = t321 * qJD(4);
t258 = -pkin(4) * t285 - pkin(9) * t282 - pkin(3);
t386 = (t225 - t255) * t264 + t258 * t247;
t212 = qJD(4) * pkin(9) + t214;
t385 = (pkin(8) * t264 + t212) * qJD(5) + t308;
t384 = qJD(3) * pkin(3);
t281 = sin(qJ(5));
t284 = cos(qJ(5));
t324 = t281 * qJDD(4) + t284 * t390;
t354 = qJD(5) * t282;
t334 = qJD(3) * t354;
t226 = -t281 * t334 + t324;
t381 = t226 * t281;
t380 = t238 * t282;
t352 = t284 * qJD(4);
t364 = qJD(3) * t282;
t248 = t281 * t364 - t352;
t378 = t248 * t264;
t358 = qJD(4) * t281;
t250 = t284 * t364 + t358;
t377 = t250 * t264;
t376 = t250 * t284;
t374 = t264 * t285;
t272 = t282 ^ 2;
t366 = -t285 ^ 2 + t272;
t360 = qJD(4) * t248;
t359 = qJD(4) * t250;
t357 = qJD(4) * t282;
t356 = qJD(4) * t285;
t355 = qJD(5) * t281;
t353 = qJD(5) * t284;
t339 = t275 * t361;
t338 = t264 * t358;
t337 = t264 * t352;
t335 = t286 * t351;
t222 = -t224 - t384;
t327 = -qJD(3) * t222 - t209;
t221 = qJD(3) * t258 - t224;
t318 = t223 * t282 - t239 * t285;
t326 = -qJDD(4) * pkin(9) + qJD(4) * t318 - qJD(5) * t221 - t209 * t285 - t380;
t320 = g(1) * (t218 * t285 + t233 * t282) + g(2) * (t216 * t285 + t232 * t282);
t220 = t231 * t285 + t282 * t312;
t208 = t220 * t284 + t230 * t281;
t207 = -t220 * t281 + t230 * t284;
t288 = qJD(3) ^ 2;
t317 = qJDD(3) * t286 - t283 * t288;
t244 = t279 * t282 + t285 * t372;
t316 = -t244 * t281 - t284 * t371;
t315 = -t244 * t284 + t281 * t371;
t243 = -t279 * t285 + t282 * t372;
t311 = t247 * t281 - t264 * t353;
t310 = t247 * t284 + t264 * t355;
t307 = g(1) * t218 + g(2) * t216 + g(3) * t231;
t305 = -qJD(5) * t212 - t319;
t299 = t320 + t326;
t296 = -pkin(8) * qJDD(4) + (t222 + t224 - t384) * qJD(4);
t293 = qJD(5) * t258 * t264 - t307;
t211 = -qJD(4) * pkin(4) + t318;
t292 = -pkin(9) * t247 + (-t211 + t318) * t264;
t291 = -pkin(8) * t247 + qJD(4) * t211 - t224 * t264 - t326;
t287 = qJD(4) ^ 2;
t289 = 0.2e1 * qJDD(3) * pkin(3) - pkin(8) * t287 + t391;
t268 = t284 * qJDD(4);
t235 = qJD(4) * t244 + t282 * t339;
t234 = -qJD(4) * t243 + t285 * t339;
t229 = t231 * qJD(3);
t228 = (t276 * t314 + t343) * qJD(3);
t227 = t284 * t334 - t268 + (t348 + (qJD(5) + t362) * qJD(4)) * t281;
t202 = -qJD(4) * t219 + t228 * t285;
t201 = qJD(4) * t220 + t228 * t282;
t200 = qJD(3) * t255 + qJDD(3) * t258 + t290;
t199 = t284 * t200;
t198 = t212 * t284 + t221 * t281;
t197 = -t212 * t281 + t221 * t284;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t261 * t280 - g(3) + (t274 ^ 2 + t277 ^ 2) * t276 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(4) * t201 - qJDD(4) * t219) * MDP(11) + (-qJD(4) * t202 - qJDD(4) * t220) * MDP(12) + (-(-qJD(5) * t208 - t202 * t281 + t229 * t284) * t264 + t207 * t247 + t201 * t248 + t219 * t227) * MDP(18) + ((qJD(5) * t207 + t202 * t284 + t229 * t281) * t264 - t208 * t247 + t201 * t250 + t219 * t226) * MDP(19) + (-MDP(5) * t231 + (-MDP(11) * t285 + MDP(12) * t282 - MDP(4)) * t230) * qJDD(3) + (-t229 * MDP(4) - t228 * MDP(5) + (-t229 * t285 + t230 * t357) * MDP(11) + (t229 * t282 + t230 * t356) * MDP(12)) * qJD(3); (-g(3) * t280 + (-g(1) * t383 + g(2) * t278) * t276 + t261) * MDP(2) + (-qJD(4) * t235 - qJDD(4) * t243) * MDP(11) + (-qJD(4) * t234 - qJDD(4) * t244) * MDP(12) + (-(qJD(5) * t315 - t234 * t281 + t284 * t340) * t264 + t316 * t247 + t235 * t248 + t243 * t227) * MDP(18) + ((qJD(5) * t316 + t234 * t284 + t281 * t340) * t264 + t315 * t247 + t235 * t250 + t243 * t226) * MDP(19) + (t317 * MDP(4) + (-qJDD(3) * t283 - t286 * t288) * MDP(5) + (-t282 * t335 + t285 * t317) * MDP(11) + (-t282 * t317 - t285 * t335) * MDP(12)) * t275; qJDD(3) * MDP(3) + t391 * MDP(4) + (-t261 * t372 - t313 * t349 + (-t263 * t371 - t314 * t365 + t224) * qJD(3) + t307) * MDP(5) + (qJDD(3) * t272 + 0.2e1 * t282 * t336) * MDP(6) + 0.2e1 * (t282 * t347 - t351 * t366) * MDP(7) + (qJDD(4) * t282 + t285 * t287) * MDP(8) + (qJDD(4) * t285 - t282 * t287) * MDP(9) + (t282 * t296 + t285 * t289) * MDP(11) + (-t282 * t289 + t285 * t296) * MDP(12) + (t226 * t282 * t284 + (-t281 * t354 + t285 * t352) * t250) * MDP(13) + ((-t248 * t284 - t250 * t281) * t356 + (-t381 - t227 * t284 + (t248 * t281 - t376) * qJD(5)) * t282) * MDP(14) + ((-t226 - t337) * t285 + (t310 + t359) * t282) * MDP(15) + ((t227 + t338) * t285 + (-t311 - t360) * t282) * MDP(16) + (-t247 * t285 - t264 * t357) * MDP(17) + (t386 * t284 + t293 * t281 + (pkin(8) * t360 + t291 * t281 + t284 * t385 - t199) * t285 + (t211 * t353 + t197 * qJD(4) + t196 * t281 - t224 * t248 + (t227 - t338) * pkin(8)) * t282) * MDP(18) + (-t386 * t281 + t293 * t284 + (pkin(8) * t359 + t291 * t284 + (t200 - t385) * t281) * t285 + (-t211 * t355 - t198 * qJD(4) + t196 * t284 - t224 * t250 + (t226 - t337) * pkin(8)) * t282) * MDP(19); MDP(8) * t348 + MDP(9) * t347 + qJDD(4) * MDP(10) + (t282 * t327 - t309 + t379) * MDP(11) + (g(3) * t220 + t327 * t285 + t320 - t380) * MDP(12) + (-t264 * t376 + t381) * MDP(13) + ((t226 + t378) * t284 + (-t227 + t377) * t281) * MDP(14) + ((-t250 * t282 + t284 * t374) * qJD(3) + t311) * MDP(15) + ((t248 * t282 - t281 * t374) * qJD(3) + t310) * MDP(16) + t264 * MDP(17) * t364 + (-pkin(4) * t227 - t197 * t364 - t214 * t248 + t292 * t281 + t284 * t387) * MDP(18) + (-pkin(4) * t226 + t198 * t364 - t214 * t250 - t281 * t387 + t292 * t284) * MDP(19) + (-MDP(6) * t282 * t285 + MDP(7) * t366) * t288; t250 * t248 * MDP(13) + (-t248 ^ 2 + t250 ^ 2) * MDP(14) + (t324 - t378) * MDP(15) + (t268 - t377) * MDP(16) + t247 * MDP(17) + (-g(3) * t207 - t198 * t264 - t211 * t250 + t199) * MDP(18) + (g(3) * t208 - t197 * t264 + t211 * t248) * MDP(19) + (-MDP(16) * t334 + MDP(18) * t305 + MDP(19) * t299) * t284 + (-MDP(15) * t334 - t390 * MDP(16) + t299 * MDP(18) + (-t200 - t305) * MDP(19)) * t281;];
tau = t1;
