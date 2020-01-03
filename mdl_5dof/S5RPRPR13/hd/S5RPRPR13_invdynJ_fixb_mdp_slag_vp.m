% Calculate vector of inverse dynamics joint torques for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:33:00
% EndTime: 2019-12-31 18:33:05
% DurationCPUTime: 3.39s
% Computational Cost: add. (1778->332), mult. (4129->405), div. (0->0), fcn. (3037->10), ass. (0->152)
t307 = cos(pkin(8));
t400 = cos(qJ(3));
t358 = t400 * t307;
t345 = qJD(1) * t358;
t306 = sin(pkin(8));
t310 = sin(qJ(3));
t381 = t306 * t310;
t357 = qJD(1) * t381;
t273 = -t345 + t357;
t309 = sin(qJ(5));
t312 = cos(qJ(5));
t249 = qJD(3) * t309 - t312 * t273;
t281 = t400 * t306 + t310 * t307;
t411 = t281 * qJD(1);
t416 = qJD(5) + t411;
t419 = t249 * t416;
t251 = qJD(3) * t312 + t273 * t309;
t349 = t416 * t251;
t372 = qJD(3) * t310;
t356 = t306 * t372;
t353 = qJDD(1) * t400;
t362 = qJDD(1) * t310;
t359 = qJD(3) * t345 + t306 * t353 + t307 * t362;
t243 = qJD(1) * t356 - t359;
t241 = -qJDD(5) + t243;
t235 = t312 * t241;
t351 = t309 * t416;
t418 = -t351 * t416 - t235;
t391 = qJDD(1) * pkin(1);
t311 = sin(qJ(1));
t313 = cos(qJ(1));
t413 = g(1) * t311 - g(2) * t313;
t332 = -qJDD(2) + t391 + t413;
t305 = pkin(8) + qJ(3);
t299 = sin(t305);
t300 = cos(t305);
t344 = g(1) * t313 + g(2) * t311;
t323 = -g(3) * t299 - t344 * t300;
t364 = qJD(1) * qJD(2);
t393 = pkin(6) + qJ(2);
t404 = t393 * qJDD(1) + t364;
t256 = t404 * t306;
t257 = t404 * t307;
t286 = t393 * t306;
t282 = qJD(1) * t286;
t287 = t393 * t307;
t283 = qJD(1) * t287;
t355 = qJD(3) * t400;
t347 = t310 * t256 - t400 * t257 + t282 * t355 + t283 * t372;
t417 = t323 - t347;
t375 = -t400 * t282 - t310 * t283;
t410 = qJD(4) - t375;
t412 = t307 * MDP(4) - t306 * MDP(5);
t409 = qJ(2) * qJDD(1);
t246 = -t310 * t282 + t400 * t283;
t239 = -qJD(3) * qJ(4) - t246;
t398 = pkin(4) * t273;
t224 = -t239 - t398;
t401 = pkin(3) + pkin(7);
t408 = t401 * t241 + (t224 - t246 + t398) * t416;
t295 = g(3) * t300;
t407 = -t344 * t299 + t295;
t231 = (qJD(2) * t306 + qJD(3) * t287) * t310 - qJD(2) * t358 + t286 * t355;
t248 = -t310 * t286 + t400 * t287;
t406 = qJD(3) * t231 - qJDD(3) * t248 - t299 * t413;
t232 = t281 * qJD(2) + t248 * qJD(3);
t247 = t400 * t286 + t310 * t287;
t405 = -qJD(3) * t232 - qJDD(3) * t247 + t300 * t413;
t403 = t273 ^ 2;
t402 = t411 ^ 2;
t278 = t281 * qJD(3);
t340 = t306 * t362 - t307 * t353;
t244 = qJD(1) * t278 + t340;
t399 = pkin(3) * t244;
t392 = qJ(4) * t273;
t390 = qJDD(3) * pkin(3);
t369 = qJD(5) * t312;
t360 = t312 * qJDD(3) + t309 * t244 + t273 * t369;
t363 = qJD(3) * qJD(5);
t217 = -t309 * t363 + t360;
t389 = t217 * t312;
t280 = -t358 + t381;
t297 = pkin(2) * t307 + pkin(1);
t333 = -qJ(4) * t281 - t297;
t227 = t401 * t280 + t333;
t388 = t227 * t241;
t387 = t241 * t309;
t386 = t249 * t273;
t385 = t251 * t273;
t384 = t273 * t411;
t383 = t280 * t309;
t379 = t309 * t311;
t378 = t309 * t313;
t377 = t311 * t312;
t376 = t312 * t313;
t373 = t306 ^ 2 + t307 ^ 2;
t285 = -t297 * qJD(1) + qJD(2);
t322 = -qJ(4) * t411 + t285;
t220 = t401 * t273 + t322;
t371 = qJD(5) * t220;
t370 = qJD(5) * t280;
t368 = t246 * qJD(3);
t366 = pkin(4) * t411 + t410;
t361 = qJDD(3) * qJ(4);
t352 = t373 * qJD(1) ^ 2;
t350 = t312 * t416;
t284 = -t297 * qJDD(1) + qJDD(2);
t321 = qJ(4) * t243 + t284;
t319 = -qJD(4) * t411 + t321;
t208 = t401 * t244 + t319;
t223 = -t401 * qJD(3) + t366;
t348 = qJD(5) * t223 + t208;
t346 = 0.2e1 * t373;
t342 = pkin(3) * t300 + qJ(4) * t299;
t339 = -t371 + t295;
t210 = t220 * t312 + t223 * t309;
t277 = -t307 * t355 + t356;
t337 = qJ(4) * t277 - qJD(4) * t281;
t331 = t297 + t342;
t330 = t278 * t309 + t280 * t369;
t329 = t400 * t256 + t310 * t257 - t282 * t372 + t283 * t355;
t215 = -qJD(3) * qJD(4) + t347 - t361;
t213 = -pkin(4) * t244 - t215;
t233 = t281 * pkin(4) + t247;
t328 = t213 * t280 + t224 * t278 + t233 * t241;
t327 = -t350 * t416 + t387;
t325 = qJDD(4) + t329;
t320 = t346 * t364 - t344;
t318 = t213 + (t416 * t401 + t392) * t416 + t323;
t230 = pkin(3) * t273 + t322;
t317 = t230 * t411 + t325 + t407;
t268 = -t299 * t379 + t376;
t267 = t299 * t377 + t378;
t266 = t299 * t378 + t377;
t265 = t299 * t376 - t379;
t259 = qJD(3) * t273;
t242 = pkin(3) * t280 + t333;
t240 = pkin(3) * t411 + t392;
t238 = -qJD(3) * pkin(3) + t410;
t237 = t312 * t244;
t234 = -t280 * pkin(4) + t248;
t226 = pkin(3) * t278 + t337;
t222 = -t277 * pkin(4) + t232;
t221 = -pkin(4) * t278 - t231;
t219 = t401 * t278 + t337;
t218 = t251 * qJD(5) + qJDD(3) * t309 - t237;
t216 = t325 - t390;
t214 = t319 + t399;
t212 = -t243 * pkin(4) - t401 * qJDD(3) + t325;
t211 = t312 * t212;
t209 = -t220 * t309 + t223 * t312;
t1 = [(-qJD(3) * t277 + qJDD(3) * t281) * MDP(10) + (-qJD(3) * t278 - qJDD(3) * t280) * MDP(11) + (-t244 * t297 + t278 * t285 + t280 * t284 + t405) * MDP(13) + (t243 * t297 - t277 * t285 + t281 * t284 + t406) * MDP(14) + (t215 * t280 + t216 * t281 + t231 * t273 + t232 * t411 - t238 * t277 + t239 * t278 - t243 * t247 - t244 * t248 - t344) * MDP(15) + (-t214 * t280 - t226 * t273 - t230 * t278 - t242 * t244 - t405) * MDP(16) + (-t214 * t281 - t226 * t411 + t230 * t277 + t242 * t243 - t406) * MDP(17) + (t214 * t242 - t215 * t248 + t216 * t247 + t230 * t226 + t239 * t231 + t238 * t232 + (-g(1) * t393 - g(2) * t331) * t313 + (g(1) * t331 - g(2) * t393) * t311) * MDP(18) + (t217 * t383 + t330 * t251) * MDP(19) + ((-t249 * t309 + t251 * t312) * t278 + (t389 - t218 * t309 + (-t249 * t312 - t251 * t309) * qJD(5)) * t280) * MDP(20) + (t217 * t281 - t241 * t383 - t251 * t277 + t330 * t416) * MDP(21) + (-t280 * t235 - t218 * t281 + t249 * t277 + (t278 * t312 - t309 * t370) * t416) * MDP(22) + (-t241 * t281 - t277 * t416) * MDP(23) + (-g(1) * t268 - g(2) * t266 - t209 * t277 + t211 * t281 + t234 * t218 + t221 * t249 + (-t208 * t281 - t219 * t416 + t388) * t309 + (t222 * t416 - t328) * t312 + ((-t227 * t312 - t233 * t309) * t416 - t210 * t281 + t224 * t383) * qJD(5)) * MDP(24) + (g(1) * t267 - g(2) * t265 + t210 * t277 + t234 * t217 + t221 * t251 + (-(qJD(5) * t233 + t219) * t416 + t388 - t348 * t281 + t224 * t370) * t312 + (-(-qJD(5) * t227 + t222) * t416 - (t212 - t371) * t281 + t328) * t309) * MDP(25) + qJDD(1) * MDP(1) + (t346 * t409 + t320) * MDP(6) + (t332 * pkin(1) + (t373 * t409 + t320) * qJ(2)) * MDP(7) + (-t243 * t281 - t277 * t411) * MDP(8) + (t243 * t280 - t244 * t281 + t273 * t277 - t278 * t411) * MDP(9) + t413 * MDP(2) + t344 * MDP(3) + t412 * (t332 + t391); -MDP(6) * t352 + (-qJ(2) * t352 - t332) * MDP(7) + ((-t273 - t357) * qJD(3) + t359) * MDP(14) + (-t402 - t403) * MDP(15) + (t243 + t259) * MDP(17) + (t399 - t239 * t273 + (-qJD(4) - t238) * t411 + t321 - t413) * MDP(18) + (t327 + t386) * MDP(24) + (t385 - t418) * MDP(25) - t412 * qJDD(1) + (MDP(13) - MDP(16)) * (0.2e1 * t411 * qJD(3) + t340); MDP(8) * t384 + (t402 - t403) * MDP(9) + ((t273 - t357) * qJD(3) + t359) * MDP(10) - t340 * MDP(11) + qJDD(3) * MDP(12) + (-t285 * t411 - t329 + t368 - t407) * MDP(13) + (qJD(3) * t375 + t273 * t285 - t417) * MDP(14) + (pkin(3) * t243 - qJ(4) * t244 + (-t239 - t246) * t411 + (t238 - t410) * t273) * MDP(15) + (t240 * t273 + t317 - t368 - 0.2e1 * t390) * MDP(16) + (0.2e1 * t361 - t230 * t273 + t240 * t411 + (0.2e1 * qJD(4) - t375) * qJD(3) + t417) * MDP(17) + (-t216 * pkin(3) - g(3) * t342 - t215 * qJ(4) - t230 * t240 - t238 * t246 - t410 * t239 + t344 * (pkin(3) * t299 - qJ(4) * t300)) * MDP(18) + (-t309 * t349 + t389) * MDP(19) + ((-t218 - t349) * t312 + (-t217 + t419) * t309) * MDP(20) + (t385 + t418) * MDP(21) + (t327 - t386) * MDP(22) + t416 * t273 * MDP(23) + (qJ(4) * t218 + t209 * t273 + t366 * t249 + t318 * t309 + t312 * t408) * MDP(24) + (qJ(4) * t217 - t210 * t273 + t366 * t251 - t309 * t408 + t318 * t312) * MDP(25); (-t243 + t259) * MDP(15) + (qJDD(3) - t384) * MDP(16) + (-qJD(3) ^ 2 - t402) * MDP(17) + (t239 * qJD(3) + t317 - t390) * MDP(18) + (-qJD(3) * t249 - t235) * MDP(24) + (-qJD(3) * t251 + t387) * MDP(25) + (-MDP(24) * t351 - MDP(25) * t350) * t416; t251 * t249 * MDP(19) + (-t249 ^ 2 + t251 ^ 2) * MDP(20) + (t360 + t419) * MDP(21) + (t237 + t349) * MDP(22) - t241 * MDP(23) + (-g(1) * t265 - g(2) * t267 + t210 * t416 - t224 * t251 + t211) * MDP(24) + (g(1) * t266 - g(2) * t268 + t209 * t416 + t224 * t249) * MDP(25) + (-MDP(22) * t363 + t339 * MDP(24) - t348 * MDP(25)) * t312 + (-MDP(21) * t363 + (-qJD(5) * t273 - qJDD(3)) * MDP(22) - t348 * MDP(24) + (-t212 - t339) * MDP(25)) * t309;];
tau = t1;
