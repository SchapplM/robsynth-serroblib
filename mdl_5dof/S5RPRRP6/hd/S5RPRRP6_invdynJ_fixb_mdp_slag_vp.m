% Calculate vector of inverse dynamics joint torques for
% S5RPRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:36
% EndTime: 2021-01-15 18:08:44
% DurationCPUTime: 3.22s
% Computational Cost: add. (2092->368), mult. (4270->467), div. (0->0), fcn. (2637->10), ass. (0->159)
t314 = sin(pkin(8));
t299 = pkin(1) * t314 + pkin(6);
t287 = t299 * qJDD(1);
t429 = -qJD(2) * qJD(3) - t287;
t311 = qJ(1) + pkin(8);
t304 = sin(t311);
t305 = cos(t311);
t425 = g(1) * t305 + g(2) * t304;
t289 = t299 * qJD(1);
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t256 = qJD(2) * t321 - t318 * t289;
t428 = qJD(3) * t256;
t376 = qJD(1) * qJD(3);
t361 = t321 * t376;
t373 = qJDD(1) * t318;
t427 = qJD(3) * qJD(4) + t361 + t373;
t317 = sin(qJ(4));
t320 = cos(qJ(4));
t380 = qJD(4) * t318;
t360 = qJD(1) * t380;
t338 = (-qJDD(3) + t360) * t320;
t388 = qJD(1) * t321;
t239 = t317 * ((qJD(4) + t388) * qJD(3) + t373) + t338;
t401 = qJDD(2) - g(3);
t426 = t401 * t321;
t404 = t317 * t321;
t251 = t304 * t404 + t305 * t320;
t253 = t304 * t320 - t305 * t404;
t405 = t317 * t318;
t424 = -g(1) * t253 + g(2) * t251 + g(3) * t405;
t371 = t321 * qJDD(1);
t272 = t318 * t376 + qJDD(4) - t371;
t295 = -qJD(4) + t388;
t377 = t320 * qJD(3);
t365 = t321 * t377;
t333 = -t317 * t380 + t365;
t403 = t318 * t320;
t423 = t272 * t403 - t333 * t295;
t384 = qJD(3) * t317;
t389 = qJD(1) * t318;
t278 = t320 * t389 + t384;
t421 = t278 ^ 2;
t420 = pkin(4) * t272;
t276 = t317 * t389 - t377;
t419 = pkin(4) * t276;
t418 = pkin(4) * t317;
t413 = g(3) * t321;
t412 = qJ(5) + pkin(7);
t354 = t317 * qJDD(3) + t427 * t320;
t238 = t317 * t360 - t354;
t411 = qJ(5) * t238;
t410 = qJ(5) * t239;
t409 = t276 * t295;
t408 = t278 * t295;
t407 = t295 * t320;
t406 = t299 * t317;
t402 = t320 * t321;
t257 = t318 * qJD(2) + t321 * t289;
t248 = qJD(3) * pkin(7) + t257;
t315 = cos(pkin(8));
t300 = -pkin(1) * t315 - pkin(2);
t268 = -pkin(3) * t321 - pkin(7) * t318 + t300;
t249 = t268 * qJD(1);
t226 = -t248 * t317 + t320 * t249;
t223 = -qJ(5) * t278 + t226;
t222 = -pkin(4) * t295 + t223;
t400 = -t223 + t222;
t399 = -t239 * t403 - t276 * t365;
t352 = pkin(3) * t318 - pkin(7) * t321;
t281 = t352 * qJD(1);
t398 = t320 * t256 + t317 * t281;
t282 = t352 * qJD(3);
t379 = qJD(4) * t320;
t397 = t268 * t379 + t317 * t282;
t359 = qJD(4) * t412;
t367 = t317 * t388;
t378 = qJD(5) * t320;
t396 = qJ(5) * t367 - t317 * t359 + t378 - t398;
t266 = t320 * t281;
t339 = pkin(4) * t318 - qJ(5) * t402;
t395 = -t339 * qJD(1) - t320 * t359 - t266 + (-qJD(5) + t256) * t317;
t383 = qJD(3) * t318;
t394 = t320 * t282 + t383 * t406;
t280 = t299 * t402;
t393 = t317 * t268 + t280;
t392 = t425 * t403;
t312 = t318 ^ 2;
t391 = -t321 ^ 2 + t312;
t390 = MDP(21) * t317;
t290 = qJD(1) * t300;
t386 = qJD(3) * t276;
t385 = qJD(3) * t299;
t382 = qJD(3) * t321;
t381 = qJD(4) * t317;
t372 = qJDD(2) * t321;
t369 = -t289 * t382 + t429 * t318;
t368 = pkin(6) + t418;
t366 = t295 * t384;
t247 = -qJD(3) * pkin(3) - t256;
t358 = -qJD(5) - t419;
t237 = t247 - t358;
t364 = t237 * t379;
t363 = t295 * t381;
t362 = t318 * t379;
t357 = t238 * t321 + t278 * t383;
t231 = qJDD(3) * pkin(7) + qJDD(2) * t318 + t287 * t321 + t428;
t240 = qJD(1) * t282 + t268 * qJDD(1);
t355 = t320 * t231 + t317 * t240 - t248 * t381 + t249 * t379;
t341 = -qJDD(3) * pkin(3) - t369;
t232 = t341 - t372;
t353 = -qJD(4) * pkin(7) * t295 + t232;
t351 = -g(1) * t251 - g(2) * t253;
t252 = -t304 * t402 + t305 * t317;
t254 = t304 * t317 + t305 * t402;
t350 = -g(1) * t252 - g(2) * t254;
t348 = g(1) * t304 - g(2) * t305;
t319 = sin(qJ(1));
t322 = cos(qJ(1));
t347 = g(1) * t319 - g(2) * t322;
t227 = t248 * t320 + t249 * t317;
t235 = t320 * t240;
t328 = -t227 * qJD(4) - t317 * t231 + t235;
t215 = -qJD(5) * t278 + t328 + t411 + t420;
t216 = -qJD(5) * t276 + t355 - t410;
t346 = -t215 * t317 + t216 * t320;
t224 = -qJ(5) * t276 + t227;
t345 = t222 * t320 + t224 * t317;
t344 = t222 * t317 - t224 * t320;
t303 = pkin(4) * t320 + pkin(3);
t343 = t303 * t321 + t318 * t412;
t340 = t347 * pkin(1);
t337 = pkin(2) + t343;
t336 = t425 * t318;
t335 = -t272 * t317 + t295 * t379;
t334 = -qJD(1) * t290 + t425;
t332 = -pkin(7) * t272 - t295 * t247;
t331 = 0.2e1 * t290 * qJD(3) - qJDD(3) * t299;
t330 = pkin(4) * t239 + qJDD(5) + t341;
t329 = g(1) * t254 - g(2) * t252 + g(3) * t403 - t355;
t323 = qJD(3) ^ 2;
t327 = -0.2e1 * qJDD(1) * t300 - t299 * t323 + t348;
t326 = t328 + t424;
t298 = g(3) * t404;
t292 = t412 * t320;
t291 = t412 * t317;
t286 = qJDD(3) * t321 - t318 * t323;
t285 = qJDD(3) * t318 + t321 * t323;
t271 = t276 ^ 2;
t263 = (t299 + t418) * t318;
t259 = t320 * t268;
t242 = pkin(4) * t367 + t257;
t241 = t299 * t382 + (t317 * t382 + t362) * pkin(4);
t236 = -qJ(5) * t405 + t393;
t233 = -qJ(5) * t403 + t259 + (-pkin(4) - t406) * t321;
t221 = t330 - t372;
t220 = (-qJ(5) * qJD(4) - t385) * t403 + (-qJD(5) * t318 + (-qJ(5) * qJD(3) - qJD(4) * t299) * t321) * t317 + t397;
t219 = -t318 * t378 + t339 * qJD(3) + (-t280 + (qJ(5) * t318 - t268) * t317) * qJD(4) + t394;
t1 = [qJDD(1) * MDP(1) + t347 * MDP(2) + (g(1) * t322 + g(2) * t319) * MDP(3) + ((t314 ^ 2 + t315 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t340) * MDP(4) + (qJDD(1) * t312 + 0.2e1 * t318 * t361) * MDP(5) + 0.2e1 * (t318 * t371 - t391 * t376) * MDP(6) + t285 * MDP(7) + t286 * MDP(8) + (t331 * t318 + t327 * t321) * MDP(10) + (-t327 * t318 + t331 * t321) * MDP(11) + (-t238 * t403 + t333 * t278) * MDP(12) + (-t278 * t362 + (-t278 * t382 + (qJD(4) * t276 + t238) * t318) * t317 + t399) * MDP(13) + (t357 + t423) * MDP(14) + ((t239 + t366) * t321 + (t335 - t386) * t318) * MDP(15) + (-t272 * t321 - t295 * t383) * MDP(16) + (-(-t268 * t381 + t394) * t295 + t259 * t272 + (t276 * t385 - t235 + (t295 * t299 + t248) * t379 + (qJD(3) * t247 + qJD(4) * t249 - t272 * t299 + t231) * t317) * t321 + (t226 * qJD(3) + t232 * t317 + t299 * t239 + t247 * t379) * t318 + t350) * MDP(17) + (t397 * t295 - t393 * t272 + (-t299 * t363 + (t247 * t320 + t278 * t299) * qJD(3) + t355) * t321 + (-t247 * t381 + t232 * t320 - t299 * t238 + (-t299 * t407 - t227) * qJD(3)) * t318 + t351) * MDP(18) + (-t219 * t295 + t233 * t272 + t239 * t263 + t241 * t276 + (t237 * t384 - t215) * t321 + (qJD(3) * t222 + t221 * t317 + t364) * t318 + t350) * MDP(19) + (t220 * t295 - t236 * t272 - t238 * t263 + t241 * t278 + (t237 * t377 + t216) * t321 + (-qJD(3) * t224 + t221 * t320 - t237 * t381) * t318 + t351) * MDP(20) + (-t219 * t278 - t220 * t276 + t233 * t238 - t236 * t239 - t345 * t382 + (t344 * qJD(4) - t215 * t320 - t216 * t317 + t348) * t318) * MDP(21) + (t215 * t233 + t216 * t236 + t222 * t219 + t224 * t220 + t221 * t263 + t237 * t241 + t340 + (-g(1) * t368 - g(2) * t337) * t305 + (g(1) * t337 - g(2) * t368) * t304) * MDP(22); t401 * MDP(4) + t286 * MDP(10) - t285 * MDP(11) + t399 * MDP(21) - g(3) * MDP(22) + (-t221 * MDP(22) + (-t344 * MDP(22) + t278 * t390) * qJD(3)) * t321 + (MDP(17) + MDP(19)) * ((-t239 + t366) * t321 + (t335 + t386) * t318) + (MDP(18) + MDP(20)) * (t357 - t423) + (-t238 * t390 + (qJD(3) * t237 + t346) * MDP(22) + ((t276 * t317 + t278 * t320) * MDP(21) - t345 * MDP(22)) * qJD(4)) * t318; MDP(7) * t373 + MDP(8) * t371 + qJDD(3) * MDP(9) + (qJD(3) * t257 + t334 * t318 + t369 + t426) * MDP(10) + (t428 + (qJD(3) * t289 - t401) * t318 + (t334 + t429) * t321) * MDP(11) + (-t238 * t317 - t278 * t407) * MDP(12) + ((-t238 + t409) * t320 + (-t239 + t408) * t317) * MDP(13) + ((-t278 * t318 + t295 * t402) * qJD(1) - t335) * MDP(14) + (t363 + t272 * t320 + (t276 * t318 - t295 * t404) * qJD(1)) * MDP(15) + t295 * MDP(16) * t389 + (-t226 * t389 - pkin(3) * t239 - t257 * t276 + t266 * t295 + (-t353 - t413) * t320 + (-t256 * t295 + t332) * t317 + t392) * MDP(17) + (pkin(3) * t238 - t398 * t295 + t227 * t389 - t257 * t278 + t298 + t332 * t320 + (-t336 + t353) * t317) * MDP(18) + (-t222 * t389 - t239 * t303 - t242 * t276 - t272 * t291 + (-t221 - t413) * t320 - t395 * t295 + (-t237 * t388 + (t237 + t419) * qJD(4)) * t317 + t392) * MDP(19) + (t364 + t238 * t303 - t242 * t278 - t272 * t292 + t298 + t396 * t295 + (t224 * t318 - t237 * t402) * qJD(1) + (pkin(4) * qJD(4) * t278 + t221 - t336) * t317) * MDP(20) + (-g(3) * t318 - t238 * t291 - t239 * t292 - t395 * t278 - t396 * t276 - t345 * qJD(4) + (t345 * qJD(1) - t425) * t321 + t346) * MDP(21) + (t216 * t292 - t215 * t291 - t221 * t303 - g(3) * t343 + (pkin(4) * t381 - t242) * t237 + t396 * t224 + t395 * t222 + t425 * (t303 * t318 - t321 * t412)) * MDP(22) + (-t318 * t321 * MDP(5) + t391 * MDP(6)) * qJD(1) ^ 2; t278 * t276 * MDP(12) + (-t271 + t421) * MDP(13) + (-t238 - t409) * MDP(14) + (-t239 - t408) * MDP(15) + t272 * MDP(16) + (-t227 * t295 - t247 * t278 + t326) * MDP(17) + (-t226 * t295 + t247 * t276 + t329) * MDP(18) + (0.2e1 * t420 + t411 - t224 * t295 + (-t237 + t358) * t278 + t326) * MDP(19) + (-pkin(4) * t421 + t410 - t223 * t295 + (qJD(5) + t237) * t276 + t329) * MDP(20) + (pkin(4) * t238 - t400 * t276) * MDP(21) + (t400 * t224 + (-t237 * t278 + t215 + t424) * pkin(4)) * MDP(22); (t338 - t408) * MDP(19) + (t354 + t409) * MDP(20) + (-t271 - t421) * MDP(21) + (t222 * t278 + t224 * t276 + t330 - t336 - t426) * MDP(22) + (t427 * MDP(19) - MDP(20) * t360) * t317;];
tau = t1;
