% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:37
% EndTime: 2022-01-20 10:48:41
% DurationCPUTime: 1.77s
% Computational Cost: add. (1409->236), mult. (2244->322), div. (0->0), fcn. (1448->16), ass. (0->148)
t330 = sin(qJ(5));
t331 = sin(qJ(4));
t334 = cos(qJ(5));
t335 = cos(qJ(4));
t266 = t330 * t335 + t331 * t334;
t323 = qJD(1) + qJD(2);
t252 = t266 * t323;
t336 = cos(qJ(2));
t415 = pkin(1) * t336;
t309 = qJDD(1) * t415;
t321 = qJDD(1) + qJDD(2);
t332 = sin(qJ(2));
t406 = pkin(1) * qJD(2);
t375 = qJD(1) * t406;
t253 = pkin(2) * t321 - t332 * t375 + t309;
t328 = sin(pkin(9));
t329 = cos(pkin(9));
t379 = qJDD(1) * t332;
t417 = pkin(1) * t379 + t336 * t375;
t232 = t253 * t329 - t328 * t417;
t228 = -pkin(3) * t321 - t232;
t327 = qJ(1) + qJ(2);
t312 = pkin(9) + t327;
t299 = cos(t312);
t418 = g(2) * t299 + t228;
t391 = t334 * t335;
t395 = t330 * t331;
t265 = -t391 + t395;
t315 = sin(t327);
t317 = cos(t327);
t416 = g(1) * t315 - g(2) * t317;
t407 = pkin(1) * qJD(1);
t377 = t336 * t407;
t275 = pkin(2) * t323 + t377;
t378 = t332 * t407;
t295 = t329 * t378;
t245 = t328 * t275 + t295;
t367 = t245 + (pkin(7) + pkin(8)) * t323;
t230 = t335 * qJD(3) - t367 * t331;
t322 = qJD(4) + qJD(5);
t231 = qJD(3) * t331 + t367 * t335;
t414 = pkin(2) * t329;
t413 = pkin(4) * t335;
t298 = sin(t312);
t412 = g(1) * t298;
t308 = pkin(2) + t415;
t396 = t329 * t332;
t388 = pkin(1) * t396 + t328 * t308;
t255 = pkin(7) + t388;
t409 = -pkin(8) - t255;
t300 = pkin(2) * t328 + pkin(7);
t408 = -pkin(8) - t300;
t405 = t231 * t334;
t294 = t328 * t378;
t258 = t329 * t377 - t294;
t404 = t258 * t322;
t326 = qJ(4) + qJ(5);
t314 = sin(t326);
t403 = t298 * t314;
t316 = cos(t326);
t402 = t298 * t316;
t401 = t299 * t314;
t400 = t299 * t316;
t399 = t321 * t335;
t398 = t323 * t331;
t397 = t328 * t332;
t392 = t331 * t335;
t390 = qJDD(3) - g(3);
t244 = t275 * t329 - t294;
t240 = -pkin(3) * t323 - t244;
t383 = qJD(4) * t331;
t389 = t240 * t383 + t335 * t412;
t387 = g(1) * t317 + g(2) * t315;
t324 = t331 ^ 2;
t386 = -t335 ^ 2 + t324;
t384 = qJD(4) * t323;
t382 = qJD(4) * t335;
t381 = qJD(5) * t330;
t376 = pkin(4) * t383;
t373 = t323 * t395;
t372 = t323 * t391;
t371 = t240 * t382 + t331 * t418;
t233 = t328 * t253 + t329 * t417;
t370 = -pkin(3) - t413;
t369 = t323 * t382;
t229 = pkin(7) * t321 + t233;
t368 = pkin(8) * t321 + t229;
t365 = qJD(4) * t409;
t364 = qJD(4) * t408;
t363 = -pkin(1) * t397 + t308 * t329;
t362 = qJD(1) * (-qJD(2) + t323);
t254 = -pkin(3) - t363;
t360 = t309 + t416;
t256 = t328 * t377 + t295;
t359 = -t256 + t376;
t333 = sin(qJ(1));
t337 = cos(qJ(1));
t358 = g(1) * t333 - g(2) * t337;
t357 = t265 * t321;
t227 = qJD(4) * pkin(4) + t230;
t355 = -t227 * t330 - t405;
t238 = t322 * t265;
t320 = qJDD(4) + qJDD(5);
t354 = -t238 * t322 + t266 * t320;
t242 = t409 * t331;
t318 = t335 * pkin(8);
t243 = t255 * t335 + t318;
t353 = t242 * t334 - t243 * t330;
t352 = t242 * t330 + t243 * t334;
t263 = t408 * t331;
t264 = t300 * t335 + t318;
t351 = t263 * t334 - t264 * t330;
t350 = t263 * t330 + t264 * t334;
t216 = (t323 * t383 - t399) * pkin(4) + t228;
t234 = t370 * t323 - t244;
t349 = -g(1) * t403 + g(2) * t401 + t216 * t266 - t234 * t238;
t239 = t322 * t266;
t348 = g(1) * t402 - g(2) * t400 + t216 * t265 + t234 * t239;
t257 = (t328 * t336 + t396) * t406;
t338 = qJD(4) ^ 2;
t347 = t254 * t321 + t255 * t338 + t257 * t323;
t301 = -pkin(3) - t414;
t346 = -t256 * t323 + t300 * t338 + t301 * t321;
t217 = qJD(5) * t372 + t266 * t321 - t322 * t373 + t334 * t369;
t250 = -t372 + t373;
t345 = t252 * t250 * MDP(15) + (t250 * t322 + t217) * MDP(17) - t357 * MDP(18) + (-t250 ^ 2 + t252 ^ 2) * MDP(16) + t320 * MDP(19);
t344 = g(1) * t299 + g(2) * t298 - t240 * t323 - t229;
t343 = -qJDD(4) * t300 + (t301 * t323 + t258) * qJD(4);
t218 = t239 * t323 + t357;
t222 = -t239 * t322 - t265 * t320;
t279 = qJDD(4) * t331 + t335 * t338;
t280 = qJDD(4) * t335 - t331 * t338;
t342 = (-t217 * t265 - t218 * t266 + t238 * t250 - t239 * t252) * MDP(16) + (t217 * t266 - t238 * t252) * MDP(15) + t354 * MDP(17) + t222 * MDP(18) + 0.2e1 * (t321 * t392 - t386 * t384) * MDP(9) + (t321 * t324 + 0.2e1 * t331 * t369) * MDP(8) + t279 * MDP(10) + t280 * MDP(11) + t321 * MDP(4);
t259 = (t329 * t336 - t397) * t406;
t341 = -qJD(4) * t259 - qJDD(4) * t255 + t254 * t384;
t311 = t335 * qJDD(3);
t208 = qJDD(4) * pkin(4) - t231 * qJD(4) - t368 * t331 + t311;
t340 = t231 * t381 + g(2) * t402 + g(1) * t400 + g(3) * t314 + t234 * t250 + (-t231 * t322 - t208) * t330;
t209 = t230 * qJD(4) + t331 * qJDD(3) + t368 * t335;
t339 = g(1) * t401 + g(2) * t403 - g(3) * t316 + t355 * qJD(5) + t334 * t208 - t330 * t209 - t234 * t252;
t278 = t370 - t414;
t261 = t335 * t364;
t260 = t331 * t364;
t248 = t254 - t413;
t246 = t257 + t376;
t226 = -t259 * t331 + t335 * t365;
t225 = t259 * t335 + t331 * t365;
t1 = [(t321 * t336 * MDP(5) + t358 * MDP(7) + (-qJDD(1) - t321) * MDP(6) * t332 + (MDP(5) * t332 + MDP(6) * t336) * qJD(2) * (-qJD(1) - t323)) * pkin(1) + t358 * MDP(2) + (g(1) * t337 + g(2) * t333) * MDP(3) + (t416 * pkin(2) + t232 * t363 + t233 * t388 - t244 * t257 + t245 * t259) * MDP(7) + t342 + t389 * MDP(13) + t371 * MDP(14) + (t246 * t250 + t248 * t218 + (-t352 * qJD(5) - t225 * t330 + t226 * t334) * t322 + t353 * t320 + t348) * MDP(20) + (t246 * t252 + t248 * t217 - (t353 * qJD(5) + t225 * t334 + t226 * t330) * t322 - t352 * t320 + t349) * MDP(21) + (t341 * MDP(13) + (t347 - t412) * MDP(14)) * t331 + ((-t347 - t418) * MDP(13) + t341 * MDP(14)) * t335 + t360 * MDP(5) + t387 * MDP(6) + qJDD(1) * MDP(1); (t343 * t331 + (-t346 - t418) * t335 + t389) * MDP(13) + t342 + (t343 * t335 + (t346 - t412) * t331 + t371) * MDP(14) + (t244 * t256 - t245 * t258 + (t232 * t329 + t233 * t328 + t416) * pkin(2)) * MDP(7) + ((t336 * t362 - t379) * pkin(1) + t387) * MDP(6) + (t332 * pkin(1) * t362 + t360) * MDP(5) + (t278 * t218 + (-t350 * qJD(5) - t260 * t330 + t261 * t334) * t322 + t351 * t320 + t266 * t404 + t359 * t250 + t348) * MDP(20) + (t278 * t217 - (t351 * qJD(5) + t260 * t334 + t261 * t330) * t322 - t350 * t320 - t265 * t404 + t359 * t252 + t349) * MDP(21); t280 * MDP(13) - t279 * MDP(14) + t222 * MDP(20) - t354 * MDP(21) + t390 * MDP(7); t331 * t321 * MDP(10) + MDP(11) * t399 + qJDD(4) * MDP(12) + (-g(3) * t335 + t344 * t331 + t311) * MDP(13) + (-t390 * t331 + t344 * t335) * MDP(14) + (-(-t230 * t330 - t405) * t322 + (-t250 * t398 + t320 * t334 - t322 * t381) * pkin(4) + t339) * MDP(20) + ((-qJD(5) * t227 + t230 * t322 - t209) * t334 + (-qJD(5) * t322 * t334 - t252 * t398 - t320 * t330) * pkin(4) + t340) * MDP(21) + t345 + (-MDP(8) * t392 + t386 * MDP(9)) * t323 ^ 2; (-t355 * t322 + t339) * MDP(20) + ((-t209 + (-qJD(5) + t322) * t227) * t334 + t340) * MDP(21) + t345;];
tau = t1;
