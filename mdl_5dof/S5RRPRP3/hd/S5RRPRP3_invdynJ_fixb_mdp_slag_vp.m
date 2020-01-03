% Calculate vector of inverse dynamics joint torques for
% S5RRPRP3
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:19
% EndTime: 2019-12-31 19:51:23
% DurationCPUTime: 2.27s
% Computational Cost: add. (2297->295), mult. (3329->334), div. (0->0), fcn. (2162->12), ass. (0->148)
t335 = qJDD(1) + qJDD(2);
t339 = qJD(1) + qJD(2);
t345 = sin(qJ(2));
t401 = qJDD(1) * t345;
t347 = cos(qJ(2));
t405 = qJD(2) * t347;
t273 = qJ(3) * t335 + qJD(3) * t339 + (qJD(1) * t405 + t401) * pkin(1);
t341 = sin(pkin(8));
t342 = cos(pkin(8));
t408 = t341 ^ 2 + t342 ^ 2;
t381 = t408 * t273;
t406 = qJD(2) * t345;
t397 = pkin(1) * t406;
t426 = pkin(1) * t347;
t409 = -qJD(1) * t397 + qJDD(1) * t426;
t387 = qJDD(3) - t409;
t425 = pkin(2) * t335;
t281 = t387 - t425;
t340 = qJ(1) + qJ(2);
t331 = cos(t340);
t322 = g(2) * t331;
t439 = t281 + t322;
t429 = cos(qJ(4));
t391 = t429 * t342;
t344 = sin(qJ(4));
t417 = t344 * t341;
t360 = t391 - t417;
t291 = t429 * t341 + t344 * t342;
t407 = qJD(1) * t347;
t371 = -pkin(1) * t407 + qJD(3);
t284 = t291 * t339;
t403 = qJD(4) * t344;
t389 = t341 * t403;
t388 = qJD(4) * t429;
t374 = t342 * t388;
t393 = t291 * t335 + t339 * t374;
t247 = t339 * t389 - t393;
t287 = t291 * qJD(4);
t369 = t360 * t335;
t248 = t287 * t339 - t369;
t319 = pkin(3) * t342 + pkin(2);
t265 = -t319 * t335 + t387;
t350 = pkin(4) * t248 + qJ(5) * t247 + t265;
t227 = -qJD(5) * t284 + t350;
t280 = -t319 * t339 + t371;
t396 = t339 * t417;
t282 = -t339 * t391 + t396;
t237 = pkin(4) * t282 - qJ(5) * t284 + t280;
t438 = -t227 * t360 + t237 * t287;
t286 = -t374 + t389;
t437 = -t227 * t291 + t237 * t286;
t436 = -t265 * t360 + t280 * t287;
t435 = t265 * t291 - t280 * t286;
t330 = sin(t340);
t434 = g(1) * t331 + g(2) * t330;
t323 = g(1) * t330;
t433 = t322 - t323;
t338 = pkin(8) + qJ(4);
t328 = sin(t338);
t329 = cos(t338);
t385 = pkin(7) * t335 + t273;
t256 = t385 * t341;
t257 = t385 * t342;
t428 = pkin(1) * t345;
t399 = qJD(1) * t428;
t293 = qJ(3) * t339 + t399;
t384 = pkin(7) * t339 + t293;
t274 = t384 * t341;
t395 = -t344 * t256 + t429 * t257 - t274 * t388;
t432 = -g(3) * t328 - t434 * t329 + t395;
t343 = -pkin(7) - qJ(3);
t304 = t343 * t341;
t332 = t342 * pkin(7);
t305 = qJ(3) * t342 + t332;
t272 = t344 * t304 + t429 * t305;
t421 = t328 * t331;
t422 = t328 * t330;
t412 = g(1) * t422 - g(2) * t421;
t361 = t429 * t304 - t344 * t305;
t415 = t361 * qJD(4) + t360 * t371;
t431 = t415 * qJD(4) + qJDD(4) * t272 + t412;
t430 = t284 ^ 2;
t346 = sin(qJ(1));
t427 = pkin(1) * t346;
t424 = qJDD(4) * pkin(4);
t423 = t282 * t284;
t420 = t329 * t331;
t419 = t331 * t343;
t275 = t384 * t342;
t418 = t344 * t275;
t414 = t272 * qJD(4) + t291 * t371;
t413 = t439 * t341;
t411 = t331 * pkin(2) + t330 * qJ(3);
t242 = -t344 * t274 + t429 * t275;
t404 = qJD(4) * t242;
t241 = -t429 * t274 - t418;
t402 = qJD(5) - t241;
t400 = qJDD(4) * qJ(5);
t394 = pkin(4) * t420 + qJ(5) * t421 + t331 * t319;
t390 = t339 * t406;
t386 = -pkin(2) * t330 + t331 * qJ(3);
t312 = pkin(1) * t405 + qJD(3);
t382 = t312 * t408;
t380 = t408 * t335;
t379 = t241 + t418;
t378 = t339 * t399;
t377 = t429 * t256 + t344 * t257 - t274 * t403 + t275 * t388;
t376 = -t434 + t381;
t375 = -t409 + t433;
t373 = -g(2) * t420 + t329 * t323;
t243 = pkin(4) * t287 + qJ(5) * t286 - qJD(5) * t291;
t372 = -t243 + t399;
t368 = pkin(4) * t329 + qJ(5) * t328;
t365 = (-t247 * t360 - t248 * t291 + t282 * t286 - t284 * t287) * MDP(12) + (-t247 * t291 - t284 * t286) * MDP(11) + (-qJD(4) * t286 + qJDD(4) * t291) * MDP(13) + (-qJD(4) * t287 + qJDD(4) * t360) * MDP(14) + t335 * MDP(4);
t225 = t400 + (qJD(5) - t418) * qJD(4) + t395;
t226 = qJDD(5) + t377 - t424;
t238 = -qJD(4) * pkin(4) + t402;
t239 = qJD(4) * qJ(5) + t242;
t363 = t225 * t360 + t226 * t291 - t238 * t286 - t239 * t287 - t434;
t318 = qJ(3) + t428;
t288 = (-pkin(7) - t318) * t341;
t289 = t318 * t342 + t332;
t362 = t429 * t288 - t344 * t289;
t261 = t344 * t288 + t429 * t289;
t359 = -t378 - t425;
t325 = -pkin(2) - t426;
t358 = pkin(1) * t390 + t325 * t335;
t233 = t362 * qJD(4) + t312 * t360;
t357 = qJD(4) * t233 + qJDD(4) * t261 + t412;
t264 = -pkin(4) * t360 - qJ(5) * t291 - t319;
t356 = g(1) * t421 + g(2) * t422 - g(3) * t329 - t377;
t234 = t261 * qJD(4) + t291 * t312;
t355 = -qJD(4) * t234 + qJDD(4) * t362 + t373;
t354 = t371 * t408;
t353 = -t414 * qJD(4) + qJDD(4) * t361 + t373;
t352 = (-g(1) * (-t319 - t368) + g(2) * t343) * t330;
t351 = t237 * t284 + qJDD(5) - t356;
t348 = cos(qJ(1));
t333 = t348 * pkin(1);
t311 = t342 * t323;
t303 = -t319 - t426;
t292 = -pkin(2) * t339 + t371;
t279 = t282 ^ 2;
t255 = t264 - t426;
t246 = pkin(4) * t284 + qJ(5) * t282;
t240 = t243 + t397;
t236 = (t282 - t396) * qJD(4) + t393;
t1 = [(t225 * t261 + t239 * t233 + t227 * t255 + t237 * t240 - t226 * t362 + t238 * t234 - g(1) * (-t419 - t427) - g(2) * (t333 + t394) + t352) * MDP(21) + (t318 * t380 + t339 * t382 + t376) * MDP(9) + qJDD(1) * MDP(1) + (-t233 * t282 + t234 * t284 + t247 * t362 - t248 * t261 + t363) * MDP(19) + (g(1) * t346 - g(2) * t348) * MDP(2) + (g(1) * t348 + g(2) * t346) * MDP(3) + t365 + (t281 * t325 + t292 * t397 - g(1) * (t386 - t427) - g(2) * (t333 + t411) + t293 * t382 + t318 * t381) * MDP(10) + (t248 * t303 + t282 * t397 + t355 + t436) * MDP(16) + (t240 * t282 + t248 * t255 + t355 + t438) * MDP(18) + (-t247 * t303 + t284 * t397 - t357 + t435) * MDP(17) + (-t240 * t284 + t247 * t255 + t357 + t437) * MDP(20) + (t311 + (-t358 - t439) * t342) * MDP(7) + ((t358 - t323) * t341 + t413) * MDP(8) + ((t335 * t347 - t390) * pkin(1) - t375) * MDP(5) + (((-qJDD(1) - t335) * t345 + (-qJD(1) - t339) * t405) * pkin(1) + t434) * MDP(6); (-t375 + t378) * MDP(5) + ((-t401 + (-qJD(2) + t339) * t407) * pkin(1) + t434) * MDP(6) + (t311 + (-t359 - t439) * t342) * MDP(7) + ((t359 - t323) * t341 + t413) * MDP(8) + (qJ(3) * t380 + t339 * t354 + t376) * MDP(9) + (-t281 * pkin(2) - g(1) * t386 - g(2) * t411 + qJ(3) * t381 - t292 * t399 + t354 * t293) * MDP(10) + (-t248 * t319 - t282 * t399 + t353 + t436) * MDP(16) + (t247 * t319 - t284 * t399 - t431 + t435) * MDP(17) + (t248 * t264 - t282 * t372 + t353 + t438) * MDP(18) + (t247 * t361 - t248 * t272 - t415 * t282 + t414 * t284 + t363) * MDP(19) + (t247 * t264 + t372 * t284 + t431 + t437) * MDP(20) + (g(1) * t419 - g(2) * t394 + t225 * t272 - t226 * t361 + t227 * t264 - t372 * t237 + t414 * t238 + t415 * t239 + t352) * MDP(21) + t365; -t408 * MDP(9) * t339 ^ 2 + (-t408 * t339 * t293 + t281 + t433) * MDP(10) + (-t279 - t430) * MDP(19) + (t239 * t282 + (-qJD(5) - t238) * t284 + t350 + t433) * MDP(21) + (MDP(16) + MDP(18)) * (0.2e1 * t284 * qJD(4) - t369) + (-t342 * MDP(7) + t341 * MDP(8)) * t335 + (-MDP(17) + MDP(20)) * ((t282 + t396) * qJD(4) - t393); MDP(11) * t423 + (-t279 + t430) * MDP(12) + t236 * MDP(13) + t369 * MDP(14) + qJDD(4) * MDP(15) + (-t280 * t284 + t356 + t404) * MDP(16) + (t379 * qJD(4) + t280 * t282 - t432) * MDP(17) + (-t246 * t282 - t351 + t404 + 0.2e1 * t424) * MDP(18) + (pkin(4) * t247 - qJ(5) * t248 + (t239 - t242) * t284 + (t238 - t402) * t282) * MDP(19) + (0.2e1 * t400 - t237 * t282 + t246 * t284 + (0.2e1 * qJD(5) - t379) * qJD(4) + t432) * MDP(20) + (-t226 * pkin(4) - g(3) * t368 + t225 * qJ(5) - t237 * t246 - t238 * t242 + t402 * t239 + t434 * (pkin(4) * t328 - qJ(5) * t329)) * MDP(21); (-qJDD(4) + t423) * MDP(18) + t236 * MDP(19) + (-qJD(4) ^ 2 - t430) * MDP(20) + (-qJD(4) * t239 + t351 - t424) * MDP(21);];
tau = t1;
