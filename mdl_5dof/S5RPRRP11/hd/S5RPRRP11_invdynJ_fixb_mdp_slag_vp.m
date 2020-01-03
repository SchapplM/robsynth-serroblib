% Calculate vector of inverse dynamics joint torques for
% S5RPRRP11
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:42
% EndTime: 2019-12-31 18:54:49
% DurationCPUTime: 5.20s
% Computational Cost: add. (3366->399), mult. (7831->493), div. (0->0), fcn. (5757->10), ass. (0->153)
t338 = cos(pkin(8));
t344 = cos(qJ(3));
t402 = t344 * t338;
t337 = sin(pkin(8));
t341 = sin(qJ(3));
t406 = t337 * t341;
t440 = t402 - t406;
t445 = t440 * qJD(1);
t387 = qJD(1) * qJD(2);
t418 = pkin(6) + qJ(2);
t428 = qJDD(1) * t418 + t387;
t291 = t428 * t337;
t292 = t428 * t338;
t317 = t418 * t337;
t310 = qJD(1) * t317;
t318 = t418 * t338;
t311 = qJD(1) * t318;
t280 = -t341 * t310 + t344 * t311;
t443 = qJD(3) * t280;
t367 = -t291 * t344 - t341 * t292 - t443;
t248 = -qJDD(3) * pkin(3) - t367;
t340 = sin(qJ(4));
t343 = cos(qJ(4));
t306 = t440 * qJD(3);
t357 = qJD(1) * t306;
t309 = t337 * t344 + t338 * t341;
t358 = t309 * qJDD(1);
t348 = t358 + t357;
t389 = t343 * qJD(3);
t391 = qJD(4) * t340;
t431 = t309 * qJD(1);
t252 = -qJD(4) * t389 - t340 * qJDD(3) - t343 * t348 + t391 * t431;
t390 = qJD(4) * t343;
t392 = qJD(3) * t340;
t253 = -t343 * qJDD(3) + t340 * t358 + t390 * t431 + t392 * (qJD(4) + t445);
t287 = t343 * t431 + t392;
t229 = pkin(4) * t253 + qJ(5) * t252 - qJD(5) * t287 + t248;
t336 = pkin(8) + qJ(3);
t329 = cos(t336);
t419 = g(3) * t329;
t328 = sin(t336);
t342 = sin(qJ(1));
t345 = cos(qJ(1));
t379 = g(1) * t345 + g(2) * t342;
t436 = t328 * t379;
t355 = -t419 + t436;
t444 = -t229 + t355;
t295 = qJD(4) - t445;
t373 = -t291 * t341 + t292 * t344;
t434 = -t310 * t344 - t341 * t311;
t247 = qJDD(3) * pkin(7) + qJD(3) * t434 + t373;
t307 = t309 * qJD(3);
t320 = qJDD(1) * t402;
t375 = qJDD(1) * t406 - t320;
t278 = qJD(1) * t307 + t375;
t323 = t338 * pkin(2) + pkin(1);
t433 = -pkin(7) * t309 - t323;
t249 = t278 * pkin(3) - pkin(7) * t357 + qJDD(1) * t433 + qJDD(2);
t315 = -qJD(1) * t323 + qJD(2);
t259 = -pkin(3) * t445 - pkin(7) * t431 + t315;
t274 = qJD(3) * pkin(7) + t280;
t362 = t343 * t247 + t340 * t249 + t259 * t390 - t274 * t391;
t272 = qJDD(4) + t278;
t416 = qJ(5) * t272;
t227 = qJD(5) * t295 + t362 + t416;
t238 = t259 * t343 - t274 * t340;
t388 = qJD(5) - t238;
t234 = -pkin(4) * t295 + t388;
t442 = t295 * t234 + t227;
t383 = t340 * t247 - t343 * t249 + t259 * t391 + t274 * t390;
t426 = pkin(4) * t272;
t228 = qJDD(5) + t383 - t426;
t239 = t259 * t340 + t274 * t343;
t235 = qJ(5) * t295 + t239;
t414 = t235 * t295;
t441 = -t228 + t414;
t432 = g(1) * t342 - g(2) * t345;
t438 = qJDD(2) - t432;
t435 = t432 * t328;
t277 = -pkin(3) * t440 + t433;
t284 = -t317 * t341 + t318 * t344;
t396 = t340 * t277 + t343 * t284;
t430 = qJ(2) * qJDD(1);
t273 = -qJD(3) * pkin(3) - t434;
t285 = t340 * t431 - t389;
t240 = pkin(4) * t285 - qJ(5) * t287 + t273;
t425 = pkin(7) * t272;
t429 = t240 * t295 - t425;
t420 = g(3) * t328;
t354 = -t379 * t329 - t420;
t427 = t287 ^ 2;
t417 = pkin(7) * qJD(4);
t415 = qJDD(1) * pkin(1);
t413 = t239 * t295;
t412 = t252 * t340;
t411 = t285 * t287;
t410 = t285 * t445;
t384 = t287 * t295;
t409 = t295 * t340;
t408 = t309 * t343;
t264 = t340 * t272;
t405 = t340 * t342;
t404 = t342 * t343;
t265 = t343 * t272;
t403 = t343 * t345;
t401 = t345 * t340;
t400 = -t340 * t253 - t285 * t390;
t399 = t295 * t390 + t264;
t398 = t409 * t445 + t265;
t275 = pkin(3) * t431 - pkin(7) * t445;
t397 = t340 * t275 + t343 * t434;
t376 = pkin(4) * t340 - qJ(5) * t343;
t395 = -qJD(5) * t340 + t295 * t376 - t280;
t394 = (g(1) * t403 + g(2) * t404) * t328;
t393 = t337 ^ 2 + t338 ^ 2;
t386 = t295 * t417;
t382 = 0.2e1 * t393;
t296 = t329 * t405 + t403;
t298 = t329 * t401 - t404;
t381 = -g(1) * t296 + g(2) * t298;
t297 = t329 * t404 - t401;
t299 = t329 * t403 + t405;
t380 = g(1) * t297 - g(2) * t299;
t377 = pkin(4) * t343 + qJ(5) * t340;
t374 = t234 * t343 - t235 * t340;
t371 = -t317 * t344 - t318 * t341;
t370 = pkin(3) + t377;
t369 = pkin(3) * t329 + pkin(7) * t328 + t323;
t368 = t415 - t438;
t366 = t386 + t419;
t365 = t306 * t340 + t309 * t390;
t364 = -t306 * t343 + t309 * t391;
t363 = t273 * t295 - t425;
t260 = qJD(2) * t440 + qJD(3) * t371;
t276 = pkin(3) * t307 - pkin(7) * t306;
t361 = t343 * t260 + t340 * t276 + t277 * t390 - t284 * t391;
t353 = g(1) * t298 + g(2) * t296 + t340 * t420 - t383;
t352 = t382 * t387 - t379;
t350 = t240 * t287 + qJDD(5) - t353;
t349 = -g(1) * t299 - g(2) * t297 - t343 * t420 + t362;
t261 = qJD(2) * t309 + qJD(3) * t284;
t314 = -qJDD(1) * t323 + qJDD(2);
t255 = pkin(4) * t287 + qJ(5) * t285;
t254 = t309 * t376 - t371;
t246 = pkin(4) * t440 - t277 * t343 + t284 * t340;
t245 = -qJ(5) * t440 + t396;
t237 = -pkin(4) * t431 - t275 * t343 + t340 * t434;
t236 = qJ(5) * t431 + t397;
t233 = t285 * t295 - t252;
t232 = t376 * t306 + (qJD(4) * t377 - qJD(5) * t343) * t309 + t261;
t231 = -pkin(4) * t307 + qJD(4) * t396 + t260 * t340 - t276 * t343;
t230 = qJ(5) * t307 - qJD(5) * t440 + t361;
t1 = [(qJD(3) * t306 + qJDD(3) * t309) * MDP(10) + (-qJD(3) * t307 + qJDD(3) * t440) * MDP(11) + (-qJD(3) * t261 + qJDD(3) * t371 - t278 * t323 + t307 * t315 - t314 * t440 + t329 * t432) * MDP(13) + (-t284 * qJDD(3) + t315 * t306 + t314 * t309 - t435 - t323 * t358 + (-t323 * t445 - t260) * qJD(3)) * MDP(14) + (-t252 * t408 - t287 * t364) * MDP(15) + ((-t285 * t343 - t287 * t340) * t306 + (t412 - t253 * t343 + (t285 * t340 - t287 * t343) * qJD(4)) * t309) * MDP(16) + (t252 * t440 + t265 * t309 + t287 * t307 - t295 * t364) * MDP(17) + (t253 * t440 - t264 * t309 - t285 * t307 - t295 * t365) * MDP(18) + (-t272 * t440 + t295 * t307) * MDP(19) + (t383 * t440 + t238 * t307 + t261 * t285 - t371 * t253 + ((-qJD(4) * t284 + t276) * t295 + t277 * t272 + t273 * qJD(4) * t309) * t343 + ((-qJD(4) * t277 - t260) * t295 - t284 * t272 + t248 * t309 + t273 * t306) * t340 + t380) * MDP(20) + (-t239 * t307 + t248 * t408 + t252 * t371 + t261 * t287 - t272 * t396 - t273 * t364 - t295 * t361 + t362 * t440 + t381) * MDP(21) + (t229 * t309 * t340 + t228 * t440 - t231 * t295 + t232 * t285 - t234 * t307 + t240 * t365 - t246 * t272 + t253 * t254 + t380) * MDP(22) + (-t230 * t285 + t231 * t287 - t245 * t253 - t246 * t252 + t435 + t374 * t306 + (-t227 * t340 + t228 * t343 + (-t234 * t340 - t235 * t343) * qJD(4)) * t309) * MDP(23) + (-t227 * t440 - t229 * t408 + t230 * t295 - t232 * t287 + t235 * t307 + t240 * t364 + t245 * t272 + t252 * t254 - t381) * MDP(24) + (t227 * t245 + t235 * t230 + t229 * t254 + t240 * t232 + t228 * t246 + t234 * t231 - g(1) * (-pkin(4) * t297 - qJ(5) * t296) - g(2) * (pkin(4) * t299 + qJ(5) * t298) + (-g(1) * t418 - g(2) * t369) * t345 + (g(1) * t369 - g(2) * t418) * t342) * MDP(25) + qJDD(1) * MDP(1) + (t382 * t430 + t352) * MDP(6) + (pkin(1) * t368 + (t393 * t430 + t352) * qJ(2)) * MDP(7) + (t306 * t431 + t309 * t348) * MDP(8) + (-t309 * t278 + t306 * t445 - t307 * t431 + t348 * t440) * MDP(9) + t432 * MDP(2) + t379 * MDP(3) + (t338 * MDP(4) - t337 * MDP(5)) * (t368 + t415); t438 * MDP(7) - t320 * MDP(13) + t398 * MDP(20) + t400 * MDP(23) + t399 * MDP(24) - t432 * MDP(25) + (-t240 * MDP(25) + (-MDP(21) + MDP(24)) * t287 + (-MDP(20) - MDP(22)) * t285) * t431 + (t272 * MDP(22) + (t252 + t410) * MDP(23) + t441 * MDP(25) + (-MDP(21) * t295 - MDP(24) * t445) * t295) * t343 + (-t272 * MDP(21) + t442 * MDP(25) + MDP(23) * t384 + (-qJD(4) * MDP(20) - MDP(22) * t295) * t295) * t340 + (-pkin(1) * MDP(7) + (MDP(14) * t341 - MDP(4)) * t338 + (MDP(13) * t341 + MDP(14) * t344 + MDP(5)) * t337) * qJDD(1) + (t431 * MDP(13) + t445 * MDP(14) + (MDP(13) * t309 + MDP(14) * t440) * qJD(1)) * qJD(3) + (-qJ(2) * MDP(7) - MDP(6)) * qJD(1) ^ 2 * t393; -t445 ^ 2 * MDP(9) + t358 * MDP(10) - t375 * MDP(11) + qJDD(3) * MDP(12) + (t355 + t367 + t443) * MDP(13) + (-t315 * t445 - t354 - t373) * MDP(14) + (t343 * t384 - t412) * MDP(15) + ((-t252 + t410) * t343 - t287 * t409 + t400) * MDP(16) + (-t295 * t343 * t445 + t399) * MDP(17) + (-t295 * t391 + t398) * MDP(18) + (-pkin(3) * t253 - t280 * t285 + (-t419 - t248 + (-t275 - t417) * t295) * t343 + (t295 * t434 + t363) * t340 + t394) * MDP(20) + (pkin(3) * t252 + t397 * t295 - t280 * t287 + t363 * t343 + (t248 + t366 - t436) * t340) * MDP(21) + (t237 * t295 - t253 * t370 + t395 * t285 + (-t229 - t366) * t343 + t429 * t340 + t394) * MDP(22) + (t236 * t285 - t237 * t287 + ((qJD(4) * t287 - t253) * pkin(7) + t442) * t343 + ((qJD(4) * t285 - t252) * pkin(7) - t441) * t340 + t354) * MDP(23) + (-t236 * t295 - t252 * t370 - t395 * t287 - t429 * t343 + (-t386 + t444) * t340) * MDP(24) + (-t234 * t237 - t235 * t236 + t395 * t240 + (qJD(4) * t374 + t227 * t343 + t228 * t340 + t354) * pkin(7) + t444 * t370) * MDP(25) + (-t315 * MDP(13) - t287 * MDP(17) + t285 * MDP(18) - t295 * MDP(19) - t238 * MDP(20) + t239 * MDP(21) + t234 * MDP(22) - t235 * MDP(24) - MDP(8) * t445 + MDP(9) * t431) * t431; MDP(15) * t411 + (-t285 ^ 2 + t427) * MDP(16) + t233 * MDP(17) + (-t253 + t384) * MDP(18) + t272 * MDP(19) + (-t273 * t287 + t353 + t413) * MDP(20) + (t238 * t295 + t273 * t285 - t349) * MDP(21) + (-t255 * t285 - t350 + t413 + 0.2e1 * t426) * MDP(22) + (pkin(4) * t252 - qJ(5) * t253 + (t235 - t239) * t287 + (t234 - t388) * t285) * MDP(23) + (0.2e1 * t416 - t240 * t285 + t255 * t287 + (0.2e1 * qJD(5) - t238) * t295 + t349) * MDP(24) + (t227 * qJ(5) - t228 * pkin(4) - t240 * t255 - t234 * t239 - g(1) * (-pkin(4) * t298 + qJ(5) * t299) - g(2) * (-pkin(4) * t296 + qJ(5) * t297) + t376 * t420 + t388 * t235) * MDP(25); (-qJD(3) * t431 - qJDD(4) - t375 + t411) * MDP(22) + t233 * MDP(23) + (-t295 ^ 2 - t427) * MDP(24) + (t350 - t414 - t426) * MDP(25);];
tau = t1;
