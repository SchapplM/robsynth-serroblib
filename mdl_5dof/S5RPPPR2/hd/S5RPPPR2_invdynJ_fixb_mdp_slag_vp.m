% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:23:08
% EndTime: 2020-01-03 11:23:14
% DurationCPUTime: 3.37s
% Computational Cost: add. (1461->339), mult. (3651->498), div. (0->0), fcn. (2888->10), ass. (0->161)
t432 = qJDD(1) * pkin(1);
t329 = qJDD(2) - t432;
t348 = sin(qJ(1));
t350 = cos(qJ(1));
t435 = -g(2) * t350 - g(3) * t348;
t438 = t435 - t329;
t341 = sin(pkin(9));
t344 = cos(pkin(9));
t346 = cos(pkin(7));
t343 = sin(pkin(7));
t345 = cos(pkin(8));
t423 = t343 * t345;
t289 = t341 * t423 + t344 * t346;
t281 = t289 * qJD(1);
t276 = qJD(5) + t281;
t342 = sin(pkin(8));
t418 = t346 * t350;
t296 = t342 * t348 + t345 * t418;
t421 = t343 * t350;
t263 = t296 * t344 + t341 * t421;
t415 = t348 * t345;
t295 = t342 * t418 - t415;
t347 = sin(qJ(5));
t349 = cos(qJ(5));
t437 = t263 * t347 - t295 * t349;
t436 = t263 * t349 + t295 * t347;
t399 = qJ(2) * qJDD(1);
t434 = g(2) * t348;
t433 = g(3) * t350;
t408 = qJD(1) * t343;
t385 = t344 * t408;
t406 = qJD(1) * t346;
t386 = t341 * t406;
t284 = t345 * t385 - t386;
t431 = t284 * t347;
t427 = t341 * t342;
t426 = t342 * t343;
t425 = t342 * t347;
t424 = t342 * t349;
t422 = t343 * t348;
t420 = t345 * t346;
t419 = t346 * t348;
t351 = qJD(1) ^ 2;
t417 = t346 * t351;
t279 = t289 * qJDD(1);
t275 = qJDD(5) + t279;
t416 = t347 * t275;
t414 = t349 * t275;
t300 = -pkin(2) * t346 - qJ(3) * t343 - pkin(1);
t404 = qJD(3) * t343;
t270 = -qJD(1) * t404 + qJDD(1) * t300 + qJDD(2);
t393 = qJDD(1) * t346;
t379 = t345 * t393;
t398 = qJD(1) * qJD(2);
t384 = t346 * t398;
t243 = qJ(2) * t379 + t270 * t342 + t345 * t384;
t397 = qJD(1) * qJD(4);
t236 = (-qJ(4) * qJDD(1) - t397) * t346 + t243;
t395 = qJDD(1) * t343;
t297 = qJ(2) * t395 + t343 * t398 + qJDD(3);
t369 = pkin(3) * t342 - qJ(4) * t345;
t247 = (qJDD(1) * t369 - t345 * t397) * t343 + t297;
t222 = t236 * t344 + t247 * t341;
t288 = qJD(1) * t300 + qJD(2);
t390 = qJ(2) * t406;
t258 = t288 * t342 + t345 * t390;
t250 = -qJ(4) * t406 + t258;
t311 = qJ(2) * t408 + qJD(3);
t268 = t369 * t408 + t311;
t229 = t250 * t344 + t268 * t341;
t273 = qJ(2) * t420 + t300 * t342;
t265 = -qJ(4) * t346 + t273;
t277 = (qJ(2) + t369) * t343;
t238 = t265 * t344 + t277 * t341;
t413 = g(1) * t346 + g(3) * t421;
t412 = pkin(1) * t350 + qJ(2) * t348;
t338 = t343 ^ 2;
t410 = t346 ^ 2 + t338;
t409 = qJD(1) * t342;
t407 = qJD(1) * t345;
t405 = qJD(2) * t346;
t403 = qJD(5) * t276;
t402 = qJD(5) * t347;
t401 = t275 * MDP(20);
t396 = qJDD(1) * t342;
t394 = qJDD(1) * t345;
t312 = t341 * t393;
t380 = t344 * t394;
t280 = t343 * t380 - t312;
t387 = t349 * t409;
t372 = t343 * t387;
t382 = t342 * t395;
t391 = qJD(5) * t372 + t280 * t349 + t347 * t382;
t389 = t342 * t408;
t388 = t347 * t409;
t337 = t342 ^ 2;
t383 = t337 * t395;
t381 = t342 * t393;
t378 = t410 * t351;
t220 = pkin(6) * t382 + t222;
t242 = -qJ(2) * t381 + t270 * t345 - t342 * t384;
t239 = pkin(3) * t393 + qJDD(4) - t242;
t224 = pkin(4) * t279 - pkin(6) * t280 + t239;
t377 = -t347 * t220 + t224 * t349;
t376 = t280 * t347 - t349 * t382;
t257 = t288 * t345 - t342 * t390;
t272 = -qJ(2) * t342 * t346 + t300 * t345;
t375 = (-t345 ^ 2 - t337) * MDP(10);
t374 = pkin(2) * t418 + qJ(3) * t421 + t412;
t373 = 0.2e1 * t410;
t293 = t342 * t419 + t345 * t350;
t371 = g(2) * t295 + g(3) * t293;
t267 = pkin(3) * t346 - t272;
t292 = -t342 * t404 + t345 * t405;
t368 = t349 * t220 + t347 * t224;
t249 = pkin(3) * t406 + qJD(4) - t257;
t225 = pkin(4) * t281 - pkin(6) * t284 + t249;
t227 = pkin(6) * t389 + t229;
t367 = t225 * t349 - t227 * t347;
t366 = -t225 * t347 - t227 * t349;
t290 = -t341 * t346 + t344 * t423;
t232 = pkin(4) * t289 - pkin(6) * t290 + t267;
t234 = pkin(6) * t426 + t238;
t365 = t232 * t349 - t234 * t347;
t364 = t232 * t347 + t234 * t349;
t221 = -t236 * t341 + t247 * t344;
t228 = -t250 * t341 + t268 * t344;
t237 = -t265 * t341 + t277 * t344;
t363 = (-qJD(4) * t345 + qJD(2)) * t343;
t333 = t348 * pkin(1);
t362 = pkin(2) * t419 - qJ(2) * t350 + qJ(3) * t422 + t333;
t361 = -t290 * t347 + t343 * t424;
t260 = t290 * t349 + t343 * t425;
t360 = t344 * t424 - t345 * t347;
t359 = t344 * t425 + t345 * t349;
t291 = t342 * t405 + t345 * t404;
t358 = qJD(1) * t291 - qJDD(1) * t272 - t242;
t357 = qJD(1) * t292 + qJDD(1) * t273 + t243;
t253 = t284 * t349 + t343 * t388;
t356 = t360 * t276;
t355 = t432 + t438;
t354 = t373 * t398 + t433;
t353 = t297 * t343 + (t398 + t399) * t338;
t294 = -t342 * t350 + t346 * t415;
t286 = (t341 * t343 + t344 * t420) * qJD(1);
t283 = t345 * t386 - t385;
t274 = -qJD(4) * t346 + t292;
t262 = t294 * t344 + t341 * t422;
t255 = t260 * qJD(5);
t254 = t361 * qJD(5);
t251 = -t372 + t431;
t246 = t344 * t274 + t341 * t363;
t245 = t274 * t341 - t344 * t363;
t241 = t262 * t349 + t293 * t347;
t240 = -t262 * t347 + t293 * t349;
t233 = -pkin(4) * t426 - t237;
t231 = qJD(5) * t253 + t376;
t230 = -t284 * t402 + t391;
t226 = -pkin(4) * t389 - t228;
t219 = -pkin(4) * t382 - t221;
t1 = [qJDD(1) * MDP(1) + t435 * MDP(2) + (-t433 + t434) * MDP(3) + t355 * t346 * MDP(4) + (t373 * t399 + t354 - t434) * MDP(6) + (-t329 * pkin(1) - g(2) * t412 - g(3) * t333 + (t399 * t410 + t354) * qJ(2)) * MDP(7) + (-g(2) * t296 - g(3) * t294 + t342 * t353 + t346 * t358) * MDP(8) + (t345 * t353 + t346 * t357 + t371) * MDP(9) + (-g(2) * t374 - g(3) * t362 + t242 * t272 + t243 * t273 - t257 * t291 + t258 * t292) * MDP(11) + (-g(2) * t263 - g(3) * t262 + t239 * t289 + t267 * t279 + t281 * t291 + (-qJD(1) * t245 + qJDD(1) * t237 + t221) * t426) * MDP(12) + (t291 * t284 + t267 * t280 + t239 * t290 - g(2) * (-t296 * t341 + t344 * t421) - g(3) * (-t294 * t341 + t344 * t422) + (-qJD(1) * t246 - qJDD(1) * t238 - t222) * t426) * MDP(13) + (-t221 * t290 - t222 * t289 - t237 * t280 - t238 * t279 + t245 * t284 - t246 * t281 - t371) * MDP(14) + (t222 * t238 + t229 * t246 + t221 * t237 - t228 * t245 + t239 * t267 + t249 * t291 - g(2) * (pkin(3) * t296 + qJ(4) * t295 + t374) - g(3) * (pkin(3) * t294 + qJ(4) * t293 + t362)) * MDP(15) + (t230 * t260 + t253 * t254) * MDP(16) + (t230 * t361 - t231 * t260 - t251 * t254 - t253 * t255) * MDP(17) + (t230 * t289 + t254 * t276 + t260 * t275) * MDP(18) + (-t231 * t289 - t255 * t276 + t275 * t361) * MDP(19) + t289 * t401 + ((-t246 * t347 + t291 * t349) * t276 + t365 * t275 + t377 * t289 + t245 * t251 + t233 * t231 - t219 * t361 + t226 * t255 - g(2) * t436 - g(3) * t241 + (-t276 * t364 + t289 * t366) * qJD(5)) * MDP(21) + (-(t246 * t349 + t291 * t347) * t276 - t364 * t275 - t368 * t289 + t245 * t253 + t233 * t230 + t219 * t260 + t226 * t254 + g(2) * t437 - g(3) * t240 + (-t276 * t365 - t289 * t367) * qJD(5)) * MDP(22) + (-t355 * MDP(5) + (-t342 * t357 + t345 * t358 + t435) * MDP(10) + (qJ(2) * t297 + qJD(2) * t311) * MDP(11)) * t343; -MDP(4) * t393 - MDP(6) * t378 + (-qJ(2) * t378 - t438) * MDP(7) + (-t342 * t378 - t379) * MDP(8) + (-t345 * t378 + t381) * MDP(9) + (t242 * t345 + t243 * t342 + (-t311 * t343 + (t257 * t342 - t258 * t345) * t346) * qJD(1) - t435) * MDP(11) + (-t341 * t383 - t279 * t345 + (-t281 * t346 + t283 * t343) * t409) * MDP(12) + (-t344 * t383 - t280 * t345 + (-t284 * t346 + t286 * t343) * t409) * MDP(13) + (t281 * t286 - t283 * t284 + (-t279 * t344 + t280 * t341) * t342) * MDP(14) + (t228 * t283 - t229 * t286 - t239 * t345 + (-t221 * t341 + t222 * t344 - t249 * t406) * t342 - t435) * MDP(15) + (-t359 * t275 + t231 * t427 - (-t286 * t347 + t346 * t387) * t276 - t283 * t251 - qJD(5) * t356) * MDP(21) + (-t360 * t275 + t230 * t427 + (t286 * t349 + t346 * t388) * t276 - t283 * t253 + t359 * t403) * MDP(22) + (MDP(5) + t375) * t395; (t297 + t413) * MDP(11) + (-t279 * t341 - t280 * t344 + (-t281 * t344 + t284 * t341) * t389) * MDP(14) + (t221 * t344 + t222 * t341 + t413) * MDP(15) + (-t344 * t231 + (-t349 * t403 - t416) * t341 + (t251 * t427 - t276 * t359) * t408) * MDP(21) + (-t344 * t230 + (t276 * t402 - t414) * t341 + (t253 * t427 - t356) * t408) * MDP(22) + ((-t345 * t417 + t396) * MDP(8) + (t342 * t417 + t394) * MDP(9) + (-t434 + (t257 * t345 + t258 * t342) * qJD(1)) * MDP(11) + (-t281 * t407 + t344 * t396) * MDP(12) + (-t284 * t407 - t341 * t396) * MDP(13) + (-t434 + (-t249 * t345 + (-t228 * t341 + t229 * t344) * t342) * qJD(1)) * MDP(15)) * t343 + (t375 + (-MDP(12) * t341 - MDP(13) * t344) * t337) * t338 * t351; t344 * MDP(12) * t393 - t312 * MDP(13) + (-t281 ^ 2 - t284 ^ 2) * MDP(14) + (-g(2) * t293 + g(3) * t295 + t228 * t284 + t229 * t281 + t239) * MDP(15) + (-t251 * t284 + t414) * MDP(21) + (-t253 * t284 - t416) * MDP(22) + ((t284 * t409 + t341 * t394) * MDP(12) + (-t281 * t409 + t380) * MDP(13) - g(1) * t342 * MDP(15)) * t343 - (MDP(21) * t347 + MDP(22) * t349) * t276 ^ 2; t253 * t251 * MDP(16) + (-t251 ^ 2 + t253 ^ 2) * MDP(17) + (t251 * t276 + t391) * MDP(18) + (t253 * t276 - t376) * MDP(19) + t401 + (-g(1) * t361 - g(2) * t240 - g(3) * t437 - t226 * t253 - t366 * t276 + t377) * MDP(21) + (g(1) * t260 + g(2) * t241 - g(3) * t436 + t226 * t251 + t367 * t276 - t368) * MDP(22) + (-MDP(18) * t431 - MDP(19) * t253 + MDP(21) * t366 - MDP(22) * t367) * qJD(5);];
tau = t1;
