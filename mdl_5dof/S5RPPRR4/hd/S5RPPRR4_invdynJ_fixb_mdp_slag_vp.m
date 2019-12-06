% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:45:10
% EndTime: 2019-12-05 17:45:18
% DurationCPUTime: 4.33s
% Computational Cost: add. (2093->340), mult. (5244->459), div. (0->0), fcn. (4071->14), ass. (0->171)
t368 = sin(pkin(9));
t370 = cos(pkin(9));
t369 = sin(pkin(8));
t431 = qJD(3) * t369;
t371 = cos(pkin(8));
t432 = qJD(2) * t371;
t321 = -t368 * t432 - t370 * t431;
t322 = -t368 * t431 + t370 * t432;
t373 = sin(qJ(4));
t376 = cos(qJ(4));
t333 = -pkin(2) * t371 - qJ(3) * t369 - pkin(1);
t327 = t370 * t333;
t384 = -pkin(6) * t369 * t370 + (-qJ(2) * t368 - pkin(3)) * t371;
t279 = t327 + t384;
t295 = t370 * t371 * qJ(2) + t368 * t333;
t448 = t368 * t369;
t284 = -pkin(6) * t448 + t295;
t441 = t373 * t279 + t376 * t284;
t469 = t441 * qJD(4) - t376 * t321 + t322 * t373;
t329 = t368 * t376 + t370 * t373;
t387 = qJD(1) * t329;
t301 = t369 * t387;
t375 = cos(qJ(5));
t435 = qJD(1) * t369;
t418 = t368 * t435;
t403 = t373 * t418;
t445 = t370 * t376;
t421 = t369 * t445;
t304 = qJD(1) * t421 - t403;
t372 = sin(qJ(5));
t449 = t304 * t372;
t258 = t375 * t301 + t449;
t434 = qJD(1) * t371;
t465 = qJD(4) - t434;
t346 = -qJD(5) - t465;
t452 = t258 * t346;
t324 = t329 * qJD(4);
t334 = qJDD(1) * t421;
t446 = t368 * t373;
t269 = t334 + (-qJD(1) * t324 - qJDD(1) * t446) * t369;
t422 = qJDD(1) * t371;
t349 = -qJDD(4) + t422;
t292 = -qJD(1) * t431 + t333 * qJDD(1) + qJDD(2);
t287 = t370 * t292;
t424 = qJD(1) * qJD(2);
t416 = t371 * t424;
t254 = t384 * qJDD(1) - t368 * t416 + t287;
t414 = t370 * t422;
t268 = qJ(2) * t414 + t368 * t292 + t370 * t416;
t423 = qJDD(1) * t369;
t415 = t368 * t423;
t256 = -pkin(6) * t415 + t268;
t319 = t333 * qJD(1) + qJD(2);
t309 = t370 * t319;
t271 = t384 * qJD(1) + t309;
t419 = qJ(2) * t434;
t283 = t368 * t319 + t370 * t419;
t276 = -pkin(6) * t418 + t283;
t396 = -t271 * t373 - t276 * t376;
t381 = t396 * qJD(4) + t376 * t254 - t373 * t256;
t237 = -pkin(4) * t349 - pkin(7) * t269 + t381;
t468 = t372 * t237;
t409 = t376 * t271 - t276 * t373;
t245 = -pkin(7) * t304 + t409;
t243 = pkin(4) * t465 + t245;
t246 = -pkin(7) * t301 - t396;
t454 = t246 * t375;
t397 = -t243 * t372 - t454;
t466 = t397 * qJD(5);
t428 = qJD(5) * t372;
t244 = t246 * t428;
t347 = qJ(2) * t435 + qJD(3);
t320 = pkin(3) * t418 + t347;
t274 = pkin(4) * t301 + t320;
t367 = pkin(9) + qJ(4);
t361 = qJ(5) + t367;
t355 = sin(t361);
t377 = cos(qJ(1));
t356 = cos(t361);
t374 = sin(qJ(1));
t442 = t374 * t356;
t297 = -t355 * t377 + t371 * t442;
t443 = t371 * t377;
t299 = -t355 * t374 - t356 * t443;
t459 = g(1) * t369;
t464 = -g(2) * t297 - g(3) * t299 + t274 * t258 + t356 * t459 + t244;
t331 = qJD(4) * t403;
t429 = qJD(4) * t376;
t417 = t370 * t429;
t270 = t369 * (qJD(1) * t417 + t329 * qJDD(1)) - t331;
t343 = -qJDD(5) + t349;
t332 = t343 * MDP(23);
t394 = -t301 * t372 + t375 * t304;
t463 = t258 * MDP(19) * t394 + (-t258 ^ 2 + t394 ^ 2) * MDP(20) - t332;
t453 = t394 * t346;
t425 = qJ(2) * qJDD(1);
t444 = t371 * t374;
t296 = t355 * t444 + t356 * t377;
t298 = t355 * t443 - t442;
t430 = qJD(4) * t373;
t389 = -t373 * t254 - t376 * t256 - t271 * t429 + t276 * t430;
t238 = -pkin(7) * t270 - t389;
t412 = t375 * t237 - t372 * t238;
t461 = -g(2) * t296 + g(3) * t298 - t274 * t394 + t355 * t459 + t412;
t410 = t269 * t372 + t375 * t270;
t240 = t394 * qJD(5) + t410;
t458 = g(2) * t374;
t457 = g(3) * t377;
t456 = qJDD(1) * pkin(1);
t455 = t243 * t375;
t451 = t301 * t465;
t450 = t304 * t465;
t447 = t368 * t371;
t440 = -t371 * t387 + t324;
t328 = t445 - t446;
t439 = t465 * t328;
t330 = pkin(3) * t448 + t369 * qJ(2);
t364 = t369 ^ 2;
t438 = t371 ^ 2 + t364;
t437 = MDP(17) * t376;
t433 = qJD(2) * t369;
t427 = qJD(5) * t375;
t426 = t349 * MDP(16);
t420 = t375 * t269 - t372 * t270 - t301 * t427;
t325 = qJ(2) * t423 + t369 * t424 + qJDD(3);
t378 = qJD(1) ^ 2;
t413 = t438 * t378;
t408 = t376 * t279 - t284 * t373;
t406 = MDP(10) * (-t368 ^ 2 - t370 ^ 2);
t405 = qJD(5) * t243 + t238;
t404 = 0.2e1 * t438;
t293 = pkin(3) * t415 + t325;
t401 = g(2) * t377 + g(3) * t374;
t400 = -t457 + t458;
t399 = qJD(5) * t329 + t440;
t398 = qJD(5) * t328 + t439;
t317 = t329 * t369;
t318 = t328 * t369;
t272 = t375 * t317 + t318 * t372;
t273 = -t317 * t372 + t318 * t375;
t393 = t424 + t425;
t267 = -t393 * t447 + t287;
t294 = -qJ(2) * t447 + t327;
t392 = -t321 * qJD(1) - t294 * qJDD(1) - t267;
t391 = t322 * qJD(1) + t295 * qJDD(1) + t268;
t388 = t279 * t429 - t284 * t430 + t373 * t321 + t376 * t322;
t239 = -t304 * t428 + t420;
t386 = -t401 - t456;
t357 = qJDD(2) - t456;
t385 = -t357 - t386;
t382 = t404 * t424 + t458;
t362 = t377 * qJ(2);
t359 = cos(t367);
t358 = sin(t367);
t313 = -t358 * t374 - t359 * t443;
t312 = t358 * t443 - t359 * t374;
t311 = -t358 * t377 + t359 * t444;
t310 = t358 * t444 + t359 * t377;
t307 = t369 * t417 - t430 * t448;
t306 = t369 * t324;
t285 = pkin(4) * t307 + t433;
t282 = -t368 * t419 + t309;
t281 = pkin(4) * t317 + t330;
t251 = pkin(4) * t270 + t293;
t250 = -pkin(7) * t317 + t441;
t249 = -pkin(4) * t371 - pkin(7) * t318 + t408;
t248 = t273 * qJD(5) - t306 * t372 + t375 * t307;
t247 = -t272 * qJD(5) - t306 * t375 - t307 * t372;
t242 = pkin(7) * t306 - t469;
t241 = -pkin(7) * t307 + t388;
t1 = [(-g(2) * t313 + g(3) * t311 + t330 * t270 + t293 * t317 + t301 * t433 + t320 * t307 - t408 * t349 - t465 * t469) * MDP(17) + (-g(2) * t312 - g(3) * t310 + t330 * t269 + t293 * t318 + t304 * t433 - t320 * t306 + t441 * t349 - t388 * t465) * MDP(18) + (t268 * t295 + t283 * t322 + t267 * t294 + t282 * t321 - g(2) * (-pkin(1) * t377 - pkin(2) * t443 - qJ(2) * t374) - g(3) * (-pkin(1) * t374 - pkin(2) * t444 + t362)) * MDP(11) + (-g(3) * t362 + (-t357 + t401) * pkin(1) + (t425 * t438 + t382) * qJ(2)) * MDP(7) + (-(t249 * t375 - t250 * t372) * t343 + t285 * t258 + t281 * t240 + t251 * t272 + t274 * t248 - g(2) * t299 + g(3) * t297 + (t241 * t372 - t242 * t375 - (-t249 * t372 - t250 * t375) * qJD(5)) * t346) * MDP(24) + qJDD(1) * MDP(1) - t400 * MDP(3) + t401 * MDP(2) + (t248 * t346 + t272 * t343) * MDP(22) + (-t247 * t346 - t273 * t343) * MDP(21) + (-t307 * t465 + t317 * t349) * MDP(15) + (-t306 * t465 - t318 * t349) * MDP(14) + (t404 * t425 + t382 - t457) * MDP(6) + (t269 * t318 - t304 * t306) * MDP(12) + (-t269 * t317 - t270 * t318 + t301 * t306 - t304 * t307) * MDP(13) + (-g(2) * t298 - g(3) * t296 + t281 * t239 + t274 * t247 + t251 * t273 + t285 * t394 + ((-qJD(5) * t250 + t242) * t346 + t249 * t343) * t372 + ((qJD(5) * t249 + t241) * t346 + t250 * t343) * t375) * MDP(25) + (t239 * t273 + t247 * t394) * MDP(19) + (-t239 * t272 - t240 * t273 - t247 * t258 - t248 * t394) * MDP(20) + (t368 * MDP(8) + t370 * MDP(9)) * (t325 * t369 + t393 * t364 + t400) + ((t325 * qJ(2) + t401 * qJ(3) + t347 * qJD(2)) * MDP(11) - t385 * MDP(5) + (-t391 * t368 + t392 * t370 + t401) * MDP(10)) * t369 + (-t381 * MDP(17) - t389 * MDP(18) + t385 * MDP(4) + (t405 * t375 - t244 + t468) * MDP(25) + (-t412 - t466) * MDP(24) + (t401 * t370 + t392) * MDP(8) + (-t401 * t368 + t391) * MDP(9) + t240 * MDP(22) - t239 * MDP(21) + t270 * MDP(15) - t269 * MDP(14) + t426 + t332) * t371; -MDP(4) * t422 - MDP(6) * t413 + (-qJ(2) * t413 + qJDD(2) + t386) * MDP(7) + (-t368 * t413 - t414) * MDP(8) + (t368 * t422 - t370 * t413) * MDP(9) + (t267 * t370 + t268 * t368 + (-t347 * t369 + (t282 * t368 - t283 * t370) * t371) * qJD(1) - t401) * MDP(11) + (-t301 * t435 - t328 * t349 - t440 * t465) * MDP(17) + (-t304 * t435 + t329 * t349 - t439 * t465) * MDP(18) + (-(t328 * t375 - t329 * t372) * t343 - t258 * t435 + (t398 * t372 + t399 * t375) * t346) * MDP(24) + ((t328 * t372 + t329 * t375) * t343 - t394 * t435 + (-t399 * t372 + t398 * t375) * t346) * MDP(25) + (MDP(5) + t406) * t423; (g(1) * t371 + t325) * MDP(11) + (-t331 + t450) * MDP(17) + (t334 - t451) * MDP(18) + (t240 - t453) * MDP(24) + (t239 + t452) * MDP(25) + t378 * t364 * t406 + (t400 * MDP(11) + (-MDP(8) * t370 + MDP(9) * t368) * t378 * t371 + ((MDP(17) * t373 + MDP(9)) * t370 + (-MDP(18) * t373 + MDP(8) + t437) * t368) * qJDD(1) + ((t282 * t370 + t283 * t368) * MDP(11) + (-t329 * MDP(18) + t370 * t437) * qJD(4)) * qJD(1)) * t369; t304 * t301 * MDP(12) + (-t301 ^ 2 + t304 ^ 2) * MDP(13) + (t269 + t451) * MDP(14) + (-t270 + t450) * MDP(15) - t426 + (-g(2) * t310 + g(3) * t312 - t320 * t304 + t358 * t459 - t396 * t465 + t381) * MDP(17) + (-g(2) * t311 - g(3) * t313 + t320 * t301 + t359 * t459 + t409 * t465 + t389) * MDP(18) + (t239 - t452) * MDP(21) + (-t240 - t453) * MDP(22) + ((-t245 * t372 - t454) * t346 + t466 + (-t258 * t304 - t343 * t375 + t346 * t428) * pkin(4) + t461) * MDP(24) + ((t246 * t346 - t237) * t372 + (-t245 * t346 - t405) * t375 + (-t304 * t394 + t343 * t372 + t346 * t427) * pkin(4) + t464) * MDP(25) + t463; (t420 - t452) * MDP(21) + (-t410 - t453) * MDP(22) + (t397 * t346 + t461) * MDP(24) + (-t375 * t238 - t468 - (-t246 * t372 + t455) * t346 + t464) * MDP(25) + (-MDP(21) * t449 - t394 * MDP(22) + t397 * MDP(24) - MDP(25) * t455) * qJD(5) + t463;];
tau = t1;
