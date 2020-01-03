% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:34
% EndTime: 2020-01-03 11:50:40
% DurationCPUTime: 3.73s
% Computational Cost: add. (2167->336), mult. (5184->449), div. (0->0), fcn. (3628->10), ass. (0->179)
t375 = sin(pkin(8));
t377 = sin(qJ(4));
t381 = cos(qJ(3));
t427 = qJDD(1) * t381;
t378 = sin(qJ(3));
t428 = qJDD(1) * t378;
t380 = cos(qJ(4));
t493 = t375 * t380;
t324 = t377 * t381 + t378 * t380;
t426 = qJD(3) + qJD(4);
t494 = t426 * t324;
t264 = (qJD(1) * t494 + t377 * t428) * t375 - t427 * t493;
t370 = t375 ^ 2;
t485 = 0.2e1 * t370;
t376 = cos(pkin(8));
t431 = qJDD(1) * qJ(2);
t433 = qJD(1) * qJD(2);
t397 = t431 + t433;
t492 = t376 * t397;
t331 = -pkin(2) * t376 - pkin(6) * t375 - pkin(1);
t321 = t381 * t331;
t469 = t375 * t381;
t425 = pkin(7) * t469;
t477 = qJ(2) * t378;
t287 = -t425 + t321 + (-pkin(3) - t477) * t376;
t476 = qJ(2) * t381;
t350 = t376 * t476;
t452 = t378 * t331 + t350;
t470 = t375 * t378;
t292 = -pkin(7) * t470 + t452;
t457 = t377 * t287 + t380 * t292;
t443 = qJD(1) * t375;
t422 = t378 * t443;
t404 = t377 * t422;
t441 = qJD(1) * t381;
t421 = t375 * t441;
t299 = t380 * t421 - t404;
t293 = t299 * qJ(5);
t312 = qJD(1) * t331 + qJD(2);
t304 = t381 * t312;
t394 = -t376 * t477 - t425;
t283 = qJD(1) * t394 + t304;
t442 = qJD(1) * t376;
t351 = -qJD(3) + t442;
t269 = -pkin(3) * t351 + t283;
t284 = -pkin(7) * t422 + qJD(1) * t350 + t378 * t312;
t276 = t377 * t284;
t414 = t380 * t269 - t276;
t491 = t293 - t414;
t379 = sin(qJ(1));
t382 = cos(qJ(1));
t490 = -g(2) * t382 - g(3) * t379;
t489 = -MDP(13) * t378 - MDP(14) * t381;
t488 = t376 * MDP(4) - MDP(5) * t375;
t374 = qJ(3) + qJ(4);
t361 = sin(t374);
t362 = cos(t374);
t467 = t376 * t379;
t305 = -t361 * t467 - t362 * t382;
t460 = t382 * t361;
t307 = -t362 * t379 + t376 * t460;
t483 = g(1) * t375;
t487 = -g(2) * t305 - g(3) * t307 + t361 * t483;
t486 = t299 ^ 2;
t371 = t376 ^ 2;
t484 = pkin(7) * t375;
t481 = g(2) * t379;
t479 = g(3) * t382;
t395 = qJD(1) * t324;
t296 = t375 * t395;
t475 = qJ(5) * t296;
t474 = qJDD(1) * pkin(1);
t313 = pkin(3) * t422 + qJ(2) * t443;
t281 = pkin(4) * t296 + qJD(5) + t313;
t473 = t281 * t299;
t472 = (-qJ(5) - pkin(7) - pkin(6)) * t375;
t383 = qJD(1) ^ 2;
t471 = t370 * t383;
t466 = t376 * t382;
t465 = t377 * t378;
t464 = t378 * t379;
t463 = t378 * t382;
t462 = t379 * t381;
t278 = t380 * t284;
t461 = t381 * t382;
t343 = -qJD(4) + t351;
t250 = -pkin(4) * t343 - t491;
t459 = t250 + t491;
t458 = t380 * t283 - t276;
t323 = t380 * t381 - t465;
t456 = (t426 - t442) * t323;
t455 = t376 * t395 - t494;
t438 = qJD(3) * t381;
t440 = qJD(2) * t376;
t454 = t331 * t438 + t381 * t440;
t453 = t426 * t404;
t330 = t381 * pkin(3) + pkin(4) * t362;
t346 = t375 * pkin(3) * t438;
t319 = t375 * qJD(2) + t346;
t352 = pkin(3) * t470;
t322 = t375 * qJ(2) + t352;
t451 = t382 * pkin(1) + t379 * qJ(2);
t449 = t370 + t371;
t373 = t381 ^ 2;
t448 = t378 ^ 2 - t373;
t447 = MDP(10) * t375;
t446 = MDP(11) * t375;
t439 = qJD(3) * t312;
t437 = qJD(4) * t377;
t436 = qJD(4) * t380;
t429 = qJDD(1) * t376;
t349 = -qJDD(3) + t429;
t337 = -qJDD(4) + t349;
t328 = t337 * MDP(19);
t435 = t349 * MDP(12);
t434 = qJD(3) + t351;
t432 = qJD(1) * qJD(3);
t430 = qJDD(1) * t375;
t423 = qJ(2) * qJD(3) * t376;
t420 = qJ(2) * t429;
t419 = t381 * t433;
t418 = t381 * t432;
t416 = t377 * t427;
t311 = qJDD(1) * t331 + qJDD(2);
t303 = t381 * t311;
t403 = qJD(1) * t423;
t259 = -pkin(3) * t349 + t303 + (-pkin(7) * t430 - t403) * t381 + (-t420 - t439 + (qJD(3) * t484 - t440) * qJD(1)) * t378;
t406 = t378 * t311 + t312 * t438 + t376 * t419 + t381 * t420;
t389 = -t378 * t403 + t406;
t263 = (-t418 - t428) * t484 + t389;
t415 = t380 * t259 - t377 * t263;
t279 = qJD(3) * t394 + t454;
t280 = -t378 * t440 + (-t350 + (-t331 + t484) * t378) * qJD(3);
t413 = -t279 * t377 + t380 * t280;
t412 = -t283 * t377 - t278;
t411 = t380 * t287 - t292 * t377;
t410 = qJD(1) * t434;
t409 = t426 * t381;
t408 = t349 + t429;
t407 = -t377 * t259 - t380 * t263 - t269 * t436 + t284 * t437;
t289 = qJ(2) * t430 + qJD(1) * t346 + qJDD(1) * t352 + t375 * t433;
t405 = 0.2e1 * t449;
t401 = -t479 + t481;
t400 = -t269 * t377 - t278;
t398 = qJD(3) * (t351 + t442);
t396 = t380 * t279 + t377 * t280 + t287 * t436 - t292 * t437;
t265 = (t416 + (qJD(1) * t409 + t428) * t380) * t375 - t453;
t393 = pkin(4) * t265 + qJDD(5) + t289;
t294 = t296 ^ 2;
t392 = t299 * t296 * MDP(15) + (-t296 * t343 - t264) * MDP(17) + (-t299 * t343 + (-t416 + (-t426 * t441 - t428) * t380) * t375 + t453) * MDP(18) + (-t294 + t486) * MDP(16) - t328;
t390 = t405 * t433 + t479;
t388 = t400 * qJD(4) + t415;
t306 = t362 * t467 - t460;
t308 = t361 * t379 + t362 * t466;
t387 = g(2) * t306 - g(3) * t308 + t313 * t296 + t362 * t483 + t407;
t385 = -t313 * t299 + t388 + t487;
t364 = t379 * pkin(1);
t358 = qJDD(2) - t474;
t357 = pkin(3) * t380 + pkin(4);
t329 = pkin(3) * t378 + pkin(4) * t361;
t325 = pkin(2) + t330;
t317 = t376 * t461 + t464;
t316 = t376 * t463 - t462;
t315 = t376 * t462 - t463;
t314 = -t376 * t464 - t461;
t310 = t323 * t375;
t309 = t324 * t375;
t275 = -t426 * t375 * t465 + t409 * t493;
t274 = t494 * t375;
t262 = -qJ(5) * t309 + t457;
t260 = -pkin(4) * t376 - qJ(5) * t310 + t411;
t254 = -t293 + t458;
t253 = t412 + t475;
t252 = -t400 - t475;
t249 = qJ(5) * t274 - qJD(4) * t457 - qJD(5) * t310 + t413;
t248 = -qJ(5) * t275 - qJD(5) * t309 + t396;
t247 = -qJ(5) * t265 - qJD(5) * t296 - t407;
t246 = -pkin(4) * t337 + qJ(5) * t264 - qJD(5) * t299 + t388;
t1 = [qJDD(1) * MDP(1) + t490 * MDP(2) + t401 * MDP(3) + (t405 * t431 + t390 - t481) * MDP(6) + (-t358 * pkin(1) - g(2) * t451 - g(3) * t364 + (t449 * t431 + t390) * qJ(2)) * MDP(7) + (qJDD(1) * t373 - 0.2e1 * t378 * t418) * t370 * MDP(8) + (-t378 * t427 + t432 * t448) * MDP(9) * t485 + (t378 * t398 - t381 * t408) * t447 + (t378 * t408 + t381 * t398) * t446 + t376 * t435 + (-g(2) * t317 - g(3) * t315 - t303 * t376 - t321 * t349 + (t351 * t376 + (t485 + t371) * qJD(1)) * qJ(2) * t438 + (qJD(3) * t331 * t351 + t397 * t485 + (qJ(2) * t349 + qJD(2) * t351 + t439 + t492) * t376) * t378) * MDP(13) + ((-t378 * t423 + t454) * t351 + t452 * t349 + t389 * t376 + g(2) * t316 - g(3) * t314 + (t419 + (-t378 * t432 + t427) * qJ(2)) * t485) * MDP(14) + (-t264 * t310 - t274 * t299) * MDP(15) + (t264 * t309 - t265 * t310 + t274 * t296 - t275 * t299) * MDP(16) + (t264 * t376 + t274 * t343 - t310 * t337) * MDP(17) + (t265 * t376 + t275 * t343 + t309 * t337) * MDP(18) + t376 * t328 + (-t413 * t343 - t411 * t337 - t415 * t376 + t319 * t296 + t322 * t265 + t289 * t309 + t313 * t275 - g(2) * t308 - g(3) * t306 + (t343 * t457 - t376 * t400) * qJD(4)) * MDP(20) + (g(2) * t307 - g(3) * t305 - t322 * t264 - t313 * t274 + t289 * t310 + t319 * t299 + t337 * t457 + t343 * t396 - t376 * t407) * MDP(21) + (-t246 * t310 - t247 * t309 - t248 * t296 - t249 * t299 + t250 * t274 - t252 * t275 + t260 * t264 - t262 * t265 + t375 * t490) * MDP(22) + (t247 * t262 + t252 * t248 + t246 * t260 + t250 * t249 + t393 * (pkin(4) * t309 + t322) + t281 * (pkin(4) * t275 + t319) - g(2) * (t329 * t379 + t451) - g(3) * (t325 * t467 - t379 * t472 + t364) + (-g(2) * (t325 * t376 - t472) - g(3) * (-qJ(2) - t329)) * t382) * MDP(23) + t488 * (-t358 + t490 + t474); (qJDD(2) - t490) * MDP(7) + (-t296 * t443 - t323 * t337) * MDP(20) + (-t299 * t443 + t324 * t337) * MDP(21) + (t264 * t323 - t265 * t324 - t296 * t456 - t299 * t455) * MDP(22) + (t246 * t323 + t247 * t324 + t250 * t455 + t252 * t456 - t281 * t443 - t490) * MDP(23) + (-MDP(13) * t381 + MDP(14) * t378) * t349 + (-MDP(20) * t455 + MDP(21) * t456) * t343 + (-pkin(1) * MDP(7) - t488) * qJDD(1) + (-t371 * MDP(6) + (-MDP(6) + t489) * t370 - t449 * MDP(7) * qJ(2)) * t383 + t489 * t351 ^ 2; t381 * t378 * MDP(8) * t471 - t448 * MDP(9) * t471 + (-t378 * t410 + t427) * t447 + (-t381 * t410 - t428) * t446 - t435 + (-g(2) * t314 - g(3) * t316 + t303 + (-t376 * t410 - t471) * t476 + (-t312 * t434 + t483 - t492) * t378) * MDP(13) + (g(1) * t469 + g(2) * t315 - g(3) * t317 - t304 * t351 + (t434 * t442 + t471) * t477 - t406) * MDP(14) + (t412 * t343 + (-t296 * t421 - t337 * t380 + t343 * t437) * pkin(3) + t385) * MDP(20) + (-t458 * t343 + (-t299 * t421 + t377 * t337 + t343 * t436) * pkin(3) + t387) * MDP(21) + (t264 * t357 + (t252 + t253) * t299 + (-t250 + t254) * t296 + (-t265 * t377 + (-t296 * t380 + t299 * t377) * qJD(4)) * pkin(3)) * MDP(22) + (t246 * t357 - t252 * t254 - t250 * t253 - pkin(4) * t473 + t329 * t483 - g(2) * (-t329 * t467 - t330 * t382) - g(3) * (t329 * t466 - t330 * t379) + (-t281 * t421 + t247 * t377 + (-t250 * t377 + t252 * t380) * qJD(4)) * pkin(3)) * MDP(23) + t392; (t343 * t400 + t385) * MDP(20) + (-t343 * t414 + t387) * MDP(21) + (pkin(4) * t264 - t296 * t459) * MDP(22) + (t459 * t252 + (t246 - t473 + t487) * pkin(4)) * MDP(23) + t392; (-t294 - t486) * MDP(22) + (g(1) * t376 + t250 * t299 + t252 * t296 - t375 * t401 + t393) * MDP(23);];
tau = t1;
