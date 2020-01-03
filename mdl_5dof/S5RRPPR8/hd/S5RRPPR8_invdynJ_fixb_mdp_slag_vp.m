% Calculate vector of inverse dynamics joint torques for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:24
% EndTime: 2019-12-31 19:39:29
% DurationCPUTime: 4.55s
% Computational Cost: add. (1732->392), mult. (3754->506), div. (0->0), fcn. (2526->10), ass. (0->183)
t416 = sin(pkin(8));
t417 = cos(pkin(8));
t420 = sin(qJ(2));
t423 = cos(qJ(2));
t350 = t416 * t420 + t417 * t423;
t475 = qJD(1) * qJD(2);
t466 = t420 * t475;
t367 = t417 * t466;
t465 = t423 * t475;
t298 = qJDD(1) * t350 + t416 * t465 - t367;
t343 = t350 * qJD(2);
t473 = qJDD(1) * t423;
t474 = qJDD(1) * t420;
t452 = -t416 * t473 + t417 * t474;
t299 = qJD(1) * t343 + t452;
t338 = t350 * qJD(1);
t482 = qJD(1) * t423;
t483 = qJD(1) * t420;
t340 = -t416 * t482 + t417 * t483;
t419 = sin(qJ(5));
t422 = cos(qJ(5));
t477 = qJD(5) * t422;
t478 = qJD(5) * t419;
t261 = t422 * t298 + t419 * t299 - t338 * t478 + t340 * t477;
t291 = -t338 * t419 + t340 * t422;
t407 = qJDD(2) - qJDD(5);
t409 = qJD(2) - qJD(5);
t496 = t291 * t409;
t495 = t340 * t419;
t506 = t422 * t338 + t495;
t516 = -(t261 + t496) * MDP(22) + MDP(19) * t291 * t506 + (t291 ^ 2 - t506 ^ 2) * MDP(20) - t407 * MDP(23);
t497 = t506 * t409;
t513 = t465 + t474;
t386 = pkin(6) * t474;
t464 = pkin(6) * t465 + qJDD(3) + t386;
t479 = qJD(4) * t420;
t504 = pkin(2) + pkin(3);
t283 = -qJ(4) * t513 - qJD(1) * t479 - t504 * qJDD(2) + t464;
t481 = qJD(2) * t420;
t500 = pkin(6) - qJ(4);
t334 = -qJD(4) * t423 - t500 * t481;
t387 = pkin(6) * t473;
t410 = qJDD(2) * qJ(3);
t411 = qJD(2) * qJD(3);
t470 = t387 + t410 + t411;
t285 = -qJ(4) * t473 + qJD(1) * t334 + t470;
t461 = -t417 * t283 + t285 * t416;
t258 = -qJDD(2) * pkin(4) - pkin(7) * t299 - t461;
t263 = t416 * t283 + t417 * t285;
t259 = -pkin(7) * t298 + t263;
t348 = -qJD(1) * pkin(1) - pkin(2) * t482 - qJ(3) * t483;
t314 = pkin(3) * t482 + qJD(4) - t348;
t286 = pkin(4) * t338 + t314;
t421 = sin(qJ(1));
t408 = pkin(8) + qJ(5);
t394 = sin(t408);
t395 = cos(t408);
t442 = t394 * t423 - t395 * t420;
t310 = t442 * t421;
t424 = cos(qJ(1));
t489 = t423 * t424;
t493 = t420 * t424;
t312 = t394 * t489 - t395 * t493;
t441 = t394 * t420 + t395 * t423;
t511 = -g(1) * t312 - g(2) * t310 - g(3) * t441 - t422 * t258 + t419 * t259 + t286 * t291;
t311 = t441 * t421;
t313 = t441 * t424;
t510 = -g(1) * t313 - g(2) * t311 + g(3) * t442 + t419 * t258 + t422 * t259 - t286 * t506;
t487 = t423 * pkin(2) + t420 * qJ(3);
t508 = -pkin(1) - t487;
t507 = g(1) * t424 + g(2) * t421;
t499 = pkin(6) * qJDD(2);
t505 = (qJD(1) * t508 + t348) * qJD(2) - t499;
t503 = pkin(7) * t338;
t502 = pkin(7) * t340;
t405 = g(1) * t421;
t501 = g(2) * t424;
t399 = t423 * pkin(3);
t415 = qJDD(1) * pkin(1);
t498 = qJDD(2) * pkin(2);
t494 = t420 * t421;
t427 = qJD(1) ^ 2;
t492 = t420 * t427;
t491 = t421 * t423;
t366 = t500 * t423;
t335 = qJD(2) * t366 - t479;
t280 = t417 * t334 + t416 * t335;
t390 = pkin(6) * t483;
t356 = qJ(4) * t483 - t390;
t467 = t504 * qJD(2);
t322 = qJD(3) - t467 - t356;
t391 = pkin(6) * t482;
t358 = -qJ(4) * t482 + t391;
t412 = qJD(2) * qJ(3);
t347 = t358 + t412;
t278 = t416 * t322 + t417 * t347;
t303 = t417 * t356 + t416 * t358;
t365 = t500 * t420;
t306 = t416 * t365 + t417 * t366;
t396 = t420 * qJD(3);
t480 = qJD(2) * t423;
t488 = qJ(3) * t480 + t396;
t413 = t420 ^ 2;
t414 = t423 ^ 2;
t485 = t413 - t414;
t472 = t423 * t492;
t471 = -t419 * t298 + t422 * t299 - t338 * t477;
t469 = t399 + t487;
t468 = -g(1) * t493 - g(2) * t494 + g(3) * t423;
t463 = t405 - t501;
t462 = -qJD(2) * pkin(2) + qJD(3);
t359 = -qJ(3) * t416 - t417 * t504;
t277 = t417 * t322 - t347 * t416;
t279 = -t334 * t416 + t417 * t335;
t302 = -t356 * t416 + t417 * t358;
t305 = t417 * t365 - t366 * t416;
t349 = pkin(1) + t469;
t459 = qJD(3) * t416 + t302;
t458 = qJD(3) * t417 - t303;
t457 = t424 * pkin(1) + pkin(2) * t489 + t421 * pkin(6) + qJ(3) * t493;
t456 = -t386 - t468;
t455 = t420 * t467;
t426 = qJD(2) ^ 2;
t454 = pkin(6) * t426 + t501;
t451 = pkin(2) * t473 + qJ(3) * t513 + qJD(1) * t396 + t415;
t267 = -qJD(2) * pkin(4) + t277 - t502;
t268 = t278 - t503;
t450 = t422 * t267 - t419 * t268;
t449 = -t419 * t267 - t422 * t268;
t448 = t277 * t416 - t278 * t417;
t438 = t416 * t423 - t417 * t420;
t281 = pkin(7) * t438 + t305;
t282 = -pkin(7) * t350 + t306;
t447 = t281 * t422 - t282 * t419;
t446 = t281 * t419 + t282 * t422;
t300 = t422 * t350 - t419 * t438;
t301 = -t350 * t419 - t422 * t438;
t355 = -pkin(4) + t359;
t360 = qJ(3) * t417 - t416 * t504;
t445 = t355 * t422 - t360 * t419;
t444 = t355 * t419 + t360 * t422;
t362 = t390 + t462;
t364 = t391 + t412;
t443 = t362 * t423 - t364 * t420;
t440 = t416 * t422 + t417 * t419;
t439 = t416 * t419 - t417 * t422;
t384 = qJ(3) * t482;
t331 = -t504 * t483 + t384;
t330 = t464 - t498;
t436 = -0.2e1 * pkin(1) * t475 - t499;
t435 = t409 * t440;
t434 = t409 * t439;
t309 = -t455 + t488;
t260 = -t340 * t478 + t471;
t433 = -t454 + 0.2e1 * t415;
t432 = pkin(3) * t473 + qJDD(4) + t451;
t296 = pkin(2) * t466 - t451;
t336 = pkin(2) * t481 - t488;
t431 = -qJD(1) * t336 - qJDD(1) * t508 - t296 - t454;
t321 = -pkin(6) * t466 + t470;
t428 = qJD(2) * t443 + t321 * t423 + t330 * t420;
t271 = -qJD(1) * t455 + t432;
t401 = t424 * pkin(6);
t378 = g(1) * t491;
t374 = qJ(3) * t489;
t372 = qJ(3) * t491;
t357 = pkin(2) * t483 - t384;
t342 = t416 * t480 - t417 * t481;
t329 = t350 * t424;
t328 = t438 * t424;
t327 = t350 * t421;
t326 = t438 * t421;
t304 = pkin(4) * t350 + t349;
t297 = -pkin(4) * t340 + t331;
t284 = pkin(4) * t342 + t309;
t273 = t303 + t502;
t272 = t302 - t503;
t270 = -pkin(7) * t342 + t280;
t269 = -pkin(7) * t343 + t279;
t266 = qJD(5) * t301 + t422 * t342 + t343 * t419;
t265 = -qJD(5) * t300 - t342 * t419 + t343 * t422;
t264 = pkin(4) * t298 + t271;
t1 = [qJDD(1) * MDP(1) + t463 * MDP(2) + t507 * MDP(3) + (qJDD(1) * t413 + 0.2e1 * t420 * t465) * MDP(4) + 0.2e1 * (t420 * t473 - t485 * t475) * MDP(5) + (qJDD(2) * t420 + t423 * t426) * MDP(6) + (qJDD(2) * t423 - t420 * t426) * MDP(7) + (t420 * t436 + t423 * t433 + t378) * MDP(9) + (t436 * t423 + (-t433 - t405) * t420) * MDP(10) + (t505 * t420 + t431 * t423 + t378) * MDP(11) + ((t413 + t414) * qJDD(1) * pkin(6) + t428 - t507) * MDP(12) + (-t505 * t423 + (t431 + t405) * t420) * MDP(13) + (pkin(6) * t428 - g(1) * t401 - g(2) * t457 + t348 * t336 + (t296 - t405) * t508) * MDP(14) + (g(1) * t327 - g(2) * t329 - qJD(2) * t279 - qJDD(2) * t305 + t271 * t350 + t298 * t349 + t309 * t338 + t314 * t342) * MDP(15) + (-g(1) * t326 + g(2) * t328 + qJD(2) * t280 + qJDD(2) * t306 - t271 * t438 + t299 * t349 + t309 * t340 + t314 * t343) * MDP(16) + (-t263 * t350 - t277 * t343 - t278 * t342 - t279 * t340 - t280 * t338 - t298 * t306 - t299 * t305 - t438 * t461 + t507) * MDP(17) + (t263 * t306 + t278 * t280 - t461 * t305 + t277 * t279 + t271 * t349 + t314 * t309 - g(1) * (-qJ(4) * t424 + t401) - g(2) * (pkin(3) * t489 + t457) + (-g(1) * (t508 - t399) + g(2) * qJ(4)) * t421) * MDP(18) + (t260 * t301 + t265 * t291) * MDP(19) + (-t260 * t300 - t261 * t301 - t265 * t506 - t266 * t291) * MDP(20) + (-t265 * t409 - t301 * t407) * MDP(21) + (t266 * t409 + t300 * t407) * MDP(22) + (t284 * t506 + t304 * t261 + t264 * t300 + t286 * t266 - (-qJD(5) * t446 + t269 * t422 - t270 * t419) * t409 - t447 * t407 + g(1) * t311 - g(2) * t313) * MDP(24) + (t284 * t291 + t304 * t260 + t264 * t301 + t286 * t265 + (qJD(5) * t447 + t269 * t419 + t270 * t422) * t409 + t446 * t407 - g(1) * t310 + g(2) * t312) * MDP(25); (t507 * t420 * t504 - g(1) * t374 - g(2) * t372 - g(3) * t469 - t448 * qJD(3) + t263 * t360 - t277 * t302 - t278 * t303 - t314 * t331 - t359 * t461) * MDP(18) + (-t445 * t407 - t297 * t506 + (t272 * t422 - t273 * t419) * t409 + qJD(3) * t435 + (t409 * t444 - t449) * qJD(5) + t511) * MDP(24) + (-t260 + t497) * MDP(21) + (t444 * t407 - t297 * t291 - (t272 * t419 + t273 * t422) * t409 - qJD(3) * t434 + (t409 * t445 + t450) * qJD(5) + t510) * MDP(25) + (-g(1) * t329 - g(2) * t327 + g(3) * t438 + qJD(2) * t458 + qJDD(2) * t360 - t314 * t338 - t331 * t340 + t263) * MDP(16) + (g(3) * t420 - t387 + (pkin(1) * t427 + t507) * t423) * MDP(10) + (t387 + 0.2e1 * t410 + 0.2e1 * t411 + (qJD(1) * t357 - g(3)) * t420 + (qJD(1) * t348 - t507) * t423) * MDP(13) + (-g(1) * t328 - g(2) * t326 - g(3) * t350 + qJD(2) * t459 - qJDD(2) * t359 + t314 * t340 - t331 * t338 + t461) * MDP(15) + (-t298 * t360 - t299 * t359 + (-t278 + t459) * t340 + (t277 - t458) * t338) * MDP(17) + (pkin(1) * t492 + t456) * MDP(9) + (t321 * qJ(3) + t364 * qJD(3) - t330 * pkin(2) - t348 * t357 - g(1) * (-pkin(2) * t493 + t374) - g(2) * (-pkin(2) * t494 + t372) - g(3) * t487 - t443 * qJD(1) * pkin(6)) * MDP(14) + MDP(7) * t473 + MDP(6) * t474 + (0.2e1 * t498 - qJDD(3) + (-t348 * t420 + t357 * t423) * qJD(1) + t456) * MDP(11) + t485 * MDP(5) * t427 - MDP(4) * t472 + ((-pkin(2) * t420 + qJ(3) * t423) * qJDD(1) + ((t364 - t412) * t420 + (-t362 + t462) * t423) * qJD(1)) * MDP(12) + qJDD(2) * MDP(8) - t516; (-qJDD(2) - t472) * MDP(11) + MDP(12) * t474 + (-t413 * t427 - t426) * MDP(13) + (-qJD(2) * t364 + t348 * t483 + t330 + t468) * MDP(14) + (-qJDD(2) * t417 - t338 * t483 - t416 * t426) * MDP(15) + (qJDD(2) * t416 - t340 * t483 - t417 * t426) * MDP(16) + (-t298 * t416 - t299 * t417 + (t338 * t417 - t340 * t416) * qJD(2)) * MDP(17) + (qJD(2) * t448 + t263 * t416 - t314 * t483 - t417 * t461 + t468) * MDP(18) + (t407 * t439 - t409 * t435 - t483 * t506) * MDP(24) + (-t291 * t483 + t407 * t440 + t409 * t434) * MDP(25); (t416 * t474 + t417 * t473 - t367) * MDP(15) + t452 * MDP(16) + (-t338 ^ 2 - t340 ^ 2) * MDP(17) + (t277 * t340 + t278 * t338 + t432 + t463) * MDP(18) + (t261 - t496) * MDP(24) + (t260 + t497) * MDP(25) + (-t340 * MDP(15) + t338 * MDP(16) + ((MDP(15) * t416 + MDP(16) * t417) * t423 + (t416 * MDP(16) - t504 * MDP(18)) * t420) * qJD(1)) * qJD(2); (t471 - t497) * MDP(21) + (t409 * t449 - t511) * MDP(24) + (-t409 * t450 - t510) * MDP(25) + (-MDP(21) * t495 + MDP(24) * t449 - MDP(25) * t450) * qJD(5) + t516;];
tau = t1;
