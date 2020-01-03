% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:36
% EndTime: 2020-01-03 12:09:41
% DurationCPUTime: 2.61s
% Computational Cost: add. (2809->316), mult. (4325->398), div. (0->0), fcn. (3045->14), ass. (0->173)
t433 = sin(qJ(2));
t484 = qJD(2) * t433;
t414 = pkin(1) * t484;
t437 = cos(qJ(2));
t513 = pkin(1) * t437;
t488 = -qJD(1) * t414 + qJDD(1) * t513;
t422 = qJDD(1) + qJDD(2);
t512 = pkin(2) * t422;
t362 = -t488 - t512;
t427 = qJ(1) + qJ(2);
t417 = sin(t427);
t418 = cos(t427);
t462 = g(2) * t418 + g(3) * t417;
t522 = t362 + t462;
t521 = qJ(4) + pkin(7);
t428 = sin(pkin(9));
t429 = cos(pkin(9));
t432 = sin(qJ(3));
t436 = cos(qJ(3));
t378 = -t428 * t432 + t429 * t436;
t424 = qJD(1) + qJD(2);
t358 = t378 * t424;
t435 = cos(qJ(5));
t345 = t435 * t358;
t379 = t428 * t436 + t429 * t432;
t359 = t379 * t424;
t431 = sin(qJ(5));
t499 = t359 * t431;
t308 = t345 - t499;
t423 = qJD(3) + qJD(5);
t500 = t308 * t423;
t517 = g(2) * t417 - g(3) * t418;
t479 = qJDD(1) * t433;
t483 = qJD(2) * t437;
t363 = pkin(7) * t422 + (qJD(1) * t483 + t479) * pkin(1);
t449 = qJ(4) * t422 + qJD(4) * t424 + t363;
t514 = pkin(1) * t433;
t478 = qJD(1) * t514;
t470 = t521 * t424 + t478;
t453 = qJD(3) * t470;
t293 = qJDD(3) * pkin(3) - t432 * t449 - t436 * t453;
t296 = -t432 * t453 + t436 * t449;
t270 = t429 * t293 - t296 * t428;
t481 = qJD(3) * t436;
t472 = t424 * t481;
t482 = qJD(3) * t432;
t473 = t424 * t482;
t314 = t379 * t422 - t428 * t473 + t429 * t472;
t264 = qJDD(3) * pkin(4) - pkin(8) * t314 + t270;
t271 = t428 * t293 + t429 * t296;
t371 = t379 * qJD(3);
t313 = -t371 * t424 + t378 * t422;
t265 = pkin(8) * t313 + t271;
t347 = t470 * t432;
t341 = qJD(3) * pkin(3) - t347;
t348 = t470 * t436;
t496 = t429 * t348;
t300 = t428 * t341 + t496;
t508 = pkin(8) * t358;
t283 = t300 + t508;
t409 = pkin(3) * t436 + pkin(2);
t485 = qJD(1) * t437;
t477 = pkin(1) * t485;
t357 = -t409 * t424 + qJD(4) - t477;
t315 = -pkin(4) * t358 + t357;
t415 = qJ(3) + pkin(9) + qJ(5);
t402 = sin(t415);
t403 = cos(t415);
t480 = qJD(5) * t431;
t520 = g(1) * t402 - t431 * t264 - t435 * t265 + t283 * t480 - t315 * t308 + t517 * t403;
t421 = qJDD(3) + qJDD(5);
t455 = t358 * t431 + t435 * t359;
t519 = t421 * MDP(20) + (-t308 ^ 2 + t455 ^ 2) * MDP(17) - t308 * MDP(16) * t455;
t501 = t455 * t423;
t416 = t436 * qJD(4);
t471 = qJD(3) * t521;
t368 = -t432 * t471 + t416;
t369 = -qJD(4) * t432 - t436 * t471;
t492 = -t368 * t428 + t429 * t369 + t379 * t477;
t491 = t429 * t368 + t428 * t369 - t378 * t477;
t439 = qJD(3) ^ 2;
t516 = pkin(7) * t439 - t512;
t497 = t402 * t418;
t498 = t402 * t417;
t515 = -g(1) * t403 + g(2) * t498 - g(3) * t497 + t435 * t264 - t431 * t265 - t315 * t455;
t469 = -t435 * t313 + t314 * t431;
t269 = qJD(5) * t455 + t469;
t511 = pkin(2) * t424;
t510 = pkin(3) * t428;
t507 = pkin(8) * t359;
t372 = t378 * qJD(3);
t506 = pkin(8) * t372;
t505 = pkin(8) * t379;
t504 = g(1) * t436;
t335 = t428 * t348;
t495 = t432 * t436;
t299 = t429 * t341 - t335;
t281 = qJD(3) * pkin(4) + t299 - t507;
t494 = t435 * t281;
t408 = pkin(7) + t514;
t493 = -qJ(4) - t408;
t467 = qJD(3) * t493;
t476 = pkin(1) * t483;
t330 = t432 * t467 + t436 * t476 + t416;
t331 = (-qJD(4) - t476) * t432 + t436 * t467;
t295 = t429 * t330 + t428 * t331;
t302 = -t429 * t347 - t335;
t376 = t493 * t432;
t419 = t436 * qJ(4);
t377 = t408 * t436 + t419;
t322 = t428 * t376 + t429 * t377;
t394 = t521 * t432;
t395 = pkin(7) * t436 + t419;
t340 = -t428 * t394 + t429 * t395;
t490 = t417 * t409 - t418 * t521;
t425 = t432 ^ 2;
t487 = -t436 ^ 2 + t425;
t413 = pkin(3) * t482;
t475 = qJD(5) * t345 + t431 * t313 + t435 * t314;
t474 = t424 * t484;
t344 = pkin(4) * t371 + t413;
t294 = -t330 * t428 + t429 * t331;
t301 = t347 * t428 - t496;
t321 = t429 * t376 - t377 * t428;
t339 = -t429 * t394 - t395 * t428;
t468 = t418 * t409 + t417 * t521;
t466 = t424 * t478;
t320 = pkin(3) * t473 - t409 * t422 + qJDD(4) - t488;
t282 = -pkin(4) * t313 + t320;
t454 = t435 * t378 - t379 * t431;
t286 = qJD(5) * t454 - t371 * t431 + t372 * t435;
t324 = t378 * t431 + t379 * t435;
t465 = g(2) * t497 + g(3) * t498 + t282 * t324 + t315 * t286;
t390 = -t477 - t511;
t464 = t390 * t481 + t522 * t432;
t463 = t344 - t478;
t375 = t378 * pkin(8);
t319 = t375 + t340;
t460 = qJD(5) * t319 - t492 + t506;
t318 = t339 - t505;
t367 = t371 * pkin(8);
t459 = -qJD(5) * t318 + t367 - t491;
t458 = -t281 * t431 - t283 * t435;
t303 = t321 - t505;
t304 = t375 + t322;
t457 = t303 * t435 - t304 * t431;
t456 = t303 * t431 + t304 * t435;
t352 = -pkin(4) * t378 - t409;
t404 = pkin(3) * t429 + pkin(4);
t452 = t404 * t431 + t435 * t510;
t451 = t404 * t435 - t431 * t510;
t450 = -t270 * t379 + t271 * t378 - t299 * t372 - t300 * t371 - t517;
t268 = -t359 * t480 + t475;
t448 = -t390 * t424 - t363 + t517;
t287 = qJD(5) * t324 + t435 * t371 + t372 * t431;
t447 = -t282 * t454 + t315 * t287 - t403 * t462;
t446 = (t268 * t454 - t269 * t324 + t286 * t308 - t287 * t455) * MDP(17) + (t268 * t324 + t286 * t455) * MDP(16) + (t286 * t423 + t324 * t421) * MDP(18) + (-t287 * t423 + t421 * t454) * MDP(19) + 0.2e1 * (-t487 * t424 * qJD(3) + t422 * t495) * MDP(8) + (t422 * t425 + 0.2e1 * t432 * t472) * MDP(7) + (qJDD(3) * t436 - t432 * t439) * MDP(10) + (qJDD(3) * t432 + t436 * t439) * MDP(9) + t422 * MDP(4);
t445 = -t462 + t466;
t410 = -pkin(2) - t513;
t444 = pkin(1) * t474 + t408 * t439 + t410 * t422;
t442 = -pkin(7) * qJDD(3) + (t477 - t511) * qJD(3);
t441 = -qJDD(3) * t408 + (t410 * t424 - t476) * qJD(3);
t438 = cos(qJ(1));
t434 = sin(qJ(1));
t373 = t390 * t482;
t343 = t352 - t513;
t334 = t344 + t414;
t329 = pkin(3) * t424 * t432 + pkin(4) * t359;
t285 = t302 - t507;
t284 = t301 - t508;
t279 = -t367 + t295;
t278 = t294 - t506;
t1 = [(((-qJDD(1) - t422) * t433 + (-qJD(1) - t424) * t483) * pkin(1) + t517) * MDP(6) + qJDD(1) * MDP(1) + (-t334 * t308 + t343 * t269 + (-qJD(5) * t456 + t278 * t435 - t279 * t431) * t423 + t457 * t421 + t447) * MDP(21) + (t334 * t455 + t343 * t268 - (qJD(5) * t457 + t278 * t431 + t279 * t435) * t423 - t456 * t421 + t465) * MDP(22) + ((t422 * t437 - t474) * pkin(1) - t462 + t488) * MDP(5) + (t373 + t441 * t432 + (-t444 - t522) * t436) * MDP(12) + (t432 * t444 + t436 * t441 + t464) * MDP(13) + (-g(2) * t438 - g(3) * t434) * MDP(2) + (g(2) * t434 - g(3) * t438) * MDP(3) + (t271 * t322 + t300 * t295 + t270 * t321 + t299 * t294 + t320 * (-t409 - t513) + t357 * (t414 + t413) - g(2) * (pkin(1) * t438 + t468) - g(3) * (pkin(1) * t434 + t490)) * MDP(15) + (-t294 * t359 + t295 * t358 + t313 * t322 - t314 * t321 + t450) * MDP(14) + t446; (t442 * t436 + (-t466 + t516) * t432 + t464) * MDP(13) + ((-t479 + (-qJD(2) + t424) * t485) * pkin(1) + t517) * MDP(6) + (t373 + t442 * t432 + (-t362 + t445 - t516) * t436) * MDP(12) + (t445 + t488) * MDP(5) + (t313 * t340 - t314 * t339 + t491 * t358 - t492 * t359 + t450) * MDP(14) + t446 + (t352 * t269 + (t318 * t435 - t319 * t431) * t421 + (t431 * t459 - t435 * t460) * t423 - t463 * t308 + t447) * MDP(21) + (t352 * t268 - (t318 * t431 + t319 * t435) * t421 + (t431 * t460 + t435 * t459) * t423 + t463 * t455 + t465) * MDP(22) + (t271 * t340 + t270 * t339 - t320 * t409 - g(2) * t468 - g(3) * t490 + (t413 - t478) * t357 + t491 * t300 + t492 * t299) * MDP(15); qJDD(3) * MDP(11) + (t432 * t448 - t504) * MDP(12) + (g(1) * t432 + t436 * t448) * MDP(13) + ((t300 + t301) * t359 + (t299 - t302) * t358 + (t313 * t428 - t314 * t429) * pkin(3)) * MDP(14) + (-t299 * t301 - t300 * t302 + (-t504 + t270 * t429 + t271 * t428 + (-t357 * t424 + t517) * t432) * pkin(3)) * MDP(15) + (t268 - t500) * MDP(18) + (-t269 + t501) * MDP(19) + (t451 * t421 + t329 * t308 - (t284 * t435 - t285 * t431) * t423 + (-t423 * t452 + t458) * qJD(5) + t515) * MDP(21) + (-t452 * t421 - t329 * t455 + (t284 * t431 + t285 * t435) * t423 + (-t423 * t451 - t494) * qJD(5) + t520) * MDP(22) + (t436 * MDP(10) + t432 * MDP(9)) * t422 + (-MDP(7) * t495 + t487 * MDP(8)) * t424 ^ 2 + t519; (-t358 ^ 2 - t359 ^ 2) * MDP(14) + (t299 * t359 - t300 * t358 + t320 + t462) * MDP(15) + (t269 + t501) * MDP(21) + (t268 + t500) * MDP(22); (t475 - t500) * MDP(18) + (-t469 + t501) * MDP(19) + (-t423 * t458 + t515) * MDP(21) + ((-t283 * t431 + t494) * t423 + t520) * MDP(22) + (-MDP(18) * t499 - MDP(19) * t455 + MDP(21) * t458 - MDP(22) * t494) * qJD(5) + t519;];
tau = t1;
