% Calculate vector of inverse dynamics joint torques for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:52
% EndTime: 2020-01-03 12:11:56
% DurationCPUTime: 2.40s
% Computational Cost: add. (2551->301), mult. (3847->364), div. (0->0), fcn. (2445->12), ass. (0->169)
t407 = sin(qJ(2));
t469 = qJD(2) * t407;
t388 = pkin(1) * t469;
t410 = cos(qJ(2));
t505 = pkin(1) * t410;
t472 = -qJD(1) * t388 + qJDD(1) * t505;
t397 = qJDD(1) + qJDD(2);
t504 = pkin(2) * t397;
t323 = -t472 - t504;
t404 = qJ(1) + qJ(2);
t390 = sin(t404);
t392 = cos(t404);
t436 = g(2) * t392 + g(3) * t390;
t517 = t323 + t436;
t409 = cos(qJ(3));
t400 = qJD(1) + qJD(2);
t506 = pkin(1) * t407;
t463 = qJD(1) * t506;
t352 = pkin(7) * t400 + t463;
t450 = pkin(8) * t400 + t352;
t311 = t450 * t409;
t405 = sin(qJ(4));
t304 = t405 * t311;
t406 = sin(qJ(3));
t310 = t450 * t406;
t307 = qJD(3) * pkin(3) - t310;
t507 = cos(qJ(4));
t448 = t507 * t307 - t304;
t342 = t405 * t409 + t507 * t406;
t320 = t342 * t400;
t493 = t320 * qJ(5);
t516 = t493 - t448;
t452 = qJD(4) * t507;
t456 = t507 * t409;
t515 = -qJD(3) * t456 - t409 * t452;
t485 = t405 * t406;
t429 = t456 - t485;
t394 = t409 * pkin(3);
t498 = pkin(2) + t394;
t508 = -pkin(7) - pkin(8);
t458 = qJD(3) * t508;
t349 = t406 * t458;
t350 = t409 * t458;
t470 = qJD(1) * t410;
t462 = pkin(1) * t470;
t366 = t508 * t406;
t393 = t409 * pkin(8);
t367 = pkin(7) * t409 + t393;
t475 = t405 * t366 + t507 * t367;
t514 = -t475 * qJD(4) + t342 * t462 - t405 * t349 + t507 * t350;
t465 = qJD(4) * t405;
t513 = -t507 * t349 - t405 * t350 - t366 * t452 + t367 * t465 + t429 * t462;
t381 = pkin(7) + t506;
t497 = -pkin(8) - t381;
t338 = t497 * t406;
t339 = t381 * t409 + t393;
t477 = t405 * t338 + t507 * t339;
t512 = g(2) * t390 - g(3) * t392;
t467 = qJD(3) * t406;
t387 = pkin(3) * t467;
t511 = t387 - t463;
t412 = qJD(3) ^ 2;
t510 = pkin(7) * t412 - t504;
t399 = qJD(3) + qJD(4);
t509 = t320 ^ 2;
t503 = pkin(2) * t400;
t403 = qJ(3) + qJ(4);
t391 = cos(t403);
t501 = g(1) * t391;
t496 = qJ(5) * t342;
t460 = t400 * t485;
t318 = -t400 * t456 + t460;
t322 = -t400 * t498 - t462;
t291 = pkin(4) * t318 + qJD(5) + t322;
t495 = t291 * t320;
t494 = t318 * qJ(5);
t389 = sin(t403);
t492 = t389 * t390;
t491 = t389 * t392;
t490 = t390 * t391;
t489 = t391 * t392;
t488 = t397 * t406;
t487 = t397 * t409;
t486 = t400 * t406;
t483 = t406 * t409;
t298 = t399 * t342;
t479 = -t298 * qJ(5) + qJD(5) * t429;
t482 = t479 - t513;
t433 = t399 * t485;
t297 = t433 + t515;
t432 = t297 * qJ(5) - t342 * qJD(5);
t481 = t432 + t514;
t272 = pkin(4) * t399 - t516;
t480 = t272 + t516;
t478 = -t507 * t310 - t304;
t473 = pkin(4) * t391 + t394;
t348 = pkin(2) + t473;
t398 = -qJ(5) + t508;
t476 = t390 * t348 + t392 * t398;
t401 = t406 ^ 2;
t471 = -t409 ^ 2 + t401;
t468 = qJD(2) * t410;
t466 = qJD(3) * t409;
t464 = qJDD(1) * t407;
t461 = pkin(1) * t468;
t383 = -pkin(2) - t505;
t306 = t507 * t311;
t455 = t400 * t469;
t454 = t400 * t467;
t453 = t400 * t466;
t451 = pkin(4) * t298 + t387;
t449 = qJD(3) * t497;
t447 = t310 * t405 - t306;
t446 = t507 * t338 - t339 * t405;
t445 = t392 * t348 - t390 * t398;
t444 = t507 * t366 - t367 * t405;
t443 = t400 * t463;
t294 = pkin(3) * t454 - t397 * t498 - t472;
t442 = g(2) * t491 + g(3) * t492 + t294 * t342 - t322 * t297;
t353 = -t462 - t503;
t441 = t353 * t466 + t517 * t406;
t440 = -t342 * t397 + t515 * t400;
t439 = -pkin(4) * t429 - t498;
t434 = t429 * t397;
t277 = t400 * t433 + t440;
t396 = qJDD(3) + qJDD(4);
t324 = pkin(7) * t397 + (qJD(1) * t468 + t464) * pkin(1);
t281 = -t352 * t466 + qJDD(3) * pkin(3) - t324 * t406 + (-t453 - t488) * pkin(8);
t282 = -t352 * t467 + t324 * t409 + (-t454 + t487) * pkin(8);
t430 = -t405 * t307 - t306;
t418 = t430 * qJD(4) + t507 * t281 - t405 * t282;
t259 = t396 * pkin(4) + t277 * qJ(5) - t320 * qJD(5) + t418;
t278 = t298 * t400 - t434;
t415 = t405 * t281 + t507 * t282 + t307 * t452 - t311 * t465;
t260 = -t278 * qJ(5) - t318 * qJD(5) + t415;
t274 = -t430 - t494;
t431 = -t259 * t342 + t260 * t429 + t272 * t297 - t274 * t298 - t512;
t308 = t406 * t449 + t409 * t461;
t309 = -t406 * t461 + t409 * t449;
t428 = t507 * t308 + t405 * t309 + t338 * t452 - t339 * t465;
t317 = t318 ^ 2;
t426 = t320 * t318 * MDP(14) + (-t440 + (t318 - t460) * t399) * MDP(16) + t434 * MDP(17) + (-t317 + t509) * MDP(15) + t396 * MDP(18);
t425 = -t353 * t400 - t324 + t512;
t424 = -g(2) * t489 - g(3) * t490 - t294 * t429 + t322 * t298;
t423 = (-t277 * t429 - t278 * t342 + t297 * t318 - t298 * t320) * MDP(15) + (-t277 * t342 - t297 * t320) * MDP(14) + (-t297 * t399 + t342 * t396) * MDP(16) + (-t298 * t399 + t396 * t429) * MDP(17) + 0.2e1 * (-t471 * t400 * qJD(3) + t397 * t483) * MDP(8) + (t397 * t401 + 0.2e1 * t406 * t453) * MDP(7) + (qJDD(3) * t409 - t406 * t412) * MDP(10) + (qJDD(3) * t406 + t409 * t412) * MDP(9) + t397 * MDP(4);
t422 = -t436 + t443;
t421 = pkin(1) * t455 + t381 * t412 + t383 * t397;
t420 = -pkin(7) * qJDD(3) + (t462 - t503) * qJD(3);
t419 = -qJDD(3) * t381 + (t383 * t400 - t461) * qJD(3);
t269 = pkin(4) * t278 + qJDD(5) + t294;
t417 = -t477 * qJD(4) - t405 * t308 + t507 * t309;
t414 = g(1) * t389 + g(2) * t490 - g(3) * t489 + t322 * t318 - t415;
t413 = g(2) * t492 - g(3) * t491 - t322 * t320 + t418 - t501;
t411 = cos(qJ(1));
t408 = sin(qJ(1));
t382 = t507 * pkin(3) + pkin(4);
t361 = t383 - t394;
t351 = t388 + t387;
t337 = t429 * qJ(5);
t331 = t353 * t467;
t293 = t337 + t475;
t292 = t444 - t496;
t286 = t337 + t477;
t285 = t446 - t496;
t276 = -t493 + t478;
t275 = t447 + t494;
t263 = t417 + t432;
t262 = t428 + t479;
t1 = [(t331 + t419 * t406 + (-t421 - t517) * t409) * MDP(12) + (t421 * t406 + t419 * t409 + t441) * MDP(13) + t423 + qJDD(1) * MDP(1) + ((t397 * t410 - t455) * pkin(1) - t436 + t472) * MDP(5) + (((-qJDD(1) - t397) * t407 + (-qJD(1) - t400) * t468) * pkin(1) + t512) * MDP(6) + (t361 * t278 + t351 * t318 + t446 * t396 + t417 * t399 + t424) * MDP(19) + (-g(2) * t411 - g(3) * t408) * MDP(2) + (g(2) * t408 - g(3) * t411) * MDP(3) + (t260 * t286 + t274 * t262 + t259 * t285 + t272 * t263 + t269 * (t439 - t505) + t291 * (t388 + t451) - g(2) * (pkin(1) * t411 + t445) - g(3) * (pkin(1) * t408 + t476)) * MDP(22) + (-t361 * t277 + t351 * t320 - t477 * t396 - t428 * t399 + t442) * MDP(20) + (-t262 * t318 - t263 * t320 + t277 * t285 - t278 * t286 + t431) * MDP(21); (t331 + t420 * t406 + (-t323 + t422 - t510) * t409) * MDP(12) + (t422 + t472) * MDP(5) + t423 + ((-t464 + (-qJD(2) + t400) * t470) * pkin(1) + t512) * MDP(6) + (t260 * t293 + t259 * t292 + t269 * t439 - g(2) * t445 - g(3) * t476 + (t451 - t463) * t291 + t482 * t274 + t481 * t272) * MDP(22) + (-t278 * t498 + t511 * t318 + t444 * t396 + t514 * t399 + t424) * MDP(19) + (t420 * t409 + (-t443 + t510) * t406 + t441) * MDP(13) + (t277 * t498 + t511 * t320 - t475 * t396 + t513 * t399 + t442) * MDP(20) + (t277 * t292 - t278 * t293 - t482 * t318 - t481 * t320 + t431) * MDP(21); MDP(9) * t488 + MDP(10) * t487 + qJDD(3) * MDP(11) + (-g(1) * t409 + t425 * t406) * MDP(12) + (g(1) * t406 + t425 * t409) * MDP(13) + (-t447 * t399 + (-t318 * t486 + t507 * t396 - t399 * t465) * pkin(3) + t413) * MDP(19) + (t478 * t399 + (-t320 * t486 - t405 * t396 - t399 * t452) * pkin(3) + t414) * MDP(20) + (t382 * t277 + (t274 + t275) * t320 + (-t272 + t276) * t318 + (-t278 * t405 + (-t507 * t318 + t320 * t405) * qJD(4)) * pkin(3)) * MDP(21) + (t259 * t382 - t274 * t276 - t272 * t275 - pkin(4) * t495 - g(1) * t473 - t512 * (-pkin(3) * t406 - pkin(4) * t389) + (-t291 * t486 + t260 * t405 + (-t272 * t405 + t507 * t274) * qJD(4)) * pkin(3)) * MDP(22) + t426 + (-MDP(7) * t483 + t471 * MDP(8)) * t400 ^ 2; (-t430 * t399 + t413) * MDP(19) + (t448 * t399 + t414) * MDP(20) + (pkin(4) * t277 - t480 * t318) * MDP(21) + (t480 * t274 + (t389 * t512 + t259 - t495 - t501) * pkin(4)) * MDP(22) + t426; (-t317 - t509) * MDP(21) + (t272 * t320 + t274 * t318 + t269 + t436) * MDP(22);];
tau = t1;
