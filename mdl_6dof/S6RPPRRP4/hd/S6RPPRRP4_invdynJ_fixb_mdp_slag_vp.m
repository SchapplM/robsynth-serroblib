% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:29
% EndTime: 2019-03-09 02:06:36
% DurationCPUTime: 4.40s
% Computational Cost: add. (3530->469), mult. (6248->588), div. (0->0), fcn. (3822->8), ass. (0->176)
t397 = cos(pkin(9));
t462 = qJD(1) * qJD(2);
t373 = t397 * t462;
t403 = -pkin(1) - pkin(2);
t365 = qJDD(1) * t403 + qJDD(2);
t396 = sin(pkin(9));
t463 = qJ(2) * qJDD(1);
t481 = t396 * t365 + t397 * t463;
t320 = t373 + t481;
t315 = -qJDD(1) * pkin(7) + t320;
t525 = -qJD(3) * qJD(4) - t315;
t402 = cos(qJ(4));
t461 = qJD(1) * qJD(4);
t444 = t402 * t461;
t400 = sin(qJ(4));
t459 = qJDD(1) * t400;
t524 = t444 + t459;
t399 = sin(qJ(5));
t401 = cos(qJ(5));
t470 = qJD(4) * t401;
t476 = qJD(1) * t400;
t351 = t399 * t476 + t470;
t468 = qJD(5) * t351;
t297 = qJDD(4) * t399 - t401 * t524 + t468;
t465 = t399 * qJD(4);
t352 = t401 * t476 - t465;
t469 = qJD(4) * t402;
t523 = -t400 * t297 + t352 * t469;
t512 = sin(qJ(1));
t513 = cos(qJ(1));
t349 = -t396 * t512 - t397 * t513;
t350 = t396 * t513 - t397 * t512;
t433 = g(1) * t349 + g(2) * t350;
t522 = g(3) * t402 - t400 * t433;
t366 = qJD(1) * t403 + qJD(2);
t477 = qJ(2) * qJD(1);
t337 = t396 * t366 + t397 * t477;
t326 = -qJD(1) * pkin(7) + t337;
t311 = qJD(3) * t402 - t400 * t326;
t521 = qJD(4) * t311;
t368 = qJD(1) * t402 + qJD(5);
t458 = MDP(22) + MDP(24);
t493 = t399 * t402;
t345 = t396 * t493 + t397 * t401;
t471 = qJD(4) * t400;
t448 = t396 * t471;
t489 = t401 * t402;
t519 = -qJD(5) * t345 - t401 * t448 - (t396 * t399 + t397 * t489) * qJD(1);
t422 = t396 * t489 - t397 * t399;
t518 = qJD(5) * t422 - t399 * t448 - (-t396 * t401 + t397 * t493) * qJD(1);
t384 = t402 * qJDD(1);
t517 = -t400 * t461 + t384;
t457 = MDP(23) - MDP(26);
t306 = -qJD(4) * pkin(4) - t311;
t282 = -pkin(5) * t351 + qJ(6) * t352 + t306;
t348 = -qJDD(5) - t517;
t510 = pkin(8) * t348;
t516 = t282 * t368 + t510;
t505 = pkin(8) * qJD(5);
t515 = -t368 * t505 + t522;
t514 = t352 ^ 2;
t511 = pkin(5) * t348;
t507 = g(3) * t400;
t504 = pkin(1) * qJDD(1);
t503 = qJ(6) * t348;
t453 = -t326 * t469 + t400 * t525;
t278 = -qJDD(4) * pkin(4) - qJDD(3) * t402 - t453;
t466 = qJD(5) * t401;
t446 = t400 * t466;
t298 = t399 * t459 + qJDD(4) * t401 - qJD(5) * t465 + (t402 * t465 + t446) * qJD(1);
t270 = -pkin(5) * t298 - qJ(6) * t297 + qJD(6) * t352 + t278;
t502 = t270 * t401;
t312 = t400 * qJD(3) + t402 * t326;
t307 = qJD(4) * pkin(8) + t312;
t336 = t366 * t397 - t396 * t477;
t325 = qJD(1) * pkin(3) - t336;
t438 = pkin(4) * t402 + pkin(8) * t400;
t308 = qJD(1) * t438 + t325;
t280 = t307 * t401 + t308 * t399;
t274 = qJ(6) * t368 + t280;
t501 = t274 * t368;
t500 = t280 * t368;
t499 = t351 * t352;
t498 = t351 * t368;
t497 = t352 * t368;
t358 = t397 * qJ(2) + t396 * t403;
t354 = -pkin(7) + t358;
t496 = t354 * t399;
t495 = t368 * t401;
t494 = t399 * t348;
t491 = t400 * t401;
t490 = t401 * t348;
t405 = qJD(1) ^ 2;
t488 = t402 * t405;
t487 = qJDD(3) + g(3);
t437 = -pkin(4) * t400 + pkin(8) * t402;
t356 = t437 * qJD(1);
t486 = t401 * t311 + t399 * t356;
t357 = -t396 * qJ(2) + t397 * t403;
t426 = pkin(3) + t438;
t327 = -t357 + t426;
t483 = t399 * t327 + t354 * t489;
t430 = pkin(5) * t399 - qJ(6) * t401;
t482 = -qJD(6) * t399 + t368 * t430 - t312;
t480 = t513 * pkin(1) + t512 * qJ(2);
t479 = g(1) * t512 - g(2) * t513;
t393 = t400 ^ 2;
t478 = -t402 ^ 2 + t393;
t474 = qJD(2) * t397;
t472 = qJD(4) * t352;
t467 = qJD(5) * t399;
t279 = -t307 * t399 + t308 * t401;
t464 = qJD(6) - t279;
t456 = 0.2e1 * t462;
t277 = qJDD(4) * pkin(8) + qJDD(3) * t400 + t315 * t402 + t521;
t371 = t396 * t462;
t443 = t365 * t397 - t396 * t463;
t319 = -t371 + t443;
t314 = qJDD(1) * pkin(3) - t319;
t288 = t517 * pkin(4) + pkin(8) * t524 + t314;
t455 = -t401 * t277 - t399 * t288 - t308 * t466;
t454 = -t352 * t446 - t399 * t523;
t341 = qJD(2) * t396 + qJD(4) * t437;
t449 = t402 * t474;
t452 = t327 * t466 + t399 * t341 + t401 * t449;
t451 = g(3) * t489 - t433 * t491;
t450 = t513 * pkin(2) + t480;
t442 = t458 * t351;
t441 = qJDD(2) - t504;
t440 = -t399 * t277 + t401 * t288 - t307 * t466 - t308 * t467;
t439 = -pkin(1) * t512 + t513 * qJ(2);
t299 = -t349 * t401 + t350 * t493;
t303 = -t349 * t493 - t350 * t401;
t436 = g(1) * t299 + g(2) * t303;
t300 = t349 * t399 + t350 * t489;
t304 = -t349 * t489 + t350 * t399;
t435 = -g(1) * t300 - g(2) * t304;
t434 = -g(1) * t350 + g(2) * t349;
t431 = -pkin(5) * t401 - qJ(6) * t399;
t418 = -t307 * t467 - t455;
t268 = qJD(6) * t368 + t418 - t503;
t269 = qJDD(6) - t440 + t511;
t429 = t268 * t401 + t269 * t399;
t273 = -pkin(5) * t368 + t464;
t428 = t273 * t401 - t274 * t399;
t427 = t273 * t399 + t274 * t401;
t404 = qJD(4) ^ 2;
t362 = qJDD(4) * t402 - t400 * t404;
t361 = -qJDD(4) * t400 - t402 * t404;
t425 = pkin(4) - t431;
t423 = t354 - t430;
t420 = t368 * t466 - t494;
t419 = t368 * t467 + t490;
t417 = g(1) * t513 + g(2) * t512;
t416 = qJD(1) * t325 - t433;
t415 = -t400 * t467 + t401 * t469;
t414 = -pkin(2) * t512 + t439;
t413 = t306 * t368 + t510;
t412 = t402 * t433 + t507;
t353 = pkin(3) - t357;
t410 = -qJDD(4) * t354 + (-qJD(1) * t353 - t325 - t474) * qJD(4);
t409 = g(1) * t303 - g(2) * t299 - t399 * t507 + t440;
t408 = -g(1) * t304 + g(2) * t300 + g(3) * t491 + t418;
t407 = qJDD(1) * t353 - t354 * t404 + t314 + t371 + t434;
t406 = -t282 * t352 + qJDD(6) - t409;
t313 = -pkin(5) * t352 - qJ(6) * t351;
t310 = t423 * t400;
t290 = -t327 * t401 + (-pkin(5) + t496) * t402;
t289 = qJ(6) * t402 + t483;
t287 = pkin(5) * t476 + t311 * t399 - t356 * t401;
t286 = -qJ(6) * t476 + t486;
t285 = t297 - t498;
t281 = t423 * t469 + (qJD(5) * t431 + qJD(6) * t401 + t474) * t400;
t272 = pkin(5) * t471 + (qJD(5) * t354 * t402 - t341) * t401 + (qJD(5) * t327 - t354 * t471 + t449) * t399;
t271 = (-t354 * t467 + qJD(6)) * t402 + (-t354 * t401 - qJ(6)) * t471 + t452;
t1 = [(t271 * t368 + t281 * t352 - t289 * t348 - t297 * t310 + (t282 * t470 + t268) * t402 + (-qJD(4) * t274 - t282 * t467 + t502) * t400 - t436) * MDP(26) + (-qJDD(2) + t479 + 0.2e1 * t504) * MDP(4) + ((t368 * t465 + t298) * t402 + (-qJD(4) * t351 + t420) * t400) * MDP(20) + (-t272 * t368 - t281 * t351 + t290 * t348 - t298 * t310 + (-t282 * t465 - t269) * t402 + (qJD(4) * t273 - t270 * t399 - t282 * t466) * t400 + t435) * MDP(24) + (-t417 + t456 + 0.2e1 * t463) * MDP(5) + (qJDD(1) * t393 + 0.2e1 * t400 * t444) * MDP(10) + (t268 * t289 + t274 * t271 + t270 * t310 + t282 * t281 + t269 * t290 + t273 * t272 - g(1) * (t300 * pkin(5) + t299 * qJ(6) + t414) - g(2) * (pkin(5) * t304 + qJ(6) * t303 + t450) + (-g(2) * pkin(7) - g(1) * t426) * t350 + (-g(1) * pkin(7) + g(2) * t426) * t349) * MDP(27) + (t320 * t358 + t319 * t357 - g(1) * t414 - g(2) * t450 + (-t336 * t396 + t337 * t397) * qJD(2)) * MDP(9) + (-qJDD(1) * t357 + 0.2e1 * t371 + t434 - t443) * MDP(7) + (-t452 * t368 + t483 * t348 + ((t354 * t368 + t307) * t467 + (-t306 * t401 - t352 * t354) * qJD(4) + t455) * t402 + (-t352 * t474 + t306 * t467 - t278 * t401 + t354 * t297 + (t354 * t495 + t280) * qJD(4)) * t400 + t436) * MDP(23) + ((-t327 * t467 + t341 * t401) * t368 - t327 * t490 + ((-t354 * t466 - t399 * t474) * t368 + t354 * t494 + (-t306 * t399 - t351 * t354) * qJD(4) + t440) * t402 + (-t351 * t474 - t306 * t466 - t278 * t399 - t354 * t298 + (t368 * t496 - t279) * qJD(4)) * t400 + t435) * MDP(22) + t417 * MDP(3) + (t400 * t410 + t402 * t407) * MDP(15) + (-t400 * t407 + t402 * t410) * MDP(16) + (-t298 * t491 - t351 * t415 + t454) * MDP(18) + (-t297 * t491 + t352 * t415) * MDP(17) + t479 * MDP(2) + (-t441 * pkin(1) - g(1) * t439 - g(2) * t480 + (t456 + t463) * qJ(2)) * MDP(6) + (qJDD(1) * t358 + 0.2e1 * t373 + t433 + t481) * MDP(8) + 0.2e1 * (t384 * t400 - t461 * t478) * MDP(11) + (t271 * t351 - t272 * t352 + t289 * t298 + t290 * t297 - t428 * t469 + (qJD(5) * t427 + t268 * t399 - t269 * t401 + t434) * t400) * MDP(25) + (-t348 * t402 - t368 * t471) * MDP(21) + ((-t368 * t470 + t297) * t402 + (t419 + t472) * t400) * MDP(19) + t361 * MDP(12) - t362 * MDP(13) + qJDD(1) * MDP(1); -qJDD(1) * MDP(4) - t405 * MDP(5) + (-qJ(2) * t405 + t441 - t479) * MDP(6) - t479 * MDP(9) + (t297 * t345 + t351 * t519 - t352 * t518) * MDP(25) + (t269 * t345 + t518 * t273 + t519 * t274 - t479) * MDP(27) + (t298 * MDP(25) + t268 * MDP(27) + t348 * t457) * t422 + (-t457 * t519 - t458 * t518) * t368 + (-t405 * MDP(8) + t319 * MDP(9) + (-MDP(15) * t402 + MDP(16) * t400 - MDP(7)) * qJDD(1) + (0.2e1 * MDP(16) * t469 - t337 * MDP(9) + (0.2e1 * MDP(15) * qJD(4) - t282 * MDP(27) + t352 * t457 + t442) * t400) * qJD(1)) * t397 + (-t405 * MDP(7) + qJDD(1) * MDP(8) + (qJD(1) * t336 + t320) * MDP(9) + (t361 - t488) * MDP(15) + (t400 * t405 - t362) * MDP(16) + (t270 * t400 + t282 * t469) * MDP(27) - t457 * t523) * t396 + t458 * (t345 * t348 + (-t298 * t400 - t351 * t469) * t396); t487 * MDP(9) + t362 * MDP(15) + t361 * MDP(16) + t454 * MDP(25) + g(3) * MDP(27) - t442 * t471 + (-t270 * MDP(27) + t458 * t298 - t457 * t297 + (t351 * t401 * MDP(25) + t427 * MDP(27) + (-t399 * t458 - t401 * t457) * t368) * qJD(4)) * t402 + ((t298 * t401 - t351 * t467) * MDP(25) + (qJD(4) * t282 + t273 * t466 - t274 * t467 + t429) * MDP(27) - t458 * t420 + t457 * (t419 - t472)) * t400; -t400 * MDP(10) * t488 + t478 * MDP(11) * t405 - MDP(12) * t459 - MDP(13) * t384 + qJDD(4) * MDP(14) + (qJD(4) * t312 + t400 * t416 + t402 * t487 + t453) * MDP(15) + (t521 + (qJD(4) * t326 - t487) * t400 + (t416 + t525) * t402) * MDP(16) + (t297 * t399 - t352 * t495) * MDP(17) + ((t297 + t498) * t401 + (t298 + t497) * t399) * MDP(18) + ((-t352 * t400 + t368 * t489) * qJD(1) + t420) * MDP(19) + ((t351 * t400 - t368 * t493) * qJD(1) - t419) * MDP(20) + t368 * MDP(21) * t476 + (t279 * t476 + pkin(4) * t298 + t312 * t351 + (-t278 + (-t356 - t505) * t368) * t401 + (t311 * t368 + t413) * t399 + t451) * MDP(22) + (-pkin(4) * t297 + t486 * t368 - t280 * t476 + t312 * t352 + t413 * t401 + (t278 - t515) * t399) * MDP(23) + (-t273 * t476 - t502 + t298 * t425 + (-pkin(8) * t466 + t287) * t368 - t482 * t351 + t516 * t399 + t451) * MDP(24) + (-t286 * t351 + t287 * t352 + (t268 + t368 * t273 + (-qJD(5) * t352 + t298) * pkin(8)) * t401 + (t269 - t501 + (t297 - t468) * pkin(8)) * t399 + t412) * MDP(25) + (t274 * t476 - t286 * t368 + t297 * t425 + t482 * t352 - t516 * t401 + (-t270 + t515) * t399) * MDP(26) + (-t273 * t287 - t274 * t286 + t482 * t282 + (qJD(5) * t428 + t412 + t429) * pkin(8) + (-t270 + t522) * t425) * MDP(27); MDP(17) * t499 + (-t351 ^ 2 + t514) * MDP(18) + t285 * MDP(19) + (t298 - t497) * MDP(20) - t348 * MDP(21) + (t306 * t352 + t409 + t500) * MDP(22) + (t279 * t368 - t306 * t351 - t408) * MDP(23) + (t313 * t351 - t406 + t500 - 0.2e1 * t511) * MDP(24) + (-pkin(5) * t297 + qJ(6) * t298 + (-t274 + t280) * t352 + (-t273 + t464) * t351) * MDP(25) + (-0.2e1 * t503 + t282 * t351 - t313 * t352 + (0.2e1 * qJD(6) - t279) * t368 + t408) * MDP(26) + (t268 * qJ(6) - t269 * pkin(5) - t282 * t313 - t273 * t280 - g(1) * (-pkin(5) * t303 + qJ(6) * t304) - g(2) * (pkin(5) * t299 - qJ(6) * t300) - t430 * t507 + t464 * t274) * MDP(27); (t348 + t499) * MDP(24) + t285 * MDP(25) + (-t368 ^ 2 - t514) * MDP(26) + (t406 - t501 + t511) * MDP(27);];
tau  = t1;
