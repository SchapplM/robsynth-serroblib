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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
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
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:59
% EndTime: 2022-01-23 09:33:04
% DurationCPUTime: 4.48s
% Computational Cost: add. (2862->387), mult. (6825->495), div. (0->0), fcn. (4829->10), ass. (0->197)
t402 = cos(qJ(3));
t455 = qJDD(1) * t402;
t396 = sin(pkin(8));
t401 = cos(qJ(4));
t522 = t396 * t401;
t355 = t455 * t522;
t398 = sin(qJ(4));
t399 = sin(qJ(3));
t456 = qJDD(1) * t399;
t491 = t398 * t402;
t346 = t399 * t401 + t491;
t454 = qJD(3) + qJD(4);
t523 = t454 * t346;
t281 = (qJD(1) * t523 + t398 * t456) * t396 - t355;
t391 = t396 ^ 2;
t514 = 0.2e1 * t391;
t400 = sin(qJ(1));
t403 = cos(qJ(1));
t443 = g(1) * t400 - g(2) * t403;
t397 = cos(pkin(8));
t460 = qJD(1) * qJD(2);
t461 = qJ(2) * qJDD(1);
t419 = t460 + t461;
t521 = t419 * t397;
t353 = -pkin(2) * t397 - pkin(6) * t396 - pkin(1);
t342 = t402 * t353;
t495 = t396 * t402;
t453 = pkin(7) * t495;
t506 = qJ(2) * t399;
t305 = -t453 + t342 + (-pkin(3) - t506) * t397;
t505 = qJ(2) * t402;
t372 = t397 * t505;
t479 = t399 * t353 + t372;
t496 = t396 * t399;
t312 = -pkin(7) * t496 + t479;
t484 = t398 * t305 + t401 * t312;
t472 = qJD(1) * t396;
t450 = t399 * t472;
t428 = t398 * t450;
t470 = qJD(1) * t402;
t449 = t396 * t470;
t319 = t401 * t449 - t428;
t313 = t319 * qJ(5);
t332 = t353 * qJD(1) + qJD(2);
t324 = t402 * t332;
t415 = -t397 * t506 - t453;
t301 = t415 * qJD(1) + t324;
t471 = qJD(1) * t397;
t373 = -qJD(3) + t471;
t287 = -pkin(3) * t373 + t301;
t302 = -pkin(7) * t450 + qJD(1) * t372 + t332 * t399;
t294 = t398 * t302;
t440 = t401 * t287 - t294;
t267 = -t313 + t440;
t457 = qJDD(1) * t397;
t371 = -qJDD(3) + t457;
t360 = -qJDD(4) + t371;
t465 = qJD(4) * t401;
t366 = -qJD(4) + t373;
t513 = pkin(3) * t366;
t520 = t398 * pkin(3) * t360 + t465 * t513;
t501 = qJDD(1) * pkin(1);
t519 = -t501 - t443;
t518 = -MDP(12) * t399 - MDP(13) * t402;
t354 = t360 * pkin(4);
t504 = qJ(5) * t281;
t517 = t504 - t354;
t395 = qJ(3) + qJ(4);
t383 = sin(t395);
t384 = cos(t395);
t494 = t397 * t400;
t325 = t383 * t494 + t384 * t403;
t493 = t397 * t403;
t327 = -t383 * t493 + t384 * t400;
t508 = g(3) * t396;
t516 = -g(1) * t327 + g(2) * t325 + t383 * t508;
t515 = t319 ^ 2;
t392 = t397 ^ 2;
t512 = pkin(7) * t396;
t507 = MDP(9) * t396;
t435 = t454 * t402;
t444 = t398 * t455;
t480 = t454 * t428;
t282 = (t444 + (qJD(1) * t435 + t456) * t401) * t396 - t480;
t503 = qJ(5) * t282;
t416 = qJD(1) * t346;
t316 = t396 * t416;
t502 = qJ(5) * t316;
t500 = t316 * t366;
t499 = t319 * t366;
t498 = (-qJ(5) - pkin(7) - pkin(6)) * t396;
t404 = qJD(1) ^ 2;
t497 = t391 * t404;
t492 = t398 * t399;
t490 = t399 * t400;
t489 = t399 * t403;
t488 = t400 * t402;
t296 = t401 * t302;
t487 = t402 * t403;
t264 = -pkin(4) * t366 + t267;
t486 = -t267 + t264;
t485 = t401 * t301 - t294;
t345 = t401 * t402 - t492;
t483 = (t454 - t471) * t345;
t482 = t397 * t416 - t523;
t467 = qJD(3) * t402;
t469 = qJD(2) * t397;
t481 = t353 * t467 + t402 * t469;
t333 = pkin(3) * t450 + qJ(2) * t472;
t352 = t402 * pkin(3) + pkin(4) * t384;
t369 = t396 * pkin(3) * t467;
t339 = t396 * qJD(2) + t369;
t374 = pkin(3) * t496;
t344 = t396 * qJ(2) + t374;
t478 = t403 * pkin(1) + t400 * qJ(2);
t477 = t391 + t392;
t394 = t402 ^ 2;
t476 = t399 ^ 2 - t394;
t475 = MDP(10) * t396;
t468 = qJD(3) * t332;
t466 = qJD(4) * t398;
t350 = t360 * MDP(18);
t464 = t371 * MDP(11);
t463 = qJD(3) + t373;
t442 = -pkin(4) * t316 - qJD(5);
t299 = -t442 + t333;
t462 = qJD(5) + t299;
t459 = qJD(1) * qJD(3);
t458 = qJDD(1) * t396;
t451 = qJ(2) * qJD(3) * t397;
t448 = qJ(2) * t457;
t447 = t402 * t460;
t446 = t402 * t459;
t331 = t353 * qJDD(1) + qJDD(2);
t323 = t402 * t331;
t427 = qJD(1) * t451;
t276 = -pkin(3) * t371 + t323 + (-pkin(7) * t458 - t427) * t402 + (-t448 - t468 + (qJD(3) * t512 - t469) * qJD(1)) * t399;
t430 = t399 * t331 + t332 * t467 + t397 * t447 + t402 * t448;
t412 = -t399 * t427 + t430;
t280 = (-t446 - t456) * t512 + t412;
t441 = t401 * t276 - t398 * t280;
t297 = t415 * qJD(3) + t481;
t298 = -t399 * t469 + (-t372 + (-t353 + t512) * t399) * qJD(3);
t439 = -t297 * t398 + t401 * t298;
t438 = -t301 * t398 - t296;
t437 = t401 * t305 - t312 * t398;
t436 = qJD(1) * t463;
t434 = pkin(3) * t449;
t433 = t371 + t457;
t432 = MDP(22) * t454;
t431 = -t398 * t276 - t401 * t280 - t287 * t465 + t302 * t466;
t308 = qJ(2) * t458 + qJD(1) * t369 + qJDD(1) * t374 + t396 * t460;
t429 = 0.2e1 * t477;
t426 = -g(1) * t325 - g(2) * t327;
t326 = t383 * t403 - t384 * t494;
t328 = t383 * t400 + t384 * t493;
t425 = -g(1) * t326 - g(2) * t328;
t424 = g(1) * t403 + g(2) * t400;
t422 = -t287 * t398 - t296;
t420 = qJD(3) * (t373 + t471);
t418 = t429 * t460;
t417 = t401 * t297 + t398 * t298 + t305 * t465 - t312 * t466;
t271 = pkin(4) * t282 + qJDD(5) + t308;
t314 = t316 ^ 2;
t414 = t319 * t316 * MDP(14) + (-t281 - t500) * MDP(16) + (-t499 + (-t444 + (-t454 * t470 - t456) * t401) * t396 + t480) * MDP(17) + (-t314 + t515) * MDP(15) - t350;
t413 = g(1) * t328 - g(2) * t326 + t384 * t508 + t431;
t411 = t422 * qJD(4) + t441;
t410 = t333 * t316 + t413;
t408 = t462 * t316 + t413 + t503;
t407 = t411 + t516;
t406 = -t333 * t319 + t407;
t386 = t403 * qJ(2);
t380 = qJDD(2) - t501;
t379 = pkin(3) * t401 + pkin(4);
t351 = pkin(3) * t399 + pkin(4) * t383;
t347 = pkin(2) + t352;
t337 = t397 * t487 + t490;
t336 = -t397 * t489 + t488;
t335 = -t397 * t488 + t489;
t334 = t397 * t490 + t487;
t330 = t345 * t396;
t329 = t346 * t396;
t309 = pkin(4) * t319 + t434;
t306 = pkin(4) * t329 + t344;
t293 = -t454 * t396 * t492 + t435 * t522;
t292 = t523 * t396;
t284 = pkin(4) * t293 + t339;
t279 = -qJ(5) * t329 + t484;
t277 = -pkin(4) * t397 - qJ(5) * t330 + t437;
t270 = -t313 + t485;
t269 = t438 + t502;
t268 = -t422 - t502;
t263 = qJ(5) * t292 - qJD(4) * t484 - qJD(5) * t330 + t439;
t262 = -qJ(5) * t293 - qJD(5) * t329 + t417;
t261 = -qJD(5) * t316 - t431 - t503;
t260 = -qJD(5) * t319 + t411 + t517;
t1 = [qJDD(1) * MDP(1) + t443 * MDP(2) + t424 * MDP(3) + (-t380 - t519) * t397 * MDP(4) + (t429 * t461 + t418 - t424) * MDP(5) + (-t380 * pkin(1) - g(1) * (-pkin(1) * t400 + t386) - g(2) * t478 + (t477 * t461 + t418) * qJ(2)) * MDP(6) + (qJDD(1) * t394 - 0.2e1 * t399 * t446) * t391 * MDP(7) + (-t399 * t455 + t476 * t459) * MDP(8) * t514 + (t399 * t420 - t433 * t402) * t507 + (t433 * t399 + t402 * t420) * t475 + t397 * t464 + (-g(1) * t335 - g(2) * t337 - t323 * t397 - t342 * t371 + (t373 * t397 + (t514 + t392) * qJD(1)) * qJ(2) * t467 + (qJD(3) * t353 * t373 + t419 * t514 + (qJ(2) * t371 + qJD(2) * t373 + t468 + t521) * t397) * t399) * MDP(12) + ((-t399 * t451 + t481) * t373 + t479 * t371 + t412 * t397 - g(1) * t334 - g(2) * t336 + (t447 + (-t399 * t459 + t455) * qJ(2)) * t514) * MDP(13) + (-t281 * t330 - t292 * t319) * MDP(14) + (t281 * t329 - t282 * t330 + t292 * t316 - t293 * t319) * MDP(15) + (t281 * t397 + t292 * t366 - t330 * t360) * MDP(16) + (t282 * t397 + t293 * t366 + t329 * t360) * MDP(17) + t397 * t350 + (-t439 * t366 - t437 * t360 - t441 * t397 + t339 * t316 + t344 * t282 + t308 * t329 + t333 * t293 + (t366 * t484 - t422 * t397) * qJD(4) + t425) * MDP(19) + (-t344 * t281 - t333 * t292 + t308 * t330 + t339 * t319 + t484 * t360 + t417 * t366 - t431 * t397 + t426) * MDP(20) + (-t260 * t397 - t263 * t366 + t271 * t329 - t277 * t360 + t282 * t306 + t284 * t316 + t293 * t299 + t425) * MDP(21) + (t261 * t397 + t262 * t366 + t271 * t330 + t279 * t360 - t281 * t306 + t284 * t319 - t292 * t299 + t426) * MDP(22) + (-t260 * t330 - t261 * t329 - t262 * t316 - t263 * t319 + t264 * t292 - t268 * t293 + t277 * t281 - t279 * t282 + t443 * t396) * MDP(23) + (t261 * t279 + t268 * t262 + t260 * t277 + t264 * t263 + t271 * t306 + t299 * t284 - g(1) * (t351 * t403 + t386) - g(2) * (t347 * t493 - t403 * t498 + t478) + (-g(1) * (-t347 * t397 - pkin(1) + t498) - g(2) * t351) * t400) * MDP(24); -MDP(4) * t457 + (qJDD(2) + t519) * MDP(6) + (t281 * t345 - t282 * t346 - t483 * t316 - t482 * t319) * MDP(23) + (t260 * t345 + t261 * t346 + t482 * t264 + t483 * t268 - t299 * t472 - t443) * MDP(24) + (-MDP(12) * t402 + MDP(13) * t399) * t371 + (MDP(19) + MDP(21)) * (-t316 * t472 - t345 * t360 - t482 * t366) + (MDP(20) + MDP(22)) * (-t319 * t472 + t346 * t360 + t483 * t366) + (-t392 * MDP(5) + (-MDP(5) + t518) * t391 - t477 * MDP(6) * qJ(2)) * t404 + t518 * t373 ^ 2; t402 * t399 * MDP(7) * t497 - t476 * MDP(8) * t497 + (-t399 * t436 + t455) * t507 + (-t402 * t436 - t456) * t475 - t464 + (-g(1) * t336 + g(2) * t334 + t323 + (-t397 * t436 - t497) * t505 + (-t463 * t332 + t508 - t521) * t399) * MDP(12) + (g(3) * t495 + g(1) * t337 - g(2) * t335 - t324 * t373 + (t463 * t471 + t497) * t506 - t430) * MDP(13) + (t438 * t366 + (-t316 * t449 - t360 * t401 + t366 * t466) * pkin(3) + t406) * MDP(19) + (-t319 * t434 - t485 * t366 + t410 + t520) * MDP(20) + (t269 * t366 - t309 * t316 - t360 * t379 - t462 * t319 + (-t296 + (-t287 + t513) * t398) * qJD(4) + t441 + t516 + t517) * MDP(21) + (-t270 * t366 - t309 * t319 + t408 + t520) * MDP(22) + (t281 * t379 + (t268 + t269) * t319 + (-t264 + t270) * t316 + (-t282 * t398 + (-t316 * t401 + t319 * t398) * qJD(4)) * pkin(3)) * MDP(23) + (t260 * t379 - t268 * t270 - t264 * t269 - t299 * t309 - g(1) * (-t351 * t493 + t352 * t400) - g(2) * (-t351 * t494 - t352 * t403) + t351 * t508 + (t261 * t398 + (-t264 * t398 + t268 * t401) * qJD(4)) * pkin(3)) * MDP(24) + t414; (t422 * t366 + t406) * MDP(19) + (-t440 * t366 + t410) * MDP(20) + (t504 - t268 * t366 - 0.2e1 * t354 + (-t299 + t442) * t319 + t407) * MDP(21) + (-pkin(4) * t515 - t267 * t366 + t408) * MDP(22) + (pkin(4) * t281 - t486 * t316) * MDP(23) + (t486 * t268 + (-t299 * t319 + t260 + t516) * pkin(4)) * MDP(24) + t414; (-t480 - t499) * MDP(21) + (t500 + t355) * MDP(22) + (-t314 - t515) * MDP(23) + (g(3) * t397 + t264 * t319 + t268 * t316 + t271) * MDP(24) + (-t424 * MDP(24) + (t346 * MDP(21) - MDP(22) * t492) * qJDD(1) + (-t432 * t491 + (MDP(21) * t435 - t399 * t432) * t401) * qJD(1)) * t396;];
tau = t1;
