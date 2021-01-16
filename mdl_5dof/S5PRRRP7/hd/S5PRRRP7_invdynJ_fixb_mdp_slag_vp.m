% Calculate vector of inverse dynamics joint torques for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:38
% EndTime: 2021-01-15 16:45:51
% DurationCPUTime: 5.16s
% Computational Cost: add. (2280->449), mult. (5315->597), div. (0->0), fcn. (4021->10), ass. (0->194)
t379 = sin(pkin(5));
t383 = sin(qJ(4));
t385 = sin(qJ(2));
t386 = cos(qJ(4));
t387 = cos(qJ(3));
t388 = cos(qJ(2));
t469 = t387 * t388;
t312 = (t383 * t385 + t386 * t469) * t379;
t384 = sin(qJ(3));
t415 = pkin(3) * t384 - pkin(8) * t387;
t349 = t415 * qJD(3);
t353 = -pkin(3) * t387 - pkin(8) * t384 - pkin(2);
t445 = qJD(4) * t386;
t510 = -qJD(1) * t312 + t383 * t349 + t353 * t445;
t441 = qJD(2) * qJD(3);
t425 = t387 * t441;
t438 = qJDD(2) * t384;
t401 = -t425 - t438;
t509 = qJD(3) * qJD(4) - t401;
t311 = (-t383 * t469 + t385 * t386) * t379;
t449 = qJD(3) * t384;
t500 = pkin(7) * t383;
t508 = -qJD(1) * t311 + t386 * t349 + t449 * t500;
t446 = qJD(4) * t384;
t424 = qJD(2) * t446;
t408 = (-qJDD(3) + t424) * t386;
t453 = qJD(2) * t387;
t279 = ((qJD(4) + t453) * qJD(3) + t438) * t383 + t408;
t470 = t386 * t387;
t367 = pkin(7) * t470;
t409 = pkin(4) * t384 - qJ(5) * t470;
t444 = qJD(5) * t386;
t492 = qJ(5) * t384;
t507 = -t384 * t444 + t409 * qJD(3) + (-t367 + (-t353 + t492) * t383) * qJD(4) + t508;
t473 = t384 * t386;
t466 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t473 + (-qJD(5) * t384 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t387) * t383 + t510;
t457 = qJD(1) * t379;
t351 = qJD(2) * pkin(7) + t385 * t457;
t381 = cos(pkin(5));
t456 = qJD(1) * t381;
t506 = -t384 * t351 + t387 * t456;
t389 = qJD(3) ^ 2;
t442 = qJD(1) * qJD(2);
t426 = t385 * t442;
t478 = t379 * t388;
t412 = -qJDD(1) * t478 + t379 * t426;
t378 = sin(pkin(9));
t380 = cos(pkin(9));
t475 = t381 * t388;
t325 = t378 * t385 - t380 * t475;
t327 = t378 * t475 + t380 * t385;
t413 = g(1) * t327 + g(2) * t325;
t505 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t389 + (-g(3) * t388 + t426) * t379 - t412 + t413;
t476 = t381 * t385;
t324 = t378 * t476 - t380 * t388;
t480 = t379 * t384;
t283 = t324 * t387 - t378 * t480;
t326 = t378 * t388 + t380 * t476;
t285 = t326 * t387 - t380 * t480;
t316 = t325 * t386;
t319 = t327 * t386;
t471 = t385 * t387;
t331 = t379 * t471 + t381 * t384;
t294 = -t331 * t383 - t386 * t478;
t499 = g(3) * t294;
t503 = -t499 - g(1) * (t283 * t383 + t319) - g(2) * (-t285 * t383 + t316);
t450 = qJD(3) * t383;
t454 = qJD(2) * t384;
t346 = t386 * t454 + t450;
t502 = t346 ^ 2;
t443 = t386 * qJD(3);
t344 = t383 * t454 - t443;
t501 = pkin(4) * t344;
t406 = -t331 * t386 + t383 * t478;
t498 = g(3) * t406;
t497 = g(3) * t311;
t496 = g(3) * t312;
t474 = t384 * t385;
t330 = t379 * t474 - t381 * t387;
t495 = g(3) * t330;
t437 = t387 * qJDD(2);
t400 = -t384 * t441 + t437;
t341 = qJDD(4) - t400;
t494 = t341 * pkin(4);
t382 = -qJ(5) - pkin(8);
t493 = qJD(2) * pkin(2);
t364 = t384 * t456;
t306 = t387 * t351 + t364;
t300 = qJD(3) * pkin(8) + t306;
t435 = t388 * t457;
t308 = t353 * qJD(2) - t435;
t272 = t300 * t386 + t308 * t383;
t266 = -qJ(5) * t344 + t272;
t366 = -qJD(4) + t453;
t490 = t266 * t366;
t418 = t383 * qJDD(3) + t509 * t386;
t278 = t383 * t424 - t418;
t489 = t278 * t383;
t320 = t324 * t383;
t486 = t325 * t383;
t317 = t326 * t383;
t485 = t327 * t383;
t484 = t344 * t366;
t483 = t346 * t366;
t482 = t346 * t386;
t481 = t378 * t381;
t479 = t379 * t387;
t477 = t380 * t381;
t472 = t384 * t388;
t468 = qJDD(1) - g(3);
t271 = -t300 * t383 + t386 * t308;
t265 = -qJ(5) * t346 + t271;
t260 = -pkin(4) * t366 + t265;
t465 = -t265 + t260;
t348 = t415 * qJD(2);
t464 = t383 * t348 + t386 * t506;
t422 = qJD(4) * t382;
t431 = t383 * t453;
t463 = qJ(5) * t431 + t383 * t422 + t444 - t464;
t335 = t386 * t348;
t462 = -t409 * qJD(2) + t386 * t422 - t335 + (-qJD(5) + t506) * t383;
t459 = t383 * t353 + t367;
t376 = t384 ^ 2;
t458 = -t387 ^ 2 + t376;
t455 = qJD(2) * t379;
t452 = qJD(3) * t344;
t451 = qJD(3) * t346;
t448 = qJD(3) * t387;
t447 = qJD(4) * t383;
t439 = qJDD(1) * t381;
t436 = pkin(4) * t383 + pkin(7);
t433 = t385 * t455;
t432 = t388 * t455;
t430 = t366 * t443;
t299 = -qJD(3) * pkin(3) - t506;
t421 = -qJD(5) - t501;
t277 = t299 - t421;
t429 = t277 * t445;
t428 = t366 * t447;
t427 = t366 * t445;
t423 = t384 * t439;
t314 = qJDD(2) * pkin(7) + (qJDD(1) * t385 + t388 * t442) * t379;
t268 = qJDD(3) * pkin(8) + qJD(3) * t506 + t314 * t387 + t423;
t276 = qJD(2) * t349 + t353 * qJDD(2) + t412;
t419 = -t386 * t268 - t383 * t276 + t300 * t447 - t308 * t445;
t417 = t344 * t435;
t416 = t346 * t435;
t284 = t326 * t384 + t380 * t479;
t286 = t324 * t384 + t378 * t479;
t414 = g(1) * t286 - g(2) * t284;
t371 = pkin(4) * t386 + pkin(3);
t411 = -t371 * t387 + t382 * t384;
t397 = qJD(3) * t364 + t384 * t314 + t351 * t448 - t387 * t439;
t269 = -qJDD(3) * pkin(3) + t397;
t259 = t279 * pkin(4) + qJDD(5) + t269;
t410 = -t259 - t414;
t404 = t341 * t383 - t427;
t403 = t341 * t386 + t428;
t402 = t279 * qJ(5) + t419;
t399 = -g(1) * t283 + g(2) * t285 + g(3) * t331;
t398 = -t414 + t495;
t396 = -pkin(8) * t341 - t366 * t299;
t329 = -t381 * t474 - t479;
t395 = -pkin(8) * qJD(4) * t366 - g(1) * (t329 * t378 + t380 * t472) - g(2) * (-t329 * t380 + t378 * t472) + t269;
t275 = t386 * t276;
t393 = -t272 * qJD(4) - t383 * t268 + t275;
t352 = -t435 - t493;
t392 = -pkin(7) * qJDD(3) + (t352 + t435 - t493) * qJD(3);
t391 = t278 * qJ(5) + t393;
t390 = qJD(2) ^ 2;
t355 = t382 * t386;
t354 = t382 * t383;
t350 = t436 * t384;
t343 = t386 * t353;
t340 = t344 ^ 2;
t332 = t381 * t471 - t480;
t321 = t324 * t386;
t318 = t326 * t386;
t315 = t383 * t495;
t310 = t327 * t387;
t309 = t325 * t387;
t307 = pkin(7) * t448 + (t383 * t448 + t384 * t445) * pkin(4);
t296 = -t383 * t492 + t459;
t293 = t331 * qJD(3) + t384 * t432;
t292 = -t330 * qJD(3) + t387 * t432;
t291 = -t332 * t378 + t380 * t469;
t289 = t332 * t380 + t378 * t469;
t281 = pkin(4) * t431 + t306;
t280 = -qJ(5) * t473 + t343 + (-pkin(4) - t500) * t387;
t263 = t294 * qJD(4) + t292 * t386 + t383 * t433;
t262 = t406 * qJD(4) - t292 * t383 + t386 * t433;
t256 = -t344 * qJD(5) - t402;
t255 = -t346 * qJD(5) + t391 + t494;
t1 = [t468 * MDP(1) + (-qJD(3) * t293 - qJDD(3) * t330) * MDP(10) + (-qJD(3) * t292 - qJDD(3) * t331) * MDP(11) + (-t262 * t346 - t263 * t344 + t278 * t294 + t279 * t406) * MDP(21) + (t255 * t294 - t256 * t406 + t259 * t330 + t260 * t262 + t263 * t266 + t277 * t293 - g(3)) * MDP(22) + (MDP(18) + MDP(20)) * (t263 * t366 - t278 * t330 + t293 * t346 + t341 * t406) + (MDP(17) + MDP(19)) * (-t262 * t366 + t279 * t330 + t293 * t344 + t294 * t341) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t387 + MDP(11) * t384 - MDP(3)) * t390) * t385 + (t400 * MDP(10) + t401 * MDP(11) + qJDD(2) * MDP(3) - t390 * MDP(4)) * t388) * t379; qJDD(2) * MDP(2) + (t468 * t478 + t413) * MDP(3) + (-t468 * t385 * t379 - g(1) * t324 + g(2) * t326) * MDP(4) + (qJDD(2) * t376 + 0.2e1 * t384 * t425) * MDP(5) + 0.2e1 * (t384 * t437 - t458 * t441) * MDP(6) + (qJDD(3) * t384 + t387 * t389) * MDP(7) + (qJDD(3) * t387 - t384 * t389) * MDP(8) + (t392 * t384 + t505 * t387) * MDP(10) + (-t505 * t384 + t392 * t387) * MDP(11) + (-t278 * t473 + (-t383 * t446 + t387 * t443) * t346) * MDP(12) + ((-t344 * t386 - t346 * t383) * t448 + (t489 - t279 * t386 + (t344 * t383 - t482) * qJD(4)) * t384) * MDP(13) + ((t278 - t430) * t387 + (t403 + t451) * t384) * MDP(14) + ((t366 * t450 + t279) * t387 + (-t404 - t452) * t384) * MDP(15) + (-t341 * t387 - t366 * t449) * MDP(16) + (t343 * t341 - g(1) * (-t310 * t386 - t320) - g(2) * (-t309 * t386 + t317) - t496 + (t353 * t447 - t508) * t366 + (t300 * t445 - t275 + (t427 + t452) * pkin(7) + (-pkin(7) * t341 + qJD(3) * t299 + qJD(4) * t308 + t268) * t383) * t387 + (pkin(7) * t279 + qJD(3) * t271 + t269 * t383 + t299 * t445 - t417) * t384) * MDP(17) + (-t459 * t341 - g(1) * (t310 * t383 - t321) - g(2) * (t309 * t383 + t318) - t497 + t510 * t366 + (t299 * t443 + (-t428 + t451) * pkin(7) - t419) * t387 + (-t416 - t299 * t447 - t272 * qJD(3) + t269 * t386 + (-t278 - t430) * pkin(7)) * t384) * MDP(18) + (g(1) * t320 - g(2) * t317 - t496 + t350 * t279 + t280 * t341 + t307 * t344 + (t277 * t450 + t413 * t386 - t255) * t387 - t507 * t366 + (qJD(3) * t260 + t259 * t383 - t417 + t429) * t384) * MDP(19) + (g(1) * t321 - g(2) * t318 - t497 - t350 * t278 - t296 * t341 + t307 * t346 + (t277 * t443 - t413 * t383 + t256) * t387 + t466 * t366 + (-qJD(3) * t266 + t259 * t386 - t277 * t447 - t416) * t384) * MDP(20) + (t278 * t280 - t279 * t296 - t507 * t346 - t466 * t344 + (-t260 * t386 - t266 * t383) * t448 + (-g(3) * t478 - t255 * t386 - t256 * t383 + (t260 * t383 - t266 * t386) * qJD(4) + t413) * t384) * MDP(21) + (t256 * t296 + t255 * t280 + t259 * t350 - g(1) * (-pkin(4) * t320 - (pkin(2) * t380 + pkin(7) * t481) * t385 + (-pkin(2) * t481 + pkin(7) * t380) * t388 + t411 * t327) - g(2) * (pkin(4) * t317 - (pkin(2) * t378 - pkin(7) * t477) * t385 + (pkin(2) * t477 + pkin(7) * t378) * t388 + t411 * t325) - g(3) * (t436 * t385 + (pkin(2) - t411) * t388) * t379 + (-t384 * t435 + t307) * t277 + t466 * t266 + t507 * t260) * MDP(22); MDP(7) * t438 + MDP(8) * t437 + qJDD(3) * MDP(9) + (qJD(3) * t306 - t352 * t454 - t397 + t398) * MDP(10) + (-t423 + (-qJD(2) * t352 - t314) * t387 + t399) * MDP(11) + (-t366 * t482 - t489) * MDP(12) + ((-t278 + t484) * t386 + (-t279 + t483) * t383) * MDP(13) + ((-t346 * t384 + t366 * t470) * qJD(2) + t404) * MDP(14) + ((-t366 * t383 * t387 + t344 * t384) * qJD(2) + t403) * MDP(15) + t366 * MDP(16) * t454 + (-t271 * t454 - pkin(3) * t279 - t306 * t344 + t335 * t366 + (-t366 * t506 + t396) * t383 + (-t395 + t495) * t386) * MDP(17) + (pkin(3) * t278 + t272 * t454 - t306 * t346 - t366 * t464 + t383 * t395 + t386 * t396 - t315) * MDP(18) + (-t260 * t454 - t279 * t371 - t281 * t344 + t341 * t354 - t462 * t366 + (-t277 * t453 + (t277 + t501) * qJD(4)) * t383 + (-t259 + t398) * t386) * MDP(19) + (t429 + t278 * t371 - t281 * t346 + t341 * t355 - t315 + t463 * t366 + (t266 * t384 - t277 * t470) * qJD(2) + (pkin(4) * qJD(4) * t346 - t410) * t383) * MDP(20) + (t278 * t354 + t279 * t355 - t462 * t346 - t463 * t344 + (t260 * t366 + t256) * t386 + (-t255 + t490) * t383 - t399) * MDP(21) + (-t256 * t355 + t255 * t354 - t259 * t371 - g(1) * (t283 * t382 + t286 * t371) - g(2) * (-t284 * t371 - t285 * t382) - g(3) * (-t330 * t371 - t331 * t382) + (pkin(4) * t447 - t281) * t277 + t463 * t266 + t462 * t260) * MDP(22) + (-t384 * t387 * MDP(5) + t458 * MDP(6)) * t390; t346 * t344 * MDP(12) + (-t340 + t502) * MDP(13) + (-t278 - t484) * MDP(14) + (-t279 - t483) * MDP(15) + t341 * MDP(16) + (-t272 * t366 - t299 * t346 - g(1) * (-t291 * t383 + t319) - g(2) * (-t289 * t383 + t316) - t499 + t393) * MDP(17) + (-t271 * t366 + t299 * t344 - g(1) * (-t291 * t386 - t485) - g(2) * (-t289 * t386 - t486) - t498 + t419) * MDP(18) + (0.2e1 * t494 - t490 + (-t277 + t421) * t346 + t391 + t503) * MDP(19) + (-t265 * t366 - t502 * pkin(4) - g(1) * (t283 * t386 - t485) - g(2) * (-t285 * t386 - t486) - t498 + (qJD(5) + t277) * t344 + t402) * MDP(20) + (pkin(4) * t278 - t344 * t465) * MDP(21) + (t465 * t266 + (-t277 * t346 + t255 + t503) * pkin(4)) * MDP(22); (t408 - t483) * MDP(19) + (t418 + t484) * MDP(20) + (-t340 - t502) * MDP(21) + (t260 * t346 + t266 * t344 - t410 - t495) * MDP(22) + (t509 * MDP(19) - MDP(20) * t424) * t383;];
tau = t1;
