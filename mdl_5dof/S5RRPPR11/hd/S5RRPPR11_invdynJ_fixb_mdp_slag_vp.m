% Calculate vector of inverse dynamics joint torques for
% S5RRPPR11
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
%   see S5RRPPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:08
% EndTime: 2019-12-31 19:48:14
% DurationCPUTime: 4.61s
% Computational Cost: add. (2017->426), mult. (4278->546), div. (0->0), fcn. (2689->10), ass. (0->200)
t511 = pkin(3) + pkin(6);
t408 = sin(qJ(2));
t484 = qJD(1) * t408;
t374 = qJD(5) + t484;
t405 = cos(pkin(8));
t404 = sin(pkin(8));
t482 = qJD(2) * t404;
t411 = cos(qJ(2));
t483 = qJD(1) * t411;
t343 = t405 * t483 + t482;
t481 = qJD(2) * t405;
t345 = -t404 * t483 + t481;
t407 = sin(qJ(5));
t410 = cos(qJ(5));
t437 = t343 * t407 - t345 * t410;
t516 = t437 * t374;
t348 = t404 * t410 + t405 * t407;
t335 = t348 * qJD(5);
t425 = t348 * t408;
t488 = -qJD(1) * t425 - t335;
t472 = qJD(1) * qJD(2);
t459 = t411 * t472;
t471 = qJDD(1) * t408;
t424 = t459 + t471;
t347 = qJDD(5) + t424;
t462 = t405 * t484;
t476 = qJD(5) * t410;
t477 = qJD(5) * t407;
t495 = t404 * t407;
t487 = -t404 * t477 + t405 * t476 + t410 * t462 - t484 * t495;
t515 = -t347 * t348 - t487 * t374;
t409 = sin(qJ(1));
t412 = cos(qJ(1));
t447 = g(1) * t412 + g(2) * t409;
t383 = pkin(6) * t484;
t513 = qJD(3) + t383;
t384 = pkin(6) * t483;
t355 = pkin(3) * t483 + t384;
t401 = qJD(2) * qJ(3);
t333 = qJD(4) + t401 + t355;
t479 = qJD(2) * t411;
t505 = pkin(2) + qJ(4);
t512 = t505 * t479 + (qJD(4) - t333) * t408;
t510 = g(1) * t409;
t507 = g(2) * t412;
t506 = g(3) * t408;
t397 = g(3) * t411;
t394 = t411 * pkin(2);
t504 = -pkin(7) - t505;
t503 = pkin(6) * qJDD(2);
t502 = qJ(3) * t411;
t501 = qJ(4) * t411;
t500 = qJDD(2) * pkin(2);
t498 = t345 * t407;
t294 = t410 * t343 + t498;
t499 = t294 * t374;
t402 = t408 ^ 2;
t414 = qJD(1) ^ 2;
t496 = t402 * t414;
t494 = t405 * t411;
t390 = t408 * qJ(3);
t493 = t408 * t409;
t492 = t408 * t412;
t491 = t408 * t414;
t490 = t409 * t411;
t489 = t411 * t412;
t460 = t408 * t472;
t373 = pkin(2) * t460;
t441 = qJ(4) * t408 - t502;
t478 = qJD(3) * t408;
t415 = t441 * qJD(2) - qJD(4) * t411 - t478;
t456 = -pkin(1) - t390;
t423 = -t505 * t411 + t456;
t276 = t415 * qJD(1) + t423 * qJDD(1) + t373;
t372 = pkin(6) * t459;
t380 = pkin(6) * t471;
t457 = qJDD(3) + t372 + t380;
t291 = t424 * pkin(3) - qJD(2) * qJD(4) - t505 * qJDD(2) + t457;
t267 = t405 * t276 + t404 * t291;
t480 = qJD(2) * t408;
t386 = pkin(2) * t480;
t301 = t386 + t415;
t356 = t511 * t479;
t281 = t405 * t301 + t404 * t356;
t315 = t423 * qJD(1);
t474 = pkin(3) * t484 + t513;
t324 = -t505 * qJD(2) + t474;
t279 = t405 * t315 + t404 * t324;
t387 = pkin(2) * t484;
t329 = t441 * qJD(1) + t387;
t290 = t405 * t329 + t404 * t355;
t486 = t394 + t390;
t445 = t486 + t501;
t342 = -pkin(1) - t445;
t363 = t511 * t408;
t298 = t405 * t342 + t404 * t363;
t364 = t511 * t411;
t403 = t411 ^ 2;
t485 = t402 - t403;
t465 = -pkin(4) * t405 - pkin(3);
t475 = -t465 * t484 + t513;
t470 = qJDD(1) * t411;
t469 = t411 * t491;
t310 = qJDD(2) * t404 + (-t460 + t470) * t405;
t452 = qJDD(2) * t405 - t404 * t470;
t311 = t404 * t460 + t452;
t468 = -t407 * t310 + t410 * t311 - t343 * t476;
t467 = g(1) * t492 + g(2) * t493 - t397;
t381 = pkin(6) * t470;
t399 = qJDD(2) * qJ(3);
t400 = qJD(2) * qJD(3);
t466 = t381 + t399 + t400;
t464 = t511 * qJD(2);
t463 = t345 * t484;
t458 = t505 * t471;
t455 = -qJD(2) * pkin(2) + qJD(3);
t266 = -t276 * t404 + t405 * t291;
t264 = t424 * pkin(4) - pkin(7) * t311 + t266;
t265 = -pkin(7) * t310 + t267;
t454 = t410 * t264 - t265 * t407;
t280 = -t301 * t404 + t405 * t356;
t453 = t410 * t310 + t407 * t311;
t278 = -t315 * t404 + t405 * t324;
t289 = -t329 * t404 + t405 * t355;
t451 = t412 * pkin(1) + pkin(2) * t489 + t409 * pkin(6) + qJ(3) * t492;
t450 = -t380 + t467;
t449 = qJD(1) * t464;
t413 = qJD(2) ^ 2;
t448 = pkin(6) * t413 + t507;
t446 = -t507 + t510;
t435 = -t405 * t410 + t495;
t444 = -t435 * t347 + t374 * t488;
t443 = t279 * t484 + t266;
t442 = -t278 * t484 + t267;
t440 = t264 * t407 + t265 * t410;
t270 = pkin(4) * t484 - pkin(7) * t345 + t278;
t271 = -pkin(7) * t343 + t279;
t262 = t270 * t410 - t271 * t407;
t263 = t270 * t407 + t271 * t410;
t351 = t405 * t363;
t283 = pkin(4) * t408 + t351 + (pkin(7) * t411 - t342) * t404;
t288 = -pkin(7) * t494 + t298;
t439 = t283 * t410 - t288 * t407;
t438 = t283 * t407 + t288 * t410;
t359 = t383 + t455;
t362 = -t384 - t401;
t436 = t359 * t411 + t362 * t408;
t434 = pkin(3) * t470 + qJDD(4) + t466;
t433 = -pkin(7) * t404 * t408 + pkin(4) * t411;
t432 = t456 - t394;
t431 = t447 * t411;
t430 = -0.2e1 * pkin(1) * t472 - t503;
t341 = t432 * qJD(1);
t429 = t341 * t484 + qJDD(3) - t450;
t428 = -qJ(3) * t479 - t478;
t357 = t504 * t404;
t427 = t433 * qJD(1) + qJD(4) * t405 + qJD(5) * t357 + t289;
t358 = t504 * t405;
t426 = pkin(7) * t462 + qJD(4) * t404 - qJD(5) * t358 + t290;
t325 = t435 * t411;
t268 = -t345 * t477 + t468;
t422 = 0.2e1 * qJDD(1) * pkin(1) - t448;
t293 = -t408 * t449 + t434;
t421 = -t293 * t411 + t447;
t360 = -pkin(1) - t486;
t420 = t503 + (-qJD(1) * t360 - t341) * qJD(2);
t419 = -t431 - t506;
t418 = t293 + t419;
t269 = -t437 * qJD(5) + t453;
t292 = t428 * qJD(1) + t432 * qJDD(1) + t373;
t331 = t386 + t428;
t417 = qJD(1) * t331 + qJDD(1) * t360 + t292 + t448;
t316 = pkin(6) * t460 - t466;
t328 = t457 - t500;
t416 = t436 * qJD(2) - t316 * t411 + t328 * t408;
t398 = pkin(8) + qJ(5);
t395 = t412 * pkin(6);
t389 = cos(t398);
t388 = sin(t398);
t379 = pkin(4) * t404 + qJ(3);
t377 = g(1) * t490;
t371 = qJ(3) * t489;
t369 = qJ(3) * t490;
t354 = t408 * t464;
t352 = -qJ(3) * t483 + t387;
t334 = pkin(4) * t494 + t364;
t326 = t348 * t411;
t323 = -t388 * t493 + t389 * t412;
t322 = t388 * t412 + t389 * t493;
t321 = t388 * t492 + t389 * t409;
t320 = -t388 * t409 + t389 * t492;
t318 = (-pkin(6) + t465) * t480;
t299 = pkin(4) * t343 + t333;
t297 = -t342 * t404 + t351;
t285 = t335 * t411 - t435 * t480;
t284 = qJD(2) * t425 + qJD(5) * t325;
t275 = pkin(7) * t405 * t480 + t281;
t274 = pkin(4) * t310 + t293;
t272 = t433 * qJD(2) + t280;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t402 + 0.2e1 * t408 * t459) * MDP(4) + 0.2e1 * (t408 * t470 - t485 * t472) * MDP(5) + (qJDD(2) * t408 + t411 * t413) * MDP(6) + (qJDD(2) * t411 - t408 * t413) * MDP(7) + (t430 * t408 + t422 * t411 + t377) * MDP(9) + (t430 * t411 + (-t422 - t510) * t408) * MDP(10) + ((t402 + t403) * qJDD(1) * pkin(6) + t416 - t447) * MDP(11) + (t420 * t408 + t417 * t411 - t377) * MDP(12) + (t420 * t411 + (-t417 + t510) * t408) * MDP(13) + (t416 * pkin(6) - g(1) * t395 - g(2) * t451 + t292 * t360 + t341 * t331 - t432 * t510) * MDP(14) + (t364 * t310 - t354 * t343 + (qJD(1) * t297 + t278) * t479 - t421 * t405 + (t280 * qJD(1) + t297 * qJDD(1) - t333 * t481 + t446 * t404 + t266) * t408) * MDP(15) + (t364 * t311 - t354 * t345 + (-qJD(1) * t298 - t279) * t479 + t421 * t404 + (-t281 * qJD(1) - t298 * qJDD(1) + t333 * t482 + t446 * t405 - t267) * t408) * MDP(16) + (-t280 * t345 - t281 * t343 - t297 * t311 - t298 * t310 + t377 + (-t278 * t404 + t279 * t405) * t480 + (t266 * t404 - t267 * t405 - t507) * t411) * MDP(17) + (t267 * t298 + t279 * t281 + t266 * t297 + t278 * t280 + t293 * t364 - t333 * t354 - g(1) * (pkin(3) * t412 + t395) - g(2) * (qJ(4) * t489 + t451) + (-g(1) * (t432 - t501) - g(2) * pkin(3)) * t409) * MDP(18) + (-t268 * t326 - t284 * t437) * MDP(19) + (t268 * t325 + t269 * t326 - t284 * t294 - t285 * t437) * MDP(20) + (t268 * t408 + t284 * t374 - t326 * t347 - t437 * t479) * MDP(21) + (-t269 * t408 + t285 * t374 - t294 * t479 + t325 * t347) * MDP(22) + (t347 * t408 + t374 * t479) * MDP(23) + ((t272 * t410 - t275 * t407) * t374 + t439 * t347 + t454 * t408 + t262 * t479 + t318 * t294 + t334 * t269 - t274 * t325 - t299 * t285 - g(1) * t323 - g(2) * t321 + (-t263 * t408 - t438 * t374) * qJD(5)) * MDP(24) + (-(t272 * t407 + t275 * t410) * t374 - t438 * t347 - t440 * t408 - t263 * t479 - t318 * t437 + t334 * t268 - t274 * t326 + t299 * t284 + g(1) * t322 - g(2) * t320 + (-t262 * t408 - t439 * t374) * qJD(5)) * MDP(25) + t446 * MDP(2) + t447 * MDP(3); -MDP(4) * t469 + t485 * MDP(5) * t414 + MDP(6) * t471 + MDP(7) * t470 + qJDD(2) * MDP(8) + (pkin(1) * t491 + t450) * MDP(9) + (t506 - t381 + (pkin(1) * t414 + t447) * t411) * MDP(10) + ((-pkin(2) * t408 + t502) * qJDD(1) + ((-t362 - t401) * t408 + (-t359 + t455) * t411) * qJD(1)) * MDP(11) + (-t352 * t483 + t429 - 0.2e1 * t500) * MDP(12) + (t381 + 0.2e1 * t399 + 0.2e1 * t400 + (qJD(1) * t352 - g(3)) * t408 + (qJD(1) * t341 - t447) * t411) * MDP(13) + (-t316 * qJ(3) - t362 * qJD(3) - t328 * pkin(2) - t341 * t352 - g(1) * (-pkin(2) * t492 + t371) - g(2) * (-pkin(2) * t493 + t369) - g(3) * t486 - t436 * qJD(1) * pkin(6)) * MDP(14) + (-t405 * t458 + qJ(3) * t310 + t474 * t343 + t418 * t404 + (-t278 * t411 - t289 * t408 - t405 * t512) * qJD(1)) * MDP(15) + (t404 * t458 + qJ(3) * t311 + t474 * t345 + t418 * t405 + (t279 * t411 + t290 * t408 + t404 * t512) * qJD(1)) * MDP(16) + (t289 * t345 + t290 * t343 + (qJD(4) * t345 + t311 * t505 - t443) * t405 + (qJD(4) * t343 + t310 * t505 - t442) * t404 + t467) * MDP(17) + (t293 * qJ(3) - t279 * t290 - t278 * t289 - g(1) * t371 - g(2) * t369 - g(3) * t445 + t474 * t333 + (-t278 * t405 - t279 * t404) * qJD(4) + (-t266 * t405 - t267 * t404 + t447 * t408) * t505) * MDP(18) + (-t268 * t435 - t437 * t488) * MDP(19) + (-t268 * t348 + t269 * t435 - t488 * t294 + t437 * t487) * MDP(20) + (t437 * t483 + t444) * MDP(21) + (t294 * t483 + t515) * MDP(22) - t374 * MDP(23) * t483 + ((-t357 * t407 + t358 * t410) * t347 + t379 * t269 + t274 * t348 - t262 * t483 + (t426 * t407 - t427 * t410) * t374 + t487 * t299 + t475 * t294 + t419 * t388) * MDP(24) + (-(t357 * t410 + t358 * t407) * t347 + t379 * t268 - t274 * t435 + t263 * t483 + (t427 * t407 + t426 * t410) * t374 + t488 * t299 - t475 * t437 + t419 * t389) * MDP(25); MDP(11) * t471 + (qJDD(2) + t469) * MDP(12) + (-t413 - t496) * MDP(13) + (t372 + t429 - t500) * MDP(14) - t467 * MDP(18) + t444 * MDP(24) + t515 * MDP(25) + (MDP(14) * t362 - MDP(15) * t343 - MDP(16) * t345 - MDP(18) * t333 - MDP(24) * t294 + MDP(25) * t437) * qJD(2) + (t424 * MDP(15) - MDP(16) * t496 + (-t343 * t484 - t311) * MDP(17) + t443 * MDP(18)) * t405 + (-MDP(15) * t496 - t424 * MDP(16) + (-t310 + t463) * MDP(17) + t442 * MDP(18)) * t404; (t310 + t463) * MDP(15) + ((-t343 + t482) * t484 + t452) * MDP(16) + (-t343 ^ 2 - t345 ^ 2) * MDP(17) + (t278 * t345 + t279 * t343 - t431 + (-g(3) - t449) * t408 + t434) * MDP(18) + (t269 - t516) * MDP(24) + (t268 - t499) * MDP(25); -t437 * t294 * MDP(19) + (-t294 ^ 2 + t437 ^ 2) * MDP(20) + (t468 + t499) * MDP(21) + (-t453 - t516) * MDP(22) + t347 * MDP(23) + (-g(1) * t320 - g(2) * t322 + t263 * t374 + t299 * t437 + t389 * t397 + t454) * MDP(24) + (g(1) * t321 - g(2) * t323 + t262 * t374 + t294 * t299 - t388 * t397 - t440) * MDP(25) + (-MDP(21) * t498 + t437 * MDP(22) - t263 * MDP(24) - t262 * MDP(25)) * qJD(5);];
tau = t1;
