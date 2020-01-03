% Calculate vector of inverse dynamics joint torques for
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:23
% EndTime: 2019-12-31 20:18:30
% DurationCPUTime: 4.64s
% Computational Cost: add. (3980->370), mult. (9560->487), div. (0->0), fcn. (7410->12), ass. (0->173)
t428 = sin(pkin(9));
t429 = cos(pkin(9));
t433 = sin(qJ(2));
t437 = cos(qJ(2));
t390 = -t428 * t433 + t429 * t437;
t380 = t390 * qJD(1);
t436 = cos(qJ(4));
t366 = t436 * t380;
t391 = t428 * t437 + t429 * t433;
t382 = t391 * qJD(1);
t432 = sin(qJ(4));
t338 = -t382 * t432 + t366;
t435 = cos(qJ(5));
t483 = qJD(5) * t435;
t536 = -t338 * t435 + t483;
t423 = qJ(2) + pkin(9) + qJ(4);
t416 = sin(t423);
t434 = sin(qJ(1));
t438 = cos(qJ(1));
t459 = g(1) * t438 + g(2) * t434;
t535 = t459 * t416;
t455 = t380 * t432 + t436 * t382;
t381 = t391 * qJD(2);
t347 = -qJD(1) * t381 + qJDD(1) * t390;
t479 = qJD(1) * qJD(2);
t472 = t437 * t479;
t473 = t433 * t479;
t348 = qJDD(1) * t391 - t428 * t473 + t429 * t472;
t485 = qJD(4) * t432;
t295 = qJD(4) * t366 + t432 * t347 + t436 * t348 - t382 * t485;
t424 = qJDD(2) + qJDD(4);
t425 = qJD(2) + qJD(4);
t431 = sin(qJ(5));
t474 = t435 * t295 + t431 * t424 + t425 * t483;
t484 = qJD(5) * t431;
t283 = -t455 * t484 + t474;
t282 = t283 * t435;
t324 = t425 * t431 + t435 * t455;
t413 = t435 * t424;
t284 = qJD(5) * t324 + t295 * t431 - t413;
t322 = -t435 * t425 + t431 * t455;
t534 = -t431 * t284 - t536 * t322 + t282;
t281 = t283 * t431;
t296 = qJD(4) * t455 - t436 * t347 + t348 * t432;
t294 = qJDD(5) + t296;
t288 = t431 * t294;
t528 = qJD(5) - t338;
t493 = t435 * t528;
t499 = t338 * t425;
t501 = t455 * t425;
t503 = t324 * t455;
t533 = t424 * MDP(17) + (-t296 + t501) * MDP(16) - t338 ^ 2 * MDP(14) + (-t338 * t493 + t483 * t528 + t288 - t503) * MDP(22) + (-t338 * MDP(13) + MDP(14) * t455 - MDP(24) * t528) * t455 + (t295 - t499) * MDP(15) + (t536 * t324 + t281) * MDP(20);
t508 = qJ(3) + pkin(6);
t406 = t508 * t437;
t398 = qJD(1) * t406;
t385 = t428 * t398;
t405 = t508 * t433;
t397 = qJD(1) * t405;
t507 = qJD(2) * pkin(2);
t389 = -t397 + t507;
t345 = t429 * t389 - t385;
t513 = pkin(7) * t382;
t317 = qJD(2) * pkin(3) + t345 - t513;
t497 = t429 * t398;
t346 = t428 * t389 + t497;
t514 = pkin(7) * t380;
t321 = t346 + t514;
t298 = t317 * t436 - t321 * t432;
t290 = -pkin(4) * t425 - t298;
t532 = t290 * t338;
t470 = qJD(2) * t508;
t378 = -qJD(3) * t433 - t437 * t470;
t344 = qJDD(2) * pkin(2) + qJD(1) * t378 - qJDD(1) * t405;
t377 = qJD(3) * t437 - t433 * t470;
t351 = qJD(1) * t377 + qJDD(1) * t406;
t310 = t429 * t344 - t351 * t428;
t297 = qJDD(2) * pkin(3) - pkin(7) * t348 + t310;
t299 = t317 * t432 + t321 * t436;
t311 = t428 * t344 + t429 * t351;
t300 = pkin(7) * t347 + t311;
t517 = qJD(4) * t299 - t436 * t297 + t432 * t300;
t273 = -pkin(4) * t424 + t517;
t417 = cos(t423);
t510 = g(3) * t417;
t530 = t273 + t510;
t420 = pkin(2) * t437 + pkin(1);
t399 = -qJD(1) * t420 + qJD(3);
t354 = -pkin(3) * t380 + t399;
t410 = g(3) * t416;
t518 = (qJD(4) * t317 + t300) * t436 + t432 * t297 - t321 * t485;
t527 = -t354 * t338 + t417 * t459 + t410 - t518;
t524 = pkin(4) * t455;
t504 = t322 * t455;
t523 = t528 * (pkin(8) * t528 + t524);
t418 = pkin(2) * t429 + pkin(3);
t515 = pkin(2) * t428;
t487 = t432 * t418 + t436 * t515;
t291 = pkin(8) * t425 + t299;
t301 = -pkin(4) * t338 - pkin(8) * t455 + t354;
t276 = -t291 * t431 + t301 * t435;
t522 = -t276 * t455 + t290 * t484 + t435 * t535;
t277 = t291 * t435 + t301 * t431;
t521 = t277 * t455 + t290 * t483 + t431 * t530;
t520 = -t354 * t455 - t510 - t517 + t535;
t328 = -t377 * t428 + t429 * t378;
t384 = t390 * qJD(2);
t314 = -pkin(7) * t384 + t328;
t329 = t429 * t377 + t428 * t378;
t315 = -pkin(7) * t381 + t329;
t355 = -t429 * t405 - t406 * t428;
t330 = -pkin(7) * t391 + t355;
t356 = -t428 * t405 + t429 * t406;
t331 = pkin(7) * t390 + t356;
t456 = t330 * t436 - t331 * t432;
t278 = qJD(4) * t456 + t314 * t432 + t315 * t436;
t350 = t390 * t432 + t391 * t436;
t360 = -pkin(3) * t390 - t420;
t454 = t436 * t390 - t391 * t432;
t305 = -pkin(4) * t454 - pkin(8) * t350 + t360;
t307 = t330 * t432 + t331 * t436;
t312 = qJD(4) * t454 - t381 * t432 + t384 * t436;
t272 = pkin(8) * t424 + t518;
t463 = qJD(5) * t301 + t272;
t516 = t273 * t350 + t290 * t312 - t307 * t294 - (qJD(5) * t305 + t278) * t528 + t454 * t463;
t509 = g(3) * t437;
t506 = t290 * t350;
t505 = t305 * t294;
t502 = t324 * t431;
t496 = t431 * t434;
t495 = t431 * t438;
t494 = t434 * t435;
t289 = t435 * t294;
t492 = t435 * t438;
t352 = t397 * t428 - t497;
t325 = t352 - t514;
t353 = -t429 * t397 - t385;
t326 = t353 - t513;
t452 = t418 * t436 - t432 * t515;
t489 = -t452 * qJD(4) + t325 * t432 + t326 * t436;
t488 = qJD(4) * t487 + t325 * t436 - t326 * t432;
t426 = t433 ^ 2;
t486 = -t437 ^ 2 + t426;
t478 = qJDD(1) * t437;
t422 = t433 * t507;
t358 = pkin(3) * t381 + t422;
t357 = t433 * qJD(1) * pkin(2) + pkin(3) * t382;
t465 = t528 * t431;
t447 = pkin(2) * t473 - qJDD(1) * t420 + qJDD(3);
t320 = -pkin(3) * t347 + t447;
t275 = pkin(4) * t296 - pkin(8) * t295 + t320;
t462 = qJD(5) * t291 - t275;
t376 = pkin(8) + t487;
t460 = -pkin(8) * t338 + qJD(5) * t376 + t357 + t524;
t458 = g(1) * t434 - g(2) * t438;
t457 = -t294 * t376 - t532;
t453 = t289 - (-t338 * t431 + t484) * t528;
t450 = -0.2e1 * pkin(1) * t479 - pkin(6) * qJDD(2);
t449 = t312 * t435 - t350 * t484;
t448 = -pkin(8) * t294 + t298 * t528 - t532;
t439 = qJD(2) ^ 2;
t445 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t439 + t458;
t440 = qJD(1) ^ 2;
t444 = pkin(1) * t440 - pkin(6) * qJDD(1) + t459;
t375 = -pkin(4) - t452;
t374 = t417 * t492 + t496;
t373 = -t417 * t495 + t494;
t372 = -t417 * t494 + t495;
t371 = t417 * t496 + t492;
t313 = qJD(4) * t350 + t436 * t381 + t384 * t432;
t285 = pkin(4) * t313 - pkin(8) * t312 + t358;
t279 = qJD(4) * t307 - t314 * t436 + t315 * t432;
t274 = t435 * t275;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t426 + 0.2e1 * t433 * t472) * MDP(4) + 0.2e1 * (t433 * t478 - t479 * t486) * MDP(5) + (qJDD(2) * t433 + t437 * t439) * MDP(6) + (qJDD(2) * t437 - t433 * t439) * MDP(7) + (t433 * t450 + t437 * t445) * MDP(9) + (-t433 * t445 + t437 * t450) * MDP(10) + (-t310 * t391 + t311 * t390 - t328 * t382 + t329 * t380 - t345 * t384 - t346 * t381 + t347 * t356 - t348 * t355 - t459) * MDP(11) + (t311 * t356 + t346 * t329 + t310 * t355 + t345 * t328 - t447 * t420 + t399 * t422 - g(1) * (-t420 * t434 + t438 * t508) - g(2) * (t420 * t438 + t434 * t508)) * MDP(12) + (t295 * t350 + t312 * t455) * MDP(13) + (t295 * t454 - t296 * t350 + t312 * t338 - t313 * t455) * MDP(14) + (t312 * t425 + t350 * t424) * MDP(15) + (-t313 * t425 + t424 * t454) * MDP(16) + (-t279 * t425 + t296 * t360 + t313 * t354 - t320 * t454 - t338 * t358 + t417 * t458 + t424 * t456) * MDP(18) + (-t278 * t425 + t295 * t360 - t307 * t424 + t312 * t354 + t320 * t350 + t358 * t455 - t416 * t458) * MDP(19) + (t350 * t282 + t324 * t449) * MDP(20) + ((-t322 * t435 - t502) * t312 + (-t281 - t284 * t435 + (t322 * t431 - t324 * t435) * qJD(5)) * t350) * MDP(21) + (-t283 * t454 + t350 * t289 + t313 * t324 + t449 * t528) * MDP(22) + (-t350 * t288 + t284 * t454 - t313 * t322 + (-t312 * t431 - t350 * t483) * t528) * MDP(23) + (-t294 * t454 + t313 * t528) * MDP(24) + (-g(1) * t372 - g(2) * t374 - t274 * t454 + t276 * t313 + t279 * t322 - t456 * t284 + (t285 * t528 + t505 + (t291 * t454 - t307 * t528 + t506) * qJD(5)) * t435 + t516 * t431) * MDP(25) + (-g(1) * t371 - g(2) * t373 - t277 * t313 + t279 * t324 - t456 * t283 + (-(-qJD(5) * t307 + t285) * t528 - t505 - t462 * t454 - qJD(5) * t506) * t431 + t516 * t435) * MDP(26) + t458 * MDP(2) + t459 * MDP(3); (t433 * t444 - t509) * MDP(9) + (-t345 * t352 - t346 * t353 + (-t509 + t310 * t429 + t311 * t428 + (-qJD(1) * t399 + t459) * t433) * pkin(2)) * MDP(12) + (t453 + t504) * MDP(23) + (t375 * t284 - t530 * t435 + t457 * t431 + t488 * t322 + (t431 * t489 - t435 * t460) * t528 + t522) * MDP(25) + (t375 * t283 + t457 * t435 - t431 * t535 + t488 * t324 + (t431 * t460 + t435 * t489) * t528 + t521) * MDP(26) + (t338 * t357 + t424 * t452 - t425 * t488 + t520) * MDP(18) + (-t502 * t528 + t534) * MDP(21) + ((t346 + t352) * t382 + (t345 - t353) * t380 + (t347 * t428 - t348 * t429) * pkin(2)) * MDP(11) + qJDD(2) * MDP(8) + (-t357 * t455 - t424 * t487 + t425 * t489 + t527) * MDP(19) + (g(3) * t433 + t437 * t444) * MDP(10) + t433 * qJDD(1) * MDP(6) + (-t433 * t437 * MDP(4) + t486 * MDP(5)) * t440 + MDP(7) * t478 + t533; (-t380 ^ 2 - t382 ^ 2) * MDP(11) + (t345 * t382 - t346 * t380 + t447 - t458) * MDP(12) + (t296 + t501) * MDP(18) + (t295 + t499) * MDP(19) + (t453 - t504) * MDP(25) + (-t493 * t528 - t288 - t503) * MDP(26); (t299 * t425 + t520) * MDP(18) + (t298 * t425 + t527) * MDP(19) + (-t324 * t465 + t534) * MDP(21) + (-t465 * t528 + t289 + t504) * MDP(23) + (-pkin(4) * t284 - t299 * t322 + t448 * t431 + (-t530 - t523) * t435 + t522) * MDP(25) + (-pkin(4) * t283 - t299 * t324 + t448 * t435 + (-t535 + t523) * t431 + t521) * MDP(26) + t533; t324 * t322 * MDP(20) + (-t322 ^ 2 + t324 ^ 2) * MDP(21) + (t322 * t528 + t474) * MDP(22) + (t324 * t528 + t413) * MDP(23) + t294 * MDP(24) + (-g(1) * t373 + g(2) * t371 + t277 * t528 - t290 * t324 + t274) * MDP(25) + (g(1) * t374 - g(2) * t372 + t276 * t528 + t290 * t322) * MDP(26) + ((-t272 + t410) * MDP(26) + (-MDP(23) * t455 - MDP(25) * t291 - MDP(26) * t301) * qJD(5)) * t435 + (-qJD(5) * t455 * MDP(22) + (-qJD(5) * t425 - t295) * MDP(23) + (-t463 + t410) * MDP(25) + t462 * MDP(26)) * t431;];
tau = t1;
