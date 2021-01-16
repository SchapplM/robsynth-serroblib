% Calculate vector of inverse dynamics joint torques for
% S5RRPRR2
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:46
% EndTime: 2021-01-15 21:24:02
% DurationCPUTime: 6.26s
% Computational Cost: add. (3905->376), mult. (9432->487), div. (0->0), fcn. (7361->16), ass. (0->172)
t445 = sin(pkin(9));
t446 = cos(pkin(9));
t454 = cos(qJ(2));
t501 = qJD(1) * t454;
t490 = t446 * t501;
t450 = sin(qJ(2));
t502 = qJD(1) * t450;
t387 = t445 * t502 - t490;
t399 = t445 * t454 + t446 * t450;
t390 = t399 * qJD(1);
t449 = sin(qJ(4));
t453 = cos(qJ(4));
t347 = -t453 * t387 - t390 * t449;
t448 = sin(qJ(5));
t452 = cos(qJ(5));
t469 = t387 * t449 - t453 * t390;
t299 = t347 * t448 - t452 * t469;
t515 = t469 * t448;
t303 = t347 * t452 + t515;
t440 = qJDD(2) + qJDD(4);
t435 = qJDD(5) + t440;
t548 = t435 * MDP(26) + (t299 ^ 2 - t303 ^ 2) * MDP(23) - t303 * MDP(22) * t299;
t389 = t399 * qJD(2);
t494 = qJDD(1) * t454;
t419 = t446 * t494;
t495 = qJDD(1) * t450;
t357 = qJD(1) * t389 + t445 * t495 - t419;
t496 = qJD(1) * qJD(2);
t489 = t450 * t496;
t417 = t445 * t489;
t488 = t454 * t496;
t358 = t399 * qJDD(1) + t446 * t488 - t417;
t499 = qJD(4) * t453;
t500 = qJD(4) * t449;
t293 = -t449 * t357 + t453 * t358 - t387 * t499 - t390 * t500;
t460 = qJD(4) * t469 - t453 * t357 - t358 * t449;
t497 = qJD(5) * t452;
t491 = t452 * t293 + t347 * t497 + t448 * t460;
t498 = qJD(5) * t448;
t274 = t469 * t498 + t491;
t486 = t293 * t448 - t452 * t460;
t461 = -qJD(5) * t299 - t486;
t441 = qJD(2) + qJD(4);
t516 = t347 * t441;
t517 = t469 * t441;
t438 = qJD(5) + t441;
t543 = t299 * t438;
t546 = t303 * t438;
t547 = t440 * MDP(19) + t347 * MDP(15) * t469 + (-t347 ^ 2 + t469 ^ 2) * MDP(16) + (t293 - t516) * MDP(17) + (t460 - t517) * MDP(18) + (t461 + t543) * MDP(25) + (t274 - t546) * MDP(24) + t548;
t447 = -qJ(3) - pkin(6);
t415 = t447 * t454;
t404 = qJD(1) * t415;
t393 = t445 * t404;
t414 = t447 * t450;
t403 = qJD(1) * t414;
t520 = qJD(2) * pkin(2);
t397 = t403 + t520;
t355 = t446 * t397 + t393;
t524 = pkin(7) * t390;
t321 = qJD(2) * pkin(3) + t355 - t524;
t513 = t446 * t404;
t356 = t445 * t397 - t513;
t525 = pkin(7) * t387;
t327 = t356 - t525;
t471 = -t321 * t449 - t327 * t453;
t541 = pkin(8) * t347;
t284 = -t471 + t541;
t432 = pkin(2) * t454 + pkin(1);
t409 = -qJD(1) * t432 + qJD(3);
t364 = pkin(3) * t387 + t409;
t312 = -pkin(4) * t347 + t364;
t442 = qJ(2) + pkin(9);
t439 = qJ(4) + t442;
t430 = qJ(5) + t439;
t423 = sin(t430);
t424 = cos(t430);
t451 = sin(qJ(1));
t455 = cos(qJ(1));
t476 = g(1) * t455 + g(2) * t451;
t536 = g(3) * t423 + t284 * t498 - t312 * t303 + t476 * t424;
t487 = qJD(2) * t447;
t385 = -qJD(3) * t450 + t454 * t487;
t354 = qJDD(2) * pkin(2) + qJD(1) * t385 + qJDD(1) * t414;
t384 = qJD(3) * t454 + t450 * t487;
t361 = qJD(1) * t384 - qJDD(1) * t415;
t306 = t446 * t354 - t361 * t445;
t295 = qJDD(2) * pkin(3) - pkin(7) * t358 + t306;
t307 = t445 * t354 + t446 * t361;
t296 = -pkin(7) * t357 + t307;
t463 = qJD(4) * t471 + t453 * t295 - t449 * t296;
t272 = pkin(4) * t440 - pkin(8) * t293 + t463;
t528 = (qJD(4) * t321 + t296) * t453 + t449 * t295 - t327 * t500;
t273 = pkin(8) * t460 + t528;
t535 = -g(3) * t424 + t452 * t272 - t448 * t273 - t312 * t299 + t476 * t423;
t540 = pkin(8) * t469;
t427 = sin(t439);
t428 = cos(t439);
t534 = g(3) * t427 - t364 * t347 + t476 * t428 - t528;
t533 = -g(3) * t428 + t364 * t469 + t476 * t427 + t463;
t365 = t446 * t414 + t415 * t445;
t340 = -pkin(7) * t399 + t365;
t366 = t445 * t414 - t446 * t415;
t514 = t445 * t450;
t398 = -t446 * t454 + t514;
t341 = -pkin(7) * t398 + t366;
t507 = t449 * t340 + t453 * t341;
t362 = -t403 * t445 + t513;
t328 = t362 + t525;
t363 = t446 * t403 + t393;
t329 = t363 - t524;
t429 = pkin(2) * t446 + pkin(3);
t527 = pkin(2) * t445;
t477 = t453 * t429 - t449 * t527;
t530 = -t477 * qJD(4) + t449 * t328 + t453 * t329;
t383 = t429 * t449 + t453 * t527;
t529 = -t383 * qJD(4) - t453 * t328 + t329 * t449;
t526 = pkin(2) * t450;
t521 = g(3) * t454;
t484 = t453 * t321 - t327 * t449;
t283 = t484 + t540;
t281 = pkin(4) * t441 + t283;
t512 = t452 * t281;
t511 = t452 * t284;
t510 = t529 + t541;
t509 = t530 + t540;
t339 = t446 * t384 + t445 * t385;
t443 = t450 ^ 2;
t505 = -t454 ^ 2 + t443;
t493 = pkin(2) * t489 + qJDD(3);
t434 = t450 * t520;
t368 = pkin(3) * t389 + t434;
t367 = pkin(2) * t502 + pkin(3) * t390;
t482 = t453 * t340 - t341 * t449;
t338 = -t384 * t445 + t446 * t385;
t479 = -qJD(5) * t281 - t273;
t475 = g(1) * t451 - g(2) * t455;
t474 = -t448 * t281 - t511;
t360 = -t398 * t449 + t399 * t453;
t287 = -pkin(8) * t360 + t482;
t359 = t453 * t398 + t399 * t449;
t288 = -pkin(8) * t359 + t507;
t473 = t287 * t452 - t288 * t448;
t472 = t287 * t448 + t288 * t452;
t310 = t452 * t359 + t360 * t448;
t311 = -t359 * t448 + t360 * t452;
t370 = pkin(3) * t398 - t432;
t468 = -0.2e1 * pkin(1) * t496 - pkin(6) * qJDD(2);
t392 = t398 * qJD(2);
t317 = pkin(7) * t392 + t338;
t318 = -pkin(7) * t389 + t339;
t467 = t449 * t317 + t453 * t318 + t340 * t499 - t341 * t500;
t381 = -qJDD(1) * t432 + t493;
t456 = qJD(2) ^ 2;
t465 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t456 + t475;
t457 = qJD(1) ^ 2;
t464 = pkin(1) * t457 - pkin(6) * qJDD(1) + t476;
t323 = pkin(3) * t357 + t381;
t462 = -t507 * qJD(4) + t453 * t317 - t318 * t449;
t437 = cos(t442);
t436 = sin(t442);
t382 = pkin(4) + t477;
t322 = pkin(4) * t359 + t370;
t313 = -pkin(4) * t469 + t367;
t309 = qJD(4) * t360 + t453 * t389 - t392 * t449;
t308 = -qJD(4) * t359 - t389 * t449 - t392 * t453;
t297 = pkin(4) * t309 + t368;
t280 = -pkin(4) * t460 + t323;
t279 = qJD(5) * t311 + t308 * t448 + t452 * t309;
t278 = -qJD(5) * t310 + t308 * t452 - t309 * t448;
t277 = -pkin(8) * t308 + t462;
t276 = -pkin(8) * t309 + t467;
t1 = [(t297 * t299 + t322 * t274 + t280 * t311 + t312 * t278 - (qJD(5) * t473 + t276 * t452 + t277 * t448) * t438 - t472 * t435 - t475 * t423) * MDP(28) + (t274 * t311 + t278 * t299) * MDP(22) + (t364 * t309 + t323 * t359 - t347 * t368 - t370 * t460 + t428 * t475 + t440 * t482 + t441 * t462) * MDP(20) + (-t293 * t359 + t308 * t347 + t309 * t469 + t360 * t460) * MDP(16) + (t293 * t360 - t308 * t469) * MDP(15) + (t370 * t293 + t364 * t308 + t323 * t360 - t368 * t469 - t427 * t475 - t440 * t507 - t441 * t467) * MDP(21) + (-t306 * t399 - t307 * t398 - t338 * t390 - t339 * t387 + t355 * t392 - t356 * t389 - t357 * t366 - t358 * t365 - t476) * MDP(13) + (-qJDD(2) * t366 - t358 * t432 + t381 * t399 - t392 * t409 - t475 * t436 + (t390 * t526 - t339) * qJD(2)) * MDP(12) + (qJDD(2) * t450 + t454 * t456) * MDP(6) + (qJDD(2) * t454 - t450 * t456) * MDP(7) + qJDD(1) * MDP(1) + (-t309 * t441 - t359 * t440) * MDP(18) + (t308 * t441 + t360 * t440) * MDP(17) + (-t279 * t438 - t310 * t435) * MDP(25) + (t278 * t438 + t311 * t435) * MDP(24) + (t450 * t468 + t454 * t465) * MDP(9) + (-t450 * t465 + t454 * t468) * MDP(10) + t475 * MDP(2) + t476 * MDP(3) + (qJDD(1) * t443 + 0.2e1 * t450 * t488) * MDP(4) + (t307 * t366 + t356 * t339 + t306 * t365 + t355 * t338 - t381 * t432 + t409 * t434 - g(1) * (-t432 * t451 - t447 * t455) - g(2) * (t432 * t455 - t447 * t451)) * MDP(14) + 0.2e1 * (t450 * t494 - t496 * t505) * MDP(5) + (qJDD(2) * t365 - t357 * t432 + t381 * t398 + t389 * t409 + t475 * t437 + (t387 * t526 + t338) * qJD(2)) * MDP(11) + (-t274 * t310 + t278 * t303 - t279 * t299 + t311 * t461) * MDP(23) + (-t297 * t303 - t322 * t461 + t280 * t310 + t312 * t279 + (-qJD(5) * t472 - t276 * t448 + t277 * t452) * t438 + t473 * t435 + t475 * t424) * MDP(27); ((t356 + t362) * t390 + (-t355 + t363) * t387 + (-t357 * t445 - t358 * t446) * pkin(2)) * MDP(13) + (-t313 * t299 + (-t382 * t435 - t272 + (qJD(5) * t383 - t510) * t438) * t448 + (-t383 * t435 + (-qJD(5) * t382 + t509) * t438 + t479) * t452 + t536) * MDP(28) + (-g(3) * t437 - qJD(2) * t362 - t390 * t409 + t476 * t436 + (qJDD(2) * t446 - t387 * t502) * pkin(2) + t306) * MDP(11) + (g(3) * t436 + qJD(2) * t363 + t387 * t409 + t476 * t437 + (-qJDD(2) * t445 - t390 * t502) * pkin(2) - t307) * MDP(12) + (g(3) * t450 + t454 * t464) * MDP(10) + MDP(6) * t495 + MDP(7) * t494 + (t347 * t367 + t477 * t440 + t529 * t441 + t533) * MDP(20) + (t367 * t469 - t383 * t440 + t530 * t441 + t534) * MDP(21) + qJDD(2) * MDP(8) + ((t382 * t452 - t383 * t448) * t435 + t313 * t303 + (t509 * t448 + t510 * t452) * t438 + ((-t382 * t448 - t383 * t452) * t438 + t474) * qJD(5) + t535) * MDP(27) + (t450 * t464 - t521) * MDP(9) + (-t355 * t362 - t356 * t363 + (-t521 + t306 * t446 + t307 * t445 + (-qJD(1) * t409 + t476) * t450) * pkin(2)) * MDP(14) + (-MDP(4) * t450 * t454 + MDP(5) * t505) * t457 + t547; -t419 * MDP(11) - t417 * MDP(12) + (-t387 ^ 2 - t390 ^ 2) * MDP(13) + (t355 * t390 + t356 * t387 - t475 + t493) * MDP(14) + (-t460 - t517) * MDP(20) + (t293 + t516) * MDP(21) + (-t461 + t543) * MDP(27) + (t274 + t546) * MDP(28) + (MDP(11) * t514 + t399 * MDP(12) - t432 * MDP(14)) * qJDD(1) + ((t445 * t501 + t446 * t502 + t390) * MDP(11) + (-t387 + t490) * MDP(12)) * qJD(2); (-t441 * t471 + t533) * MDP(20) + (t441 * t484 + t534) * MDP(21) + (-(-t283 * t448 - t511) * t438 + t474 * qJD(5) + (-t303 * t469 + t452 * t435 - t438 * t498) * pkin(4) + t535) * MDP(27) + ((-t284 * t438 - t272) * t448 + (t283 * t438 + t479) * t452 + (t299 * t469 - t448 * t435 - t438 * t497) * pkin(4) + t536) * MDP(28) + t547; (t491 - t546) * MDP(24) + (-t486 + t543) * MDP(25) + (-t438 * t474 + t535) * MDP(27) + (-t452 * t273 - t448 * t272 + (-t284 * t448 + t512) * t438 + t536) * MDP(28) + (MDP(24) * t515 - t299 * MDP(25) + t474 * MDP(27) - MDP(28) * t512) * qJD(5) + t548;];
tau = t1;
