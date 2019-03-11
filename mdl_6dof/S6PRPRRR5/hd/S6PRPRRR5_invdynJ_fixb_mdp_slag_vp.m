% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:46
% EndTime: 2019-03-08 20:43:52
% DurationCPUTime: 4.78s
% Computational Cost: add. (2993->433), mult. (6102->581), div. (0->0), fcn. (4733->14), ass. (0->185)
t404 = qJD(4) + qJD(5);
t416 = cos(qJ(6));
t413 = sin(qJ(5));
t417 = cos(qJ(5));
t418 = cos(qJ(4));
t489 = qJD(2) * t418;
t414 = sin(qJ(4));
t490 = qJD(2) * t414;
t356 = t413 * t490 - t417 * t489;
t412 = sin(qJ(6));
t508 = t356 * t412;
t333 = -t416 * t404 - t508;
t363 = t413 * t418 + t414 * t417;
t357 = t363 * qJD(2);
t528 = qJD(6) + t357;
t533 = t333 * t528;
t441 = t356 * t416 - t404 * t412;
t532 = t441 * t528;
t403 = qJDD(4) + qJDD(5);
t420 = -pkin(2) - pkin(8);
t409 = sin(pkin(6));
t415 = sin(qJ(2));
t503 = t409 * t415;
t468 = qJD(2) * t503;
t370 = qJD(1) * t468;
t419 = cos(qJ(2));
t480 = qJDD(1) * t409;
t463 = t419 * t480;
t438 = qJDD(3) + t370 - t463;
t336 = t420 * qJDD(2) + t438;
t411 = cos(pkin(6));
t474 = t411 * qJDD(1);
t446 = t418 * t336 - t414 * t474;
t494 = qJD(1) * t409;
t471 = t419 * t494;
t447 = qJD(3) - t471;
t355 = t420 * qJD(2) + t447;
t462 = pkin(9) * qJD(2) - t355;
t493 = qJD(1) * t411;
t469 = t418 * t493;
t477 = qJDD(2) * t418;
t288 = -pkin(9) * t477 + qJDD(4) * pkin(4) + (t414 * t462 - t469) * qJD(4) + t446;
t470 = t414 * t493;
t291 = t418 * t474 + (-pkin(9) * qJDD(2) + t336) * t414 + (-t418 * t462 - t470) * qJD(4);
t322 = -pkin(9) * t489 + t418 * t355 - t470;
t314 = qJD(4) * pkin(4) + t322;
t323 = -pkin(9) * t490 + t355 * t414 + t469;
t509 = t323 * t417;
t293 = t314 * t413 + t509;
t522 = qJD(5) * t293 - t417 * t288 + t291 * t413;
t273 = -pkin(5) * t403 + t522;
t408 = sin(pkin(11));
t410 = cos(pkin(11));
t500 = t411 * t419;
t349 = t408 * t415 - t410 * t500;
t351 = t408 * t500 + t410 * t415;
t407 = qJ(4) + qJ(5);
t400 = sin(t407);
t401 = cos(t407);
t502 = t409 * t419;
t505 = t409 * t410;
t506 = t408 * t409;
t436 = -g(3) * (-t400 * t411 - t401 * t502) - g(2) * (t349 * t401 + t400 * t505) - g(1) * (t351 * t401 - t400 * t506);
t432 = -t273 + t436;
t481 = qJD(2) * qJD(4);
t464 = t418 * t481;
t478 = qJDD(2) * t414;
t531 = t464 + t478;
t530 = -t414 * t481 + t477;
t486 = qJD(5) * t413;
t467 = t414 * t486;
t440 = -qJD(2) * t467 + t530 * t413;
t453 = t404 * t418;
t305 = (qJD(2) * t453 + t478) * t417 + t440;
t302 = qJDD(6) + t305;
t362 = t413 * t414 - t417 * t418;
t390 = t414 * pkin(4) + qJ(3);
t316 = pkin(5) * t363 + pkin(10) * t362 + t390;
t501 = t411 * t415;
t350 = t408 * t419 + t410 * t501;
t352 = -t408 * t501 + t410 * t419;
t433 = g(1) * t352 + g(2) * t350 + g(3) * t503;
t430 = t433 * t400;
t529 = t316 * t302 - t430;
t527 = t357 * t404;
t487 = qJD(4) * t418;
t448 = -pkin(4) * t487 - t447;
t450 = g(1) * t351 + g(2) * t349;
t526 = -g(3) * t502 + t450;
t524 = (qJD(5) * t314 + t291) * t417 + t288 * t413 - t323 * t486;
t272 = pkin(10) * t403 + t524;
t510 = t323 * t413;
t292 = t314 * t417 - t510;
t289 = -pkin(5) * t404 - t292;
t472 = t415 * t494;
t492 = qJD(2) * qJ(3);
t364 = t472 + t492;
t344 = pkin(4) * t490 + t364;
t303 = pkin(5) * t357 + pkin(10) * t356 + t344;
t485 = qJD(5) * t417;
t488 = qJD(4) * t414;
t324 = -t413 * t487 - t414 * t485 - t417 * t488 - t418 * t486;
t518 = pkin(9) - t420;
t366 = t518 * t414;
t367 = t518 * t418;
t330 = -t366 * t417 - t367 * t413;
t329 = -t366 * t413 + t367 * t417;
t359 = t518 * t488;
t360 = qJD(4) * t367;
t499 = qJD(5) * t329 - t359 * t413 + t360 * t417 + t363 * t472;
t525 = -(qJD(6) * t303 + t272) * t363 + t289 * t324 + (-qJD(6) * t316 + t499) * t528 - t273 * t362 - t330 * t302 + t526;
t523 = qJD(4) * (-t364 + t472 - t492) - qJDD(4) * t420;
t517 = qJDD(2) * pkin(2);
t304 = -t413 * t478 + t417 * t477 - t527;
t483 = qJD(6) * t416;
t473 = t416 * t304 + t412 * t403 + t404 * t483;
t484 = qJD(6) * t412;
t283 = t356 * t484 + t473;
t516 = t283 * t412;
t515 = t289 * t357;
t514 = t289 * t362;
t513 = t302 * t412;
t512 = t302 * t416;
t507 = t362 * t416;
t504 = t409 * t414;
t498 = qJD(5) * t330 - t359 * t417 - t360 * t413 - t362 * t472;
t497 = t324 * t404 - t362 * t403;
t406 = t418 ^ 2;
t496 = t414 ^ 2 - t406;
t421 = qJD(4) ^ 2;
t422 = qJD(2) ^ 2;
t495 = -t421 - t422;
t491 = qJD(2) * t364;
t479 = qJDD(2) * qJ(3);
t476 = qJDD(4) * t414;
t466 = t528 * t484;
t290 = pkin(10) * t404 + t293;
t444 = t290 * t412 - t303 * t416;
t461 = t289 * t484 - t356 * t444;
t458 = t304 * t412 - t416 * t403;
t457 = t416 * t528;
t315 = -pkin(5) * t356 + pkin(10) * t357;
t393 = pkin(4) * t413 + pkin(10);
t454 = pkin(4) * t489 + qJD(6) * t393 + t315;
t294 = t322 * t413 + t509;
t452 = pkin(4) * t486 - t294;
t295 = t322 * t417 - t510;
t451 = -pkin(4) * t485 + t295;
t325 = -t413 * t488 + t417 * t453 - t467;
t449 = pkin(5) * t325 - pkin(10) * t324 - t448;
t445 = -t393 * t302 + t515;
t280 = t290 * t416 + t303 * t412;
t443 = -t325 * t404 - t363 * t403;
t353 = -t411 * t414 - t418 * t502;
t439 = -t411 * t418 + t414 * t502;
t442 = t353 * t417 + t413 * t439;
t309 = t353 * t413 - t417 * t439;
t437 = t324 * t416 + t362 * t484;
t434 = -t280 * t356 + t289 * t483 - t432 * t412;
t375 = t415 * t480;
t431 = -t375 + t433;
t429 = t463 + t526;
t428 = qJDD(3) - t429;
t337 = t479 + t375 + (qJD(3) + t471) * qJD(2);
t313 = t531 * pkin(4) + t337;
t284 = -qJD(6) * t441 + t458;
t427 = ((t283 - t533) * t416 + (-t284 + t532) * t412) * MDP(23) + (-t441 * t457 + t516) * MDP(22) + (-t412 * t528 ^ 2 - t333 * t356 + t512) * MDP(25) + (-t356 * t441 + t457 * t528 + t513) * MDP(24) + (t304 + t527) * MDP(17) + (-t356 * t404 + (-t404 * t489 - t478) * t417 - t440) * MDP(18) + (t356 ^ 2 - t357 ^ 2) * MDP(16) + t403 * MDP(19) + (-MDP(15) * t357 + MDP(26) * t528) * t356;
t319 = t351 * t400 + t401 * t506;
t321 = -t349 * t400 + t401 * t505;
t342 = -t400 * t502 + t401 * t411;
t426 = g(1) * t319 - g(2) * t321 + g(3) * t342 + t344 * t357 - t524;
t425 = t344 * t356 + t436 - t522;
t424 = qJD(2) * t447 - t420 * t421 + t337 - t433 + t479;
t399 = qJDD(4) * t418;
t394 = -pkin(4) * t417 - pkin(5);
t371 = t416 * t503;
t361 = -qJD(2) * pkin(2) + t447;
t339 = t438 - t517;
t327 = qJD(4) * t439 + t418 * t468;
t326 = qJD(4) * t353 + t414 * t468;
t282 = qJD(5) * t309 + t326 * t413 - t327 * t417;
t281 = qJD(5) * t442 + t326 * t417 + t327 * t413;
t278 = pkin(5) * t305 - pkin(10) * t304 + t313;
t277 = t416 * t278;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (qJDD(1) * t411 ^ 2 - g(3)) * MDP(7) + (qJD(4) * t327 + qJDD(4) * t353) * MDP(13) + (-qJD(4) * t326 + qJDD(4) * t439) * MDP(14) + (-t282 * t404 + t403 * t442) * MDP(20) + (-t281 * t404 - t309 * t403) * MDP(21) + ((-t281 * t412 - t309 * t483) * t528 + (-t309 * t412 + t371) * t302 + t282 * t333 - t442 * t284) * MDP(27) + (-(t281 * t416 - t309 * t484) * t528 - t309 * t512 - t282 * t441 - t442 * t283) * MDP(28) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t415 + t419 * t422) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t419 + t415 * t422) + (-t339 * MDP(7) + (t414 * MDP(13) + t418 * MDP(14)) * t422 + (MDP(20) * t357 - MDP(21) * t356 + MDP(7) * t364 + (MDP(27) * t416 - MDP(28) * t412) * t528) * qJD(2)) * t419 + ((qJD(2) * t361 + t337) * MDP(7) + t531 * MDP(13) + t530 * MDP(14) + t305 * MDP(20) + t304 * MDP(21) - MDP(27) * t466 + (-t483 * t528 - t513) * MDP(28)) * t415) * t409; ((-t333 * t416 + t412 * t441) * t324 + (t516 + t284 * t416 + (-t333 * t412 - t416 * t441) * qJD(6)) * t362) * MDP(23) + (-t283 * t507 - t437 * t441) * MDP(22) + (t305 * t390 + t313 * t363 + t325 * t344 - t329 * t403 - t357 * t448 - t404 * t498 - t430) * MDP(20) + t431 * MDP(4) + t429 * MDP(3) + (t523 * t414 + t424 * t418) * MDP(14) + (t424 * t414 - t523 * t418) * MDP(13) + (t428 - 0.2e1 * t517) * MDP(5) + (-t418 * t421 - t476) * MDP(11) + (qJDD(2) * t406 - 0.2e1 * t414 * t464) * MDP(8) + t443 * MDP(18) + 0.2e1 * (-t414 * t477 + t481 * t496) * MDP(9) + t497 * MDP(17) + (t304 * t390 - t313 * t362 + t324 * t344 - t330 * t403 + t448 * t356 - t433 * t401 + t404 * t499) * MDP(21) + (t337 * qJ(3) + t364 * qJD(3) - t339 * pkin(2) - g(1) * (-pkin(2) * t351 + qJ(3) * t352) - g(2) * (-pkin(2) * t349 + qJ(3) * t350) + (-g(3) * (pkin(2) * t419 + qJ(3) * t415) + (-t361 * t415 - t364 * t419) * qJD(1)) * t409) * MDP(7) + (t283 * t363 - t302 * t507 - t325 * t441 + t437 * t528) * MDP(24) + (t362 * t513 - t284 * t363 - t325 * t333 + (-t324 * t412 + t362 * t483) * t528) * MDP(25) + (-t280 * t325 + t329 * t283 - t498 * t441 + (-(-qJD(6) * t290 + t278) * t363 + qJD(6) * t514 + (qJD(6) * t330 - t449) * t528 - t529) * t412 + t525 * t416) * MDP(28) + (t277 * t363 - t444 * t325 + t329 * t284 + t498 * t333 + (t449 * t528 + (-t290 * t363 - t330 * t528 - t514) * qJD(6) + t529) * t416 + t525 * t412) * MDP(27) + (t302 * t363 + t325 * t528) * MDP(26) + qJDD(2) * MDP(2) + (-t304 * t362 - t324 * t356) * MDP(15) + (-t304 * t363 + t305 * t362 - t324 * t357 + t325 * t356) * MDP(16) + (-t414 * t421 + t399) * MDP(10) + (0.2e1 * qJD(2) * qJD(3) - t431 + 0.2e1 * t479) * MDP(6); qJDD(2) * MDP(5) - t422 * MDP(6) + (t370 + t428 - t491 - t517) * MDP(7) + (t414 * t495 + t399) * MDP(13) + (t418 * t495 - t476) * MDP(14) + (-qJD(2) * t357 + t497) * MDP(20) + (qJD(2) * t356 + t443) * MDP(21) + (t284 * t362 - t324 * t333 - t363 * t513) * MDP(27) + (t283 * t362 + t324 * t441 - t363 * t512) * MDP(28) + ((-qJD(2) * t416 - t325 * t412 - t363 * t483) * MDP(27) + (qJD(2) * t412 - t325 * t416 + t363 * t484) * MDP(28)) * t528; (t394 * t284 + t452 * t333 + (t451 * t528 + t445) * t412 + (-t454 * t528 + t432) * t416 + t461) * MDP(27) + (t394 * t283 + t445 * t416 - t452 * t441 + (t412 * t454 + t416 * t451) * t528 + t434) * MDP(28) + (t295 * t404 + (t356 * t489 - t403 * t413 - t404 * t485) * pkin(4) + t426) * MDP(21) + (t294 * t404 + (-t357 * t489 + t403 * t417 - t404 * t486) * pkin(4) + t425) * MDP(20) + (-t364 * t489 - g(1) * (t351 * t418 - t408 * t504) - g(2) * (t349 * t418 + t410 * t504) - g(3) * t353 + t446) * MDP(13) + (-g(3) * t439 + (-t474 + (g(1) * t408 - g(2) * t410) * t409) * t418 + (-t336 + t450 + t491) * t414) * MDP(14) - MDP(11) * t478 + t427 + MDP(10) * t477 + qJDD(4) * MDP(12) + (t418 * t414 * MDP(8) - t496 * MDP(9)) * t422; (-pkin(5) * t283 + (t292 * t416 + t315 * t412) * t528 + t293 * t441 + t416 * t515 + (t466 - t512) * pkin(10) + t434) * MDP(28) + (-pkin(5) * t284 - t293 * t333 + (-pkin(10) * t302 + t292 * t528 + t515) * t412 + ((-pkin(10) * qJD(6) - t315) * t528 + t432) * t416 + t461) * MDP(27) + (t292 * t404 + t426) * MDP(21) + (t293 * t404 + t425) * MDP(20) + t427; -t441 * t333 * MDP(22) + (-t333 ^ 2 + t441 ^ 2) * MDP(23) + (t473 + t533) * MDP(24) + (-t458 - t532) * MDP(25) + t302 * MDP(26) + (-t412 * t272 + t277 + t280 * t528 + t289 * t441 - g(1) * (-t319 * t412 + t352 * t416) - g(2) * (t321 * t412 + t350 * t416) - g(3) * (-t342 * t412 + t371)) * MDP(27) + (-t416 * t272 - t412 * t278 - t444 * t528 + t289 * t333 - g(1) * (-t319 * t416 - t352 * t412) - g(2) * (t321 * t416 - t350 * t412) - g(3) * (-t342 * t416 - t412 * t503)) * MDP(28) + (MDP(24) * t508 + MDP(25) * t441 - MDP(27) * t280 + MDP(28) * t444) * qJD(6);];
tau  = t1;
