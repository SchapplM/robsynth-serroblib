% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:23
% EndTime: 2019-03-09 06:55:34
% DurationCPUTime: 4.64s
% Computational Cost: add. (5848->360), mult. (13863->486), div. (0->0), fcn. (10330->10), ass. (0->177)
t424 = cos(qJ(6));
t486 = qJD(6) * t424;
t422 = sin(qJ(4));
t423 = sin(qJ(3));
t493 = qJD(1) * t423;
t478 = t422 * t493;
t426 = cos(qJ(4));
t427 = cos(qJ(3));
t492 = qJD(1) * t427;
t479 = t426 * t492;
t387 = -t478 + t479;
t388 = -t422 * t492 - t426 * t493;
t421 = sin(qJ(5));
t425 = cos(qJ(5));
t349 = t425 * t387 + t388 * t421;
t538 = t349 * t424;
t541 = t486 - t538;
t415 = qJD(3) + qJD(4);
t393 = t422 * t427 + t423 * t426;
t529 = qJD(1) * t393;
t431 = t415 * t529;
t405 = sin(pkin(11)) * pkin(1) + pkin(7);
t522 = pkin(8) + t405;
t473 = t522 * qJD(1);
t372 = qJD(2) * t423 + t427 * t473;
t365 = t372 * qJD(3);
t490 = qJD(4) * t422;
t466 = -t422 * t365 - t372 * t490;
t371 = t427 * qJD(2) - t473 * t423;
t364 = t371 * qJD(3);
t521 = qJD(3) * pkin(3);
t369 = t371 + t521;
t530 = t426 * (qJD(4) * t369 + t364);
t300 = -pkin(9) * t431 + t466 + t530;
t385 = t388 * pkin(9);
t366 = t422 * t372;
t467 = t426 * t369 - t366;
t326 = t385 + t467;
t323 = pkin(4) * t415 + t326;
t483 = qJD(1) * qJD(3);
t477 = t427 * t483;
t354 = qJD(4) * t479 - t415 * t478 + t426 * t477;
t368 = t426 * t372;
t447 = -t369 * t422 - t368;
t468 = -t422 * t364 - t426 * t365;
t433 = qJD(4) * t447 + t468;
t301 = -pkin(9) * t354 + t433;
t523 = pkin(9) * t387;
t327 = -t447 + t523;
t489 = qJD(5) * t421;
t470 = t301 * t421 - t327 * t489;
t273 = (qJD(5) * t323 + t300) * t425 + t470;
t406 = -cos(pkin(11)) * pkin(1) - pkin(2);
t394 = -pkin(3) * t427 + t406;
t389 = t394 * qJD(1);
t358 = -pkin(4) * t387 + t389;
t512 = t349 * t358;
t540 = -t273 - t512;
t446 = t387 * t421 - t425 * t388;
t420 = sin(qJ(6));
t487 = qJD(6) * t420;
t488 = qJD(5) * t425;
t311 = t425 * t354 + t387 * t488 + t388 * t489 - t421 * t431;
t413 = qJD(5) + t415;
t498 = t424 * t311 + t413 * t486;
t289 = -t446 * t487 + t498;
t287 = t289 * t420;
t288 = t289 * t424;
t336 = t413 * t420 + t424 * t446;
t518 = t311 * t420;
t290 = qJD(6) * t336 + t518;
t312 = qJD(5) * t446 + t354 * t421 + t425 * t431;
t310 = t424 * t312;
t510 = t446 * t420;
t334 = -t424 * t413 + t510;
t484 = -qJD(6) + t349;
t308 = t420 * t312;
t499 = -t484 * t486 + t308;
t537 = t484 * t420;
t539 = -t312 * MDP(22) - t349 ^ 2 * MDP(20) + (-t349 * t413 + t311) * MDP(21) + (-t349 * MDP(19) + MDP(20) * t446 + t413 * MDP(22) + MDP(30) * t484) * t446 + (t336 * t541 + t287) * MDP(26) + (-t336 * t446 + t484 * t538 + t499) * MDP(28) + (t334 * t446 - t484 * t537 + t310) * MDP(29) + (-t420 * t290 - t334 * t541 + t336 * t537 + t288) * MDP(27);
t516 = t327 * t421;
t298 = t323 * t425 - t516;
t293 = -pkin(5) * t413 - t298;
t520 = t293 * t349;
t515 = t327 * t425;
t299 = t323 * t421 + t515;
t471 = t300 * t421 - t425 * t301;
t274 = qJD(5) * t299 + t471;
t514 = t446 * t358;
t536 = -t274 - t514;
t324 = pkin(5) * t446 - pkin(10) * t349;
t532 = MDP(5) * t423;
t531 = MDP(6) * (t423 ^ 2 - t427 ^ 2);
t294 = pkin(10) * t413 + t299;
t313 = -pkin(5) * t349 - pkin(10) * t446 + t358;
t449 = t294 * t420 - t313 * t424;
t443 = -t274 * t424 + t293 * t487 + t446 * t449;
t282 = t294 * t424 + t313 * t420;
t454 = t274 * t420 + t282 * t446 + t293 * t486;
t438 = t393 * qJD(4);
t362 = qJD(3) * t393 + t438;
t474 = qJD(3) * t522;
t383 = t423 * t474;
t384 = t427 * t474;
t391 = t522 * t427;
t390 = t522 * t423;
t508 = t390 * t426;
t439 = -qJD(4) * t508 - t426 * t383 - t422 * t384 - t391 * t490;
t314 = -pkin(9) * t362 + t439;
t392 = t422 * t423 - t426 * t427;
t361 = t415 * t392;
t445 = t390 * t422 - t391 * t426;
t432 = qJD(4) * t445 + t383 * t422 - t426 * t384;
t315 = pkin(9) * t361 + t432;
t338 = -pkin(9) * t393 - t391 * t422 - t508;
t339 = -pkin(9) * t392 - t445;
t448 = t338 * t425 - t339 * t421;
t275 = qJD(5) * t448 + t314 * t425 + t315 * t421;
t317 = t338 * t421 + t339 * t425;
t356 = t425 * t392 + t393 * t421;
t319 = -qJD(5) * t356 - t361 * t425 - t362 * t421;
t357 = -t392 * t421 + t393 * t425;
t370 = pkin(4) * t392 + t394;
t322 = pkin(5) * t356 - pkin(10) * t357 + t370;
t526 = t274 * t357 + t293 * t319 - t317 * t312 + (qJD(6) * t322 + t275) * t484 - (qJD(6) * t313 + t273) * t356;
t524 = pkin(4) * t388;
t519 = t293 * t357;
t517 = t322 * t312;
t509 = t389 * t388;
t507 = t421 * t422;
t506 = t422 * t425;
t428 = qJD(3) ^ 2;
t505 = t423 * t428;
t504 = t427 * t428;
t320 = qJD(5) * t357 - t361 * t421 + t425 * t362;
t502 = t289 * t356 + t336 * t320;
t465 = -t371 * t422 - t368;
t328 = t465 - t523;
t496 = t426 * t371 - t366;
t329 = t385 + t496;
t410 = pkin(3) * t426 + pkin(4);
t500 = t328 * t421 + t329 * t425 - t410 * t488 - (-t422 * t489 + (t425 * t426 - t507) * qJD(4)) * pkin(3);
t497 = t328 * t425 - t329 * t421 + t410 * t489 + (t422 * t488 + (t421 * t426 + t506) * qJD(4)) * pkin(3);
t396 = qJD(1) * t406;
t412 = t423 * t521;
t411 = pkin(3) * t493;
t482 = t357 * t308;
t481 = t357 * t310;
t476 = -pkin(3) * t415 - t369;
t475 = -pkin(4) * t413 - t323;
t351 = pkin(4) * t362 + t412;
t321 = t324 - t524;
t382 = pkin(3) * t506 + t410 * t421 + pkin(10);
t458 = qJD(6) * t382 + t321 + t411;
t408 = pkin(4) * t421 + pkin(10);
t457 = qJD(6) * t408 + t321;
t302 = t326 * t421 + t515;
t456 = pkin(4) * t489 - t302;
t303 = t326 * t425 - t516;
t455 = -pkin(4) * t488 + t303;
t452 = -t290 * t356 - t320 * t334;
t451 = -t312 * t382 - t520;
t450 = -t312 * t408 - t520;
t444 = 0.2e1 * qJD(3) * t396;
t442 = -t389 * t387 - t466;
t441 = -t319 * t420 - t357 * t486;
t440 = -t319 * t424 + t357 * t487;
t340 = pkin(4) * t431 + qJD(3) * t411;
t430 = t388 * t387 * MDP(12) + (-t387 * t415 + t354) * MDP(14) + (-t388 * t415 - t431) * MDP(15) + (-t387 ^ 2 + t388 ^ 2) * MDP(13) + t539;
t409 = -pkin(4) * t425 - pkin(5);
t381 = pkin(3) * t507 - t410 * t425 - pkin(5);
t373 = t411 - t524;
t284 = pkin(5) * t320 - pkin(10) * t319 + t351;
t283 = t312 * pkin(5) - t311 * pkin(10) + t340;
t280 = t424 * t283;
t276 = qJD(5) * t317 + t314 * t421 - t315 * t425;
t1 = [0.2e1 * t477 * t532 + MDP(7) * t504 + (-t354 * t392 - t361 * t387 + t388 * t362 - t393 * t431) * MDP(13) + (t312 * t370 + t320 * t358 + t340 * t356 - t349 * t351) * MDP(24) + (t311 * t370 + t319 * t358 + t340 * t357 + t351 * t446) * MDP(25) + (t354 * t393 + t361 * t388) * MDP(12) + (t312 * t356 - t320 * t484) * MDP(30) + (-t311 * t356 - t312 * t357 + t319 * t349 - t320 * t446) * MDP(20) + (t311 * t357 + t319 * t446) * MDP(19) + (t288 * t357 - t336 * t440) * MDP(26) + ((-t334 * t424 - t336 * t420) * t319 + (-t287 - t290 * t424 + (t334 * t420 - t336 * t424) * qJD(6)) * t357) * MDP(27) + (t276 * t336 - t282 * t320 - t448 * t289 + ((-qJD(6) * t317 + t284) * t484 - t517 - (-qJD(6) * t294 + t283) * t356 - qJD(6) * t519) * t420 + t526 * t424) * MDP(32) + (t276 * t334 + t280 * t356 - t449 * t320 - t448 * t290 + (-t284 * t484 + t517 + (-t294 * t356 + t317 * t484 + t519) * qJD(6)) * t424 + t526 * t420) * MDP(31) - MDP(8) * t505 + (t405 * t505 + t427 * t444) * MDP(11) + (t440 * t484 + t481 + t502) * MDP(28) + (-t405 * t504 + t423 * t444) * MDP(10) + (t394 * t354 - t389 * t361 + (-t388 + t529) * t412) * MDP(18) + (-t387 * t412 + t389 * t362 + (t394 * t438 + (t423 * pkin(3) * t392 + t393 * t394) * qJD(3)) * qJD(1)) * MDP(17) + (-t441 * t484 + t452 - t482) * MDP(29) - 0.2e1 * t483 * t531 + (-t361 * MDP(14) - t362 * MDP(15) + MDP(17) * t432 - MDP(18) * t439) * t415 + (MDP(21) * t319 - MDP(22) * t320 - MDP(24) * t276 - MDP(25) * t275) * t413; (-t452 - t482) * MDP(31) + (-t481 + t502) * MDP(32) + (-MDP(10) * t423 - MDP(11) * t427) * t428 + (-t362 * MDP(17) + t361 * MDP(18)) * t415 + (-MDP(24) * t320 - MDP(25) * t319) * t413 - (MDP(31) * t441 + MDP(32) * t440) * t484; (t381 * t289 + t451 * t424 + t497 * t336 - (t420 * t458 + t424 * t500) * t484 + t454) * MDP(32) + (t349 * t373 - t413 * t497 + t536) * MDP(24) + t430 + (-t373 * t446 + t413 * t500 + t540) * MDP(25) + (t381 * t290 + t451 * t420 + t497 * t334 - (t420 * t500 - t424 * t458) * t484 + t443) * MDP(31) + (t388 * t411 + t496 * t415 + (qJD(4) * t476 - t364) * t426 + t442) * MDP(18) + (t387 * t411 + t509 - t465 * t415 + (t422 * t476 - t368) * qJD(4) + t468) * MDP(17) + (-t427 * t532 + t531) * qJD(1) ^ 2 + (-MDP(10) * t493 - MDP(11) * t492) * t396; (t409 * t289 + t450 * t424 + t456 * t336 - (t420 * t457 + t424 * t455) * t484 + t454) * MDP(32) + (-t415 * t447 + t433 + t509) * MDP(17) + (t415 * t467 + t442 - t530) * MDP(18) + (t446 * t524 + t303 * t413 - t512 + (qJD(5) * t475 - t300) * t425 - t470) * MDP(25) + (-t349 * t524 + t302 * t413 - t514 + (t421 * t475 - t515) * qJD(5) - t471) * MDP(24) + (t409 * t290 + t450 * t420 + t456 * t334 - (t420 * t455 - t424 * t457) * t484 + t443) * MDP(31) + t430; (t299 * t413 + t536) * MDP(24) + (t298 * t413 + t540) * MDP(25) + (-pkin(5) * t290 + (-t298 * t420 + t324 * t424) * t484 - t299 * t334 - t420 * t520 - t499 * pkin(10) + t443) * MDP(31) + (-pkin(5) * t289 - (t298 * t424 + t324 * t420) * t484 - t299 * t336 - t293 * t538 + (-t484 * t487 - t310) * pkin(10) + t454) * MDP(32) + t539; t336 * t334 * MDP(26) + (-t334 ^ 2 + t336 ^ 2) * MDP(27) + (-t334 * t484 + t498) * MDP(28) + (-t336 * t484 - t518) * MDP(29) + t312 * MDP(30) + (-t273 * t420 - t282 * t484 - t293 * t336 + t280) * MDP(31) + (-t273 * t424 - t283 * t420 + t293 * t334 + t449 * t484) * MDP(32) + (-MDP(28) * t510 - MDP(29) * t336 - MDP(31) * t282 + MDP(32) * t449) * qJD(6);];
tauc  = t1;
