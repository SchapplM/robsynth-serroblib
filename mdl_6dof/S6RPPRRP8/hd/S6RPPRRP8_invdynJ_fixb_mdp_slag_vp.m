% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP8
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:35
% EndTime: 2019-03-09 02:16:43
% DurationCPUTime: 5.90s
% Computational Cost: add. (4685->471), mult. (9207->560), div. (0->0), fcn. (6467->10), ass. (0->186)
t393 = sin(pkin(9));
t394 = cos(pkin(9));
t398 = sin(qJ(4));
t513 = cos(qJ(4));
t423 = -t393 * t513 - t398 * t394;
t471 = qJD(1) * t393;
t455 = t398 * t471;
t409 = qJD(4) * t455 + qJDD(1) * t423;
t456 = t513 * t394;
t446 = qJD(1) * t456;
t318 = qJD(4) * t446 - t409;
t315 = qJDD(5) + t318;
t396 = -pkin(1) - qJ(3);
t364 = qJD(1) * t396 + qJD(2);
t452 = -pkin(7) * qJD(1) + t364;
t336 = t452 * t393;
t337 = t452 * t394;
t304 = t513 * t336 + t398 * t337;
t301 = qJD(4) * pkin(8) + t304;
t345 = t423 * qJD(1);
t347 = t446 - t455;
t375 = qJD(1) * qJ(2) + qJD(3);
t358 = pkin(3) * t471 + t375;
t302 = -pkin(4) * t345 - pkin(8) * t347 + t358;
t397 = sin(qJ(5));
t400 = cos(qJ(5));
t279 = t301 * t400 + t302 * t397;
t525 = qJD(5) - t345;
t270 = qJ(6) * t525 + t279;
t528 = t270 * t525;
t470 = qJD(4) * t397;
t326 = t347 * t400 + t470;
t449 = t326 * t525;
t422 = -t398 * t393 + t456;
t413 = t422 * qJDD(1);
t520 = qJD(4) * t345;
t403 = t413 + t520;
t527 = -qJD(4) * qJD(5) - t403;
t486 = t394 * MDP(8);
t526 = t393 * MDP(7) + t486;
t389 = pkin(9) + qJ(4);
t376 = sin(t389);
t377 = cos(t389);
t399 = sin(qJ(1));
t401 = cos(qJ(1));
t522 = g(1) * t399 - g(2) * t401;
t412 = g(3) * t376 - t377 * t522;
t381 = t393 * pkin(3);
t369 = qJ(2) + t381;
t316 = -pkin(4) * t423 - pkin(8) * t422 + t369;
t506 = -pkin(7) + t396;
t355 = t506 * t393;
t356 = t506 * t394;
t321 = t355 * t513 + t398 * t356;
t478 = t397 * t316 + t400 * t321;
t524 = -t398 * t336 + t513 * t337;
t390 = qJDD(1) * qJ(2);
t391 = qJD(1) * qJD(2);
t521 = t390 + t391;
t360 = qJDD(3) + t521;
t510 = g(1) * t401;
t442 = g(2) * t399 + t510;
t523 = t360 - t442;
t472 = t393 ^ 2 + t394 ^ 2;
t519 = MDP(23) + MDP(25);
t460 = MDP(24) - MDP(27);
t516 = -qJD(1) * qJD(3) + qJDD(1) * t396;
t354 = qJDD(2) + t516;
t450 = -pkin(7) * qJDD(1) + t354;
t330 = t450 * t393;
t331 = t450 * t394;
t426 = t330 * t513 + t398 * t331;
t276 = qJDD(4) * pkin(8) + qJD(4) * t524 + t426;
t462 = qJDD(1) * t393;
t351 = pkin(3) * t462 + t360;
t285 = t318 * pkin(4) - pkin(8) * t403 + t351;
t467 = qJD(5) * t400;
t468 = qJD(5) * t397;
t447 = t397 * t276 - t400 * t285 + t301 * t467 + t302 * t468;
t512 = pkin(5) * t315;
t264 = qJDD(6) + t447 - t512;
t518 = -t264 + t528;
t300 = -qJD(4) * pkin(4) - t524;
t324 = -t400 * qJD(4) + t347 * t397;
t280 = t324 * pkin(5) - t326 * qJ(6) + t300;
t511 = pkin(8) * t315;
t517 = t280 * t525 - t511;
t508 = g(3) * t377;
t411 = -t376 * t522 - t508;
t515 = t326 ^ 2;
t514 = 0.2e1 * t391;
t507 = g(3) * t400;
t505 = pkin(8) * qJD(5);
t504 = pkin(1) * qJDD(1);
t503 = qJ(6) * t315;
t454 = qJD(4) * t513;
t469 = qJD(4) * t398;
t418 = t398 * t330 - t331 * t513 + t336 * t454 + t337 * t469;
t277 = -qJDD(4) * pkin(4) + t418;
t289 = -t397 * qJDD(4) + t347 * t468 + t400 * t527;
t440 = -t400 * qJDD(4) + t347 * t467;
t290 = t397 * t413 + (qJD(5) + t345) * t470 + t440;
t265 = t290 * pkin(5) + t289 * qJ(6) - t326 * qJD(6) + t277;
t502 = t265 * t422;
t501 = t279 * t525;
t500 = t289 * t397;
t499 = t290 * t400;
t498 = t324 * t345;
t497 = t324 * t347;
t496 = t324 * t397;
t495 = t324 * t400;
t494 = t326 * t324;
t493 = t326 * t347;
t492 = t326 * t397;
t491 = t326 * t400;
t490 = t525 * t400;
t489 = t422 * t400;
t488 = t377 * t399;
t308 = t397 * t315;
t485 = t397 * t399;
t484 = t397 * t401;
t483 = t399 * t400;
t310 = t400 * t315;
t482 = t401 * t400;
t480 = -t397 * t290 - t324 * t467;
t317 = pkin(4) * t347 - pkin(8) * t345;
t479 = t397 * t317 + t400 * t524;
t348 = -t393 * t454 - t394 * t469;
t477 = t348 * qJD(4) + qJDD(4) * t422;
t438 = pkin(5) * t397 - qJ(6) * t400;
t476 = -qJD(6) * t397 + t438 * t525 - t304;
t475 = g(2) * t377 * t482 + t376 * t507;
t474 = t401 * pkin(1) + t399 * qJ(2);
t278 = -t301 * t397 + t302 * t400;
t465 = qJD(6) - t278;
t459 = g(1) * t488;
t458 = t525 * t505;
t453 = qJDD(2) - t522;
t451 = t472 * MDP(9);
t448 = t397 * t525;
t445 = -t265 - t458;
t340 = t376 * t485 - t482;
t342 = t376 * t484 + t483;
t444 = g(1) * t342 + g(2) * t340;
t341 = t376 * t483 + t484;
t343 = t376 * t482 - t485;
t443 = -g(1) * t343 - g(2) * t341;
t439 = -pkin(5) * t400 - qJ(6) * t397;
t416 = t400 * t276 + t397 * t285 - t301 * t468 + t302 * t467;
t263 = qJD(6) * t525 + t416 + t503;
t437 = t263 * t400 + t264 * t397;
t269 = -pkin(5) * t525 + t465;
t436 = t269 * t400 - t270 * t397;
t435 = t269 * t397 + t270 * t400;
t433 = t269 * t525 + t263;
t429 = pkin(4) - t439;
t428 = -t345 * t490 + t467 * t525 + t308;
t427 = t310 + (t345 * t397 - t468) * t525;
t424 = -t398 * t355 + t356 * t513;
t421 = -t348 * t397 - t422 * t467;
t420 = t348 * t400 - t422 * t468;
t419 = t300 * t525 - t511;
t297 = qJD(3) * t423 + qJD(4) * t424;
t349 = -t393 * t469 + t394 * t454;
t313 = pkin(4) * t349 - pkin(8) * t348 + qJD(2);
t415 = t400 * t297 + t397 * t313 + t316 * t467 - t321 * t468;
t414 = pkin(4) * t376 - pkin(8) * t377 + t381;
t410 = g(1) * t340 - g(2) * t342 + t397 * t508 - t447;
t406 = t280 * t326 + qJDD(6) - t410;
t405 = -g(1) * t341 + g(2) * t343 - t377 * t507 + t416;
t404 = (t491 + t496) * MDP(26) + t436 * MDP(28);
t298 = qJD(3) * t422 + qJD(4) * t321;
t402 = qJD(1) ^ 2;
t395 = -pkin(7) - qJ(3);
t383 = t401 * qJ(2);
t292 = pkin(5) * t326 + qJ(6) * t324;
t291 = t422 * t438 - t424;
t284 = pkin(5) * t423 - t316 * t400 + t321 * t397;
t283 = -qJ(6) * t423 + t478;
t275 = -pkin(5) * t347 - t317 * t400 + t397 * t524;
t274 = qJ(6) * t347 + t479;
t271 = t324 * t525 - t289;
t268 = t438 * t348 - (qJD(5) * t439 + qJD(6) * t400) * t422 + t298;
t267 = -pkin(5) * t349 + qJD(5) * t478 + t297 * t397 - t313 * t400;
t266 = qJ(6) * t349 - qJD(6) * t423 + t415;
t1 = [(t453 - 0.2e1 * t504) * MDP(4) + (t263 * t283 + t270 * t266 + t265 * t291 + t280 * t268 + t264 * t284 + t269 * t267 - g(1) * (pkin(5) * t343 + qJ(6) * t342 + t383) - g(2) * (pkin(5) * t341 + qJ(6) * t340 + t474) + (-g(1) * t414 + g(2) * t395) * t401 + (-g(1) * (-pkin(1) + t395) - g(2) * t414) * t399) * MDP(28) + t477 * MDP(13) + (-t289 * t489 + t326 * t420) * MDP(18) + t442 * MDP(3) + (t360 * qJ(2) + t375 * qJD(2) - g(1) * (t396 * t399 + t383) - g(2) * (qJ(3) * t401 + t474) + (-t364 * qJD(3) + t396 * t354) * t472) * MDP(10) + (-g(2) * t488 + qJD(2) * t347 - t297 * qJD(4) - t321 * qJDD(4) + t358 * t348 + t351 * t422 + t369 * t403 - t377 * t510) * MDP(17) + ((-t492 - t495) * t348 - (-t500 + t499 + (t491 - t496) * qJD(5)) * t422) * MDP(19) + (-t266 * t324 + t267 * t326 - t283 * t290 - t284 * t289 + t442 * t377 + t436 * t348 - (qJD(5) * t435 + t263 * t397 - t264 * t400) * t422) * MDP(26) + (t347 * t348 + t403 * t422) * MDP(11) + (-qJD(4) * t349 + qJDD(4) * t423) * MDP(14) + (-t318 * t422 + t348 * t345 - t347 * t349 + t403 * t423) * MDP(12) + (-qJD(2) * t345 - qJD(4) * t298 + qJDD(4) * t424 + t318 * t369 + t349 * t358 - t351 * t423 - t442 * t376) * MDP(16) + t522 * MDP(2) + (t522 + t472 * (-t354 - t516)) * MDP(9) + t526 * (t521 + t523) + (-(qJDD(2) - t504) * pkin(1) - g(1) * (-pkin(1) * t399 + t383) - g(2) * t474 + (t390 + t514) * qJ(2)) * MDP(6) + (0.2e1 * t390 + t514 - t442) * MDP(5) + (-t263 * t423 - t265 * t489 + t266 * t525 - t268 * t326 + t270 * t349 - t280 * t420 + t283 * t315 + t289 * t291 - t444) * MDP(27) + (t264 * t423 - t267 * t525 + t268 * t324 - t269 * t349 - t280 * t421 - t284 * t315 + t290 * t291 + t397 * t502 + t443) * MDP(25) + (-t315 * t423 + t349 * t525) * MDP(22) + (t290 * t423 - t308 * t422 - t324 * t349 + t421 * t525) * MDP(21) + (t289 * t423 + t310 * t422 + t326 * t349 + t420 * t525) * MDP(20) + (t447 * t423 + t278 * t349 + t298 * t324 - t424 * t290 + ((-qJD(5) * t321 + t313) * t525 + t316 * t315 + t300 * qJD(5) * t422) * t400 + ((-qJD(5) * t316 - t297) * t525 - t321 * t315 + t277 * t422 + t300 * t348) * t397 + t443) * MDP(23) + (t277 * t489 - t279 * t349 + t289 * t424 + t298 * t326 + t300 * t420 - t315 * t478 - t415 * t525 + t416 * t423 + t444) * MDP(24) + qJDD(1) * MDP(1); t453 * MDP(6) + (t354 * t472 - t522) * MDP(10) + t477 * MDP(16) + (-t280 * t348 - t502 - t522) * MDP(28) + (-qJD(4) * MDP(17) + (t492 - t495) * MDP(26) + t435 * MDP(28)) * t349 + (-t375 * MDP(10) + t345 * MDP(16) - t347 * MDP(17) + t404) * qJD(1) + (-MDP(6) * qJ(2) - MDP(5) - t526) * t402 + (-pkin(1) * MDP(6) + MDP(4) - t451) * qJDD(1) + (t519 * (-qJD(1) * t400 - t349 * t397) + t460 * (qJD(1) * t397 - t349 * t400)) * t525 - (-qJDD(4) * MDP(17) + (-t499 - t500) * MDP(26) + t437 * MDP(28) + (-t397 * t519 - t400 * t460) * t315 + ((t397 * t460 - t400 * t519) * t525 + t404) * qJD(5)) * t423 + t519 * (-t290 * t422 - t348 * t324) - t460 * (-t289 * t422 + t326 * t348); MDP(7) * t462 + qJDD(1) * t486 - t402 * t451 + (qJD(1) * t364 * t472 + t523) * MDP(10) + ((t446 + t347) * qJD(4) - t409) * MDP(16) + (t413 + 0.2e1 * t520) * MDP(17) + (t427 - t497) * MDP(23) + (-t490 * t525 - t308 - t493) * MDP(24) + (-t448 * t525 + t310 - t497) * MDP(25) + ((t289 + t498) * t400 + t397 * t449 + t480) * MDP(26) + (t428 + t493) * MDP(27) + (-t280 * t347 + t433 * t397 + t400 * t518 - t442) * MDP(28); -t347 * t345 * MDP(11) + (-t345 ^ 2 + t347 ^ 2) * MDP(12) + t413 * MDP(13) + ((-t446 + t347) * qJD(4) + t409) * MDP(14) + qJDD(4) * MDP(15) + (t304 * qJD(4) - t358 * t347 + t412 - t418) * MDP(16) + (-t358 * t345 - t411 - t426) * MDP(17) + (t400 * t449 - t500) * MDP(18) + ((-t289 + t498) * t400 - t326 * t448 + t480) * MDP(19) + (t428 - t493) * MDP(20) + (t427 + t497) * MDP(21) - t525 * t347 * MDP(22) + (-pkin(4) * t290 - t278 * t347 - t304 * t324 + (-t459 - t277 + (-t317 - t505) * t525) * t400 + (t524 * t525 + t419) * t397 + t475) * MDP(23) + (pkin(4) * t289 + t479 * t525 + t279 * t347 - t304 * t326 + t419 * t400 + (t277 - t412 + t458) * t397) * MDP(24) + (t269 * t347 + t275 * t525 - t290 * t429 + t476 * t324 + (t445 - t459) * t400 + t517 * t397 + t475) * MDP(25) + (t274 * t324 - t275 * t326 + ((qJD(5) * t326 - t290) * pkin(8) + t433) * t400 + ((qJD(5) * t324 - t289) * pkin(8) - t518) * t397 + t411) * MDP(26) + (-t270 * t347 - t274 * t525 - t289 * t429 - t476 * t326 - t517 * t400 + (t412 + t445) * t397) * MDP(27) + (-t269 * t275 - t270 * t274 + t476 * t280 + (qJD(5) * t436 + t411 + t437) * pkin(8) + (-t265 + t412) * t429) * MDP(28); MDP(18) * t494 + (-t324 ^ 2 + t515) * MDP(19) + t271 * MDP(20) + (t397 * t527 - t440 + t449) * MDP(21) + t315 * MDP(22) + (-t300 * t326 + t410 + t501) * MDP(23) + (t278 * t525 + t300 * t324 - t405) * MDP(24) + (-t292 * t324 - t406 + t501 + 0.2e1 * t512) * MDP(25) + (pkin(5) * t289 - qJ(6) * t290 + (t270 - t279) * t326 + (t269 - t465) * t324) * MDP(26) + (0.2e1 * t503 - t280 * t324 + t292 * t326 + (0.2e1 * qJD(6) - t278) * t525 + t405) * MDP(27) + (t263 * qJ(6) - t264 * pkin(5) - t280 * t292 - t269 * t279 - g(1) * (-pkin(5) * t340 + qJ(6) * t341) - g(2) * (pkin(5) * t342 - qJ(6) * t343) + t438 * t508 + t465 * t270) * MDP(28); (t494 - t315) * MDP(25) + t271 * MDP(26) + (-t525 ^ 2 - t515) * MDP(27) + (t406 - t512 - t528) * MDP(28);];
tau  = t1;
