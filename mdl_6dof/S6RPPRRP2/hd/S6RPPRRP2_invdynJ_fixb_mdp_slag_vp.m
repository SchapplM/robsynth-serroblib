% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:27
% EndTime: 2019-03-09 02:01:34
% DurationCPUTime: 5.42s
% Computational Cost: add. (4777->463), mult. (10047->560), div. (0->0), fcn. (7310->14), ass. (0->188)
t408 = sin(pkin(9));
t385 = pkin(1) * t408 + qJ(3);
t372 = qJD(1) * qJD(3) + qJDD(1) * t385;
t409 = cos(pkin(10));
t392 = t409 * qJDD(2);
t407 = sin(pkin(10));
t338 = t392 + (-pkin(7) * qJDD(1) - t372) * t407;
t351 = t407 * qJDD(2) + t409 * t372;
t464 = qJDD(1) * t409;
t339 = pkin(7) * t464 + t351;
t413 = sin(qJ(4));
t416 = cos(qJ(4));
t379 = t385 * qJD(1);
t394 = t409 * qJD(2);
t344 = t394 + (-pkin(7) * qJD(1) - t379) * t407;
t358 = t407 * qJD(2) + t409 * t379;
t472 = qJD(1) * t409;
t345 = pkin(7) * t472 + t358;
t305 = t413 * t344 + t416 * t345;
t528 = qJD(4) * t305;
t438 = t338 * t416 - t413 * t339 - t528;
t279 = -qJDD(4) * pkin(4) - t438;
t412 = sin(qJ(5));
t415 = cos(qJ(5));
t485 = t407 * t413;
t374 = -t416 * t409 + t485;
t366 = t374 * qJD(4);
t428 = qJD(1) * t366;
t375 = t407 * t416 + t409 * t413;
t429 = t375 * qJDD(1);
t420 = t429 - t428;
t468 = t415 * qJD(4);
t470 = qJD(5) * t412;
t518 = t375 * qJD(1);
t297 = -qJD(5) * t468 - t412 * qJDD(4) - t415 * t420 + t470 * t518;
t430 = t374 * qJD(1);
t469 = qJD(5) * t415;
t471 = qJD(4) * t412;
t298 = -t415 * qJDD(4) + t412 * t429 + t469 * t518 + t471 * (qJD(5) - t430);
t348 = t415 * t518 + t471;
t270 = pkin(5) * t298 + qJ(6) * t297 - qJD(6) * t348 + t279;
t405 = pkin(10) + qJ(4);
t397 = cos(t405);
t506 = g(3) * t397;
t395 = sin(t405);
t406 = qJ(1) + pkin(9);
t396 = sin(t406);
t398 = cos(t406);
t453 = g(1) * t398 + g(2) * t396;
t523 = t395 * t453;
t425 = -t506 + t523;
t529 = -t270 + t425;
t473 = qJD(1) * t407;
t364 = -t413 * t473 + t416 * t472;
t360 = qJD(5) - t364;
t445 = t338 * t413 + t339 * t416;
t520 = t344 * t416 - t413 * t345;
t278 = qJDD(4) * pkin(8) + qJD(4) * t520 + t445;
t367 = t375 * qJD(4);
t383 = t416 * t464;
t448 = qJDD(1) * t485 - t383;
t332 = qJD(1) * t367 + t448;
t388 = t409 * pkin(3) + pkin(2);
t410 = cos(pkin(9));
t505 = t410 * pkin(1);
t377 = -t388 - t505;
t519 = -pkin(8) * t375 + t377;
t291 = t332 * pkin(4) + pkin(8) * t428 + qJDD(1) * t519 + qJDD(3);
t302 = qJD(4) * pkin(8) + t305;
t362 = qJD(1) * t377 + qJD(3);
t312 = -pkin(4) * t364 - pkin(8) * t518 + t362;
t433 = t415 * t278 + t412 * t291 - t302 * t470 + t312 * t469;
t329 = qJDD(5) + t332;
t502 = qJ(6) * t329;
t268 = qJD(6) * t360 + t433 + t502;
t280 = -t302 * t412 + t312 * t415;
t467 = qJD(6) - t280;
t274 = -pkin(5) * t360 + t467;
t527 = t360 * t274 + t268;
t456 = t412 * t278 - t415 * t291 + t302 * t469 + t312 * t470;
t513 = pkin(5) * t329;
t269 = qJDD(6) + t456 - t513;
t281 = t302 * t415 + t312 * t412;
t275 = qJ(6) * t360 + t281;
t501 = t275 * t360;
t526 = -t269 + t501;
t458 = -g(1) * t396 + g(2) * t398;
t522 = t458 * t395;
t346 = t412 * t518 - t468;
t521 = t374 * t298 + t367 * t346;
t480 = -t374 * t297 + t348 * t367;
t319 = pkin(4) * t374 + t519;
t504 = pkin(7) + t385;
t368 = t504 * t407;
t369 = t504 * t409;
t325 = -t368 * t413 + t369 * t416;
t479 = t412 * t319 + t415 * t325;
t517 = MDP(21) + MDP(23);
t301 = -qJD(4) * pkin(4) - t520;
t284 = pkin(5) * t346 - qJ(6) * t348 + t301;
t512 = pkin(8) * t329;
t516 = t284 * t360 - t512;
t507 = g(3) * t395;
t426 = -t453 * t397 - t507;
t515 = t348 ^ 2;
t414 = sin(qJ(1));
t514 = pkin(1) * t414;
t503 = pkin(8) * qJD(5);
t500 = t281 * t360;
t499 = t297 * t412;
t497 = t346 * t348;
t496 = t346 * t364;
t495 = t346 * t412;
t457 = t348 * t360;
t494 = t348 * t415;
t493 = t360 * t412;
t492 = t366 * t412;
t491 = t366 * t415;
t490 = t375 * t415;
t489 = t396 * t412;
t488 = t396 * t415;
t487 = t398 * t412;
t486 = t398 * t415;
t320 = t412 * t329;
t321 = t415 * t329;
t484 = -t298 * t490 + t346 * t491;
t483 = -t412 * t298 - t346 * t469;
t330 = pkin(4) * t518 - pkin(8) * t364;
t482 = t412 * t330 + t415 * t520;
t328 = t360 * t491;
t481 = t375 * t321 - t328;
t478 = t360 * t469 + t320;
t477 = t364 * t493 + t321;
t449 = pkin(5) * t412 - qJ(6) * t415;
t476 = -qJD(6) * t412 + t360 * t449 - t305;
t475 = (g(1) * t486 + g(2) * t488) * t395;
t474 = t407 ^ 2 + t409 ^ 2;
t389 = -pkin(2) - t505;
t465 = qJDD(1) * t389;
t462 = MDP(22) - MDP(25);
t461 = t360 * t503;
t460 = t348 * t492;
t459 = t375 * t470;
t353 = t397 * t489 + t486;
t355 = t397 * t487 - t488;
t455 = -g(1) * t353 + g(2) * t355;
t354 = t397 * t488 - t487;
t356 = t397 * t486 + t489;
t454 = g(1) * t354 - g(2) * t356;
t417 = cos(qJ(1));
t451 = g(1) * t414 - g(2) * t417;
t450 = pkin(5) * t415 + qJ(6) * t412;
t447 = t268 * t415 + t269 * t412;
t446 = t274 * t415 - t275 * t412;
t350 = -t372 * t407 + t392;
t443 = -t350 * t407 + t351 * t409;
t442 = -t368 * t416 - t369 * t413;
t441 = pkin(4) + t450;
t440 = pkin(4) * t397 + pkin(8) * t395 + t388;
t437 = t461 + t506;
t436 = t375 * t469 - t492;
t435 = t459 + t491;
t434 = t301 * t360 - t512;
t308 = -qJD(3) * t374 + qJD(4) * t442;
t331 = pkin(4) * t367 + pkin(8) * t366;
t432 = t415 * t308 + t319 * t469 - t325 * t470 + t412 * t331;
t424 = g(1) * t355 + g(2) * t353 + t412 * t507 - t456;
t422 = t284 * t348 + qJDD(6) - t424;
t421 = -g(1) * t356 - g(2) * t354 - t415 * t507 + t433;
t309 = qJD(3) * t375 + qJD(4) * t325;
t411 = -pkin(7) - qJ(3);
t402 = t417 * pkin(1);
t376 = qJDD(3) + t465;
t359 = qJDD(1) * t377 + qJDD(3);
t357 = -t379 * t407 + t394;
t334 = -qJD(4) * t367 - qJDD(4) * t374;
t333 = -qJD(4) * t366 + qJDD(4) * t375;
t310 = pkin(5) * t348 + qJ(6) * t346;
t296 = t375 * t449 - t442;
t287 = -pkin(5) * t374 - t319 * t415 + t325 * t412;
t286 = qJ(6) * t374 + t479;
t285 = t346 * t360 - t297;
t283 = -pkin(5) * t518 - t330 * t415 + t412 * t520;
t282 = qJ(6) * t518 + t482;
t273 = -t449 * t366 + (qJD(5) * t450 - qJD(6) * t415) * t375 + t309;
t272 = -pkin(5) * t367 + qJD(5) * t479 + t308 * t412 - t331 * t415;
t271 = qJ(6) * t367 + qJD(6) * t374 + t432;
t1 = [(t268 * t286 + t275 * t271 + t270 * t296 + t284 * t273 + t269 * t287 + t274 * t272 - g(1) * (-pkin(5) * t354 - qJ(6) * t353 - t514) - g(2) * (pkin(5) * t356 + qJ(6) * t355 + t402) + (g(1) * t411 - g(2) * t440) * t398 + (g(1) * t440 + g(2) * t411) * t396) * MDP(26) + (t376 * t389 - g(1) * (-pkin(2) * t396 + qJ(3) * t398 - t514) - g(2) * (pkin(2) * t398 + qJ(3) * t396 + t402) + t443 * t385 + (-t357 * t407 + t358 * t409) * qJD(3)) * MDP(8) + (t460 + (t499 + (-t494 + t495) * qJD(5)) * t375 + t484) * MDP(17) + (-t297 * t490 - t348 * t435) * MDP(16) + (t268 * t374 - t270 * t490 + t271 * t360 - t273 * t348 + t275 * t367 + t284 * t435 + t286 * t329 + t296 * t297 - t455) * MDP(25) + (-t360 * t459 + t480 + t481) * MDP(18) + (t279 * t490 - t281 * t367 + t297 * t442 - t301 * t435 + t309 * t348 - t329 * t479 - t360 * t432 - t374 * t433 + t455) * MDP(22) + (-qJD(4) * t309 + qJDD(4) * t442 + t332 * t377 + t359 * t374 + t362 * t367 - t397 * t458) * MDP(14) + (-t456 * t374 + t280 * t367 + t309 * t346 - t442 * t298 + ((-qJD(5) * t325 + t331) * t360 + t319 * t329 + t301 * qJD(5) * t375) * t415 + ((-qJD(5) * t319 - t308) * t360 - t325 * t329 + t279 * t375 - t301 * t366) * t412 + t454) * MDP(21) + (-t375 * t332 - t366 * t364 - t367 * t518 - t374 * t420) * MDP(10) + (-t366 * t518 + t375 * t420) * MDP(9) + (t372 * t474 + t443 - t453) * MDP(7) + t451 * MDP(2) + (t270 * t375 * t412 - t269 * t374 - t272 * t360 + t273 * t346 - t274 * t367 + t284 * t436 - t287 * t329 + t296 * t298 + t454) * MDP(23) + qJDD(1) * MDP(1) + (t451 + (t408 ^ 2 + t410 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (g(1) * t417 + g(2) * t414) * MDP(3) + (t329 * t374 + t360 * t367) * MDP(20) + t333 * MDP(11) + t334 * MDP(12) + (-t320 * t375 - t360 * t436 - t521) * MDP(19) + (-t325 * qJDD(4) + t359 * t375 - t362 * t366 + t522 + t377 * t429 + (-t377 * t430 - t308) * qJD(4)) * MDP(15) + (-t271 * t346 + t272 * t348 - t286 * t298 - t287 * t297 - t522 - t446 * t366 + (-t268 * t412 + t269 * t415 + (-t274 * t412 - t275 * t415) * qJD(5)) * t375) * MDP(24) + (t409 * MDP(5) - t407 * MDP(6)) * (-t376 - t458 - t465); (qJDD(2) - g(3)) * MDP(4) + (t350 * t409 + t351 * t407 - g(3)) * MDP(8) + t334 * MDP(14) - t333 * MDP(15) + (t328 + t480) * MDP(22) + (-t460 + t484) * MDP(24) + (t481 - t480) * MDP(25) + (t270 * t374 - t274 * t492 - t275 * t491 + t284 * t367 - g(3)) * MDP(26) + (-MDP(24) * t499 + t447 * MDP(26) + (-t415 * MDP(22) - t412 * t517) * t329 + ((t494 + t495) * MDP(24) + t446 * MDP(26) + (t412 * t462 - t415 * t517) * t360) * qJD(5)) * t375 + t517 * (t360 * t492 + t521); (t357 * t473 - t358 * t472 + qJDD(3) + t458) * MDP(8) - t383 * MDP(14) + t477 * MDP(21) + t483 * MDP(24) + t478 * MDP(25) + t458 * MDP(26) - t474 * MDP(7) * qJD(1) ^ 2 + (-t284 * MDP(26) - t346 * t517 - t348 * t462) * t518 + (t329 * MDP(23) + (t297 + t496) * MDP(24) + t526 * MDP(26) + (-MDP(22) * t360 - t364 * MDP(25)) * t360) * t415 + (-t329 * MDP(22) + t527 * MDP(26) + MDP(24) * t457 + (-qJD(5) * MDP(21) - t360 * MDP(23)) * t360) * t412 + (t389 * MDP(8) + (MDP(15) * t413 - MDP(5)) * t409 + (MDP(14) * t413 + MDP(15) * t416 + MDP(6)) * t407) * qJDD(1) + (t518 * MDP(14) + t364 * MDP(15) + (MDP(14) * t375 - MDP(15) * t374) * qJD(1)) * qJD(4); -t364 ^ 2 * MDP(10) + (t429 + (-t364 - t430) * qJD(4)) * MDP(11) - t448 * MDP(12) + qJDD(4) * MDP(13) + (t425 + t438 + t528) * MDP(14) + (-t362 * t364 - t426 - t445) * MDP(15) + (t415 * t457 - t499) * MDP(16) + ((-t297 + t496) * t415 - t348 * t493 + t483) * MDP(17) + (-t360 * t364 * t415 + t478) * MDP(18) + (-t360 * t470 + t477) * MDP(19) + (-pkin(4) * t298 - t305 * t346 + (-t506 - t279 + (-t330 - t503) * t360) * t415 + (t360 * t520 + t434) * t412 + t475) * MDP(21) + (pkin(4) * t297 + t482 * t360 - t305 * t348 + t434 * t415 + (t279 + t437 - t523) * t412) * MDP(22) + (t283 * t360 - t298 * t441 + t476 * t346 + (-t270 - t437) * t415 + t516 * t412 + t475) * MDP(23) + (t282 * t346 - t283 * t348 + ((qJD(5) * t348 - t298) * pkin(8) + t527) * t415 + ((qJD(5) * t346 - t297) * pkin(8) - t526) * t412 + t426) * MDP(24) + (-t282 * t360 - t297 * t441 - t476 * t348 - t516 * t415 + (-t461 + t529) * t412) * MDP(25) + (-t274 * t283 - t275 * t282 + t476 * t284 + (qJD(5) * t446 + t426 + t447) * pkin(8) + t529 * t441) * MDP(26) + (MDP(10) * t518 - t362 * MDP(14) - t348 * MDP(18) + t346 * MDP(19) - t360 * MDP(20) - t280 * MDP(21) + t281 * MDP(22) + t274 * MDP(23) - t275 * MDP(25) - t364 * MDP(9)) * t518; MDP(16) * t497 + (-t346 ^ 2 + t515) * MDP(17) + t285 * MDP(18) + (-t298 + t457) * MDP(19) + t329 * MDP(20) + (-t301 * t348 + t424 + t500) * MDP(21) + (t280 * t360 + t301 * t346 - t421) * MDP(22) + (-t310 * t346 - t422 + t500 + 0.2e1 * t513) * MDP(23) + (pkin(5) * t297 - qJ(6) * t298 + (t275 - t281) * t348 + (t274 - t467) * t346) * MDP(24) + (0.2e1 * t502 - t284 * t346 + t310 * t348 + (0.2e1 * qJD(6) - t280) * t360 + t421) * MDP(25) + (t268 * qJ(6) - t269 * pkin(5) - t284 * t310 - t274 * t281 - g(1) * (-pkin(5) * t355 + qJ(6) * t356) - g(2) * (-pkin(5) * t353 + qJ(6) * t354) + t449 * t507 + t467 * t275) * MDP(26); (-qJD(4) * t518 - qJDD(5) - t448 + t497) * MDP(23) + t285 * MDP(24) + (-t360 ^ 2 - t515) * MDP(25) + (t422 - t501 - t513) * MDP(26);];
tau  = t1;
