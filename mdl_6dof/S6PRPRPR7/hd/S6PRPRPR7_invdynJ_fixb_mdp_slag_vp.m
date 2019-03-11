% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:56
% EndTime: 2019-03-08 19:54:02
% DurationCPUTime: 4.14s
% Computational Cost: add. (1601->415), mult. (3262->533), div. (0->0), fcn. (2360->10), ass. (0->188)
t371 = cos(qJ(4));
t453 = qJD(2) * qJD(4);
t428 = t371 * t453;
t368 = sin(qJ(4));
t449 = qJDD(2) * t368;
t395 = t428 + t449;
t369 = sin(qJ(2));
t364 = sin(pkin(6));
t471 = qJD(1) * t364;
t438 = t369 * t471;
t468 = qJD(2) * t368;
t417 = pkin(4) * t468 + t438;
t494 = qJ(5) * t371;
t424 = qJ(3) - t494;
t291 = qJD(2) * t424 + t417;
t469 = qJD(2) * qJ(3);
t324 = t438 + t469;
t467 = qJD(2) * t371;
t351 = qJD(6) + t467;
t367 = sin(qJ(6));
t370 = cos(qJ(6));
t472 = MDP(25) * t370;
t508 = (t351 * (MDP(24) * t367 + t472) - t291 * MDP(18) - t324 * MDP(7)) * qJD(2);
t374 = -pkin(2) - pkin(8);
t372 = cos(qJ(2));
t437 = t372 * t471;
t415 = qJD(3) - t437;
t312 = qJD(2) * t374 + t415;
t366 = cos(pkin(6));
t470 = qJD(1) * t366;
t286 = -t371 * t312 + t368 * t470;
t507 = qJD(5) + t286;
t443 = MDP(13) - MDP(16);
t442 = MDP(14) - MDP(17);
t441 = MDP(15) * qJDD(2);
t458 = qJD(4) * t371;
t460 = qJD(4) * t368;
t506 = pkin(4) * t458 + qJ(5) * t460;
t273 = -qJD(4) * pkin(4) + t507;
t429 = t368 * t453;
t444 = t371 * qJDD(2);
t394 = t429 - t444;
t318 = -qJDD(6) + t394;
t303 = t370 * t318;
t456 = qJD(6) * t367;
t396 = -t351 * t456 - t303;
t287 = t368 * t312 + t371 * t470;
t466 = qJD(4) * qJ(5);
t274 = -t466 - t287;
t363 = sin(pkin(10));
t365 = cos(pkin(10));
t482 = t366 * t369;
t307 = t363 * t372 + t365 * t482;
t309 = -t363 * t482 + t365 * t372;
t505 = -g(1) * t309 - g(2) * t307;
t454 = pkin(5) * t467 + t507;
t504 = qJDD(2) * t374;
t485 = t364 * t369;
t436 = qJD(2) * t485;
t331 = qJD(1) * t436;
t452 = qJDD(1) * t364;
t427 = t372 * t452;
t398 = qJDD(3) + t331 - t427;
t288 = t398 + t504;
t503 = -t288 * t371 + qJDD(5);
t272 = -pkin(5) * t468 + t287;
t270 = t272 + t466;
t373 = -pkin(4) - pkin(9);
t502 = -t373 * t318 + (t270 - t272) * t351;
t360 = t368 * pkin(4);
t325 = t360 + t424;
t445 = qJDD(4) * t374;
t500 = qJD(4) * (-qJD(2) * t325 - t291 + t438) - t445;
t499 = qJD(4) * (-t324 + t438 - t469) - t445;
t430 = qJD(4) * t470;
t451 = qJDD(1) * t366;
t439 = t312 * t460 + t368 * t451 + t371 * t430;
t391 = t439 + t503;
t492 = qJDD(4) * pkin(4);
t264 = t391 - t492;
t422 = -t312 * t458 - t371 * t451 + (-t288 + t430) * t368;
t448 = qJDD(4) * qJ(5);
t385 = qJD(4) * qJD(5) - t422 + t448;
t498 = qJD(4) * (t273 * t368 - t274 * t371) + t385 * t368 - t264 * t371;
t495 = pkin(5) - t374;
t493 = qJDD(2) * pkin(2);
t260 = -pkin(5) * t395 + t385;
t491 = t260 * t368;
t432 = t370 * t468;
t421 = qJD(6) * t432 + t370 * qJDD(4) + t395 * t367;
t268 = -qJD(4) * t456 + t421;
t490 = t268 * t370;
t461 = qJD(4) * t367;
t320 = -t432 + t461;
t488 = t320 * t351;
t459 = qJD(4) * t370;
t322 = t367 * t468 + t459;
t487 = t322 * t351;
t486 = t364 * t368;
t484 = t364 * t371;
t483 = t364 * t372;
t481 = t366 * t372;
t480 = t367 * t318;
t479 = t367 * t368;
t478 = t367 * t371;
t477 = t368 * t369;
t476 = t370 * t371;
t323 = pkin(4) * t467 + qJ(5) * t468;
t361 = t368 ^ 2;
t362 = t371 ^ 2;
t475 = t361 - t362;
t474 = t361 + t362;
t375 = qJD(4) ^ 2;
t376 = qJD(2) ^ 2;
t473 = t375 + t376;
t465 = qJD(4) * t274;
t464 = qJD(4) * t287;
t463 = qJD(4) * t320;
t462 = qJD(4) * t322;
t457 = qJD(5) * t371;
t455 = qJD(6) * t370;
t450 = qJDD(2) * qJ(3);
t447 = qJDD(4) * t368;
t446 = qJDD(4) * t371;
t440 = t368 * t371 * t376;
t435 = t351 * t461;
t434 = t351 * t459;
t431 = g(3) * (pkin(2) * t483 + qJ(3) * t485);
t426 = qJD(4) * t495;
t259 = -pkin(5) * t394 + qJDD(4) * t373 + t391;
t393 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t371;
t401 = pkin(9) * t368 + t424;
t339 = t369 * t452;
t420 = t395 * pkin(4) + qJ(5) * t429 + t339;
t265 = t401 * qJDD(2) + (t393 + t437) * qJD(2) + t420;
t423 = t370 * t259 - t367 * t265;
t306 = t363 * t369 - t365 * t481;
t308 = t363 * t481 + t365 * t369;
t419 = -g(1) * t308 - g(2) * t306;
t416 = -t494 + t360;
t414 = qJD(3) + t437;
t413 = qJDD(4) * t367 - t395 * t370;
t267 = qJD(4) * t373 + t454;
t281 = qJD(2) * t401 + t417;
t261 = t267 * t370 - t281 * t367;
t262 = t267 * t367 + t281 * t370;
t408 = t351 * t367;
t403 = -g(3) * t483 - t419;
t310 = t366 * t368 + t371 * t483;
t311 = t366 * t371 - t368 * t483;
t400 = t367 * t372 + t369 * t476;
t399 = -t369 * t478 + t370 * t372;
t397 = t351 * t455 - t480;
t276 = t308 * t368 + t363 * t484;
t278 = -t306 * t368 + t365 * t484;
t392 = -g(1) * t276 + g(2) * t278 - g(3) * t311;
t388 = g(3) * t485 - t505;
t387 = -t339 + t388;
t386 = t403 + t427;
t275 = -t308 * t371 + t363 * t486;
t277 = t306 * t371 + t365 * t486;
t384 = g(1) * t275 - g(2) * t277 + g(3) * t310 - t439;
t383 = -t374 * t375 - t388;
t382 = qJDD(3) - t386;
t381 = t392 - t422;
t380 = t260 + (pkin(9) * t467 - qJD(6) * t373 + t323) * t351 + t392;
t379 = t291 * t467 - t384 + t503;
t266 = t424 * qJDD(2) + (t414 - t457) * qJD(2) + t420;
t299 = qJD(3) - t457 + t506;
t378 = -qJDD(2) * t325 - t266 + (-t299 + t437) * qJD(2) - t383;
t289 = qJD(2) * t414 + t339 + t450;
t377 = qJD(2) * t415 + t289 + t383 + t450;
t327 = t495 * t371;
t326 = t495 * t368;
t319 = -qJD(2) * pkin(2) + t415;
t317 = t371 * t426;
t316 = t368 * t426;
t315 = t360 + t401;
t302 = t308 * pkin(2);
t301 = t306 * pkin(2);
t292 = t398 - t493;
t290 = t393 + t506;
t283 = t310 * t367 + t370 * t485;
t282 = t310 * t370 - t367 * t485;
t280 = qJD(4) * t311 - t371 * t436;
t279 = qJD(4) * t310 - t368 * t436;
t269 = qJD(6) * t322 + t413;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (qJDD(1) * t366 ^ 2 - g(3)) * MDP(7) + (t264 * t310 + t273 * t280 + t274 * t279 + t311 * t385 - g(3)) * MDP(18) + ((t280 * t370 - t310 * t456) * t351 - t282 * t318 - t279 * t320 + t311 * t269) * MDP(24) + (-(t280 * t367 + t310 * t455) * t351 + t283 * t318 - t279 * t322 + t311 * t268) * MDP(25) + (t310 * t371 - t311 * t368) * t441 + (t279 * t368 + t280 * t371 + (-t310 * t368 - t311 * t371) * qJD(4)) * MDP(15) * qJD(2) + t442 * (qJD(4) * t279 - qJDD(4) * t311) - t443 * (qJD(4) * t280 + qJDD(4) * t310) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t369 + t372 * t376) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t372 + t369 * t376) + (-t292 * MDP(7) - t508 + (t368 * t443 + t371 * t442) * t376) * t372 + ((qJD(2) * t319 + t289) * MDP(7) + t266 * MDP(18) + (-MDP(24) * t370 + MDP(25) * t367) * t351 * qJD(6) + t443 * t395 - t442 * t394) * t369) * t364; qJDD(2) * MDP(2) + t386 * MDP(3) + t387 * MDP(4) + (t382 - 0.2e1 * t493) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) - t387 + 0.2e1 * t450) * MDP(6) + (t289 * qJ(3) + t324 * qJD(3) - t292 * pkin(2) - g(1) * (qJ(3) * t309 - t302) - g(2) * (qJ(3) * t307 - t301) - t431 + (-t319 * t369 - t324 * t372) * t471) * MDP(7) + (qJDD(2) * t362 - 0.2e1 * t368 * t428) * MDP(8) + 0.2e1 * (-t368 * t444 + t453 * t475) * MDP(9) + (-t368 * t375 + t446) * MDP(10) + (-t371 * t375 - t447) * MDP(11) + (t377 * t368 - t371 * t499) * MDP(13) + (t368 * t499 + t377 * t371) * MDP(14) + (t403 + (t331 - t504) * t474 - t498) * MDP(15) + (t378 * t368 + t371 * t500) * MDP(16) + (-t368 * t500 + t378 * t371) * MDP(17) + (t266 * t325 + t291 * t299 - g(1) * (-pkin(8) * t308 - t302) - g(2) * (-pkin(8) * t306 - t301) - t431 + t498 * t374 + (-g(3) * (pkin(8) * t372 + t369 * t416) + (-t291 * t372 + (t273 * t371 + t274 * t368) * t369) * qJD(1)) * t364 + t505 * (qJ(3) + t416)) * MDP(18) + (t268 * t479 + (t367 * t458 + t368 * t455) * t322) * MDP(19) + ((-t320 * t367 + t322 * t370) * t458 + (t490 - t269 * t367 + (-t320 * t370 - t322 * t367) * qJD(6)) * t368) * MDP(20) + ((t268 + t435) * t371 + (t397 - t462) * t368) * MDP(21) + ((-t269 + t434) * t371 + (t396 + t463) * t368) * MDP(22) + (-t318 * t371 - t351 * t460) * MDP(23) + ((-t290 * t367 - t316 * t370) * t351 - (-t315 * t367 + t327 * t370) * t318 + t423 * t371 - t317 * t320 - t326 * t269 - t370 * t491 - g(1) * (-t308 * t370 - t309 * t478) - g(2) * (-t306 * t370 - t307 * t478) + (-t261 * t368 - t270 * t476) * qJD(4) + ((-t315 * t370 - t327 * t367) * t351 - t262 * t371 + t270 * t479) * qJD(6) + (-g(3) * t399 + (-t320 * t477 + t351 * t400) * qJD(1)) * t364) * MDP(24) + (t262 * t460 - t326 * t268 - t317 * t322 + (-(qJD(6) * t327 + t290) * t351 + t315 * t318 + t270 * qJD(6) * t368 + (-qJD(6) * t267 - t265 - t505) * t371) * t370 + (-(-qJD(6) * t315 - t316) * t351 + t327 * t318 + t491 + (qJD(4) * t270 + qJD(6) * t281 - t259) * t371 + t419) * t367 + (g(3) * t400 + (-t322 * t477 + t351 * t399) * qJD(1)) * t364) * MDP(25); -t376 * MDP(6) + (t331 + t382) * MDP(7) - t403 * MDP(18) + t443 * (-t368 * t473 + t446) - t442 * (t371 * t473 + t447) + t508 + ((qJD(4) * t273 + t385) * MDP(18) + (t269 + t434) * MDP(24) + (t268 - t435) * MDP(25)) * t368 + (-MDP(15) * t474 - pkin(2) * MDP(7) + MDP(5)) * qJDD(2) + ((-t264 - t465) * MDP(18) + (-t396 + t463) * MDP(24) + (t397 + t462) * MDP(25)) * t371; MDP(8) * t440 - t475 * t376 * MDP(9) + MDP(10) * t444 - MDP(11) * t449 + qJDD(4) * MDP(12) + (t464 + (-qJD(2) * t324 + t288) * t371 + t384) * MDP(13) + (-qJD(4) * t286 + t324 * t468 - t381) * MDP(14) + (-pkin(4) * t371 - qJ(5) * t368) * t441 + (t323 * t468 + t379 - t464 - 0.2e1 * t492) * MDP(16) + (0.2e1 * t448 + (0.2e1 * qJD(5) + t286) * qJD(4) + (-t291 * t368 + t323 * t371) * qJD(2) + t381) * MDP(17) + (t385 * qJ(5) - t264 * pkin(4) - t291 * t323 - t273 * t287 - g(1) * (-pkin(4) * t275 + qJ(5) * t276) - g(2) * (pkin(4) * t277 - qJ(5) * t278) - g(3) * (-pkin(4) * t310 + qJ(5) * t311) - t507 * t274) * MDP(18) + (-t322 * t408 + t490) * MDP(19) + ((-t269 - t487) * t370 + (-t268 + t488) * t367) * MDP(20) + ((t322 * t368 - t351 * t478) * qJD(2) + t396) * MDP(21) + ((-t320 * t368 - t351 * t476) * qJD(2) - t397) * MDP(22) + t351 * MDP(23) * t468 + (qJ(5) * t269 + t261 * t468 + t454 * t320 + t380 * t367 + t370 * t502) * MDP(24) + (qJ(5) * t268 - t262 * t468 + t454 * t322 - t367 * t502 + t380 * t370) * MDP(25); t371 * t441 + (qJDD(4) - t440) * MDP(16) + (-t362 * t376 - t375) * MDP(17) + (t379 + t465 - t492) * MDP(18) + (-t303 - t463) * MDP(24) + (-t462 + t480) * MDP(25) + (-MDP(24) * t408 - t351 * t472) * t351; t322 * t320 * MDP(19) + (-t320 ^ 2 + t322 ^ 2) * MDP(20) + (t421 + t488) * MDP(21) + (-t413 + t487) * MDP(22) - t318 * MDP(23) + (t262 * t351 - t270 * t322 - g(1) * (t275 * t370 - t309 * t367) - g(2) * (-t277 * t370 - t307 * t367) - g(3) * t282 + t423) * MDP(24) + (-t370 * t265 - t367 * t259 + t261 * t351 + t270 * t320 - g(1) * (-t275 * t367 - t309 * t370) - g(2) * (t277 * t367 - t307 * t370) + g(3) * t283) * MDP(25) + (-MDP(21) * t461 - MDP(22) * t322 - MDP(24) * t262 - MDP(25) * t261) * qJD(6);];
tau  = t1;
