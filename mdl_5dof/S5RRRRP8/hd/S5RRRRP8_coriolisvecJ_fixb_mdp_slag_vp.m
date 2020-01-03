% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:16
% EndTime: 2019-12-31 22:02:25
% DurationCPUTime: 4.11s
% Computational Cost: add. (3111->366), mult. (7694->509), div. (0->0), fcn. (5144->6), ass. (0->168)
t408 = cos(qJ(3));
t407 = sin(qJ(2));
t409 = cos(qJ(2));
t483 = t408 * t409;
t423 = pkin(3) * t407 - pkin(8) * t483;
t498 = pkin(7) + pkin(8);
t449 = qJD(3) * t498;
t425 = pkin(2) * t407 - pkin(7) * t409;
t363 = t425 * qJD(1);
t406 = sin(qJ(3));
t465 = qJD(1) * t407;
t444 = t406 * t465;
t468 = pkin(6) * t444 + t408 * t363;
t519 = qJD(1) * t423 + t408 * t449 + t468;
t343 = t406 * t363;
t485 = t407 * t408;
t486 = t406 * t409;
t518 = -t343 - (-t485 * pkin(6) - pkin(8) * t486) * qJD(1) - t406 * t449;
t458 = qJD(3) * t407;
t439 = qJD(1) * t458;
t453 = qJD(1) * qJD(2);
t440 = t409 * t453;
t515 = qJD(2) * qJD(3) + t440;
t328 = -t406 * t439 + t515 * t408;
t461 = qJD(2) * t408;
t356 = -t444 + t461;
t463 = qJD(2) * t406;
t357 = t408 * t465 + t463;
t405 = sin(qJ(4));
t497 = cos(qJ(4));
t442 = t497 * qJD(4);
t450 = t515 * t406 + t408 * t439;
t456 = qJD(4) * t405;
t278 = -t497 * t328 - t356 * t442 + t357 * t456 + t405 * t450;
t421 = t405 * t356 + t357 * t497;
t279 = qJD(4) * t421 + t405 * t328 + t497 * t450;
t312 = -t497 * t356 + t357 * t405;
t310 = t312 ^ 2;
t464 = qJD(1) * t409;
t391 = -qJD(3) + t464;
t380 = -qJD(4) + t391;
t437 = MDP(22) * t465;
t499 = t421 ^ 2;
t517 = qJD(2) * t437 + (-t380 * t421 - t279) * MDP(21) + t312 * t421 * MDP(18) + (-t312 * t380 - t278) * MDP(20) + (-t310 + t499) * MDP(19);
t516 = t312 * qJ(5);
t489 = t405 * t406;
t420 = t408 * t497 - t489;
t502 = qJD(3) + qJD(4);
t503 = t497 * qJD(3) + t442;
t474 = -t503 * t408 + t420 * t464 + t502 * t489;
t359 = t405 * t408 + t406 * t497;
t323 = t502 * t359;
t473 = -t359 * t464 + t323;
t400 = pkin(6) * t464;
t459 = qJD(3) * t406;
t506 = -t400 + (-t406 * t464 + t459) * pkin(3);
t460 = qJD(2) * t409;
t445 = t406 * t460;
t457 = qJD(3) * t408;
t446 = t407 * t457;
t514 = t445 + t446;
t373 = -qJD(2) * pkin(2) + pkin(6) * t465;
t329 = -pkin(3) * t356 + t373;
t374 = qJD(2) * pkin(7) + t400;
t368 = -pkin(2) * t409 - pkin(7) * t407 - pkin(1);
t349 = t368 * qJD(1);
t488 = t406 * t349;
t320 = t374 * t408 + t488;
t366 = t425 * qJD(2);
t350 = qJD(1) * t366;
t441 = t407 * t453;
t428 = pkin(6) * t441;
t472 = -t408 * t350 - t406 * t428;
t417 = -qJD(3) * t320 - t472;
t276 = pkin(3) * t441 - pkin(8) * t328 + t417;
t422 = t349 * t457 + t406 * t350 - t374 * t459;
t416 = -t408 * t428 + t422;
t281 = -pkin(8) * t450 + t416;
t319 = t408 * t349 - t374 * t406;
t298 = -pkin(8) * t357 + t319;
t293 = -pkin(3) * t391 + t298;
t299 = pkin(8) * t356 + t320;
t427 = -t405 * t276 - t497 * t281 - t293 * t442 + t299 * t456;
t513 = t312 * t329 + t427;
t511 = -0.2e1 * t453;
t510 = MDP(4) * t407;
t403 = t407 ^ 2;
t509 = MDP(5) * (-t409 ^ 2 + t403);
t508 = qJ(5) * t421;
t290 = pkin(4) * t312 + qJD(5) + t329;
t507 = t290 * t421;
t355 = t408 * t368;
t496 = pkin(6) * t406;
t318 = -pkin(8) * t485 + t355 + (-pkin(3) - t496) * t409;
t393 = pkin(6) * t483;
t467 = t406 * t368 + t393;
t487 = t406 * t407;
t324 = -pkin(8) * t487 + t467;
t475 = t405 * t318 + t497 * t324;
t375 = t498 * t406;
t376 = t498 * t408;
t469 = -t405 * t375 + t497 * t376;
t505 = t469 * qJD(4) + t518 * t405 + t519 * t497;
t504 = -t375 * t442 - t376 * t456 - t519 * t405 + t518 * t497;
t297 = t497 * t299;
t269 = t405 * t293 + t297;
t414 = -qJD(4) * t269 + t497 * t276 - t405 * t281;
t501 = -t329 * t421 + t414;
t495 = t328 * t406;
t494 = t356 * t391;
t493 = t357 * t391;
t492 = t373 * t406;
t491 = t373 * t408;
t490 = t391 * t408;
t295 = t405 * t299;
t410 = qJD(2) ^ 2;
t484 = t407 * t410;
t482 = t409 * t410;
t411 = qJD(1) ^ 2;
t481 = t409 * t411;
t268 = t497 * t293 - t295;
t263 = t268 - t508;
t262 = -pkin(4) * t380 + t263;
t480 = t262 - t263;
t479 = -t473 * qJ(5) + qJD(5) * t420 + t504;
t478 = -pkin(4) * t465 + t474 * qJ(5) - t359 * qJD(5) - t505;
t477 = t497 * t298 - t295;
t471 = t406 * t366 + t368 * t457;
t462 = qJD(2) * t407;
t470 = t408 * t366 + t462 * t496;
t367 = pkin(3) * t487 + t407 * pkin(6);
t454 = t407 * MDP(15);
t330 = t514 * pkin(3) + pkin(6) * t460;
t398 = -pkin(3) * t408 - pkin(2);
t448 = t406 * t458;
t447 = t409 * t459;
t438 = qJD(2) * t454;
t436 = pkin(1) * t511;
t435 = -t298 * t405 - t297;
t433 = t497 * t318 - t324 * t405;
t431 = -t497 * t375 - t376 * t405;
t430 = -t356 + t461;
t429 = -t357 + t463;
t426 = t497 * t460;
t424 = qJD(1) * t403 - t391 * t409;
t309 = pkin(3) * t450 + pkin(6) * t440;
t285 = t423 * qJD(2) + (-t393 + (pkin(8) * t407 - t368) * t406) * qJD(3) + t470;
t287 = -t514 * pkin(8) + (-t407 * t461 - t447) * pkin(6) + t471;
t419 = t405 * t285 + t497 * t287 + t318 * t442 - t324 * t456;
t267 = t279 * pkin(4) + t309;
t413 = -t475 * qJD(4) + t497 * t285 - t405 * t287;
t397 = pkin(3) * t497 + pkin(4);
t338 = t420 * t407;
t337 = t359 * t407;
t303 = qJ(5) * t420 + t469;
t302 = -qJ(5) * t359 + t431;
t289 = t406 * t426 - t405 * t448 - t456 * t487 + (t405 * t460 + t503 * t407) * t408;
t288 = t323 * t407 + t405 * t445 - t408 * t426;
t282 = -qJ(5) * t337 + t475;
t280 = -pkin(4) * t409 - qJ(5) * t338 + t433;
t266 = t477 - t508;
t265 = t435 + t516;
t264 = t269 - t516;
t261 = -qJ(5) * t289 - qJD(5) * t337 + t419;
t260 = pkin(4) * t462 + t288 * qJ(5) - t338 * qJD(5) + t413;
t259 = -qJ(5) * t279 - qJD(5) * t312 - t427;
t258 = pkin(4) * t441 + t278 * qJ(5) - qJD(5) * t421 + t414;
t1 = [0.2e1 * t440 * t510 + t509 * t511 + MDP(6) * t482 - MDP(7) * t484 + (-pkin(6) * t482 + t407 * t436) * MDP(9) + (pkin(6) * t484 + t409 * t436) * MDP(10) + (t328 * t485 + (t408 * t460 - t448) * t357) * MDP(11) + ((t356 * t408 - t357 * t406) * t460 + (-t408 * t450 - t495 + (-t406 * t356 - t357 * t408) * qJD(3)) * t407) * MDP(12) + (t391 * t448 - t328 * t409 + (t357 * t407 + t408 * t424) * qJD(2)) * MDP(13) + (t391 * t446 + t450 * t409 + (t356 * t407 - t406 * t424) * qJD(2)) * MDP(14) + (-t391 - t464) * t438 + (-(-t368 * t459 + t470) * t391 + (pkin(6) * t450 + t373 * t457 + (t355 * qJD(1) + t319) * qJD(2)) * t407 + ((-pkin(6) * t356 + t492) * qJD(2) + (t488 + (pkin(6) * t391 + t374) * t408) * qJD(3) + t472) * t409) * MDP(16) + ((-pkin(6) * t447 + t471) * t391 + t422 * t409 + (pkin(6) * t328 - t373 * t459) * t407 + ((pkin(6) * t357 + t491) * t409 + (-pkin(6) * t490 - qJD(1) * t467 - t320) * t407) * qJD(2)) * MDP(17) + (-t278 * t338 - t288 * t421) * MDP(18) + (t278 * t337 - t279 * t338 + t288 * t312 - t289 * t421) * MDP(19) + (t278 * t409 + t288 * t380 + (qJD(1) * t338 + t421) * t462) * MDP(20) + (t279 * t409 + t289 * t380 + (-qJD(1) * t337 - t312) * t462) * MDP(21) + (-t380 - t464) * MDP(22) * t462 + (t268 * t462 + t367 * t279 + t329 * t289 + t309 * t337 + t330 * t312 - t380 * t413 - t409 * t414 + t433 * t441) * MDP(23) + (t419 * t380 - t427 * t409 + t330 * t421 - t367 * t278 + t309 * t338 - t329 * t288 + (-qJD(1) * t475 - t269) * t462) * MDP(24) + (-t258 * t338 - t259 * t337 - t260 * t421 - t261 * t312 + t262 * t288 - t264 * t289 + t278 * t280 - t279 * t282) * MDP(25) + (t259 * t282 + t264 * t261 + t258 * t280 + t262 * t260 + t267 * (pkin(4) * t337 + t367) + t290 * (pkin(4) * t289 + t330)) * MDP(26); -t481 * t510 + t411 * t509 + (-t357 * t490 + t495) * MDP(11) + ((t328 - t494) * t408 + (-t450 + t493) * t406) * MDP(12) + (-t391 * t457 + (t391 * t483 + t407 * t429) * qJD(1)) * MDP(13) + (t391 * t459 + (-t391 * t486 + t407 * t430) * qJD(1)) * MDP(14) + t391 * qJD(1) * t454 + (-pkin(2) * t450 + t468 * t391 + (pkin(7) * t490 + t492) * qJD(3) + ((-pkin(7) * t463 - t319) * t407 + (-pkin(6) * t430 - t492) * t409) * qJD(1)) * MDP(16) + (-pkin(2) * t328 - t343 * t391 + (-pkin(7) * t391 * t406 + t491) * qJD(3) + (-t373 * t483 + (-pkin(7) * t461 + t320) * t407 + (t391 * t485 + t409 * t429) * pkin(6)) * qJD(1)) * MDP(17) + (-t278 * t359 - t421 * t474) * MDP(18) + (-t278 * t420 - t279 * t359 + t312 * t474 - t421 * t473) * MDP(19) + (t398 * t279 - t309 * t420 + t506 * t312 + t473 * t329 + t431 * t441) * MDP(23) + (-t398 * t278 + t309 * t359 - t474 * t329 + t506 * t421) * MDP(24) + (-t258 * t359 + t259 * t420 + t262 * t474 - t264 * t473 + t278 * t302 - t279 * t303 - t312 * t479 - t421 * t478) * MDP(25) + (t259 * t303 + t258 * t302 + t267 * (-pkin(4) * t420 + t398) + (pkin(4) * t473 + t506) * t290 + t479 * t264 + t478 * t262) * MDP(26) + ((qJD(2) * t359 - t421) * MDP(20) + (qJD(2) * t420 + t312) * MDP(21) - t268 * MDP(23) + (-qJD(2) * t469 + t269) * MDP(24)) * t465 + (t474 * MDP(20) + t473 * MDP(21) + t505 * MDP(23) + t504 * MDP(24) + t437) * t380 + (MDP(9) * t407 * t411 + MDP(10) * t481) * pkin(1); -t357 * t356 * MDP(11) + (-t356 ^ 2 + t357 ^ 2) * MDP(12) + (t328 + t494) * MDP(13) + (-t450 - t493) * MDP(14) + qJD(1) * t438 + (-t320 * t391 - t357 * t373 + t417) * MDP(16) + (-t319 * t391 - t356 * t373 - t416) * MDP(17) + (t435 * t380 + (-t357 * t312 + t380 * t456 + t441 * t497) * pkin(3) + t501) * MDP(23) + (-t477 * t380 + (-t357 * t421 + t380 * t442 - t405 * t441) * pkin(3) + t513) * MDP(24) + (-t262 * t312 + t264 * t421 + t265 * t421 + t266 * t312 + t397 * t278 + (-t279 * t405 + (-t312 * t497 + t405 * t421) * qJD(4)) * pkin(3)) * MDP(25) + (-pkin(4) * t507 + t258 * t397 - t262 * t265 - t264 * t266 + (t259 * t405 - t290 * t357 + (-t262 * t405 + t264 * t497) * qJD(4)) * pkin(3)) * MDP(26) + t517; (-t269 * t380 + t501) * MDP(23) + (-t268 * t380 + t513) * MDP(24) + (pkin(4) * t278 - t312 * t480) * MDP(25) + (t480 * t264 + (t258 - t507) * pkin(4)) * MDP(26) + t517; (-t310 - t499) * MDP(25) + (t262 * t421 + t264 * t312 + t267) * MDP(26);];
tauc = t1;
