% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:53
% EndTime: 2019-03-09 08:20:01
% DurationCPUTime: 4.29s
% Computational Cost: add. (2251->429), mult. (5241->571), div. (0->0), fcn. (3113->6), ass. (0->192)
t506 = pkin(3) + pkin(7);
t411 = sin(qJ(2));
t478 = qJD(1) * t411;
t392 = pkin(7) * t478;
t465 = pkin(3) * t478 + qJD(3) + t392;
t463 = MDP(17) + MDP(20);
t464 = qJD(1) * qJD(2);
t514 = -0.2e1 * t464;
t405 = t411 ^ 2;
t413 = cos(qJ(2));
t513 = MDP(5) * (-t413 ^ 2 + t405);
t408 = cos(pkin(9));
t410 = sin(qJ(6));
t407 = sin(pkin(9));
t412 = cos(qJ(6));
t489 = t412 * t407;
t435 = -t408 * t410 + t489;
t334 = t435 * t413;
t511 = t408 * qJD(5) - t465;
t404 = qJD(2) * qJ(3);
t510 = qJD(4) + t404;
t384 = -qJD(6) + t478;
t509 = qJD(6) - t384;
t508 = MDP(15) + MDP(19);
t507 = MDP(16) - MDP(21);
t502 = qJ(5) * t407;
t505 = pkin(4) + pkin(5);
t430 = t408 * t505 + t502;
t477 = qJD(1) * t413;
t456 = t407 * t477;
t475 = qJD(2) * t408;
t355 = -t456 + t475;
t352 = t355 ^ 2;
t409 = -pkin(2) - qJ(4);
t504 = pkin(8) + t409;
t503 = qJD(2) * pkin(2);
t452 = -qJ(3) * t411 - pkin(1);
t350 = t409 * t413 + t452;
t328 = t350 * qJD(1);
t332 = qJD(2) * t409 + t465;
t286 = -t407 * t328 + t332 * t408;
t442 = qJD(5) - t286;
t281 = -pkin(4) * t478 + t442;
t501 = t281 * t413;
t287 = t408 * t328 + t407 * t332;
t282 = qJ(5) * t478 + t287;
t500 = t282 * t413;
t499 = t286 * t413;
t498 = t287 * t413;
t415 = qJD(1) ^ 2;
t497 = t405 * t415;
t496 = t407 * t411;
t495 = t407 * t413;
t493 = t408 * t411;
t492 = t408 * t413;
t414 = qJD(2) ^ 2;
t491 = t411 * t414;
t459 = t408 * t477;
t476 = qJD(2) * t407;
t353 = t459 + t476;
t490 = t412 * t353;
t488 = t413 * t414;
t487 = t413 * t415;
t455 = t411 * t464;
t383 = pkin(2) * t455;
t438 = -qJ(3) * t413 + qJ(4) * t411;
t472 = qJD(3) * t411;
t417 = qJD(2) * t438 - qJD(4) * t413 - t472;
t306 = qJD(1) * t417 + t383;
t454 = t413 * t464;
t382 = pkin(7) * t454;
t395 = pkin(3) * t477;
t338 = t382 + (-qJD(4) + t395) * qJD(2);
t276 = t408 * t306 + t407 * t338;
t474 = qJD(2) * t411;
t396 = pkin(2) * t474;
t318 = t396 + t417;
t473 = qJD(2) * t413;
t364 = t506 * t473;
t290 = t408 * t318 + t407 * t364;
t441 = pkin(4) * t408 + t502;
t486 = -t441 * t478 + t511;
t397 = pkin(2) * t478;
t339 = qJD(1) * t438 + t397;
t394 = pkin(7) * t477;
t363 = t394 + t395;
t301 = t408 * t339 + t407 * t363;
t376 = t506 * t411;
t312 = t408 * t350 + t407 * t376;
t447 = t407 * t455;
t467 = qJD(6) * t412;
t485 = t353 * t467 + t412 * t447;
t460 = t408 * t478;
t484 = t435 * qJD(6) + t410 * t460 - t478 * t489;
t434 = t410 * t407 + t412 * t408;
t345 = t434 * qJD(6);
t425 = t434 * t411;
t483 = -qJD(1) * t425 + t345;
t482 = -t430 * t478 + t511;
t377 = t506 * t413;
t403 = qJD(2) * qJD(3);
t480 = -pkin(7) * t455 + t403;
t372 = -pkin(2) * t413 + t452;
t349 = qJD(1) * t372;
t471 = qJD(4) * t353;
t470 = qJD(4) * t355;
t469 = qJD(5) * t355;
t269 = pkin(8) * t353 + t282;
t468 = qJD(6) * t269;
t399 = t411 * qJD(5);
t462 = t505 * t411;
t461 = t411 * t487;
t293 = qJ(5) * t477 + t301;
t304 = t411 * qJ(5) + t312;
t458 = t409 * t473;
t457 = qJD(5) * t495;
t453 = MDP(27) * t477;
t451 = pkin(1) * t514;
t450 = qJD(3) - t503;
t449 = qJ(5) * t408 - qJ(3);
t275 = -t407 * t306 + t338 * t408;
t289 = -t407 * t318 + t364 * t408;
t300 = -t407 * t339 + t363 * t408;
t311 = -t407 * t350 + t376 * t408;
t268 = qJ(5) * t454 + qJD(1) * t399 + t276;
t446 = t408 * t455;
t266 = -pkin(8) * t446 + t268;
t267 = -pkin(8) * t355 - qJD(1) * t462 + t442;
t448 = -qJD(6) * t267 - t266;
t279 = qJ(5) * t473 + t290 + t399;
t343 = t363 + t510;
t270 = -pkin(4) * t454 - t275;
t445 = t282 * t478 - t270;
t444 = t287 * t478 + t275;
t440 = -t343 + t510;
t424 = qJ(5) * t355 - t343;
t291 = pkin(4) * t353 - t424;
t369 = pkin(4) * t407 - t449;
t439 = -qJD(2) * t369 - qJD(4) + t291;
t262 = t267 * t412 - t269 * t410;
t263 = t267 * t410 + t269 * t412;
t288 = pkin(8) * t495 - t311 - t462;
t294 = pkin(8) * t492 + t304;
t437 = t288 * t412 - t294 * t410;
t436 = t288 * t410 + t294 * t412;
t309 = t353 * t410 + t355 * t412;
t433 = -0.2e1 * qJD(2) * t349;
t432 = -pkin(3) - t441;
t428 = -qJ(3) * t473 - t472;
t329 = qJD(1) * t428 + t383;
t341 = t396 + t428;
t431 = pkin(7) * t414 + qJD(1) * t341 + t329;
t366 = t504 * t407;
t423 = -pkin(8) * t496 - t413 * t505;
t427 = -qJD(1) * t423 + qJD(4) * t408 - qJD(6) * t366 + t300;
t367 = t504 * t408;
t426 = pkin(8) * t460 - qJD(4) * t407 - qJD(6) * t367 - t293;
t337 = -pkin(3) * t455 + t480;
t422 = pkin(3) + t430;
t421 = -qJD(6) * t355 - t446;
t420 = qJD(2) * t425;
t419 = t423 * qJD(2);
t418 = t432 * t455 + t480;
t371 = t392 + t450;
t375 = -t394 - t404;
t416 = t480 * t413 + (t371 * t413 + (t375 + t394) * t411) * qJD(2);
t365 = t408 * t409 * t454;
t362 = t506 * t474;
t360 = -qJ(3) * t477 + t397;
t346 = -t407 * t505 + t449;
t333 = t434 * t413;
t331 = t349 * t478;
t321 = t413 * t441 + t377;
t315 = -t413 * t430 - t377;
t307 = t355 * t410 - t490;
t305 = -pkin(4) * t411 - t311;
t303 = t457 + (-pkin(7) + t432) * t474;
t297 = t345 * t413 + t435 * t474;
t296 = -qJD(6) * t334 + t420;
t295 = -pkin(4) * t477 - t300;
t292 = -t457 + (pkin(7) + t422) * t474;
t284 = -pkin(4) * t473 - t289;
t283 = t418 - t469;
t278 = qJD(1) * t420 + qJD(6) * t309;
t277 = t410 * t421 + t485;
t274 = -t353 * t505 + t424;
t273 = t422 * t455 + t469 - t480;
t272 = -pkin(8) * t408 * t474 + t279;
t271 = -t289 + t419;
t265 = qJD(1) * t419 - t275;
t264 = t412 * t265;
t1 = [(-t411 * t431 + t413 * t433) * MDP(13) + (t411 * t433 + t413 * t431) * MDP(12) + (pkin(7) * t416 + t329 * t372 + t341 * t349) * MDP(14) + t416 * MDP(11) + (t384 + t478) * MDP(27) * t473 + (t275 * t311 + t276 * t312 + t286 * t289 + t287 * t290 + t337 * t377 - t343 * t362) * MDP(18) + (-t337 * t495 - t355 * t362 + (-qJD(1) * t290 - t276) * t411 + (t343 * t496 - t498 + (-t312 * t413 + t377 * t496) * qJD(1)) * qJD(2)) * MDP(16) + (t337 * t492 - t353 * t362 + (qJD(1) * t289 + t275) * t411 + (-t343 * t493 + t499 + (t311 * t413 - t377 * t493) * qJD(1)) * qJD(2)) * MDP(15) + 0.2e1 * t411 * MDP(4) * t454 + (-t277 * t334 + t297 * t309) * MDP(23) + (t277 * t333 + t278 * t334 - t296 * t309 - t297 * t307) * MDP(24) + (-t277 * t411 - t297 * t384 + (qJD(1) * t334 - t309) * t473) * MDP(25) + ((t271 * t410 + t272 * t412) * t384 + (t265 * t410 + t266 * t412) * t411 + t292 * t309 + t315 * t277 - t273 * t334 + t274 * t297 + (t262 * t411 + t384 * t437) * qJD(6) + (qJD(1) * t436 + t263) * t473) * MDP(29) + (t268 * t304 + t270 * t305 + t279 * t282 + t281 * t284 + t283 * t321 + t291 * t303) * MDP(22) + MDP(6) * t488 + t513 * t514 + (t278 * t411 + t296 * t384 + (-qJD(1) * t333 + t307) * t473) * MDP(26) + (-(t271 * t412 - t272 * t410) * t384 - (-t266 * t410 + t264) * t411 + t292 * t307 + t315 * t278 - t273 * t333 + t274 * t296 + (t263 * t411 + t384 * t436) * qJD(6) + (-qJD(1) * t437 - t262) * t473) * MDP(28) + (-t289 * t355 - t290 * t353 + (t275 * t407 - t276 * t408) * t413 + (-t286 * t407 + t287 * t408 + (-t311 * t407 + t312 * t408) * qJD(1)) * t474) * MDP(17) + (-t279 * t353 + t284 * t355 + (-t268 * t408 - t270 * t407) * t413 + (t281 * t407 + t282 * t408 + (t304 * t408 + t305 * t407) * qJD(1)) * t474) * MDP(20) + (-pkin(7) * t488 + t411 * t451) * MDP(9) - MDP(7) * t491 + (pkin(7) * t491 + t413 * t451) * MDP(10) + (t283 * t495 - t303 * t355 + (qJD(1) * t279 + t268) * t411 + (-t291 * t496 + t500 + (t304 * t413 - t321 * t496) * qJD(1)) * qJD(2)) * MDP(21) + (t283 * t492 + t303 * t353 + (-qJD(1) * t284 - t270) * t411 + (-t291 * t493 - t501 + (-t305 * t413 - t321 * t493) * qJD(1)) * qJD(2)) * MDP(19); -MDP(4) * t461 + t415 * t513 + ((-t375 - t404) * t411 + (-t371 + t450) * t413) * qJD(1) * MDP(11) + (-t360 * t477 + t331) * MDP(12) + (0.2e1 * t403 + (t349 * t413 + t360 * t411) * qJD(1)) * MDP(13) + (qJ(3) * t480 - qJD(3) * t375 - t349 * t360 + (-t375 * t411 + (-t371 - t503) * t413) * qJD(1) * pkin(7)) * MDP(14) + (t337 * t407 + t365 + t465 * t353 + (-t499 + (-t408 * t440 - t300) * t411) * qJD(1)) * MDP(15) + (t337 * t408 + t465 * t355 + (t498 + t301 * t411 + (t411 * t440 - t458) * t407) * qJD(1)) * MDP(16) + (t300 * t355 + t301 * t353 + (-t444 + t470) * t408 + (t286 * t478 - t276 + t471) * t407) * MDP(17) + (qJ(3) * t337 - t286 * t300 - t287 * t301 + (t275 * t408 + t276 * t407) * t409 + t465 * t343 + (-t286 * t408 - t287 * t407) * qJD(4)) * MDP(18) + (t283 * t407 + t365 - t486 * t353 + (t501 + (t408 * t439 + t295) * t411) * qJD(1)) * MDP(19) + (t293 * t353 - t295 * t355 + (-t445 + t470) * t408 + (-t281 * t478 - t268 + t471) * t407) * MDP(20) + (-t283 * t408 + t486 * t355 + (-t500 - t293 * t411 + (t411 * t439 + t458) * t407) * qJD(1)) * MDP(21) + (-t281 * t295 - t282 * t293 + t283 * t369 + (t268 * t407 - t270 * t408) * t409 - t486 * t291 + (t281 * t408 - t282 * t407) * qJD(4)) * MDP(22) + (t277 * t434 + t309 * t484) * MDP(23) + (t277 * t435 - t278 * t434 - t307 * t484 - t309 * t483) * MDP(24) + (-t484 * t384 + (-qJD(2) * t434 + t309) * t477) * MDP(25) + (t483 * t384 + (-qJD(2) * t435 - t307) * t477) * MDP(26) - t384 * t453 + (-t273 * t435 + t346 * t278 + (t410 * t426 - t412 * t427) * t384 + t482 * t307 + t483 * t274 + (-(-t366 * t410 - t367 * t412) * qJD(2) + t262) * t477) * MDP(28) + (t273 * t434 + t346 * t277 + (t410 * t427 + t412 * t426) * t384 + t482 * t309 + t484 * t274 + ((t366 * t412 - t367 * t410) * qJD(2) - t263) * t477) * MDP(29) + (MDP(9) * t411 * t415 + MDP(10) * t487) * pkin(1); MDP(12) * t461 + (-t414 - t497) * MDP(13) + (t331 + t382) * MDP(14) - t463 * t353 * t460 + t508 * (-t407 * t497 + (-t353 + t459) * qJD(2)) + (t276 * MDP(18) + t268 * MDP(22) + (-t286 * MDP(18) + MDP(22) * t281 + t355 * t463) * t478) * t407 + (t444 * MDP(18) + t445 * MDP(22) - t497 * t507) * t408 + (MDP(28) * t484 - MDP(29) * t483) * t384 + (t375 * MDP(14) - t343 * MDP(18) - t291 * MDP(22) + (t434 * t477 + t307) * MDP(28) + (t435 * t477 + t309) * MDP(29) - t507 * (t355 + t456)) * qJD(2); (t286 * t355 + t287 * t353 + t337) * MDP(18) + (t282 * t353 + (-qJD(5) - t281) * t355 + t418) * MDP(22) + (-t309 * t509 - t434 * t455) * MDP(28) + (t384 * t490 + (t355 * t509 + t446) * t410 - t485) * MDP(29) + t463 * (-t353 ^ 2 - t352) + (t508 * (t355 - t475) + t507 * (-t353 + t476)) * t478; t355 * t353 * MDP(19) + (-t352 - t497) * MDP(21) + (t291 * t355 - t275) * MDP(22) + (qJD(6) * t384 * t410 - t307 * t355) * MDP(28) + (-t309 * t355 + t384 * t467) * MDP(29) + ((-MDP(22) * pkin(4) - t412 * MDP(28) + t410 * MDP(29) - MDP(19)) * t473 + ((t353 + t476) * MDP(20) - t282 * MDP(22) + (-t410 * MDP(28) - t412 * MDP(29)) * t384) * t411) * qJD(1); -t307 ^ 2 * MDP(24) + (-t307 * t384 + t485) * MDP(25) - qJD(2) * t453 + (-t263 * t384 + t264) * MDP(28) + (-t262 * t384 + t274 * t307) * MDP(29) + (MDP(23) * t307 + MDP(24) * t309 - MDP(26) * t384 - MDP(28) * t274) * t309 + (MDP(26) * t421 - MDP(28) * t468 + MDP(29) * t448) * t412 + (t421 * MDP(25) + (-qJD(6) * t353 - t447) * MDP(26) + t448 * MDP(28) + (-t265 + t468) * MDP(29)) * t410;];
tauc  = t1;
