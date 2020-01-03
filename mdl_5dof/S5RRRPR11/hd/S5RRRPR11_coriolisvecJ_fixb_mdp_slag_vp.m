% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:23
% EndTime: 2019-12-31 21:35:34
% DurationCPUTime: 5.54s
% Computational Cost: add. (2245->428), mult. (5433->584), div. (0->0), fcn. (3419->6), ass. (0->194)
t386 = sin(qJ(3));
t389 = cos(qJ(3));
t455 = t389 * qJD(2);
t387 = sin(qJ(2));
t470 = qJD(1) * t387;
t340 = t386 * t470 - t455;
t441 = t389 * t470;
t468 = qJD(2) * t386;
t342 = t441 + t468;
t385 = sin(qJ(5));
t388 = cos(qJ(5));
t297 = t340 * t385 + t342 * t388;
t408 = -t388 * t340 + t342 * t385;
t456 = t387 * MDP(26);
t431 = qJD(2) * t456;
t524 = MDP(22) * t297 * t408 + (t297 ^ 2 - t408 ^ 2) * MDP(23) - qJD(1) * t431;
t390 = cos(qJ(2));
t469 = qJD(1) * t390;
t522 = qJD(3) - t469;
t452 = -qJD(5) + t522;
t523 = t297 * t452;
t511 = qJD(3) - qJD(5);
t359 = -qJD(2) * pkin(2) + pkin(6) * t470;
t403 = qJ(4) * t342 - t359;
t506 = pkin(3) + pkin(4);
t281 = -t506 * t340 + t403;
t363 = t522 * qJD(4);
t451 = qJD(1) * qJD(2);
t435 = t387 * t451;
t369 = qJ(4) * t435;
t353 = -pkin(2) * t390 - pkin(7) * t387 - pkin(1);
t333 = t353 * qJD(1);
t419 = pkin(2) * t387 - pkin(7) * t390;
t351 = t419 * qJD(2);
t334 = qJD(1) * t351;
t380 = pkin(6) * t469;
t360 = qJD(2) * pkin(7) + t380;
t464 = qJD(3) * t389;
t465 = qJD(3) * t386;
t404 = t333 * t464 + t386 * t334 - t360 * t465;
t423 = pkin(6) * t435;
t395 = t389 * t423 - t404;
t268 = t363 + t369 - t395;
t436 = t387 * t464;
t466 = qJD(2) * t390;
t398 = t386 * t466 + t436;
t450 = qJD(2) * qJD(3);
t314 = t398 * qJD(1) + t386 * t450;
t263 = pkin(8) * t314 + t268;
t434 = t390 * t451;
t438 = t387 * t465;
t313 = qJD(1) * t438 + (-t434 - t450) * t389;
t422 = -t333 * t465 + t389 * t334 - t360 * t464 + t386 * t423;
t264 = pkin(8) * t313 - t506 * t435 - t422;
t429 = t385 * t263 - t388 * t264;
t520 = t281 * t297 + t429;
t518 = -0.2e1 * t451;
t517 = MDP(4) * t387;
t516 = t452 * t408;
t383 = t387 ^ 2;
t515 = (-t390 ^ 2 + t383) * MDP(5);
t424 = pkin(3) * t435;
t272 = -t422 - t424;
t300 = t386 * t333 + t389 * t360;
t365 = t522 * qJ(4);
t289 = t365 + t300;
t514 = -t289 * t522 + t272;
t513 = qJD(4) * t386 + t380;
t512 = t390 * t455 - t438;
t485 = t386 * t388;
t407 = t385 * t389 - t485;
t299 = t389 * t333 - t386 * t360;
t453 = qJD(4) - t299;
t454 = pkin(8) * t342 - t453;
t274 = -t506 * t522 - t454;
t460 = qJD(5) * t388;
t445 = t388 * t263 + t385 * t264 + t274 * t460;
t510 = -t281 * t408 + t445;
t502 = qJ(4) * t386;
t509 = -t506 * t389 - t502;
t507 = t342 ^ 2;
t505 = pkin(7) - pkin(8);
t504 = pkin(6) * t386;
t503 = qJ(4) * t340;
t501 = qJ(4) * t389;
t396 = -pkin(6) * t434 - qJ(4) * t313 + qJD(4) * t342;
t270 = pkin(3) * t314 - t396;
t500 = t270 * t386;
t499 = t270 * t389;
t291 = pkin(3) * t340 - t403;
t497 = t291 * t342;
t496 = t313 * t386;
t495 = t340 * t342;
t494 = t340 * t522;
t493 = t342 * t522;
t348 = t419 * qJD(1);
t492 = t348 * t389;
t491 = t359 * t386;
t490 = t359 * t389;
t489 = t522 * t389;
t284 = pkin(8) * t340 + t300;
t279 = t284 + t365;
t488 = t385 * t279;
t486 = t386 * t387;
t484 = t386 * t390;
t483 = t387 * t389;
t392 = qJD(2) ^ 2;
t482 = t387 * t392;
t481 = t389 * t390;
t480 = t390 * t392;
t393 = qJD(1) ^ 2;
t479 = t390 * t393;
t344 = t385 * t386 + t388 * t389;
t400 = t344 * t390;
t478 = -qJD(1) * t400 + t511 * t344;
t461 = qJD(5) * t385;
t477 = t385 * t464 + t386 * t460 - t388 * t465 - t389 * t461 - t407 * t469;
t414 = pkin(3) * t386 - t501;
t476 = -t522 * t414 + t513;
t402 = -t506 * t386 + t501;
t475 = t522 * t402 + t513;
t331 = t386 * t348;
t474 = qJ(4) * t470 + t331;
t473 = t386 * t351 + t353 * t464;
t375 = pkin(6) * t481;
t472 = t386 * t353 + t375;
t467 = qJD(2) * t387;
t462 = qJD(4) * t389;
t458 = t359 * qJD(3);
t457 = t387 * MDP(15);
t449 = pkin(7) * t522 * t386;
t448 = pkin(7) * t489;
t447 = pkin(7) * t467;
t446 = pkin(7) * t455;
t362 = t505 * t389;
t444 = -t388 * t313 + t385 * t314 + t340 * t460;
t443 = qJ(4) * t467 + t473;
t442 = -pkin(3) - t504;
t437 = t390 * t465;
t432 = qJD(2) * t457;
t430 = pkin(1) * t518;
t428 = -t313 * t385 - t388 * t314;
t374 = pkin(6) * t484;
t427 = t353 * t389 - t374;
t426 = t340 + t455;
t425 = -t342 + t468;
t420 = t442 * t387;
t311 = -qJ(4) * t390 + t472;
t418 = -qJD(3) * t375 + t351 * t389 - t353 * t465;
t397 = -pkin(8) * t481 + (-pkin(4) + t442) * t387;
t417 = t397 * qJD(1) - t511 * t362 - t492;
t361 = t505 * t386;
t416 = -qJD(5) * t361 + (-pkin(6) * t483 + pkin(8) * t484) * qJD(1) + t474 + t505 * t465;
t415 = pkin(3) * t389 + t502;
t413 = qJ(4) * t388 - t385 * t506;
t412 = qJ(4) * t385 + t388 * t506;
t259 = t385 * t274 + t388 * t279;
t288 = -pkin(3) * t522 + t453;
t411 = t288 * t389 - t289 * t386;
t382 = t390 * pkin(3);
t292 = pkin(4) * t390 + t374 + t382 + (-pkin(8) * t387 - t353) * t389;
t298 = pkin(8) * t486 + t311;
t410 = t292 * t388 - t298 * t385;
t409 = t292 * t385 + t298 * t388;
t406 = qJD(1) * t383 + t390 * t522;
t405 = pkin(6) + t414;
t401 = t300 * t522 + t422;
t327 = t344 * t387;
t266 = -t342 * t461 + t444;
t399 = -pkin(6) + t402;
t267 = t297 * qJD(5) + t428;
t394 = t299 * t522 + t395;
t352 = -pkin(2) - t415;
t337 = pkin(2) - t509;
t326 = t385 * t483 - t387 * t485;
t321 = t405 * t387;
t312 = t382 - t427;
t310 = t399 * t387;
t306 = pkin(3) * t342 + t503;
t305 = qJD(1) * t420 - t492;
t304 = -pkin(6) * t441 + t474;
t286 = -t506 * t342 - t503;
t285 = -t313 + t494;
t282 = (t415 * qJD(3) - t462) * t387 + t405 * t466;
t280 = qJD(2) * t420 - t418;
t278 = -qJD(4) * t390 + (-t387 * t455 - t437) * pkin(6) + t443;
t277 = t511 * t387 * t407 + qJD(2) * t400;
t276 = qJD(5) * t327 + t512 * t385 - t398 * t388;
t275 = (t509 * qJD(3) + t462) * t387 + t399 * t466;
t271 = (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t483 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t386) * t390 + t443;
t269 = pkin(8) * t438 + t397 * qJD(2) - t418;
t265 = -t506 * t314 + t396;
t258 = t274 * t388 - t488;
t1 = [t515 * t518 + (-pkin(6) * t480 + t387 * t430) * MDP(9) + (pkin(6) * t482 + t390 * t430) * MDP(10) + (-t313 * t483 + t512 * t342) * MDP(11) + ((-t340 * t389 - t342 * t386) * t466 + (t496 - t314 * t389 + (t340 * t386 - t342 * t389) * qJD(3)) * t387) * MDP(12) + (-t522 * t438 + t313 * t390 + (t342 * t387 + t406 * t389) * qJD(2)) * MDP(13) + (-t522 * t436 + t314 * t390 + (-t340 * t387 - t406 * t386) * qJD(2)) * MDP(14) + (t522 - t469) * t432 + (t418 * t522 - t422 * t390 + (pkin(6) * t314 + t389 * t458) * t387 + ((pkin(6) * t340 + t491) * t390 + (t427 * qJD(1) + t299 + (t522 + t469) * t504) * t387) * qJD(2)) * MDP(16) + (-(-pkin(6) * t437 + t473) * t522 + t404 * t390 + (-pkin(6) * t313 - t386 * t458) * t387 + ((pkin(6) * t342 + t490) * t390 + (pkin(6) * t489 - t472 * qJD(1) - t300) * t387) * qJD(2)) * MDP(17) + (-t280 * t522 + t282 * t340 + t314 * t321 + (t291 * t468 + t272) * t390 + (t291 * t464 + t500 + (-qJD(1) * t312 - t288) * qJD(2)) * t387) * MDP(18) + (-t278 * t340 + t280 * t342 - t311 * t314 - t312 * t313 + t411 * t466 + (-t268 * t386 + t272 * t389 + (-t288 * t386 - t289 * t389) * qJD(3)) * t387) * MDP(19) + (t278 * t522 - t282 * t342 + t313 * t321 + (-t291 * t455 - t268) * t390 + (t291 * t465 - t499 + (qJD(1) * t311 + t289) * qJD(2)) * t387) * MDP(20) + (t268 * t311 + t270 * t321 + t272 * t312 + t278 * t289 + t280 * t288 + t282 * t291) * MDP(21) + (t266 * t327 + t277 * t297) * MDP(22) + (-t266 * t326 - t267 * t327 - t276 * t297 - t277 * t408) * MDP(23) + (t266 * t390 - t277 * t452 + (-qJD(1) * t327 - t297) * t467) * MDP(24) + (-t267 * t390 + t276 * t452 + (qJD(1) * t326 + t408) * t467) * MDP(25) + (t452 - t469) * t431 + (-(t269 * t388 - t271 * t385) * t452 - t429 * t390 + t275 * t408 + t310 * t267 + t265 * t326 + t281 * t276 + (-t259 * t390 + t409 * t452) * qJD(5) + (-t410 * qJD(1) - t258) * t467) * MDP(27) + ((t410 * qJD(5) + t269 * t385 + t271 * t388) * t452 - (-t279 * t461 + t445) * t390 + t275 * t297 + t310 * t266 + t265 * t327 + t281 * t277 + (t409 * qJD(1) + t259) * t467) * MDP(28) + MDP(6) * t480 - MDP(7) * t482 + 0.2e1 * t434 * t517; -t479 * t517 + t393 * t515 + (t342 * t489 - t496) * MDP(11) + ((-t313 - t494) * t389 + (-t314 - t493) * t386) * MDP(12) + t522 * t464 * MDP(13) - t522 * t465 * MDP(14) + (-t348 * t489 - pkin(2) * t314 + (-t448 + t491) * qJD(3)) * MDP(16) + (pkin(2) * t313 + t331 * t522 + (t449 + t490) * qJD(3)) * MDP(17) + (-t499 + t305 * t522 + t314 * t352 - t476 * t340 + (t291 * t386 - t448) * qJD(3)) * MDP(18) + (t304 * t340 - t305 * t342 + (t268 + t522 * t288 + (qJD(3) * t342 - t314) * pkin(7)) * t389 + ((qJD(3) * t340 - t313) * pkin(7) + t514) * t386) * MDP(19) + (-t500 - t304 * t522 + t313 * t352 + t476 * t342 + (-t291 * t389 - t449) * qJD(3)) * MDP(20) + (t270 * t352 - t288 * t305 - t289 * t304 - t476 * t291 + (t411 * qJD(3) + t268 * t389 + t272 * t386) * pkin(7)) * MDP(21) + (-t266 * t407 + t478 * t297) * MDP(22) + (-t266 * t344 + t267 * t407 - t477 * t297 - t408 * t478) * MDP(23) + (-t478 * t452 + (qJD(2) * t407 + t297) * t470) * MDP(24) + (t477 * t452 + (qJD(2) * t344 - t408) * t470) * MDP(25) + (t265 * t344 + t337 * t267 - (t416 * t385 - t417 * t388) * t452 + t475 * t408 + t477 * t281 + (-(t361 * t388 - t362 * t385) * qJD(2) + t258) * t470) * MDP(27) + (-t265 * t407 + t337 * t266 - (t417 * t385 + t416 * t388) * t452 + t475 * t297 + t478 * t281 + ((t361 * t385 + t362 * t388) * qJD(2) - t259) * t470) * MDP(28) + ((t425 * t387 - t481 * t522) * MDP(13) + (t426 * t387 + t484 * t522) * MDP(14) - t522 * t457 + (-t299 * t387 + (-t359 * t390 - t447) * t386 + (-t426 * t390 - t486 * t522) * pkin(6)) * MDP(16) + (-t359 * t481 + (t300 - t446) * t387 + (t425 * t390 - t483 * t522) * pkin(6)) * MDP(17) + (t288 * t387 + (-t291 * t390 - t447) * t386) * MDP(18) + (t291 * t481 + (-t289 + t446) * t387) * MDP(20) - t452 * t456) * qJD(1) + (t393 * t387 * MDP(9) + MDP(10) * t479) * pkin(1); MDP(11) * t495 + (-t340 ^ 2 + t507) * MDP(12) + t285 * MDP(13) + (-t314 + t493) * MDP(14) + qJD(1) * t432 + (-t342 * t359 + t401) * MDP(16) + (t340 * t359 + t394) * MDP(17) + (-t306 * t340 + t401 + 0.2e1 * t424 - t497) * MDP(18) + (pkin(3) * t313 - qJ(4) * t314 + (t289 - t300) * t342 + (t288 - t453) * t340) * MDP(19) + (-t291 * t340 + t306 * t342 + 0.2e1 * t363 + 0.2e1 * t369 - t394) * MDP(20) + (-pkin(3) * t272 + qJ(4) * t268 - t288 * t300 + t453 * t289 - t291 * t306) * MDP(21) + (-t266 + t516) * MDP(24) + (t267 + t523) * MDP(25) + (t412 * t435 - t286 * t408 - (-t388 * t284 + t454 * t385) * t452 + (t413 * t452 + t259) * qJD(5) + t520) * MDP(27) + (t413 * t435 - t286 * t297 - (t385 * t284 + t454 * t388) * t452 + (-t412 * t452 - t488) * qJD(5) + t510) * MDP(28) - t524; (-t435 + t495) * MDP(18) + t285 * MDP(19) + (-t522 ^ 2 - t507) * MDP(20) + (t497 + t514) * MDP(21) + (-t342 * t408 - t388 * t435) * MDP(27) + (-t297 * t342 + t385 * t435) * MDP(28) - (MDP(27) * t385 + MDP(28) * t388) * t452 ^ 2; (t444 - t516) * MDP(24) + (-t428 - t523) * MDP(25) + (-t259 * t452 - t520) * MDP(27) + (-t258 * t452 - t510) * MDP(28) + ((-MDP(25) * t342 - MDP(27) * t279) * t388 + (-t342 * MDP(24) - t340 * MDP(25) - MDP(27) * t274 + MDP(28) * t279) * t385) * qJD(5) + t524;];
tauc = t1;
