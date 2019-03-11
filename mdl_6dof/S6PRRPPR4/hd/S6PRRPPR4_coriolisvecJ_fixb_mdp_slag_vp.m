% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:35
% EndTime: 2019-03-08 21:16:44
% DurationCPUTime: 4.65s
% Computational Cost: add. (2325->440), mult. (6059->600), div. (0->0), fcn. (4388->10), ass. (0->200)
t531 = MDP(14) + MDP(17);
t417 = sin(qJ(3));
t538 = MDP(5) * t417;
t420 = cos(qJ(3));
t411 = t420 ^ 2;
t537 = MDP(6) * (t417 ^ 2 - t411);
t445 = pkin(3) * t417 - qJ(4) * t420;
t377 = t445 * qJD(2);
t412 = sin(pkin(11));
t414 = cos(pkin(11));
t418 = sin(qJ(2));
t413 = sin(pkin(6));
t496 = qJD(1) * t413;
t471 = t418 * t496;
t382 = qJD(2) * pkin(8) + t471;
t415 = cos(pkin(6));
t495 = qJD(1) * t415;
t534 = -t417 * t382 + t420 * t495;
t305 = t412 * t377 + t414 * t534;
t493 = qJD(2) * t417;
t296 = qJ(5) * t493 + t305;
t487 = qJD(4) * t414;
t536 = -t296 + t487;
t357 = qJD(3) * t445 - qJD(4) * t417;
t421 = cos(qJ(2));
t470 = t421 * t496;
t514 = t412 * t420;
t535 = t470 * t514 + (t357 - t471) * t414;
t532 = MDP(12) + MDP(16);
t479 = MDP(13) - MDP(18);
t509 = t414 * t421;
t342 = (t412 * t418 + t420 * t509) * t496;
t484 = qJD(5) * t420;
t490 = qJD(3) * t417;
t530 = qJ(5) * t490 - t342 - t484;
t481 = t414 * qJD(3);
t371 = t412 * t493 - t481;
t491 = qJD(3) * t412;
t373 = t414 * t493 + t491;
t416 = sin(qJ(6));
t419 = cos(qJ(6));
t318 = t371 * t416 + t373 * t419;
t480 = qJD(2) * qJD(3);
t464 = t420 * t480;
t394 = t412 * t464;
t449 = t414 * t464;
t428 = -t419 * t394 + t416 * t449;
t289 = qJD(6) * t318 + t428;
t369 = t373 ^ 2;
t529 = pkin(4) + pkin(5);
t528 = pkin(9) * t417;
t527 = -pkin(9) + qJ(4);
t526 = qJD(2) * pkin(2);
t525 = qJ(4) * t417;
t524 = qJ(5) * t414;
t397 = t417 * t495;
t494 = qJD(2) * t413;
t467 = t421 * t494;
t450 = qJD(1) * t467;
t489 = qJD(3) * t420;
t313 = qJD(3) * t397 + t382 * t489 + t417 * t450;
t454 = -pkin(4) * t394 - t313;
t424 = -qJ(5) * t449 - t454;
t486 = qJD(5) * t373;
t283 = t424 - t486;
t523 = t283 * t412;
t522 = t283 * t414;
t521 = t313 * t412;
t520 = t313 * t414;
t516 = t373 * t416;
t316 = -t419 * t371 + t516;
t492 = qJD(2) * t420;
t403 = qJD(6) + t492;
t519 = t316 * t403;
t518 = t318 * t403;
t515 = t412 * t419;
t513 = t413 * t418;
t512 = t413 * t421;
t511 = t414 * t417;
t510 = t414 * t420;
t422 = qJD(3) ^ 2;
t508 = t417 * t422;
t507 = t420 * t422;
t473 = -pkin(8) * t412 - pkin(4);
t478 = pkin(9) * t510;
t506 = (-t478 + (-pkin(5) + t473) * t417) * qJD(3) - t535;
t350 = t412 * t357;
t505 = -t350 - (-pkin(8) * t511 + pkin(9) * t514) * qJD(3) - t530;
t309 = t420 * t450 + (qJD(4) + t534) * qJD(3);
t322 = (t357 + t471) * qJD(2);
t281 = t414 * t309 + t412 * t322;
t477 = pkin(8) * t490;
t455 = t414 * t477;
t326 = t350 - t455;
t503 = t326 + t530;
t502 = t473 * t490 - t535;
t344 = t420 * t382 + t397;
t338 = qJD(3) * qJ(4) + t344;
t446 = -pkin(3) * t420 - t525;
t387 = -pkin(2) + t446;
t345 = qJD(2) * t387 - t470;
t294 = t414 * t338 + t412 * t345;
t456 = t412 * t477;
t501 = t456 + t535;
t500 = t326 - t342;
t375 = t412 * t416 + t414 * t419;
t431 = t375 * t420;
t499 = -qJD(2) * t431 - t375 * qJD(6);
t469 = t412 * t492;
t476 = t416 * t510;
t482 = qJD(6) * t419;
t483 = qJD(6) * t416;
t498 = -qJD(2) * t476 + t412 * t482 - t414 * t483 + t419 * t469;
t347 = pkin(8) * t510 + t412 * t387;
t488 = qJD(4) * t373;
t485 = qJD(5) * t412;
t465 = t417 * t480;
t475 = qJ(5) * t465 + t281;
t474 = t371 * t482 + t416 * t394 + t419 * t449;
t472 = qJ(4) * t490;
t468 = t418 * t494;
t466 = qJD(5) * t511;
t463 = MDP(24) * t493;
t462 = qJ(5) * t412 + pkin(3);
t461 = qJD(3) * pkin(3) - qJD(4);
t280 = -t412 * t309 + t322 * t414;
t426 = (-t417 * t529 - t478) * qJD(2);
t274 = qJD(3) * t426 - t280;
t275 = (pkin(9) * t491 - qJD(5)) * t492 + t475;
t460 = t419 * t274 - t275 * t416;
t293 = -t412 * t338 + t345 * t414;
t304 = t377 * t414 - t412 * t534;
t401 = pkin(8) * t514;
t346 = t387 * t414 - t401;
t436 = -t412 * t529 + t524;
t458 = -t436 * t492 + t344 + t485;
t444 = pkin(4) * t412 - t524;
t311 = t444 * t492 + t344;
t457 = t311 + t485;
t453 = t371 * t470;
t452 = t373 * t470;
t451 = t417 * t470;
t337 = -qJ(5) * t420 + t347;
t383 = -t470 - t526;
t448 = -t383 - t470;
t443 = t274 * t416 + t275 * t419;
t285 = pkin(4) * t492 + qJD(5) - t293;
t278 = pkin(5) * t492 - pkin(9) * t373 + t285;
t286 = -qJ(5) * t492 + t294;
t282 = pkin(9) * t371 + t286;
t271 = t278 * t419 - t282 * t416;
t272 = t278 * t416 + t282 * t419;
t409 = t420 * pkin(4);
t314 = pkin(5) * t420 + t401 + t409 + (-t387 - t528) * t414;
t320 = t412 * t528 + t337;
t442 = t314 * t419 - t320 * t416;
t441 = t314 * t416 + t320 * t419;
t376 = -t414 * t416 + t515;
t439 = -t416 * MDP(25) - t419 * MDP(26);
t438 = pkin(8) + t444;
t437 = t461 + t534;
t363 = t415 * t417 + t420 * t513;
t362 = -t415 * t420 + t417 * t513;
t391 = t527 * t414;
t434 = qJD(4) * t412 - qJD(6) * t391 + t304 - t426;
t390 = t527 * t412;
t433 = pkin(9) * t469 - qJD(6) * t390 - t536;
t432 = -t419 * MDP(25) + t416 * MDP(26) - MDP(16);
t353 = t375 * t417;
t430 = t373 * t483 - t474;
t429 = -pkin(8) + t436;
t427 = qJ(5) * t373 + t437;
t425 = qJD(3) * (-t448 - t526);
t423 = qJD(2) ^ 2;
t393 = qJD(4) * t469;
t384 = -pkin(4) * t414 - t462;
t365 = t414 * t529 + t462;
t355 = t371 * t492;
t354 = t371 * t487;
t352 = t416 * t511 - t417 * t515;
t351 = t438 * t417;
t339 = -t346 + t409;
t333 = t429 * t417;
t330 = qJD(3) * t363 + t417 * t467;
t329 = -qJD(3) * t362 + t420 * t467;
t328 = t363 * t414 - t412 * t512;
t327 = t363 * t412 + t413 * t509;
t321 = t438 * t489 - t466;
t312 = t429 * t489 + t466;
t308 = qJD(3) * t476 + qJD(6) * t353 - t489 * t515;
t307 = qJD(6) * t376 * t417 + qJD(3) * t431;
t303 = t329 * t414 + t412 * t468;
t302 = t329 * t412 - t414 * t468;
t299 = -pkin(4) * t493 - t304;
t291 = pkin(4) * t371 - t427;
t284 = -t371 * t529 + t427;
t279 = t486 + (-pkin(5) * t412 + t524) * t464 + t454;
t277 = -pkin(4) * t465 - t280;
t276 = -qJD(2) * t484 + t475;
t1 = [-t329 * qJD(3) * MDP(11) + (-t280 * t327 + t281 * t328 - t293 * t302 + t294 * t303 + t313 * t362) * MDP(15) + (t276 * t328 + t277 * t327 + t283 * t362 + t285 * t302 + t286 * t303) * MDP(19) + ((t302 * t419 - t303 * t416 + (-t327 * t416 - t328 * t419) * qJD(6)) * t403 - t362 * t289) * MDP(25) + (-(t302 * t416 + t303 * t419 + (t327 * t419 - t328 * t416) * qJD(6)) * t403 + t362 * t430) * MDP(26) + (-MDP(4) * t421 + (-MDP(10) * t420 + MDP(11) * t417 - MDP(3)) * t418) * t423 * t413 + (-qJD(3) * MDP(10) - MDP(15) * t437 + t291 * MDP(19) - t316 * MDP(25) - t318 * MDP(26) + t373 * t479) * t330 + ((t532 * t302 + t479 * t303) * t420 + ((t479 * t414 * t362 - MDP(11) * t512 + t531 * (t327 * t414 - t328 * t412)) * t420 + (-MDP(10) * t512 + (-t439 - t479) * t328 + (-MDP(12) + t432) * t327) * t417) * qJD(3)) * qJD(2) + t532 * (t330 * t371 + t362 * t394) + t531 * (t302 * t373 - t303 * t371); 0.2e1 * t464 * t538 - 0.2e1 * t480 * t537 + MDP(7) * t507 - MDP(8) * t508 + (-pkin(8) * t507 + t417 * t425) * MDP(10) + (pkin(8) * t508 + t420 * t425) * MDP(11) + ((-t453 + t521 + (qJD(2) * t346 + t293) * qJD(3)) * t417 + (-t280 + (pkin(8) * t371 - t412 * t437) * qJD(3) + (t456 - t501) * qJD(2)) * t420) * MDP(12) + ((-t452 + t520 + (-qJD(2) * t347 - t294) * qJD(3)) * t417 + (t281 + (pkin(8) * t373 - t414 * t437) * qJD(3) + (t455 + t500) * qJD(2)) * t420) * MDP(13) + ((-t280 * t414 - t281 * t412) * t417 - t501 * t373 - t500 * t371 + (-t293 * t414 - t294 * t412 + (-t346 * t414 - t347 * t412) * qJD(2)) * t489) * MDP(14) + (t437 * t451 + t280 * t346 + t281 * t347 + t500 * t294 + t501 * t293 + (t313 * t417 - t437 * t489) * pkin(8)) * MDP(15) + (t321 * t371 + (-t453 + t523 + (-qJD(2) * t339 - t285) * qJD(3)) * t417 + (t291 * t491 + t277 + (t351 * t491 + t502) * qJD(2)) * t420) * MDP(16) + ((-t276 * t412 + t277 * t414) * t417 + t502 * t373 - t503 * t371 + (t285 * t414 - t286 * t412 + (-t337 * t412 + t339 * t414) * qJD(2)) * t489) * MDP(17) + (-t321 * t373 + (t452 - t522 + (qJD(2) * t337 + t286) * qJD(3)) * t417 + (-t291 * t481 - t276 + (-t351 * t481 - t503) * qJD(2)) * t420) * MDP(18) + (t276 * t337 + t277 * t339 + t283 * t351 + (t321 - t451) * t291 + t503 * t286 + t502 * t285) * MDP(19) + (t307 * t318 - t353 * t430) * MDP(20) + (-t289 * t353 - t307 * t316 - t308 * t318 + t352 * t430) * MDP(21) + (-t430 * t420 + t307 * t403 + (-qJD(2) * t353 - t318) * t490) * MDP(22) + (-t289 * t420 - t308 * t403 + (qJD(2) * t352 + t316) * t490) * MDP(23) + (-t403 - t492) * MDP(24) * t490 + (t460 * t420 + t312 * t316 + t333 * t289 + t279 * t352 + t284 * t308 + (t416 * t505 + t419 * t506) * t403 + (-t272 * t420 - t403 * t441) * qJD(6) + (t316 * t470 + (-qJD(2) * t442 - t271) * qJD(3)) * t417) * MDP(25) + (-t443 * t420 + t312 * t318 - t333 * t430 + t279 * t353 + t284 * t307 + (-t416 * t506 + t419 * t505) * t403 + (-t271 * t420 - t403 * t442) * qJD(6) + (t318 * t470 + (qJD(2) * t441 + t272) * qJD(3)) * t417) * MDP(26); (qJD(3) * t344 - t383 * t493 - t313) * MDP(10) + t448 * t492 * MDP(11) + (-t520 - t344 * t371 + t393 + (-t293 * t417 + t304 * t420 + (qJD(3) * t446 + t420 * t437) * t412) * qJD(2)) * MDP(12) + (t521 - t344 * t373 + (t294 * t417 - t305 * t420 + (-t472 + (t437 - t461) * t420) * t414) * qJD(2)) * MDP(13) + (t304 * t373 + t305 * t371 - t354 + (t293 * t492 + t281) * t414 + (t294 * t492 - t280 + t488) * t412) * MDP(14) + (-pkin(3) * t313 - t293 * t304 - t294 * t305 + t437 * t344 + (-t293 * t412 + t294 * t414) * qJD(4) + (-t280 * t412 + t281 * t414) * qJ(4)) * MDP(15) + (-t522 + t393 - t457 * t371 + (t285 * t417 - t299 * t420 + (-t291 * t420 + (t384 * t420 - t525) * qJD(3)) * t412) * qJD(2)) * MDP(16) + (t296 * t371 - t299 * t373 - t354 + (-t285 * t492 + t276) * t414 + (t286 * t492 + t277 + t488) * t412) * MDP(17) + (-t523 + t457 * t373 + (-t286 * t417 + t296 * t420 + (t472 + (-qJD(3) * t384 - qJD(4) + t291) * t420) * t414) * qJD(2)) * MDP(18) + (qJ(4) * t276 * t414 + t283 * t384 - t285 * t299 - t291 * t311 + t536 * t286 + (qJ(4) * t277 + qJD(4) * t285 - qJD(5) * t291) * t412) * MDP(19) + (t318 * t499 - t376 * t430) * MDP(20) + (-t289 * t376 - t316 * t499 - t318 * t498 + t375 * t430) * MDP(21) + (t499 * t403 + (-qJD(3) * t376 + t318) * t493) * MDP(22) + (-t498 * t403 + (qJD(3) * t375 - t316) * t493) * MDP(23) + t403 * t463 + (t279 * t375 + t365 * t289 + (t416 * t433 + t419 * t434) * t403 + t458 * t316 + t498 * t284 + (-(t390 * t419 - t391 * t416) * qJD(3) + t271) * t493) * MDP(25) + (t279 * t376 - t365 * t430 + (-t416 * t434 + t419 * t433) * t403 + t458 * t318 + t499 * t284 + ((t390 * t416 + t391 * t419) * qJD(3) - t272) * t493) * MDP(26) + (-t420 * t538 + t537) * t423; (t293 * t373 + t294 * t371 + t313) * MDP(15) + (t286 * t371 + (-qJD(5) - t285) * t373 + t424) * MDP(19) + (-t289 - t518) * MDP(25) + (t430 + t519) * MDP(26) + t532 * (-t373 * t492 + t394) + t479 * (t355 + t449) + t531 * (-t371 ^ 2 - t369); t373 * t371 * MDP(16) - t355 * MDP(17) + (-t411 * t423 - t369) * MDP(18) + (t291 * t373 - t280) * MDP(19) + (-t316 * t373 - t403 * t483) * MDP(25) + (-t318 * t373 - t403 * t482) * MDP(26) + ((t286 * MDP(19) + t403 * t439) * t420 + (MDP(17) * t510 + (-MDP(19) * pkin(4) + t432) * t417) * qJD(3)) * qJD(2); t318 * t316 * MDP(20) + (-t316 ^ 2 + t318 ^ 2) * MDP(21) + (t474 + t519) * MDP(22) + (-t428 + t518) * MDP(23) - qJD(3) * t463 + (t272 * t403 - t284 * t318 + t460) * MDP(25) + (t271 * t403 + t284 * t316 - t443) * MDP(26) + (-MDP(22) * t516 - MDP(23) * t318 - MDP(25) * t272 - MDP(26) * t271) * qJD(6);];
tauc  = t1;
