% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:51
% EndTime: 2019-03-09 04:45:02
% DurationCPUTime: 6.51s
% Computational Cost: add. (5011->443), mult. (12733->528), div. (0->0), fcn. (9252->6), ass. (0->169)
t413 = sin(pkin(9));
t414 = cos(pkin(9));
t416 = sin(qJ(3));
t418 = cos(qJ(3));
t524 = -t413 * t416 + t418 * t414;
t533 = t524 * qJD(1);
t503 = pkin(7) + qJ(2);
t400 = t503 * t413;
t396 = qJD(1) * t400;
t401 = t503 * t414;
t397 = qJD(1) * t401;
t346 = -t416 * t396 + t418 * t397;
t532 = qJD(3) * t346;
t531 = MDP(23) - MDP(28);
t395 = t413 * t418 + t414 * t416;
t387 = t395 * qJD(1);
t417 = cos(qJ(4));
t470 = qJD(4) * t417;
t415 = sin(qJ(4));
t472 = qJD(3) * t415;
t318 = t387 * t470 + (qJD(4) + t533) * t472;
t389 = t395 * qJD(3);
t372 = qJD(1) * t389;
t459 = MDP(22) + MDP(26);
t368 = t372 * pkin(4);
t429 = t524 * qJD(2);
t512 = -t418 * t396 - t416 * t397;
t308 = qJD(1) * t429 + qJD(3) * t512;
t408 = -pkin(2) * t414 - pkin(1);
t398 = t408 * qJD(1) + qJD(2);
t323 = -pkin(3) * t533 - pkin(8) * t387 + t398;
t388 = t524 * qJD(3);
t425 = qJD(1) * t388;
t328 = t372 * pkin(3) - pkin(8) * t425;
t341 = qJD(3) * pkin(8) + t346;
t471 = qJD(4) * t415;
t446 = t415 * t308 + t323 * t471 - t417 * t328 + t341 * t470;
t276 = -t368 + t446;
t298 = t415 * t323 + t417 * t341;
t378 = qJD(4) - t533;
t371 = t378 * qJ(5);
t291 = t371 + t298;
t517 = -t291 * t378 + t276;
t530 = -t517 * MDP(25) + t459 * t372;
t529 = (t413 ^ 2 + t414 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t358 = t387 * t417 + t472;
t466 = t417 * qJD(3);
t317 = -qJD(4) * t466 + t387 * t471 - t417 * t425;
t426 = -qJ(6) * t317 - t276;
t422 = -pkin(5) * t372 - t426;
t270 = -qJD(6) * t358 + t422;
t356 = t387 * t415 - t466;
t288 = qJ(6) * t356 + t298;
t283 = t288 + t371;
t496 = t283 * t378;
t528 = -t270 + t496;
t432 = -t417 * t308 - t323 * t470 - t415 * t328 + t341 * t471;
t525 = t372 * qJ(5) + t378 * qJD(5);
t275 = -t432 + t525;
t516 = t318 * qJ(6) + t356 * qJD(6);
t271 = t275 + t516;
t297 = t417 * t323 - t415 * t341;
t287 = qJ(6) * t358 + t297;
t465 = qJD(5) - t287;
t505 = pkin(4) + pkin(5);
t281 = -t505 * t378 + t465;
t527 = t378 * t281 + t271;
t464 = qJD(5) - t297;
t290 = -pkin(4) * t378 + t464;
t526 = t378 * t290 + t275;
t522 = MDP(24) + MDP(27);
t515 = -qJD(5) * t415 - t346;
t514 = 0.2e1 * t525;
t513 = t389 * qJ(5) - qJD(5) * t524;
t351 = t400 * t418 + t416 * t401;
t352 = -t400 * t416 + t401 * t418;
t510 = t352 * qJD(3);
t445 = qJD(3) * pkin(3) + t512;
t427 = qJ(5) * t358 + t445;
t299 = pkin(4) * t356 - t427;
t504 = pkin(8) * t372;
t508 = t378 * t299 - t504;
t498 = qJ(5) * t415;
t507 = -t505 * t417 - t498;
t506 = t356 ^ 2;
t355 = t358 ^ 2;
t502 = pkin(8) - qJ(6);
t500 = qJ(5) * t318;
t499 = qJ(5) * t356;
t497 = qJ(5) * t417;
t494 = t317 * t415;
t493 = t356 * t378;
t492 = t356 * t533;
t491 = t358 * t378;
t448 = t378 * t417;
t490 = t533 * t415;
t489 = t395 * t415;
t488 = t395 * t417;
t363 = t415 * t372;
t365 = t417 * t372;
t338 = t415 * t512;
t342 = pkin(3) * t387 - pkin(8) * t533;
t403 = t502 * t417;
t482 = t338 + (-qJ(6) * t533 - t342) * t417 - t505 * t387 - qJD(4) * t403 + qJD(6) * t415;
t476 = t415 * t342 + t417 * t512;
t293 = t387 * qJ(5) + t476;
t481 = -qJ(6) * t490 - qJD(6) * t417 - t502 * t471 - t293;
t437 = -t505 * t415 + t497;
t480 = -t378 * t437 + t515;
t443 = pkin(4) * t415 - t497;
t479 = -t378 * t443 - t515;
t478 = -t415 * t318 - t356 * t470;
t477 = t378 * t490 + t365;
t344 = -pkin(3) * t524 - pkin(8) * t395 + t408;
t475 = t415 * t344 + t417 * t352;
t474 = t378 * t470 + t363;
t468 = qJD(5) * t417;
t467 = t299 * MDP(25);
t286 = -t505 * t356 + qJD(6) + t427;
t463 = qJD(6) + t286;
t461 = qJD(1) * qJD(2);
t458 = pkin(8) * t471;
t324 = -t351 * qJD(3) + t429;
t457 = t415 * t324 + t344 * t471 + t352 * t470;
t343 = pkin(3) * t389 - pkin(8) * t388;
t456 = t417 * t324 + t415 * t343 + t344 * t470;
t300 = -qJ(5) * t524 + t475;
t455 = t395 * t471;
t454 = t395 * t470;
t453 = t378 * t471;
t451 = MDP(21) - t522;
t348 = t415 * t352;
t450 = t344 * t417 - t348;
t449 = t378 * t415;
t447 = qJD(4) * t356 - t317;
t309 = t395 * t461 + t532;
t444 = pkin(4) * t417 + t498;
t442 = t290 * t417 - t291 * t415;
t440 = -qJ(6) * t388 - qJD(6) * t395;
t438 = t343 * t417 - t457;
t436 = t388 * t415 + t454;
t435 = -t388 * t417 + t455;
t434 = -t378 * t445 - t504;
t433 = t298 * t378 - t446;
t431 = -t352 * t471 + t456;
t278 = t318 * pkin(4) + t317 * qJ(5) - t358 * qJD(5) + t309;
t424 = t317 - t493;
t423 = t297 * t378 + t432;
t274 = -pkin(5) * t318 - t278;
t325 = t395 * qJD(2) + t510;
t402 = t502 * t415;
t399 = -pkin(3) - t444;
t392 = pkin(3) - t507;
t373 = t378 ^ 2;
t310 = pkin(4) * t358 + t499;
t305 = t443 * t395 + t351;
t303 = -t505 * t358 - t499;
t302 = t437 * t395 - t351;
t301 = pkin(4) * t524 - t450;
t295 = -pkin(4) * t387 - t342 * t417 + t338;
t292 = qJ(6) * t489 + t300;
t285 = t348 + (-qJ(6) * t395 - t344) * t417 + t505 * t524;
t282 = t443 * t388 + (t444 * qJD(4) - t468) * t395 + t325;
t280 = -pkin(4) * t389 - t438;
t279 = -t510 + t437 * t388 + (qJD(4) * t507 - qJD(2) + t468) * t395;
t277 = t431 + t513;
t273 = qJ(6) * t454 + (-qJD(4) * t352 - t440) * t415 + t456 + t513;
t272 = qJ(6) * t455 - t505 * t389 + (-t343 + t440) * t417 + t457;
t1 = [(t387 * t388 + t395 * t425) * MDP(8) + (-t395 * t372 - t387 * t389 + t388 * t533 + t425 * t524) * MDP(9) + (t372 * t408 + t389 * t398) * MDP(13) + t398 * t388 * MDP(14) + (-t317 * t488 - t435 * t358) * MDP(15) + ((-t356 * t417 - t358 * t415) * t388 + (t494 - t318 * t417 + (t356 * t415 - t358 * t417) * qJD(4)) * t395) * MDP(16) + (t317 * t524 + t358 * t389 + t395 * t365 - t435 * t378) * MDP(17) + (t318 * t524 - t356 * t389 - t395 * t363 - t436 * t378) * MDP(18) + (-t372 * t524 + t378 * t389) * MDP(19) + (t297 * t389 + t309 * t489 + t351 * t318 + t325 * t356 + t450 * t372 + t438 * t378 - t436 * t445 + t446 * t524) * MDP(20) + (-t298 * t389 + t309 * t488 - t351 * t317 + t325 * t358 - t475 * t372 - t431 * t378 - t432 * t524 + t435 * t445) * MDP(21) + (t276 * t524 + t278 * t489 - t280 * t378 + t282 * t356 - t290 * t389 + t436 * t299 - t301 * t372 + t305 * t318) * MDP(22) + (-t277 * t356 + t280 * t358 - t300 * t318 - t301 * t317 + t442 * t388 + (-t275 * t415 + t276 * t417 + (-t290 * t415 - t291 * t417) * qJD(4)) * t395) * MDP(23) + (-t275 * t524 + t277 * t378 - t278 * t488 - t282 * t358 + t291 * t389 + t435 * t299 + t300 * t372 + t305 * t317) * MDP(24) + (t275 * t300 + t276 * t301 + t277 * t291 + t278 * t305 + t280 * t290 + t282 * t299) * MDP(25) + (t270 * t524 - t272 * t378 - t274 * t489 - t279 * t356 - t281 * t389 - t285 * t372 - t436 * t286 - t302 * t318) * MDP(26) + (-t271 * t524 + t273 * t378 + t274 * t488 + t279 * t358 + t283 * t389 - t435 * t286 + t292 * t372 - t302 * t317) * MDP(27) + (-t272 * t358 + t273 * t356 + t285 * t317 + t292 * t318 + (-t281 * t417 + t283 * t415) * t388 + (-t270 * t417 + t271 * t415 + (t281 * t415 + t283 * t417) * qJD(4)) * t395) * MDP(28) + (t270 * t285 + t271 * t292 + t272 * t281 + t273 * t283 + t274 * t302 + t279 * t286) * MDP(29) + 0.2e1 * t461 * t529 + (t388 * MDP(10) - t389 * MDP(11) - t325 * MDP(13) + (t408 * t533 - t324) * MDP(14)) * qJD(3); t477 * MDP(20) + t478 * MDP(23) - t459 * t453 + (-t467 + t286 * MDP(29) - t451 * t358 + (-MDP(20) - t459) * t356) * t387 - qJD(1) ^ 2 * t529 + (t387 * MDP(13) + t533 * MDP(14) + (t395 * MDP(13) + MDP(14) * t524) * qJD(1)) * qJD(3) + (-t372 * MDP(21) + t526 * MDP(25) + t318 * MDP(28) + t527 * MDP(29) + (-qJD(4) * MDP(20) + t531 * t358 + t459 * t533) * t378) * t415 + ((t317 + t492) * MDP(23) + (t447 - t492) * MDP(28) + t528 * MDP(29) + (-qJD(4) * MDP(21) + t451 * t533) * t378 + t530) * t417 + t522 * t474; (-t309 + t532) * MDP(13) + (t358 * t448 - t494) * MDP(15) + ((-t317 + t492) * t417 - t358 * t449 + t478) * MDP(16) + t474 * MDP(17) + (-t453 + t477) * MDP(18) + (-pkin(3) * t318 - t309 * t417 - t346 * t356 + (t338 + (-pkin(8) * qJD(4) - t342) * t417) * t378 + t434 * t415) * MDP(20) + (pkin(3) * t317 + t309 * t415 - t346 * t358 + (t458 + t476) * t378 + t434 * t417) * MDP(21) + (-t278 * t417 + t318 * t399 + (-pkin(8) * t470 + t295) * t378 - t479 * t356 + t508 * t415) * MDP(22) + (t293 * t356 - t295 * t358 + ((qJD(4) * t358 - t318) * pkin(8) + t526) * t417 + (t447 * pkin(8) + t517) * t415) * MDP(23) + (-t278 * t415 + t317 * t399 + (-t293 - t458) * t378 + t479 * t358 - t508 * t417) * MDP(24) + (t278 * t399 - t290 * t295 - t291 * t293 - t479 * t299 + (t442 * qJD(4) + t275 * t417 + t276 * t415) * pkin(8)) * MDP(25) + (t274 * t417 - t286 * t449 - t318 * t392 + t480 * t356 - t372 * t402 + t482 * t378) * MDP(26) + (t274 * t415 + t286 * t448 - t317 * t392 - t480 * t358 + t372 * t403 + t481 * t378) * MDP(27) + (t317 * t402 + t318 * t403 + t481 * t356 + t482 * t358 + t528 * t415 - t527 * t417) * MDP(28) + (t270 * t402 + t271 * t403 + t274 * t392 - t482 * t281 + t481 * t283 - t480 * t286) * MDP(29) + (-t398 * MDP(13) - t358 * MDP(17) + t356 * MDP(18) - t378 * MDP(19) - t297 * MDP(20) + t298 * MDP(21) + t290 * MDP(22) - t291 * MDP(24) + t281 * MDP(26) - t283 * MDP(27) + MDP(9) * t387) * t387 + ((-qJD(2) - t398) * MDP(14) - t448 * MDP(17) - t387 * MDP(8) - MDP(9) * t533) * t533; t358 * t356 * MDP(15) + (t355 - t506) * MDP(16) - t424 * MDP(17) + (t491 - t318) * MDP(18) + t372 * MDP(19) + (t358 * t445 + t433) * MDP(20) + (-t356 * t445 + t423) * MDP(21) + (-t299 * t358 - t310 * t356 + 0.2e1 * t368 + t433) * MDP(22) + (pkin(4) * t317 - t500 + (t291 - t298) * t358 + (t290 - t464) * t356) * MDP(23) + (-t299 * t356 + t310 * t358 - t423 + t514) * MDP(24) + (-pkin(4) * t276 + qJ(5) * t275 - t290 * t298 + t464 * t291 - t299 * t310) * MDP(25) + (t288 * t378 + t303 * t356 + (pkin(5) + t505) * t372 + t463 * t358 + t426) * MDP(26) + (t286 * t356 - t287 * t378 - t303 * t358 - t432 + t514 + t516) * MDP(27) + (t500 - t317 * t505 + (-t283 + t288) * t358 + (-t281 + t465) * t356) * MDP(28) + (qJ(5) * t271 - t270 * t505 - t281 * t288 + t465 * t283 - t286 * t303) * MDP(29); (-t373 - t355) * MDP(24) - t373 * MDP(27) + (t422 - t496) * MDP(29) + (-MDP(27) * t358 - t463 * MDP(29) + t459 * t356 + t467) * t358 - t531 * t424 - t530; (-t317 - t493) * MDP(27) + (-t355 - t506) * MDP(28) + (t281 * t358 - t283 * t356 + t274) * MDP(29) + (-t491 - t318) * MDP(26);];
tauc  = t1;
