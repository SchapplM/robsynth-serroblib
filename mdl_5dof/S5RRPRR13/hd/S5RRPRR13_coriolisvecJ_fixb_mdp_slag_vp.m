% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:34:01
% EndTime: 2019-12-31 20:34:11
% DurationCPUTime: 5.88s
% Computational Cost: add. (3322->401), mult. (8493->557), div. (0->0), fcn. (6337->8), ass. (0->186)
t439 = cos(qJ(2));
t497 = qJD(1) * t439;
t422 = -qJD(4) + t497;
t416 = -qJD(5) + t422;
t432 = sin(pkin(9));
t436 = sin(qJ(2));
t498 = qJD(1) * t436;
t478 = t432 * t498;
t433 = cos(pkin(9));
t486 = t433 * qJD(2);
t391 = t478 - t486;
t477 = t433 * t498;
t496 = qJD(2) * t432;
t393 = t477 + t496;
t435 = sin(qJ(4));
t438 = cos(qJ(4));
t340 = t438 * t391 + t393 * t435;
t437 = cos(qJ(5));
t339 = t391 * t435 - t393 * t438;
t434 = sin(qJ(5));
t516 = t339 * t434;
t534 = -t437 * t340 + t516;
t533 = t416 * t534;
t456 = pkin(2) * t436 - qJ(3) * t439;
t401 = t456 * qJD(1);
t382 = t432 * t401;
t512 = t433 * t436;
t513 = t432 * t439;
t450 = -pkin(6) * t512 - pkin(7) * t513;
t346 = qJD(1) * t450 + t382;
t535 = qJD(3) * t433 - t346;
t532 = t339 * t422;
t531 = t340 * t422;
t454 = t339 * t437 + t340 * t434;
t530 = t416 * t454;
t509 = t438 * t433;
t515 = t432 * t435;
t398 = -t509 + t515;
t448 = t398 * t439;
t502 = qJD(1) * t448 - t398 * qJD(4);
t399 = t432 * t438 + t433 * t435;
t449 = t399 * t439;
t501 = -qJD(1) * t449 + t399 * qJD(4);
t406 = -pkin(2) * t439 - qJ(3) * t436 - pkin(1);
t384 = t406 * qJD(1);
t426 = pkin(6) * t497;
t411 = qJD(2) * qJ(3) + t426;
t348 = t433 * t384 - t411 * t432;
t481 = pkin(3) * t497;
t309 = -pkin(7) * t393 + t348 - t481;
t349 = t432 * t384 + t433 * t411;
t313 = -pkin(7) * t391 + t349;
t282 = t309 * t435 + t313 * t438;
t277 = -pkin(8) * t340 + t282;
t489 = qJD(5) * t434;
t275 = t277 * t489;
t425 = pkin(6) * t498;
t519 = qJD(2) * pkin(2);
t470 = qJD(3) - t519;
t405 = t425 + t470;
t355 = pkin(3) * t391 + t405;
t299 = pkin(4) * t340 + t355;
t528 = -t299 * t534 + t275;
t483 = qJD(1) * qJD(2);
t475 = t439 * t483;
t460 = t432 * t475;
t490 = qJD(4) * t438;
t302 = -t391 * t490 + t475 * t509 + (-qJD(4) * t393 - t460) * t435;
t377 = qJD(2) * t456 - qJD(3) * t436;
t366 = t377 * qJD(1);
t404 = (qJD(3) - t425) * qJD(2);
t327 = t433 * t366 - t404 * t432;
t511 = t433 * t439;
t451 = pkin(3) * t436 - pkin(7) * t511;
t445 = t451 * qJD(2);
t306 = qJD(1) * t445 + t327;
t328 = t432 * t366 + t433 * t404;
t310 = -pkin(7) * t460 + t328;
t467 = t438 * t306 - t310 * t435;
t442 = -t282 * qJD(4) + t467;
t476 = t436 * t483;
t268 = pkin(4) * t476 - pkin(8) * t302 + t442;
t444 = qJD(2) * t449;
t521 = qJD(4) * t339;
t303 = qJD(1) * t444 - t521;
t492 = qJD(4) * t435;
t447 = t435 * t306 + t309 * t490 + t438 * t310 - t313 * t492;
t269 = -pkin(8) * t303 + t447;
t469 = t437 * t268 - t434 * t269;
t527 = t299 * t454 + t469;
t474 = MDP(26) * t498;
t526 = qJD(2) * t474 + (t454 ^ 2 - t534 ^ 2) * MDP(23) + t534 * t454 * MDP(22);
t525 = -0.2e1 * t483;
t524 = MDP(4) * t436;
t523 = MDP(5) * (t436 ^ 2 - t439 ^ 2);
t390 = t433 * t406;
t347 = -pkin(7) * t512 + t390 + (-pkin(6) * t432 - pkin(3)) * t439;
t420 = pkin(6) * t511;
t362 = t432 * t406 + t420;
t514 = t432 * t436;
t354 = -pkin(7) * t514 + t362;
t503 = t435 * t347 + t438 * t354;
t520 = pkin(7) + qJ(3);
t408 = t520 * t432;
t409 = t520 * t433;
t500 = -t435 * t408 + t438 * t409;
t356 = pkin(6) * t478 + t433 * t401;
t329 = qJD(1) * t451 + t356;
t452 = qJD(3) * t432 + qJD(4) * t409;
t522 = -t408 * t490 + t535 * t438 + (-t329 - t452) * t435;
t468 = t434 * t302 + t437 * t303;
t272 = -qJD(5) * t454 + t468;
t281 = t438 * t309 - t313 * t435;
t276 = pkin(8) * t339 + t281;
t274 = -pkin(4) * t422 + t276;
t518 = t274 * t437;
t517 = t277 * t437;
t440 = qJD(2) ^ 2;
t510 = t436 * t440;
t508 = t439 * t440;
t441 = qJD(1) ^ 2;
t507 = t439 * t441;
t344 = t437 * t398 + t399 * t434;
t506 = -qJD(5) * t344 - t434 * t501 + t437 * t502;
t345 = -t398 * t434 + t399 * t437;
t505 = qJD(5) * t345 + t434 * t502 + t437 * t501;
t495 = qJD(2) * t436;
t480 = pkin(6) * t495;
t352 = t433 * t377 + t432 * t480;
t421 = pkin(6) * t475;
t376 = pkin(3) * t460 + t421;
t385 = t432 * t481 + t426;
t494 = qJD(2) * t439;
t427 = pkin(6) * t494;
t386 = t432 * pkin(3) * t494 + t427;
t402 = pkin(3) * t514 + t436 * pkin(6);
t491 = qJD(4) * t436;
t488 = qJD(5) * t437;
t485 = MDP(11) * qJD(1);
t484 = MDP(12) * qJD(1);
t482 = pkin(6) * t513;
t479 = t437 * t302 - t434 * t303 - t340 * t488;
t424 = -pkin(3) * t433 - pkin(2);
t473 = MDP(19) * t495;
t472 = pkin(4) * t501 - t385;
t471 = pkin(1) * t525;
t319 = t445 + t352;
t368 = t432 * t377;
t331 = qJD(2) * t450 + t368;
t466 = t438 * t319 - t331 * t435;
t465 = t438 * t347 - t354 * t435;
t464 = -t438 * t408 - t409 * t435;
t463 = t391 + t486;
t462 = -t393 + t496;
t461 = qJD(5) * t274 + t269;
t459 = -t405 + t470;
t321 = t438 * t329;
t323 = -pkin(8) * t398 + t500;
t458 = pkin(4) * t498 + pkin(8) * t502 + t399 * qJD(3) + qJD(4) * t500 + qJD(5) * t323 - t435 * t346 + t321;
t322 = -pkin(8) * t399 + t464;
t457 = -pkin(8) * t501 + qJD(5) * t322 + t522;
t266 = t274 * t434 + t517;
t374 = t398 * t436;
t284 = -pkin(4) * t439 + pkin(8) * t374 + t465;
t373 = t399 * t436;
t285 = -pkin(8) * t373 + t503;
t455 = t284 * t434 + t285 * t437;
t317 = t437 * t373 - t374 * t434;
t318 = -t373 * t434 - t374 * t437;
t446 = t435 * t319 + t438 * t331 + t347 * t490 - t354 * t492;
t271 = t339 * t489 + t479;
t367 = pkin(4) * t398 + t424;
t361 = t390 - t482;
t357 = -pkin(6) * t477 + t382;
t353 = -t433 * t480 + t368;
t351 = pkin(4) * t373 + t402;
t325 = t490 * t512 - t491 * t515 + t444;
t324 = -qJD(2) * t448 - t399 * t491;
t304 = pkin(4) * t325 + t386;
t286 = pkin(4) * t303 + t376;
t279 = qJD(5) * t318 + t324 * t434 + t437 * t325;
t278 = -qJD(5) * t317 + t324 * t437 - t325 * t434;
t273 = -pkin(8) * t325 + t446;
t270 = pkin(4) * t495 - pkin(8) * t324 - qJD(4) * t503 + t466;
t265 = -t277 * t434 + t518;
t1 = [MDP(6) * t508 - MDP(7) * t510 + t523 * t525 + (-pkin(6) * t508 + t436 * t471) * MDP(9) + (pkin(6) * t510 + t439 * t471) * MDP(10) + ((-qJD(1) * t352 - t327) * t439 + ((pkin(6) * t391 + t405 * t432) * t439 + (t348 + (t361 + 0.2e1 * t482) * qJD(1)) * t436) * qJD(2)) * MDP(11) + ((qJD(1) * t353 + t328) * t439 + ((pkin(6) * t393 + t405 * t433) * t439 + (-t349 + (-t362 + 0.2e1 * t420) * qJD(1)) * t436) * qJD(2)) * MDP(12) + (-t352 * t393 - t353 * t391 + (-t327 * t433 - t328 * t432) * t436 + (-t348 * t433 - t349 * t432 + (-t361 * t433 - t362 * t432) * qJD(1)) * t494) * MDP(13) + (t327 * t361 + t328 * t362 + t348 * t352 + t349 * t353 + (t405 + t425) * t427) * MDP(14) + (-t302 * t374 - t324 * t339) * MDP(15) + (-t302 * t373 + t303 * t374 - t324 * t340 + t325 * t339) * MDP(16) + (-t302 * t439 - t324 * t422 + (-qJD(1) * t374 - t339) * t495) * MDP(17) + (t303 * t439 + t325 * t422 + (-qJD(1) * t373 - t340) * t495) * MDP(18) + (-t422 - t497) * t473 + (-t466 * t422 - t467 * t439 + t386 * t340 + t402 * t303 + t376 * t373 + t355 * t325 + (t282 * t439 + t422 * t503) * qJD(4) + (qJD(1) * t465 + t281) * t495) * MDP(20) + (t446 * t422 + t447 * t439 - t386 * t339 + t402 * t302 - t376 * t374 + t355 * t324 + (-qJD(1) * t503 - t282) * t495) * MDP(21) + (t271 * t318 - t278 * t454) * MDP(22) + (-t271 * t317 - t272 * t318 + t278 * t534 + t279 * t454) * MDP(23) + (-t271 * t439 - t278 * t416 + (qJD(1) * t318 - t454) * t495) * MDP(24) + (t272 * t439 + t279 * t416 + (-qJD(1) * t317 + t534) * t495) * MDP(25) + (-t416 - t497) * MDP(26) * t495 + (-(t270 * t437 - t273 * t434) * t416 - t469 * t439 - t304 * t534 + t351 * t272 + t286 * t317 + t299 * t279 + (t266 * t439 + t416 * t455) * qJD(5) + ((t284 * t437 - t285 * t434) * qJD(1) + t265) * t495) * MDP(27) + (t351 * t271 - t275 * t439 + t299 * t278 + t286 * t318 - t304 * t454 + ((-qJD(5) * t285 + t270) * t416 + t268 * t439) * t434 + ((qJD(5) * t284 + t273) * t416 + t461 * t439) * t437 + (-qJD(1) * t455 - t266) * t495) * MDP(28) + 0.2e1 * t475 * t524; -t507 * t524 + t441 * t523 + ((-qJ(3) * t496 - t348) * t436 + (-pkin(6) * t463 + t432 * t459 + t356) * t439) * t485 + ((-qJ(3) * t486 + t349) * t436 + (pkin(6) * t462 + t433 * t459 - t357) * t439) * t484 + (t356 * t393 + t357 * t391 + (-qJD(3) * t391 + t348 * t497 + t328) * t433 + (qJD(3) * t393 + t349 * t497 - t327) * t432) * MDP(13) + (-t348 * t356 - t349 * t357 + (-t348 * t432 + t349 * t433) * qJD(3) + (-t327 * t432 + t328 * t433) * qJ(3) + (-t405 - t519) * t426) * MDP(14) + (t302 * t399 - t339 * t502) * MDP(15) + (-t302 * t398 - t303 * t399 + t339 * t501 - t340 * t502) * MDP(16) + (-t502 * t422 + (qJD(2) * t399 + t339) * t498) * MDP(17) + (t501 * t422 + (-qJD(2) * t398 + t340) * t498) * MDP(18) + t422 * MDP(19) * t498 + (t424 * t303 - t385 * t340 + t376 * t398 + (t321 + t452 * t438 + (-qJD(4) * t408 + t535) * t435) * t422 + t501 * t355 + (qJD(2) * t464 - t281) * t498) * MDP(20) + (t424 * t302 + t385 * t339 + t376 * t399 + t522 * t422 + t502 * t355 + (-qJD(2) * t500 + t282) * t498) * MDP(21) + (t271 * t345 - t454 * t506) * MDP(22) + (-t271 * t344 - t272 * t345 + t454 * t505 + t506 * t534) * MDP(23) + (-t506 * t416 + (qJD(2) * t345 + t454) * t498) * MDP(24) + (t505 * t416 + (-qJD(2) * t344 - t534) * t498) * MDP(25) + t416 * t474 + (t367 * t272 + t286 * t344 + (t434 * t457 + t437 * t458) * t416 + t505 * t299 - t472 * t534 + ((t322 * t437 - t323 * t434) * qJD(2) - t265) * t498) * MDP(27) + (t367 * t271 + t286 * t345 + (-t434 * t458 + t437 * t457) * t416 + t506 * t299 - t472 * t454 + (-(t322 * t434 + t323 * t437) * qJD(2) + t266) * t498) * MDP(28) + (MDP(9) * t436 * t441 + MDP(10) * t507) * pkin(1); (-t391 ^ 2 - t393 ^ 2) * MDP(13) + (t348 * t393 + t349 * t391 + t421) * MDP(14) + (t303 + t532) * MDP(20) + (t302 + t531) * MDP(21) + (t272 + t530) * MDP(27) + (t271 - t533) * MDP(28) + (t462 * t485 + t463 * t484) * t439; -t339 * t340 * MDP(15) + (t339 ^ 2 - t340 ^ 2) * MDP(16) + (t302 - t531) * MDP(17) + (-t399 * t475 + t521 + t532) * MDP(18) + qJD(1) * t473 + (-t282 * t422 + t339 * t355 + t442) * MDP(20) + (-t281 * t422 + t340 * t355 - t447) * MDP(21) + (t271 + t533) * MDP(24) + (-t272 + t530) * MDP(25) + ((-t276 * t434 - t517) * t416 - t266 * qJD(5) + (-t339 * t534 + t416 * t489 + t437 * t476) * pkin(4) + t527) * MDP(27) + ((t277 * t416 - t268) * t434 + (-t276 * t416 - t461) * t437 + (-t339 * t454 + t416 * t488 - t434 * t476) * pkin(4) + t528) * MDP(28) + t526; (t479 + t533) * MDP(24) + (-t468 + t530) * MDP(25) + (-t266 * t416 + t527) * MDP(27) + (-t265 * t416 - t434 * t268 - t437 * t269 + t528) * MDP(28) + (MDP(24) * t516 + MDP(25) * t454 - MDP(27) * t266 - MDP(28) * t518) * qJD(5) + t526;];
tauc = t1;
