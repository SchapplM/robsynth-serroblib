% Calculate Coriolis joint torque vector for
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:32
% EndTime: 2021-01-16 00:21:45
% DurationCPUTime: 5.13s
% Computational Cost: add. (4140->426), mult. (10213->563), div. (0->0), fcn. (6884->6), ass. (0->185)
t435 = cos(qJ(3));
t434 = sin(qJ(2));
t436 = cos(qJ(2));
t516 = t435 * t436;
t452 = pkin(3) * t434 - pkin(8) * t516;
t535 = pkin(7) + pkin(8);
t481 = qJD(3) * t535;
t454 = pkin(2) * t434 - pkin(7) * t436;
t389 = t454 * qJD(1);
t433 = sin(qJ(3));
t496 = qJD(1) * t434;
t476 = t433 * t496;
t500 = pkin(6) * t476 + t435 * t389;
t557 = qJD(1) * t452 + t435 * t481 + t500;
t368 = t433 * t389;
t518 = t434 * t435;
t519 = t433 * t436;
t556 = -t368 - (-pkin(6) * t518 - pkin(8) * t519) * qJD(1) - t433 * t481;
t489 = qJD(3) * t434;
t471 = qJD(1) * t489;
t485 = qJD(1) * qJD(2);
t472 = t436 * t485;
t553 = qJD(2) * qJD(3) + t472;
t352 = -t433 * t471 + t553 * t435;
t492 = qJD(2) * t435;
t382 = -t476 + t492;
t494 = qJD(2) * t433;
t383 = t435 * t496 + t494;
t432 = sin(qJ(4));
t534 = cos(qJ(4));
t474 = t534 * qJD(4);
t482 = t553 * t433 + t435 * t471;
t487 = qJD(4) * t432;
t297 = -t534 * t352 - t382 * t474 + t383 * t487 + t432 * t482;
t450 = t432 * t382 + t383 * t534;
t298 = qJD(4) * t450 + t432 * t352 + t534 * t482;
t335 = -t534 * t382 + t383 * t432;
t333 = t335 ^ 2;
t469 = MDP(22) * t496;
t495 = qJD(1) * t436;
t418 = -qJD(3) + t495;
t406 = -qJD(4) + t418;
t529 = t450 * t406;
t530 = t335 * t406;
t536 = t450 ^ 2;
t555 = qJD(2) * t469 + t335 * MDP(18) * t450 + (-t298 - t529) * MDP(21) + (-t297 - t530) * MDP(20) + (-t333 + t536) * MDP(19);
t554 = t335 * qJ(5);
t522 = t432 * t433;
t449 = t435 * t534 - t522;
t538 = qJD(3) + qJD(4);
t540 = t534 * qJD(3) + t474;
t506 = -t435 * t540 + t449 * t495 + t522 * t538;
t385 = t432 * t435 + t433 * t534;
t347 = t538 * t385;
t505 = -t385 * t495 + t347;
t427 = pkin(6) * t495;
t490 = qJD(3) * t433;
t544 = -t427 + (-t433 * t495 + t490) * pkin(3);
t488 = qJD(3) * t435;
t477 = t434 * t488;
t491 = qJD(2) * t436;
t480 = t433 * t491;
t552 = t477 + t480;
t399 = -qJD(2) * pkin(2) + pkin(6) * t496;
t353 = -pkin(3) * t382 + t399;
t400 = qJD(2) * pkin(7) + t427;
t394 = -pkin(2) * t436 - pkin(7) * t434 - pkin(1);
t374 = t394 * qJD(1);
t521 = t433 * t374;
t343 = t400 * t435 + t521;
t392 = t454 * qJD(2);
t375 = qJD(1) * t392;
t473 = t434 * t485;
t458 = pkin(6) * t473;
t504 = -t435 * t375 - t433 * t458;
t445 = -qJD(3) * t343 - t504;
t295 = pkin(3) * t473 - pkin(8) * t352 + t445;
t451 = t374 * t488 + t433 * t375 - t400 * t490;
t444 = -t435 * t458 + t451;
t300 = -pkin(8) * t482 + t444;
t342 = t435 * t374 - t400 * t433;
t319 = -pkin(8) * t383 + t342;
t313 = -pkin(3) * t418 + t319;
t320 = pkin(8) * t382 + t343;
t457 = -t432 * t295 - t534 * t300 - t313 * t474 + t320 * t487;
t551 = t335 * t353 + t457;
t549 = -0.2e1 * t485;
t548 = MDP(4) * t434;
t430 = t434 ^ 2;
t547 = MDP(5) * (-t436 ^ 2 + t430);
t546 = qJ(5) * t450;
t468 = -pkin(4) * t335 - qJD(5);
t310 = t353 - t468;
t545 = t310 * t450;
t381 = t435 * t394;
t532 = pkin(6) * t433;
t341 = -pkin(8) * t518 + t381 + (-pkin(3) - t532) * t436;
t420 = pkin(6) * t516;
t499 = t433 * t394 + t420;
t520 = t433 * t434;
t348 = -pkin(8) * t520 + t499;
t507 = t432 * t341 + t534 * t348;
t401 = t535 * t433;
t402 = t535 * t435;
t501 = -t432 * t401 + t534 * t402;
t416 = pkin(4) * t473;
t543 = -t450 * qJD(5) + t416;
t542 = qJD(4) * t501 + t556 * t432 + t557 * t534;
t541 = -t401 * t474 - t402 * t487 - t557 * t432 + t556 * t534;
t316 = t432 * t320;
t287 = t534 * t313 - t316;
t539 = t287 * MDP(23);
t318 = t534 * t320;
t288 = t432 * t313 + t318;
t466 = t534 * t295 - t432 * t300;
t442 = -qJD(4) * t288 + t466;
t537 = -t353 * t450 + t442;
t533 = pkin(3) * t406;
t531 = t297 * qJ(5);
t528 = t352 * t433;
t527 = t382 * t418;
t526 = t383 * t418;
t525 = t399 * t433;
t524 = t399 * t435;
t523 = t418 * t435;
t437 = qJD(2) ^ 2;
t517 = t434 * t437;
t515 = t436 * t437;
t438 = qJD(1) ^ 2;
t514 = t436 * t438;
t282 = t287 - t546;
t281 = -pkin(4) * t406 + t282;
t513 = t281 - t282;
t512 = pkin(4) * t496 - t506 * qJ(5) + t385 * qJD(5) + t542;
t511 = -t505 * qJ(5) + qJD(5) * t449 + t541;
t510 = t534 * t319 - t316;
t509 = t505 * pkin(4) + t544;
t503 = t433 * t392 + t394 * t488;
t493 = qJD(2) * t434;
t502 = t435 * t392 + t493 * t532;
t393 = pkin(3) * t520 + t434 * pkin(6);
t354 = t552 * pkin(3) + pkin(6) * t491;
t425 = -pkin(3) * t435 - pkin(2);
t479 = t433 * t489;
t478 = t436 * t490;
t470 = MDP(15) * t496;
t467 = pkin(1) * t549;
t465 = -t319 * t432 - t318;
t463 = t534 * t341 - t348 * t432;
t461 = -t534 * t401 - t402 * t432;
t460 = -t382 + t492;
t459 = -t383 + t494;
t456 = t534 * t491;
t455 = t432 * t473;
t453 = qJD(1) * t430 - t418 * t436;
t332 = pkin(3) * t482 + pkin(6) * t472;
t448 = qJ(5) * t298 + t457;
t305 = t452 * qJD(2) + (-t420 + (pkin(8) * t434 - t394) * t433) * qJD(3) + t502;
t307 = -t552 * pkin(8) + (-t434 * t492 - t478) * pkin(6) + t503;
t447 = t432 * t305 + t534 * t307 + t341 * t474 - t348 * t487;
t286 = t298 * pkin(4) + t332;
t278 = -qJD(5) * t335 - t448;
t441 = -qJD(4) * t507 + t534 * t305 - t432 * t307;
t439 = t442 + t531;
t424 = pkin(3) * t534 + pkin(4);
t378 = t474 * t533;
t363 = t449 * t434;
t362 = t385 * t434;
t357 = -pkin(4) * t449 + t425;
t344 = pkin(4) * t362 + t393;
t324 = qJ(5) * t449 + t501;
t323 = -qJ(5) * t385 + t461;
t315 = pkin(3) * t383 + pkin(4) * t450;
t309 = t433 * t456 - t432 * t479 - t487 * t520 + (t432 * t491 + t434 * t540) * t435;
t308 = t347 * t434 + t432 * t480 - t435 * t456;
t302 = pkin(4) * t309 + t354;
t301 = -qJ(5) * t362 + t507;
t299 = -pkin(4) * t436 - qJ(5) * t363 + t463;
t285 = t510 - t546;
t284 = t465 + t554;
t283 = t288 - t554;
t280 = -qJ(5) * t309 - qJD(5) * t362 + t447;
t279 = pkin(4) * t493 + t308 * qJ(5) - t363 * qJD(5) + t441;
t277 = t439 + t543;
t1 = [t547 * t549 + (-pkin(6) * t515 + t434 * t467) * MDP(9) + (pkin(6) * t517 + t436 * t467) * MDP(10) + (t352 * t518 + (t435 * t491 - t479) * t383) * MDP(11) + ((t382 * t435 - t383 * t433) * t491 + (-t435 * t482 - t528 + (-t433 * t382 - t383 * t435) * qJD(3)) * t434) * MDP(12) + (t418 * t479 - t352 * t436 + (t383 * t434 + t435 * t453) * qJD(2)) * MDP(13) + (t418 * t477 + t482 * t436 + (t382 * t434 - t433 * t453) * qJD(2)) * MDP(14) + (-(-t394 * t490 + t502) * t418 + (pkin(6) * t482 + t399 * t488 + (t381 * qJD(1) + t342) * qJD(2)) * t434 + ((-pkin(6) * t382 + t525) * qJD(2) + (t521 + (pkin(6) * t418 + t400) * t435) * qJD(3) + t504) * t436) * MDP(16) + ((-pkin(6) * t478 + t503) * t418 + t451 * t436 + (pkin(6) * t352 - t399 * t490) * t434 + ((pkin(6) * t383 + t524) * t436 + (-pkin(6) * t523 - qJD(1) * t499 - t343) * t434) * qJD(2)) * MDP(17) + (-t297 * t363 - t308 * t450) * MDP(18) + (t297 * t362 - t298 * t363 + t308 * t335 - t309 * t450) * MDP(19) + (t297 * t436 + t308 * t406) * MDP(20) + (t298 * t436 + t309 * t406) * MDP(21) + (t393 * t298 + t353 * t309 + t332 * t362 + t354 * t335 - t406 * t441 - t436 * t442 + t463 * t473) * MDP(23) + (-t393 * t297 - t353 * t308 + t332 * t363 + t354 * t450 + t447 * t406 - t457 * t436) * MDP(24) + (-t277 * t436 - t279 * t406 + t286 * t362 + t298 * t344 + t302 * t335 + t309 * t310) * MDP(25) + (t278 * t436 + t280 * t406 + t286 * t363 - t297 * t344 + t302 * t450 - t308 * t310) * MDP(26) + (-t277 * t363 - t278 * t362 - t279 * t450 - t280 * t335 + t281 * t308 - t283 * t309 + t297 * t299 - t298 * t301) * MDP(27) + (t277 * t299 + t278 * t301 + t279 * t281 + t280 * t283 + t286 * t344 + t302 * t310) * MDP(28) + 0.2e1 * t472 * t548 + MDP(6) * t515 - MDP(7) * t517 + ((-t418 - t495) * MDP(15) + (qJD(1) * t363 + t450) * MDP(20) + (-qJD(1) * t362 - t335) * MDP(21) + (-t406 - t495) * MDP(22) + t539 + (-qJD(1) * t507 - t288) * MDP(24) + (qJD(1) * t299 + t281) * MDP(25) + (-qJD(1) * t301 - t283) * MDP(26)) * t493; -t514 * t548 + t438 * t547 + (-t383 * t523 + t528) * MDP(11) + ((t352 - t527) * t435 + (-t482 + t526) * t433) * MDP(12) + (-t418 * t488 + (t418 * t516 + t434 * t459) * qJD(1)) * MDP(13) + (t418 * t490 + (-t418 * t519 + t434 * t460) * qJD(1)) * MDP(14) + t418 * t470 + (-pkin(2) * t482 + t500 * t418 + (pkin(7) * t523 + t525) * qJD(3) + ((-pkin(7) * t494 - t342) * t434 + (-pkin(6) * t460 - t525) * t436) * qJD(1)) * MDP(16) + (-pkin(2) * t352 - t368 * t418 + (-pkin(7) * t418 * t433 + t524) * qJD(3) + (-t399 * t516 + (-pkin(7) * t492 + t343) * t434 + (t418 * t518 + t436 * t459) * pkin(6)) * qJD(1)) * MDP(17) + (-t297 * t385 - t450 * t506) * MDP(18) + (-t297 * t449 - t298 * t385 + t335 * t506 - t450 * t505) * MDP(19) + (t425 * t298 - t332 * t449 + t544 * t335 + t505 * t353 + t461 * t473) * MDP(23) + (-t425 * t297 + t332 * t385 - t506 * t353 + t544 * t450) * MDP(24) + (-t286 * t449 + t298 * t357 + t505 * t310 + t509 * t335) * MDP(25) + (t286 * t385 - t297 * t357 - t506 * t310 + t450 * t509) * MDP(26) + (-t277 * t385 + t278 * t449 + t281 * t506 - t283 * t505 + t297 * t323 - t298 * t324 - t335 * t511 + t450 * t512) * MDP(27) + (t277 * t323 + t278 * t324 - t281 * t512 + t283 * t511 + t286 * t357 + t310 * t509) * MDP(28) + ((qJD(2) * t385 - t450) * MDP(20) + (qJD(2) * t449 + t335) * MDP(21) - t539 + (-qJD(2) * t501 + t288) * MDP(24) + (qJD(2) * t323 - t281) * MDP(25) + (-qJD(2) * t324 + t283) * MDP(26)) * t496 + (t506 * MDP(20) + t505 * MDP(21) + MDP(23) * t542 + MDP(24) * t541 + t512 * MDP(25) + t511 * MDP(26) + t469) * t406 + (MDP(9) * t434 * t438 + MDP(10) * t514) * pkin(1); -t383 * t382 * MDP(11) + (-t382 ^ 2 + t383 ^ 2) * MDP(12) + (t352 + t527) * MDP(13) + (-t482 - t526) * MDP(14) + qJD(2) * t470 + (-t343 * t418 - t383 * t399 + t445) * MDP(16) + (-t342 * t418 - t382 * t399 - t444) * MDP(17) + (t465 * t406 + (-t383 * t335 + t406 * t487 + t473 * t534) * pkin(3) + t537) * MDP(23) + (t378 - t510 * t406 + (-t383 * t450 - t455) * pkin(3) + t551) * MDP(24) + (t424 * t473 + t531 + t284 * t406 - t545 - t315 * t335 + (-t318 + (-t313 + t533) * t432) * qJD(4) + t466 + t543) * MDP(25) + (-pkin(3) * t455 - t285 * t406 + t310 * t335 - t315 * t450 - t278 + t378) * MDP(26) + (-t281 * t335 + t283 * t450 + t284 * t450 + t285 * t335 + t424 * t297 + (-t298 * t432 + (-t335 * t534 + t432 * t450) * qJD(4)) * pkin(3)) * MDP(27) + (t277 * t424 - t281 * t284 - t283 * t285 - t310 * t315 + (t278 * t432 + (-t281 * t432 + t283 * t534) * qJD(4)) * pkin(3)) * MDP(28) + t555; (-t288 * t406 + t537) * MDP(23) + (-t287 * t406 + t551) * MDP(24) + (-t283 * t406 + 0.2e1 * t416 + (-t310 + t468) * t450 + t439) * MDP(25) + (-pkin(4) * t536 - t282 * t406 + (qJD(5) + t310) * t335 + t448) * MDP(26) + (pkin(4) * t297 - t335 * t513) * MDP(27) + (t513 * t283 + (t277 - t545) * pkin(4)) * MDP(28) + t555; (t298 - t529) * MDP(25) + (-t297 + t530) * MDP(26) + (-t333 - t536) * MDP(27) + (t281 * t450 + t283 * t335 + t286) * MDP(28);];
tauc = t1;
