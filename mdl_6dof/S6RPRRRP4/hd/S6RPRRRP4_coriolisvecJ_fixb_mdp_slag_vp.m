% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:54
% EndTime: 2019-03-09 06:09:03
% DurationCPUTime: 4.97s
% Computational Cost: add. (7225->393), mult. (18961->504), div. (0->0), fcn. (14926->8), ass. (0->180)
t423 = sin(qJ(5));
t422 = cos(pkin(10));
t428 = cos(qJ(3));
t496 = t422 * t428;
t421 = sin(pkin(10));
t425 = sin(qJ(3));
t497 = t421 * t425;
t443 = -t496 + t497;
t394 = t443 * qJD(1);
t402 = t421 * t428 + t422 * t425;
t395 = t402 * qJD(1);
t424 = sin(qJ(4));
t427 = cos(qJ(4));
t371 = t427 * t394 + t395 * t424;
t535 = qJD(5) + t371;
t544 = t535 * t423;
t548 = pkin(5) * t544;
t426 = cos(qJ(5));
t447 = -t394 * t424 + t427 * t395;
t470 = qJD(3) + qJD(4);
t363 = t423 * t470 + t426 * t447;
t547 = t363 * t544;
t481 = qJD(3) * t428;
t410 = t422 * qJD(1) * t481;
t465 = qJD(1) * t497;
t388 = -qJD(3) * t465 + t410;
t397 = t402 * qJD(3);
t389 = qJD(1) * t397;
t479 = qJD(4) * t427;
t480 = qJD(4) * t424;
t437 = -t427 * t388 + t424 * t389 + t394 * t479 + t395 * t480;
t433 = t423 * t437;
t309 = t363 * qJD(5) - t433;
t451 = t426 * t470;
t361 = t423 * t447 - t451;
t514 = pkin(7) + qJ(2);
t406 = t514 * t421;
t403 = qJD(1) * t406;
t407 = t514 * t422;
t404 = qJD(1) * t407;
t527 = -t428 * t403 - t404 * t425;
t345 = -pkin(8) * t389 - qJD(2) * t394 + qJD(3) * t527;
t359 = -pkin(8) * t395 + t527;
t357 = qJD(3) * pkin(3) + t359;
t452 = qJD(4) * t357 + t345;
t435 = t402 * qJD(2);
t434 = qJD(1) * t435;
t445 = t403 * t425 - t404 * t428;
t346 = -pkin(8) * t388 + qJD(3) * t445 - t434;
t360 = -pkin(8) * t394 - t445;
t458 = t424 * t346 - t360 * t480;
t287 = t452 * t427 + t458;
t457 = t424 * t388 + t427 * t389;
t523 = t447 * qJD(4);
t336 = t457 + t523;
t300 = t389 * pkin(3) + t336 * pkin(4) + pkin(9) * t437;
t356 = t427 * t360;
t319 = t424 * t357 + t356;
t313 = pkin(9) * t470 + t319;
t414 = -pkin(2) * t422 - pkin(1);
t405 = qJD(1) * t414 + qJD(2);
t379 = pkin(3) * t394 + t405;
t320 = pkin(4) * t371 - pkin(9) * t447 + t379;
t477 = qJD(5) * t426;
t478 = qJD(5) * t423;
t436 = t426 * t287 + t423 * t300 - t313 * t478 + t320 * t477;
t276 = -qJ(6) * t309 - qJD(6) * t361 + t436;
t308 = -qJD(5) * t451 + t426 * t437 + t447 * t478;
t293 = t313 * t426 + t320 * t423;
t299 = t426 * t300;
t431 = -qJD(5) * t293 - t287 * t423 + t299;
t274 = pkin(5) * t336 + qJ(6) * t308 - qJD(6) * t363 + t431;
t283 = -qJ(6) * t361 + t293;
t521 = t283 * t535 + t274;
t546 = t276 * t426 - t423 * t521;
t355 = t424 * t360;
t322 = t359 * t427 - t355;
t537 = -pkin(3) * t479 + t322;
t540 = t371 * t423;
t543 = -qJ(6) * t540 + t426 * qJD(6);
t307 = t426 * t308;
t492 = -t423 * t309 - t361 * t477;
t539 = t371 * t426;
t542 = -t361 * t539 - t307 + t492;
t306 = t308 * t423;
t475 = t447 * qJD(3);
t333 = t423 * t336;
t488 = t477 * t535 + t333;
t507 = t363 * t447;
t529 = t371 * t470;
t541 = (t475 - t457) * MDP(18) - t371 ^ 2 * MDP(16) + (t371 * MDP(15) + MDP(16) * t447 - t535 * MDP(26)) * t447 + (-t437 + t529) * MDP(17) + (-t306 + (t477 + t539) * t363) * MDP(22) + (t535 * t539 + t488 - t507) * MDP(24);
t538 = t379 * t371;
t340 = pkin(4) * t447 + pkin(9) * t371;
t418 = t426 * qJ(6);
t534 = -pkin(5) * t447 - t371 * t418;
t508 = t361 * t447;
t367 = -pkin(8) * t402 - t406 * t428 - t407 * t425;
t444 = t406 * t425 - t407 * t428;
t368 = -pkin(8) * t443 - t444;
t330 = -t367 * t427 + t424 * t368;
t516 = pkin(3) * t395;
t325 = t340 + t516;
t526 = t423 * t325 + t537 * t426;
t335 = t426 * t336;
t525 = t478 * t535 - t335;
t522 = (t421 ^ 2 + t422 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t288 = t424 * t345 - t427 * t346 + t357 * t480 + t360 * t479;
t520 = -t379 * t447 - t288;
t292 = -t313 * t423 + t426 * t320;
t318 = t427 * t357 - t355;
t312 = -pkin(4) * t470 - t318;
t509 = t288 * t426;
t519 = -t292 * t447 + t312 * t478 - t509;
t311 = t312 * t477;
t518 = t288 * t423 + t293 * t447 + t311;
t517 = t363 ^ 2;
t515 = pkin(3) * t427;
t513 = -qJ(6) - pkin(9);
t512 = pkin(3) * qJD(4);
t282 = -qJ(6) * t363 + t292;
t280 = pkin(5) * t535 + t282;
t510 = t280 * t426;
t506 = t363 * t423;
t376 = t402 * t427 - t424 * t443;
t499 = t376 * t423;
t331 = t367 * t424 + t368 * t427;
t327 = t426 * t331;
t415 = pkin(3) * t424 + pkin(9);
t495 = -qJ(6) - t415;
t494 = t280 - t282;
t491 = t426 * t318 + t423 * t340;
t383 = pkin(3) * t443 + t414;
t446 = -t402 * t424 - t427 * t443;
t332 = -pkin(4) * t446 - pkin(9) * t376 + t383;
t489 = t423 * t332 + t327;
t456 = qJD(5) * t495;
t486 = t423 * t456 - t526 + t543;
t324 = t426 * t325;
t485 = t426 * t456 - t324 + (-qJD(6) + t537) * t423 + t534;
t462 = qJD(5) * t513;
t484 = t423 * t462 - t491 + t543;
t459 = -t318 * t423 + t426 * t340;
t483 = -qJD(6) * t423 + t426 * t462 - t459 + t534;
t472 = qJD(1) * qJD(2);
t432 = -t406 * t481 + qJD(2) * t496 + (-qJD(2) * t421 - qJD(3) * t407) * t425;
t347 = -pkin(8) * t397 + t432;
t396 = t443 * qJD(3);
t430 = qJD(3) * t444 - t435;
t348 = pkin(8) * t396 + t430;
t296 = -qJD(4) * t330 + t347 * t427 + t348 * t424;
t343 = qJD(4) * t446 - t396 * t427 - t397 * t424;
t344 = qJD(4) * t376 - t396 * t424 + t427 * t397;
t303 = t397 * pkin(3) + pkin(4) * t344 - pkin(9) * t343;
t468 = t426 * t296 + t423 * t303 + t332 * t477;
t466 = -pkin(5) * t426 - pkin(4);
t464 = t376 * t477;
t454 = t535 * t426;
t321 = t359 * t424 + t356;
t449 = pkin(3) * t480 - t321;
t448 = t312 * t371 - t336 * t415;
t442 = -qJ(6) * t343 - qJD(6) * t376;
t440 = -t535 * t540 - t525;
t439 = t343 * t423 + t464;
t438 = t343 * t426 - t376 * t478;
t279 = pkin(5) * t309 + t288;
t297 = qJD(4) * t331 + t347 * t424 - t348 * t427;
t416 = -pkin(4) - t515;
t409 = pkin(9) * t426 + t418;
t408 = t513 * t423;
t400 = t415 * t426 + t418;
t399 = t495 * t423;
t358 = t361 ^ 2;
t329 = t426 * t332;
t304 = t361 * pkin(5) + qJD(6) + t312;
t302 = t426 * t303;
t295 = -qJ(6) * t499 + t489;
t289 = -pkin(5) * t446 - t331 * t423 - t376 * t418 + t329;
t278 = -qJ(6) * t464 + (-qJD(5) * t331 + t442) * t423 + t468;
t277 = pkin(5) * t344 - t296 * t423 + t302 + t442 * t426 + (-t327 + (qJ(6) * t376 - t332) * t423) * qJD(5);
t1 = [(t388 * t402 - t395 * t396) * MDP(8) + (-t388 * t443 - t389 * t402 + t394 * t396 - t395 * t397) * MDP(9) + (t414 * t389 + t405 * t397) * MDP(13) + (t414 * t388 - t405 * t396) * MDP(14) + (t343 * t447 - t376 * t437) * MDP(15) + (-t376 * t336 - t343 * t371 - t344 * t447 - t437 * t446) * MDP(16) + (t383 * t336 + t379 * t344 + (t397 * t371 - t389 * t446) * pkin(3)) * MDP(20) + (-t383 * t437 + t379 * t343 + (t376 * t389 + t397 * t447) * pkin(3)) * MDP(21) + (-t307 * t376 + t363 * t438) * MDP(22) + ((-t361 * t426 - t506) * t343 + (t306 - t309 * t426 + (t361 * t423 - t363 * t426) * qJD(5)) * t376) * MDP(23) + (t308 * t446 + t335 * t376 + t344 * t363 + t438 * t535) * MDP(24) + (t309 * t446 - t333 * t376 - t344 * t361 - t439 * t535) * MDP(25) + (-t336 * t446 + t344 * t535) * MDP(26) + ((-t331 * t477 + t302) * t535 + t329 * t336 - (-t313 * t477 + t299) * t446 + t292 * t344 + t297 * t361 + t330 * t309 + t376 * t311 + ((-qJD(5) * t332 - t296) * t535 - t331 * t336 - (-qJD(5) * t320 - t287) * t446 + t288 * t376 + t312 * t343) * t423) * MDP(27) + (-(-t331 * t478 + t468) * t535 - t489 * t336 + t436 * t446 - t293 * t344 + t297 * t363 - t330 * t308 + t376 * t509 + t438 * t312) * MDP(28) + (-t277 * t363 - t278 * t361 + t289 * t308 - t295 * t309 + (-t283 * t423 - t510) * t343 + (-t274 * t426 - t276 * t423 + (t280 * t423 - t283 * t426) * qJD(5)) * t376) * MDP(29) + (t276 * t295 + t283 * t278 + t274 * t289 + t280 * t277 + t279 * (pkin(5) * t499 + t330) + t304 * (pkin(5) * t439 + t297)) * MDP(30) + (t343 * MDP(17) - t344 * MDP(18) - t297 * MDP(20) - t296 * MDP(21)) * t470 + 0.2e1 * t472 * t522 + (-MDP(10) * t396 - MDP(11) * t397 + MDP(13) * t430 - MDP(14) * t432) * qJD(3); 0.2e1 * t395 * qJD(3) * MDP(13) + (t410 + (-t394 - t465) * qJD(3)) * MDP(14) + (t457 + t475 + 0.2e1 * t523) * MDP(20) + (-t437 - t529) * MDP(21) + (t440 - t508) * MDP(27) + (-t454 * t535 - t333 - t507) * MDP(28) + ((-t361 * t371 + t308) * t426 + t547 + t492) * MDP(29) + (-t304 * t447 + t521 * t426 + (-t280 * t535 + t276) * t423) * MDP(30) - qJD(1) ^ 2 * t522; (-t447 * t516 + t538 + t322 * t470 + (-t470 * t512 - t452) * t427 - t458) * MDP(21) + (-t405 * t395 - t434) * MDP(13) + (-t416 * t308 + t448 * t426 + t449 * t363 + (t415 * t478 + t526) * t535 + t518) * MDP(28) + (t416 * t309 + t448 * t423 + t449 * t361 + (-t415 * t477 + t423 * t537 - t324) * t535 + t519) * MDP(27) + (t321 * t470 + (-t395 * t371 - t470 * t480) * pkin(3) + t520) * MDP(20) + (-t280 * t454 + t308 * t399 - t309 * t400 - t486 * t361 - t485 * t363 + t546) * MDP(29) + (t405 * t394 + t443 * t472) * MDP(14) + (-t506 * t535 + t542) * MDP(23) + (t410 + (t394 - t465) * qJD(3)) * MDP(10) + (t276 * t400 + t274 * t399 + t279 * (t466 - t515) + (-t356 + (-t359 + t512) * t424 + t548) * t304 + t486 * t283 + t485 * t280) * MDP(30) + (t440 + t508) * MDP(25) + (-t394 ^ 2 + t395 ^ 2) * MDP(9) + t394 * t395 * MDP(8) + t541; (t319 * t470 + t520) * MDP(20) + (t318 * t470 - t287 + t538) * MDP(21) + (t542 - t547) * MDP(23) + (-t535 * t544 + t335 + t508) * MDP(25) + (-pkin(4) * t309 - pkin(9) * t488 + t312 * t540 - t319 * t361 - t459 * t535 + t519) * MDP(27) + (pkin(4) * t308 + pkin(9) * t525 + t312 * t539 - t319 * t363 + t491 * t535 + t518) * MDP(28) + (t308 * t408 - t309 * t409 - t484 * t361 - t483 * t363 - t535 * t510 + t546) * MDP(29) + (t276 * t409 + t274 * t408 + t279 * t466 + (-t319 + t548) * t304 + t484 * t283 + t483 * t280) * MDP(30) + t541; t363 * t361 * MDP(22) + (-t358 + t517) * MDP(23) + (t361 * t535 - t308) * MDP(24) + (t433 + (-qJD(5) + t535) * t363) * MDP(25) + t336 * MDP(26) + (t293 * t535 - t312 * t363 + t431) * MDP(27) + (t292 * t535 + t312 * t361 - t436) * MDP(28) + (pkin(5) * t308 - t361 * t494) * MDP(29) + (t494 * t283 + (-t304 * t363 + t274) * pkin(5)) * MDP(30); (-t358 - t517) * MDP(29) + (t280 * t363 + t283 * t361 + t279) * MDP(30);];
tauc  = t1;
