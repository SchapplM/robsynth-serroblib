% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:11
% EndTime: 2019-03-09 03:56:19
% DurationCPUTime: 3.98s
% Computational Cost: add. (4128->336), mult. (9244->461), div. (0->0), fcn. (6918->8), ass. (0->154)
t419 = cos(qJ(6));
t462 = qJD(6) * t419;
t414 = sin(pkin(10));
t415 = cos(pkin(10));
t418 = sin(qJ(3));
t421 = cos(qJ(3));
t436 = t414 * t421 + t415 * t418;
t381 = t436 * qJD(1);
t420 = cos(qJ(5));
t368 = t420 * t381;
t469 = qJD(1) * t421;
t470 = qJD(1) * t418;
t384 = -t414 * t470 + t415 * t469;
t417 = sin(qJ(5));
t332 = t384 * t417 + t368;
t515 = t332 * t419;
t519 = t462 + t515;
t438 = -t381 * t417 + t420 * t384;
t416 = sin(qJ(6));
t463 = qJD(6) * t416;
t435 = t414 * t418 - t415 * t421;
t457 = qJD(1) * qJD(3);
t372 = t435 * t457;
t373 = t436 * t457;
t464 = qJD(5) * t417;
t296 = -qJD(5) * t368 + t417 * t372 - t420 * t373 - t384 * t464;
t409 = qJD(3) + qJD(5);
t476 = t419 * t296 + t409 * t462;
t282 = -t438 * t463 + t476;
t279 = t282 * t419;
t320 = t409 * t416 + t419 * t438;
t494 = t296 * t416;
t283 = t320 * qJD(6) + t494;
t486 = t438 * t416;
t318 = -t419 * t409 + t486;
t518 = -t416 * t283 - t519 * t318 + t279;
t278 = t282 * t416;
t297 = qJD(5) * t438 - t420 * t372 - t373 * t417;
t293 = t416 * t297;
t513 = -qJD(6) - t332;
t477 = -t462 * t513 + t293;
t488 = t332 * t409;
t490 = t438 * t409;
t492 = t320 * t438;
t517 = (-t297 + t490) * MDP(19) - t332 ^ 2 * MDP(17) + (t332 * MDP(16) + MDP(17) * t438 + MDP(27) * t513) * t438 + (t296 + t488) * MDP(18) + (t519 * t320 + t278) * MDP(23) + (-t513 * t515 + t477 - t492) * MDP(25);
t422 = -pkin(1) - pkin(7);
t396 = qJD(1) * t422 + qJD(2);
t378 = -qJ(4) * t470 + t396 * t418;
t359 = t414 * t378;
t379 = -qJ(4) * t469 + t421 * t396;
t363 = qJD(3) * pkin(3) + t379;
t321 = t415 * t363 - t359;
t498 = pkin(8) * t384;
t308 = qJD(3) * pkin(4) + t321 - t498;
t485 = t415 * t378;
t322 = t414 * t363 + t485;
t499 = pkin(8) * t381;
t309 = t322 - t499;
t284 = t308 * t420 - t309 * t417;
t280 = -pkin(5) * t409 - t284;
t516 = t280 * t332;
t446 = t513 * t416;
t300 = pkin(5) * t438 + pkin(9) * t332;
t465 = qJD(4) * t421;
t468 = qJD(3) * t418;
t347 = -t396 * t468 + (qJ(4) * t468 - t465) * qJD(1);
t466 = qJD(4) * t418;
t467 = qJD(3) * t421;
t348 = t396 * t467 + (-qJ(4) * t467 - t466) * qJD(1);
t312 = t415 * t347 - t348 * t414;
t301 = pkin(8) * t373 + t312;
t313 = t414 * t347 + t415 * t348;
t302 = pkin(8) * t372 + t313;
t266 = t420 * (qJD(5) * t308 + t302) + t301 * t417 - t309 * t464;
t392 = pkin(3) * t470 + qJD(1) * qJ(2) + qJD(4);
t346 = pkin(4) * t381 + t392;
t512 = t332 * t346 - t266;
t508 = MDP(8) * (t418 ^ 2 - t421 ^ 2);
t493 = t318 * t438;
t506 = qJ(2) * MDP(6) + MDP(5);
t295 = t419 * t297;
t505 = -t463 * t513 - t295;
t285 = t308 * t417 + t309 * t420;
t267 = qJD(5) * t285 - t420 * t301 + t302 * t417;
t281 = pkin(9) * t409 + t285;
t288 = pkin(5) * t332 - pkin(9) * t438 + t346;
t440 = t281 * t416 - t288 * t419;
t504 = -t267 * t419 + t280 * t463 + t438 * t440;
t269 = t281 * t419 + t288 * t416;
t503 = t267 * t416 + t269 * t438 + t280 * t462;
t502 = -t346 * t438 - t267;
t382 = t414 * t468 - t415 * t467;
t383 = t436 * qJD(3);
t437 = -t417 * t435 + t420 * t436;
t306 = qJD(5) * t437 - t382 * t417 + t420 * t383;
t481 = qJ(4) - t422;
t374 = t468 * t481 - t465;
t394 = t481 * t421;
t375 = -qJD(3) * t394 - t466;
t324 = t415 * t374 - t375 * t414;
t310 = pkin(8) * t383 + t324;
t325 = t414 * t374 + t415 * t375;
t311 = pkin(8) * t382 + t325;
t393 = t481 * t418;
t344 = t393 * t414 - t415 * t394;
t326 = pkin(8) * t435 + t344;
t345 = -t415 * t393 - t414 * t394;
t327 = -pkin(8) * t436 + t345;
t439 = t326 * t420 - t327 * t417;
t270 = qJD(5) * t439 + t310 * t417 + t311 * t420;
t291 = t326 * t417 + t327 * t420;
t338 = -t417 * t436 - t420 * t435;
t482 = t418 * pkin(3) + qJ(2);
t364 = pkin(4) * t436 + t482;
t292 = pkin(5) * t437 - pkin(9) * t338 + t364;
t501 = t267 * t338 - t280 * t306 - t291 * t297 + (qJD(6) * t292 + t270) * t513 - t437 * (qJD(6) * t288 + t266);
t500 = pkin(3) * t414;
t496 = t280 * t338;
t495 = t292 * t297;
t491 = t320 * t416;
t424 = qJD(1) ^ 2;
t484 = t421 * t424;
t423 = qJD(3) ^ 2;
t483 = t422 * t423;
t328 = -t379 * t414 - t485;
t314 = t328 + t499;
t329 = t415 * t379 - t359;
t315 = t329 - t498;
t403 = pkin(3) * t415 + pkin(4);
t432 = t403 * t420 - t417 * t500;
t478 = -t432 * qJD(5) + t314 * t417 + t315 * t420;
t433 = t403 * t417 + t420 * t500;
t475 = t433 * qJD(5) + t314 * t420 - t315 * t417;
t410 = qJD(1) * qJD(2);
t452 = t421 * t457;
t474 = pkin(3) * t452 + t410;
t460 = pkin(3) * t467 + qJD(2);
t456 = 0.2e1 * qJD(1);
t352 = pkin(3) * t469 + pkin(4) * t384;
t377 = pkin(9) + t433;
t442 = qJD(6) * t377 + t300 + t352;
t340 = -pkin(4) * t372 + t474;
t349 = -pkin(4) * t382 + t460;
t441 = -t297 * t377 + t516;
t434 = t332 * t446 - t505;
t431 = -t306 * t419 - t338 * t463;
t428 = -t312 * t435 + t313 * t436 - t321 * t383 - t322 * t382;
t425 = qJD(5) * t338 - t420 * t382 - t383 * t417;
t376 = -pkin(5) - t432;
t274 = pkin(5) * t425 + pkin(9) * t306 + t349;
t273 = pkin(5) * t297 - pkin(9) * t296 + t340;
t272 = t419 * t273;
t271 = qJD(5) * t291 - t310 * t420 + t311 * t417;
t1 = [0.2e1 * t457 * t508 - t423 * t421 * MDP(10) + qJ(2) * t467 * t456 * MDP(12) + (-t421 * t483 + (-qJ(2) * t468 + qJD(2) * t421) * t456) * MDP(13) + (-t324 * t384 - t325 * t381 + t344 * t373 + t345 * t372 - t428) * MDP(14) + (t312 * t344 + t313 * t345 + t321 * t324 + t322 * t325 + t392 * t460 + t474 * t482) * MDP(15) + (t296 * t338 - t306 * t438) * MDP(16) + (-t296 * t437 - t297 * t338 + t306 * t332 - t425 * t438) * MDP(17) + (t364 * t297 + t332 * t349 + t340 * t437 + t346 * t425) * MDP(21) + (t296 * t364 - t306 * t346 + t338 * t340 + t349 * t438) * MDP(22) + (t279 * t338 + t320 * t431) * MDP(23) + (-(-t318 * t419 - t491) * t306 + (-t278 - t283 * t419 + (t318 * t416 - t320 * t419) * qJD(6)) * t338) * MDP(24) + (t282 * t437 + t295 * t338 + t320 * t425 - t431 * t513) * MDP(25) + (-t338 * t293 - t283 * t437 - t425 * t318 - (t306 * t416 - t338 * t462) * t513) * MDP(26) + (t297 * t437 - t425 * t513) * MDP(27) + (-t440 * t425 + t271 * t318 + t272 * t437 - t439 * t283 + (-t274 * t513 + t495 + (-t281 * t437 + t291 * t513 + t496) * qJD(6)) * t419 + t501 * t416) * MDP(28) + (-t269 * t425 + t271 * t320 - t439 * t282 + ((-qJD(6) * t291 + t274) * t513 - t495 - (-qJD(6) * t281 + t273) * t437 - qJD(6) * t496) * t416 + t501 * t419) * MDP(29) + 0.2e1 * t506 * t410 + (-0.2e1 * MDP(7) * t452 - t423 * MDP(9) + (qJD(2) * t456 - t483) * MDP(12)) * t418 + (-MDP(18) * t306 - MDP(19) * t425 - t271 * MDP(21) - t270 * MDP(22)) * t409; (t372 * t436 - t373 * t435 + t381 * t382 + t383 * t384) * MDP(14) + (-qJD(1) * t392 + t428) * MDP(15) + (-qJD(1) * t332 - t306 * t409) * MDP(21) + (-qJD(1) * t438 - t409 * t425) * MDP(22) + (-t283 * t338 - t293 * t437 + t306 * t318) * MDP(28) + (-t282 * t338 - t295 * t437 + t306 * t320) * MDP(29) - t506 * t424 - ((-qJD(1) * t419 - t416 * t425 - t437 * t462) * MDP(28) + (qJD(1) * t416 - t419 * t425 + t437 * t463) * MDP(29)) * t513 + (t418 * MDP(12) + t421 * MDP(13)) * (-t423 - t424); t418 * MDP(7) * t484 - t424 * t508 + ((t322 + t328) * t384 - (t321 - t329) * t381 + (t372 * t414 + t373 * t415) * pkin(3)) * MDP(14) + (-t321 * t328 - t322 * t329 + (t312 * t415 + t313 * t414 - t392 * t469) * pkin(3)) * MDP(15) + (-t332 * t352 - t409 * t475 + t502) * MDP(21) + (-t352 * t438 + t409 * t478 + t512) * MDP(22) + (t491 * t513 + t518) * MDP(24) + (t434 + t493) * MDP(26) + (t376 * t283 + t441 * t416 + t475 * t318 - (t416 * t478 - t419 * t442) * t513 + t504) * MDP(28) + (t376 * t282 + t441 * t419 + t475 * t320 - (t416 * t442 + t419 * t478) * t513 + t503) * MDP(29) + (MDP(13) * t418 * t424 - MDP(12) * t484) * qJ(2) + t517; (-t381 ^ 2 - t384 ^ 2) * MDP(14) + (t321 * t384 + t322 * t381 + t474) * MDP(15) + (t297 + t490) * MDP(21) + (t296 - t488) * MDP(22) + (t434 - t493) * MDP(28) + (-t419 * t513 ^ 2 - t293 - t492) * MDP(29); (t285 * t409 + t502) * MDP(21) + (t284 * t409 + t512) * MDP(22) + (t320 * t446 + t518) * MDP(24) + (-t446 * t513 + t295 + t493) * MDP(26) + (-pkin(5) * t283 + (-t284 * t416 + t300 * t419) * t513 - t285 * t318 + t416 * t516 - t477 * pkin(9) + t504) * MDP(28) + (-pkin(5) * t282 - (t284 * t419 + t300 * t416) * t513 - t285 * t320 + t280 * t515 + t505 * pkin(9) + t503) * MDP(29) + t517; t320 * t318 * MDP(23) + (-t318 ^ 2 + t320 ^ 2) * MDP(24) + (-t318 * t513 + t476) * MDP(25) + (-t320 * t513 - t494) * MDP(26) + t297 * MDP(27) + (-t266 * t416 - t269 * t513 - t280 * t320 + t272) * MDP(28) + (-t266 * t419 - t273 * t416 + t280 * t318 + t440 * t513) * MDP(29) + (-MDP(25) * t486 - MDP(26) * t320 - MDP(28) * t269 + MDP(29) * t440) * qJD(6);];
tauc  = t1;
