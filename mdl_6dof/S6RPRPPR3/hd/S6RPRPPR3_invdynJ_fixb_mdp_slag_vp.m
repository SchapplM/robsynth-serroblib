% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:34
% EndTime: 2019-03-09 02:45:38
% DurationCPUTime: 4.13s
% Computational Cost: add. (1686->411), mult. (3100->491), div. (0->0), fcn. (1721->10), ass. (0->186)
t387 = cos(pkin(9));
t344 = -t387 * pkin(1) - pkin(2);
t323 = qJD(1) * t344;
t386 = sin(pkin(9));
t343 = pkin(1) * t386 + pkin(7);
t322 = t343 * qJD(1);
t391 = sin(qJ(3));
t394 = cos(qJ(3));
t480 = -t394 * qJD(2) + t391 * t322;
t518 = qJD(4) + t480;
t321 = qJDD(1) * t344;
t507 = pkin(3) + pkin(4);
t517 = t391 * t507;
t320 = t343 * qJDD(1);
t516 = qJD(2) * qJD(3) + t320;
t515 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t287 = -qJD(3) * pkin(3) + t518;
t462 = qJD(1) * qJD(3);
t446 = t394 * t462;
t459 = qJDD(1) * t391;
t514 = t446 + t459;
t310 = qJDD(6) + t514;
t390 = sin(qJ(6));
t393 = cos(qJ(6));
t473 = qJD(1) * t391;
t339 = qJD(6) + t473;
t468 = qJD(6) * t339;
t513 = t390 * t310 + t393 * t468;
t450 = t507 * qJD(3);
t285 = -qJ(5) * t473 + t480;
t465 = -qJD(4) - t285;
t279 = -t450 - t465;
t472 = qJD(1) * t394;
t512 = t507 * qJDD(3);
t483 = -qJ(5) + t343;
t300 = t391 * qJD(2) + t394 * t322;
t286 = -qJ(5) * t472 + t300;
t382 = qJD(3) * qJ(4);
t283 = -t286 - t382;
t511 = 0.2e1 * t515;
t369 = t391 * qJ(4);
t476 = t394 * pkin(3) + t369;
t291 = t382 + t300;
t378 = qJ(1) + pkin(9);
t360 = sin(t378);
t361 = cos(t378);
t431 = g(1) * t361 + g(2) * t360;
t470 = qJD(3) * t391;
t435 = t391 * qJDD(2) - t322 * t470 + t394 * t516;
t270 = t435 + t515;
t447 = t391 * t462;
t337 = qJ(5) * t447;
t461 = qJD(1) * qJD(5);
t463 = qJ(5) * qJDD(1);
t268 = (t461 + t463) * t394 - t270 - t337;
t266 = qJDD(3) * pkin(5) - t268;
t280 = qJD(3) * pkin(5) - t283;
t304 = -t476 + t344;
t374 = t394 * pkin(4);
t298 = t374 - t304;
t432 = pkin(5) * t391 + pkin(8) * t394;
t284 = t298 + t432;
t469 = qJD(3) * t394;
t290 = -qJD(5) * t391 + t483 * t469;
t305 = t483 * t391;
t381 = -pkin(8) - t507;
t436 = -t394 * qJDD(2) + t322 * t469 + t391 * t516;
t427 = -qJDD(4) - t436;
t400 = -qJ(5) * t514 - t391 * t461 - t427;
t265 = t381 * qJDD(3) + t400;
t292 = -pkin(3) * t472 - qJ(4) * t473 + t323;
t282 = pkin(4) * t472 + qJD(5) - t292;
t274 = t432 * qJD(1) + t282;
t439 = -qJD(6) * t274 - t265;
t510 = -(qJD(6) * t284 + t290) * t339 + (qJD(3) * t280 + t439) * t391 - t266 * t394 - t305 * t310;
t457 = qJDD(3) * t343;
t509 = (qJD(1) * t304 + t292) * qJD(3) - t457;
t353 = qJ(4) * t472;
t419 = pkin(5) * t394 + t381 * t391;
t502 = g(3) * t391;
t508 = (t419 * qJD(1) + qJD(6) * t381 + t353) * t339 + t431 * t394 - t266 + t502;
t506 = g(1) * t360;
t503 = g(2) * t361;
t377 = g(3) * t394;
t501 = qJ(4) * t394;
t500 = qJDD(3) * pkin(3);
t458 = qJDD(1) * t394;
t424 = qJD(3) * qJD(6) + t458;
t467 = qJD(6) * t394;
t449 = t390 * t467;
t478 = qJD(1) * t449 + t393 * t447;
t275 = -qJDD(3) * t390 - t424 * t393 + t478;
t499 = t275 * t390;
t316 = qJD(3) * t390 + t393 * t472;
t479 = -t393 * qJDD(3) - t390 * t447;
t276 = t316 * qJD(6) + t390 * t458 + t479;
t498 = t276 * t391;
t497 = t284 * t310;
t466 = t393 * qJD(3);
t315 = t390 * t472 - t466;
t496 = t315 * t339;
t495 = t315 * t394;
t494 = t316 * t339;
t493 = t316 * t393;
t492 = t316 * t394;
t491 = t360 * t391;
t490 = t360 * t394;
t489 = t361 * t391;
t488 = t361 * t394;
t486 = t390 * t391;
t485 = t391 * t393;
t484 = t393 * t394;
t482 = t275 * t391 - t316 * t469;
t481 = t513 * t394;
t365 = t391 * qJD(4);
t477 = qJ(4) * t469 + t365;
t383 = t391 ^ 2;
t384 = t394 ^ 2;
t475 = t383 - t384;
t474 = t383 + t384;
t471 = qJD(3) * t283;
t464 = -qJD(5) - t282;
t456 = -MDP(12) + MDP(17);
t455 = MDP(14) + MDP(16);
t454 = t310 * t484;
t453 = t339 * t486;
t452 = t339 * t485;
t398 = qJD(1) ^ 2;
t451 = t391 * t394 * t398;
t392 = sin(qJ(1));
t445 = -pkin(1) * t392 + t361 * pkin(7);
t442 = qJD(1) * t298 + t282;
t423 = pkin(3) * t458 + qJ(4) * t514 + qJD(1) * t365 - t321;
t412 = pkin(4) * t458 + qJDD(5) + t423;
t413 = t419 * qJD(3);
t262 = qJD(1) * t413 + t432 * qJDD(1) + t412;
t277 = t381 * qJD(3) - t465;
t440 = qJD(6) * t277 - t262;
t434 = t391 * t450;
t395 = cos(qJ(1));
t430 = g(1) * t392 - g(2) * t395;
t397 = qJD(3) ^ 2;
t429 = t343 * t397 + t503;
t428 = t395 * pkin(1) + pkin(3) * t488 + t360 * pkin(7) + (pkin(2) + t369) * t361;
t426 = -pkin(2) - t476;
t422 = t453 - t495;
t420 = -t393 * t310 + t390 * t468;
t417 = -qJD(3) * t480 - t435;
t415 = -g(1) * t489 - g(2) * t491 + t377 + t436;
t414 = t391 * t466 + t449;
t411 = -qJDD(4) - t415;
t410 = 0.2e1 * t321 + t429;
t409 = 0.2e1 * t323 * qJD(3) - t457;
t269 = -qJD(1) * t434 + t412;
t288 = -t434 + t477;
t408 = -qJD(1) * t288 - qJDD(1) * t298 - t269 + t503;
t407 = qJD(3) * t300 - t415;
t406 = -t411 - t500;
t405 = -t381 * t310 + (-t280 + t286) * t339;
t272 = pkin(3) * t447 - t423;
t303 = pkin(3) * t470 - t477;
t404 = -qJD(1) * t303 - qJDD(1) * t304 - t272 - t429;
t271 = -t427 - t500;
t403 = t270 * t394 + t271 * t391 + (t287 * t394 - t291 * t391) * qJD(3);
t401 = (t279 * t391 - t283 * t394) * MDP(19) + (-t453 - t495) * MDP(25);
t389 = qJ(4) + pkin(5);
t331 = g(1) * t490;
t330 = g(1) * t491;
t327 = qJ(4) * t488;
t325 = qJ(4) * t490;
t319 = qJDD(3) * t394 - t391 * t397;
t318 = qJDD(3) * t391 + t394 * t397;
t317 = pkin(3) * t473 - t353;
t306 = t483 * t394;
t301 = -t507 * t473 + t353;
t296 = -t360 * t390 + t361 * t485;
t295 = -t360 * t393 - t361 * t486;
t294 = -t360 * t485 - t361 * t390;
t293 = t360 * t486 - t361 * t393;
t289 = qJD(5) * t394 + t470 * t483;
t278 = t413 + t477;
t267 = t400 - t512;
t264 = t274 * t390 + t277 * t393;
t263 = t274 * t393 - t277 * t390;
t261 = t393 * t262;
t1 = [(t320 * t474 + t403 - t431) * MDP(13) + (-t264 * t469 - g(1) * t293 - g(2) * t295 + t306 * t275 + t289 * t316 + (-(-qJD(6) * t305 + t278) * t339 - t497 + t440 * t391 + t280 * t467) * t390 + t510 * t393) * MDP(26) + (t263 * t469 - g(1) * t294 - g(2) * t296 + t261 * t391 - t306 * t276 + t289 * t315 + (t278 * t339 + t497 + (-t277 * t391 - t280 * t394 - t305 * t339) * qJD(6)) * t393 + t510 * t390) * MDP(25) + (t391 * t509 + t404 * t394 + t331) * MDP(12) + (t404 * t391 - t394 * t509 + t330) * MDP(14) + t430 * MDP(2) + (qJDD(3) * t305 - t331 + (t442 * t391 + t290) * qJD(3) + t408 * t394) * MDP(17) + (qJDD(3) * t306 + t330 + (t442 * t394 - t289) * qJD(3) - t408 * t391) * MDP(16) + (g(1) * t395 + g(2) * t392) * MDP(3) + t318 * MDP(7) + t319 * MDP(8) + (t310 * t391 + t339 * t469) * MDP(24) + ((-qJD(3) * t279 - qJDD(1) * t306 + t268 + (-qJD(3) * t305 + t289) * qJD(1)) * t394 + (-t471 - qJDD(1) * t305 - t267 + (qJD(3) * t306 - t290) * qJD(1)) * t391 + t431) * MDP(18) + (t410 * t391 + t409 * t394 - t330) * MDP(11) + (t409 * t391 - t410 * t394 + t331) * MDP(10) + (qJDD(1) * t383 + 0.2e1 * t391 * t446) * MDP(5) + (t414 * t339 - t454 + t482) * MDP(22) + (-t275 * t484 - t414 * t316) * MDP(20) + qJDD(1) * MDP(1) + 0.2e1 * (t391 * t458 - t475 * t462) * MDP(6) + (-t422 * qJD(3) + t481 + t498) * MDP(23) + ((t315 * t393 + t316 * t390) * t470 + (t499 - t276 * t393 + (t315 * t390 - t493) * qJD(6)) * t394) * MDP(21) + (t430 + (t386 ^ 2 + t387 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t267 * t305 + t279 * t290 - t268 * t306 + t283 * t289 + t269 * t298 + t282 * t288 - g(1) * (-qJ(5) * t361 + t445) - g(2) * (pkin(4) * t488 + t428) + (-g(1) * (t426 - t374) + g(2) * qJ(5)) * t360) * MDP(19) + (-g(1) * t445 - g(2) * t428 + t272 * t304 + t292 * t303 + t403 * t343 - t426 * t506) * MDP(15); (qJDD(2) - g(3)) * MDP(4) + (t270 * t391 - t271 * t394 - g(3)) * MDP(15) + (-t267 * t394 - t268 * t391 - g(3)) * MDP(19) + (t481 - t498) * MDP(25) + (-t339 * t449 + t454 + t482) * MDP(26) + (MDP(10) - t456) * t319 + (-MDP(11) + t455) * t318 + ((t287 * t391 + t291 * t394) * MDP(15) - MDP(26) * t452 + t401) * qJD(3); -MDP(5) * t451 + t475 * MDP(6) * t398 + MDP(7) * t459 + MDP(8) * t458 + qJDD(3) * MDP(9) + (-t323 * t473 + t407) * MDP(10) + (t502 + (-qJD(1) * t323 + t431) * t394 + t417) * MDP(11) + (0.2e1 * t500 - qJDD(4) + (-t292 * t391 + t317 * t394) * qJD(1) + t407) * MDP(12) + ((qJD(1) * t317 - g(3)) * t391 + (qJD(1) * t292 - t431) * t394 - t417 + t511) * MDP(14) + (t270 * qJ(4) - t271 * pkin(3) - t292 * t317 - t287 * t300 - g(1) * (-pkin(3) * t489 + t327) - g(2) * (-pkin(3) * t491 + t325) - g(3) * t476 + t518 * t291) * MDP(15) + (qJD(3) * t285 + t337 + (-qJD(1) * t301 - g(3)) * t391 + (t464 * qJD(1) - t431 - t463) * t394 + t435 + t511) * MDP(16) + (-qJ(5) * t459 - qJD(3) * t286 - 0.2e1 * t512 + ((-qJ(5) * qJD(3) + t301) * t394 + t464 * t391) * qJD(1) - t411) * MDP(17) + (-t267 * t507 - t268 * qJ(4) - t279 * t286 - t282 * t301 - g(1) * t327 - g(2) * t325 - g(3) * (t374 + t476) + t465 * t283 + t431 * t517) * MDP(19) + (t339 * t493 - t499) * MDP(20) + ((-t275 - t496) * t393 + (-t276 - t494) * t390) * MDP(21) + ((-t452 + t492) * qJD(1) - t513) * MDP(22) + (t422 * qJD(1) + t420) * MDP(23) - t339 * MDP(24) * t472 + (-t263 * t472 - t389 * t276 + t465 * t315 + t405 * t390 - t393 * t508) * MDP(25) + (t264 * t472 + t389 * t275 + t465 * t316 + t390 * t508 + t405 * t393) * MDP(26) + ((-pkin(3) * t391 + t501) * MDP(13) + (-t501 + t517) * MDP(18)) * qJDD(1); (-qJD(3) * t291 + t406) * MDP(15) + (-qJDD(3) * pkin(4) - qJ(5) * t446 + t406 + t471) * MDP(19) + (qJD(3) * t315 - t513) * MDP(25) + (qJD(3) * t316 + t420) * MDP(26) + t455 * (-t383 * t398 - t397) + t456 * (qJDD(3) + t451) + ((-MDP(19) * qJ(5) + MDP(13) - MDP(18)) * qJDD(1) + (t292 * MDP(15) + t464 * MDP(19) + (-MDP(25) * t393 + MDP(26) * t390) * t339) * qJD(1)) * t391; (t412 - t503 + t506) * MDP(19) - t420 * MDP(25) - t513 * MDP(26) - t474 * MDP(18) * t398 + (t391 * MDP(16) - MDP(17) * t394) * qJDD(1) + ((-t452 - t492) * MDP(26) + (0.2e1 * t394 * MDP(16) + (-t507 * MDP(19) + 0.2e1 * MDP(17)) * t391) * qJD(3) + t401) * qJD(1); t316 * t315 * MDP(20) + (-t315 ^ 2 + t316 ^ 2) * MDP(21) + (t478 - t496) * MDP(22) + (t479 - t494) * MDP(23) + t310 * MDP(24) + (-g(1) * t295 + g(2) * t293 + t264 * t339 + t280 * t316 + t261) * MDP(25) + (g(1) * t296 - g(2) * t294 + t263 * t339 - t280 * t315) * MDP(26) + (-MDP(22) * t458 + (-t265 - t377) * MDP(26) + (-MDP(22) * qJD(3) + MDP(23) * t472 - MDP(25) * t277 - MDP(26) * t274) * qJD(6)) * t393 + (-qJDD(3) * MDP(22) + t424 * MDP(23) + (t439 - t377) * MDP(25) + t440 * MDP(26)) * t390;];
tau  = t1;
