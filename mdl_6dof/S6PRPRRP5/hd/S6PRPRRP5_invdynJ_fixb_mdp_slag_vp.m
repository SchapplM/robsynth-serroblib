% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:53
% EndTime: 2019-03-08 20:16:58
% DurationCPUTime: 3.83s
% Computational Cost: add. (2334->413), mult. (4814->546), div. (0->0), fcn. (3472->10), ass. (0->186)
t360 = sin(qJ(4));
t363 = cos(qJ(4));
t420 = qJDD(2) * t363;
t425 = qJD(2) * qJD(4);
t494 = -t360 * t425 + t420;
t359 = sin(qJ(5));
t483 = pkin(5) * t359;
t493 = pkin(8) + t483;
t397 = pkin(4) * t363 + pkin(9) * t360;
t319 = t397 * qJD(4) + qJD(3);
t362 = cos(qJ(5));
t364 = cos(qJ(2));
t355 = sin(pkin(6));
t442 = qJD(1) * t355;
t361 = sin(qJ(2));
t457 = t359 * t361;
t492 = -(-t360 * t457 + t362 * t364) * t442 + t362 * t319;
t327 = pkin(4) * t360 - pkin(9) * t363 + qJ(3);
t365 = -pkin(2) - pkin(8);
t426 = t362 * qJD(4);
t429 = qJD(5) * t362;
t454 = t361 * t362;
t491 = (t359 * t364 + t360 * t454) * t442 - t363 * t365 * t426 - t359 * t319 - t327 * t429;
t455 = t360 * t362;
t445 = t359 * t327 + t365 * t455;
t412 = t364 * t442;
t394 = qJD(3) - t412;
t311 = t365 * qJD(2) + t394;
t357 = cos(pkin(6));
t460 = t357 * t360;
t490 = -qJD(1) * t460 + t311 * t363;
t354 = sin(pkin(10));
t356 = cos(pkin(10));
t458 = t357 * t364;
t307 = t354 * t458 + t356 * t361;
t462 = t355 * t363;
t276 = t307 * t360 + t354 * t462;
t305 = t354 * t361 - t356 * t458;
t278 = -t305 * t360 + t356 * t462;
t461 = t355 * t364;
t310 = t357 * t363 - t360 * t461;
t281 = -t310 * t359 + t355 * t454;
t459 = t357 * t361;
t306 = t354 * t364 + t356 * t459;
t308 = -t354 * t459 + t356 * t364;
t489 = -g(1) * (-t276 * t359 + t308 * t362) - g(2) * (t278 * t359 + t306 * t362) - g(3) * t281;
t377 = g(1) * t307 + g(2) * t305 - g(3) * t461;
t488 = -g(1) * t308 - g(2) * t306;
t452 = qJDD(1) - g(3);
t463 = t355 * t361;
t487 = -t452 * t463 - t488;
t409 = t360 * t426;
t430 = qJD(5) * t359;
t379 = t363 * t430 + t409;
t268 = t379 * qJD(2) - qJD(5) * t426 - t359 * qJDD(4) - t362 * t420;
t413 = t361 * t442;
t440 = qJD(2) * qJ(3);
t326 = t413 + t440;
t486 = (-t326 + t413 - t440) * qJD(4) - qJDD(4) * t365;
t433 = qJD(4) * t359;
t436 = qJD(2) * t363;
t323 = t362 * t436 + t433;
t269 = t323 * qJD(5) - t362 * qJDD(4) + t494 * t359;
t441 = qJD(1) * t363;
t340 = t357 * t441;
t285 = t360 * t311 + t340;
t273 = qJD(4) * pkin(9) + t285;
t437 = qJD(2) * t360;
t345 = qJD(5) + t437;
t376 = g(3) * t463 - t488;
t485 = (t345 * t365 + t273) * qJD(5) + t376;
t484 = t323 ^ 2;
t479 = qJ(6) + pkin(9);
t478 = qJ(6) * t359;
t477 = qJDD(2) * pkin(2);
t476 = t268 * t359;
t475 = t269 * t362;
t406 = t363 * t425;
t421 = qJDD(2) * t360;
t383 = t406 + t421;
t318 = qJDD(5) + t383;
t473 = t318 * t359;
t472 = t318 * t362;
t321 = t359 * t436 - t426;
t471 = t321 * t345;
t470 = t321 * t359;
t469 = t321 * t362;
t468 = t323 * t345;
t467 = t323 * t359;
t466 = t323 * t362;
t465 = t345 * t359;
t464 = t355 * t360;
t456 = t359 * t365;
t453 = t362 * t363;
t402 = pkin(5) - t456;
t428 = qJD(6) * t362;
t451 = qJ(6) * t409 - t445 * qJD(5) + (qJ(6) * t430 + t402 * qJD(4) - t428) * t363 + t492;
t418 = qJ(6) * t453;
t450 = -qJD(5) * t418 + (-qJD(6) * t363 + (qJ(6) * qJD(4) - qJD(5) * t365) * t360) * t359 - t491;
t291 = t327 * qJD(2) + t413;
t261 = -t273 * t359 + t362 * t291;
t253 = -qJ(6) * t323 + t261;
t250 = pkin(5) * t345 + t253;
t449 = -t253 + t250;
t325 = t397 * qJD(2);
t448 = t359 * t325 + t362 * t490;
t400 = qJD(5) * t479;
t447 = -t359 * t400 - t437 * t478 + t428 - t448;
t313 = t362 * t325;
t446 = -t362 * t400 - t313 - (pkin(5) * t363 + qJ(6) * t455) * qJD(2) + (-qJD(6) + t490) * t359;
t353 = t363 ^ 2;
t444 = t360 ^ 2 - t353;
t366 = qJD(4) ^ 2;
t367 = qJD(2) ^ 2;
t443 = -t366 - t367;
t439 = qJD(2) * t326;
t438 = qJD(2) * t355;
t435 = qJD(4) * t321;
t434 = qJD(4) * t323;
t432 = qJD(4) * t360;
t431 = qJD(5) * t345;
t272 = -qJD(4) * pkin(4) - t490;
t267 = pkin(5) * t321 + qJD(6) + t272;
t427 = t267 * qJD(4);
t424 = qJDD(1) * t355;
t423 = qJDD(1) * t357;
t422 = qJDD(2) * qJ(3);
t411 = t361 * t438;
t333 = qJD(1) * t411;
t404 = t364 * t424;
t387 = qJDD(3) + t333 - t404;
t287 = t365 * qJDD(2) + t387;
t403 = t363 * t423;
t258 = qJDD(4) * pkin(9) + qJD(4) * t490 + t287 * t360 + t403;
t405 = t361 * t424;
t266 = t405 + t327 * qJDD(2) + (t319 + t412) * qJD(2);
t417 = -t362 * t258 - t359 * t266 - t291 * t429;
t416 = -qJD(4) * t340 - t311 * t432 - t360 * t423;
t410 = t364 * t438;
t408 = g(3) * (pkin(2) * t461 + qJ(3) * t463);
t401 = -t365 + t483;
t398 = -t287 + t439;
t262 = t273 * t362 + t291 * t359;
t254 = -qJ(6) * t321 + t262;
t393 = t250 * t362 + t254 * t359;
t392 = t250 * t359 - t254 * t362;
t347 = pkin(5) * t362 + pkin(4);
t391 = t347 * t360 - t363 * t479;
t282 = t310 * t362 + t355 * t457;
t309 = t363 * t461 + t460;
t385 = t345 * t429 + t473;
t384 = -t345 * t430 + t472;
t382 = t273 * t430 + t417;
t275 = -t307 * t363 + t354 * t464;
t277 = t305 * t363 + t356 * t464;
t381 = g(1) * t275 - g(2) * t277 + g(3) * t309;
t380 = g(1) * t276 - g(2) * t278 + g(3) * t310;
t259 = -qJDD(4) * pkin(4) - t287 * t363 - t416;
t375 = -pkin(9) * t318 + t345 * t272;
t374 = t377 + t404;
t373 = qJDD(3) - t374;
t265 = t362 * t266;
t372 = -t262 * qJD(5) - t359 * t258 + t265;
t371 = pkin(9) * t431 + t259 - t381;
t249 = pkin(5) * t269 + qJDD(6) + t259;
t288 = t405 + t422 + (qJD(3) + t412) * qJD(2);
t369 = t394 * qJD(2) - t365 * t366 + t288 - t376 + t422;
t368 = (t466 + t470) * MDP(22) - t393 * MDP(23) + (-t362 * MDP(20) + t359 * MDP(21)) * t345;
t350 = qJDD(4) * t363;
t330 = t479 * t362;
t329 = t479 * t359;
t320 = -qJD(2) * pkin(2) + t394;
t317 = t321 ^ 2;
t316 = t362 * t327;
t298 = t307 * pkin(2);
t297 = t305 * pkin(2);
t292 = t387 - t477;
t280 = t310 * qJD(4) - t363 * t411;
t279 = -t309 * qJD(4) + t360 * t411;
t274 = -t363 * t478 + t445;
t270 = t402 * t360 + t316 - t418;
t256 = t281 * qJD(5) + t279 * t362 + t359 * t410;
t255 = -t282 * qJD(5) - t279 * t359 + t362 * t410;
t248 = -qJ(6) * t269 - qJD(6) * t321 - t382;
t247 = pkin(5) * t318 + qJ(6) * t268 - qJD(6) * t323 + t372;
t1 = [t452 * MDP(1) + (qJDD(1) * t357 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t280 - qJDD(4) * t309) * MDP(13) + (-qJD(4) * t279 - qJDD(4) * t310) * MDP(14) + (t255 * t345 + t269 * t309 + t280 * t321 + t281 * t318) * MDP(20) + (-t256 * t345 - t268 * t309 + t280 * t323 - t282 * t318) * MDP(21) + (-t255 * t323 - t256 * t321 + t268 * t281 - t269 * t282) * MDP(22) + (t247 * t281 + t248 * t282 + t249 * t309 + t250 * t255 + t254 * t256 + t267 * t280 - g(3)) * MDP(23) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t361 + t364 * t367) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t364 + t361 * t367) + ((-t292 + t439) * MDP(7) + (MDP(13) * t360 + MDP(14) * t363) * t367) * t364 + ((qJD(2) * t320 + t288) * MDP(7) + t383 * MDP(13) + t494 * MDP(14)) * t361) * t355; qJDD(2) * MDP(2) + t374 * MDP(3) + t487 * MDP(4) + (t373 - 0.2e1 * t477) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t422 - t487) * MDP(6) + (t288 * qJ(3) + t326 * qJD(3) - t292 * pkin(2) - g(1) * (qJ(3) * t308 - t298) - g(2) * (qJ(3) * t306 - t297) - t408 + (-t320 * t361 - t326 * t364) * t442) * MDP(7) + (qJDD(2) * t353 - 0.2e1 * t360 * t406) * MDP(8) + 0.2e1 * (-t360 * t420 + t444 * t425) * MDP(9) + (-t360 * t366 + t350) * MDP(10) + (-qJDD(4) * t360 - t363 * t366) * MDP(11) + (t369 * t360 - t363 * t486) * MDP(13) + (t360 * t486 + t369 * t363) * MDP(14) + (-t268 * t453 - t379 * t323) * MDP(15) + ((t467 + t469) * t432 + (t476 - t475 + (-t466 + t470) * qJD(5)) * t363) * MDP(16) + ((-t345 * t426 - t268) * t360 + (t384 + t434) * t363) * MDP(17) + ((t345 * t433 - t269) * t360 + (-t385 - t435) * t363) * MDP(18) + (qJD(4) * t345 * t363 + t318 * t360) * MDP(19) + (t316 * t318 + t492 * t345 + (-t327 * t431 + t377) * t359 + (t365 * t435 + t265 + (-qJD(4) * t272 - qJD(5) * t291 - t318 * t365 - t258) * t359 - t485 * t362) * t360 + (t321 * t413 + t272 * t429 + t259 * t359 - t365 * t269 + (-t345 * t456 + t261) * qJD(4)) * t363) * MDP(20) + (-t445 * t318 + t491 * t345 + t377 * t362 + ((-t272 * t362 + t323 * t365) * qJD(4) + t485 * t359 + t417) * t360 + (-qJD(4) * t262 + t259 * t362 + t268 * t365 - t272 * t430 + t323 * t413) * t363) * MDP(21) + (t268 * t270 - t269 * t274 - t451 * t323 - t450 * t321 + t393 * t432 + (t392 * qJD(5) - t247 * t362 - t248 * t359 + t376) * t363) * MDP(22) + (t248 * t274 + t247 * t270 - g(1) * (-t493 * t307 - t298) - g(2) * (-t493 * t305 - t297) - t408 - t401 * t360 * t427 + (t267 * pkin(5) * t429 + t249 * t401) * t363 + t450 * t254 + t451 * t250 + (-g(3) * t493 * t364 + (-g(3) * t391 + t267 * t441) * t361) * t355 + t488 * (qJ(3) + t391)) * MDP(23); qJDD(2) * MDP(5) - t367 * MDP(6) + (t333 + t373 - t477) * MDP(7) + t350 * MDP(13) - t377 * MDP(23) + (-t326 * MDP(7) + t368) * qJD(2) + (t443 * MDP(14) - t269 * MDP(20) + t268 * MDP(21) - t249 * MDP(23) + ((t467 - t469) * MDP(22) - t392 * MDP(23) + (-t359 * MDP(20) - t362 * MDP(21)) * t345) * qJD(4)) * t363 + (t443 * MDP(13) - qJDD(4) * MDP(14) + (t435 - t473) * MDP(20) + (t434 - t472) * MDP(21) + (-t475 - t476) * MDP(22) + (-t247 * t359 + t248 * t362 + t427) * MDP(23) + t368 * qJD(5)) * t360; MDP(10) * t420 - MDP(11) * t421 + qJDD(4) * MDP(12) + (qJD(4) * t285 - t363 * t398 + t381 + t416) * MDP(13) + (t398 * t360 + t380 - t403) * MDP(14) + (t345 * t466 - t476) * MDP(15) + ((-t268 - t471) * t362 + (-t269 - t468) * t359) * MDP(16) + ((-t323 * t363 + t345 * t455) * qJD(2) + t385) * MDP(17) + ((t321 * t363 - t360 * t465) * qJD(2) + t384) * MDP(18) - t345 * MDP(19) * t436 + (-t261 * t436 - pkin(4) * t269 - t285 * t321 - t313 * t345 + (t345 * t490 + t375) * t359 - t371 * t362) * MDP(20) + (pkin(4) * t268 + t262 * t436 - t285 * t323 + t345 * t448 + t359 * t371 + t362 * t375) * MDP(21) + (-t268 * t329 - t269 * t330 - t446 * t323 - t447 * t321 + (-t250 * t345 + t248) * t362 + (-t254 * t345 - t247) * t359 - t380) * MDP(22) + (t248 * t330 - t247 * t329 - t249 * t347 - g(1) * (-t275 * t347 + t276 * t479) - g(2) * (t277 * t347 - t278 * t479) - g(3) * (-t309 * t347 + t310 * t479) + (pkin(5) * t465 - t285) * t267 + t447 * t254 + t446 * t250) * MDP(23) + (t363 * t360 * MDP(8) - t444 * MDP(9)) * t367; t323 * t321 * MDP(15) + (-t317 + t484) * MDP(16) + (-t268 + t471) * MDP(17) + (-t269 + t468) * MDP(18) + t318 * MDP(19) + (t262 * t345 - t272 * t323 + t372 + t489) * MDP(20) + (t261 * t345 + t272 * t321 - g(1) * (-t276 * t362 - t308 * t359) - g(2) * (t278 * t362 - t306 * t359) + g(3) * t282 + t382) * MDP(21) + (pkin(5) * t268 - t321 * t449) * MDP(22) + (t449 * t254 + (-t267 * t323 + t247 + t489) * pkin(5)) * MDP(23); (-t317 - t484) * MDP(22) + (t250 * t323 + t254 * t321 + t249 - t381) * MDP(23);];
tau  = t1;
