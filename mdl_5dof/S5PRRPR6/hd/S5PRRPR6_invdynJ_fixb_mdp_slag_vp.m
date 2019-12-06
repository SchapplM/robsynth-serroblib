% Calculate vector of inverse dynamics joint torques for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:33:02
% EndTime: 2019-12-05 16:33:12
% DurationCPUTime: 5.27s
% Computational Cost: add. (2085->417), mult. (4974->600), div. (0->0), fcn. (3956->14), ass. (0->181)
t382 = cos(qJ(3));
t449 = qJD(2) * t382;
t362 = -qJD(5) + t449;
t374 = sin(pkin(10));
t376 = cos(pkin(10));
t443 = t376 * qJD(3);
t379 = sin(qJ(3));
t450 = qJD(2) * t379;
t338 = t374 * t450 - t443;
t448 = qJD(3) * t374;
t340 = t376 * t450 + t448;
t378 = sin(qJ(5));
t381 = cos(qJ(5));
t406 = t338 * t378 - t340 * t381;
t487 = t362 * t406;
t413 = pkin(3) * t379 - qJ(4) * t382;
t320 = t413 * qJD(3) - qJD(4) * t379;
t380 = sin(qJ(2));
t447 = qJD(3) * t379;
t436 = pkin(7) * t447;
t375 = sin(pkin(5));
t454 = qJD(1) * t375;
t383 = cos(qJ(2));
t463 = t382 * t383;
t459 = t376 * t320 + t374 * t436 - (-t374 * t463 + t376 * t380) * t454;
t486 = t374 * t320 - (t374 * t380 + t376 * t463) * t454;
t347 = qJD(2) * pkin(7) + t380 * t454;
t377 = cos(pkin(5));
t464 = t377 * t382;
t485 = qJD(1) * t464 - t379 * t347;
t484 = -qJDD(3) * pkin(3) + qJDD(4);
t384 = qJD(3) ^ 2;
t442 = qJD(1) * qJD(2);
t426 = t380 * t442;
t467 = t375 * t383;
t412 = -qJDD(1) * t467 + t375 * t426;
t477 = sin(pkin(9));
t419 = t477 * t380;
t478 = cos(pkin(9));
t420 = t478 * t383;
t322 = -t377 * t420 + t419;
t481 = g(2) * t322;
t418 = t477 * t383;
t421 = t478 * t380;
t324 = t377 * t418 + t421;
t482 = g(1) * t324;
t415 = t481 + t482;
t483 = (2 * qJDD(2) * pkin(2)) - pkin(7) * t384 + (-g(3) * t383 + t426) * t375 - t412 + t415;
t480 = pkin(8) + qJ(4);
t479 = qJD(2) * pkin(2);
t473 = t340 * t378;
t276 = t381 * t338 + t473;
t474 = t276 * t362;
t371 = pkin(10) + qJ(5);
t368 = sin(t371);
t472 = t368 * t382;
t369 = cos(t371);
t471 = t369 * t382;
t470 = t374 * t378;
t469 = t374 * t382;
t468 = t375 * t380;
t466 = t376 * t379;
t465 = t376 * t382;
t462 = qJDD(1) - g(3);
t314 = qJDD(2) * pkin(7) + (qJDD(1) * t380 + t383 * t442) * t375;
t440 = qJDD(1) * t377;
t424 = t379 * t440;
t258 = t424 + qJDD(3) * qJ(4) + t314 * t382 + (qJD(4) + t485) * qJD(3);
t404 = pkin(3) * t382 + qJ(4) * t379 + pkin(2);
t265 = t320 * qJD(2) - t404 * qJDD(2) + t412;
t251 = t376 * t258 + t374 * t265;
t405 = pkin(4) * t379 - pkin(8) * t465;
t461 = t405 * qJD(3) + t459;
t460 = (-pkin(7) * t466 - pkin(8) * t469) * qJD(3) + t486;
t453 = qJD(1) * t379;
t360 = t377 * t453;
t305 = t382 * t347 + t360;
t297 = qJD(3) * qJ(4) + t305;
t452 = qJD(1) * t383;
t433 = t375 * t452;
t306 = -t404 * qJD(2) - t433;
t262 = t376 * t297 + t374 * t306;
t458 = -t376 * t436 + t486;
t345 = t413 * qJD(2);
t271 = t374 * t345 + t376 * t485;
t343 = -t381 * t376 + t470;
t398 = t343 * t382;
t457 = qJD(2) * t398 - t343 * qJD(5);
t344 = t374 * t381 + t376 * t378;
t399 = t344 * t382;
t456 = -qJD(2) * t399 + t344 * qJD(5);
t310 = pkin(7) * t465 - t374 * t404;
t372 = t379 ^ 2;
t455 = -t382 ^ 2 + t372;
t451 = qJD(2) * t375;
t446 = qJD(3) * t382;
t445 = qJD(5) * t379;
t444 = qJD(5) * t381;
t441 = qJD(2) * qJD(3);
t439 = qJDD(2) * t379;
t438 = qJDD(3) * t374;
t437 = t382 * qJDD(2);
t366 = t376 * qJDD(3);
t425 = t382 * t441;
t397 = t425 + t439;
t307 = t397 * t374 - t366;
t308 = t397 * t376 + t438;
t435 = -t378 * t307 + t381 * t308 - t338 * t444;
t434 = pkin(4) * t374 + pkin(7);
t431 = t379 * t452;
t430 = t374 * t449;
t429 = t380 * t451;
t428 = t383 * t451;
t427 = qJ(4) * t437;
t423 = t375 * t478;
t422 = t375 * t477;
t250 = -t258 * t374 + t376 * t265;
t396 = t379 * t441 - t437;
t246 = t396 * pkin(4) - pkin(8) * t308 + t250;
t249 = -pkin(8) * t307 + t251;
t417 = t381 * t246 - t378 * t249;
t261 = -t297 * t374 + t376 * t306;
t416 = t381 * t307 + t378 * t308;
t270 = t376 * t345 - t374 * t485;
t323 = t377 * t421 + t418;
t325 = -t377 * t419 + t420;
t414 = g(1) * t325 + g(2) * t323;
t411 = t378 * t246 + t381 * t249;
t255 = -pkin(4) * t449 - pkin(8) * t340 + t261;
t256 = -pkin(8) * t338 + t262;
t247 = t255 * t381 - t256 * t378;
t248 = t255 * t378 + t256 * t381;
t337 = t376 * t404;
t279 = -pkin(8) * t466 - t337 + (-pkin(7) * t374 - pkin(4)) * t382;
t292 = -pkin(8) * t374 * t379 + t310;
t410 = t279 * t381 - t292 * t378;
t409 = t279 * t378 + t292 * t381;
t330 = t377 * t379 + t382 * t468;
t284 = -t330 * t374 - t376 * t467;
t285 = t330 * t376 - t374 * t467;
t408 = t284 * t381 - t285 * t378;
t407 = t284 * t378 + t285 * t381;
t329 = t379 * t468 - t464;
t352 = t480 * t376;
t401 = t405 * qJD(2) + qJD(4) * t374 + qJD(5) * t352 + t270;
t351 = t480 * t374;
t400 = pkin(8) * t430 + qJD(4) * t376 - qJD(5) * t351 - t271;
t252 = -qJD(5) * t473 + t435;
t395 = g(1) * (t325 * t379 - t382 * t422) + g(2) * (t323 * t379 + t382 * t423) + g(3) * t329;
t287 = t323 * t382 - t379 * t423;
t289 = t325 * t382 + t379 * t422;
t394 = g(1) * t289 + g(2) * t287 + g(3) * t330;
t393 = qJD(3) * t360 + t379 * t314 + t347 * t446 - t382 * t440;
t392 = g(3) * t467 - t415;
t391 = -g(3) * t468 - t414;
t259 = t393 + t484;
t390 = -t259 + t395;
t293 = -qJD(3) * pkin(3) + qJD(4) - t485;
t389 = -qJ(4) * t447 + (qJD(4) - t293) * t382;
t253 = -t406 * qJD(5) + t416;
t348 = -t433 - t479;
t387 = -pkin(7) * qJDD(3) + (t348 + t433 - t479) * qJD(3);
t386 = -t393 + t395;
t385 = qJD(2) ^ 2;
t365 = -pkin(4) * t376 - pkin(3);
t346 = t434 * t379;
t342 = qJDD(5) + t396;
t333 = t434 * t446;
t318 = t343 * t379;
t317 = t344 * t379;
t309 = -pkin(7) * t469 - t337;
t291 = t330 * qJD(3) + t379 * t428;
t290 = -t329 * qJD(3) + t382 * t428;
t280 = pkin(4) * t430 + t305;
t274 = qJD(3) * t399 + t444 * t466 - t445 * t470;
t273 = -qJD(3) * t398 - t344 * t445;
t272 = pkin(4) * t338 + t293;
t269 = t290 * t376 + t374 * t429;
t268 = -t290 * t374 + t376 * t429;
t254 = pkin(4) * t307 + t259;
t1 = [t462 * MDP(1) + (-qJD(3) * t291 - qJDD(3) * t329) * MDP(10) + (-qJD(3) * t290 - qJDD(3) * t330) * MDP(11) + (-t284 * t437 + t291 * t338 + t307 * t329) * MDP(12) + (t285 * t437 + t291 * t340 + t308 * t329) * MDP(13) + (-t268 * t340 - t269 * t338 - t284 * t308 - t285 * t307) * MDP(14) + (t250 * t284 + t251 * t285 + t259 * t329 + t261 * t268 + t262 * t269 + t291 * t293 - g(3)) * MDP(15) + (-(-t407 * qJD(5) + t268 * t381 - t269 * t378) * t362 + t408 * t342 + t291 * t276 + t329 * t253) * MDP(21) + (-t291 * t406 + t329 * t252 + (t408 * qJD(5) + t268 * t378 + t269 * t381) * t362 - t407 * t342) * MDP(22) + ((-t268 * t382 + t284 * t447) * MDP(12) + (t269 * t382 - t285 * t447) * MDP(13)) * qJD(2) + ((-(qJDD(2) * MDP(4)) + (-MDP(10) * t382 + MDP(11) * t379 - MDP(3)) * t385) * t380 + (-t396 * MDP(10) - t397 * MDP(11) + qJDD(2) * MDP(3) - t385 * MDP(4)) * t383) * t375; (qJDD(2) * MDP(2)) + (t462 * t467 + t415) * MDP(3) + (-t462 * t468 + t414) * MDP(4) + (qJDD(2) * t372 + 0.2e1 * t379 * t425) * MDP(5) + 0.2e1 * (t379 * t437 - t455 * t441) * MDP(6) + (qJDD(3) * t379 + t382 * t384) * MDP(7) + (qJDD(3) * t382 - t379 * t384) * MDP(8) + (t387 * t379 + t483 * t382) * MDP(10) + (-t483 * t379 + t387 * t382) * MDP(11) + (t391 * t374 + (-t338 * t433 + pkin(7) * t307 + t259 * t374 + (qJD(2) * t309 + t261) * qJD(3)) * t379 + (-t309 * qJDD(2) - t250 + (pkin(7) * t338 + t293 * t374) * qJD(3) - t459 * qJD(2) - t392 * t376) * t382) * MDP(12) + (t391 * t376 + (-t340 * t433 + pkin(7) * t308 + t259 * t376 + (-qJD(2) * t310 - t262) * qJD(3)) * t379 + (t310 * qJDD(2) + t251 + (pkin(7) * t340 + t293 * t376) * qJD(3) + t458 * qJD(2) + t392 * t374) * t382) * MDP(13) + (-t307 * t310 - t308 * t309 - t459 * t340 - t458 * t338 + (-t261 * t376 - t262 * t374) * t446 + (-t250 * t376 - t251 * t374 - t392) * t379) * MDP(14) + (t250 * t309 + t251 * t310 + t458 * t262 + t459 * t261 + t404 * t482 + t404 * t481 + (t259 * t379 + t293 * t446 - t414) * pkin(7) + (-g(3) * pkin(7) * t380 + (-g(3) * t404 - t293 * t453) * t383) * t375) * MDP(15) + (-t252 * t318 - t273 * t406) * MDP(16) + (-t252 * t317 + t253 * t318 - t273 * t276 + t274 * t406) * MDP(17) + (-t252 * t382 - t273 * t362 - t318 * t342 - t406 * t447) * MDP(18) + (t253 * t382 + t274 * t362 - t276 * t447 - t317 * t342) * MDP(19) + (-t342 * t382 - t362 * t447) * MDP(20) + (t410 * t342 - t417 * t382 + t247 * t447 + t333 * t276 + t346 * t253 + t254 * t317 + t272 * t274 - g(1) * (-t324 * t471 + t325 * t368) - g(2) * (-t322 * t471 + t323 * t368) + (t460 * t378 - t461 * t381) * t362 + (t248 * t382 + t409 * t362) * qJD(5) + (-t276 * t431 - g(3) * (t368 * t380 + t369 * t463)) * t375) * MDP(21) + (-t333 * t406 + t346 * t252 - t254 * t318 + t272 * t273 - t409 * t342 + t411 * t382 - t248 * t447 - g(1) * (t324 * t472 + t325 * t369) - g(2) * (t322 * t472 + t323 * t369) + (t461 * t378 + t460 * t381) * t362 + (t247 * t382 + t410 * t362) * qJD(5) + (t406 * t431 - g(3) * (-t368 * t463 + t369 * t380)) * t375) * MDP(22); MDP(7) * t439 + MDP(8) * t437 + qJDD(3) * MDP(9) + (qJD(3) * t305 - t348 * t450 + t386) * MDP(10) + (-t424 + (-qJD(2) * t348 - t314) * t382 + t394) * MDP(11) + (t374 * t427 - pkin(3) * t307 - t305 * t338 + t390 * t376 + (-t261 * t379 + t270 * t382 + t374 * t389) * qJD(2)) * MDP(12) + (t376 * t427 - pkin(3) * t308 - t305 * t340 - t390 * t374 + (t262 * t379 - t271 * t382 + t376 * t389) * qJD(2)) * MDP(13) + (t270 * t340 + t271 * t338 + (-qJ(4) * t307 - qJD(4) * t338 + t261 * t449 + t251) * t376 + (qJ(4) * t308 + qJD(4) * t340 + t262 * t449 - t250) * t374 - t394) * MDP(14) + (-t261 * t270 - t262 * t271 - t293 * t305 + (-t261 * t374 + t262 * t376) * qJD(4) + t390 * pkin(3) + (-t250 * t374 + t251 * t376 - t394) * qJ(4)) * MDP(15) + (t252 * t344 - t406 * t457) * MDP(16) + (-t252 * t343 - t253 * t344 - t276 * t457 + t406 * t456) * MDP(17) + (t342 * t344 - t362 * t457 + t406 * t450) * MDP(18) + (t276 * t450 - t342 * t343 + t362 * t456) * MDP(19) + t362 * MDP(20) * t450 + ((-t351 * t381 - t352 * t378) * t342 + t365 * t253 + t254 * t343 - t247 * t450 - t280 * t276 + (t378 * t400 + t381 * t401) * t362 + t456 * t272 + t395 * t369) * MDP(21) + (t365 * t252 + t254 * t344 - (-t351 * t378 + t352 * t381) * t342 + t280 * t406 + t248 * t450 + (-t378 * t401 + t381 * t400) * t362 + t457 * t272 - t395 * t368) * MDP(22) + (-t379 * t382 * MDP(5) + t455 * MDP(6)) * t385; (t374 * t439 - t366 + (-t340 + t448) * t449) * MDP(12) + (t376 * t439 + t438 + (t338 + t443) * t449) * MDP(13) + (-t338 ^ 2 - t340 ^ 2) * MDP(14) + (t261 * t340 + t262 * t338 - t386 + t484) * MDP(15) + (t253 + t487) * MDP(21) + (t252 + t474) * MDP(22); -t406 * t276 * MDP(16) + (-t276 ^ 2 + t406 ^ 2) * MDP(17) + (t435 - t474) * MDP(18) + (-t416 + t487) * MDP(19) + t342 * MDP(20) + (-t248 * t362 + t272 * t406 - g(1) * (-t289 * t368 + t324 * t369) - g(2) * (-t287 * t368 + t322 * t369) - g(3) * (-t330 * t368 - t369 * t467) + t417) * MDP(21) + (t272 * t276 - t247 * t362 - g(1) * (-t289 * t369 - t324 * t368) - g(2) * (-t287 * t369 - t322 * t368) - g(3) * (-t330 * t369 + t368 * t467) - t411) * MDP(22) + (-MDP(18) * t473 + MDP(19) * t406 - MDP(21) * t248 - MDP(22) * t247) * qJD(5);];
tau = t1;
