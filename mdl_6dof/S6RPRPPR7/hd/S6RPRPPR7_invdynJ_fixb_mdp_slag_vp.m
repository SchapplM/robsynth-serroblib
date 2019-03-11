% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:32
% EndTime: 2019-03-09 02:57:40
% DurationCPUTime: 4.26s
% Computational Cost: add. (2693->419), mult. (5209->506), div. (0->0), fcn. (3457->10), ass. (0->187)
t385 = sin(qJ(3));
t457 = qJD(1) * t385;
t337 = pkin(3) * t457 + qJD(1) * qJ(2) + qJD(4);
t386 = sin(qJ(1));
t389 = cos(qJ(1));
t461 = g(1) * t386 - g(2) * t389;
t506 = -qJD(1) * t337 - t461;
t392 = qJD(1) ^ 2;
t405 = -qJ(2) * t392 - t461;
t505 = MDP(14) + MDP(16);
t489 = sin(pkin(9));
t432 = t489 * t385;
t347 = qJD(1) * t432;
t388 = cos(qJ(3));
t490 = cos(pkin(9));
t433 = t490 * t388;
t423 = qJD(1) * t433;
t330 = t423 - t347;
t320 = qJD(6) + t330;
t504 = t320 ^ 2;
t390 = -pkin(1) - pkin(7);
t343 = t390 * qJD(1) + qJD(2);
t318 = -qJ(4) * t457 + t343 * t385;
t309 = t489 * t318;
t456 = qJD(1) * t388;
t319 = -qJ(4) * t456 + t388 * t343;
t293 = t490 * t319 - t309;
t448 = -qJD(5) + t293;
t378 = qJDD(1) * qJ(2);
t422 = g(1) * t389 + g(2) * t386;
t379 = qJD(1) * qJD(2);
t440 = 0.2e1 * t379;
t502 = 0.2e1 * t378 + t440 - t422;
t312 = qJD(3) * pkin(3) + t319;
t434 = t490 * t318;
t285 = t489 * t312 + t434;
t280 = -qJD(3) * qJ(5) - t285;
t401 = t490 * t385 + t489 * t388;
t326 = t401 * qJD(1);
t492 = t326 * pkin(5);
t269 = -t280 - t492;
t292 = t489 * t319 + t434;
t427 = qJDD(1) * t489;
t428 = qJDD(1) * t490;
t410 = t385 * t427 - t388 * t428;
t446 = qJD(1) * qJD(3);
t301 = t401 * t446 + t410;
t298 = -qJDD(6) + t301;
t359 = -t490 * pkin(3) - pkin(4);
t351 = -pkin(8) + t359;
t501 = -t351 * t298 + (t269 - t292 + t492) * t320;
t342 = t390 * qJDD(1) + qJDD(2);
t335 = t388 * t342;
t442 = qJDD(1) * t388;
t445 = qJD(1) * qJD(4);
t454 = qJD(3) * t385;
t283 = -t388 * t445 - t343 * t454 + qJDD(3) * pkin(3) + t335 + (t385 * t446 - t442) * qJ(4);
t453 = qJD(3) * t388;
t289 = (-qJ(4) * qJD(1) + t343) * t453 + (-qJ(4) * qJDD(1) + t342 - t445) * t385;
t263 = t490 * t283 - t489 * t289;
t264 = t489 * t283 + t490 * t289;
t284 = t490 * t312 - t309;
t430 = qJD(3) * t489;
t431 = qJD(3) * t490;
t328 = t385 * t430 - t388 * t431;
t329 = -t385 * t431 - t388 * t430;
t333 = -t432 + t433;
t500 = -t263 * t333 - t264 * t401 - t284 * t329 + t285 * t328;
t439 = qJDD(3) * qJ(5) + t264;
t261 = -qJD(3) * qJD(5) - t439;
t409 = qJDD(5) - t263;
t485 = qJDD(3) * pkin(4);
t262 = t409 - t485;
t412 = qJD(5) - t284;
t278 = -qJD(3) * pkin(4) + t412;
t499 = t261 * t401 + t262 * t333 + t278 * t329 - t280 * t328;
t325 = t330 ^ 2;
t496 = pkin(4) + pkin(8);
t437 = -qJD(3) * t423 - t385 * t428 - t388 * t427;
t300 = qJD(3) * t347 + t437;
t495 = pkin(4) * t300;
t376 = qJ(3) + pkin(9);
t366 = sin(t376);
t494 = g(3) * t366;
t493 = g(3) * t385;
t491 = t330 * pkin(5);
t372 = t385 * pkin(3);
t488 = pkin(1) * qJDD(1);
t486 = qJ(5) * t301;
t384 = sin(qJ(6));
t387 = cos(qJ(6));
t450 = qJD(6) * t387;
t438 = t387 * qJDD(3) - t384 * t300 + t326 * t450;
t444 = qJD(3) * qJD(6);
t267 = -t384 * t444 + t438;
t480 = t267 * t387;
t464 = qJ(2) + t372;
t420 = -qJ(5) * t333 + t464;
t279 = t401 * t496 + t420;
t478 = t279 * t298;
t474 = t298 * t384;
t304 = qJD(3) * t384 - t387 * t326;
t473 = t304 * t320;
t472 = t304 * t326;
t306 = qJD(3) * t387 + t326 * t384;
t471 = t306 * t320;
t470 = t306 * t326;
t469 = t401 * t384;
t468 = t384 * t386;
t467 = t384 * t389;
t466 = t386 * t387;
t295 = t387 * t298;
t465 = t387 * t389;
t463 = qJ(4) - t390;
t462 = t389 * pkin(1) + t386 * qJ(2);
t382 = t388 ^ 2;
t460 = t385 ^ 2 - t382;
t391 = qJD(3) ^ 2;
t459 = -t391 - t392;
t413 = -qJ(5) * t330 + t337;
t270 = t496 * t326 + t413;
t452 = qJD(6) * t270;
t451 = qJD(6) * t384;
t449 = pkin(3) * t453 + qJD(2);
t447 = t491 - t448;
t443 = qJDD(1) * t385;
t441 = qJDD(3) * t385;
t436 = t388 * t446;
t435 = -pkin(1) * t386 + t389 * qJ(2);
t339 = t463 * t388;
t429 = pkin(3) * t456 + qJ(5) * t326;
t426 = t320 * t384;
t416 = qJDD(4) + t378 + t379 + (t436 + t443) * pkin(3);
t397 = -qJD(5) * t330 + t416 + t486;
t254 = -t496 * t300 + t397;
t266 = -t496 * qJD(3) + t412 + t491;
t425 = qJD(6) * t266 + t254;
t424 = qJDD(2) - t488;
t367 = cos(t376);
t419 = pkin(4) * t366 - qJ(5) * t367;
t315 = -qJD(4) * t388 + t463 * t454;
t316 = -qJD(3) * t339 - qJD(4) * t385;
t287 = -t490 * t315 + t489 * t316;
t338 = t463 * t385;
t302 = -t489 * t338 + t490 * t339;
t418 = -t452 - t494;
t259 = t266 * t384 + t270 * t387;
t383 = -qJ(4) - pkin(7);
t415 = t389 * t372 + t386 * t383 + t435;
t414 = t386 * t372 - t383 * t389 + t462;
t411 = -t320 * t426 - t295;
t408 = -t328 * t384 + t401 * t450;
t407 = 0.2e1 * qJ(2) * t446 + qJDD(3) * t390;
t404 = -qJ(5) * t329 - qJD(5) * t333 + t449;
t257 = pkin(5) * t300 - t261;
t290 = t333 * pkin(5) + t302;
t403 = t257 * t401 - t269 * t328 + t290 * t298;
t402 = -t387 * t504 + t474;
t288 = t489 * t315 + t490 * t316;
t303 = -t490 * t338 - t489 * t339;
t400 = -g(3) * t367 - t366 * t461;
t399 = t416 - t422;
t396 = t287 * t330 - t288 * t326 + t300 * t303 - t301 * t302 + t461;
t395 = -t390 * t391 + t502;
t286 = pkin(4) * t326 + t413;
t394 = t286 * t330 + t461 * t367 + t409 - t494;
t393 = t257 + (-qJD(6) * t351 + t496 * t330 + t429) * t320 + t400;
t369 = qJDD(3) * t388;
t355 = t489 * pkin(3) + qJ(5);
t324 = -t367 * t468 + t465;
t323 = -t367 * t466 - t467;
t322 = -t367 * t467 - t466;
t321 = -t367 * t465 + t468;
t299 = pkin(4) * t401 + t420;
t297 = t387 * t300;
t294 = pkin(4) * t330 + t429;
t291 = -pkin(5) * t401 + t303;
t276 = -pkin(4) * t328 + t404;
t272 = t328 * pkin(5) + t288;
t271 = t329 * pkin(5) + t287;
t268 = t306 * qJD(6) + qJDD(3) * t384 + t297;
t265 = -t496 * t328 + t404;
t260 = t397 - t495;
t258 = t266 * t387 - t270 * t384;
t256 = -t301 * pkin(5) - t496 * qJDD(3) + t409;
t255 = t387 * t256;
t1 = [(t395 * t385 + t407 * t388) * MDP(12) + (-t407 * t385 + t395 * t388) * MDP(13) + (-t401 * t295 - t268 * t333 - t304 * t329 + (-t328 * t387 - t401 * t451) * t320) * MDP(23) + (qJD(3) * t287 + qJDD(3) * t302 - t260 * t401 - t276 * t326 + t286 * t328 + t299 * t300 + t422 * t366) * MDP(17) + (-g(1) * t321 - g(2) * t323 - t259 * t329 + t291 * t267 + t272 * t306 + (-(qJD(6) * t290 + t265) * t320 + t478 - t425 * t333 + t269 * qJD(6) * t401) * t387 + (-(-qJD(6) * t279 + t271) * t320 - (t256 - t452) * t333 + t403) * t384) * MDP(26) + ((t304 * t384 - t306 * t387) * t328 - (-t480 + t268 * t384 + (t304 * t387 + t306 * t384) * qJD(6)) * t401) * MDP(21) + (qJDD(1) * t382 - 0.2e1 * t385 * t436) * MDP(7) + (t260 * t299 + t286 * t276 - t261 * t303 - t280 * t288 + t262 * t302 + t278 * t287 - g(1) * (t419 * t389 + t415) - g(2) * (t419 * t386 + t414)) * MDP(19) + 0.2e1 * (-t385 * t442 + t460 * t446) * MDP(8) + t461 * MDP(2) + (-t424 * pkin(1) - g(1) * t435 - g(2) * t462 + (t440 + t378) * qJ(2)) * MDP(6) + qJDD(1) * MDP(1) + (-g(1) * t415 - g(2) * t414 - t263 * t302 + t264 * t303 - t284 * t287 + t285 * t288 + t337 * t449 + t416 * t464) * MDP(15) + (-t388 * t391 - t441) * MDP(10) + (qJD(3) * t288 + qJDD(3) * t303 - t260 * t333 - t276 * t330 - t286 * t329 + t299 * t301 + t422 * t367) * MDP(18) + t422 * MDP(3) + (qJDD(2) - t461 - 0.2e1 * t488) * MDP(4) + (t267 * t469 + t408 * t306) * MDP(20) + (t267 * t333 - t298 * t469 + t306 * t329 + t408 * t320) * MDP(22) + (-g(1) * t322 - g(2) * t324 + t255 * t333 + t258 * t329 + t291 * t268 + t272 * t304 + (-t254 * t333 - t265 * t320 + t478) * t384 + (t271 * t320 - t403) * t387 + ((-t279 * t387 - t290 * t384) * t320 - t259 * t333 + t269 * t469) * qJD(6)) * MDP(25) + t502 * MDP(5) + (-t385 * t391 + t369) * MDP(9) + (t396 + t499) * MDP(16) + (t396 + t500) * MDP(14) + (-t298 * t333 + t320 * t329) * MDP(24); qJDD(1) * MDP(4) - t392 * MDP(5) + (t424 + t405) * MDP(6) + (t459 * t385 + t369) * MDP(12) + (t459 * t388 - t441) * MDP(13) + (-t500 + t506) * MDP(15) + (qJD(1) * t326 - qJD(3) * t329 - qJDD(3) * t333) * MDP(17) + (qJD(1) * t330 - qJD(3) * t328 + qJDD(3) * t401) * MDP(18) + (-qJD(1) * t286 - t461 - t499) * MDP(19) + (t268 * t401 + t333 * t295 - t304 * t328) * MDP(25) + (t267 * t401 - t306 * t328 - t333 * t474) * MDP(26) + ((qJD(1) * t384 - t329 * t387 + t333 * t451) * MDP(25) + (qJD(1) * t387 + t329 * t384 + t333 * t450) * MDP(26)) * t320 + t505 * (t300 * t401 + t301 * t333 + t326 * t328 - t329 * t330); MDP(9) * t442 - MDP(10) * t443 + qJDD(3) * MDP(11) + (t405 * t388 + t335 + t493) * MDP(12) + (g(3) * t388 + (-t342 - t405) * t385) * MDP(13) + ((t285 - t292) * t330 + (-t284 + t293) * t326 + (t489 * t300 + t490 * t301) * pkin(3)) * MDP(14) + (t284 * t292 - t285 * t293 + (t490 * t263 + t489 * t264 + t506 * t388 + t493) * pkin(3)) * MDP(15) + (t300 * t355 - t301 * t359 + (-t280 - t292) * t330 + (t278 + t448) * t326) * MDP(16) + (-t292 * qJD(3) + t294 * t326 + (-pkin(4) + t359) * qJDD(3) + t394) * MDP(17) + (qJDD(3) * t355 - t286 * t326 + t294 * t330 + (0.2e1 * qJD(5) - t293) * qJD(3) + t400 + t439) * MDP(18) + (-t261 * t355 + t262 * t359 - t286 * t294 - t278 * t292 - g(3) * (-t372 - t419) + t448 * t280 - t461 * (pkin(3) * t388 + pkin(4) * t367 + qJ(5) * t366)) * MDP(19) + (-t306 * t426 + t480) * MDP(20) + ((-t268 - t471) * t387 + (-t267 + t473) * t384) * MDP(21) + (t411 + t470) * MDP(22) + (t402 - t472) * MDP(23) + t320 * t326 * MDP(24) + (t258 * t326 + t355 * t268 + t447 * t304 + t393 * t384 + t501 * t387) * MDP(25) + (-t259 * t326 + t355 * t267 + t447 * t306 - t384 * t501 + t393 * t387) * MDP(26) + (MDP(7) * t385 * t388 - MDP(8) * t460) * t392; (t284 * t330 + t285 * t326 + t399) * MDP(15) + ((t347 - t330) * qJD(3) + t437) * MDP(17) + (qJD(3) * t326 + t301) * MDP(18) + (-t495 + t486 - t280 * t326 + (-qJD(5) - t278) * t330 + t399) * MDP(19) + (t402 + t472) * MDP(25) + (t384 * t504 + t295 + t470) * MDP(26) + t505 * (-t326 ^ 2 - t325); -t410 * MDP(16) + (-t326 * t330 + qJDD(3)) * MDP(17) + (-t325 - t391) * MDP(18) + (t280 * qJD(3) + t394 - t485) * MDP(19) + (-qJD(3) * t304 + t411) * MDP(25) + (-qJD(3) * t306 + t402) * MDP(26); t306 * t304 * MDP(20) + (-t304 ^ 2 + t306 ^ 2) * MDP(21) + (t438 + t473) * MDP(22) + (-t297 + t471) * MDP(23) - t298 * MDP(24) + (-g(1) * t323 + g(2) * t321 + t259 * t320 - t269 * t306 + t255) * MDP(25) + (g(1) * t324 - g(2) * t322 + t258 * t320 + t269 * t304) * MDP(26) + (-MDP(23) * t444 + t418 * MDP(25) - t425 * MDP(26)) * t387 + (-MDP(22) * t444 + (-qJD(6) * t326 - qJDD(3)) * MDP(23) - t425 * MDP(25) + (-t256 - t418) * MDP(26)) * t384;];
tau  = t1;
