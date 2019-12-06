% Calculate vector of inverse dynamics joint torques for
% S5PRRRR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:10
% EndTime: 2019-12-05 17:13:15
% DurationCPUTime: 3.33s
% Computational Cost: add. (1976->324), mult. (4333->442), div. (0->0), fcn. (3275->14), ass. (0->155)
t363 = sin(qJ(4));
t367 = cos(qJ(4));
t368 = cos(qJ(3));
t424 = qJDD(2) * t368;
t364 = sin(qJ(3));
t425 = qJDD(2) * t364;
t319 = t363 * t368 + t364 * t367;
t356 = qJD(3) + qJD(4);
t464 = t356 * t319;
t262 = qJD(2) * t464 + t363 * t425 - t367 * t424;
t434 = qJD(2) * t368;
t418 = t367 * t434;
t428 = t364 * qJD(2);
t419 = t363 * t428;
t308 = -t418 + t419;
t366 = cos(qJ(5));
t310 = -t363 * t434 - t367 * t428;
t362 = sin(qJ(5));
t449 = t310 * t362;
t271 = -t366 * t308 + t449;
t423 = -qJD(4) - qJD(5);
t351 = qJD(3) - t423;
t466 = t271 * t351;
t396 = t308 * t362 + t366 * t310;
t465 = t351 * t396;
t365 = sin(qJ(2));
t360 = sin(pkin(9));
t361 = cos(pkin(9));
t405 = g(1) * t361 + g(2) * t360;
t393 = t405 * t365;
t369 = cos(qJ(2));
t453 = g(3) * t369;
t384 = t393 - t453;
t426 = qJD(2) * qJD(3);
t416 = t368 * t426;
t463 = t416 + t425;
t436 = qJD(1) * t365;
t333 = qJD(2) * pkin(6) + t436;
t415 = pkin(7) * qJD(2) + t333;
t303 = t415 * t368;
t290 = t367 * t303;
t302 = t415 * t364;
t291 = qJD(3) * pkin(3) - t302;
t399 = -t363 * t291 - t290;
t455 = pkin(8) * t308;
t252 = -t399 - t455;
t349 = -pkin(3) * t368 - pkin(2);
t435 = qJD(1) * t369;
t316 = qJD(2) * t349 - t435;
t277 = pkin(4) * t308 + t316;
t359 = qJ(3) + qJ(4);
t354 = qJ(5) + t359;
t346 = sin(t354);
t347 = cos(t354);
t430 = qJD(5) * t362;
t447 = t361 * t369;
t448 = t360 * t369;
t454 = g(3) * t365;
t381 = -t277 * t271 - g(1) * (-t346 * t360 - t347 * t447) - g(2) * (t346 * t361 - t347 * t448) + t252 * t430 + t347 * t454;
t261 = qJD(4) * t418 - t356 * t419 + t363 * t424 + t463 * t367;
t355 = qJDD(3) + qJDD(4);
t427 = qJD(1) * qJD(2);
t315 = qJDD(2) * pkin(6) + qJDD(1) * t365 + t369 * t427;
t264 = -qJD(3) * t333 * t368 + qJDD(3) * pkin(3) - t463 * pkin(7) - t364 * t315;
t417 = t364 * t426;
t433 = qJD(3) * t364;
t265 = -t333 * t433 + t368 * t315 + (-t417 + t424) * pkin(7);
t379 = qJD(4) * t399 + t367 * t264 - t363 * t265;
t235 = pkin(4) * t355 - pkin(8) * t261 + t379;
t432 = qJD(4) * t363;
t458 = -(qJD(4) * t291 + t265) * t367 - t363 * t264 + t303 * t432;
t236 = -pkin(8) * t262 - t458;
t377 = t277 * t396 - g(1) * (-t346 * t447 + t347 * t360) - g(2) * (-t346 * t448 - t347 * t361) + t366 * t235 - t362 * t236 + t346 * t454;
t350 = qJDD(5) + t355;
t462 = t350 * MDP(23) + t271 * MDP(19) * t396 + (-t271 ^ 2 + t396 ^ 2) * MDP(20);
t318 = t363 * t364 - t367 * t368;
t305 = t318 * t365;
t457 = pkin(6) + pkin(7);
t420 = qJD(3) * t457;
t323 = t364 * t420;
t324 = t368 * t420;
t389 = t319 * t369;
t327 = t457 * t364;
t328 = t457 * t368;
t439 = -t363 * t327 + t367 * t328;
t461 = qJD(1) * t389 - t439 * qJD(4) + t323 * t363 - t367 * t324;
t388 = t318 * t369;
t431 = qJD(4) * t367;
t460 = -qJD(1) * t388 + t367 * t323 + t363 * t324 + t327 * t431 + t328 * t432;
t306 = t310 * pkin(8);
t288 = t363 * t303;
t412 = t367 * t291 - t288;
t251 = t306 + t412;
t391 = pkin(3) * t433 - t436;
t370 = qJD(3) ^ 2;
t406 = -qJDD(1) * t369 + t365 * t427;
t459 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t370 + t365 * (t405 + t427) - t406 - t453;
t414 = t261 * t362 + t366 * t262;
t241 = -qJD(5) * t396 + t414;
t456 = pkin(3) * t363;
t452 = qJD(2) * pkin(2);
t450 = t252 * t366;
t446 = t362 * t350;
t445 = t363 * t366;
t249 = pkin(4) * t356 + t251;
t444 = t366 * t249;
t443 = t366 * t350;
t371 = qJD(2) ^ 2;
t442 = t368 * t371;
t441 = qJDD(1) - g(3);
t440 = -t367 * t302 - t288;
t357 = t364 ^ 2;
t438 = -t368 ^ 2 + t357;
t429 = qJD(5) * t366;
t421 = t366 * t261 - t362 * t262 - t308 * t429;
t411 = t302 * t363 - t290;
t410 = -t367 * t327 - t328 * t363;
t409 = pkin(4) * t464 + t391;
t407 = -qJD(5) * t249 - t236;
t404 = g(1) * t360 - g(2) * t361;
t266 = -pkin(8) * t319 + t410;
t403 = pkin(8) * t464 - qJD(5) * t266 + t460;
t267 = -pkin(8) * t318 + t439;
t278 = t356 * t318;
t402 = -pkin(8) * t278 + qJD(5) * t267 - t461;
t400 = -t249 * t362 - t450;
t304 = t319 * t365;
t398 = -t304 * t366 + t305 * t362;
t397 = -t304 * t362 - t305 * t366;
t275 = t366 * t318 + t319 * t362;
t276 = -t318 * t362 + t319 * t366;
t395 = qJDD(3) * t368 - t364 * t370;
t394 = qJDD(3) * t364 + t368 * t370;
t392 = t405 * t369;
t240 = t310 * t430 + t421;
t280 = pkin(3) * t417 + qJDD(2) * t349 + t406;
t334 = -t435 - t452;
t380 = -pkin(6) * qJDD(3) + (t334 + t435 - t452) * qJD(3);
t376 = -t334 * qJD(2) - t315 + t392 + t454;
t374 = -t310 * t308 * MDP(12) + (t240 - t466) * MDP(21) + (-t241 - t465) * MDP(22) + (t308 * t356 + t261) * MDP(14) + (-t310 * t356 - t262) * MDP(15) + (-t308 ^ 2 + t310 ^ 2) * MDP(13) + t355 * MDP(16) + t462;
t352 = sin(t359);
t353 = cos(t359);
t373 = -g(1) * (-t352 * t360 - t353 * t447) - g(2) * (t352 * t361 - t353 * t448) + t316 * t308 + t353 * t454 + t458;
t372 = -g(1) * (-t352 * t447 + t353 * t360) - g(2) * (-t352 * t448 - t353 * t361) + t316 * t310 + t352 * t454 + t379;
t348 = pkin(3) * t367 + pkin(4);
t295 = pkin(4) * t318 + t349;
t282 = pkin(3) * t428 - pkin(4) * t310;
t256 = t306 + t440;
t255 = t411 + t455;
t254 = -qJD(2) * t389 + t356 * t305;
t253 = -qJD(2) * t388 - t365 * t464;
t246 = pkin(4) * t262 + t280;
t243 = qJD(5) * t276 - t278 * t362 + t366 * t464;
t242 = -qJD(5) * t275 - t278 * t366 - t362 * t464;
t1 = [t441 * MDP(1) + (t254 * t356 - t304 * t355) * MDP(17) + (-t253 * t356 + t305 * t355) * MDP(18) + ((-qJD(5) * t397 - t253 * t362 + t254 * t366) * t351 + t398 * t350) * MDP(24) + (-(qJD(5) * t398 + t253 * t366 + t254 * t362) * t351 - t397 * t350) * MDP(25) + (qJDD(2) * MDP(3) - t371 * MDP(4) + (-0.2e1 * t417 + t424) * MDP(10) + (-0.2e1 * t416 - t425) * MDP(11) - t262 * MDP(17) - t261 * MDP(18) - t241 * MDP(24) - t240 * MDP(25)) * t369 + (-t371 * MDP(3) - qJDD(2) * MDP(4) + (-t394 - t442) * MDP(10) + (t364 * t371 - t395) * MDP(11) + (MDP(17) * t308 - MDP(18) * t310 - MDP(24) * t271 - MDP(25) * t396) * qJD(2)) * t365; qJDD(2) * MDP(2) + (t441 * t369 + t393) * MDP(3) + (-t441 * t365 + t392) * MDP(4) + (qJDD(2) * t357 + 0.2e1 * t364 * t416) * MDP(5) + 0.2e1 * (t364 * t424 - t426 * t438) * MDP(6) + t394 * MDP(7) + t395 * MDP(8) + (t380 * t364 + t459 * t368) * MDP(10) + (-t459 * t364 + t380 * t368) * MDP(11) + (t261 * t319 + t278 * t310) * MDP(12) + (-t261 * t318 - t262 * t319 + t278 * t308 + t310 * t464) * MDP(13) + (-t278 * t356 + t319 * t355) * MDP(14) + (-t318 * t355 - t356 * t464) * MDP(15) + (t349 * t262 + t280 * t318 + t391 * t308 + t316 * t464 + t384 * t353 + t410 * t355 + t461 * t356) * MDP(17) + (t349 * t261 - t316 * t278 + t280 * t319 - t391 * t310 - t352 * t384 - t439 * t355 + t460 * t356) * MDP(18) + (t240 * t276 - t242 * t396) * MDP(19) + (-t240 * t275 - t241 * t276 + t242 * t271 + t243 * t396) * MDP(20) + (t242 * t351 + t276 * t350) * MDP(21) + (-t243 * t351 - t275 * t350) * MDP(22) + ((t266 * t366 - t267 * t362) * t350 + t295 * t241 + t246 * t275 + t277 * t243 + (t362 * t403 - t366 * t402) * t351 - t409 * t271 + t384 * t347) * MDP(24) + (-(t266 * t362 + t267 * t366) * t350 + t295 * t240 + t246 * t276 + t277 * t242 + (t362 * t402 + t366 * t403) * t351 - t409 * t396 - t384 * t346) * MDP(25); (-t411 * t356 + (-t308 * t428 + t355 * t367 - t356 * t432) * pkin(3) + t372) * MDP(17) + (t440 * t356 + (t310 * t428 - t355 * t363 - t356 * t431) * pkin(3) + t373) * MDP(18) + (t282 * t396 + (-t348 * t350 - t235 + (-t423 * t456 + t255) * t351) * t362 + (-t350 * t456 + (-pkin(3) * t431 - qJD(5) * t348 + t256) * t351 + t407) * t366 + t381) * MDP(25) - t364 * MDP(5) * t442 + MDP(8) * t424 + t438 * MDP(6) * t371 + (t364 * t404 + t368 * t376) * MDP(11) + t374 + (t348 * t443 - (t255 * t366 - t256 * t362) * t351 + t282 * t271 + (-t363 * t446 + (-t362 * t367 - t445) * t351 * qJD(4)) * pkin(3) + ((-pkin(3) * t445 - t348 * t362) * t351 + t400) * qJD(5) + t377) * MDP(24) + (t364 * t376 - t368 * t404) * MDP(10) + qJDD(3) * MDP(9) + MDP(7) * t425; (-t356 * t399 + t372) * MDP(17) + (t356 * t412 + t373) * MDP(18) + t374 + (-(-t251 * t362 - t450) * t351 + t400 * qJD(5) + (-t271 * t310 - t351 * t430 + t443) * pkin(4) + t377) * MDP(24) + ((-t252 * t351 - t235) * t362 + (t251 * t351 + t407) * t366 + (-t310 * t396 - t351 * t429 - t446) * pkin(4) + t381) * MDP(25); (t421 - t466) * MDP(21) + (-t414 - t465) * MDP(22) + (-t351 * t400 + t377) * MDP(24) + (-t366 * t236 - t362 * t235 + (-t252 * t362 + t444) * t351 + t381) * MDP(25) + (MDP(21) * t449 + MDP(22) * t396 + MDP(24) * t400 - MDP(25) * t444) * qJD(5) + t462;];
tau = t1;
