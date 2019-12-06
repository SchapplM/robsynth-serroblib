% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:07:04
% EndTime: 2019-12-05 18:07:11
% DurationCPUTime: 3.56s
% Computational Cost: add. (2167->332), mult. (5184->443), div. (0->0), fcn. (3628->10), ass. (0->178)
t365 = sin(pkin(8));
t367 = sin(qJ(4));
t371 = cos(qJ(3));
t419 = qJDD(1) * t371;
t368 = sin(qJ(3));
t420 = qJDD(1) * t368;
t370 = cos(qJ(4));
t483 = t365 * t370;
t318 = t367 * t371 + t368 * t370;
t418 = qJD(3) + qJD(4);
t484 = t418 * t318;
t258 = (qJD(1) * t484 + t367 * t420) * t365 - t419 * t483;
t360 = t365 ^ 2;
t475 = 0.2e1 * t360;
t372 = cos(qJ(1));
t468 = g(3) * t372;
t366 = cos(pkin(8));
t423 = qJDD(1) * qJ(2);
t425 = qJD(1) * qJD(2);
t388 = t423 + t425;
t482 = t366 * t388;
t325 = -pkin(2) * t366 - pkin(6) * t365 - pkin(1);
t315 = t371 * t325;
t459 = t365 * t371;
t417 = pkin(7) * t459;
t466 = qJ(2) * t368;
t281 = -t417 + t315 + (-pkin(3) - t466) * t366;
t465 = qJ(2) * t371;
t344 = t366 * t465;
t442 = t368 * t325 + t344;
t460 = t365 * t368;
t286 = -pkin(7) * t460 + t442;
t447 = t367 * t281 + t370 * t286;
t435 = qJD(1) * t365;
t414 = t368 * t435;
t396 = t367 * t414;
t433 = qJD(1) * t371;
t413 = t365 * t433;
t293 = t370 * t413 - t396;
t287 = t293 * qJ(5);
t306 = qJD(1) * t325 + qJD(2);
t298 = t371 * t306;
t384 = -t366 * t466 - t417;
t277 = qJD(1) * t384 + t298;
t434 = qJD(1) * t366;
t345 = -qJD(3) + t434;
t263 = -pkin(3) * t345 + t277;
t278 = -pkin(7) * t414 + qJD(1) * t344 + t368 * t306;
t270 = t367 * t278;
t406 = t370 * t263 - t270;
t481 = t287 - t406;
t480 = -MDP(13) * t368 - MDP(14) * t371;
t479 = t366 * MDP(4) - MDP(5) * t365;
t369 = sin(qJ(1));
t394 = g(2) * t372 + g(3) * t369;
t478 = qJDD(2) - t394;
t364 = qJ(3) + qJ(4);
t355 = sin(t364);
t356 = cos(t364);
t457 = t366 * t369;
t299 = t355 * t457 + t356 * t372;
t452 = t369 * t356;
t456 = t366 * t372;
t301 = t355 * t456 - t452;
t473 = g(1) * t365;
t477 = -g(2) * t299 + g(3) * t301 + t355 * t473;
t476 = t293 ^ 2;
t361 = t366 ^ 2;
t474 = pkin(7) * t365;
t471 = g(2) * t369;
t469 = qJ(2) * t468;
t385 = qJD(1) * t318;
t290 = t365 * t385;
t464 = qJ(5) * t290;
t463 = qJDD(1) * pkin(1);
t307 = pkin(3) * t414 + qJ(2) * t435;
t275 = pkin(4) * t290 + qJD(5) + t307;
t462 = t275 * t293;
t373 = qJD(1) ^ 2;
t461 = t360 * t373;
t455 = t367 * t368;
t454 = t368 * t369;
t453 = t368 * t372;
t451 = t369 * t371;
t272 = t370 * t278;
t450 = t371 * t372;
t337 = -qJD(4) + t345;
t244 = -pkin(4) * t337 - t481;
t449 = t244 + t481;
t448 = t370 * t277 - t270;
t317 = t370 * t371 - t455;
t446 = (t418 - t434) * t317;
t445 = t366 * t385 - t484;
t430 = qJD(3) * t371;
t432 = qJD(2) * t366;
t444 = t325 * t430 + t371 * t432;
t443 = t418 * t396;
t324 = t371 * pkin(3) + pkin(4) * t356;
t340 = t365 * pkin(3) * t430;
t313 = t365 * qJD(2) + t340;
t346 = pkin(3) * t460;
t316 = t365 * qJ(2) + t346;
t441 = t360 + t361;
t363 = t371 ^ 2;
t440 = t368 ^ 2 - t363;
t439 = MDP(10) * t365;
t438 = MDP(11) * t365;
t431 = qJD(3) * t306;
t429 = qJD(4) * t367;
t428 = qJD(4) * t370;
t421 = qJDD(1) * t366;
t343 = -qJDD(3) + t421;
t331 = -qJDD(4) + t343;
t322 = t331 * MDP(19);
t427 = t343 * MDP(12);
t426 = qJD(3) + t345;
t424 = qJD(1) * qJD(3);
t422 = qJDD(1) * t365;
t415 = qJ(2) * qJD(3) * t366;
t412 = qJ(2) * t421;
t411 = t371 * t425;
t410 = t371 * t424;
t408 = t367 * t419;
t305 = qJDD(1) * t325 + qJDD(2);
t297 = t371 * t305;
t395 = qJD(1) * t415;
t253 = -pkin(3) * t343 + t297 + (-pkin(7) * t422 - t395) * t371 + (-t412 - t431 + (qJD(3) * t474 - t432) * qJD(1)) * t368;
t398 = t368 * t305 + t306 * t430 + t366 * t411 + t371 * t412;
t379 = -t368 * t395 + t398;
t257 = (-t410 - t420) * t474 + t379;
t407 = t370 * t253 - t367 * t257;
t273 = qJD(3) * t384 + t444;
t274 = -t368 * t432 + (-t344 + (-t325 + t474) * t368) * qJD(3);
t405 = -t273 * t367 + t370 * t274;
t404 = -t277 * t367 - t272;
t403 = t370 * t281 - t286 * t367;
t402 = qJD(1) * t426;
t401 = t418 * t371;
t400 = t343 + t421;
t399 = -t367 * t253 - t370 * t257 - t263 * t428 + t278 * t429;
t283 = qJ(2) * t422 + qJD(1) * t340 + qJDD(1) * t346 + t365 * t425;
t397 = 0.2e1 * t441;
t393 = -t468 + t471;
t392 = -t263 * t367 - t272;
t390 = qJD(3) * (t345 + t434);
t389 = t463 - t478;
t387 = (pkin(2) + t324) * t366 - (-qJ(5) - pkin(7) - pkin(6)) * t365 + pkin(1);
t386 = t370 * t273 + t367 * t274 + t281 * t428 - t286 * t429;
t259 = (t408 + (qJD(1) * t401 + t420) * t370) * t365 - t443;
t383 = pkin(4) * t259 + qJDD(5) + t283;
t288 = t290 ^ 2;
t382 = t293 * t290 * MDP(15) + (-t290 * t337 - t258) * MDP(17) + (-t293 * t337 + (-t408 + (-t418 * t433 - t420) * t370) * t365 + t443) * MDP(18) + (-t288 + t476) * MDP(16) - t322;
t380 = t397 * t425 + t471;
t378 = t392 * qJD(4) + t407;
t300 = -t355 * t372 + t366 * t452;
t302 = -t355 * t369 - t356 * t456;
t377 = -g(2) * t300 - g(3) * t302 + t307 * t290 + t356 * t473 + t399;
t375 = -t307 * t293 + t378 + t477;
t351 = pkin(3) * t370 + pkin(4);
t323 = pkin(3) * t368 + pkin(4) * t355;
t311 = -t366 * t450 - t454;
t310 = t366 * t453 - t451;
t309 = t366 * t451 - t453;
t308 = t366 * t454 + t450;
t304 = t317 * t365;
t303 = t318 * t365;
t269 = -t418 * t365 * t455 + t401 * t483;
t268 = t484 * t365;
t256 = -qJ(5) * t303 + t447;
t254 = -pkin(4) * t366 - qJ(5) * t304 + t403;
t248 = -t287 + t448;
t247 = t404 + t464;
t246 = -t392 - t464;
t243 = qJ(5) * t268 - qJD(4) * t447 - qJD(5) * t304 + t405;
t242 = -qJ(5) * t269 - qJD(5) * t303 + t386;
t241 = -qJ(5) * t259 - qJD(5) * t290 - t399;
t240 = -pkin(4) * t331 + qJ(5) * t258 - qJD(5) * t293 + t378;
t1 = [qJDD(1) * MDP(1) + t394 * MDP(2) - t393 * MDP(3) + (t397 * t423 + t380 - t468) * MDP(6) + (pkin(1) * t389 - t469 + (t441 * t423 + t380) * qJ(2)) * MDP(7) + (qJDD(1) * t363 - 0.2e1 * t368 * t410) * t360 * MDP(8) + (-t368 * t419 + t424 * t440) * MDP(9) * t475 + (t368 * t390 - t371 * t400) * t439 + (t368 * t400 + t371 * t390) * t438 + t366 * t427 + (-g(2) * t311 + g(3) * t309 - t297 * t366 - t315 * t343 + (t345 * t366 + (t475 + t361) * qJD(1)) * qJ(2) * t430 + (qJD(3) * t325 * t345 + t388 * t475 + (qJ(2) * t343 + qJD(2) * t345 + t431 + t482) * t366) * t368) * MDP(13) + ((-t368 * t415 + t444) * t345 + t442 * t343 + t379 * t366 - g(2) * t310 - g(3) * t308 + (t411 + (-t368 * t424 + t419) * qJ(2)) * t475) * MDP(14) + (-t258 * t304 - t268 * t293) * MDP(15) + (t258 * t303 - t259 * t304 + t268 * t290 - t269 * t293) * MDP(16) + (t258 * t366 + t268 * t337 - t304 * t331) * MDP(17) + (t259 * t366 + t269 * t337 + t303 * t331) * MDP(18) + t366 * t322 + (-t405 * t337 - t403 * t331 - t407 * t366 + t313 * t290 + t316 * t259 + t283 * t303 + t307 * t269 - g(2) * t302 + g(3) * t300 + (t337 * t447 - t366 * t392) * qJD(4)) * MDP(20) + (-g(2) * t301 - g(3) * t299 - t316 * t258 - t307 * t268 + t283 * t304 + t313 * t293 + t331 * t447 + t337 * t386 - t366 * t399) * MDP(21) + (-t240 * t304 - t241 * t303 - t242 * t290 - t243 * t293 + t244 * t268 - t246 * t269 + t254 * t258 - t256 * t259 + t365 * t394) * MDP(22) + (t241 * t256 + t246 * t242 + t240 * t254 + t244 * t243 + t383 * (pkin(4) * t303 + t316) + t275 * (pkin(4) * t269 + t313) - t469 + (g(2) * t387 - g(3) * t323) * t372 + (-g(2) * (-qJ(2) - t323) + g(3) * t387) * t369) * MDP(23) + t479 * (t389 + t463); t478 * MDP(7) + (-t290 * t435 - t317 * t331) * MDP(20) + (-t293 * t435 + t318 * t331) * MDP(21) + (t258 * t317 - t259 * t318 - t290 * t446 - t293 * t445) * MDP(22) + (t240 * t317 + t241 * t318 + t244 * t445 + t246 * t446 - t275 * t435 - t394) * MDP(23) + (-MDP(13) * t371 + MDP(14) * t368) * t343 + (-MDP(20) * t445 + MDP(21) * t446) * t337 + (-pkin(1) * MDP(7) - t479) * qJDD(1) + (-t361 * MDP(6) + (-MDP(6) + t480) * t360 - t441 * MDP(7) * qJ(2)) * t373 + t480 * t345 ^ 2; t371 * t368 * MDP(8) * t461 - t440 * MDP(9) * t461 + (-t368 * t402 + t419) * t439 + (-t371 * t402 - t420) * t438 - t427 + (-g(2) * t308 + g(3) * t310 + t297 + (-t366 * t402 - t461) * t465 + (-t306 * t426 + t473 - t482) * t368) * MDP(13) + (g(1) * t459 - g(2) * t309 - g(3) * t311 - t298 * t345 + (t426 * t434 + t461) * t466 - t398) * MDP(14) + (t404 * t337 + (-t290 * t413 - t331 * t370 + t337 * t429) * pkin(3) + t375) * MDP(20) + (-t448 * t337 + (-t293 * t413 + t367 * t331 + t337 * t428) * pkin(3) + t377) * MDP(21) + (t258 * t351 + (t246 + t247) * t293 + (-t244 + t248) * t290 + (-t259 * t367 + (-t290 * t370 + t293 * t367) * qJD(4)) * pkin(3)) * MDP(22) + (t240 * t351 - t246 * t248 - t244 * t247 - pkin(4) * t462 + t323 * t473 - g(2) * (t323 * t457 + t324 * t372) - g(3) * (-t323 * t456 + t324 * t369) + (-t275 * t413 + t241 * t367 + (-t244 * t367 + t246 * t370) * qJD(4)) * pkin(3)) * MDP(23) + t382; (t337 * t392 + t375) * MDP(20) + (-t337 * t406 + t377) * MDP(21) + (pkin(4) * t258 - t290 * t449) * MDP(22) + (t449 * t246 + (t240 - t462 + t477) * pkin(4)) * MDP(23) + t382; (-t288 - t476) * MDP(22) + (g(1) * t366 + t244 * t293 + t246 * t290 + t365 * t393 + t383) * MDP(23);];
tau = t1;
