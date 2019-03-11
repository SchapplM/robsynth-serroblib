% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:17
% EndTime: 2019-03-09 01:47:22
% DurationCPUTime: 3.43s
% Computational Cost: add. (2690->374), mult. (5086->477), div. (0->0), fcn. (3406->12), ass. (0->180)
t367 = sin(pkin(10));
t373 = sin(qJ(4));
t375 = cos(qJ(4));
t464 = cos(pkin(10));
t311 = t367 * t375 + t464 * t373;
t307 = t311 * qJD(1);
t372 = sin(qJ(6));
t374 = cos(qJ(6));
t289 = -t374 * qJD(4) - t307 * t372;
t418 = t464 * t375;
t327 = qJD(1) * t418;
t444 = qJD(1) * t373;
t305 = t367 * t444 - t327;
t437 = qJD(6) - t305;
t476 = t437 * t289;
t291 = qJD(4) * t372 - t307 * t374;
t475 = t437 * t291;
t434 = qJD(1) * qJD(4);
t420 = t373 * t434;
t283 = qJD(4) * t327 + t311 * qJDD(1) - t367 * t420;
t474 = -qJD(4) * qJD(6) + t283;
t414 = t437 * t374;
t306 = t311 * qJD(4);
t430 = qJDD(1) * t373;
t282 = qJD(1) * t306 - qJDD(1) * t418 + t367 * t430;
t281 = -qJDD(6) + t282;
t456 = t372 * t281;
t473 = t437 * t414 - t456;
t442 = qJD(4) * t373;
t308 = qJD(4) * t418 - t367 * t442;
t353 = t375 * qJDD(3);
t369 = cos(pkin(9));
t435 = qJD(1) * qJD(2);
t331 = t369 * t435;
t376 = -pkin(1) - pkin(2);
t323 = t376 * qJDD(1) + qJDD(2);
t368 = sin(pkin(9));
t431 = qJDD(1) * t369;
t451 = qJ(2) * t431 + t368 * t323;
t288 = t331 + t451;
t285 = -qJDD(1) * pkin(7) + t288;
t383 = qJ(5) * qJDD(1) + qJD(1) * qJD(5) - qJD(3) * qJD(4) - t285;
t324 = t376 * qJD(1) + qJD(2);
t445 = qJD(1) * t369;
t303 = qJ(2) * t445 + t368 * t324;
t298 = -qJD(1) * pkin(7) + t303;
t413 = qJ(5) * qJD(1) - t298;
t395 = t413 * qJD(4);
t249 = qJDD(4) * pkin(4) + t383 * t373 + t375 * t395 + t353;
t250 = (qJDD(3) + t395) * t373 - t383 * t375;
t241 = t464 * t249 - t367 * t250;
t239 = -qJDD(4) * pkin(5) - t241;
t279 = qJ(5) * t444 + t375 * qJD(3) - t298 * t373;
t273 = qJD(4) * pkin(4) + t279;
t280 = qJD(3) * t373 - t413 * t375;
t457 = t367 * t280;
t253 = t464 * t273 - t457;
t251 = -qJD(4) * pkin(5) - t253;
t443 = qJD(2) * t369;
t407 = -qJD(5) + t443;
t317 = t369 * qJ(2) + t368 * t376;
t314 = -pkin(7) + t317;
t455 = qJ(5) - t314;
t415 = qJD(4) * t455;
t278 = t373 * t415 + t407 * t375;
t381 = -t407 * t373 + t375 * t415;
t257 = t464 * t278 + t367 * t381;
t299 = t455 * t375;
t417 = t455 * t373;
t263 = -t464 * t299 + t367 * t417;
t389 = -t367 * t373 + t418;
t316 = -t368 * qJ(2) + t369 * t376;
t313 = pkin(3) - t316;
t358 = t375 * pkin(4);
t396 = t313 + t358;
t265 = pkin(5) * t389 + pkin(8) * t311 + t396;
t242 = t367 * t249 + t464 * t250;
t446 = qJD(1) * t368;
t302 = -qJ(2) * t446 + t324 * t369;
t297 = qJD(1) * pkin(3) - t302;
t286 = qJD(1) * t358 + qJD(5) + t297;
t261 = -pkin(5) * t305 + pkin(8) * t307 + t286;
t412 = qJDD(4) * pkin(8) + qJD(6) * t261 + t242;
t468 = sin(qJ(1));
t469 = cos(qJ(1));
t310 = -t468 * t368 - t469 * t369;
t467 = g(1) * t310;
t471 = -t239 * t311 - t251 * t308 + t263 * t281 - (qJD(6) * t265 + t257) * t437 - t412 * t389 - t467;
t338 = pkin(4) * t367 + pkin(8);
t363 = qJ(4) + pkin(10);
t347 = sin(t363);
t348 = cos(t363);
t312 = t469 * t368 - t468 * t369;
t402 = g(2) * t312 + t467;
t470 = (-pkin(4) * t444 - pkin(5) * t307 - pkin(8) * t305 + qJD(6) * t338) * t437 + t402 * t347 - g(3) * t348 + t239;
t466 = g(2) * t310;
t465 = g(3) * t375;
t463 = pkin(1) * qJDD(1);
t423 = t372 * qJDD(4) - t374 * t474;
t438 = qJD(6) * t372;
t259 = t307 * t438 + t423;
t462 = t259 * t372;
t461 = t289 * t307;
t460 = t291 * t307;
t459 = t310 * t348;
t458 = t312 * t348;
t272 = t374 * t281;
t454 = qJDD(3) + g(3);
t270 = t464 * t280;
t254 = t367 * t273 + t270;
t453 = -t369 * t307 + t308 * t368;
t452 = -t368 * t306 - t389 * t445;
t450 = t469 * pkin(1) + t468 * qJ(2);
t449 = g(1) * t468 - g(2) * t469;
t364 = t373 ^ 2;
t448 = -t375 ^ 2 + t364;
t377 = qJD(4) ^ 2;
t378 = qJD(1) ^ 2;
t447 = t377 + t378;
t252 = qJD(4) * pkin(8) + t254;
t441 = qJD(6) * t252;
t440 = qJD(6) * t307;
t439 = qJD(6) * t311;
t436 = qJ(2) * qJDD(1);
t432 = qJDD(1) * t368;
t429 = qJDD(1) * t375;
t428 = qJDD(4) * t373;
t427 = qJDD(4) * t375;
t426 = 0.2e1 * t435;
t425 = t311 * t456;
t424 = t311 * t272;
t422 = t469 * pkin(2) + t450;
t329 = t368 * t435;
t416 = -qJ(2) * t432 + t323 * t369;
t409 = 0.2e1 * t375 * t434;
t408 = qJDD(2) - t463;
t406 = -t468 * pkin(1) + t469 * qJ(2);
t405 = -pkin(4) * t442 + t368 * qJD(2);
t403 = g(1) * t312 - t466;
t287 = -t329 + t416;
t401 = -qJD(6) * t369 + t452;
t400 = -t441 + t466;
t399 = t259 * t389 - t291 * t306;
t352 = t374 * qJDD(4);
t260 = t291 * qJD(6) - t283 * t372 - t352;
t398 = -t260 * t389 + t289 * t306;
t397 = t302 * t368 - t303 * t369;
t301 = t389 * t368;
t394 = qJD(6) * t301 + t446;
t284 = qJDD(1) * pkin(3) - t287;
t393 = -t272 + (t305 * t372 - t438) * t437;
t392 = -g(1) * t458 - t265 * t281;
t391 = t308 * t372 + t374 * t439;
t390 = -t308 * t374 + t311 * t438;
t388 = g(1) * t469 + g(2) * t468;
t387 = -t468 * pkin(2) + t406;
t386 = qJD(1) * t297 - t285 - t402;
t385 = -g(2) * t458 - g(3) * t347 - t412;
t258 = t464 * t279 - t457;
t384 = t338 * t281 + (t251 + t258) * t437;
t380 = -qJDD(4) * t314 + (-qJD(1) * t313 - t297 - t443) * qJD(4);
t266 = qJDD(5) + t284 + (-t420 + t429) * pkin(4);
t379 = qJDD(1) * t313 - t314 * t377 + t284 + t329 - t403;
t371 = -qJ(5) - pkin(7);
t343 = t358 + pkin(3);
t339 = -t464 * pkin(4) - pkin(5);
t321 = -t373 * t377 + t427;
t320 = -t375 * t377 - t428;
t300 = t311 * t368;
t275 = t312 * t372 - t374 * t459;
t274 = t312 * t374 + t372 * t459;
t264 = -pkin(5) * t306 + pkin(8) * t308 + t405;
t262 = -t299 * t367 - t464 * t417;
t256 = t279 * t367 + t270;
t255 = t278 * t367 - t464 * t381;
t246 = -pkin(5) * t282 + pkin(8) * t283 + t266;
t245 = t374 * t246;
t244 = t252 * t374 + t261 * t372;
t243 = -t252 * t372 + t261 * t374;
t1 = [qJDD(1) * MDP(1) + (-qJDD(2) + t449 + 0.2e1 * t463) * MDP(4) + (-t388 + t426 + 0.2e1 * t436) * MDP(5) + (-t408 * pkin(1) - g(1) * t406 - g(2) * t450 + (t426 + t436) * qJ(2)) * MDP(6) + (-qJDD(1) * t316 + 0.2e1 * t329 - t403 - t416) * MDP(7) + (qJDD(1) * t317 + 0.2e1 * t331 + t402 + t451) * MDP(8) + (-g(1) * t387 - g(2) * t422 - t397 * qJD(2) + t287 * t316 + t288 * t317) * MDP(9) + (qJDD(1) * t364 + t373 * t409) * MDP(10) + 0.2e1 * (t373 * t429 - t448 * t434) * MDP(11) + t320 * MDP(12) - t321 * MDP(13) + (t380 * t373 + t379 * t375) * MDP(15) + (-t379 * t373 + t380 * t375) * MDP(16) + (t241 * t311 - t242 * t389 + t253 * t308 + t254 * t306 - t255 * t307 + t257 * t305 - t262 * t283 + t263 * t282 - t402) * MDP(17) + (t242 * t263 + t254 * t257 - t241 * t262 - t253 * t255 + t266 * t396 + t286 * t405 - g(1) * (-t310 * t371 + t312 * t343 + t387) - g(2) * (-t310 * t343 - t312 * t371 + t422)) * MDP(18) + (-t259 * t311 * t374 + t390 * t291) * MDP(19) + ((t289 * t374 + t291 * t372) * t308 + (t462 + t260 * t374 + (-t289 * t372 + t291 * t374) * qJD(6)) * t311) * MDP(20) + (t390 * t437 + t399 + t424) * MDP(21) + (t391 * t437 + t398 - t425) * MDP(22) + (-t281 * t389 - t306 * t437) * MDP(23) + (-g(2) * t275 - t243 * t306 + t245 * t389 + t255 * t289 + t262 * t260 + (t264 * t437 + (-t251 * t311 - t252 * t389 - t263 * t437) * qJD(6) + t392) * t374 + t471 * t372) * MDP(24) + (-g(2) * t274 + t244 * t306 + t255 * t291 + t262 * t259 + (-(-qJD(6) * t263 + t264) * t437 - (t246 - t441) * t389 + t251 * t439 - t392) * t372 + t471 * t374) * MDP(25) + t449 * MDP(2) + t388 * MDP(3); -qJDD(1) * MDP(4) - t378 * MDP(5) + (-qJ(2) * t378 + t408 - t449) * MDP(6) + (-t368 * t378 - t431) * MDP(7) + (-t369 * t378 + t432) * MDP(8) + (t397 * qJD(1) + t287 * t369 + t288 * t368 - t449) * MDP(9) + ((0.2e1 * t420 - t429) * t369 + (-t447 * t375 - t428) * t368) * MDP(15) + ((t409 + t430) * t369 + (t447 * t373 - t427) * t368) * MDP(16) + (t282 * t301 - t283 * t300 + t452 * t305 - t453 * t307) * MDP(17) + (-t241 * t300 + t242 * t301 - t453 * t253 + t452 * t254 - t266 * t369 - t286 * t446 - t449) * MDP(18) + (-(-t301 * t372 - t369 * t374) * t281 + t300 * t260 - (t401 * t372 + t394 * t374) * t437 + t453 * t289) * MDP(24) + ((t301 * t374 - t369 * t372) * t281 + t300 * t259 - (-t394 * t372 + t401 * t374) * t437 + t453 * t291) * MDP(25); t454 * MDP(9) + t321 * MDP(15) + t320 * MDP(16) + (t282 * t311 + t283 * t389 + t305 * t308 - t306 * t307) * MDP(17) + (t241 * t389 + t242 * t311 - t253 * t306 + t254 * t308 + g(3)) * MDP(18) + (t398 + t425) * MDP(24) + (-t399 + t424) * MDP(25) - (t391 * MDP(24) - t390 * MDP(25)) * t437; -MDP(12) * t430 - MDP(13) * t429 + qJDD(4) * MDP(14) + (t386 * t373 + t353 + t465) * MDP(15) + (-t454 * t373 + t386 * t375) * MDP(16) + ((-t254 + t256) * t307 + (t253 - t258) * t305 + (t282 * t367 + t464 * t283) * pkin(4)) * MDP(17) + (t253 * t256 - t254 * t258 + (t464 * t241 + t465 + t242 * t367 + (qJD(1) * t286 - t402) * t373) * pkin(4)) * MDP(18) + (t291 * t414 + t462) * MDP(19) + ((t259 - t476) * t374 + (-t260 - t475) * t372) * MDP(20) + (t460 + t473) * MDP(21) + (t393 - t461) * MDP(22) + t437 * t307 * MDP(23) + (t243 * t307 - t256 * t289 + t339 * t260 + t384 * t372 - t470 * t374) * MDP(24) + (-t244 * t307 - t256 * t291 + t339 * t259 + t470 * t372 + t384 * t374) * MDP(25) + (-t373 * t375 * MDP(10) + t448 * MDP(11)) * t378; (-t305 ^ 2 - t307 ^ 2) * MDP(17) + (-t253 * t307 - t254 * t305 + t266 - t403) * MDP(18) + (t393 + t461) * MDP(24) + (t460 - t473) * MDP(25); t291 * t289 * MDP(19) + (-t289 ^ 2 + t291 ^ 2) * MDP(20) + (t423 + t476) * MDP(21) + (t352 + t475) * MDP(22) - t281 * MDP(23) + (-g(1) * t274 + t244 * t437 - t251 * t291 + t245) * MDP(24) + (g(1) * t275 + t243 * t437 + t251 * t289) * MDP(25) + (MDP(22) * t440 + t400 * MDP(24) + t385 * MDP(25)) * t374 + (MDP(21) * t440 + t474 * MDP(22) + t385 * MDP(24) + (-t246 - t400) * MDP(25)) * t372;];
tau  = t1;
