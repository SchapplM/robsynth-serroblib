% Calculate vector of inverse dynamics joint torques for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:31:00
% EndTime: 2019-12-31 17:31:07
% DurationCPUTime: 4.21s
% Computational Cost: add. (2544->401), mult. (6512->573), div. (0->0), fcn. (5110->10), ass. (0->177)
t370 = cos(qJ(2));
t362 = sin(pkin(4));
t366 = sin(qJ(2));
t429 = t366 * qJDD(1);
t414 = t362 * t429;
t432 = qJD(1) * qJD(2);
t415 = t362 * t432;
t482 = -t370 * t415 - t414;
t442 = qJD(1) * t362;
t420 = t370 * t442;
t343 = -qJD(3) + t420;
t369 = cos(qJ(3));
t365 = sin(qJ(3));
t421 = t366 * t442;
t402 = t365 * t421;
t363 = cos(pkin(4));
t441 = qJD(1) * t363;
t407 = qJD(2) + t441;
t304 = -t369 * t407 + t402;
t299 = qJD(4) + t304;
t367 = sin(qJ(1));
t371 = cos(qJ(1));
t453 = t366 * t371;
t324 = t363 * t453 + t367 * t370;
t460 = t362 * t365;
t288 = t324 * t369 - t371 * t460;
t448 = t370 * t371;
t454 = t366 * t367;
t323 = -t363 * t448 + t454;
t364 = sin(qJ(4));
t368 = cos(qJ(4));
t481 = t288 * t364 - t323 * t368;
t480 = t288 * t368 + t323 * t364;
t431 = qJDD(1) * t363;
t354 = qJDD(2) + t431;
t388 = qJD(3) * t407;
t439 = qJD(2) * t370;
t418 = t365 * t439;
t435 = qJD(3) * t369;
t264 = -t369 * t354 + ((t366 * t435 + t418) * qJD(1) + t365 * t429) * t362 + t365 * t388;
t457 = t362 * t370;
t474 = pkin(1) * t366;
t444 = pkin(6) * t457 + t363 * t474;
t312 = pkin(7) * t363 + t444;
t389 = -pkin(2) * t370 - pkin(7) * t366 - pkin(1);
t313 = t389 * t362;
t479 = t369 * t312 + t365 * t313;
t430 = qJDD(1) * t370;
t353 = t362 * t430;
t401 = t366 * t415;
t314 = qJDD(3) - t353 + t401;
t427 = pkin(1) * t441;
t318 = pkin(6) * t420 + t366 * t427;
t294 = t407 * pkin(7) + t318;
t298 = qJD(1) * t313;
t262 = t294 * t369 + t298 * t365;
t404 = qJD(2) * t427;
t426 = pkin(1) * t431;
t422 = -pkin(6) * t353 - t366 * t426 - t370 * t404;
t379 = -pkin(6) * t401 - t422;
t275 = pkin(7) * t354 + t379;
t399 = pkin(2) * t366 - pkin(7) * t370;
t387 = t399 * qJD(2);
t278 = (qJD(1) * t387 + t389 * qJDD(1)) * t362;
t375 = -t262 * qJD(3) - t275 * t365 + t369 * t278;
t244 = -pkin(3) * t314 - t375;
t306 = t365 * t407 + t369 * t421;
t326 = -t363 * t454 + t448;
t458 = t362 * t369;
t290 = -t326 * t365 + t367 * t458;
t459 = t362 * t366;
t321 = -t363 * t369 + t365 * t459;
t450 = t369 * t371;
t464 = t324 * t365;
t384 = g(1) * t290 + g(2) * (-t362 * t450 - t464) - g(3) * t321;
t478 = (pkin(3) * t306 + t299 * pkin(8)) * t299 + t244 + t384;
t260 = qJDD(4) + t264;
t339 = -pkin(3) * t369 - pkin(8) * t365 - pkin(2);
t477 = (t318 + t343 * (pkin(3) * t365 - pkin(8) * t369)) * t299 - t339 * t260;
t317 = t362 * t387;
t355 = pkin(6) * t459;
t449 = t370 * t363;
t319 = (pkin(1) * t449 - t355) * qJD(2);
t475 = -qJD(3) * t479 + t317 * t369 - t319 * t365;
t372 = qJD(1) ^ 2;
t473 = pkin(7) * qJD(3);
t263 = -qJD(3) * t402 + t365 * t354 + (t388 - t482) * t369;
t433 = qJD(4) * t368;
t423 = t368 * t263 + t364 * t314 - t343 * t433;
t434 = qJD(4) * t364;
t249 = -t306 * t434 + t423;
t472 = t249 * t364;
t467 = t306 * t364;
t280 = t368 * t343 + t467;
t471 = t280 * t299;
t282 = t306 * t368 - t343 * t364;
t470 = t282 * t299;
t469 = t304 * t343;
t468 = t306 * t343;
t462 = t354 * MDP(8);
t359 = t362 ^ 2;
t461 = t359 * t372;
t456 = t364 * t260;
t455 = t365 * t343;
t452 = t368 * t260;
t451 = t369 * t370;
t315 = -pkin(6) * t421 + t370 * t427;
t316 = t399 * t442;
t445 = t369 * t315 + t365 * t316;
t360 = t366 ^ 2;
t443 = -t370 ^ 2 + t360;
t440 = qJD(2) * t366;
t438 = qJD(3) * t364;
t437 = qJD(3) * t365;
t436 = qJD(3) * t368;
t428 = 0.2e1 * t359;
t425 = t370 * t461;
t424 = t364 * t457;
t419 = t362 * t440;
t417 = t362 * t363 * t372;
t416 = t370 * t432;
t386 = t369 * t275 + t365 * t278 - t294 * t437 + t298 * t435;
t243 = pkin(8) * t314 + t386;
t403 = t482 * pkin(6) - t366 * t404 + t370 * t426;
t276 = -pkin(2) * t354 - t403;
t246 = pkin(3) * t264 - pkin(8) * t263 + t276;
t411 = -t364 * t243 + t368 * t246;
t410 = t263 * t364 - t368 * t314;
t408 = t368 * t299;
t406 = qJD(2) + 0.2e1 * t441;
t405 = t354 + t431;
t397 = -g(1) * t326 - g(2) * t324;
t297 = (t364 * t366 + t368 * t451) * t442;
t396 = t368 * t435 - t297;
t395 = t368 * t243 + t364 * t246;
t293 = -t407 * pkin(2) - t315;
t256 = t304 * pkin(3) - t306 * pkin(8) + t293;
t258 = -pkin(8) * t343 + t262;
t247 = t256 * t368 - t258 * t364;
t248 = t256 * t364 + t258 * t368;
t311 = t355 + (-pkin(1) * t370 - pkin(2)) * t363;
t322 = t363 * t365 + t366 * t458;
t265 = pkin(3) * t321 - pkin(8) * t322 + t311;
t267 = -pkin(8) * t457 + t479;
t394 = t265 * t368 - t267 * t364;
t393 = t265 * t364 + t267 * t368;
t261 = -t294 * t365 + t298 * t369;
t391 = -t312 * t365 + t313 * t369;
t285 = t322 * t364 + t368 * t457;
t385 = -t312 * t437 + t313 * t435 + t365 * t317 + t369 * t319;
t325 = t367 * t449 + t453;
t382 = -g(1) * t325 - g(2) * t323 + g(3) * t457;
t320 = t444 * qJD(2);
t378 = -t276 - t382;
t257 = pkin(3) * t343 - t261;
t377 = -pkin(8) * t260 + (t257 + t261) * t299;
t376 = -pkin(7) * t314 - t343 * t293;
t374 = pkin(7) * qJD(4) * t299 + t382;
t373 = -g(3) * t459 + (pkin(8) * t421 - qJD(4) * t339 + t445) * t299 + t397;
t296 = t364 * t369 * t420 - t368 * t421;
t291 = t326 * t369 + t367 * t460;
t286 = t322 * t368 - t424;
t284 = -t321 * qJD(3) + t439 * t458;
t283 = t322 * qJD(3) + t362 * t418;
t271 = t291 * t368 + t325 * t364;
t270 = -t291 * t364 + t325 * t368;
t268 = -pkin(3) * t421 + t315 * t365 - t316 * t369;
t266 = pkin(3) * t457 - t391;
t255 = -t285 * qJD(4) + t284 * t368 + t364 * t419;
t254 = -qJD(4) * t424 + t284 * t364 + t322 * t433 - t368 * t419;
t253 = pkin(3) * t283 - pkin(8) * t284 + t320;
t252 = -pkin(3) * t419 - t475;
t251 = pkin(8) * t419 + t385;
t250 = t282 * qJD(4) + t410;
t242 = -t248 * qJD(4) + t411;
t241 = t247 * qJD(4) + t395;
t1 = [qJDD(1) * MDP(1) + (g(1) * t367 - g(2) * t371) * MDP(2) + (g(1) * t371 + g(2) * t367) * MDP(3) + t363 * t462 + (-t320 * t407 - t355 * t354 + t403 * t363 + g(1) * t324 - g(2) * t326 + (t354 * t449 + (-t366 * t432 + t430) * t428) * pkin(1)) * MDP(9) + (-t319 * t407 - t444 * t354 - t379 * t363 - g(1) * t323 + g(2) * t325 + (-t416 - t429) * pkin(1) * t428) * MDP(10) + (t263 * t322 + t284 * t306) * MDP(11) + (-t263 * t321 - t264 * t322 - t283 * t306 - t284 * t304) * MDP(12) + (-t284 * t343 + t314 * t322) * MDP(13) + (t283 * t343 - t314 * t321) * MDP(14) + (g(1) * t288 - g(2) * t291 + t311 * t264 + t276 * t321 + t293 * t283 + t320 * t304 + t391 * t314 - t475 * t343) * MDP(16) + (-g(1) * t464 - g(2) * t290 + t311 * t263 + t276 * t322 + t293 * t284 + t320 * t306 - t314 * t479 + t385 * t343) * MDP(17) + (t249 * t286 + t255 * t282) * MDP(18) + (-t249 * t285 - t250 * t286 - t254 * t282 - t255 * t280) * MDP(19) + (t249 * t321 + t255 * t299 + t260 * t286 + t282 * t283) * MDP(20) + (-t250 * t321 - t254 * t299 - t260 * t285 - t280 * t283) * MDP(21) + (t260 * t321 + t283 * t299) * MDP(22) + ((-t393 * qJD(4) - t251 * t364 + t253 * t368) * t299 + t394 * t260 + t242 * t321 + t247 * t283 + t252 * t280 + t266 * t250 + t244 * t285 + t257 * t254 + g(1) * t480 - g(2) * t271) * MDP(23) + (-(t394 * qJD(4) + t251 * t368 + t253 * t364) * t299 - t393 * t260 - t241 * t321 - t248 * t283 + t252 * t282 + t266 * t249 + t244 * t286 + t257 * t255 - g(1) * t481 - g(2) * t270) * MDP(24) + ((qJDD(1) * t360 + 0.2e1 * t366 * t416) * MDP(4) + 0.2e1 * (t370 * t429 - t443 * t432) * MDP(5)) * t359 + ((t405 * t366 + t406 * t439) * MDP(6) + (t405 * t370 - t406 * t440) * MDP(7) + (-t263 * t370 + t306 * t440) * MDP(13) + (t264 * t370 - t304 * t440) * MDP(14) + (-t314 * t370 - t343 * t440) * MDP(15) + (t261 * t440 - t375 * t370) * MDP(16) + (-g(1) * t450 - t262 * t440 + t386 * t370) * MDP(17)) * t362; -t366 * MDP(4) * t425 + t443 * MDP(5) * t461 + (-t370 * t417 + t414) * MDP(6) + (t366 * t417 + t353) * MDP(7) + t462 + (t318 * t407 + t461 * t474 - t382 + t403) * MDP(9) + (t315 * t407 + pkin(1) * t425 + (pkin(6) * t432 + g(3)) * t459 - t397 + t422) * MDP(10) + (t263 * t365 - t369 * t468) * MDP(11) + ((t263 + t469) * t369 + (-t264 + t468) * t365) * MDP(12) + (-t343 * t435 + t365 * t314 + (-t306 * t366 + t343 * t451) * t442) * MDP(13) + (t343 * t437 + t314 * t369 + (t304 * t366 - t370 * t455) * t442) * MDP(14) + t343 * MDP(15) * t421 + (-t261 * t421 - pkin(2) * t264 - t318 * t304 + (-t315 * t343 + t376) * t365 + ((t316 + t473) * t343 + t378) * t369) * MDP(16) + (-pkin(2) * t263 - t445 * t343 + t262 * t421 - t318 * t306 + t376 * t369 + (-t343 * t473 - t378) * t365) * MDP(17) + (t249 * t365 * t368 + (-t365 * t434 + t396) * t282) * MDP(18) + (t280 * t297 + t282 * t296 + (-t280 * t368 - t282 * t364) * t435 + (-t472 - t250 * t368 + (t280 * t364 - t282 * t368) * qJD(4)) * t365) * MDP(19) + (-t249 * t369 + t396 * t299 + (-t343 * t282 - t299 * t434 + t452) * t365) * MDP(20) + (t250 * t369 + (-t364 * t435 + t296) * t299 + (t343 * t280 - t299 * t433 - t456) * t365) * MDP(21) + (-t260 * t369 - t299 * t455) * MDP(22) + (-t257 * t296 - t268 * t280 - t477 * t368 + t373 * t364 + (t257 * t438 - t242 + (qJD(3) * t280 - t456) * pkin(7) - t374 * t368) * t369 + (t257 * t433 + t244 * t364 - t343 * t247 + (t299 * t438 + t250) * pkin(7)) * t365) * MDP(23) + (-t257 * t297 - t268 * t282 + t477 * t364 + t373 * t368 + (t257 * t436 + t241 + (qJD(3) * t282 - t452) * pkin(7) + t374 * t364) * t369 + (-t257 * t434 + t244 * t368 + t343 * t248 + (t299 * t436 + t249) * pkin(7)) * t365) * MDP(24); -t304 ^ 2 * MDP(12) + (t263 - t469) * MDP(13) + (-t264 - t468) * MDP(14) + t314 * MDP(15) + (-t262 * t343 + t375 - t384) * MDP(16) + (g(1) * t291 + g(2) * t288 + g(3) * t322 - t261 * t343 + t293 * t304 - t386) * MDP(17) + (t282 * t408 + t472) * MDP(18) + ((t249 - t471) * t368 + (-t250 - t470) * t364) * MDP(19) + (t299 * t408 + t456) * MDP(20) + (-t299 ^ 2 * t364 + t452) * MDP(21) + (-pkin(3) * t250 - t262 * t280 + t377 * t364 - t478 * t368) * MDP(23) + (-pkin(3) * t249 - t262 * t282 + t478 * t364 + t377 * t368) * MDP(24) + (t304 * MDP(11) + t306 * MDP(12) - t293 * MDP(16) - t282 * MDP(20) + t280 * MDP(21) - t299 * MDP(22) - t247 * MDP(23) + t248 * MDP(24)) * t306; t282 * t280 * MDP(18) + (-t280 ^ 2 + t282 ^ 2) * MDP(19) + (t423 + t471) * MDP(20) + (-t410 + t470) * MDP(21) + t260 * MDP(22) + (-g(1) * t270 + g(2) * t481 + g(3) * t285 + t248 * t299 - t257 * t282 + t411) * MDP(23) + (g(1) * t271 + g(2) * t480 + g(3) * t286 + t247 * t299 + t257 * t280 - t395) * MDP(24) + (-MDP(20) * t467 - t282 * MDP(21) - t248 * MDP(23) - t247 * MDP(24)) * qJD(4);];
tau = t1;
