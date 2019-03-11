% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:27
% EndTime: 2019-03-09 02:36:36
% DurationCPUTime: 4.33s
% Computational Cost: add. (3467->388), mult. (7960->516), div. (0->0), fcn. (6018->8), ass. (0->166)
t373 = cos(pkin(10));
t456 = cos(qJ(4));
t412 = t456 * t373;
t399 = qJD(1) * t412;
t377 = sin(qJ(4));
t372 = sin(pkin(10));
t429 = qJD(1) * t372;
t408 = t377 * t429;
t337 = t399 - t408;
t376 = sin(qJ(5));
t379 = cos(qJ(5));
t419 = t379 * qJD(4);
t316 = t337 * t376 - t419;
t391 = -t456 * t372 - t377 * t373;
t461 = t391 * qJD(1);
t472 = qJD(5) - t461;
t477 = t316 * t472;
t378 = cos(qJ(6));
t318 = qJD(4) * t376 + t337 * t379;
t375 = sin(qJ(6));
t446 = t318 * t375;
t266 = t378 * t316 + t446;
t329 = qJD(6) + t472;
t476 = t266 * t329;
t394 = t316 * t375 - t378 * t318;
t475 = t329 * t394;
t344 = t375 * t379 + t376 * t378;
t459 = qJD(5) + qJD(6);
t312 = t459 * t344;
t431 = -t344 * t461 + t312;
t425 = qJD(5) * t376;
t442 = t461 * t376;
t474 = t425 - t442;
t402 = t379 * t472;
t353 = qJD(4) * t408;
t331 = qJD(4) * t399 - t353;
t437 = t376 * t331;
t473 = -t402 * t472 - t437;
t374 = -pkin(1) - qJ(3);
t462 = t374 * qJD(1);
t354 = qJD(2) + t462;
t406 = -pkin(7) * qJD(1) + t354;
t332 = t406 * t372;
t333 = t406 * t373;
t292 = t456 * t332 + t377 * t333;
t286 = qJD(4) * pkin(8) + t292;
t371 = qJD(1) * qJ(2);
t365 = qJD(3) + t371;
t349 = pkin(3) * t429 + t365;
t287 = -pkin(4) * t461 - pkin(8) * t337 + t349;
t257 = t286 * t379 + t287 * t376;
t249 = -pkin(9) * t316 + t257;
t423 = qJD(6) * t375;
t247 = t249 * t423;
t466 = -t377 * t332 + t456 * t333;
t285 = -qJD(4) * pkin(4) - t466;
t260 = t316 * pkin(5) + t285;
t471 = t260 * t266 + t247;
t330 = t461 * qJD(4);
t276 = qJD(5) * t419 + t379 * t330 - t337 * t425;
t382 = t391 * qJD(3);
t262 = qJD(1) * t382 + qJD(4) * t466;
t370 = qJD(1) * qJD(2);
t288 = pkin(4) * t331 - pkin(8) * t330 + t370;
t282 = t379 * t288;
t381 = -qJD(5) * t257 - t262 * t376 + t282;
t238 = pkin(5) * t331 - pkin(9) * t276 + t381;
t444 = t330 * t376;
t277 = qJD(5) * t318 + t444;
t424 = qJD(5) * t379;
t385 = t379 * t262 - t286 * t425 + t287 * t424 + t376 * t288;
t241 = -pkin(9) * t277 + t385;
t405 = t378 * t238 - t375 * t241;
t470 = t260 * t394 + t405;
t469 = qJ(2) * MDP(6) + t372 * MDP(7) + t373 * MDP(8) + MDP(5);
t468 = t331 * MDP(29) + (-t266 ^ 2 + t394 ^ 2) * MDP(26) - t266 * t394 * MDP(25);
t465 = t375 * t425 + t376 * t423;
t455 = -pkin(7) + t374;
t347 = t455 * t372;
t348 = t455 * t373;
t464 = -t377 * t347 + t456 * t348;
t463 = -qJD(6) * t379 - t424;
t430 = t372 ^ 2 + t373 ^ 2;
t460 = t430 * qJD(3);
t343 = t375 * t376 - t378 * t379;
t432 = (t461 - t459) * t343;
t443 = t331 * t344;
t458 = -t432 * t329 - t443;
t404 = t276 * t375 + t378 * t277;
t245 = -qJD(6) * t394 + t404;
t457 = pkin(8) + pkin(9);
t256 = -t286 * t376 + t379 * t287;
t248 = -pkin(9) * t318 + t256;
t246 = pkin(5) * t472 + t248;
t453 = t246 * t378;
t452 = t249 * t378;
t451 = t266 * t337;
t450 = t394 * t337;
t449 = t276 * t376;
t448 = t316 * t337;
t447 = t318 * t337;
t407 = qJD(4) * t456;
t426 = qJD(4) * t377;
t339 = -t372 * t426 + t373 * t407;
t445 = t329 * t339;
t390 = -t377 * t372 + t412;
t441 = t390 * t376;
t440 = t390 * t379;
t308 = t343 * t331;
t338 = -t372 * t407 - t373 * t426;
t436 = t376 * t338;
t310 = t456 * t347 + t377 * t348;
t306 = t379 * t310;
t322 = t379 * t331;
t435 = t379 * t338;
t359 = t372 * pkin(3) + qJ(2);
t305 = pkin(4) * t337 - pkin(8) * t461;
t434 = t376 * t305 + t379 * t466;
t304 = -pkin(4) * t391 - pkin(8) * t390 + t359;
t433 = t376 * t304 + t306;
t428 = qJD(4) * t338;
t427 = qJD(4) * t339;
t422 = qJD(6) * t378;
t416 = t378 * t276 - t375 * t277 - t316 * t422;
t415 = qJD(5) * t457;
t411 = t390 * t424;
t403 = qJD(1) * t430;
t401 = qJD(6) * t246 + t241;
t400 = -qJD(5) * t391 + qJD(1);
t398 = pkin(5) * t474 - t292;
t397 = -t329 * t431 - t308;
t299 = t379 * t305;
t351 = t457 * t379;
t396 = pkin(5) * t337 + qJD(6) * t351 - t466 * t376 + t299 + (-pkin(9) * t461 + t415) * t379;
t350 = t457 * t376;
t395 = -pkin(9) * t442 + qJD(6) * t350 + t376 * t415 + t434;
t240 = t246 * t375 + t452;
t393 = -t472 * t474 + t322;
t389 = -t411 - t436;
t388 = -t390 * t425 + t435;
t387 = -pkin(8) * t331 + t285 * t472;
t386 = t329 * t343;
t279 = qJD(4) * t464 + t382;
t302 = pkin(4) * t339 - pkin(8) * t338 + qJD(2);
t384 = t379 * t279 + t376 * t302 + t304 * t424 - t310 * t425;
t244 = -t318 * t423 + t416;
t263 = qJD(3) * t337 + t332 * t407 + t333 * t426;
t280 = qJD(3) * t390 + qJD(4) * t310;
t380 = qJD(1) ^ 2;
t364 = -pkin(5) * t379 - pkin(4);
t307 = t331 * t391;
t301 = t343 * t390;
t300 = t344 * t390;
t297 = t379 * t304;
t295 = t379 * t302;
t283 = pkin(5) * t441 - t464;
t259 = -pkin(5) * t389 + t280;
t258 = -pkin(9) * t441 + t433;
t255 = -pkin(5) * t391 - pkin(9) * t440 - t310 * t376 + t297;
t253 = pkin(5) * t277 + t263;
t251 = t375 * t435 + (t440 * t459 + t436) * t378 - t465 * t390;
t250 = -t312 * t390 - t338 * t343;
t243 = pkin(9) * t389 + t384;
t242 = -pkin(9) * t435 + pkin(5) * t339 - t279 * t376 + t295 + (-t306 + (pkin(9) * t390 - t304) * t376) * qJD(5);
t239 = -t249 * t375 + t453;
t1 = [(t339 * t472 - t307) * MDP(22) + (t276 * t440 + t318 * t388) * MDP(18) + (-t257 * t339 + t263 * t440 - t276 * t464 + t280 * t318 + t388 * t285 - t433 * t331 - t384 * t472 + t385 * t391) * MDP(24) + (-t244 * t391 + t250 * t329 - t301 * t331 - t339 * t394) * MDP(27) + (-t240 * t339 + t283 * t244 - t247 * t391 + t260 * t250 - t253 * t301 - t259 * t394 + (-(-qJD(6) * t258 + t242) * t329 - t255 * t331 + t238 * t391) * t375 + (-(qJD(6) * t255 + t243) * t329 - t258 * t331 + t401 * t391) * t378) * MDP(31) + (t245 * t391 - t251 * t329 - t266 * t339 - t300 * t331) * MDP(28) + ((t242 * t378 - t243 * t375) * t329 + (t255 * t378 - t258 * t375) * t331 - t405 * t391 + t239 * t339 + t259 * t266 + t283 * t245 + t253 * t300 + t260 * t251 + ((-t255 * t375 - t258 * t378) * t329 + t240 * t391) * qJD(6)) * MDP(30) + (t277 * t391 - t316 * t339 + t389 * t472 - t390 * t437) * MDP(21) + (-t276 * t391 + t318 * t339 + t322 * t390 + t388 * t472) * MDP(20) + (t330 * t391 - t331 * t390 - t337 * t339 + t338 * t461) * MDP(12) + (-0.2e1 * qJD(2) * t461 - qJD(4) * t280 + t331 * t359 + t339 * t349) * MDP(16) + ((-t310 * t424 + t295) * t472 + t297 * t331 - (-t286 * t424 + t282) * t391 + t256 * t339 + t280 * t316 - t464 * t277 + t285 * t411 + ((-qJD(5) * t304 - t279) * t472 - t310 * t331 - (-qJD(5) * t287 - t262) * t391 + t263 * t390 + t285 * t338) * t376) * MDP(23) + (-t244 * t301 - t250 * t394) * MDP(25) + (-t244 * t300 + t245 * t301 - t250 * t266 + t251 * t394) * MDP(26) + (t330 * t390 + t337 * t338) * MDP(11) + (-qJD(4) * t279 + t330 * t359 + t338 * t349 + (qJD(1) * t390 + t337) * qJD(2)) * MDP(17) + ((-t316 * t379 - t318 * t376) * t338 - (t449 + t277 * t379 + (-t316 * t376 + t318 * t379) * qJD(5)) * t390) * MDP(19) + ((t365 + t371) * qJD(2) + (-t354 - t462) * t460) * MDP(10) + 0.2e1 * qJD(3) * MDP(9) * t403 - MDP(14) * t427 + (-t307 + t445) * MDP(29) + MDP(13) * t428 + 0.2e1 * t469 * t370; (-t365 - t460) * qJD(1) * MDP(10) + (qJD(1) * t461 + t428) * MDP(16) + (-qJD(1) * t337 - t427) * MDP(17) + (t391 * t437 - t277 * t390 - t316 * t338 + (-t339 * t376 - t379 * t400) * t472) * MDP(23) + (t391 * t322 - t276 * t390 - t318 * t338 + (-t339 * t379 + t376 * t400) * t472) * MDP(24) + (-t390 * t245 - t338 * t266 - t344 * t445 + qJD(1) * t386 - ((t378 * t463 + t465) * t329 - t443) * t391) * MDP(30) + (-t390 * t244 + t338 * t394 + t339 * t386 + t344 * t329 * qJD(1) - (-(t375 * t463 - t376 * t422 - t378 * t425) * t329 + t308) * t391) * MDP(31) - t469 * t380; -t430 * t380 * MDP(9) + (t354 * t403 + t370) * MDP(10) + (-t353 + (t399 + t337) * qJD(4)) * MDP(16) + 0.2e1 * MDP(17) * t330 + (t393 - t448) * MDP(23) + (-t447 + t473) * MDP(24) + (t397 - t451) * MDP(30) + (t450 + t458) * MDP(31); -t331 * MDP(14) + (qJD(4) * t292 - t263) * MDP(16) + (t318 * t402 + t449) * MDP(18) + ((t276 - t477) * t379 + (-t318 * t472 - t277) * t376) * MDP(19) + (-t447 - t473) * MDP(20) + (t393 + t448) * MDP(21) + (-pkin(4) * t277 - t263 * t379 - t292 * t316 + (-pkin(8) * t424 - t299) * t472 + (t466 * t472 + t387) * t376) * MDP(23) + (-pkin(4) * t276 + t263 * t376 - t292 * t318 + (pkin(8) * t425 + t434) * t472 + t387 * t379) * MDP(24) + (t244 * t344 - t394 * t432) * MDP(25) + (-t244 * t343 - t245 * t344 - t432 * t266 + t394 * t431) * MDP(26) + (t450 - t458) * MDP(27) + (t397 + t451) * MDP(28) + ((-t350 * t378 - t351 * t375) * t331 + t364 * t245 + t253 * t343 + (t375 * t395 - t378 * t396) * t329 + t398 * t266 + t431 * t260) * MDP(30) + (-(-t350 * t375 + t351 * t378) * t331 + t364 * t244 + t253 * t344 + (t375 * t396 + t378 * t395) * t329 - t398 * t394 + t432 * t260) * MDP(31) + (MDP(12) * t337 + qJD(4) * MDP(14) - t349 * MDP(16) - MDP(22) * t472 - t256 * MDP(23) + t257 * MDP(24) - t329 * MDP(29) - t239 * MDP(30) + t240 * MDP(31)) * t337 + ((-qJD(3) - t349) * MDP(17) - t337 * MDP(11) - MDP(12) * t461) * t461; t318 * t316 * MDP(18) + (-t316 ^ 2 + t318 ^ 2) * MDP(19) + (t276 + t477) * MDP(20) + (-t444 + (-qJD(5) + t472) * t318) * MDP(21) + t331 * MDP(22) + (t257 * t472 - t285 * t318 + t381) * MDP(23) + (t256 * t472 + t285 * t316 - t385) * MDP(24) + (t244 + t476) * MDP(27) + (-t245 - t475) * MDP(28) + (-(-t248 * t375 - t452) * t329 - t240 * qJD(6) + (-t266 * t318 - t329 * t423 + t378 * t331) * pkin(5) + t470) * MDP(30) + ((-t249 * t329 - t238) * t375 + (t248 * t329 - t401) * t378 + (t318 * t394 - t329 * t422 - t375 * t331) * pkin(5) + t471) * MDP(31) + t468; (t416 + t476) * MDP(27) + (-t404 - t475) * MDP(28) + (t240 * t329 + t470) * MDP(30) + (-t375 * t238 + t239 * t329 - t378 * t241 + t471) * MDP(31) + (-MDP(27) * t446 + MDP(28) * t394 - MDP(30) * t240 - MDP(31) * t453) * qJD(6) + t468;];
tauc  = t1;
