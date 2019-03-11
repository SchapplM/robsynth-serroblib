% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:51
% EndTime: 2019-03-09 01:58:57
% DurationCPUTime: 3.98s
% Computational Cost: add. (3610->394), mult. (7703->495), div. (0->0), fcn. (5584->14), ass. (0->178)
t378 = sin(pkin(9));
t354 = pkin(1) * t378 + qJ(3);
t340 = qJD(1) * qJD(3) + qJDD(1) * t354;
t379 = cos(pkin(10));
t362 = t379 * qJDD(2);
t377 = sin(pkin(10));
t468 = pkin(7) * qJDD(1);
t307 = t362 + (-t340 - t468) * t377;
t321 = t377 * qJDD(2) + t379 * t340;
t308 = t379 * t468 + t321;
t384 = sin(qJ(4));
t387 = cos(qJ(4));
t346 = t354 * qJD(1);
t364 = t379 * qJD(2);
t314 = t364 + (-pkin(7) * qJD(1) - t346) * t377;
t327 = t377 * qJD(2) + t379 * t346;
t437 = qJD(1) * t379;
t315 = pkin(7) * t437 + t327;
t276 = t384 * t314 + t387 * t315;
t490 = qJD(4) * t276;
t401 = t307 * t387 - t384 * t308 - t490;
t256 = -qJDD(4) * pkin(4) - t401;
t375 = pkin(10) + qJ(4);
t365 = sin(t375);
t367 = cos(t375);
t376 = qJ(1) + pkin(9);
t366 = sin(t376);
t368 = cos(t376);
t415 = g(1) * t368 + g(2) * t366;
t392 = -g(3) * t367 + t365 * t415;
t435 = qJD(1) * t387;
t353 = t379 * t435;
t436 = qJD(1) * t384;
t422 = t377 * t436;
t334 = t353 - t422;
t329 = qJD(5) - t334;
t434 = qJD(5) * t329;
t491 = -pkin(8) * t434 - t256 + t392;
t408 = t307 * t384 + t308 * t387;
t483 = t314 * t387 - t384 * t315;
t255 = qJDD(4) * pkin(8) + qJD(4) * t483 + t408;
t273 = qJD(4) * pkin(8) + t276;
t357 = pkin(3) * t379 + pkin(2);
t380 = cos(pkin(9));
t479 = pkin(1) * t380;
t345 = -t357 - t479;
t331 = qJD(1) * t345 + qJD(3);
t343 = t377 * t387 + t379 * t384;
t335 = t343 * qJD(1);
t281 = -pkin(4) * t334 - pkin(8) * t335 + t331;
t383 = sin(qJ(5));
t386 = cos(qJ(5));
t259 = t273 * t386 + t281 * t383;
t337 = t343 * qJD(4);
t427 = qJDD(1) * t387;
t352 = t379 * t427;
t428 = qJDD(1) * t384;
t411 = -t377 * t428 + t352;
t301 = qJD(1) * t337 - t411;
t328 = qJDD(1) * t345 + qJDD(3);
t423 = qJD(4) * t353 + t377 * t427 + t379 * t428;
t395 = qJD(4) * t422 - t423;
t264 = t301 * pkin(4) + pkin(8) * t395 + t328;
t263 = t386 * t264;
t431 = t386 * qJD(4);
t433 = qJD(5) * t383;
t269 = -qJD(5) * t431 - t383 * qJDD(4) + t335 * t433 + t386 * t395;
t298 = qJDD(5) + t301;
t318 = qJD(4) * t383 + t335 * t386;
t245 = pkin(5) * t298 + qJ(6) * t269 - qJD(5) * t259 - qJD(6) * t318 - t255 * t383 + t263;
t316 = t335 * t383 - t431;
t252 = -qJ(6) * t316 + t259;
t489 = t252 * t329 + t245;
t270 = qJD(5) * t318 - t386 * qJDD(4) - t383 * t395;
t432 = qJD(5) * t386;
t396 = t386 * t255 + t383 * t264 - t273 * t433 + t281 * t432;
t246 = -qJ(6) * t270 - qJD(6) * t316 + t396;
t258 = -t273 * t383 + t386 * t281;
t251 = -qJ(6) * t318 + t258;
t250 = pkin(5) * t329 + t251;
t488 = -t250 * t329 + t246;
t419 = -g(1) * t366 + g(2) * t368;
t485 = t419 * t365;
t449 = t368 * t386;
t452 = t366 * t383;
t322 = t367 * t452 + t449;
t450 = t368 * t383;
t451 = t366 * t386;
t324 = -t367 * t450 + t451;
t482 = -g(1) * t324 + g(2) * t322;
t472 = g(3) * t365;
t481 = t415 * t367 + t472;
t480 = t318 ^ 2;
t385 = sin(qJ(1));
t475 = g(1) * t385;
t470 = pkin(7) + t354;
t469 = qJ(6) + pkin(8);
t467 = t269 * t383;
t465 = t316 * t334;
t464 = t316 * t335;
t463 = t316 * t383;
t462 = t318 * t329;
t461 = t318 * t335;
t460 = t318 * t386;
t459 = (-t346 * t377 + t364) * t377;
t458 = t329 * t383;
t457 = t334 * t383;
t342 = t377 * t384 - t387 * t379;
t336 = t342 * qJD(4);
t456 = t336 * t383;
t455 = t336 * t386;
t454 = t343 * t383;
t453 = t343 * t386;
t448 = t379 * MDP(5);
t338 = t470 * t377;
t339 = t470 * t379;
t295 = -t338 * t384 + t339 * t387;
t287 = t386 * t295;
t289 = t386 * t298;
t447 = -t251 + t250;
t446 = -t270 * t453 + t316 * t455;
t445 = -t383 * t270 - t316 * t432;
t299 = pkin(4) * t335 - pkin(8) * t334;
t444 = t383 * t299 + t386 * t483;
t443 = -t269 * t342 + t318 * t337;
t288 = pkin(4) * t342 - pkin(8) * t343 + t345;
t442 = t383 * t288 + t287;
t441 = t329 * t457 + t289;
t418 = qJD(5) * t469;
t440 = qJ(6) * t457 + qJD(6) * t386 - t383 * t418 - t444;
t292 = t386 * t299;
t439 = -pkin(5) * t335 - t292 + (qJ(6) * t334 - t418) * t386 + (-qJD(6) + t483) * t383;
t438 = t377 ^ 2 + t379 ^ 2;
t358 = -pkin(2) - t479;
t429 = qJDD(1) * t358;
t425 = t318 * t456;
t405 = -t338 * t387 - t339 * t384;
t278 = -qJD(3) * t342 + qJD(4) * t405;
t300 = pkin(4) * t337 + pkin(8) * t336;
t424 = t386 * t278 + t288 * t432 + t383 * t300;
t421 = t343 * t432;
t420 = pkin(5) * t383 + pkin(7) + qJ(3);
t417 = t329 * t386;
t416 = -qJD(5) * t281 - t255;
t388 = cos(qJ(1));
t413 = -g(2) * t388 + t475;
t412 = -t273 * t432 + t263;
t410 = -t250 * t386 - t252 * t383;
t409 = -t270 * t342 - t316 * t337;
t320 = -t340 * t377 + t362;
t406 = -t320 * t377 + t321 * t379;
t360 = pkin(5) * t386 + pkin(4);
t404 = t360 * t367 + t365 * t469;
t402 = qJ(6) * t336 - qJD(6) * t343;
t272 = -qJD(4) * pkin(4) - t483;
t400 = t357 + t404;
t399 = t421 - t456;
t398 = -t343 * t433 - t455;
t397 = -pkin(8) * t298 + t272 * t329;
t249 = pkin(5) * t270 + qJDD(6) + t256;
t279 = qJD(3) * t343 + qJD(4) * t295;
t372 = t388 * pkin(1);
t348 = t469 * t386;
t347 = t469 * t383;
t344 = qJDD(3) + t429;
t325 = t367 * t449 + t452;
t323 = -t367 * t451 + t450;
t313 = t316 ^ 2;
t303 = -qJD(4) * t337 - qJDD(4) * t342;
t302 = -qJD(4) * t336 + qJDD(4) * t343;
t293 = t386 * t300;
t285 = t386 * t288;
t265 = pkin(5) * t316 + qJD(6) + t272;
t261 = -qJ(6) * t454 + t442;
t260 = pkin(5) * t342 - qJ(6) * t453 - t295 * t383 + t285;
t248 = -qJ(6) * t421 + (-qJD(5) * t295 + t402) * t383 + t424;
t247 = pkin(5) * t337 - t278 * t383 + t293 + t402 * t386 + (-t287 + (qJ(6) * t343 - t288) * t383) * qJD(5);
t1 = [qJDD(1) * MDP(1) + t413 * MDP(2) + (g(1) * t388 + g(2) * t385) * MDP(3) + (t413 + (t378 ^ 2 + t380 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t340 * t438 + t406 - t415) * MDP(7) + (t344 * t358 - g(1) * (-pkin(1) * t385 - pkin(2) * t366 + qJ(3) * t368) - g(2) * (pkin(2) * t368 + qJ(3) * t366 + t372) + t406 * t354 + (t327 * t379 - t459) * qJD(3)) * MDP(8) + (-t335 * t336 - t343 * t395) * MDP(9) + (-t343 * t301 - t336 * t334 - t335 * t337 + t342 * t395) * MDP(10) + t302 * MDP(11) + t303 * MDP(12) + (-qJD(4) * t279 + qJDD(4) * t405 + t301 * t345 + t328 * t342 + t331 * t337 - t367 * t419) * MDP(14) + (-t278 * qJD(4) - t295 * qJDD(4) + t328 * t343 - t331 * t336 - t345 * t395 + t485) * MDP(15) + (-t269 * t453 + t318 * t398) * MDP(16) + (t425 + (t467 + (-t460 + t463) * qJD(5)) * t343 + t446) * MDP(17) + (t289 * t343 + t329 * t398 + t443) * MDP(18) + (-t298 * t454 - t329 * t399 + t409) * MDP(19) + (t298 * t342 + t329 * t337) * MDP(20) + ((-t295 * t432 + t293) * t329 + t285 * t298 + t412 * t342 + t258 * t337 + t279 * t316 - t405 * t270 + t272 * t421 - g(1) * t323 - g(2) * t325 + ((-qJD(5) * t288 - t278) * t329 - t295 * t298 + t416 * t342 + t256 * t343 - t272 * t336) * t383) * MDP(21) + (-(-t295 * t433 + t424) * t329 - t442 * t298 - t396 * t342 - t259 * t337 + t279 * t318 + t405 * t269 + t256 * t453 - g(1) * t322 - g(2) * t324 + t398 * t272) * MDP(22) + (-t247 * t318 - t248 * t316 + t260 * t269 - t261 * t270 - t485 - t410 * t336 + (-t245 * t386 - t246 * t383 + (t250 * t383 - t252 * t386) * qJD(5)) * t343) * MDP(23) + (t246 * t261 + t252 * t248 + t245 * t260 + t250 * t247 + t249 * (pkin(5) * t454 - t405) + t265 * (pkin(5) * t399 + t279) + pkin(1) * t475 - g(2) * t372 + (-g(1) * t420 - g(2) * t400) * t368 + (g(1) * t400 - g(2) * t420) * t366) * MDP(24) + (-t377 * MDP(6) + t448) * (-t344 - t419 - t429); (qJDD(2) - g(3)) * MDP(4) + (t320 * t379 + t321 * t377 - g(3)) * MDP(8) + t303 * MDP(14) - t302 * MDP(15) + (t329 * t456 - t409) * MDP(21) + (t329 * t455 + t443) * MDP(22) + (-t425 + t446) * MDP(23) + (t249 * t342 + t250 * t456 - t252 * t455 + t265 * t337 - g(3)) * MDP(24) + (-MDP(23) * t467 + (-t245 * t383 + t246 * t386) * MDP(24) + (-t383 * MDP(21) - t386 * MDP(22)) * t298 + ((t460 + t463) * MDP(23) + t410 * MDP(24) + (-t386 * MDP(21) + t383 * MDP(22)) * t329) * qJD(5)) * t343; (qJD(1) * t459 - t327 * t437 + qJDD(3) + t419) * MDP(8) - t352 * MDP(14) + t423 * MDP(15) + (t441 - t464) * MDP(21) - MDP(22) * t461 + t445 * MDP(23) + (-t265 * t335 + t419) * MDP(24) - t438 * MDP(7) * qJD(1) ^ 2 + (-MDP(21) * t434 - t298 * MDP(22) + MDP(23) * t462 + MDP(24) * t488) * t383 + (-t448 + t358 * MDP(8) + (MDP(14) * t384 + MDP(6)) * t377) * qJDD(1) + ((t377 * t435 + t379 * t436 + t335) * MDP(14) + (t334 - t422) * MDP(15)) * qJD(4) + ((t269 + t465) * MDP(23) + t489 * MDP(24) - t329 ^ 2 * MDP(22)) * t386; -t334 ^ 2 * MDP(10) + ((-t334 - t422) * qJD(4) + t423) * MDP(11) + t411 * MDP(12) + qJDD(4) * MDP(13) + (t392 + t401 + t490) * MDP(14) + (-t331 * t334 - t408 + t481) * MDP(15) + (t318 * t417 - t467) * MDP(16) + ((-t269 + t465) * t386 - t318 * t458 + t445) * MDP(17) + (t298 * t383 + t329 * t417 - t461) * MDP(18) + (-t329 * t433 + t441 + t464) * MDP(19) + (-pkin(4) * t270 - t276 * t316 - t292 * t329 + (t329 * t483 + t397) * t383 + t491 * t386) * MDP(21) + (pkin(4) * t269 - t276 * t318 + t444 * t329 - t383 * t491 + t397 * t386) * MDP(22) + (-t269 * t347 - t270 * t348 - t440 * t316 - t439 * t318 - t383 * t489 + t488 * t386 - t481) * MDP(23) + (t246 * t348 - t245 * t347 - t249 * t360 - g(3) * t404 + (pkin(5) * t458 - t276) * t265 + t440 * t252 + t439 * t250 + t415 * (t360 * t365 - t367 * t469)) * MDP(24) + (t335 * MDP(10) - t331 * MDP(14) - t329 * MDP(20) - t258 * MDP(21) + t259 * MDP(22) - MDP(9) * t334) * t335; t318 * t316 * MDP(16) + (-t313 + t480) * MDP(17) + (t316 * t329 - t269) * MDP(18) + (-t270 + t462) * MDP(19) + t298 * MDP(20) + (t259 * t329 - t272 * t318 + (t416 + t472) * t383 + t412 + t482) * MDP(21) + (g(1) * t325 - g(2) * t323 + t258 * t329 + t272 * t316 + t386 * t472 - t396) * MDP(22) + (pkin(5) * t269 - t316 * t447) * MDP(23) + (t447 * t252 + (-t265 * t318 + t383 * t472 + t245 + t482) * pkin(5)) * MDP(24); (-t313 - t480) * MDP(23) + (t250 * t318 + t252 * t316 + t249 - t392) * MDP(24);];
tau  = t1;
