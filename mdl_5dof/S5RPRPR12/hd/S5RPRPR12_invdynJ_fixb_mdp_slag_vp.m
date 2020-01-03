% Calculate vector of inverse dynamics joint torques for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:24
% EndTime: 2019-12-31 18:30:32
% DurationCPUTime: 5.86s
% Computational Cost: add. (3024->371), mult. (7226->487), div. (0->0), fcn. (5586->14), ass. (0->163)
t374 = cos(pkin(8));
t452 = cos(qJ(3));
t418 = t452 * t374;
t354 = qJD(1) * t418;
t372 = sin(pkin(8));
t377 = sin(qJ(3));
t430 = t377 * t372;
t417 = qJD(1) * t430;
t331 = -t354 + t417;
t326 = qJD(5) + t331;
t341 = t372 * t452 + t377 * t374;
t333 = t341 * qJD(1);
t371 = sin(pkin(9));
t373 = cos(pkin(9));
t312 = -t373 * qJD(3) + t333 * t371;
t379 = cos(qJ(5));
t314 = qJD(3) * t371 + t333 * t373;
t376 = sin(qJ(5));
t437 = t314 * t376;
t462 = -t379 * t312 - t437;
t463 = t462 * t326;
t423 = qJD(1) * qJD(2);
t447 = pkin(6) + qJ(2);
t453 = qJDD(1) * t447 + t423;
t322 = t453 * t372;
t323 = t453 * t374;
t347 = t447 * t372;
t342 = qJD(1) * t347;
t349 = t447 * t374;
t343 = qJD(1) * t349;
t306 = -t377 * t342 + t452 * t343;
t461 = t306 * qJD(3);
t391 = t322 * t452 + t377 * t323 + t461;
t260 = -qJDD(3) * pkin(3) + qJDD(4) + t391;
t370 = pkin(8) + qJ(3);
t363 = sin(t370);
t365 = cos(t370);
t378 = sin(qJ(1));
t380 = cos(qJ(1));
t410 = g(1) * t380 + g(2) * t378;
t389 = -g(3) * t365 + t363 * t410;
t384 = t260 - t389;
t401 = t312 * t376 - t314 * t379;
t460 = t326 * t401;
t340 = t371 * t379 + t373 * t376;
t335 = t340 * qJD(5);
t428 = t340 * t331 + t335;
t443 = qJDD(1) * pkin(1);
t457 = g(1) * t378 - g(2) * t380;
t400 = -qJDD(2) + t443 + t457;
t458 = t457 * t363;
t398 = t457 * t365;
t396 = t452 * t342 + t377 * t343;
t456 = MDP(4) * t374 - MDP(5) * t372;
t455 = qJ(2) * qJDD(1);
t337 = t341 * qJD(3);
t415 = qJDD(1) * t452;
t422 = qJDD(1) * t377;
t407 = t372 * t422 - t374 * t415;
t303 = qJD(1) * t337 + t407;
t300 = qJDD(5) + t303;
t338 = t371 * t376 - t379 * t373;
t429 = t326 * t338;
t454 = -t300 * t340 + t326 * t429;
t449 = g(3) * t363;
t387 = -t410 * t365 - t449;
t327 = t331 ^ 2;
t451 = pkin(7) * t373;
t446 = pkin(7) + qJ(4);
t442 = t462 * t333;
t441 = t401 * t333;
t439 = t303 * t371;
t438 = t303 * t373;
t436 = t331 * t371;
t394 = t418 - t430;
t336 = t394 * qJD(3);
t435 = t336 * t371;
t434 = t341 * t371;
t433 = t341 * t373;
t432 = t365 * t378;
t431 = t365 * t380;
t420 = qJD(3) * t354 + t372 * t415 + t374 * t422;
t302 = -qJD(3) * t417 + t420;
t357 = pkin(2) * t374 + pkin(1);
t344 = -qJDD(1) * t357 + qJDD(2);
t251 = pkin(3) * t303 - qJ(4) * t302 - qJD(4) * t333 + t344;
t397 = -t377 * t322 + t323 * t452;
t257 = qJDD(3) * qJ(4) + (qJD(4) - t396) * qJD(3) + t397;
t237 = t371 * t251 + t373 * t257;
t273 = pkin(3) * t337 - qJ(4) * t336 - qJD(4) * t341;
t395 = -t347 * t452 - t377 * t349;
t282 = qJD(2) * t394 + qJD(3) * t395;
t249 = t371 * t273 + t373 * t282;
t345 = -qJD(1) * t357 + qJD(2);
t281 = pkin(3) * t331 - qJ(4) * t333 + t345;
t298 = qJD(3) * qJ(4) + t306;
t255 = t371 * t281 + t373 * t298;
t299 = pkin(3) * t333 + qJ(4) * t331;
t262 = t371 * t299 - t373 * t396;
t301 = -pkin(3) * t394 - qJ(4) * t341 - t357;
t311 = -t377 * t347 + t349 * t452;
t264 = t371 * t301 + t373 * t311;
t427 = t372 ^ 2 + t374 ^ 2;
t426 = qJD(5) * t376;
t425 = qJD(5) * t379;
t296 = -qJD(3) * pkin(3) + qJD(4) + t396;
t424 = -qJD(4) + t296;
t287 = -t373 * qJDD(3) + t302 * t371;
t288 = qJDD(3) * t371 + t302 * t373;
t421 = -t376 * t287 + t379 * t288 - t312 * t425;
t414 = t427 * qJD(1) ^ 2;
t236 = t373 * t251 - t257 * t371;
t232 = pkin(4) * t303 - pkin(7) * t288 + t236;
t235 = -pkin(7) * t287 + t237;
t413 = t379 * t232 - t235 * t376;
t248 = t373 * t273 - t282 * t371;
t254 = t373 * t281 - t298 * t371;
t412 = t379 * t287 + t376 * t288;
t261 = t373 * t299 + t371 * t396;
t263 = t373 * t301 - t311 * t371;
t411 = 0.2e1 * t427;
t408 = -t338 * t300 - t428 * t326;
t406 = t232 * t376 + t235 * t379;
t405 = -t236 * t373 - t237 * t371;
t241 = pkin(4) * t331 - pkin(7) * t314 + t254;
t244 = -pkin(7) * t312 + t255;
t233 = t241 * t379 - t244 * t376;
t234 = t241 * t376 + t244 * t379;
t250 = -pkin(4) * t394 - pkin(7) * t433 + t263;
t256 = -pkin(7) * t434 + t264;
t404 = t250 * t379 - t256 * t376;
t403 = t250 * t376 + t256 * t379;
t402 = -t254 * t371 + t255 * t373;
t399 = pkin(3) * t365 + qJ(4) * t363 + t357;
t348 = t446 * t373;
t393 = pkin(4) * t333 + qJD(4) * t371 + qJD(5) * t348 + t331 * t451 + t261;
t346 = t446 * t371;
t392 = pkin(7) * t436 - qJD(4) * t373 + qJD(5) * t346 + t262;
t238 = -t314 * t426 + t421;
t386 = t260 * t341 + t296 * t336 - t410;
t239 = -qJD(5) * t401 + t412;
t383 = t411 * t423 - t410;
t283 = qJD(2) * t341 + qJD(3) * t311;
t369 = pkin(9) + qJ(5);
t364 = cos(t369);
t362 = sin(t369);
t358 = -pkin(4) * t373 - pkin(3);
t320 = t362 * t378 + t364 * t431;
t319 = -t362 * t431 + t364 * t378;
t318 = t362 * t380 - t364 * t432;
t317 = t362 * t432 + t364 * t380;
t293 = t338 * t341;
t292 = t340 * t341;
t284 = pkin(4) * t434 - t395;
t272 = -pkin(4) * t436 + t306;
t266 = t312 * pkin(4) + t296;
t265 = pkin(4) * t435 + t283;
t259 = t336 * t340 + t425 * t433 - t426 * t434;
t258 = -t335 * t341 - t336 * t338;
t243 = -pkin(7) * t435 + t249;
t242 = t287 * pkin(4) + t260;
t240 = pkin(4) * t337 - t336 * t451 + t248;
t1 = [qJDD(1) * MDP(1) + (t411 * t455 + t383) * MDP(6) + (pkin(1) * t400 + (t427 * t455 + t383) * qJ(2)) * MDP(7) + (t302 * t341 + t333 * t336) * MDP(8) + (t302 * t394 - t303 * t341 - t331 * t336 - t333 * t337) * MDP(9) + (qJD(3) * t336 + qJDD(3) * t341) * MDP(10) + (-qJD(3) * t337 + qJDD(3) * t394) * MDP(11) + (-qJD(3) * t283 + qJDD(3) * t395 - t303 * t357 + t337 * t345 - t344 * t394 + t398) * MDP(13) + (-qJD(3) * t282 - qJDD(3) * t311 - t302 * t357 + t336 * t345 + t341 * t344 - t458) * MDP(14) + (-t236 * t394 + t248 * t331 + t254 * t337 + t263 * t303 + t283 * t312 - t287 * t395 + t371 * t386 + t373 * t398) * MDP(15) + (t237 * t394 - t249 * t331 - t255 * t337 - t264 * t303 + t283 * t314 - t288 * t395 - t371 * t398 + t373 * t386) * MDP(16) + (-t248 * t314 - t249 * t312 - t263 * t288 - t264 * t287 + t458 + t405 * t341 + (-t254 * t373 - t255 * t371) * t336) * MDP(17) + (t236 * t263 + t237 * t264 + t254 * t248 + t255 * t249 - t260 * t395 + t296 * t283 + (-g(1) * t447 - g(2) * t399) * t380 + (g(1) * t399 - g(2) * t447) * t378) * MDP(18) + (-t238 * t293 - t258 * t401) * MDP(19) + (-t238 * t292 + t239 * t293 + t258 * t462 + t259 * t401) * MDP(20) + (-t238 * t394 + t258 * t326 - t293 * t300 - t337 * t401) * MDP(21) + (t239 * t394 - t259 * t326 - t292 * t300 + t337 * t462) * MDP(22) + (-t300 * t394 + t326 * t337) * MDP(23) + ((t240 * t379 - t243 * t376) * t326 + t404 * t300 - t413 * t394 + t233 * t337 - t265 * t462 + t284 * t239 + t242 * t292 + t266 * t259 - g(1) * t318 - g(2) * t320 + (t234 * t394 - t326 * t403) * qJD(5)) * MDP(24) + (-(t240 * t376 + t243 * t379) * t326 - t403 * t300 + t406 * t394 - t234 * t337 - t265 * t401 + t284 * t238 - t242 * t293 + t266 * t258 - g(1) * t317 - g(2) * t319 + (t233 * t394 - t326 * t404) * qJD(5)) * MDP(25) + t457 * MDP(2) + t410 * MDP(3) + t456 * (t400 + t443); -MDP(6) * t414 + (-qJ(2) * t414 - t400) * MDP(7) + (0.2e1 * qJD(3) * t333 + t407) * MDP(13) + ((-t331 - t417) * qJD(3) + t420) * MDP(14) + (-t312 * t333 - t327 * t371 + t438) * MDP(15) + (-t314 * t333 - t327 * t373 - t439) * MDP(16) + (-t287 * t371 - t288 * t373 + (-t312 * t373 + t314 * t371) * t331) * MDP(17) + (-t296 * t333 + t331 * t402 - t405 - t457) * MDP(18) + (t408 + t442) * MDP(24) + (t441 + t454) * MDP(25) - t456 * qJDD(1); -t327 * MDP(9) + ((t331 - t417) * qJD(3) + t420) * MDP(10) - t407 * MDP(11) + qJDD(3) * MDP(12) + (t389 - t391 + t461) * MDP(13) + (t345 * t331 - t387 - t397) * MDP(14) + (-qJ(4) * t439 - pkin(3) * t287 - t306 * t312 + (t371 * t424 - t261) * t331 - t384 * t373) * MDP(15) + (-qJ(4) * t438 - pkin(3) * t288 - t306 * t314 + (t373 * t424 + t262) * t331 + t384 * t371) * MDP(16) + (t261 * t314 + t262 * t312 + (-qJ(4) * t287 - qJD(4) * t312 - t254 * t331 + t237) * t373 + (qJ(4) * t288 + qJD(4) * t314 - t255 * t331 - t236) * t371 + t387) * MDP(17) + (-t254 * t261 - t255 * t262 - t296 * t306 + t402 * qJD(4) - t384 * pkin(3) + (-t236 * t371 + t237 * t373 + t387) * qJ(4)) * MDP(18) + (t238 * t340 + t401 * t429) * MDP(19) + (-t238 * t338 - t239 * t340 + t401 * t428 - t429 * t462) * MDP(20) + (t441 - t454) * MDP(21) + (t408 - t442) * MDP(22) + ((-t346 * t379 - t348 * t376) * t300 + t358 * t239 + t242 * t338 + t272 * t462 + (t376 * t392 - t379 * t393) * t326 + t428 * t266 + t389 * t364) * MDP(24) + (-(-t346 * t376 + t348 * t379) * t300 + t358 * t238 + t242 * t340 + t272 * t401 + (t376 * t393 + t379 * t392) * t326 - t429 * t266 - t389 * t362) * MDP(25) + (-t345 * MDP(13) - t254 * MDP(15) + t255 * MDP(16) - t326 * MDP(23) - t233 * MDP(24) + t234 * MDP(25) + MDP(8) * t331 + t333 * MDP(9)) * t333; (t314 * t331 + t287) * MDP(15) + (-t312 * t331 + t288) * MDP(16) + (-t312 ^ 2 - t314 ^ 2) * MDP(17) + (t254 * t314 + t255 * t312 + t384) * MDP(18) + (t239 - t460) * MDP(24) + (t238 + t463) * MDP(25); t401 * t462 * MDP(19) + (t401 ^ 2 - t462 ^ 2) * MDP(20) + (t421 - t463) * MDP(21) + (-t412 - t460) * MDP(22) + t300 * MDP(23) + (-g(1) * t319 + g(2) * t317 + t234 * t326 + t266 * t401 + t362 * t449 + t413) * MDP(24) + (g(1) * t320 - g(2) * t318 + t233 * t326 - t266 * t462 + t364 * t449 - t406) * MDP(25) + (-MDP(21) * t437 + MDP(22) * t401 - MDP(24) * t234 - MDP(25) * t233) * qJD(5);];
tau = t1;
