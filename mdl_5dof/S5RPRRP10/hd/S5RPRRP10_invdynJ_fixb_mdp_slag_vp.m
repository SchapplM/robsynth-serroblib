% Calculate vector of inverse dynamics joint torques for
% S5RPRRP10
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:15:06
% EndTime: 2021-01-15 19:15:24
% DurationCPUTime: 5.13s
% Computational Cost: add. (3316->390), mult. (7807->474), div. (0->0), fcn. (5766->10), ass. (0->167)
t351 = cos(pkin(8));
t358 = cos(qJ(3));
t415 = t358 * t351;
t350 = sin(pkin(8));
t355 = sin(qJ(3));
t423 = t350 * t355;
t460 = t415 - t423;
t463 = t460 * qJD(1);
t353 = pkin(6) + qJ(2);
t327 = t353 * t350;
t321 = qJD(1) * t327;
t328 = t353 * t351;
t322 = qJD(1) * t328;
t290 = -t355 * t321 + t358 * t322;
t462 = qJD(3) * t290;
t320 = t350 * t358 + t351 * t355;
t316 = t320 * qJD(1);
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t369 = t320 * qJDD(1);
t404 = qJD(4) * t357;
t406 = qJD(3) * t354;
t260 = -t357 * qJDD(3) + t316 * t404 + t354 * t369 + (qJD(4) + t463) * t406;
t298 = t316 * t357 + t406;
t337 = t351 * pkin(2) + pkin(1);
t326 = -t337 * qJD(1) + qJD(2);
t266 = -pkin(3) * t463 - pkin(7) * t316 + t326;
t284 = qJD(3) * pkin(7) + t290;
t249 = t266 * t354 + t284 * t357;
t318 = t320 * qJD(3);
t334 = qJDD(1) * t415;
t386 = -qJDD(1) * t423 + t334;
t288 = qJD(1) * t318 - t386;
t317 = t460 * qJD(3);
t368 = qJD(1) * t317;
t455 = -pkin(7) * t320 - t337;
t257 = t288 * pkin(3) - pkin(7) * t368 + qJDD(1) * t455 + qJDD(2);
t253 = t357 * t257;
t402 = qJD(1) * qJD(2);
t450 = t353 * qJDD(1) + t402;
t301 = t450 * t350;
t302 = t450 * t351;
t385 = -t301 * t355 + t302 * t358;
t456 = -t321 * t358 - t355 * t322;
t254 = qJDD(3) * pkin(7) + qJD(3) * t456 + t385;
t366 = -t249 * qJD(4) - t254 * t354 + t253;
t363 = t369 + t368;
t403 = t357 * qJD(3);
t405 = qJD(4) * t354;
t259 = -qJD(4) * t403 - t354 * qJDD(3) + t316 * t405 - t357 * t363;
t437 = qJ(5) * t259;
t282 = qJDD(4) + t288;
t448 = pkin(4) * t282;
t235 = -qJD(5) * t298 + t366 + t437 + t448;
t296 = t316 * t354 - t403;
t244 = -qJ(5) * t296 + t249;
t305 = qJD(4) - t463;
t434 = t244 * t305;
t461 = t235 + t434;
t356 = sin(qJ(1));
t359 = cos(qJ(1));
t454 = g(1) * t356 - g(2) * t359;
t458 = qJDD(2) - t454;
t349 = pkin(8) + qJ(3);
t341 = sin(t349);
t457 = t454 * t341;
t388 = g(1) * t359 + g(2) * t356;
t375 = t388 * t341;
t453 = qJ(2) * qJDD(1);
t342 = cos(t349);
t416 = t357 * t359;
t419 = t354 * t356;
t306 = t342 * t419 + t416;
t417 = t356 * t357;
t418 = t354 * t359;
t308 = -t342 * t418 + t417;
t438 = g(3) * t354;
t452 = -g(1) * t308 + g(2) * t306 + t341 * t438;
t440 = g(3) * t341;
t451 = t388 * t342 + t440;
t449 = t298 ^ 2;
t447 = pkin(4) * t296;
t439 = g(3) * t342;
t352 = -qJ(5) - pkin(7);
t436 = qJ(5) * t260;
t435 = qJDD(1) * pkin(1);
t433 = t259 * t354;
t432 = t296 * t305;
t431 = t296 * t463;
t430 = t296 * t316;
t429 = t298 * t305;
t428 = t298 * t316;
t427 = t463 * t354;
t426 = t320 * t354;
t425 = t320 * t357;
t422 = t351 * MDP(4);
t421 = t354 * t282;
t420 = t354 * t305;
t272 = t357 * t282;
t294 = -t327 * t355 + t328 * t358;
t291 = t357 * t294;
t248 = t357 * t266 - t284 * t354;
t243 = -qJ(5) * t298 + t248;
t242 = pkin(4) * t305 + t243;
t414 = -t243 + t242;
t413 = -t354 * t260 - t296 * t404;
t285 = pkin(3) * t316 - pkin(7) * t463;
t412 = t354 * t285 + t357 * t456;
t287 = -pkin(3) * t460 + t455;
t411 = t354 * t287 + t291;
t396 = qJD(4) * t352;
t410 = qJ(5) * t427 + qJD(5) * t357 + t354 * t396 - t412;
t275 = t357 * t285;
t409 = -pkin(4) * t316 - t275 + (qJ(5) * t463 + t396) * t357 + (-qJD(5) + t456) * t354;
t408 = (g(1) * t416 + g(2) * t417) * t341;
t407 = t350 ^ 2 + t351 ^ 2;
t383 = -t327 * t358 - t328 * t355;
t267 = qJD(2) * t460 + t383 * qJD(3);
t286 = pkin(3) * t318 - pkin(7) * t317;
t401 = t357 * t267 + t354 * t286 + t287 * t404;
t400 = t320 * t404;
t399 = pkin(4) * t354 + t353;
t377 = -t301 * t358 - t355 * t302 - t462;
t255 = -qJDD(3) * pkin(3) - t377;
t239 = pkin(4) * t260 + qJDD(5) + t255;
t397 = -t239 - t439;
t395 = -qJD(5) - t447;
t394 = t305 * t357;
t393 = t357 * t254 + t354 * t257 + t266 * t404 - t284 * t405;
t392 = 0.2e1 * t407;
t391 = qJD(4) * pkin(7) * t305 + t255;
t390 = -g(1) * t306 - g(2) * t308;
t307 = -t342 * t417 + t418;
t309 = t342 * t416 + t419;
t389 = -g(1) * t307 - g(2) * t309;
t339 = pkin(4) * t357 + pkin(3);
t382 = t339 * t342 - t341 * t352;
t236 = -qJD(5) * t296 + t393 - t436;
t380 = -t305 * t242 + t236;
t379 = -qJ(5) * t317 - qJD(5) * t320;
t378 = t435 - t458;
t376 = -t305 * t405 + t420 * t463 + t272;
t283 = -qJD(3) * pkin(3) - t456;
t374 = t337 + t382;
t373 = t317 * t354 + t400;
t372 = t317 * t357 - t320 * t405;
t371 = -pkin(7) * t282 + t305 * t283;
t367 = g(1) * t309 - g(2) * t307 + t357 * t440 - t393;
t365 = t392 * t402 - t388;
t268 = t320 * qJD(2) + t294 * qJD(3);
t362 = t366 + t452;
t333 = t342 * t438;
t330 = t352 * t357;
t329 = t352 * t354;
t325 = -t337 * qJDD(1) + qJDD(2);
t295 = t296 ^ 2;
t278 = t357 * t287;
t276 = t357 * t286;
t269 = pkin(4) * t426 - t383;
t262 = pkin(4) * t427 + t290;
t261 = t283 - t395;
t256 = t373 * pkin(4) + t268;
t250 = -qJ(5) * t426 + t411;
t246 = -pkin(4) * t460 - qJ(5) * t425 - t294 * t354 + t278;
t238 = -qJ(5) * t400 + (-qJD(4) * t294 + t379) * t354 + t401;
t237 = pkin(4) * t318 - t267 * t354 + t276 + t379 * t357 + (-t291 + (qJ(5) * t320 - t287) * t354) * qJD(4);
t1 = [qJDD(1) * MDP(1) + t454 * MDP(2) + t388 * MDP(3) + (t378 + t435) * t422 + (t392 * t453 + t365) * MDP(5) + (t378 * pkin(1) + (t407 * t453 + t365) * qJ(2)) * MDP(6) + (t316 * t317 + t363 * t320) * MDP(7) + (-t320 * t288 - t316 * t318 + t317 * t463 + t363 * t460) * MDP(8) + (qJD(3) * t317 + qJDD(3) * t320) * MDP(9) + (-qJD(3) * t318 + qJDD(3) * t460) * MDP(10) + (-qJD(3) * t268 + qJDD(3) * t383 - t288 * t337 + t318 * t326 - t325 * t460 + t342 * t454) * MDP(12) + (-t294 * qJDD(3) + t326 * t317 + t325 * t320 - t457 - t337 * t369 + (-t337 * t463 - t267) * qJD(3)) * MDP(13) + (-t259 * t425 + t372 * t298) * MDP(14) + ((-t296 * t357 - t298 * t354) * t317 + (t433 - t260 * t357 + (t296 * t354 - t298 * t357) * qJD(4)) * t320) * MDP(15) + (t259 * t460 + t320 * t272 + t298 * t318 + t372 * t305) * MDP(16) + (t260 * t460 - t296 * t318 - t373 * t305 - t320 * t421) * MDP(17) + (-t282 * t460 + t305 * t318) * MDP(18) + ((-t294 * t404 + t276) * t305 + t278 * t282 - (-t284 * t404 + t253) * t460 + t248 * t318 + t268 * t296 - t383 * t260 + t283 * t400 + ((-qJD(4) * t287 - t267) * t305 - t294 * t282 - (-qJD(4) * t266 - t254) * t460 + t255 * t320 + t283 * t317) * t354 + t389) * MDP(19) + (-(-t294 * t405 + t401) * t305 - t411 * t282 + t393 * t460 - t249 * t318 + t268 * t298 + t383 * t259 + t255 * t425 + t372 * t283 + t390) * MDP(20) + (-t235 * t460 + t237 * t305 + t239 * t426 + t242 * t318 + t246 * t282 + t256 * t296 + t260 * t269 + t373 * t261 + t389) * MDP(21) + (t236 * t460 - t238 * t305 + t239 * t425 - t244 * t318 - t250 * t282 + t256 * t298 - t259 * t269 + t372 * t261 + t390) * MDP(22) + (-t237 * t298 - t238 * t296 + t246 * t259 - t250 * t260 + t457 + (-t242 * t357 - t244 * t354) * t317 + (-t235 * t357 - t236 * t354 + (t242 * t354 - t244 * t357) * qJD(4)) * t320) * MDP(23) + (t235 * t246 + t236 * t250 + t242 * t237 + t244 * t238 + t239 * t269 + t261 * t256 + (-g(1) * t399 - g(2) * t374) * t359 + (g(1) * t374 - g(2) * t399) * t356) * MDP(24); t458 * MDP(6) - t334 * MDP(12) + t413 * MDP(23) + (-t261 * t316 - t454) * MDP(24) + (MDP(20) + MDP(22)) * (-t305 ^ 2 * t357 - t421 - t428) + (MDP(19) + MDP(21)) * (t376 - t430) + ((t259 + t431) * MDP(23) + t461 * MDP(24)) * t357 + (MDP(23) * t429 + t380 * MDP(24)) * t354 + (MDP(12) * t423 + t320 * MDP(13) - pkin(1) * MDP(6) - t422) * qJDD(1) + (t316 * MDP(12) + t463 * MDP(13) + (t320 * MDP(12) + MDP(13) * t460) * qJD(1)) * qJD(3) + (-MDP(6) * qJ(2) - MDP(5)) * qJD(1) ^ 2 * t407; -t463 ^ 2 * MDP(8) + t369 * MDP(9) + t386 * MDP(10) + qJDD(3) * MDP(11) + (t375 + t377 - t439 + t462) * MDP(12) + (-t326 * t463 - t385 + t451) * MDP(13) + (t298 * t394 - t433) * MDP(14) + ((-t259 + t431) * t357 - t298 * t420 + t413) * MDP(15) + (t305 * t394 + t421 - t428) * MDP(16) + (t376 + t430) * MDP(17) + (-pkin(3) * t260 - t275 * t305 - t290 * t296 + (-t391 - t439) * t357 + (t305 * t456 + t371) * t354 + t408) * MDP(19) + (pkin(3) * t259 + t412 * t305 - t290 * t298 + t333 + t371 * t357 + (-t375 + t391) * t354) * MDP(20) + (-t260 * t339 - t262 * t296 + t282 * t329 + t397 * t357 + t409 * t305 + (-t261 * t463 + (t261 + t447) * qJD(4)) * t354 + t408) * MDP(21) + (t259 * t339 - t262 * t298 + t282 * t330 + t333 - t410 * t305 + t261 * t394 + (pkin(4) * qJD(4) * t298 + t239 - t375) * t354) * MDP(22) + (t259 * t329 + t260 * t330 - t410 * t296 - t409 * t298 - t354 * t461 + t380 * t357 - t451) * MDP(23) + (-t236 * t330 + t235 * t329 - t239 * t339 - g(3) * t382 + (pkin(4) * t405 - t262) * t261 + t410 * t244 + t409 * t242 + t388 * (t339 * t341 + t342 * t352)) * MDP(24) + (-t326 * MDP(12) - t305 * MDP(18) - t248 * MDP(19) + t249 * MDP(20) - t242 * MDP(21) + t244 * MDP(22) - MDP(7) * t463 + MDP(8) * t316) * t316; t298 * t296 * MDP(14) + (-t295 + t449) * MDP(15) + (-t259 + t432) * MDP(16) + (-t260 + t429) * MDP(17) + t282 * MDP(18) + (t249 * t305 - t283 * t298 + t362) * MDP(19) + (t248 * t305 + t283 * t296 + t367) * MDP(20) + (0.2e1 * t448 + t437 + t434 + (-t261 + t395) * t298 + t362) * MDP(21) + (-pkin(4) * t449 + t436 + t243 * t305 + (qJD(5) + t261) * t296 + t367) * MDP(22) + (pkin(4) * t259 - t414 * t296) * MDP(23) + (t414 * t244 + (-t261 * t298 + t235 + t452) * pkin(4)) * MDP(24); (t260 + t429) * MDP(21) + (-t259 - t432) * MDP(22) + (-t295 - t449) * MDP(23) + (t242 * t298 + t244 * t296 - t375 - t397) * MDP(24);];
tau = t1;
