% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:19
% EndTime: 2019-03-09 01:56:24
% DurationCPUTime: 4.17s
% Computational Cost: add. (2425->404), mult. (4731->476), div. (0->0), fcn. (3318->10), ass. (0->177)
t351 = sin(qJ(4));
t346 = sin(pkin(9));
t419 = qJD(1) * t346;
t322 = t351 * t419;
t347 = cos(pkin(9));
t354 = cos(qJ(4));
t426 = t354 * t347;
t402 = qJD(1) * t426;
t305 = -t322 + t402;
t297 = qJD(6) + t305;
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t464 = g(1) * t352 - g(2) * t355;
t342 = pkin(9) + qJ(4);
t331 = sin(t342);
t332 = cos(t342);
t364 = -g(3) * t332 - t464 * t331;
t349 = -pkin(1) - qJ(3);
t456 = -qJD(1) * qJD(3) + qJDD(1) * t349;
t313 = qJDD(2) + t456;
t395 = -pkin(7) * qJDD(1) + t313;
t285 = t395 * t346;
t286 = t395 * t347;
t320 = qJD(1) * t349 + qJD(2);
t398 = -pkin(7) * qJD(1) + t320;
t293 = t398 * t346;
t294 = t398 * t347;
t415 = qJD(4) * t354;
t416 = qJD(4) * t351;
t390 = -t354 * t285 - t351 * t286 + t293 * t416 - t294 * t415;
t467 = t364 - t390;
t321 = qJDD(1) * t426;
t407 = qJDD(1) * t346;
t382 = t351 * t407 - t321;
t432 = t346 * t354;
t312 = t347 * t351 + t432;
t303 = t312 * qJD(1);
t417 = qJD(4) * t303;
t271 = t382 + t417;
t425 = -t351 * t293 + t354 * t294;
t466 = qJD(5) - t425;
t405 = MDP(16) - MDP(19);
t350 = sin(qJ(6));
t394 = t297 * t350;
t421 = t346 ^ 2 + t347 ^ 2;
t465 = t320 * t421;
t343 = qJDD(1) * qJ(2);
t344 = qJD(1) * qJD(2);
t463 = t343 + t344;
t317 = qJDD(3) + t463;
t388 = g(1) * t355 + g(2) * t352;
t367 = -t388 + t317;
t461 = t346 * MDP(7) + t347 * MDP(8);
t410 = pkin(5) * t305 + t466;
t404 = -MDP(17) + MDP(20);
t264 = t354 * t293 + t351 * t294;
t259 = -qJD(4) * qJ(5) - t264;
t450 = pkin(5) * t303;
t249 = -t259 - t450;
t270 = -qJDD(6) + t271;
t452 = pkin(4) + pkin(8);
t460 = t452 * t270 + (t249 - t264 + t450) * t297;
t447 = -pkin(7) + t349;
t314 = t447 * t346;
t315 = t447 * t347;
t274 = t314 * t354 + t315 * t351;
t311 = t346 * t351 - t426;
t258 = -qJD(3) * t311 + qJD(4) * t274;
t273 = t314 * t351 - t315 * t354;
t459 = -qJD(4) * t258 - qJDD(4) * t273 - t331 * t388;
t257 = qJD(3) * t312 + t314 * t416 - t315 * t415;
t458 = qJD(4) * t257 - qJDD(4) * t274 - t332 * t388;
t318 = qJD(4) * t322;
t457 = -t312 * qJDD(1) + t318;
t455 = t303 ^ 2;
t454 = t305 ^ 2;
t453 = 0.2e1 * t344;
t401 = t347 * t415;
t272 = qJD(1) * t401 - t457;
t451 = pkin(4) * t272;
t449 = g(3) * t331;
t334 = t346 * pkin(3);
t446 = pkin(1) * qJDD(1);
t445 = qJ(5) * t271;
t444 = qJ(5) * t303;
t443 = qJDD(4) * pkin(4);
t353 = cos(qJ(6));
t412 = qJD(6) * t353;
t403 = t353 * qJDD(4) + t350 * t272 + t303 * t412;
t408 = qJD(4) * qJD(6);
t244 = -t350 * t408 + t403;
t442 = t244 * t353;
t326 = qJ(2) + t334;
t386 = qJ(5) * t311 + t326;
t255 = t312 * t452 + t386;
t441 = t255 * t270;
t440 = t270 * t350;
t278 = qJD(4) * t350 - t353 * t303;
t439 = t278 * t297;
t438 = t278 * t303;
t280 = qJD(4) * t353 + t303 * t350;
t437 = t280 * t297;
t436 = t280 * t303;
t435 = t303 * t305;
t434 = t312 * t350;
t430 = t350 * t352;
t429 = t350 * t355;
t428 = t352 * t353;
t265 = t353 * t270;
t427 = t353 * t355;
t307 = -t346 * t415 - t347 * t416;
t424 = t307 * qJD(4) - t311 * qJDD(4);
t423 = t355 * pkin(1) + t352 * qJ(2);
t418 = qJD(4) * t264;
t330 = qJD(1) * qJ(2) + qJD(3);
t316 = pkin(3) * t419 + t330;
t374 = -qJ(5) * t305 + t316;
t250 = t303 * t452 + t374;
t414 = qJD(6) * t250;
t413 = qJD(6) * t350;
t406 = qJDD(4) * qJ(5);
t400 = g(2) * t423;
t399 = qJDD(2) - t464;
t397 = t421 * MDP(9);
t396 = t421 * t313;
t393 = t297 * t353;
t328 = pkin(3) * t407;
t310 = t328 + t317;
t363 = -qJD(5) * t305 + t310 + t445;
t237 = t272 * t452 + t363;
t246 = -qJD(4) * t452 + t410;
t392 = qJD(6) * t246 + t237;
t391 = MDP(27) * t297;
t383 = -pkin(4) * t331 + qJ(5) * t332;
t381 = -t414 - t449;
t239 = t246 * t350 + t250 * t353;
t379 = t271 * t311 + t305 * t307;
t308 = -t346 * t416 + t401;
t376 = -qJD(4) * t308 - qJDD(4) * t312;
t375 = t351 * t285 - t286 * t354 + t293 * t415 + t294 * t416;
t373 = t308 * t350 + t312 * t412;
t372 = -qJ(5) * t307 + qJD(5) * t311 + qJD(2);
t370 = qJDD(5) + t375;
t369 = -t383 + t334;
t240 = -qJD(4) * qJD(5) + t390 - t406;
t236 = -pkin(5) * t272 - t240;
t261 = -pkin(5) * t311 + t273;
t368 = t236 * t312 + t249 * t308 + t261 * t270;
t365 = t351 * t404 + t354 * t405;
t361 = -t464 * t332 - t375 + t449;
t242 = t370 - t443;
t256 = -qJD(4) * pkin(4) + t466;
t360 = -t240 * t312 + t242 * t311 - t256 * t307 - t259 * t308 - t464;
t359 = t236 + (t297 * t452 + t444) * t297 + t364;
t260 = pkin(4) * t303 + t374;
t358 = t260 * t305 + qJDD(5) - t361;
t357 = qJD(1) ^ 2;
t348 = -pkin(7) - qJ(3);
t336 = t355 * qJ(2);
t301 = -t332 * t430 + t427;
t300 = -t332 * t428 - t429;
t299 = -t332 * t429 - t428;
t298 = -t332 * t427 + t430;
t269 = pkin(4) * t305 + t444;
t268 = pkin(4) * t312 + t386;
t267 = t353 * t272;
t262 = -pkin(5) * t312 + t274;
t254 = pkin(4) * t308 + t372;
t248 = pkin(5) * t307 + t258;
t247 = -pkin(5) * t308 - t257;
t245 = t280 * qJD(6) + qJDD(4) * t350 - t267;
t243 = t308 * t452 + t372;
t241 = t363 + t451;
t238 = t246 * t353 - t250 * t350;
t235 = -pkin(5) * t271 - qJDD(4) * t452 + t370;
t234 = t353 * t235;
t1 = [(-(qJDD(2) - t446) * pkin(1) - g(1) * (-pkin(1) * t352 + t336) - t400 + (t343 + t453) * qJ(2)) * MDP(6) + (0.2e1 * t343 + t453 - t388) * MDP(5) + ((-t278 * t350 + t280 * t353) * t308 + (t442 - t245 * t350 + (-t278 * t353 - t280 * t350) * qJD(6)) * t312) * MDP(23) + (t399 - 0.2e1 * t446) * MDP(4) + (-g(1) * t299 - g(2) * t301 - t234 * t311 + t238 * t307 + t262 * t245 + t247 * t278 + (t237 * t311 - t243 * t297 + t441) * t350 + (t248 * t297 - t368) * t353 + ((-t255 * t353 - t261 * t350) * t297 + t239 * t311 + t249 * t434) * qJD(6)) * MDP(27) + (-g(1) * t298 - g(2) * t300 - t239 * t307 + t262 * t244 + t247 * t280 + (-(qJD(6) * t261 + t243) * t297 + t441 + t392 * t311 + t249 * qJD(6) * t312) * t353 + (-(-qJD(6) * t255 + t248) * t297 + (t235 - t414) * t311 + t368) * t350) * MDP(28) + (-t244 * t311 - t270 * t434 + t280 * t307 + t297 * t373) * MDP(24) + (t244 * t434 + t280 * t373) * MDP(22) + (-t312 * t265 + t245 * t311 - t278 * t307 + (t308 * t353 - t312 * t413) * t297) * MDP(25) + t424 * MDP(13) + (t241 * t268 + t260 * t254 - t240 * t274 + t259 * t257 + t242 * t273 + t256 * t258 - g(1) * t336 - t400 + (-g(1) * t369 + g(2) * t348) * t355 + (-g(1) * (-pkin(1) + t348) - g(2) * t369) * t352) * MDP(21) + t388 * MDP(3) + t379 * MDP(11) + t376 * MDP(14) + t464 * MDP(2) + (t464 + t421 * (-t313 - t456)) * MDP(9) + (t257 * t303 + t258 * t305 - t271 * t273 - t272 * t274 - t360) * MDP(18) + (qJD(2) * t305 - t271 * t326 + t307 * t316 - t310 * t311 + t458) * MDP(17) + (t241 * t311 - t254 * t305 - t260 * t307 + t268 * t271 - t458) * MDP(20) + (qJD(2) * t303 + t272 * t326 + t308 * t316 + t310 * t312 + t459) * MDP(16) + (-t241 * t312 - t254 * t303 - t260 * t308 - t268 * t272 - t459) * MDP(19) + t461 * (t367 + t463) + (t317 * qJ(2) + t330 * qJD(2) - g(1) * (t349 * t352 + t336) - g(2) * (qJ(3) * t355 + t423) + t349 * t396 - qJD(3) * t465) * MDP(10) + qJDD(1) * MDP(1) + (t270 * t311 + t297 * t307) * MDP(26) + (t271 * t312 + t272 * t311 - t303 * t307 - t305 * t308) * MDP(12); t399 * MDP(6) + (-qJD(1) * t330 + t396 - t464) * MDP(10) + (-t272 * t312 - t303 * t308 - t379) * MDP(18) + (-qJD(1) * t260 + t360) * MDP(21) + (t245 * t312 - t265 * t311 + t278 * t308) * MDP(27) + (t244 * t312 + t280 * t308 + t311 * t440) * MDP(28) + (-MDP(6) * qJ(2) - MDP(5) - t461) * t357 + (-pkin(1) * MDP(6) + MDP(4) - t397) * qJDD(1) + ((qJD(1) * t350 - t307 * t353 - t311 * t413) * MDP(27) + (qJD(1) * t353 + t307 * t350 - t311 * t412) * MDP(28)) * t297 + t404 * (qJD(1) * t305 - t376) - t405 * (qJD(1) * t303 - t424); (qJD(1) * t465 + t367) * MDP(10) + t321 * MDP(17) + (-t454 - t455) * MDP(18) + (t417 - t321) * MDP(20) + (t451 + t445 - t259 * t303 + t328 + (-qJD(5) - t256) * t305 + t367) * MDP(21) + (t438 + t440) * MDP(27) + (t265 + t436) * MDP(28) - t357 * t397 - t405 * t318 + (MDP(28) * t394 - t353 * t391) * t297 + ((t351 * t405 + MDP(8)) * t347 + (MDP(7) + t365) * t346) * qJDD(1) + (-t303 * MDP(17) + t405 * t305 + (t347 * t365 + t404 * t432) * qJD(1)) * qJD(4); MDP(11) * t435 + (t454 - t455) * MDP(12) - t382 * MDP(13) + ((t305 - t402) * qJD(4) + t457) * MDP(14) + qJDD(4) * MDP(15) + (-t305 * t316 + t361 + t418) * MDP(16) + (qJD(4) * t425 + t303 * t316 - t467) * MDP(17) + (pkin(4) * t271 - qJ(5) * t272 + (-t259 - t264) * t305 + (t256 - t466) * t303) * MDP(18) + (t269 * t303 + t358 - t418 - 0.2e1 * t443) * MDP(19) + (0.2e1 * t406 - t260 * t303 + t269 * t305 + (0.2e1 * qJD(5) - t425) * qJD(4) + t467) * MDP(20) + (-t242 * pkin(4) - g(3) * t383 - t240 * qJ(5) - t256 * t264 - t259 * t466 - t260 * t269 - t464 * (pkin(4) * t332 + qJ(5) * t331)) * MDP(21) + (-t280 * t394 + t442) * MDP(22) + ((-t245 - t437) * t353 + (-t244 + t439) * t350) * MDP(23) + (-t297 * t394 - t265 + t436) * MDP(24) + (-t297 * t393 - t438 + t440) * MDP(25) + t297 * t303 * MDP(26) + (qJ(5) * t245 + t238 * t303 + t410 * t278 + t359 * t350 + t460 * t353) * MDP(27) + (qJ(5) * t244 - t239 * t303 + t410 * t280 - t460 * t350 + t359 * t353) * MDP(28); (t417 - t271) * MDP(18) + (qJDD(4) - t435) * MDP(19) + (-qJD(4) ^ 2 - t454) * MDP(20) + (qJD(4) * t259 + t358 - t443) * MDP(21) + (-qJD(4) * t278 - t265) * MDP(27) + (-qJD(4) * t280 + t440) * MDP(28) + (-MDP(28) * t393 - t350 * t391) * t297; t280 * t278 * MDP(22) + (-t278 ^ 2 + t280 ^ 2) * MDP(23) + (t403 + t439) * MDP(24) + (t267 + t437) * MDP(25) - t270 * MDP(26) + (-g(1) * t300 + g(2) * t298 + t239 * t297 - t249 * t280 + t234) * MDP(27) + (g(1) * t301 - g(2) * t299 + t238 * t297 + t249 * t278) * MDP(28) + (-MDP(25) * t408 + MDP(27) * t381 - MDP(28) * t392) * t353 + (-MDP(24) * t408 + (-qJD(6) * t303 - qJDD(4)) * MDP(25) - t392 * MDP(27) + (-t235 - t381) * MDP(28)) * t350;];
tau  = t1;
