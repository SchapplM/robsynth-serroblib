% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:10
% EndTime: 2019-03-09 02:29:15
% DurationCPUTime: 3.64s
% Computational Cost: add. (2407->389), mult. (4511->487), div. (0->0), fcn. (3032->10), ass. (0->169)
t370 = sin(qJ(5));
t374 = cos(qJ(5));
t375 = cos(qJ(4));
t432 = qJD(1) * t375;
t371 = sin(qJ(4));
t433 = qJD(1) * t371;
t299 = t370 * t433 - t374 * t432;
t357 = qJD(4) + qJD(5);
t369 = sin(qJ(6));
t373 = cos(qJ(6));
t282 = -t299 * t369 - t373 * t357;
t304 = t370 * t375 + t371 * t374;
t300 = t304 * qJD(1);
t468 = qJD(6) + t300;
t472 = t282 * t468;
t284 = -t299 * t373 + t357 * t369;
t471 = t284 * t468;
t421 = qJDD(1) * t375;
t446 = t300 * t357;
t470 = -t374 * t421 + t446;
t424 = qJD(1) * qJD(4);
t414 = t375 * t424;
t422 = qJDD(1) * t371;
t469 = t414 + t422;
t372 = sin(qJ(1));
t376 = cos(qJ(1));
t397 = g(1) * t376 + g(2) * t372;
t467 = t369 * t468;
t466 = g(1) * t372 - g(2) * t376;
t368 = pkin(1) + qJ(3);
t465 = qJD(1) * t368;
t360 = qJD(1) * qJD(2);
t361 = qJ(2) * qJDD(1);
t412 = qJDD(3) + t360 + t361;
t305 = -pkin(7) * qJDD(1) + t412;
t298 = t375 * t305;
t344 = qJ(2) * qJD(1) + qJD(3);
t317 = -pkin(7) * qJD(1) + t344;
t415 = t371 * t424;
t431 = qJD(4) * t371;
t270 = -t317 * t431 + qJDD(4) * pkin(4) + t298 + (t415 - t421) * pkin(8);
t430 = qJD(4) * t375;
t273 = -t469 * pkin(8) + t305 * t371 + t317 * t430;
t291 = -pkin(8) * t432 + t375 * t317;
t287 = qJD(4) * pkin(4) + t291;
t290 = -pkin(8) * t433 + t317 * t371;
t429 = qJD(5) * t370;
t464 = (qJD(5) * t287 + t273) * t374 + t270 * t370 - t290 * t429;
t318 = -qJD(2) + t465;
t367 = -pkin(7) + qJ(2);
t463 = qJD(4) * (qJD(2) + t318 + t465) + qJDD(4) * t367;
t448 = t290 * t374;
t269 = t287 * t370 + t448;
t462 = qJD(5) * t269 - t374 * t270 + t273 * t370;
t356 = qJDD(4) + qJDD(5);
t245 = -pkin(5) * t356 + t462;
t458 = pkin(8) - t367;
t310 = t458 * t371;
t311 = t458 * t375;
t280 = -t310 * t370 + t311 * t374;
t288 = qJD(2) * t375 + t431 * t458;
t289 = qJD(2) * t371 - qJD(4) * t311;
t252 = -qJD(5) * t280 + t288 * t370 + t289 * t374;
t416 = t371 * t429;
t439 = qJD(1) * t416 + t370 * t415;
t390 = t370 * t421 - t439;
t403 = t375 * t357;
t264 = (qJD(1) * t403 + t422) * t374 + t390;
t260 = qJDD(6) + t264;
t449 = t290 * t370;
t268 = t287 * t374 - t449;
t261 = -pkin(5) * t357 - t268;
t303 = t370 * t371 - t374 * t375;
t331 = t371 * pkin(4) + t368;
t275 = pkin(5) * t304 + pkin(9) * t303 + t331;
t428 = qJD(5) * t374;
t278 = -t370 * t430 - t371 * t428 - t374 * t431 - t375 * t429;
t281 = -t310 * t374 - t311 * t370;
t244 = pkin(9) * t356 + t464;
t302 = pkin(4) * t433 + t318;
t267 = pkin(5) * t300 + pkin(9) * t299 + t302;
t407 = qJD(6) * t267 + t244;
t461 = -t245 * t303 - t281 * t260 + t261 * t278 - (qJD(6) * t275 + t252) * t468 - t304 * t407;
t350 = 0.2e1 * t360;
t365 = qJ(4) + qJ(5);
t346 = sin(t365);
t334 = g(3) * t346;
t347 = cos(t365);
t335 = g(3) * t347;
t263 = -t370 * t422 - t470;
t426 = qJD(6) * t373;
t417 = t373 * t263 + t369 * t356 + t357 * t426;
t427 = qJD(6) * t369;
t250 = t299 * t427 + t417;
t457 = t250 * t369;
t456 = t260 * t369;
t455 = t260 * t373;
t454 = t261 * t300;
t453 = t261 * t303;
t452 = t275 * t260;
t451 = t282 * t299;
t450 = t284 * t299;
t447 = t299 * t357;
t445 = t303 * t373;
t444 = t369 * t372;
t443 = t369 * t376;
t442 = t372 * t373;
t441 = t373 * t376;
t440 = t278 * t357 - t303 * t356;
t438 = t376 * pkin(1) + t372 * qJ(2);
t364 = t375 ^ 2;
t436 = t371 ^ 2 - t364;
t377 = qJD(4) ^ 2;
t378 = qJD(1) ^ 2;
t435 = -t377 - t378;
t434 = MDP(22) * t374;
t320 = pkin(4) * t430 + qJD(3);
t359 = qJD(3) * qJD(1);
t423 = qJDD(1) * t368;
t419 = qJDD(4) * t371;
t366 = qJDD(1) * pkin(1);
t418 = t366 - qJDD(2);
t413 = qJDD(2) - t466;
t409 = t373 * t468;
t358 = qJDD(1) * qJ(3);
t400 = -t358 - t359 - t418;
t285 = t469 * pkin(4) - t400;
t247 = pkin(5) * t264 - pkin(9) * t263 + t285;
t262 = pkin(9) * t357 + t269;
t406 = qJD(6) * t262 - t247;
t276 = -pkin(5) * t299 + pkin(9) * t300;
t338 = pkin(4) * t370 + pkin(9);
t404 = pkin(4) * t432 + qJD(6) * t338 + t276;
t402 = MDP(23) * t357;
t401 = -t366 + t413;
t271 = t291 * t370 + t448;
t399 = pkin(4) * t429 - t271;
t272 = t291 * t374 - t449;
t398 = -pkin(4) * t428 + t272;
t248 = -t262 * t369 + t267 * t373;
t395 = t248 * t299 + t261 * t427 + t373 * t334;
t393 = -t338 * t260 + t454;
t279 = -t370 * t431 + t374 * t403 - t416;
t392 = -t279 * t357 - t304 * t356;
t391 = -t358 + t401;
t389 = t397 * t347;
t388 = t278 * t373 + t303 * t427;
t387 = t350 + 0.2e1 * t361 - t397;
t386 = qJD(1) * t318 + t397;
t385 = -t245 - t389;
t249 = t262 * t373 + t267 * t369;
t383 = -t249 * t299 + t261 * t426 + (t245 - t334) * t369 + (g(1) * t443 + g(2) * t444) * t347;
t382 = -t367 * t377 + t359 - t400 + t423 + t466;
t381 = t300 * t302 + t397 * t346 + t335 - t464;
t330 = t373 * t356;
t251 = qJD(6) * t284 + t263 * t369 - t330;
t380 = ((t250 - t472) * t373 + (-t251 - t471) * t369) * MDP(25) + (t284 * t409 + t457) * MDP(24) + (-t467 * t468 - t451 + t455) * MDP(27) + (t409 * t468 + t450 + t456) * MDP(26) + (t263 + t446) * MDP(19) + (-t447 + (-t357 * t432 - t422) * t374 - t390) * MDP(20) + (t299 ^ 2 - t300 ^ 2) * MDP(18) + t356 * MDP(21) + (-MDP(17) * t300 + MDP(28) * t468) * t299;
t379 = t299 * t302 + t334 - t389 - t462;
t349 = t376 * qJ(2);
t345 = qJDD(4) * t375;
t339 = -pkin(4) * t374 - pkin(5);
t295 = t346 * t441 - t444;
t294 = -t346 * t443 - t442;
t293 = -t346 * t442 - t443;
t292 = t346 * t444 - t441;
t256 = pkin(5) * t279 - pkin(9) * t278 + t320;
t253 = qJD(5) * t281 - t288 * t374 + t289 * t370;
t246 = t373 * t247;
t1 = [(-g(1) * t292 - g(2) * t294 - t249 * t279 + t280 * t250 + t253 * t284 + (-(-qJD(6) * t281 + t256) * t468 - t452 + t406 * t304 + qJD(6) * t453) * t369 + t461 * t373) * MDP(30) + (-g(1) * t293 - g(2) * t295 + t246 * t304 + t248 * t279 + t280 * t251 + t253 * t282 + (t256 * t468 + t452 + (-t262 * t304 - t281 * t468 - t453) * qJD(6)) * t373 + t461 * t369) * MDP(29) + (t260 * t304 + t279 * t468) * MDP(28) + (t250 * t304 - t260 * t445 + t279 * t284 + t388 * t468) * MDP(26) + (t303 * t456 - t251 * t304 - t279 * t282 + (-t278 * t369 + t303 * t426) * t468) * MDP(27) + (-t400 * t368 + t318 * qJD(3) + t412 * qJ(2) + t344 * qJD(2) - g(1) * (-t368 * t372 + t349) - g(2) * (qJ(3) * t376 + t438)) * MDP(9) + (-t253 * t357 + t264 * t331 + t279 * t302 - t280 * t356 + t285 * t304 + t300 * t320 + t346 * t466) * MDP(22) + (-t252 * t357 + t263 * t331 + t278 * t302 - t281 * t356 - t285 * t303 - t299 * t320 + t347 * t466) * MDP(23) + t466 * MDP(2) + (-t371 * t377 + t345) * MDP(12) + ((-t282 * t373 - t284 * t369) * t278 + (t457 + t251 * t373 + (-t282 * t369 + t284 * t373) * qJD(6)) * t303) * MDP(25) + t392 * MDP(20) + (t382 * t371 + t463 * t375) * MDP(15) + (-t463 * t371 + t382 * t375) * MDP(16) + qJDD(1) * MDP(1) + (qJDD(3) + t387) * MDP(7) + t387 * MDP(5) + (-0.2e1 * t366 + t413) * MDP(4) + (t418 * pkin(1) - g(1) * (-pkin(1) * t372 + t349) - g(2) * t438 + (t350 + t361) * qJ(2)) * MDP(6) + 0.2e1 * (-t371 * t421 + t424 * t436) * MDP(11) + t440 * MDP(19) + (qJDD(1) * t364 - 0.2e1 * t371 * t414) * MDP(10) + (0.2e1 * t359 - t391 + t423) * MDP(8) + (-t375 * t377 - t419) * MDP(13) + t397 * MDP(3) + (-t263 * t303 - t278 * t299) * MDP(17) + (-t263 * t304 + t264 * t303 - t278 * t300 + t279 * t299) * MDP(18) + (-t250 * t445 + t284 * t388) * MDP(24); t401 * MDP(6) + (-t359 + t391) * MDP(9) + (t439 + t447) * MDP(22) + t470 * MDP(23) + (-t451 - t455) * MDP(29) + (-t450 + t456) * MDP(30) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t378 + (MDP(29) * t467 + MDP(30) * t409) * t468 + (MDP(4) - MDP(8) + (-MDP(22) * t370 - MDP(16)) * t375 + (MDP(23) * t370 - MDP(15) - t434) * t371) * qJDD(1) + (-t344 * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + t374 * t402) * t371 + (-0.2e1 * qJD(4) * MDP(15) - t357 * t434 + t370 * t402) * t375) * qJD(1); qJDD(1) * MDP(7) - t378 * MDP(8) + (-t386 + t412) * MDP(9) + (t371 * t435 + t345) * MDP(15) + (t375 * t435 - t419) * MDP(16) + (-qJD(1) * t300 + t440) * MDP(22) + (qJD(1) * t299 + t392) * MDP(23) + (t251 * t303 - t278 * t282 - t304 * t456) * MDP(29) + (t250 * t303 - t278 * t284 - t304 * t455) * MDP(30) + ((-qJD(1) * t373 - t279 * t369 - t304 * t426) * MDP(29) + (qJD(1) * t369 - t279 * t373 + t304 * t427) * MDP(30)) * t468; (t272 * t357 + (t299 * t432 - t356 * t370 - t357 * t428) * pkin(4) + t381) * MDP(23) + (t339 * t251 + t399 * t282 + (t398 * t468 + t393) * t369 + (-t404 * t468 + t385) * t373 + t395) * MDP(29) + (t339 * t250 + t393 * t373 + t399 * t284 + (t369 * t404 + t373 * t398) * t468 + t383) * MDP(30) + qJDD(4) * MDP(14) + (g(3) * t375 + (-t305 + t386) * t371) * MDP(16) + MDP(12) * t421 + t380 - MDP(13) * t422 + (g(3) * t371 - t375 * t386 + t298) * MDP(15) + (t271 * t357 + (-t300 * t432 + t356 * t374 - t357 * t429) * pkin(4) + t379) * MDP(22) + (MDP(10) * t371 * t375 - MDP(11) * t436) * t378; (t268 * t357 + t381) * MDP(23) + (-pkin(5) * t251 - t269 * t282 + (-pkin(9) * t260 + t268 * t468 + t454) * t369 + ((-pkin(9) * qJD(6) - t276) * t468 + t385) * t373 + t395) * MDP(29) + (-pkin(5) * t250 + (t268 * t373 + t276 * t369) * t468 - t269 * t284 + t373 * t454 + (t427 * t468 - t455) * pkin(9) + t383) * MDP(30) + t380 + (t269 * t357 + t379) * MDP(22); t284 * t282 * MDP(24) + (-t282 ^ 2 + t284 ^ 2) * MDP(25) + (t417 + t472) * MDP(26) + (t330 + t471) * MDP(27) + t260 * MDP(28) + (-g(1) * t294 + g(2) * t292 + t249 * t468 - t261 * t284 + t246) * MDP(29) + (g(1) * t295 - g(2) * t293 + t248 * t468 + t261 * t282) * MDP(30) + ((-t244 + t335) * MDP(30) + (MDP(27) * t299 - MDP(29) * t262 - MDP(30) * t267) * qJD(6)) * t373 + (qJD(6) * t299 * MDP(26) + (-qJD(6) * t357 - t263) * MDP(27) + (-t407 + t335) * MDP(29) + t406 * MDP(30)) * t369;];
tau  = t1;
