% Calculate vector of inverse dynamics joint torques for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP11_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:14
% EndTime: 2019-12-31 20:14:19
% DurationCPUTime: 4.26s
% Computational Cost: add. (2239->440), mult. (4531->533), div. (0->0), fcn. (2570->6), ass. (0->184)
t350 = sin(qJ(2));
t413 = qJD(1) * qJD(2);
t401 = t350 * t413;
t353 = cos(qJ(2));
t411 = qJDD(1) * t353;
t469 = -t401 + t411;
t459 = pkin(3) + pkin(6);
t349 = sin(qJ(4));
t352 = cos(qJ(4));
t427 = qJD(1) * t353;
t403 = t349 * t427;
t259 = -qJD(4) * t403 + qJDD(2) * t349 + (qJD(2) * qJD(4) + t469) * t352;
t424 = qJD(2) * t352;
t297 = -t403 + t424;
t428 = qJD(1) * t350;
t321 = qJD(4) + t428;
t443 = t297 * t321;
t468 = -t259 + t443;
t467 = t321 ^ 2;
t426 = qJD(2) * t349;
t295 = t352 * t427 + t426;
t420 = qJD(4) * t295;
t258 = -t352 * qJDD(2) + t349 * t469 + t420;
t425 = qJD(2) * t350;
t302 = t459 * t425;
t328 = pkin(6) * t411;
t344 = qJDD(2) * qJ(3);
t345 = qJD(2) * qJD(3);
t406 = t328 + t344 + t345;
t264 = pkin(3) * t411 - qJD(1) * t302 + t406;
t241 = pkin(4) * t259 + qJ(5) * t258 - qJD(5) * t297 + t264;
t351 = sin(qJ(1));
t354 = cos(qJ(1));
t390 = g(1) * t354 + g(2) * t351;
t452 = g(3) * t350;
t364 = -t390 * t353 - t452;
t460 = pkin(2) + pkin(7);
t416 = qJD(4) * t460;
t359 = t321 * t416 + t364;
t466 = t241 + t359;
t336 = t350 * qJ(3);
t340 = t353 * pkin(2);
t430 = t340 + t336;
t307 = -pkin(1) - t430;
t457 = pkin(7) * t353;
t292 = t307 - t457;
t309 = t459 * t350;
t431 = t352 * t292 + t349 * t309;
t400 = t353 * t413;
t412 = qJDD(1) * t350;
t371 = t400 + t412;
t294 = qJDD(4) + t371;
t283 = t352 * t294;
t419 = qJD(4) * t349;
t465 = -t321 * t419 + t283;
t330 = pkin(6) * t428;
t464 = qJD(3) + t330;
t463 = -MDP(20) - MDP(22);
t410 = -MDP(21) + MDP(24);
t438 = t351 * t352;
t440 = t350 * t354;
t287 = t349 * t440 + t438;
t436 = t352 * t354;
t441 = t350 * t351;
t289 = -t349 * t441 + t436;
t343 = g(3) * t353;
t320 = pkin(2) * t401;
t450 = qJ(3) * t353;
t385 = pkin(7) * t350 - t450;
t422 = qJD(3) * t350;
t362 = qJD(2) * t385 - t422;
t398 = -pkin(1) - t336;
t369 = -t353 * t460 + t398;
t250 = qJD(1) * t362 + qJDD(1) * t369 + t320;
t319 = pkin(6) * t400;
t327 = pkin(6) * t412;
t399 = qJDD(3) + t319 + t327;
t263 = pkin(3) * t371 - qJDD(2) * t460 + t399;
t415 = pkin(3) * t428 + t464;
t277 = -qJD(2) * t460 + t415;
t418 = qJD(4) * t352;
t408 = t250 * t352 + t263 * t349 + t277 * t418;
t274 = t369 * qJD(1);
t421 = qJD(4) * t274;
t462 = -g(1) * t287 + g(2) * t289 + (-t421 + t343) * t349 + t408;
t461 = t297 ^ 2;
t458 = pkin(4) * t294;
t456 = g(1) * t351;
t453 = g(2) * t354;
t451 = pkin(6) * qJDD(2);
t449 = qJ(5) * t294;
t448 = qJDD(2) * pkin(2);
t252 = t274 * t352 + t277 * t349;
t447 = t252 * t321;
t446 = t258 * t352;
t445 = t295 * t297;
t444 = t295 * t321;
t442 = t321 * t350;
t357 = qJD(1) ^ 2;
t439 = t350 * t357;
t437 = t351 * t353;
t435 = t353 * t354;
t434 = t460 * t294;
t334 = pkin(2) * t428;
t280 = qJD(1) * t385 + t334;
t331 = pkin(6) * t427;
t303 = pkin(3) * t427 + t331;
t433 = t280 * t352 + t303 * t349;
t387 = pkin(4) * t352 + qJ(5) * t349;
t378 = -pkin(3) - t387;
t432 = qJD(4) * t387 - qJD(5) * t352 - t378 * t428 + t464;
t310 = t459 * t353;
t347 = t350 ^ 2;
t348 = t353 ^ 2;
t429 = t347 - t348;
t346 = qJD(2) * qJ(3);
t423 = qJD(2) * t353;
t417 = qJD(4) * t353;
t251 = -t274 * t349 + t277 * t352;
t414 = qJD(5) - t251;
t409 = t353 * t439;
t407 = g(1) * t440 + g(2) * t441 - t343;
t284 = t346 + t303;
t402 = g(3) * t430;
t397 = -qJD(2) * pkin(2) + qJD(3);
t396 = t250 * t349 - t263 * t352 + t274 * t418 + t277 * t419;
t395 = pkin(1) * t354 + pkin(2) * t435 + pkin(6) * t351 + qJ(3) * t440;
t394 = -t327 + t407;
t356 = qJD(2) ^ 2;
t393 = pkin(6) * t356 + t453;
t286 = t349 * t351 - t350 * t436;
t288 = t349 * t354 + t350 * t438;
t392 = -g(1) * t288 - g(2) * t286;
t391 = -g(1) * t289 - g(2) * t287;
t239 = qJD(5) * t321 - t274 * t419 + t408 + t449;
t244 = -pkin(4) * t321 + t414;
t389 = -t244 * t428 - t239;
t240 = qJDD(5) + t396 - t458;
t246 = qJ(5) * t321 + t252;
t388 = -t246 * t428 + t240;
t386 = pkin(4) * t349 - qJ(5) * t352;
t382 = t244 * t349 + t246 * t352;
t305 = t330 + t397;
t308 = -t331 - t346;
t381 = t305 * t353 + t308 * t350;
t379 = t398 - t340;
t375 = -0.2e1 * pkin(1) * t413 - t451;
t285 = t379 * qJD(1);
t374 = t285 * t428 + qJDD(3) - t394;
t373 = -t349 * t294 - t321 * t418;
t372 = -qJ(3) * t423 - t422;
t333 = pkin(2) * t425;
t272 = t333 + t362;
t304 = t459 * t423;
t370 = t272 * t352 - t292 * t419 + t304 * t349 + t309 * t418;
t368 = 0.2e1 * qJDD(1) * pkin(1) - t393;
t366 = t451 + (-qJD(1) * t307 - t285) * qJD(2);
t253 = pkin(4) * t295 - qJ(5) * t297 + t284;
t365 = t253 * t321 - t434;
t363 = g(1) * t286 - g(2) * t288 + t343 * t352 - t396;
t260 = qJD(1) * t372 + qJDD(1) * t379 + t320;
t282 = t333 + t372;
t361 = qJD(1) * t282 + qJDD(1) * t307 + t260 + t393;
t275 = pkin(6) * t401 - t406;
t279 = t399 - t448;
t360 = qJD(2) * t381 - t275 * t353 + t279 * t350;
t358 = t253 * t297 + qJDD(5) - t363;
t341 = t354 * pkin(6);
t325 = g(1) * t437;
t318 = qJ(3) * t435;
t316 = qJ(3) * t437;
t306 = qJ(3) + t386;
t300 = -qJ(3) * t427 + t334;
t271 = t353 * t387 + t310;
t265 = pkin(4) * t297 + qJ(5) * t295;
t262 = -pkin(4) * t350 + t292 * t349 - t309 * t352;
t261 = qJ(5) * t350 + t431;
t255 = -pkin(4) * t427 + t280 * t349 - t303 * t352;
t254 = qJ(5) * t427 + t433;
t247 = (-qJD(4) * t386 + qJD(5) * t349) * t353 + (-pkin(6) + t378) * t425;
t245 = t444 - t258;
t243 = -pkin(4) * t423 + qJD(4) * t431 + t272 * t349 - t304 * t352;
t242 = qJ(5) * t423 + qJD(5) * t350 + t370;
t1 = [(t294 * t350 + t321 * t423) * MDP(19) + (-t396 * t350 + t251 * t423 - t302 * t295 + t310 * t259 + ((-qJD(4) * t309 - t272) * t321 - t292 * t294 - t284 * t417) * t349 + ((-qJD(4) * t292 + t304) * t321 + t309 * t294 + t264 * t353 - t284 * t425) * t352 + t391) * MDP(20) + (-t370 * t321 - t431 * t294 - t302 * t297 - t310 * t258 + ((qJD(2) * t284 + t421) * t349 - t408) * t350 + (-qJD(2) * t252 - t264 * t349 - t284 * t418) * t353 - t392) * MDP(21) + (-t243 * t321 + t247 * t295 + t259 * t271 - t262 * t294 + (-t253 * t424 - t240) * t350 + (-qJD(2) * t244 + t241 * t352 - t253 * t419) * t353 + t391) * MDP(22) + (-t242 * t295 + t243 * t297 - t258 * t262 - t259 * t261 + t325 + t382 * t425 + (-t453 - t239 * t352 - t240 * t349 + (-t244 * t352 + t246 * t349) * qJD(4)) * t353) * MDP(23) + (t242 * t321 - t247 * t297 + t258 * t271 + t261 * t294 + (-t253 * t426 + t239) * t350 + (qJD(2) * t246 + t241 * t349 + t253 * t418) * t353 + t392) * MDP(24) + (t239 * t261 + t246 * t242 + t241 * t271 + t253 * t247 + t240 * t262 + t244 * t243 - g(1) * (pkin(3) * t354 + pkin(4) * t289 + qJ(5) * t288 + t341) - g(2) * (pkin(4) * t287 + pkin(7) * t435 + qJ(5) * t286 + t395) + (-g(1) * (t379 - t457) - g(2) * pkin(3)) * t351) * MDP(25) + qJDD(1) * MDP(1) + (qJDD(1) * t347 + 0.2e1 * t350 * t400) * MDP(4) + 0.2e1 * (t350 * t411 - t413 * t429) * MDP(5) + (qJDD(2) * t350 + t353 * t356) * MDP(6) + (qJDD(2) * t353 - t350 * t356) * MDP(7) + (t350 * t375 + t353 * t368 + t325) * MDP(9) + (t375 * t353 + (-t368 - t456) * t350) * MDP(10) + ((t347 + t348) * qJDD(1) * pkin(6) + t360 - t390) * MDP(11) + (t350 * t366 + t353 * t361 - t325) * MDP(12) + (t366 * t353 + (-t361 + t456) * t350) * MDP(13) + (pkin(6) * t360 - g(1) * t341 - g(2) * t395 + t260 * t307 + t285 * t282 - t379 * t456) * MDP(14) + (t258 * t349 * t353 + (t349 * t425 - t352 * t417) * t297) * MDP(15) + ((-t295 * t349 + t297 * t352) * t425 + (t446 + t259 * t349 + (t295 * t352 + t297 * t349) * qJD(4)) * t353) * MDP(16) + ((t321 * t426 - t258) * t350 + (qJD(2) * t297 + t373) * t353) * MDP(17) + ((t321 * t424 - t259) * t350 + (-qJD(2) * t295 - t465) * t353) * MDP(18) + (-t453 + t456) * MDP(2) + t390 * MDP(3); -MDP(4) * t409 + t429 * MDP(5) * t357 + MDP(6) * t412 + MDP(7) * t411 + qJDD(2) * MDP(8) + (pkin(1) * t439 + t394) * MDP(9) + (t452 - t328 + (pkin(1) * t357 + t390) * t353) * MDP(10) + ((-pkin(2) * t350 + t450) * qJDD(1) + ((-t308 - t346) * t350 + (-t305 + t397) * t353) * qJD(1)) * MDP(11) + (-t300 * t427 + t374 - 0.2e1 * t448) * MDP(12) + (t328 + 0.2e1 * t344 + 0.2e1 * t345 + (qJD(1) * t300 - g(3)) * t350 + (qJD(1) * t285 - t390) * t353) * MDP(13) + (-t275 * qJ(3) - t308 * qJD(3) - t279 * pkin(2) - t285 * t300 - g(1) * (-pkin(2) * t440 + t318) - g(2) * (-pkin(2) * t441 + t316) - t402 - t381 * qJD(1) * pkin(6)) * MDP(14) + (-t349 * t443 - t446) * MDP(15) + ((-t259 - t443) * t352 + (t258 + t444) * t349) * MDP(16) + ((-t297 * t353 - t349 * t442) * qJD(1) + t465) * MDP(17) + ((t295 * t353 - t352 * t442) * qJD(1) + t373) * MDP(18) - t321 * MDP(19) * t427 + (-t251 * t427 + qJ(3) * t259 + t415 * t295 + (-t434 + (t284 - t303) * t321) * t352 + (t264 + (t280 + t416) * t321 + t364) * t349) * MDP(20) + (-qJ(3) * t258 + t433 * t321 + t252 * t427 + t415 * t297 + (-t284 * t321 + t434) * t349 + (t264 + t359) * t352) * MDP(21) + (t244 * t427 + t255 * t321 + t259 * t306 + t432 * t295 + t349 * t466 + t365 * t352) * MDP(22) + (t254 * t295 - t255 * t297 + (-t258 * t460 + (t295 * t460 - t246) * qJD(4) + t388) * t352 + (t259 * t460 + (-t297 * t460 - t244) * qJD(4) + t389) * t349 + t407) * MDP(23) + (-t246 * t427 - t254 * t321 + t258 * t306 - t432 * t297 + t365 * t349 - t352 * t466) * MDP(24) + (t241 * t306 - t246 * t254 - t244 * t255 - g(1) * t318 - g(2) * t316 - t402 + (-g(3) * pkin(7) - t386 * t390) * t353 + t432 * t253 + (-g(3) * t386 + t390 * t460) * t350 - (qJD(4) * t382 + t239 * t349 - t240 * t352) * t460) * MDP(25); MDP(11) * t412 + (qJDD(2) + t409) * MDP(12) + (-t347 * t357 - t356) * MDP(13) + (t319 + t374 - t448) * MDP(14) + t283 * MDP(20) - t407 * MDP(25) + (t308 * MDP(14) - t253 * MDP(25) + t295 * t463 + t410 * t297) * qJD(2) + (t294 * MDP(22) + (-t295 * t428 + t258 - t420) * MDP(23) + (qJD(4) * t246 - t388) * MDP(25) + t410 * t467) * t352 + (t468 * MDP(23) + (qJD(4) * t244 - t389) * MDP(25) + t410 * t294 + t463 * t467) * t349; MDP(15) * t445 + (-t295 ^ 2 + t461) * MDP(16) + t245 * MDP(17) + t468 * MDP(18) + t294 * MDP(19) + (-t284 * t297 + t363 + t447) * MDP(20) + (t251 * t321 + t284 * t295 - t462) * MDP(21) + (-t265 * t295 - t358 + t447 + 0.2e1 * t458) * MDP(22) + (pkin(4) * t258 - qJ(5) * t259 + (t246 - t252) * t297 + (t244 - t414) * t295) * MDP(23) + (0.2e1 * t449 - t253 * t295 + t265 * t297 + (0.2e1 * qJD(5) - t251) * t321 + t462) * MDP(24) + (t239 * qJ(5) - t240 * pkin(4) - t253 * t265 - t244 * t252 - g(1) * (-pkin(4) * t286 + qJ(5) * t287) - g(2) * (pkin(4) * t288 - qJ(5) * t289) + t387 * t343 + t414 * t246) * MDP(25); (-t294 + t445) * MDP(22) + t245 * MDP(23) + (-t461 - t467) * MDP(24) + (-t246 * t321 + t358 - t458) * MDP(25);];
tau = t1;
