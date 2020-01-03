% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR16_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR16_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:34
% EndTime: 2019-12-31 20:47:45
% DurationCPUTime: 5.22s
% Computational Cost: add. (2696->418), mult. (7123->583), div. (0->0), fcn. (5186->8), ass. (0->182)
t368 = cos(pkin(5));
t447 = qJD(1) * t368;
t356 = qJD(2) + t447;
t370 = sin(qJ(4));
t373 = cos(qJ(4));
t374 = cos(qJ(2));
t367 = sin(pkin(5));
t448 = qJD(1) * t367;
t426 = t374 * t448;
t310 = t356 * t370 + t373 * t426;
t309 = qJD(5) + t310;
t371 = sin(qJ(2));
t446 = qJD(1) * t371;
t427 = t367 * t446;
t434 = pkin(1) * t447;
t450 = -pkin(7) * t427 + t374 * t434;
t476 = t450 - qJD(3);
t365 = t371 ^ 2;
t480 = MDP(5) * (-t374 ^ 2 + t365);
t343 = qJD(4) + t427;
t372 = cos(qJ(5));
t409 = t370 * t426;
t312 = t356 * t373 - t409;
t369 = sin(qJ(5));
t465 = t312 * t369;
t269 = -t372 * t343 + t465;
t479 = t269 * t343;
t460 = t367 * t371;
t358 = pkin(7) * t460;
t428 = -pkin(1) * t374 - pkin(2);
t284 = pkin(3) * t460 + t358 + (-pkin(8) + t428) * t368;
t375 = -pkin(2) - pkin(8);
t420 = -qJ(3) * t371 - pkin(1);
t298 = (t374 * t375 + t420) * t367;
t478 = t370 * t284 + t373 * t298;
t472 = pkin(1) * t371;
t361 = t368 * t472;
t459 = t367 * t374;
t477 = pkin(7) * t459 + t361;
t437 = pkin(3) * t427 - t476;
t435 = qJD(1) * qJD(2);
t421 = t367 * t435;
t345 = t374 * t421;
t264 = t356 * t375 + t437;
t281 = qJD(1) * t298;
t252 = t264 * t370 + t281 * t373;
t408 = t371 * t421;
t338 = pkin(2) * t408;
t402 = pkin(8) * t371 - qJ(3) * t374;
t443 = qJD(3) * t371;
t377 = (qJD(2) * t402 - t443) * t367;
t266 = qJD(1) * t377 + t338;
t433 = pkin(1) * qJD(2) * t368;
t411 = qJD(1) * t433;
t313 = pkin(7) * t345 + t371 * t411;
t285 = pkin(3) * t345 + t313;
t378 = -qJD(4) * t252 - t370 * t266 + t373 * t285;
t241 = -pkin(4) * t345 - t378;
t475 = t309 * (pkin(4) * t312 + pkin(9) * t309) + t241;
t445 = qJD(2) * t371;
t425 = t367 * t445;
t350 = pkin(2) * t425;
t278 = t350 + t377;
t473 = pkin(3) + pkin(7);
t303 = (t459 * t473 + t361) * qJD(2);
t474 = -qJD(4) * t478 - t278 * t370 + t303 * t373;
t273 = -qJD(4) * t310 + t370 * t408;
t438 = qJD(5) * t372;
t429 = t372 * t273 + t343 * t438 + t369 * t345;
t439 = qJD(5) * t369;
t249 = -t312 * t439 + t429;
t471 = t249 * t369;
t470 = t269 * t309;
t271 = t312 * t372 + t343 * t369;
t469 = t271 * t309;
t468 = t309 * t375;
t467 = t310 * t343;
t466 = t312 * t343;
t464 = t343 * t370;
t463 = t343 * t373;
t462 = t343 * t375;
t364 = t367 ^ 2;
t461 = t364 * qJD(1) ^ 2;
t440 = qJD(4) * t373;
t274 = -qJD(4) * t409 + t356 * t440 - t373 * t408;
t458 = t369 * t274;
t457 = t371 * t372;
t456 = t372 * t274;
t455 = t373 * t249;
t405 = pkin(4) * t373 + pkin(9) * t370;
t454 = (-pkin(3) - t405) * t427 - qJD(4) * t405 + t476;
t355 = pkin(2) * t427;
t300 = t402 * t448 + t355;
t319 = pkin(7) * t426 + t371 * t434;
t302 = pkin(3) * t426 + t319;
t452 = t373 * t300 + t370 * t302;
t451 = -pkin(7) * t408 + t374 * t411;
t444 = qJD(2) * t374;
t442 = qJD(4) * t370;
t441 = qJD(4) * t372;
t432 = t371 * t463;
t431 = t374 * t461;
t430 = t370 * t459;
t289 = -t356 * qJD(3) - t451;
t314 = -t368 * qJ(3) - t477;
t424 = t367 * t444;
t423 = t343 * t440;
t422 = t364 * t435;
t382 = t264 * t440 + t373 * t266 - t281 * t442 + t370 * t285;
t240 = pkin(9) * t345 + t382;
t267 = -pkin(3) * t408 - t289;
t245 = pkin(4) * t274 - pkin(9) * t273 + t267;
t418 = -t240 * t369 + t372 * t245;
t416 = t273 * t369 - t372 * t345;
t415 = t309 * t372;
t332 = pkin(4) * t370 - pkin(9) * t373 + qJ(3);
t413 = pkin(9) * t426 - qJD(5) * t332 + t452;
t410 = t371 * t431;
t346 = t356 * qJ(3);
t276 = t346 + t302;
t297 = pkin(3) * t459 - t314;
t407 = MDP(19) * t426;
t406 = -0.2e1 * pkin(1) * t422;
t404 = t319 * t356 - t313;
t307 = (t369 * t374 + t370 * t457) * t448;
t403 = -t370 * t441 - t307;
t401 = t240 * t372 + t245 * t369;
t247 = pkin(9) * t343 + t252;
t253 = pkin(4) * t310 - pkin(9) * t312 + t276;
t239 = t247 * t372 + t253 * t369;
t400 = t247 * t369 - t253 * t372;
t257 = pkin(9) * t460 + t478;
t321 = t368 * t370 + t373 * t459;
t322 = t368 * t373 - t430;
t260 = pkin(4) * t321 - pkin(9) * t322 + t297;
t399 = t257 * t372 + t260 * t369;
t398 = -t257 * t369 + t260 * t372;
t251 = t264 * t373 - t281 * t370;
t396 = t284 * t373 - t298 * t370;
t394 = -t300 * t370 + t302 * t373;
t320 = t477 * qJD(2);
t393 = t313 * t368 + t320 * t356;
t354 = t374 * t433;
t392 = -pkin(7) * t425 + t354;
t391 = -t356 * t426 + t345;
t390 = -t322 * t369 + t367 * t457;
t295 = t322 * t372 + t369 * t460;
t389 = t276 * t371 + t375 * t444;
t388 = -t309 * t438 - t458;
t387 = t309 * t439 - t456;
t384 = t343 * t271;
t381 = t373 * t278 + t284 * t440 - t298 * t442 + t370 * t303;
t315 = (-pkin(2) * t374 + t420) * t367;
t362 = t368 * qJD(3);
t283 = -t425 * t473 + t354 + t362;
t380 = (-qJ(3) * t444 - t443) * t367;
t246 = -pkin(4) * t343 - t251;
t379 = -pkin(9) * t274 + (t246 + t251) * t309;
t330 = t373 * t345;
t317 = -qJ(3) * t426 + t355;
t316 = t368 * t428 + t358;
t308 = -t362 - t392;
t306 = t369 * t370 * t427 - t372 * t426;
t305 = qJD(1) * t315;
t304 = t350 + t380;
t299 = -t346 - t319;
t296 = -pkin(2) * t356 - t476;
t293 = -qJD(4) * t430 + t368 * t440 - t373 * t425;
t292 = -qJD(4) * t321 + t370 * t425;
t287 = qJD(1) * t380 + t338;
t282 = t305 * t427;
t258 = -pkin(4) * t426 - t394;
t256 = -pkin(4) * t460 - t396;
t255 = qJD(5) * t390 + t292 * t372 + t369 * t424;
t254 = qJD(5) * t295 + t292 * t369 - t372 * t424;
t250 = qJD(5) * t271 + t416;
t248 = pkin(4) * t293 - pkin(9) * t292 + t283;
t243 = -pkin(4) * t424 - t474;
t242 = pkin(9) * t424 + t381;
t237 = -qJD(5) * t239 + t418;
t236 = -qJD(5) * t400 + t401;
t1 = [(t371 * t406 - t393) * MDP(9) + (-t356 * t392 - t368 * t451 + t374 * t406) * MDP(10) + (-t289 * t374 + t313 * t371 + (t296 * t374 + t299 * t371) * qJD(2) + (-t308 * t374 + t320 * t371 + (t314 * t371 + t316 * t374) * qJD(2)) * qJD(1)) * t367 * MDP(11) + ((-t305 * t445 + t287 * t374 + (t304 * t374 - t315 * t445) * qJD(1)) * t367 + t393) * MDP(12) + (-t289 * t368 - t308 * t356 + (-t305 * t444 - t287 * t371 + (-t304 * t371 - t315 * t444) * qJD(1)) * t367) * MDP(13) + (t287 * t315 + t289 * t314 + t296 * t320 + t299 * t308 + t304 * t305 + t313 * t316) * MDP(14) + (t273 * t322 + t292 * t312) * MDP(15) + (-t273 * t321 - t274 * t322 - t292 * t310 - t293 * t312) * MDP(16) + (t292 * t343 + (t273 * t371 + (qJD(1) * t322 + t312) * t444) * t367) * MDP(17) + (-t293 * t343 + (-t274 * t371 + (-qJD(1) * t321 - t310) * t444) * t367) * MDP(18) + (t343 * t367 + t364 * t446) * MDP(19) * t444 + (t474 * t343 + t283 * t310 + t297 * t274 + t267 * t321 + t276 * t293 + (t378 * t371 + (qJD(1) * t396 + t251) * t444) * t367) * MDP(20) + (-t381 * t343 + t283 * t312 + t297 * t273 + t267 * t322 + t276 * t292 + (-t382 * t371 + (-qJD(1) * t478 - t252) * t444) * t367) * MDP(21) + (t249 * t295 + t255 * t271) * MDP(22) + (t249 * t390 - t250 * t295 - t254 * t271 - t255 * t269) * MDP(23) + (t249 * t321 + t255 * t309 + t271 * t293 + t274 * t295) * MDP(24) + (-t250 * t321 - t254 * t309 - t269 * t293 + t274 * t390) * MDP(25) + (t274 * t321 + t293 * t309) * MDP(26) + ((-qJD(5) * t399 - t242 * t369 + t248 * t372) * t309 + t398 * t274 + t237 * t321 - t400 * t293 + t243 * t269 + t256 * t250 - t241 * t390 + t246 * t254) * MDP(27) + (-(qJD(5) * t398 + t242 * t372 + t248 * t369) * t309 - t399 * t274 - t236 * t321 - t239 * t293 + t243 * t271 + t256 * t249 + t241 * t295 + t246 * t255) * MDP(28) + (MDP(6) * t424 - MDP(7) * t425) * (t356 + t447) + 0.2e1 * (t371 * t374 * MDP(4) - t480) * t422; -MDP(4) * t410 - t343 * t407 + t461 * t480 + t391 * MDP(6) + (-qJD(2) + t356) * MDP(7) * t427 + (t461 * t472 + t404) * MDP(9) + (pkin(1) * t431 + t356 * t450 - t451) * MDP(10) + ((-qJ(3) * qJD(2) - t299 - t319) * t371 + (-pkin(2) * qJD(2) - t296 - t476) * t374) * MDP(11) * t448 + (-t317 * t426 + t282 - t404) * MDP(12) + (-t476 * t356 + (t305 * t374 + t317 * t371) * t448 - t289) * MDP(13) + (-pkin(2) * t313 - qJ(3) * t289 - t296 * t319 + t299 * t476 - t305 * t317) * MDP(14) + (t273 * t373 - t312 * t464) * MDP(15) + ((-t274 - t466) * t373 + (-t273 + t467) * t370) * MDP(16) + (-t343 * t442 + t330 + (-t312 * t374 - t371 * t464) * t448) * MDP(17) + (-t423 + (-t432 + (-qJD(2) * t370 + t310) * t374) * t448) * MDP(18) + (qJ(3) * t274 + t267 * t370 - t394 * t343 + t437 * t310 + (t276 * t373 - t370 * t462) * qJD(4) + (-t251 * t374 + t373 * t389) * t448) * MDP(20) + (qJ(3) * t273 + t267 * t373 + t452 * t343 + t437 * t312 + (-t276 * t370 - t373 * t462) * qJD(4) + (t252 * t374 - t370 * t389) * t448) * MDP(21) + (t372 * t455 + (-t373 * t439 + t403) * t271) * MDP(22) + (t269 * t307 + t271 * t306 + (t269 * t372 + t271 * t369) * t442 + (-t471 - t250 * t372 + (t269 * t369 - t271 * t372) * qJD(5)) * t373) * MDP(23) + (t249 * t370 + t403 * t309 + (t384 - t387) * t373) * MDP(24) + (-t250 * t370 + (t369 * t442 + t306) * t309 + (t388 - t479) * t373) * MDP(25) + (t274 * t370 + t309 * t463) * MDP(26) + (t332 * t456 - t246 * t306 - t258 * t269 + (t369 * t413 - t372 * t454) * t309 + (-t246 * t369 * qJD(4) + t237 + (qJD(4) * t269 + t388) * t375) * t370 + (-t400 * t427 + t246 * t438 + t241 * t369 - t375 * t250 + (-t369 * t468 - t400) * qJD(4)) * t373) * MDP(27) + (-t332 * t458 - t246 * t307 - t258 * t271 + (t369 * t454 + t372 * t413) * t309 + (-t246 * t441 - t236 + (qJD(4) * t271 + t387) * t375) * t370 + (-t239 * t427 - t246 * t439 + t241 * t372 - t375 * t249 + (-t372 * t468 - t239) * qJD(4)) * t373) * MDP(28); t391 * MDP(11) + MDP(12) * t410 + (-t356 ^ 2 - t365 * t461) * MDP(13) + (t299 * t356 + t282 + t313) * MDP(14) + (-t310 * t356 - t343 * t464 + t330) * MDP(20) + (-t423 - t312 * t356 + (-t370 * t444 - t432) * t448) * MDP(21) + (-t373 * t250 + (-t372 * t356 - t369 * t463) * t309 + (t388 + t479) * t370) * MDP(27) + (-t455 + (t369 * t356 - t372 * t463) * t309 + (t384 + t387) * t370) * MDP(28); -t310 ^ 2 * MDP(16) + (t273 + t467) * MDP(17) + (-t274 + t466) * MDP(18) + qJD(2) * t407 + (t252 * t343 + t378) * MDP(20) + (t251 * t343 + t276 * t310 - t382) * MDP(21) + (t271 * t415 + t471) * MDP(22) + ((t249 - t470) * t372 + (-t250 - t469) * t369) * MDP(23) + (t309 * t415 + t458) * MDP(24) + (-t309 ^ 2 * t369 + t456) * MDP(25) + (-pkin(4) * t250 - t252 * t269 + t379 * t369 - t372 * t475) * MDP(27) + (-pkin(4) * t249 - t252 * t271 + t369 * t475 + t379 * t372) * MDP(28) + (MDP(15) * t310 + t312 * MDP(16) - t276 * MDP(20) - t271 * MDP(24) + t269 * MDP(25) - t309 * MDP(26) + MDP(27) * t400 + MDP(28) * t239) * t312; t271 * t269 * MDP(22) + (-t269 ^ 2 + t271 ^ 2) * MDP(23) + (t429 + t470) * MDP(24) + (-t416 + t469) * MDP(25) + t274 * MDP(26) + (t239 * t309 - t246 * t271 + t418) * MDP(27) + (t246 * t269 - t309 * t400 - t401) * MDP(28) + (-MDP(24) * t465 - MDP(25) * t271 - MDP(27) * t239 + MDP(28) * t400) * qJD(5);];
tauc = t1;
