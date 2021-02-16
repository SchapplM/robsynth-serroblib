% Calculate Coriolis joint torque vector for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:17
% EndTime: 2021-01-15 22:37:29
% DurationCPUTime: 3.98s
% Computational Cost: add. (3631->409), mult. (9064->543), div. (0->0), fcn. (5798->6), ass. (0->170)
t388 = sin(qJ(3));
t389 = sin(qJ(2));
t439 = qJD(1) * t389;
t423 = t388 * t439;
t390 = cos(qJ(3));
t428 = t390 * qJD(2);
t344 = t423 - t428;
t435 = qJD(2) * t388;
t346 = t390 * t439 + t435;
t386 = sin(pkin(8));
t387 = cos(pkin(8));
t300 = t387 * t344 + t346 * t386;
t391 = cos(qJ(2));
t438 = qJD(1) * t391;
t370 = -qJD(3) + t438;
t480 = t300 * t370;
t340 = t386 * t390 + t387 * t388;
t329 = t340 * qJD(3);
t448 = t340 * t438 - t329;
t432 = qJD(3) * t388;
t417 = t386 * t432;
t422 = t388 * t438;
t431 = qJD(3) * t390;
t463 = t387 * t390;
t447 = -t386 * t422 - t387 * t431 + t438 * t463 + t417;
t418 = t389 * t431;
t433 = qJD(2) * t391;
t421 = t388 * t433;
t474 = t418 + t421;
t401 = -t344 * t386 + t387 * t346;
t479 = t401 ^ 2;
t426 = qJD(1) * qJD(2);
t478 = -0.2e1 * t426;
t477 = MDP(4) * t389;
t358 = -qJD(2) * pkin(2) + pkin(6) * t439;
t313 = pkin(3) * t344 + qJD(4) + t358;
t262 = pkin(4) * t300 - qJ(5) * t401 + t313;
t476 = t262 * t401;
t384 = t389 ^ 2;
t475 = (-t391 ^ 2 + t384) * MDP(5);
t457 = t390 * t391;
t398 = pkin(3) * t389 - qJ(4) * t457;
t406 = pkin(2) * t389 - pkin(7) * t391;
t347 = t406 * qJD(1);
t442 = pkin(6) * t423 + t390 * t347;
t294 = t398 * qJD(1) + t442;
t331 = t388 * t347;
t459 = t389 * t390;
t460 = t388 * t391;
t304 = t331 + (-pkin(6) * t459 - qJ(4) * t460) * qJD(1);
t472 = qJ(4) + pkin(7);
t412 = qJD(3) * t472;
t430 = qJD(4) * t390;
t326 = -t388 * t412 + t430;
t396 = -qJD(4) * t388 - t390 * t412;
t451 = (t294 - t396) * t387 + (-t304 + t326) * t386;
t380 = pkin(6) * t438;
t407 = -t380 + (-t422 + t432) * pkin(3);
t473 = pkin(6) * t388;
t359 = qJD(2) * pkin(7) + t380;
t355 = -pkin(2) * t391 - pkin(7) * t389 - pkin(1);
t335 = t355 * qJD(1);
t462 = t388 * t335;
t308 = t359 * t390 + t462;
t288 = -qJ(4) * t344 + t308;
t285 = t387 * t288;
t307 = t390 * t335 - t359 * t388;
t287 = -qJ(4) * t346 + t307;
t260 = t287 * t386 + t285;
t471 = t260 * t401;
t470 = t288 * t386;
t419 = t389 * t432;
t414 = t391 * t426;
t425 = qJD(2) * qJD(3);
t443 = (t414 + t425) * t390;
t314 = -qJD(1) * t419 + t443;
t469 = t314 * t388;
t468 = t344 * t370;
t467 = t346 * t370;
t466 = t358 * t388;
t465 = t358 * t390;
t464 = t370 * t390;
t461 = t388 * t389;
t392 = qJD(2) ^ 2;
t458 = t389 * t392;
t456 = t391 * t392;
t393 = qJD(1) ^ 2;
t455 = t391 * t393;
t348 = t406 * qJD(2);
t336 = qJD(1) * t348;
t415 = t389 * t426;
t408 = pkin(6) * t415;
t446 = -t390 * t336 - t388 * t408;
t394 = -t308 * qJD(3) - t446;
t254 = pkin(3) * t415 - qJ(4) * t314 - qJD(4) * t346 + t394;
t413 = t388 * t425;
t315 = t474 * qJD(1) + t413;
t399 = -t335 * t431 - t388 * t336 + t359 * t432;
t259 = -qJ(4) * t315 - qJD(4) * t344 - t390 * t408 - t399;
t454 = -t387 * t254 + t386 * t259;
t243 = t386 * t254 + t387 * t259;
t269 = t386 * t294 + t387 * t304;
t264 = qJ(5) * t439 + t269;
t292 = t387 * t326 + t386 * t396;
t453 = t264 - t292;
t452 = pkin(4) * t439 + t451;
t450 = -t269 + t292;
t373 = pkin(6) * t457;
t434 = qJD(2) * t389;
t444 = t390 * t348 + t434 * t473;
t270 = -t389 * t430 + t398 * qJD(2) + (-t373 + (qJ(4) * t389 - t355) * t388) * qJD(3) + t444;
t445 = t388 * t348 + t355 * t431;
t275 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t459 + (-qJD(4) * t389 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t391) * t388 + t445;
t248 = t386 * t270 + t387 * t275;
t449 = t448 * pkin(4) - t447 * qJ(5) + qJD(5) * t340 - t407;
t283 = -pkin(3) * t370 + t287;
t256 = t386 * t283 + t285;
t342 = t390 * t355;
t305 = -qJ(4) * t459 + t342 + (-pkin(3) - t473) * t391;
t441 = t388 * t355 + t373;
t309 = -qJ(4) * t461 + t441;
t277 = t386 * t305 + t387 * t309;
t349 = pkin(3) * t461 + t389 * pkin(6);
t357 = t472 * t390;
t416 = t472 * t388;
t311 = t357 * t386 + t387 * t416;
t437 = qJD(2) * t311;
t312 = t387 * t357 - t386 * t416;
t436 = qJD(2) * t312;
t429 = qJD(5) * t370;
t261 = t287 * t387 - t470;
t427 = qJD(5) - t261;
t316 = t474 * pkin(3) + pkin(6) * t433;
t378 = -pkin(3) * t390 - pkin(2);
t420 = t391 * t428;
t298 = pkin(3) * t315 + pkin(6) * t414;
t411 = pkin(1) * t478;
t281 = t314 * t386 + t387 * t315;
t410 = t344 + t428;
t409 = -t346 + t435;
t405 = -t261 * t370 - t243;
t404 = -t260 * t370 - t454;
t282 = t314 * t387 - t315 * t386;
t403 = -t312 * t281 + t282 * t311 - t292 * t300;
t247 = t270 * t387 - t275 * t386;
t255 = t283 * t387 - t470;
t276 = t305 * t387 - t309 * t386;
t400 = qJD(1) * t384 - t370 * t391;
t241 = -pkin(4) * t415 + t454;
t244 = pkin(4) * t281 - qJ(5) * t282 - qJD(5) * t401 + t298;
t376 = -pkin(3) * t387 - pkin(4);
t374 = pkin(3) * t386 + qJ(5);
t368 = qJ(5) * t415;
t339 = t386 * t388 - t463;
t322 = -t386 * t461 + t387 * t459;
t321 = t340 * t389;
t296 = pkin(4) * t339 - qJ(5) * t340 + t378;
t290 = t389 * t329 + t386 * t421 - t387 * t420;
t289 = -t386 * t420 - t474 * t387 + t389 * t417;
t284 = pkin(4) * t321 - qJ(5) * t322 + t349;
t274 = pkin(4) * t391 - t276;
t273 = -qJ(5) * t391 + t277;
t263 = pkin(3) * t346 + pkin(4) * t401 + qJ(5) * t300;
t251 = -qJ(5) * t370 + t256;
t250 = pkin(4) * t370 + qJD(5) - t255;
t249 = -pkin(4) * t289 + qJ(5) * t290 - qJD(5) * t322 + t316;
t246 = -pkin(4) * t434 - t247;
t245 = qJ(5) * t434 - qJD(5) * t391 + t248;
t240 = t368 - t429 + t243;
t1 = [0.2e1 * t414 * t477 + t475 * t478 + MDP(6) * t456 - MDP(7) * t458 + (-pkin(6) * t456 + t389 * t411) * MDP(9) + (pkin(6) * t458 + t391 * t411) * MDP(10) + (t314 * t459 + (-t419 + t420) * t346) * MDP(11) + ((-t344 * t390 - t346 * t388) * t433 + (-t469 - t315 * t390 + (t344 * t388 - t346 * t390) * qJD(3)) * t389) * MDP(12) + (t370 * t419 - t314 * t391 + (t346 * t389 + t400 * t390) * qJD(2)) * MDP(13) + (t370 * t418 + t315 * t391 + (-t344 * t389 - t400 * t388) * qJD(2)) * MDP(14) + (-t370 - t438) * MDP(15) * t434 + (-(-t355 * t432 + t444) * t370 + (t358 * t431 + pkin(6) * t315 + (t342 * qJD(1) + t307) * qJD(2)) * t389 + ((pkin(6) * t344 + t466) * qJD(2) + (t462 + (pkin(6) * t370 + t359) * t390) * qJD(3) + t446) * t391) * MDP(16) + ((-pkin(6) * t391 * t432 + t445) * t370 - t399 * t391 + (pkin(6) * t314 - t358 * t432) * t389 + ((pkin(6) * t346 + t465) * t391 + (-pkin(6) * t464 - t441 * qJD(1) - t308) * t389) * qJD(2)) * MDP(17) + (t454 * t391 - t247 * t370 + t281 * t349 - t289 * t313 + t298 * t321 + t300 * t316 + (qJD(1) * t276 + t255) * t434) * MDP(18) + (t243 * t391 + t248 * t370 + t282 * t349 - t290 * t313 + t298 * t322 + t401 * t316 + (-qJD(1) * t277 - t256) * t434) * MDP(19) + (-t243 * t321 - t247 * t401 - t248 * t300 + t255 * t290 + t256 * t289 - t276 * t282 - t277 * t281 + t322 * t454) * MDP(20) + (t243 * t277 + t247 * t255 + t248 * t256 - t276 * t454 + t298 * t349 + t313 * t316) * MDP(21) + (t241 * t391 + t244 * t321 + t246 * t370 + t249 * t300 - t262 * t289 + t281 * t284 + (-qJD(1) * t274 - t250) * t434) * MDP(22) + (-t240 * t321 + t241 * t322 - t245 * t300 + t246 * t401 - t250 * t290 + t251 * t289 - t273 * t281 + t274 * t282) * MDP(23) + (-t240 * t391 - t244 * t322 - t245 * t370 - t249 * t401 + t262 * t290 - t282 * t284 + (qJD(1) * t273 + t251) * t434) * MDP(24) + (t240 * t273 + t241 * t274 + t244 * t284 + t245 * t251 + t246 * t250 + t249 * t262) * MDP(25); -t455 * t477 + t393 * t475 + (-t346 * t464 + t469) * MDP(11) + ((t314 + t468) * t390 + (-t315 + t467) * t388) * MDP(12) + (-t370 * t431 + (t370 * t457 + t409 * t389) * qJD(1)) * MDP(13) + (t370 * t432 + (-t370 * t460 + t410 * t389) * qJD(1)) * MDP(14) + t370 * MDP(15) * t439 + (-pkin(2) * t315 + t442 * t370 + (pkin(7) * t464 + t466) * qJD(3) + ((-pkin(7) * t435 - t307) * t389 + (-t410 * pkin(6) - t466) * t391) * qJD(1)) * MDP(16) + (-pkin(2) * t314 - t331 * t370 + (-pkin(7) * t370 * t388 + t465) * qJD(3) + (-t358 * t457 + (-pkin(7) * t428 + t308) * t389 + (t370 * t459 + t409 * t391) * pkin(6)) * qJD(1)) * MDP(17) + (t281 * t378 + t298 * t339 + t451 * t370 - t448 * t313 + t407 * t300 + (-t255 - t437) * t439) * MDP(18) + (t282 * t378 + t298 * t340 + t450 * t370 - t447 * t313 + t407 * t401 + (t256 - t436) * t439) * MDP(19) + (-t243 * t339 + t447 * t255 + t448 * t256 + t269 * t300 + t340 * t454 + t401 * t451 + t403) * MDP(20) + (t243 * t312 - t451 * t255 + t450 * t256 + t298 * t378 + t311 * t454 + t407 * t313) * MDP(21) + (t244 * t339 + t281 * t296 + t452 * t370 - t449 * t300 - t448 * t262 + (t250 - t437) * t439) * MDP(22) + (-t240 * t339 + t241 * t340 - t447 * t250 + t448 * t251 + t264 * t300 + t401 * t452 + t403) * MDP(23) + (-t244 * t340 - t282 * t296 + t453 * t370 + t449 * t401 + t447 * t262 + (-t251 + t436) * t439) * MDP(24) + (t240 * t312 + t241 * t311 + t244 * t296 + t452 * t250 - t453 * t251 - t449 * t262) * MDP(25) + (t393 * t389 * MDP(9) + MDP(10) * t455) * pkin(1); t346 * t344 * MDP(11) + (-t344 ^ 2 + t346 ^ 2) * MDP(12) + (t443 - t468) * MDP(13) + (-t413 - t467) * MDP(14) + (-t308 * t370 - t346 * t358 + t394) * MDP(16) + (-t307 * t370 + t344 * t358 + t399) * MDP(17) + (-t313 * t401 + t404) * MDP(18) + (t300 * t313 + t405) * MDP(19) + (t256 * t401 - t471 + (-t255 + t261) * t300) * MDP(20) + (t255 * t260 - t256 * t261) * MDP(21) + (-t263 * t300 + t404 - t476) * MDP(22) + (t251 * t401 - t281 * t374 + t282 * t376 - t471 + (t250 - t427) * t300) * MDP(23) + (-t262 * t300 + t263 * t401 + t368 - t405 - 0.2e1 * t429) * MDP(24) + (t240 * t374 + t241 * t376 - t250 * t260 + t427 * t251 - t262 * t263) * MDP(25) + (-MDP(14) * t421 + ((-t388 * MDP(13) - t390 * MDP(14)) * qJD(3) + (MDP(15) + t390 * pkin(6) * MDP(17) + (pkin(4) - t376) * MDP(22) + t374 * MDP(24)) * qJD(2)) * t389) * qJD(1) + ((-t300 * t346 + t387 * t415) * MDP(18) + (-t346 * t401 - t386 * t415) * MDP(19) + (-t281 * t386 - t282 * t387) * MDP(20) + (t243 * t386 - t313 * t346 - t387 * t454) * MDP(21)) * pkin(3); (t255 * t401 + t256 * t300 + t298) * MDP(21) + (-t250 * t401 + t251 * t300 + t244) * MDP(25) + (MDP(20) + MDP(23)) * (-t300 ^ 2 - t479) + (MDP(18) + MDP(22)) * (-t370 * t401 + t281) + (-MDP(19) + MDP(24)) * (-t282 - t480); (t300 * t401 - t415) * MDP(22) + (t282 - t480) * MDP(23) + (-t370 ^ 2 - t479) * MDP(24) + (t251 * t370 + t241 + t476) * MDP(25);];
tauc = t1;
