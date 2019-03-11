% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:47
% EndTime: 2019-03-09 03:09:52
% DurationCPUTime: 3.52s
% Computational Cost: add. (3907->396), mult. (9141->535), div. (0->0), fcn. (6115->8), ass. (0->160)
t399 = sin(pkin(10));
t401 = cos(pkin(10));
t404 = sin(qJ(3));
t460 = qJD(1) * t404;
t365 = qJD(3) * t401 - t399 * t460;
t366 = qJD(3) * t399 + t401 * t460;
t403 = sin(qJ(5));
t405 = cos(qJ(5));
t317 = -t405 * t365 + t366 * t403;
t406 = cos(qJ(3));
t459 = qJD(1) * t406;
t389 = -qJD(5) + t459;
t490 = t317 * t389;
t474 = t405 * t401;
t368 = t399 * t403 - t474;
t452 = qJD(5) * t405;
t453 = qJD(5) * t403;
t483 = -t399 * t453 + t401 * t452;
t463 = t368 * t459 + t483;
t369 = t399 * t405 + t401 * t403;
t359 = t369 * qJD(5);
t416 = t369 * t406;
t462 = -qJD(1) * t416 + t359;
t489 = MDP(5) * t404;
t397 = t404 ^ 2;
t488 = MDP(6) * (-t406 ^ 2 + t397);
t454 = qJD(3) * t406;
t433 = t454 * t474;
t441 = t399 * t454;
t305 = t359 * t404 + t403 * t441 - t433;
t349 = t368 * t404;
t448 = qJD(1) * qJD(3);
t438 = t404 * t448;
t467 = t305 * t389 - t349 * t438;
t430 = pkin(3) * t404 - qJ(4) * t406;
t373 = t430 * qJD(1);
t392 = sin(pkin(9)) * pkin(1) + pkin(7);
t381 = t392 * qJD(1);
t484 = qJD(2) * t406 - t404 * t381;
t310 = t401 * t373 - t399 * t484;
t476 = t401 * t406;
t419 = pkin(4) * t404 - pkin(8) * t476;
t287 = qJD(1) * t419 + t310;
t311 = t399 * t373 + t401 * t484;
t443 = t399 * t459;
t297 = -pkin(8) * t443 + t311;
t481 = pkin(8) + qJ(4);
t378 = t481 * t399;
t379 = t481 * t401;
t422 = -t378 * t405 - t379 * t403;
t487 = -qJD(4) * t368 + qJD(5) * t422 - t403 * t287 - t405 * t297;
t327 = -t378 * t403 + t379 * t405;
t486 = qJD(4) * t369 + qJD(5) * t327 + t287 * t405 - t297 * t403;
t393 = -cos(pkin(9)) * pkin(1) - pkin(2);
t362 = -pkin(3) * t406 - qJ(4) * t404 + t393;
t351 = t401 * t362;
t477 = t401 * t404;
t303 = -pkin(8) * t477 + t351 + (-t392 * t399 - pkin(4)) * t406;
t372 = t392 * t476;
t323 = t399 * t362 + t372;
t479 = t399 * t404;
t309 = -pkin(8) * t479 + t323;
t485 = t403 * t303 + t405 * t309;
t447 = MDP(21) + MDP(23);
t446 = MDP(22) - MDP(25);
t437 = t406 * t448;
t432 = t399 * t437;
t464 = qJD(1) * t433 + t365 * t452;
t281 = t403 * (qJD(5) * t366 + t432) - t464;
t356 = qJD(3) * t430 - qJD(4) * t404;
t455 = qJD(3) * t404;
t442 = t392 * t455;
t314 = t401 * t356 + t399 * t442;
t413 = t419 * qJD(3);
t295 = t413 + t314;
t345 = t399 * t356;
t383 = t404 * t392;
t478 = t399 * t406;
t302 = t345 + (-pkin(8) * t478 - t383 * t401) * qJD(3);
t482 = -qJD(5) * t485 + t295 * t405 - t302 * t403;
t396 = t404 * qJD(2);
t353 = t406 * t381 + t396;
t339 = qJD(3) * qJ(4) + t353;
t342 = t362 * qJD(1);
t291 = t401 * t339 + t399 * t342;
t279 = pkin(8) * t365 + t291;
t480 = t279 * t403;
t411 = qJD(3) * t416;
t306 = t404 * t483 + t411;
t294 = t306 * t389;
t407 = qJD(3) ^ 2;
t475 = t404 * t407;
t473 = t406 * t407;
t472 = qJ(6) * t460 - t487;
t471 = pkin(5) * t460 + t486;
t320 = t365 * t403 + t366 * t405;
t282 = qJD(1) * t411 + qJD(5) * t320;
t470 = t349 * t282 + t305 * t317;
t328 = pkin(4) * t443 + t353;
t469 = -pkin(5) * t462 + qJ(6) * t463 + qJD(6) * t369 + t328;
t348 = t369 * t404;
t466 = -t348 * t438 + t294;
t334 = (qJD(4) + t484) * qJD(3);
t343 = t356 * qJD(1);
t286 = t401 * t334 + t399 * t343;
t344 = qJD(3) * t396 + t381 * t454;
t347 = pkin(4) * t441 + t392 * t454;
t355 = pkin(4) * t479 + t383;
t382 = qJD(1) * t393;
t457 = qJD(3) * t422;
t456 = qJD(3) * t327;
t451 = qJD(6) * t389;
t450 = t404 * MDP(20);
t290 = -t339 * t399 + t401 * t342;
t274 = -pkin(4) * t459 - pkin(8) * t366 + t290;
t257 = t274 * t405 - t480;
t449 = qJD(6) - t257;
t445 = t392 * t478;
t285 = -t334 * t399 + t401 * t343;
t275 = qJD(1) * t413 + t285;
t280 = -pkin(8) * t432 + t286;
t444 = t274 * t452 + t403 * t275 + t405 * t280;
t324 = pkin(4) * t432 + t344;
t394 = -pkin(4) * t401 - pkin(3);
t308 = t317 * t455;
t436 = -qJD(3) * pkin(3) + qJD(4);
t435 = t447 * t403;
t434 = -t274 * t453 + t405 * t275 - t279 * t452 - t403 * t280;
t258 = t274 * t403 + t279 * t405;
t429 = -t281 * t348 + t306 * t320;
t428 = -t285 * t399 + t286 * t401;
t426 = -t290 * t399 + t291 * t401;
t424 = t303 * t405 - t309 * t403;
t420 = 0.2e1 * qJD(3) * t382;
t418 = -t257 * t389 - t444;
t415 = -t279 * t453 + t444;
t414 = t403 * t295 + t405 * t302 + t303 * t452 - t309 * t453;
t337 = -t484 + t436;
t251 = -pkin(5) * t438 - t434;
t410 = -t403 * t446 + t405 * t447;
t312 = -pkin(4) * t365 + t337;
t409 = -qJ(4) * t455 + (-t337 + t436) * t406;
t254 = pkin(5) * t282 + qJ(6) * t281 - qJD(6) * t320 + t324;
t322 = t351 - t445;
t315 = -t401 * t442 + t345;
t313 = pkin(5) * t368 - qJ(6) * t369 + t394;
t307 = t320 * t455;
t283 = pkin(5) * t348 + qJ(6) * t349 + t355;
t273 = pkin(5) * t320 + qJ(6) * t317;
t265 = pkin(5) * t406 - t424;
t264 = -qJ(6) * t406 + t485;
t263 = -t281 - t490;
t262 = pkin(5) * t317 - qJ(6) * t320 + t312;
t259 = pkin(5) * t306 + qJ(6) * t305 + qJD(6) * t349 + t347;
t256 = -qJ(6) * t389 + t258;
t255 = pkin(5) * t389 + t449;
t253 = -pkin(5) * t455 - t482;
t252 = qJ(6) * t455 - qJD(6) * t406 + t414;
t250 = qJ(6) * t438 + t415 - t451;
t1 = [0.2e1 * t437 * t489 - 0.2e1 * t448 * t488 + MDP(7) * t473 - MDP(8) * t475 + (-t392 * t473 + t404 * t420) * MDP(10) + (t392 * t475 + t406 * t420) * MDP(11) + (t344 * t479 + (-qJD(1) * t314 - t285) * t406 + ((t337 * t399 - t365 * t392) * t406 + (t290 + (t322 + t445) * qJD(1)) * t404) * qJD(3)) * MDP(12) + (t344 * t477 + (qJD(1) * t315 + t286) * t406 + ((t337 * t401 + t366 * t392) * t406 + (-t291 + (-t323 + t372) * qJD(1)) * t404) * qJD(3)) * MDP(13) + (-t314 * t366 + t315 * t365 + (-t285 * t401 - t286 * t399) * t404 + (-t290 * t401 - t291 * t399 + (-t322 * t401 - t323 * t399) * qJD(1)) * t454) * MDP(14) + (t285 * t322 + t286 * t323 + t290 * t314 + t291 * t315 + (t337 * t454 + t344 * t404) * t392) * MDP(15) + (t281 * t349 - t305 * t320) * MDP(16) + (-t429 + t470) * MDP(17) + (t281 * t406 + t307 + t467) * MDP(18) + (t282 * t406 - t308 + t466) * MDP(19) + (-t389 - t459) * qJD(3) * t450 + (-t482 * t389 - t434 * t406 + t347 * t317 + t355 * t282 + t324 * t348 + t312 * t306 + (qJD(1) * t424 + t257) * t455) * MDP(21) + (t414 * t389 + t415 * t406 + t347 * t320 - t355 * t281 - t324 * t349 - t312 * t305 + (-qJD(1) * t485 - t258) * t455) * MDP(22) + (t251 * t406 + t253 * t389 + t254 * t348 + t259 * t317 + t262 * t306 + t282 * t283 + (-qJD(1) * t265 - t255) * t455) * MDP(23) + (-t250 * t348 - t251 * t349 - t252 * t317 + t253 * t320 - t255 * t305 - t256 * t306 - t264 * t282 - t265 * t281) * MDP(24) + (-t250 * t406 - t252 * t389 + t254 * t349 - t259 * t320 + t262 * t305 + t281 * t283 + (qJD(1) * t264 + t256) * t455) * MDP(25) + (t250 * t264 + t251 * t265 + t252 * t256 + t253 * t255 + t254 * t283 + t259 * t262) * MDP(26); (t308 + t466) * MDP(21) + (t307 - t467) * MDP(22) + (t308 + t294) * MDP(23) + (t429 + t470) * MDP(24) + t467 * MDP(25) + (-t250 * t349 + t251 * t348 + t255 * t306 - t256 * t305) * MDP(26) + (-t407 * MDP(11) - t344 * MDP(15) - t254 * MDP(26) + t446 * t281 - t447 * t282) * t406 + (-t407 * MDP(10) + t428 * MDP(15)) * t404 + ((-t399 * MDP(12) - t401 * MDP(13)) * t397 * qJD(1) + ((t365 * t401 + t366 * t399) * MDP(14) + t426 * MDP(15)) * t406 + (-MDP(23) * qJD(1) * t348 - MDP(12) * t365 + MDP(13) * t366 + MDP(15) * t337 - MDP(25) * t320 + MDP(26) * t262) * t404) * qJD(3); (qJD(3) * t353 - t382 * t460 - t344) * MDP(10) - t382 * t459 * MDP(11) + (-t344 * t401 + t353 * t365) * MDP(12) + (t344 * t399 - t353 * t366) * MDP(13) + (t310 * t366 - t311 * t365 + (qJD(4) * t365 + t290 * t459 + t286) * t401 + (qJD(4) * t366 + t291 * t459 - t285) * t399) * MDP(14) + (-pkin(3) * t344 + qJ(4) * t428 + qJD(4) * t426 - t290 * t310 - t291 * t311 - t337 * t353) * MDP(15) + (-t281 * t369 + t320 * t463) * MDP(16) + (t281 * t368 - t282 * t369 - t317 * t463 - t320 * t462) * MDP(17) + (-t463 * t389 + (qJD(3) * t369 - t320) * t460) * MDP(18) + (t462 * t389 + (-qJD(3) * t368 + t317) * t460) * MDP(19) + (t394 * t282 - t328 * t317 + t324 * t368 + t486 * t389 + t462 * t312 + (-t257 + t457) * t460) * MDP(21) + (-t394 * t281 - t328 * t320 + t324 * t369 + t487 * t389 + t463 * t312 + (t258 - t456) * t460) * MDP(22) + (t254 * t368 + t282 * t313 + t471 * t389 - t469 * t317 + t462 * t262 + (t255 + t457) * t460) * MDP(23) + (-t250 * t368 + t251 * t369 + t255 * t463 - t256 * t462 + t281 * t422 - t282 * t327 + t317 * t472 + t320 * t471) * MDP(24) + (-t254 * t369 + t281 * t313 + t472 * t389 + t469 * t320 - t463 * t262 + (-t256 + t456) * t460) * MDP(25) + (t250 * t327 - t251 * t422 + t254 * t313 + t255 * t471 - t256 * t472 - t262 * t469) * MDP(26) + ((-t290 * t404 + t310 * t406 + t399 * t409) * MDP(12) + (t291 * t404 - t311 * t406 + t401 * t409) * MDP(13) + t389 * t450 + (-t406 * t489 + t488) * qJD(1)) * qJD(1); (-t365 ^ 2 - t366 ^ 2) * MDP(14) + (t290 * t366 - t291 * t365 + t344) * MDP(15) - t317 ^ 2 * MDP(24) + (t256 * t317 + t254) * MDP(26) + (-t320 * MDP(24) - t255 * MDP(26) - t389 * t447) * t320 + (t365 * t435 + t366 * t410) * qJD(5) + (-t366 * MDP(12) - t365 * MDP(13) + ((MDP(13) + t435) * t401 + (MDP(12) + t410) * t399) * qJD(3)) * t459 - t446 * (-t464 - t490); t263 * MDP(18) + t418 * MDP(22) + (pkin(5) * t281 - qJ(6) * t282) * MDP(24) + (-t418 - 0.2e1 * t451) * MDP(25) + (-pkin(5) * t251 + qJ(6) * t250 - t255 * t258 + t256 * t449 - t262 * t273) * MDP(26) + (-MDP(19) * t320 + t446 * t480) * qJD(5) + (-t389 * MDP(19) - t312 * MDP(21) - t262 * MDP(23) + (t256 - t258) * MDP(24) + t273 * MDP(25) + MDP(17) * t320) * t320 + (-MDP(19) * t416 + (0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + MDP(20)) * t404) * t448 + (t320 * MDP(16) + t312 * MDP(22) - t273 * MDP(23) + (t255 - t449) * MDP(24) - t262 * MDP(25) - MDP(17) * t317) * t317 + t447 * (-t258 * t389 + t434); (t317 * t320 - t438) * MDP(23) + t263 * MDP(24) + (-t320 ^ 2 - t389 ^ 2) * MDP(25) + (t256 * t389 + t262 * t320 + t251) * MDP(26);];
tauc  = t1;
