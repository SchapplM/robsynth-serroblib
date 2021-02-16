% Calculate vector of inverse dynamics joint torques for
% S5RRPRP8
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
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:28
% EndTime: 2021-01-15 20:52:37
% DurationCPUTime: 3.74s
% Computational Cost: add. (1881->357), mult. (3947->409), div. (0->0), fcn. (2452->6), ass. (0->173)
t379 = sin(qJ(4));
t380 = sin(qJ(2));
t382 = cos(qJ(4));
t383 = cos(qJ(2));
t405 = t380 * t379 + t383 * t382;
t297 = t405 * qJD(1);
t295 = t297 ^ 2;
t444 = qJD(1) * t383;
t445 = qJD(1) * t380;
t299 = -t379 * t444 + t382 * t445;
t369 = qJDD(2) - qJDD(4);
t473 = t299 ^ 2;
t440 = qJD(4) * t382;
t442 = qJD(2) * t383;
t395 = t379 * t442 + t380 * t440;
t439 = qJD(1) * qJD(2);
t431 = t380 * t439;
t441 = qJD(4) * t379;
t432 = t383 * t441;
t449 = qJD(1) * t432 + t382 * t431;
t254 = qJD(1) * t395 + qJDD(1) * t405 - t449;
t371 = qJD(2) - qJD(4);
t458 = t299 * t371;
t481 = t254 + t458;
t437 = qJDD(1) * t383;
t486 = t431 - t437;
t430 = t383 * t439;
t438 = qJDD(1) * t380;
t487 = t430 + t438;
t397 = -t486 * t379 - t487 * t382;
t398 = t405 * qJD(4);
t253 = qJD(1) * t398 + t397;
t459 = t297 * t371;
t489 = t253 + t459;
t490 = -t489 * MDP(17) + t299 * t297 * MDP(15) - (-t473 + t295) * MDP(16) - t481 * MDP(18) - t369 * MDP(19);
t381 = sin(qJ(1));
t384 = cos(qJ(1));
t479 = g(1) * t384 + g(2) * t381;
t355 = pkin(6) * t437;
t372 = qJDD(2) * qJ(3);
t373 = qJD(2) * qJD(3);
t282 = -pkin(6) * t431 + t355 + t372 + t373;
t354 = pkin(6) * t438;
t429 = pkin(6) * t430 + qJDD(3) + t354;
t463 = qJDD(2) * pkin(2);
t289 = t429 - t463;
t359 = pkin(6) * t445;
t424 = -qJD(2) * pkin(2) + qJD(3);
t324 = t359 + t424;
t360 = pkin(6) * t444;
t374 = qJD(2) * qJ(3);
t327 = t360 + t374;
t406 = t324 * t383 - t327 * t380;
t488 = qJD(2) * t406 + t282 * t383 + t289 * t380 - t479;
t455 = t383 * t379;
t457 = t380 * t382;
t313 = -t455 + t457;
t485 = g(3) * t405 - t479 * t313;
t319 = pkin(7) * t445 - t359;
t385 = -pkin(2) - pkin(3);
t433 = t385 * qJD(2);
t286 = qJD(3) + t433 - t319;
t321 = -pkin(7) * t444 + t360;
t301 = t321 + t374;
t407 = -t286 * t379 - t301 * t382;
t465 = qJ(5) * t297;
t250 = -t407 - t465;
t483 = t250 * t371;
t480 = -t254 * qJ(5) - t297 * qJD(5);
t472 = pkin(6) - pkin(7);
t329 = t472 * t380;
t330 = t472 * t383;
t450 = t379 * t329 + t382 * t330;
t423 = -qJ(3) * t379 + t382 * t385;
t367 = g(1) * t381;
t469 = g(2) * t384;
t426 = t367 - t469;
t476 = pkin(4) * t254 + qJDD(5);
t302 = -qJD(1) * pkin(1) - pkin(2) * t444 - qJ(3) * t445;
t281 = pkin(3) * t444 - t302;
t425 = -pkin(4) * t297 - qJD(5);
t262 = t281 - t425;
t475 = t262 * t297 - t480;
t269 = -t487 * pkin(7) + t385 * qJDD(2) + t429;
t270 = t486 * pkin(7) + t282;
t416 = t379 * t269 + t382 * t270 + t286 * t440 - t301 * t441;
t390 = g(3) * t313 + t405 * t479 - t416;
t415 = -t382 * t269 + t379 * t270 + t286 * t441 + t301 * t440;
t391 = -t415 + t485;
t410 = pkin(2) * t383 + qJ(3) * t380;
t325 = pkin(1) + t410;
t468 = pkin(6) * qJDD(2);
t474 = (qJD(1) * t325 - t302) * qJD(2) + t468;
t470 = pkin(4) * t369;
t466 = qJ(5) * t253;
t464 = qJ(5) * t299;
t377 = qJDD(1) * pkin(1);
t461 = t281 * t297;
t460 = t281 * t299;
t387 = qJD(1) ^ 2;
t456 = t380 * t387;
t422 = t382 * t286 - t301 * t379;
t249 = t422 - t464;
t248 = -pkin(4) * t371 + t249;
t454 = -t249 + t248;
t451 = t382 * t319 + t379 * t321;
t257 = t451 + t464;
t287 = qJD(3) * t382 + qJD(4) * t423;
t453 = t287 - t257;
t421 = -t319 * t379 + t382 * t321;
t256 = t421 - t465;
t323 = qJ(3) * t382 + t379 * t385;
t288 = -qJD(3) * t379 - t323 * qJD(4);
t452 = t288 - t256;
t363 = t380 * qJD(3);
t448 = qJ(3) * t442 + t363;
t375 = t380 ^ 2;
t376 = t383 ^ 2;
t446 = t375 - t376;
t443 = qJD(2) * t380;
t436 = t385 * t380;
t435 = t383 * t456;
t434 = -g(3) * t383 + t479 * t380;
t408 = pkin(2) * t437 + t487 * qJ(3) + qJD(1) * t363 + t377;
t404 = pkin(3) * t437 + t408;
t413 = t380 * t433;
t255 = qJD(1) * t413 + t404;
t247 = t255 + t476;
t428 = -t247 + t469;
t427 = -t255 + t469;
t420 = t382 * t329 - t330 * t379;
t344 = pkin(4) * t382 - t385;
t345 = pkin(4) * t379 + qJ(3);
t419 = t344 * t383 + t345 * t380;
t417 = t371 ^ 2;
t414 = -t354 + t434;
t386 = qJD(2) ^ 2;
t412 = pkin(6) * t386 + t469;
t409 = -pkin(2) * t380 + qJ(3) * t383;
t349 = qJ(3) * t444;
t290 = qJD(1) * t436 + t349;
t403 = -0.2e1 * pkin(1) * t439 - t468;
t309 = t383 * pkin(3) + t325;
t278 = t413 + t448;
t320 = t472 * t443;
t322 = qJD(2) * t330;
t399 = -t382 * t320 + t379 * t322 + t329 * t440 - t330 * t441;
t396 = -t412 + 0.2e1 * t377;
t394 = -qJD(4) * t450 + t320 * t379 + t382 * t322;
t267 = pkin(2) * t431 - t408;
t292 = pkin(2) * t443 - t448;
t393 = -qJD(1) * t292 + qJDD(1) * t325 - t267 - t412;
t389 = -t288 * t371 - t391;
t388 = t287 * t371 + t323 * t369 - t390;
t370 = qJ(5) - t472;
t342 = t383 * t367;
t318 = -pkin(4) + t423;
t317 = pkin(2) * t445 - t349;
t294 = t405 * t367;
t293 = t313 * t367;
t280 = pkin(1) + t419;
t273 = pkin(4) * t405 + t309;
t272 = qJD(2) * t405 - t398;
t271 = -t382 * t443 + t395 - t432;
t268 = -pkin(4) * t299 + t290;
t261 = -qJ(5) * t405 + t450;
t260 = -qJ(5) * t313 + t420;
t251 = pkin(4) * t271 + t278;
t246 = -qJ(5) * t272 - qJD(5) * t313 + t394;
t245 = -qJ(5) * t271 - qJD(5) * t405 + t399;
t244 = t416 + t480;
t243 = -qJD(5) * t299 - t415 + t466 - t470;
t1 = [qJDD(1) * MDP(1) + t426 * MDP(2) + t479 * MDP(3) + (qJDD(1) * t375 + 0.2e1 * t380 * t430) * MDP(4) + 0.2e1 * (t380 * t437 - t446 * t439) * MDP(5) + (qJDD(2) * t380 + t383 * t386) * MDP(6) + (qJDD(2) * t383 - t380 * t386) * MDP(7) + (t380 * t403 + t383 * t396 + t342) * MDP(9) + (t403 * t383 + (-t396 - t367) * t380) * MDP(10) + (-t380 * t474 + t393 * t383 + t342) * MDP(11) + ((t375 + t376) * qJDD(1) * pkin(6) + t488) * MDP(12) + (t474 * t383 + (t393 + t367) * t380) * MDP(13) + (t302 * t292 + (-t267 + t426) * t325 + t488 * pkin(6)) * MDP(14) + (-t253 * t313 + t272 * t299) * MDP(15) + (t253 * t405 - t254 * t313 - t271 * t299 - t272 * t297) * MDP(16) + (-t272 * t371 - t313 * t369) * MDP(17) + (t271 * t371 + t369 * t405) * MDP(18) + (t309 * t254 + t281 * t271 + t278 * t297 - t420 * t369 - t394 * t371 - t405 * t427 + t294) * MDP(20) + (-t309 * t253 + t281 * t272 + t278 * t299 - t427 * t313 + t450 * t369 + t399 * t371 + t293) * MDP(21) + (-t246 * t371 + t251 * t297 + t254 * t273 - t260 * t369 + t262 * t271 - t405 * t428 + t294) * MDP(22) + (t245 * t371 + t251 * t299 - t253 * t273 + t261 * t369 + t262 * t272 - t428 * t313 + t293) * MDP(23) + (-t243 * t313 - t244 * t405 - t245 * t297 - t246 * t299 - t248 * t272 - t250 * t271 + t253 * t260 - t254 * t261 + t479) * MDP(24) + (t244 * t261 + t250 * t245 + t243 * t260 + t248 * t246 + t247 * t273 + t262 * t251 - g(1) * (-t280 * t381 - t370 * t384) - g(2) * (t280 * t384 - t370 * t381)) * MDP(25); (-t257 * t371 - t268 * t299 + t388 - t475) * MDP(23) + (pkin(1) * t456 + t414) * MDP(9) + (t409 * qJDD(1) + ((t327 - t374) * t380 + (-t324 + t424) * t383) * qJD(1)) * MDP(12) + (g(3) * t380 - t355 + (pkin(1) * t387 + t479) * t383) * MDP(10) + (t355 + 0.2e1 * t372 + 0.2e1 * t373 + (qJD(1) * t317 - g(3)) * t380 + (qJD(1) * t302 - t479) * t383) * MDP(13) + (t244 * t323 + t243 * t318 - t262 * t268 - g(3) * t419 - t479 * (-t344 * t380 + t345 * t383) + t453 * t250 + t452 * t248) * MDP(25) + (-t406 * qJD(1) * pkin(6) - t289 * pkin(2) - g(3) * t410 + t282 * qJ(3) + t327 * qJD(3) - t302 * t317 - t409 * t479) * MDP(14) + (-t290 * t299 - t451 * t371 + t388 - t461) * MDP(21) + (t253 * t318 - t254 * t323 + (-t250 - t452) * t299 + (t248 - t453) * t297) * MDP(24) + (-t290 * t297 - t423 * t369 + t421 * t371 + t389 + t460) * MDP(20) + (-t466 + t256 * t371 - t268 * t297 + (pkin(4) - t318) * t369 + (qJD(5) + t262) * t299 + t389) * MDP(22) + t446 * MDP(5) * t387 - MDP(4) * t435 + (0.2e1 * t463 - qJDD(3) + (-t302 * t380 + t317 * t383) * qJD(1) + t414) * MDP(11) + qJDD(2) * MDP(8) + MDP(7) * t437 + MDP(6) * t438 - t490; (-qJDD(2) - t435) * MDP(11) + MDP(12) * t438 + (-t375 * t387 - t386) * MDP(13) + (-qJD(2) * t327 + t302 * t445 + t289 - t434) * MDP(14) + (-t262 * t445 - t434) * MDP(25) + (MDP(21) + MDP(23)) * (-t299 * t445 + t379 * t369 - t382 * t417) + (MDP(20) + MDP(22)) * (-t297 * t445 - t369 * t382 - t379 * t417) + (t489 * MDP(24) + (t243 - t483) * MDP(25)) * t382 + (-t481 * MDP(24) + (t248 * t371 + t244) * MDP(25)) * t379; (t371 * t407 + t391 - t460) * MDP(20) + (-t422 * t371 + t390 + t461) * MDP(21) + (-0.2e1 * t470 + t466 - t483 + (-t262 + t425) * t299 + t391) * MDP(22) + (-pkin(4) * t473 - t249 * t371 + t390 + t475) * MDP(23) + (pkin(4) * t253 - t454 * t297) * MDP(24) + (t454 * t250 + (-t262 * t299 + t243 + t485) * pkin(4)) * MDP(25) + t490; (t379 * t438 + t382 * t437 - t449 - t458) * MDP(22) + (-t397 + t459) * MDP(23) + (-t295 - t473) * MDP(24) + (t248 * t299 + t250 * t297 + t404 + t426 + t476) * MDP(25) + ((MDP(22) * t457 - MDP(23) * t405) * qJD(4) + (MDP(22) * t455 + MDP(25) * t436) * qJD(2)) * qJD(1);];
tau = t1;
