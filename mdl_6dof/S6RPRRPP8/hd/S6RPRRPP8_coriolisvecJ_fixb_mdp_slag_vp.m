% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:56:00
% EndTime: 2019-03-09 04:56:07
% DurationCPUTime: 4.06s
% Computational Cost: add. (3012->444), mult. (6047->575), div. (0->0), fcn. (3312->4), ass. (0->180)
t355 = sin(qJ(3));
t354 = sin(qJ(4));
t357 = cos(qJ(3));
t426 = qJD(4) * t357;
t406 = t354 * t426;
t356 = cos(qJ(4));
t424 = t356 * qJD(3);
t364 = t355 * t424 + t406;
t324 = pkin(3) * t355 - pkin(8) * t357 + qJ(2);
t302 = t324 * qJD(1);
t358 = -pkin(1) - pkin(7);
t339 = qJD(1) * t358 + qJD(2);
t323 = t355 * t339;
t304 = qJD(3) * pkin(8) + t323;
t267 = -t356 * t302 + t304 * t354;
t433 = qJD(3) * t354;
t434 = qJD(1) * t357;
t316 = t356 * t434 + t433;
t373 = pkin(5) * t316 + t267;
t423 = qJD(5) + t373;
t474 = MDP(23) + MDP(26);
t416 = 0.2e1 * qJD(1);
t352 = t357 ^ 2;
t481 = MDP(8) * (t355 ^ 2 - t352);
t480 = qJ(2) * MDP(6) + MDP(5);
t436 = qJD(1) * t355;
t344 = qJD(4) + t436;
t331 = t344 * qJD(5);
t420 = qJD(1) * qJD(3);
t347 = t357 * t420;
t342 = qJ(5) * t347;
t432 = qJD(3) * t355;
t410 = t354 * t432;
t332 = qJD(1) * t410;
t429 = qJD(4) * t316;
t281 = -t332 + t429;
t384 = pkin(3) * t357 + pkin(8) * t355;
t313 = qJD(3) * t384 + qJD(2);
t294 = t313 * qJD(1);
t449 = t356 * t357;
t310 = t339 * t449;
t427 = qJD(4) * t356;
t428 = qJD(4) * t354;
t388 = qJD(3) * t310 + t354 * t294 + t302 * t427 - t304 * t428;
t366 = -pkin(5) * t281 + t388;
t244 = t331 + t342 + t366;
t471 = pkin(4) + qJ(6);
t251 = -t344 * t471 + t423;
t479 = t251 * t344 + t244;
t435 = qJD(1) * t356;
t467 = qJ(5) * t355;
t477 = -pkin(4) * t428 + qJD(5) * t354 + t435 * t467 + t323;
t476 = qJD(1) * t364;
t475 = MDP(19) - MDP(22);
t417 = MDP(21) + MDP(25);
t341 = t344 ^ 2;
t473 = pkin(5) + pkin(8);
t314 = t354 * t434 - t424;
t472 = pkin(5) * t314;
t469 = qJ(5) * t281;
t468 = qJ(5) * t314;
t466 = qJ(5) * t356;
t419 = qJD(3) * qJD(4);
t348 = t356 * t419;
t280 = -t348 + t476;
t368 = qJ(5) * t280 - qJD(5) * t316 + t339 * t432;
t248 = pkin(4) * t281 + t368;
t465 = t248 * t354;
t464 = t248 * t356;
t268 = t354 * t302 + t356 * t304;
t261 = -qJ(5) * t344 - t268;
t253 = qJD(6) - t261 - t472;
t462 = t253 * t344;
t461 = t261 * t344;
t460 = t280 * t354;
t459 = t314 * t316;
t458 = t314 * t344;
t457 = t314 * t354;
t456 = t316 * t344;
t455 = t339 * t357;
t454 = t344 * t356;
t453 = t354 * t355;
t452 = t354 * t357;
t451 = t355 * t356;
t450 = t355 * t358;
t360 = qJD(1) ^ 2;
t448 = t357 * t360;
t359 = qJD(3) ^ 2;
t447 = t358 * t359;
t330 = t473 * t356;
t365 = -pkin(5) * t451 - t357 * t471;
t319 = t384 * qJD(1);
t397 = t319 * t356 - t339 * t452;
t446 = qJD(1) * t365 - qJD(4) * t330 - t397;
t441 = t354 * t319 + t310;
t445 = (pkin(5) * t453 + qJ(5) * t357) * qJD(1) + t441 + t473 * t428;
t381 = qJ(6) * t354 - t466;
t403 = t471 * t355;
t444 = -qJD(1) * t354 * t403 - qJD(4) * t381 + qJD(6) * t356 + t477;
t443 = -pkin(4) * t354 * t436 + qJ(5) * t427 + t477;
t439 = t354 * t324 + t356 * t450;
t437 = -t359 - t360;
t431 = qJD(3) * t357;
t430 = qJD(3) * t358;
t425 = qJD(4) * t358;
t422 = qJD(5) + t267;
t257 = t268 - t472;
t421 = -qJD(6) - t257;
t418 = MDP(19) + MDP(27);
t415 = pkin(8) * t344 * t354;
t414 = pkin(8) * t454;
t413 = pkin(8) * t431;
t407 = t357 * t430;
t411 = t354 * t313 + t324 * t427 + t356 * t407;
t409 = t354 * t431;
t405 = t355 * t425;
t404 = t356 * t426;
t402 = -qJ(5) * t354 - pkin(3);
t401 = qJD(3) * t471;
t400 = MDP(22) - t418;
t399 = MDP(20) - t474;
t305 = -qJD(3) * pkin(3) - t455;
t398 = -t305 + t455;
t395 = t348 + t458;
t336 = t354 * t450;
t394 = t324 * t356 - t336;
t393 = t314 + t424;
t392 = -t316 + t433;
t390 = t356 * t417;
t387 = t356 * t294 - t302 * t428 - t304 * t427 - t339 * t409;
t386 = pkin(4) * t404 + t364 * qJ(5) + t355 * t430;
t385 = t355 * t347;
t278 = -t439 - t467;
t383 = t400 * t354;
t382 = t331 + t388;
t245 = -t342 - t382;
t247 = -pkin(4) * t347 - t387;
t380 = -t245 * t356 + t247 * t354;
t379 = t251 * t356 - t253 * t354;
t378 = t251 * t354 + t253 * t356;
t260 = -pkin(4) * t344 + t422;
t377 = t260 * t356 + t261 * t354;
t376 = t260 * t354 - t261 * t356;
t375 = qJD(1) * t352 - t344 * t355;
t374 = -t324 * t428 - t354 * t407 + (t313 - t405) * t356;
t363 = -qJ(5) * t316 + t305;
t265 = pkin(4) * t314 + t363;
t372 = -t265 * t355 + t413;
t371 = t305 * t355 - t413;
t243 = qJD(6) * t314 + t281 * t471 + t368;
t254 = t314 * t471 + t363;
t370 = t243 * t354 + t254 * t427;
t369 = -t243 * t356 + t254 * t428;
t367 = pkin(5) * t280 + t387;
t362 = -qJD(6) * t344 - t367;
t361 = t377 * MDP(24) + t379 * MDP(28) + t417 * t457 + (t354 * t399 + t356 * t400) * t344;
t346 = pkin(4) * t452;
t340 = 0.2e1 * t342;
t329 = t473 * t354;
t325 = -pkin(4) * t356 + t402;
t322 = t356 * t385;
t308 = -t356 * t471 + t402;
t284 = t346 + (-t358 - t466) * t357;
t279 = -pkin(4) * t355 - t394;
t277 = t346 + (-t358 + t381) * t357;
t274 = pkin(4) * t316 + t468;
t272 = -pkin(5) * t452 - t278;
t271 = -pkin(4) * t434 - t397;
t270 = -qJ(5) * t434 - t441;
t266 = t336 + (pkin(5) * t357 - t324) * t356 - t403;
t263 = t316 * t471 + t468;
t262 = t395 - t476;
t259 = -pkin(4) * t410 - qJD(5) * t449 + t386;
t255 = -pkin(4) * t431 - t374;
t252 = -qJ(5) * t431 + (t354 * t425 - qJD(5)) * t355 - t411;
t250 = (qJ(6) * qJD(4) - qJD(5)) * t449 + (qJD(6) * t357 - t355 * t401) * t354 + t386;
t249 = (-pkin(5) * t427 + qJ(5) * qJD(3)) * t357 + (qJD(5) + (pkin(5) * qJD(3) - t425) * t354) * t355 + t411;
t246 = -pkin(5) * t406 + qJD(3) * t365 - qJD(6) * t355 - t374;
t242 = -t401 * t434 + t362;
t1 = [-0.2e1 * MDP(7) * t385 + 0.2e1 * t420 * t481 + (-t355 * t447 + (qJ(2) * t431 + qJD(2) * t355) * t416) * MDP(12) + (-t357 * t447 + (-qJ(2) * t432 + qJD(2) * t357) * t416) * MDP(13) + (-t280 * t449 - t316 * t364) * MDP(14) + ((t314 * t356 + t316 * t354) * t432 + (t460 - t281 * t356 + (-t316 * t356 + t457) * qJD(4)) * t357) * MDP(15) + (-t344 * t406 - t280 * t355 + (t316 * t357 + t356 * t375) * qJD(3)) * MDP(16) + (-t344 * t404 - t281 * t355 + (-t314 * t357 - t354 * t375) * qJD(3)) * MDP(17) + (t344 + t436) * MDP(18) * t431 + (t374 * t344 + t387 * t355 + (-t358 * t281 + t305 * t427) * t357 + ((t394 * qJD(1) - t267) * t357 + (t358 * t314 + t398 * t354) * t355) * qJD(3)) * MDP(19) + (-(-t354 * t405 + t411) * t344 - t388 * t355 + (t358 * t280 - t305 * t428) * t357 + ((-t439 * qJD(1) - t268) * t357 + (t358 * t316 + t398 * t356) * t355) * qJD(3)) * MDP(20) + (t252 * t314 + t255 * t316 + t278 * t281 - t279 * t280 - t377 * t432 + (-qJD(4) * t376 + t245 * t354 + t247 * t356) * t357) * MDP(21) + (t255 * t344 - t259 * t314 - t281 * t284 + (t265 * t433 + t247) * t355 + (-t265 * t427 - t465 + (qJD(1) * t279 + t260) * qJD(3)) * t357) * MDP(22) + (-t252 * t344 - t259 * t316 + t280 * t284 + (t265 * t424 - t245) * t355 + (t265 * t428 - t464 + (-qJD(1) * t278 - t261) * qJD(3)) * t357) * MDP(23) + (t245 * t278 + t247 * t279 + t248 * t284 + t252 * t261 + t255 * t260 + t259 * t265) * MDP(24) + (t246 * t316 - t249 * t314 - t266 * t280 - t272 * t281 - t379 * t432 + (-qJD(4) * t378 + t242 * t356 - t244 * t354) * t357) * MDP(25) + (t249 * t344 - t250 * t316 + t277 * t280 + (t254 * t424 + t244) * t355 + ((qJD(1) * t272 + t253) * qJD(3) + t369) * t357) * MDP(26) + (-t246 * t344 + t250 * t314 + t277 * t281 + (-t254 * t433 - t242) * t355 + ((-qJD(1) * t266 - t251) * qJD(3) + t370) * t357) * MDP(27) + (t242 * t266 + t243 * t277 + t244 * t272 + t246 * t251 + t249 * t253 + t250 * t254) * MDP(28) + t480 * qJD(2) * t416 + (-MDP(10) * t357 - MDP(9) * t355) * t359; -t322 * MDP(20) - t480 * t360 + t418 * t314 * t432 + t361 * qJD(1) + (t437 * MDP(13) - t248 * MDP(24) - t243 * MDP(28) + t400 * t281 + t399 * t280 + (t376 * MDP(24) + t378 * MDP(28) - t314 * t390 + (-t356 * t399 + t383) * t344) * qJD(3)) * t357 + (t437 * MDP(12) + t380 * MDP(24) + (t242 * t354 + t244 * t356) * MDP(28) - t281 * t390 + (t316 * MDP(20) - t314 * MDP(22) + t265 * MDP(24) + t254 * MDP(28) + t383 * t434) * qJD(3) + t361 * qJD(4)) * t355 + t474 * (-t316 * t432 + t322) + t417 * (-t280 * t453 + (t355 * t427 + t409 + t435) * t316); t355 * MDP(7) * t448 - t360 * t481 + (t316 * t454 - t460) * MDP(14) + ((-t280 - t458) * t356 + (-t281 - t456) * t354) * MDP(15) + (t344 * t427 + (t344 * t451 + t357 * t392) * qJD(1)) * MDP(16) + (-t344 * t428 + (-t344 * t453 + t357 * t393) * qJD(1)) * MDP(17) - t344 * MDP(18) * t434 + (-pkin(3) * t281 - t397 * t344 - t393 * t323 + (t305 * t354 - t414) * qJD(4) + (t267 * t357 + t354 * t371) * qJD(1)) * MDP(19) + (pkin(3) * t280 + t441 * t344 + t392 * t323 + (t305 * t356 + t415) * qJD(4) + (t268 * t357 + t356 * t371) * qJD(1)) * MDP(20) + (-t270 * t314 - t271 * t316 + (-t245 + t344 * t260 + (-t281 + t429) * pkin(8)) * t356 + (t247 + t461 + (qJD(4) * t314 - t280) * pkin(8)) * t354) * MDP(21) + (t464 - t271 * t344 - t281 * t325 + t443 * t314 + (-t265 * t354 + t414) * qJD(4) + (-t260 * t357 + t354 * t372) * qJD(1)) * MDP(22) + (-t465 + t270 * t344 + t280 * t325 + t443 * t316 + (-t265 * t356 - t415) * qJD(4) + (t261 * t357 + t356 * t372) * qJD(1)) * MDP(23) + (t248 * t325 - t260 * t271 - t261 * t270 - t443 * t265 + (qJD(4) * t377 + t380) * pkin(8)) * MDP(24) + (-t280 * t329 - t281 * t330 - t446 * t316 + t445 * t314 + t479 * t356 + (t242 - t462) * t354) * MDP(25) + (t280 * t308 - t445 * t344 + t444 * t316 + (-t254 * t451 + (qJD(3) * t330 - t253) * t357) * qJD(1) - t370) * MDP(26) + (t281 * t308 + t446 * t344 - t444 * t314 + (t254 * t453 + (-qJD(3) * t329 + t251) * t357) * qJD(1) + t369) * MDP(27) + (t242 * t329 + t243 * t308 + t244 * t330 - t251 * t446 - t253 * t445 - t254 * t444) * MDP(28) + (MDP(13) * t355 * t360 - MDP(12) * t448) * qJ(2); t262 * MDP(16) + (-t354 * t419 + t332) * MDP(17) - t388 * MDP(20) + (pkin(4) * t280 - t469) * MDP(21) + (t340 + t382) * MDP(23) + (-pkin(4) * t247 - qJ(5) * t245 - t260 * t268 - t261 * t422 - t265 * t274) * MDP(24) + (t280 * t471 - t469) * MDP(25) + (0.2e1 * t331 + t340 + t366) * MDP(26) + t367 * MDP(27) + (qJ(5) * t244 - t242 * t471 + t251 * t421 + t253 * t423 - t254 * t263) * MDP(28) + (-t267 * MDP(20) + t422 * MDP(23) + t373 * MDP(26) + (0.2e1 * qJD(6) + t257) * MDP(27) + t475 * t268) * t344 + (-MDP(17) * t427 + (-0.2e1 * pkin(4) * MDP(22) + 0.2e1 * t471 * MDP(27) + MDP(18)) * qJD(3)) * t434 + (t344 * MDP(17) - t305 * MDP(19) + (-t261 - t268) * MDP(21) + t265 * MDP(22) + t274 * MDP(23) + (t253 + t421) * MDP(25) + t263 * MDP(26) - t254 * MDP(27) + MDP(15) * t316) * t316 + (t316 * MDP(14) + t305 * MDP(20) + (t260 - t422) * MDP(21) + t274 * MDP(22) - t265 * MDP(23) + (t251 - t423) * MDP(25) - t254 * MDP(26) - t263 * MDP(27) - MDP(15) * t314) * t314 + t475 * t387; t395 * MDP(21) + (t265 * t316 - t387 + t461) * MDP(24) + t262 * MDP(25) + (t254 * t316 + t362 - t462) * MDP(28) + (-MDP(21) * t406 + (-MDP(21) * t451 + (-pkin(4) * MDP(24) - MDP(28) * t471) * t357) * qJD(3)) * qJD(1) + (-MDP(22) + MDP(27)) * (-t347 + t459) + t474 * (-t316 ^ 2 - t341); (t456 - t281) * MDP(25) + (t347 + t459) * MDP(26) + (-t314 ^ 2 - t341) * MDP(27) + (-t254 * t314 + t479) * MDP(28);];
tauc  = t1;
