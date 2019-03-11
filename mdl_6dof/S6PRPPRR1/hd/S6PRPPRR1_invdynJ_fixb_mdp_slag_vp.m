% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:22
% EndTime: 2019-03-08 19:16:28
% DurationCPUTime: 4.06s
% Computational Cost: add. (2351->375), mult. (5560->514), div. (0->0), fcn. (4865->16), ass. (0->185)
t382 = sin(qJ(2));
t385 = cos(qJ(2));
t375 = sin(pkin(6));
t440 = qJD(2) * t375;
t420 = qJD(1) * t440;
t432 = qJDD(1) * t375;
t477 = t382 * t432 + t385 * t420;
t372 = sin(pkin(12));
t376 = cos(pkin(12));
t452 = t376 * MDP(6);
t411 = MDP(7) * t372 - t452;
t381 = sin(qJ(5));
t384 = cos(qJ(5));
t331 = t372 * t381 - t384 * t376;
t333 = t372 * t384 + t376 * t381;
t424 = -pkin(4) * t376 - pkin(3);
t377 = cos(pkin(11));
t471 = pkin(2) * t377;
t344 = t424 - t471;
t284 = pkin(5) * t331 - pkin(9) * t333 + t344;
t327 = t333 * qJD(5);
t430 = qJDD(2) * t384;
t353 = t376 * t430;
t431 = qJDD(2) * t381;
t412 = -t372 * t431 + t353;
t295 = qJD(2) * t327 - t412;
t291 = qJDD(6) + t295;
t371 = pkin(12) + qJ(5);
t365 = cos(t371);
t374 = sin(pkin(10));
t378 = cos(pkin(10));
t373 = sin(pkin(11));
t448 = t385 * t377;
t332 = t373 * t382 - t448;
t379 = cos(pkin(6));
t396 = t332 * t379;
t403 = t373 * t385 + t377 * t382;
t279 = -t374 * t403 - t378 * t396;
t282 = t374 * t396 - t378 * t403;
t454 = t375 * t382;
t319 = t373 * t454 - t375 * t448;
t393 = g(1) * t282 + g(2) * t279 - g(3) * t319;
t391 = t393 * t365;
t476 = t284 * t291 - t391;
t441 = qJD(1) * t375;
t423 = t382 * t441;
t345 = t373 * t423;
t422 = t385 * t441;
t317 = t377 * t422 - t345;
t433 = qJD(4) - t317;
t437 = qJD(2) * t384;
t354 = t376 * t437;
t438 = qJD(2) * t381;
t421 = t372 * t438;
t323 = -t354 + t421;
t322 = qJD(6) + t323;
t442 = t372 ^ 2 + t376 ^ 2;
t475 = MDP(8) * t442;
t320 = t403 * t375;
t302 = -t320 * t372 + t376 * t379;
t303 = t320 * t376 + t372 * t379;
t258 = t302 * t381 + t303 * t384;
t383 = cos(qJ(6));
t313 = t319 * t383;
t380 = sin(qJ(6));
t474 = -t258 * t380 + t313;
t325 = t333 * qJD(2);
t337 = qJD(2) * pkin(2) + t422;
t346 = t377 * t423;
t306 = t373 * t337 + t346;
t304 = qJD(2) * qJ(4) + t306;
t357 = qJD(1) * t379 + qJD(3);
t343 = t376 * t357;
t271 = t343 + (-pkin(8) * qJD(2) - t304) * t372;
t276 = t376 * t304 + t372 * t357;
t439 = qJD(2) * t376;
t272 = pkin(8) * t439 + t276;
t248 = t271 * t381 + t272 * t384;
t351 = t385 * t432;
t318 = qJDD(2) * pkin(2) - t382 * t420 + t351;
t274 = t373 * t318 + t377 * t477;
t269 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t274;
t355 = t379 * qJDD(1) + qJDD(3);
t339 = t376 * t355;
t468 = pkin(8) * qJDD(2);
t250 = t339 + (-t269 - t468) * t372;
t254 = t376 * t269 + t372 * t355;
t251 = t376 * t468 + t254;
t409 = -t250 * t384 + t251 * t381;
t238 = -qJDD(5) * pkin(5) + qJD(5) * t248 + t409;
t321 = t403 * t379;
t280 = t321 * t378 - t332 * t374;
t281 = t374 * t321 + t332 * t378;
t364 = sin(t371);
t455 = t375 * t378;
t457 = t374 * t375;
t395 = g(1) * (t281 * t364 + t365 * t457) + g(2) * (-t280 * t364 - t365 * t455) + g(3) * (-t320 * t364 + t365 * t379);
t473 = t322 * (pkin(5) * t325 + pkin(9) * t322) + t238 + t395;
t247 = t271 * t384 - t272 * t381;
t408 = t250 * t381 + t251 * t384;
t237 = qJDD(5) * pkin(9) + qJD(5) * t247 + t408;
t245 = -qJD(5) * pkin(5) - t247;
t305 = t337 * t377 - t345;
t413 = qJD(4) - t305;
t290 = qJD(2) * t424 + t413;
t252 = pkin(5) * t323 - pkin(9) * t325 + t290;
t358 = pkin(2) * t373 + qJ(4);
t470 = pkin(8) + t358;
t328 = t470 * t372;
t329 = t470 * t376;
t287 = -t328 * t381 + t329 * t384;
t326 = t331 * qJD(5);
t394 = g(1) * t281 - g(2) * t280 - g(3) * t320;
t404 = -t328 * t384 - t329 * t381;
t445 = -qJD(5) * t404 + t331 * t433;
t472 = -(qJD(6) * t252 + t237) * t331 + t238 * t333 - t245 * t326 + (-qJD(6) * t284 + t445) * t322 - t287 * t291 + t394;
t425 = qJD(5) * t354 + t372 * t430 + t376 * t431;
t294 = -qJD(5) * t421 + t425;
t434 = t383 * qJD(5);
t426 = qJD(6) * t434 + t380 * qJDD(5) + t383 * t294;
t435 = qJD(6) * t380;
t255 = -t325 * t435 + t426;
t467 = t255 * t380;
t275 = -t304 * t372 + t343;
t465 = t275 * t372;
t458 = t325 * t380;
t307 = -t434 + t458;
t463 = t307 * t322;
t462 = t307 * t325;
t309 = qJD(5) * t380 + t325 * t383;
t461 = t309 * t322;
t460 = t309 * t325;
t459 = t319 * t380;
t456 = t374 * t382;
t453 = t375 * t385;
t451 = t379 * t382;
t450 = t379 * t385;
t449 = t380 * t291;
t285 = t383 * t291;
t447 = qJDD(1) - g(3);
t446 = t255 * t331 + t309 * t327;
t444 = qJD(5) * t287 + t333 * t433;
t314 = t373 * t422 + t346;
t443 = pkin(5) * t327 + pkin(9) * t326 - t314;
t436 = qJD(6) * t333;
t429 = t333 * t449;
t428 = t333 * t285;
t427 = t378 * t450;
t417 = -t383 * qJDD(5) + t294 * t380;
t416 = t322 * t383;
t273 = t318 * t377 - t373 * t477;
t246 = qJD(5) * pkin(9) + t248;
t240 = t246 * t383 + t252 * t380;
t410 = t246 * t380 - t252 * t383;
t256 = qJD(6) * t309 + t417;
t407 = -t256 * t331 - t307 * t327;
t406 = t258 * t383 + t459;
t405 = t302 * t384 - t303 * t381;
t402 = t285 + (-t323 * t380 - t435) * t322;
t401 = qJDD(4) - t273;
t400 = -t374 * t450 - t378 * t382;
t399 = t326 * t380 - t383 * t436;
t398 = t326 * t383 + t333 * t435;
t397 = -g(1) * t457 + g(2) * t455 - g(3) * t379;
t390 = -pkin(9) * t291 + (t245 + t247) * t322;
t389 = -g(1) * t400 - g(3) * t453;
t259 = qJDD(2) * t424 + t401;
t386 = qJD(2) ^ 2;
t361 = -pkin(3) - t471;
t347 = pkin(2) * t427;
t316 = t332 * t440;
t315 = qJD(2) * t320;
t301 = -qJD(2) * pkin(3) + t413;
t299 = t320 * t365 + t364 * t379;
t297 = -qJD(5) * t327 - qJDD(5) * t331;
t296 = -qJD(5) * t326 + qJDD(5) * t333;
t270 = -qJDD(2) * pkin(3) + t401;
t264 = -t281 * t365 + t364 * t457;
t262 = t280 * t365 - t364 * t455;
t253 = -t269 * t372 + t339;
t244 = qJD(5) * t258 - t316 * t333;
t243 = qJD(5) * t405 + t316 * t331;
t242 = pkin(5) * t295 - pkin(9) * t294 + t259;
t241 = t383 * t242;
t1 = [t447 * MDP(1) + (-t273 * t319 + t274 * t320 - t305 * t315 + t355 * t379 - g(3)) * MDP(5) + (t253 * t302 + t254 * t303 + t270 * t319 + t301 * t315 - g(3)) * MDP(9) + (-qJD(5) * t244 + qJDD(5) * t405 + t295 * t319 + t315 * t323) * MDP(15) + (-qJD(5) * t243 - qJDD(5) * t258 + t294 * t319 + t315 * t325) * MDP(16) + ((-qJD(6) * t406 - t243 * t380 + t315 * t383) * t322 + t474 * t291 + t244 * t307 - t405 * t256) * MDP(22) + (-(qJD(6) * t474 + t243 * t383 + t315 * t380) * t322 - t406 * t291 + t244 * t309 - t405 * t255) * MDP(23) + t411 * t315 * qJD(2) + ((-t302 * t372 + t303 * t376) * MDP(8) + t411 * t319) * qJDD(2) - (t306 * MDP(5) + (t276 * t376 - t465) * MDP(9) + qJD(2) * t475) * t316 + ((qJDD(2) * t385 - t382 * t386) * MDP(3) + (-qJDD(2) * t382 - t385 * t386) * MDP(4)) * t375; qJDD(2) * MDP(2) + (t351 - g(2) * (t427 - t456) + t389) * MDP(3) + (-g(1) * (t374 * t451 - t378 * t385) - g(2) * (-t374 * t385 - t378 * t451) - t447 * t454) * MDP(4) + (-g(2) * t347 + t305 * t314 - t306 * t317 + (g(2) * t456 + t273 * t377 + t274 * t373 + t389) * pkin(2)) * MDP(5) + (-t253 * t372 + t254 * t376 + t394 + (qJD(2) * t433 + t358 * qJDD(2)) * t442) * MDP(8) + (t270 * t361 - t301 * t314 - g(1) * (pkin(2) * t400 + pkin(3) * t282 - qJ(4) * t281) - g(2) * (-pkin(2) * t456 + pkin(3) * t279 + qJ(4) * t280 + t347) - g(3) * (pkin(2) * t453 - pkin(3) * t319 + qJ(4) * t320) + (t254 * t358 + t276 * t433) * t376 + (-t253 * t358 - t275 * t433) * t372) * MDP(9) + (t294 * t333 - t325 * t326) * MDP(10) + (-t294 * t331 - t295 * t333 + t323 * t326 - t325 * t327) * MDP(11) + t296 * MDP(12) + t297 * MDP(13) + (-qJD(5) * t444 + qJDD(5) * t404 + t259 * t331 + t290 * t327 + t295 * t344 - t314 * t323 - t391) * MDP(15) + (qJD(5) * t445 - qJDD(5) * t287 + t259 * t333 - t290 * t326 + t294 * t344 - t314 * t325 + t364 * t393) * MDP(16) + (t255 * t333 * t383 - t309 * t398) * MDP(17) + (-(-t307 * t383 - t309 * t380) * t326 + (-t467 - t256 * t383 + (t307 * t380 - t309 * t383) * qJD(6)) * t333) * MDP(18) + (-t322 * t398 + t428 + t446) * MDP(19) + (t322 * t399 + t407 - t429) * MDP(20) + (t291 * t331 + t322 * t327) * MDP(21) + (-t410 * t327 + t241 * t331 - t404 * t256 + t444 * t307 + (t443 * t322 + (t245 * t333 - t246 * t331 - t287 * t322) * qJD(6) + t476) * t383 + t472 * t380) * MDP(22) + (-t240 * t327 - t404 * t255 + t444 * t309 + (-(-qJD(6) * t246 + t242) * t331 - t245 * t436 + (qJD(6) * t287 - t443) * t322 - t476) * t380 + t472 * t383) * MDP(23) + t411 * (-qJD(2) * t314 + qJDD(2) * t361 + t270 + t393); (t397 + t355) * MDP(5) + (t253 * t376 + t254 * t372 + t397) * MDP(9) + t297 * MDP(15) - t296 * MDP(16) + (-t407 - t429) * MDP(22) + (-t428 + t446) * MDP(23) + (MDP(22) * t399 + MDP(23) * t398) * t322; (qJD(2) * t465 - t276 * t439 + t393 + t401) * MDP(9) - t353 * MDP(15) + t425 * MDP(16) + (t402 - t462) * MDP(22) + (-t322 ^ 2 * t383 - t449 - t460) * MDP(23) - t386 * t475 + (-t452 - pkin(3) * MDP(9) + (MDP(15) * t381 + MDP(7)) * t372) * qJDD(2) + ((t372 * t437 + t376 * t438 + t325) * MDP(15) + (-t323 - t421) * MDP(16)) * qJD(5); -t323 ^ 2 * MDP(11) + ((t323 - t421) * qJD(5) + t425) * MDP(12) + t412 * MDP(13) + qJDD(5) * MDP(14) + (-t395 - t409) * MDP(15) + (g(1) * t264 + g(2) * t262 + g(3) * t299 + t290 * t323 - t408) * MDP(16) + (t309 * t416 + t467) * MDP(17) + ((t255 - t463) * t383 + (-t256 - t461) * t380) * MDP(18) + (t322 * t416 + t449 - t460) * MDP(19) + (t402 + t462) * MDP(20) + (-pkin(5) * t256 - t248 * t307 + t390 * t380 - t383 * t473) * MDP(22) + (-pkin(5) * t255 - t248 * t309 + t380 * t473 + t390 * t383) * MDP(23) + (MDP(10) * t323 + t325 * MDP(11) - t290 * MDP(15) - t322 * MDP(21) + MDP(22) * t410 + t240 * MDP(23)) * t325; t309 * t307 * MDP(17) + (-t307 ^ 2 + t309 ^ 2) * MDP(18) + (t426 + t463) * MDP(19) + (-t417 + t461) * MDP(20) + t291 * MDP(21) + (-t380 * t237 + t241 + t240 * t322 - t245 * t309 - g(1) * (-t264 * t380 - t282 * t383) - g(2) * (-t262 * t380 - t279 * t383) - g(3) * (-t299 * t380 + t313)) * MDP(22) + (-t383 * t237 - t380 * t242 - t410 * t322 + t245 * t307 - g(1) * (-t264 * t383 + t282 * t380) - g(2) * (-t262 * t383 + t279 * t380) - g(3) * (-t299 * t383 - t459)) * MDP(23) + (-MDP(19) * t458 - MDP(20) * t309 - MDP(22) * t240 + MDP(23) * t410) * qJD(6);];
tau  = t1;
