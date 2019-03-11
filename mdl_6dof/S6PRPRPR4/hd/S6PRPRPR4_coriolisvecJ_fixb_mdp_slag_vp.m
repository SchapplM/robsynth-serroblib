% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:29
% EndTime: 2019-03-08 19:41:36
% DurationCPUTime: 4.86s
% Computational Cost: add. (3313->373), mult. (8800->529), div. (0->0), fcn. (7201->12), ass. (0->159)
t379 = cos(pkin(11));
t452 = cos(qJ(4));
t412 = t452 * t379;
t367 = qJD(2) * t412;
t376 = sin(pkin(11));
t382 = sin(qJ(4));
t434 = t382 * t376;
t409 = qJD(2) * t434;
t343 = -t367 + t409;
t339 = qJD(6) + t343;
t393 = t412 - t434;
t348 = t393 * qJD(4);
t356 = t452 * t376 + t382 * t379;
t349 = t356 * qJD(4);
t383 = sin(qJ(2));
t377 = sin(pkin(6));
t421 = qJD(1) * t377;
t411 = t383 * t421;
t462 = pkin(4) * t349 - qJ(5) * t348 - qJD(5) * t356 - t411;
t385 = cos(qJ(2));
t435 = t377 * t385;
t388 = t393 * t435;
t450 = pkin(8) + qJ(3);
t360 = t450 * t376;
t362 = t450 * t379;
t454 = -t452 * t360 - t382 * t362;
t461 = -qJD(1) * t388 + t393 * qJD(3) + t454 * qJD(4);
t358 = qJD(2) * qJ(3) + t411;
t380 = cos(pkin(6));
t420 = qJD(1) * t380;
t366 = t379 * t420;
t448 = pkin(8) * qJD(2);
t319 = t366 + (-t358 - t448) * t376;
t332 = t379 * t358 + t376 * t420;
t320 = t379 * t448 + t332;
t279 = t382 * t319 + t452 * t320;
t460 = qJD(4) * t279;
t345 = qJD(2) * t356;
t375 = sin(pkin(12));
t378 = cos(pkin(12));
t326 = -t378 * qJD(4) + t345 * t375;
t384 = cos(qJ(6));
t459 = t384 * t326;
t381 = sin(qJ(6));
t355 = t375 * t384 + t378 * t381;
t423 = t339 * t355;
t328 = qJD(4) * t375 + t345 * t378;
t458 = t326 * t381 - t328 * t384;
t410 = t385 * t421;
t352 = (qJD(3) + t410) * qJD(2);
t457 = t393 * t352;
t397 = (-t358 * t376 + t366) * t376 - t332 * t379;
t456 = t397 * t385;
t422 = t376 ^ 2 + t379 ^ 2;
t455 = t422 * MDP(7);
t429 = -t461 * t375 + t378 * t462;
t428 = t375 * t462 + t461 * t378;
t323 = -t382 * t360 + t452 * t362;
t389 = t356 * t435;
t426 = -qJD(1) * t389 + qJD(3) * t356 + qJD(4) * t323;
t395 = -t452 * t319 + t382 * t320;
t403 = qJD(3) - t410;
t337 = qJD(2) * t349;
t432 = t384 * t378;
t437 = t375 * t381;
t353 = -t432 + t437;
t424 = t339 * t353;
t453 = -t337 * t355 + t424 * t339;
t340 = t343 ^ 2;
t451 = pkin(9) * t378;
t449 = pkin(9) + qJ(5);
t447 = qJD(2) * pkin(2);
t281 = t328 * t381 + t459;
t446 = t281 * t345;
t445 = t458 * t345;
t364 = qJD(4) * t367;
t336 = -qJD(4) * t409 + t364;
t444 = t336 * t375;
t443 = t336 * t378;
t441 = t343 * t375;
t440 = t348 * t375;
t439 = t356 * t375;
t438 = t356 * t378;
t436 = t377 * t383;
t431 = pkin(5) * t349 - t348 * t451 + t429;
t430 = pkin(9) * t440 - t428;
t251 = t457 + (qJD(5) - t395) * qJD(4);
t408 = qJD(2) * t436;
t363 = qJD(1) * t408;
t277 = pkin(4) * t337 - qJ(5) * t336 - qJD(5) * t345 + t363;
t240 = t378 * t251 + t375 * t277;
t274 = qJD(4) * qJ(5) + t279;
t371 = -pkin(3) * t379 - pkin(2);
t338 = t371 * qJD(2) + t403;
t288 = pkin(4) * t343 - qJ(5) * t345 + t338;
t246 = t378 * t274 + t375 * t288;
t308 = pkin(4) * t345 + qJ(5) * t343;
t253 = t375 * t308 - t378 * t395;
t427 = pkin(5) * t440 + t426;
t309 = -pkin(4) * t393 - qJ(5) * t356 + t371;
t272 = t375 * t309 + t378 * t323;
t417 = qJD(6) * t384;
t425 = -t326 * t417 + t336 * t432;
t242 = -pkin(9) * t326 + t246;
t419 = qJD(6) * t242;
t418 = qJD(6) * t356;
t407 = t422 * t352;
t239 = -t251 * t375 + t378 * t277;
t245 = -t274 * t375 + t378 * t288;
t252 = t378 * t308 + t375 * t395;
t271 = t378 * t309 - t323 * t375;
t238 = -pkin(9) * t444 + t240;
t241 = pkin(5) * t343 - pkin(9) * t328 + t245;
t405 = -qJD(6) * t241 - t238;
t256 = t352 * t356 + t460;
t404 = -t353 * t337 - t339 * t423;
t234 = t241 * t384 - t242 * t381;
t235 = t241 * t381 + t242 * t384;
t402 = -t245 * t375 + t246 * t378;
t259 = -pkin(5) * t393 - pkin(9) * t438 + t271;
t261 = -pkin(9) * t439 + t272;
t401 = t259 * t384 - t261 * t381;
t400 = t259 * t381 + t261 * t384;
t341 = -t376 * t436 + t379 * t380;
t342 = t376 * t380 + t379 * t436;
t301 = t382 * t341 + t452 * t342;
t289 = -t301 * t375 - t378 * t435;
t290 = t301 * t378 - t375 * t435;
t399 = t289 * t384 - t290 * t381;
t398 = t289 * t381 + t290 * t384;
t396 = -qJD(6) * t328 - t444;
t394 = t452 * t341 - t382 * t342;
t361 = t449 * t378;
t392 = pkin(5) * t345 + qJD(5) * t375 + qJD(6) * t361 + t343 * t451 + t252;
t359 = t449 * t375;
t391 = pkin(9) * t441 - qJD(5) * t378 + qJD(6) * t359 + t253;
t270 = -qJD(4) * pkin(4) + qJD(5) + t395;
t390 = t256 * t356 + t270 * t348 - t336 * t454;
t387 = -pkin(4) * t336 - qJ(5) * t337 + (-qJD(5) + t270) * t343;
t255 = -qJD(6) * t458 + t336 * t355;
t386 = qJD(2) ^ 2;
t370 = -pkin(5) * t378 - pkin(4);
t357 = t403 - t447;
t305 = t353 * t356;
t304 = t355 * t356;
t295 = pkin(5) * t439 - t454;
t276 = qJD(2) * t389 + qJD(4) * t301;
t275 = qJD(2) * t388 + qJD(4) * t394;
t266 = t275 * t378 + t375 * t408;
t265 = -t275 * t375 + t378 * t408;
t264 = t348 * t355 + t417 * t438 - t418 * t437;
t263 = -t348 * t353 - t355 * t418;
t262 = -pkin(5) * t441 + t279;
t260 = t326 * pkin(5) + t270;
t254 = t381 * t396 + t425;
t247 = pkin(5) * t444 + t256;
t237 = pkin(5) * t337 - pkin(9) * t443 + t239;
t236 = t384 * t237;
t1 = [(t265 * t343 + t276 * t326 + t289 * t337 - t394 * t444) * MDP(16) + (-t266 * t343 + t276 * t328 - t290 * t337 - t394 * t443) * MDP(17) + (-t265 * t328 - t266 * t326 + (-t289 * t378 - t290 * t375) * t336) * MDP(18) + (t239 * t289 + t240 * t290 + t245 * t265 + t246 * t266 - t256 * t394 + t270 * t276) * MDP(19) + ((-qJD(6) * t398 + t265 * t384 - t266 * t381) * t339 + t399 * t337 + t276 * t281 - t394 * t255) * MDP(25) + (-(qJD(6) * t399 + t265 * t381 + t266 * t384) * t339 - t398 * t337 - t276 * t458 - t394 * t254) * MDP(26) + (-t341 * t376 + t342 * t379) * MDP(8) * t352 + (-MDP(14) * t276 - MDP(15) * t275) * qJD(4) + ((-t337 * MDP(14) - t336 * MDP(15)) * t385 + (-MDP(8) * t456 + (t343 * MDP(14) + t345 * MDP(15) + (t357 - t410) * MDP(8)) * t383) * qJD(2) + ((-MDP(4) + t455) * t385 + (-MDP(5) * t379 + MDP(6) * t376 - MDP(3)) * t383) * t386) * t377; (t403 * qJD(2) * t422 + t407) * MDP(7) + (-t397 * qJD(3) + qJ(3) * t407 + (t456 + (-t357 - t447) * t383) * t421) * MDP(8) + (t336 * t356 + t345 * t348) * MDP(9) + (t336 * t393 - t337 * t356 - t343 * t348 - t345 * t349) * MDP(10) + (t337 * t371 + t338 * t349 + (-qJD(2) * t393 - t343) * t411) * MDP(14) + (t336 * t371 + t338 * t348) * MDP(15) + (-t239 * t393 + t245 * t349 + t271 * t337 + t426 * t326 + t429 * t343 + t390 * t375) * MDP(16) + (t240 * t393 - t246 * t349 - t272 * t337 + t426 * t328 - t428 * t343 + t390 * t378) * MDP(17) + (-t429 * t328 - t428 * t326 + (-t239 * t356 - t245 * t348 - t271 * t336) * t378 + (-t240 * t356 - t246 * t348 - t272 * t336) * t375) * MDP(18) + (t239 * t271 + t240 * t272 + t429 * t245 + t428 * t246 - t256 * t454 + t426 * t270) * MDP(19) + (-t254 * t305 - t263 * t458) * MDP(20) + (-t254 * t304 + t255 * t305 - t263 * t281 + t264 * t458) * MDP(21) + (-t254 * t393 + t263 * t339 - t305 * t337 - t349 * t458) * MDP(22) + (t255 * t393 - t264 * t339 - t281 * t349 - t304 * t337) * MDP(23) + (-t337 * t393 + t339 * t349) * MDP(24) + (t401 * t337 - (-t238 * t381 + t236) * t393 + t234 * t349 + t295 * t255 + t247 * t304 + t260 * t264 + (t430 * t381 + t431 * t384) * t339 + t427 * t281 + (t235 * t393 - t339 * t400) * qJD(6)) * MDP(25) + (-t400 * t337 + (t237 * t381 + t238 * t384) * t393 - t235 * t349 + t295 * t254 - t247 * t305 + t260 * t263 + (-t431 * t381 + t430 * t384) * t339 - t427 * t458 + (t234 * t393 - t339 * t401) * qJD(6)) * MDP(26) + (t348 * MDP(11) - t349 * MDP(12) - t426 * MDP(14) - MDP(15) * t461) * qJD(4); -t386 * t455 + (t397 * qJD(2) + t363) * MDP(8) + 0.2e1 * t345 * qJD(4) * MDP(14) + (t364 + (-t343 - t409) * qJD(4)) * MDP(15) + (-t326 * t345 + t337 * t378 - t340 * t375) * MDP(16) + (-t328 * t345 - t337 * t375 - t340 * t378) * MDP(17) + ((-t326 * t378 + t328 * t375) * t343 + (-t375 ^ 2 - t378 ^ 2) * t336) * MDP(18) + (t239 * t378 + t240 * t375 - t270 * t345 + t343 * t402) * MDP(19) + (t404 - t446) * MDP(25) + (t445 + t453) * MDP(26); -t340 * MDP(10) + (t364 + (t343 - t409) * qJD(4)) * MDP(11) + (-t256 + t460) * MDP(14) + (t338 * t343 - t457) * MDP(15) + (-t252 * t343 - t256 * t378 - t279 * t326 + t375 * t387) * MDP(16) + (t253 * t343 + t256 * t375 - t279 * t328 + t378 * t387) * MDP(17) + (t252 * t328 + t253 * t326 + (-qJD(5) * t326 - t245 * t343 + t240) * t378 + (qJD(5) * t328 - t246 * t343 - t239) * t375) * MDP(18) + (-pkin(4) * t256 - t245 * t252 - t246 * t253 - t270 * t279 + t402 * qJD(5) + (-t239 * t375 + t240 * t378) * qJ(5)) * MDP(19) + (t254 * t355 + t424 * t458) * MDP(20) + (-t254 * t353 - t255 * t355 + t424 * t281 + t423 * t458) * MDP(21) + (t445 - t453) * MDP(22) + (t404 + t446) * MDP(23) + ((-t359 * t384 - t361 * t381) * t337 + t370 * t255 + t247 * t353 - t262 * t281 + (t381 * t391 - t384 * t392) * t339 + t423 * t260) * MDP(25) + (-(-t359 * t381 + t361 * t384) * t337 + t370 * t254 + t247 * t355 + t262 * t458 + (t381 * t392 + t384 * t391) * t339 - t424 * t260) * MDP(26) + (MDP(10) * t345 - t338 * MDP(14) - t245 * MDP(16) + t246 * MDP(17) - t339 * MDP(24) - t234 * MDP(25) + t235 * MDP(26) + t343 * MDP(9)) * t345; (t328 * t343 + t444) * MDP(16) + (-t326 * t343 + t443) * MDP(17) + (-t326 ^ 2 - t328 ^ 2) * MDP(18) + (t245 * t328 + t246 * t326 + t256) * MDP(19) + (-t339 * t458 + t255) * MDP(25) + (-t339 * t459 + (-t328 * t339 + t396) * t381 + t425) * MDP(26); -t281 ^ 2 * MDP(21) + (t281 * t339 + t425) * MDP(22) + t337 * MDP(24) + (t235 * t339 + t236) * MDP(25) + (t234 * t339 + t260 * t281) * MDP(26) - (MDP(20) * t281 - MDP(21) * t458 + t339 * MDP(23) - t260 * MDP(25)) * t458 + (MDP(23) * t396 - MDP(25) * t419 + MDP(26) * t405) * t384 + (t396 * MDP(22) + (qJD(6) * t326 - t443) * MDP(23) + t405 * MDP(25) + (-t237 + t419) * MDP(26)) * t381;];
tauc  = t1;
