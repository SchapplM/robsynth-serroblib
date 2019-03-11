% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:41
% EndTime: 2019-03-08 19:49:47
% DurationCPUTime: 3.76s
% Computational Cost: add. (1884->352), mult. (4399->511), div. (0->0), fcn. (3154->10), ass. (0->161)
t328 = sin(pkin(11));
t330 = cos(pkin(11));
t332 = sin(qJ(6));
t335 = cos(qJ(6));
t432 = -t328 * t332 + t335 * t330;
t429 = t432 * qJD(6);
t333 = sin(qJ(4));
t397 = qJD(2) * t333;
t322 = qJD(6) + t397;
t336 = cos(qJ(4));
t396 = qJD(2) * t336;
t383 = t328 * t396;
t389 = t330 * qJD(4);
t298 = t383 - t389;
t380 = t330 * t396;
t394 = qJD(4) * t328;
t300 = t380 + t394;
t354 = t298 * t332 - t300 * t335;
t437 = t322 * t354;
t436 = MDP(9) * (t333 ^ 2 - t336 ^ 2);
t364 = pkin(4) * t336 + qJ(5) * t333;
t284 = qJD(4) * t364 - qJD(5) * t336 + qJD(3);
t338 = -pkin(2) - pkin(8);
t391 = qJD(4) * t338;
t378 = t336 * t391;
t317 = t330 * t378;
t337 = cos(qJ(2));
t329 = sin(pkin(6));
t401 = qJD(1) * t329;
t334 = sin(qJ(2));
t416 = t333 * t334;
t410 = t328 * t284 - (t328 * t337 + t330 * t416) * t401 + t317;
t435 = t330 * t284 - (-t328 * t416 + t330 * t337) * t401;
t385 = t337 * t401;
t363 = qJD(3) - t385;
t294 = qJD(2) * t338 + t363;
t331 = cos(pkin(6));
t400 = qJD(1) * t331;
t434 = t294 * t336 - t333 * t400;
t433 = MDP(13) * t333 + MDP(14) * t336;
t392 = qJD(4) * t336;
t431 = qJD(6) * t322;
t373 = qJD(4) * pkin(4) - qJD(5);
t259 = -t434 - t373;
t430 = t259 * MDP(18);
t304 = t328 * t335 + t330 * t332;
t348 = t304 * qJD(6);
t428 = pkin(9) * t336;
t427 = pkin(9) + qJ(5);
t426 = qJD(2) * pkin(2);
t398 = qJD(2) * t329;
t382 = t334 * t398;
t368 = qJD(1) * t382;
t321 = t336 * t400;
t393 = qJD(4) * t333;
t409 = qJD(4) * t321 + t294 * t393;
t246 = -t336 * t368 + t409;
t425 = t246 * t328;
t424 = t246 * t330;
t423 = t246 * t336;
t286 = t335 * t298;
t420 = t300 * t332;
t253 = t286 + t420;
t422 = t253 * t322;
t418 = t329 * t334;
t417 = t329 * t337;
t415 = t333 * t338;
t241 = t333 * t368 + (qJD(5) + t434) * qJD(4);
t256 = (t284 + t385) * qJD(2);
t226 = t330 * t241 + t328 * t256;
t375 = -t328 * t338 + pkin(5);
t387 = pkin(9) * t330 * t333;
t413 = (t336 * t375 + t387) * qJD(4) + t435;
t379 = t328 * t393;
t371 = pkin(9) * t379;
t412 = -t371 - t410;
t268 = t333 * t294 + t321;
t262 = qJD(4) * qJ(5) + t268;
t311 = pkin(4) * t333 - qJ(5) * t336 + qJ(3);
t386 = t334 * t401;
t275 = qJD(2) * t311 + t386;
t230 = t330 * t262 + t328 * t275;
t369 = t328 * t378;
t411 = -t369 + t435;
t306 = t364 * qJD(2);
t238 = t328 * t306 + t330 * t434;
t349 = t432 * t333;
t408 = qJD(2) * t349 + t429;
t347 = t304 * qJD(2);
t407 = t333 * t347 + t348;
t388 = qJD(2) * qJD(4);
t377 = t333 * t388;
t366 = t335 * t377;
t367 = t332 * t377;
t406 = -t328 * t366 - t330 * t367;
t274 = t328 * t311 + t330 * t415;
t399 = qJD(2) * qJ(3);
t395 = qJD(4) * t322;
t384 = t328 * t397;
t381 = t337 * t398;
t376 = MDP(23) * t396;
t374 = pkin(5) * t328 - t338;
t225 = -t241 * t328 + t330 * t256;
t345 = (pkin(5) * t336 + t387) * qJD(2);
t223 = qJD(4) * t345 + t225;
t224 = qJD(2) * t371 + t226;
t372 = t335 * t223 - t224 * t332;
t229 = -t262 * t328 + t330 * t275;
t237 = t330 * t306 - t328 * t434;
t370 = t336 * t386;
t307 = t386 + t399;
t365 = -t307 + t386;
t362 = t223 * t332 + t224 * t335;
t361 = -t225 * t328 + t226 * t330;
t227 = pkin(5) * t397 - pkin(9) * t300 + t229;
t228 = -pkin(9) * t298 + t230;
t220 = t227 * t335 - t228 * t332;
t221 = t227 * t332 + t228 * t335;
t360 = -t229 * t330 - t230 * t328;
t359 = -t229 * t328 + t230 * t330;
t296 = t330 * t311;
t251 = -t330 * t428 + t333 * t375 + t296;
t261 = -t328 * t428 + t274;
t358 = t251 * t335 - t261 * t332;
t357 = t251 * t332 + t261 * t335;
t291 = t331 * t336 - t333 * t417;
t263 = -t291 * t328 + t330 * t418;
t264 = t291 * t330 + t328 * t418;
t356 = t263 * t335 - t264 * t332;
t355 = t263 * t332 + t264 * t335;
t290 = t331 * t333 + t336 * t417;
t352 = t365 * qJD(2);
t316 = t427 * t330;
t351 = qJD(5) * t328 + qJD(6) * t316 + t237 + t345;
t315 = t427 * t328;
t350 = pkin(9) * t384 - qJD(5) * t330 + qJD(6) * t315 + t238;
t346 = t432 * qJD(2);
t343 = -qJD(6) * t286 + t328 * t367 - t330 * t366;
t233 = -qJD(6) * t354 + t406;
t342 = -qJ(5) * t392 + (t259 + t373) * t333;
t302 = (qJD(3) + t385) * qJD(2);
t339 = qJD(4) ^ 2;
t341 = qJD(2) * t363 - t338 * t339 + t302;
t232 = -qJD(6) * t420 + t343;
t340 = qJD(2) ^ 2;
t324 = -pkin(5) * t330 - pkin(4);
t305 = t363 - t426;
t297 = t374 * t336;
t287 = t374 * t393;
t283 = t432 * t336;
t282 = t304 * t336;
t273 = -t328 * t415 + t296;
t265 = -qJD(4) * t290 + t333 * t382;
t250 = -pkin(5) * t384 + t268;
t245 = -t332 * t333 * t389 - t335 * t379 + t336 * t429;
t244 = -qJD(4) * t349 - t336 * t348;
t243 = t265 * t330 + t328 * t381;
t242 = -t265 * t328 + t330 * t381;
t240 = pkin(5) * t298 + t259;
t235 = (-pkin(5) * t379 - t370) * qJD(2) + t409;
t1 = [-t265 * qJD(4) * MDP(14) + (-t242 * t300 - t243 * t298) * MDP(17) + (t225 * t263 + t226 * t264 + t229 * t242 + t230 * t243 + t246 * t290) * MDP(18) + ((-qJD(6) * t355 + t242 * t335 - t243 * t332) * t322 + t290 * t233) * MDP(24) + (-(qJD(6) * t356 + t242 * t332 + t243 * t335) * t322 + t290 * t232) * MDP(25) + (-MDP(13) * qJD(4) + t298 * MDP(15) + t300 * MDP(16) + t253 * MDP(24) - MDP(25) * t354 + t430) * (qJD(4) * t291 - t336 * t382) + ((MDP(15) * t242 - MDP(16) * t243) * t333 + ((t263 * MDP(15) - t264 * MDP(16) + MDP(24) * t356 - MDP(25) * t355) * t336 + ((t263 * t330 + t264 * t328) * MDP(17) + (-t328 * MDP(15) - t330 * MDP(16)) * t290) * t333) * qJD(4)) * qJD(2) + ((t307 * qJD(2) * MDP(7) + (-MDP(4) + MDP(6) + t433) * t340) * t337 + (t302 * MDP(7) + (-MDP(3) + MDP(5)) * t340 + ((MDP(13) * t336 - MDP(14) * t333) * qJD(4) + (t305 - t385) * MDP(7)) * qJD(2)) * t334) * t329; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t302 + qJD(3) * t307 + (-t307 * t337 + (-t305 - t426) * t334) * t401) * MDP(7) + 0.2e1 * t388 * t436 + (-t411 * t300 - t410 * t298 + ((t273 * t330 + t274 * t328) * qJD(2) - t360) * t393) * MDP(17) + (t225 * t273 + t226 * t274 + t411 * t229 + t410 * t230 + t259 * t370 - t338 * t423) * MDP(18) + (t232 * t283 - t244 * t354) * MDP(19) + (-t232 * t282 - t233 * t283 - t244 * t253 + t245 * t354) * MDP(20) + (t244 * t322 + (qJD(2) * t283 - t354) * t392) * MDP(21) + (-t245 * t322 + (-qJD(2) * t282 - t253) * t392) * MDP(22) + (t322 + t397) * MDP(23) * t392 + (-t287 * t253 + t297 * t233 + t235 * t282 + t240 * t245 + (t332 * t412 + t335 * t413) * t322 - t357 * t431) * MDP(24) + (t287 * t354 + t297 * t232 + t235 * t283 + t240 * t244 + (-t332 * t413 + t335 * t412) * t322 - t358 * t431) * MDP(25) + (-MDP(13) * t392 + MDP(14) * t393) * (t365 - t399) + (-t339 * MDP(10) + t341 * MDP(13) + (t225 + (-t259 * t328 + t298 * t338) * qJD(4) + (t369 + t411) * qJD(2)) * MDP(15) + (-t226 + (-t259 * t330 + t300 * t338) * qJD(4) + (t317 - t410) * qJD(2)) * MDP(16) + t391 * t430 + t232 * MDP(21) - t233 * MDP(22) + (-qJD(6) * t221 + t372) * MDP(24) + (-qJD(6) * t220 - t362) * MDP(25)) * t333 + (-0.2e1 * MDP(8) * t377 - t339 * MDP(11) + t341 * MDP(14) + (t298 * t386 + t425 + (qJD(2) * t273 + t229) * qJD(4)) * MDP(15) + (t300 * t386 + t424 + (-qJD(2) * t274 - t230) * qJD(4)) * MDP(16) + (-t225 * t330 - t226 * t328) * MDP(17) + (t253 * t386 + (qJD(2) * t358 + t220) * qJD(4)) * MDP(24) + (-t354 * t386 + (-qJD(2) * t357 - t221) * qJD(4)) * MDP(25)) * t336; -t340 * MDP(6) + MDP(7) * t352 + ((-t298 * t330 + t300 * t328) * t392 + (t298 * t328 + t300 * t330) * qJD(2)) * MDP(17) + (t360 * qJD(2) + t359 * t392 - t423) * MDP(18) + (-t322 * t346 + (-t304 * t395 - t233) * t336) * MDP(24) + (t322 * t347 + (-t395 * t432 - t232) * t336) * MDP(25) + t433 * (-t339 - t340) + ((-t330 * t340 + (t298 - t383) * qJD(4)) * MDP(15) + (t328 * t340 + (t300 - t380) * qJD(4)) * MDP(16) + (qJD(4) * t259 + t361) * MDP(18) + (-t322 * t429 + (-t304 * t396 + t253) * qJD(4)) * MDP(24) + (t322 * t348 + (-t336 * t346 - t354) * qJD(4)) * MDP(25)) * t333; (qJD(4) * t268 + t336 * t352 - t409) * MDP(13) - t365 * t397 * MDP(14) + (-t424 - t268 * t298 + (-t229 * t336 - t237 * t333 + t328 * t342) * qJD(2)) * MDP(15) + (t425 - t268 * t300 + (t230 * t336 + t238 * t333 + t330 * t342) * qJD(2)) * MDP(16) + (t237 * t300 + t238 * t298 + (-qJD(5) * t298 - t229 * t397 + t226) * t330 + (qJD(5) * t300 - t230 * t397 - t225) * t328) * MDP(17) + (-pkin(4) * t246 + qJ(5) * t361 + qJD(5) * t359 - t229 * t237 - t230 * t238 - t259 * t268) * MDP(18) + (t232 * t304 - t354 * t408) * MDP(19) + (t232 * t432 - t233 * t304 - t253 * t408 + t354 * t407) * MDP(20) + (t408 * t322 + (qJD(4) * t304 + t354) * t396) * MDP(21) + (-t407 * t322 + (qJD(4) * t432 + t253) * t396) * MDP(22) - t322 * t376 + (t324 * t233 - t235 * t432 - t250 * t253 + (t332 * t350 - t335 * t351) * t322 + t407 * t240 + ((-t315 * t335 - t316 * t332) * qJD(4) - t220) * t396) * MDP(24) + (t324 * t232 + t235 * t304 + t250 * t354 + (t332 * t351 + t335 * t350) * t322 + t408 * t240 + (-(-t315 * t332 + t316 * t335) * qJD(4) + t221) * t396) * MDP(25) + (t336 * t333 * MDP(8) - t436) * t340; (-t298 ^ 2 - t300 ^ 2) * MDP(17) + (t229 * t300 + t230 * t298 + t246) * MDP(18) + (t233 - t437) * MDP(24) + (t232 - t422) * MDP(25) + ((t300 - t394) * MDP(15) + (-t298 - t389) * MDP(16)) * t397; -t354 * t253 * MDP(19) + (-t253 ^ 2 + t354 ^ 2) * MDP(20) + (t343 + t422) * MDP(21) + (-t406 - t437) * MDP(22) + qJD(4) * t376 + (t221 * t322 + t240 * t354 + t372) * MDP(24) + (t220 * t322 + t240 * t253 - t362) * MDP(25) + (-MDP(21) * t420 + MDP(22) * t354 - t221 * MDP(24) - t220 * MDP(25)) * qJD(6);];
tauc  = t1;
