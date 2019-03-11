% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:14
% EndTime: 2019-03-09 03:26:20
% DurationCPUTime: 3.19s
% Computational Cost: add. (3938->380), mult. (8348->503), div. (0->0), fcn. (5516->6), ass. (0->148)
t334 = sin(qJ(3));
t336 = cos(qJ(3));
t421 = sin(pkin(9));
t422 = cos(pkin(9));
t341 = t334 * t422 + t336 * t421;
t296 = t341 * qJD(1);
t430 = qJD(5) + t296;
t429 = MDP(8) * (t334 ^ 2 - t336 ^ 2);
t428 = MDP(6) * qJ(2) + MDP(5);
t303 = -t334 * t421 + t336 * t422;
t401 = pkin(3) * t334 + qJ(2);
t268 = pkin(4) * t341 - pkin(8) * t303 + t401;
t337 = -pkin(1) - pkin(7);
t400 = qJ(4) - t337;
t307 = t400 * t334;
t308 = t400 * t336;
t272 = -t307 * t422 - t308 * t421;
t333 = sin(qJ(5));
t335 = cos(qJ(5));
t394 = t268 * t333 + t272 * t335;
t375 = MDP(21) + MDP(23);
t427 = MDP(22) - MDP(25);
t312 = qJD(1) * t337 + qJD(2);
t388 = qJD(1) * t334;
t292 = -qJ(4) * t388 + t312 * t334;
t285 = t421 * t292;
t387 = qJD(1) * t336;
t293 = -qJ(4) * t387 + t312 * t336;
t288 = qJD(3) * pkin(3) + t293;
t253 = t288 * t422 - t285;
t247 = -qJD(3) * pkin(4) - t253;
t299 = t303 * qJD(1);
t380 = t335 * qJD(3);
t275 = t299 * t333 - t380;
t277 = qJD(3) * t333 + t299 * t335;
t228 = pkin(5) * t275 - qJ(6) * t277 + t247;
t365 = qJD(3) * t421;
t358 = qJD(1) * t365;
t366 = qJD(3) * t422;
t359 = qJD(1) * t366;
t289 = t334 * t358 - t336 * t359;
t321 = pkin(3) * t421 + pkin(8);
t404 = t321 * t289;
t426 = t228 * t430 + t404;
t369 = t422 * t292;
t254 = t288 * t421 + t369;
t248 = qJD(3) * pkin(8) + t254;
t306 = pkin(3) * t388 + qJD(1) * qJ(2) + qJD(4);
t255 = pkin(4) * t296 - pkin(8) * t299 + t306;
t226 = -t248 * t333 + t255 * t335;
t377 = qJD(6) - t226;
t222 = -pkin(5) * t430 + t377;
t227 = t248 * t335 + t255 * t333;
t223 = qJ(6) * t430 + t227;
t354 = t222 * t335 - t223 * t333;
t407 = t277 * t335;
t411 = t275 * t333;
t425 = -(t407 + t411) * MDP(24) - t354 * MDP(26);
t424 = t277 ^ 2;
t423 = pkin(5) * t289;
t419 = qJ(6) * t289;
t384 = qJD(4) * t334;
t385 = qJD(3) * t336;
t274 = t312 * t385 + (-qJ(4) * t385 - t384) * qJD(1);
t383 = qJD(4) * t336;
t386 = qJD(3) * t334;
t340 = -t312 * t386 + (qJ(4) * t386 - t383) * qJD(1);
t237 = t274 * t421 - t340 * t422;
t382 = qJD(5) * t333;
t392 = t334 * t359 + t336 * t358;
t249 = -qJD(5) * t380 + t299 * t382 + t335 * t392;
t364 = t333 * t392;
t250 = qJD(5) * t277 - t364;
t219 = pkin(5) * t250 + qJ(6) * t249 - qJD(6) * t277 + t237;
t418 = t219 * t333;
t417 = t227 * t430;
t416 = t237 * t303;
t415 = t249 * t333;
t414 = t250 * t335;
t298 = -t334 * t366 - t336 * t365;
t413 = t253 * t298;
t412 = t275 * t296;
t410 = t275 * t335;
t409 = t277 * t275;
t408 = t277 * t333;
t406 = t430 * t296;
t405 = t303 * t335;
t281 = t333 * t289;
t283 = t335 * t289;
t339 = qJD(1) ^ 2;
t403 = t336 * t339;
t338 = qJD(3) ^ 2;
t402 = t337 * t338;
t260 = t293 * t421 + t369;
t356 = pkin(5) * t333 - qJ(6) * t335;
t399 = qJD(6) * t333 - t356 * t430 + t260;
t381 = qJD(5) * t335;
t397 = -t250 * t333 - t275 * t381;
t261 = t293 * t422 - t285;
t263 = pkin(3) * t387 + pkin(4) * t299 + pkin(8) * t296;
t396 = t261 * t335 + t263 * t333;
t395 = -t333 * t406 - t283;
t393 = t381 * t430 - t281;
t329 = qJD(1) * qJD(2);
t376 = qJD(1) * qJD(3);
t370 = t336 * t376;
t391 = pkin(3) * t370 + t329;
t379 = pkin(3) * t385 + qJD(2);
t373 = 0.2e1 * qJD(1);
t371 = t321 * t382;
t238 = t274 * t422 + t340 * t421;
t244 = -pkin(4) * t289 + pkin(8) * t392 + t391;
t344 = t238 * t335 + t244 * t333 - t248 * t382 + t255 * t381;
t216 = qJD(6) * t430 + t344 - t419;
t363 = t222 * t296 + t216;
t360 = t238 * t333 - t244 * t335 + t248 * t381 + t255 * t382;
t217 = t360 + t423;
t362 = -t223 * t296 + t217;
t361 = t277 * t430;
t322 = -pkin(3) * t422 - pkin(4);
t357 = pkin(5) * t335 + qJ(6) * t333;
t290 = t386 * t400 - t383;
t291 = -qJD(3) * t308 - t384;
t257 = -t290 * t422 + t291 * t421;
t271 = -t307 * t421 + t308 * t422;
t355 = -t216 * t335 - t217 * t333;
t353 = -t222 * t333 - t223 * t335;
t352 = t272 * t289 + t416;
t347 = t298 * t333 + t303 * t381;
t346 = -t298 * t335 + t303 * t382;
t345 = t228 * t277 + t360;
t258 = t290 * t421 + t291 * t422;
t297 = t334 * t365 - t336 * t366;
t262 = -pkin(4) * t297 - pkin(8) * t298 + t379;
t343 = t258 * t335 + t262 * t333 + t268 * t381 - t272 * t382;
t342 = t247 * t430 + t404;
t301 = -t357 + t322;
t242 = pkin(5) * t277 + qJ(6) * t275;
t233 = t303 * t356 + t271;
t231 = -pkin(5) * t341 - t268 * t335 + t272 * t333;
t230 = qJ(6) * t341 + t394;
t229 = t275 * t430 - t249;
t225 = -pkin(5) * t299 + t261 * t333 - t263 * t335;
t224 = qJ(6) * t299 + t396;
t221 = t356 * t298 + (qJD(5) * t357 - qJD(6) * t335) * t303 + t257;
t220 = pkin(5) * t297 + qJD(5) * t394 + t258 * t333 - t262 * t335;
t218 = -qJ(6) * t297 + qJD(6) * t341 + t343;
t1 = [0.2e1 * t376 * t429 - t338 * t336 * MDP(10) + qJ(2) * t385 * t373 * MDP(12) + (-t336 * t402 + (-qJ(2) * t386 + qJD(2) * t336) * t373) * MDP(13) + (-t238 * t341 + t254 * t297 + t257 * t299 - t258 * t296 - t271 * t392 + t352 - t413) * MDP(14) + (t237 * t271 + t238 * t272 - t253 * t257 + t254 * t258 + t306 * t379 + t391 * t401) * MDP(15) + (-t249 * t405 - t277 * t346) * MDP(16) + ((-t408 - t410) * t298 + (t415 - t414 + (-t407 + t411) * qJD(5)) * t303) * MDP(17) + (-t249 * t341 - t277 * t297 - t283 * t303 - t346 * t430) * MDP(18) + (-t250 * t341 + t275 * t297 + t281 * t303 - t347 * t430) * MDP(19) + (-t289 * t341 - t297 * t430) * MDP(20) + (-t360 * t341 - t226 * t297 + t257 * t275 + t271 * t250 + ((-qJD(5) * t272 + t262) * t430 - t268 * t289 + t247 * qJD(5) * t303) * t335 + ((-qJD(5) * t268 - t258) * t430 + t247 * t298 + t352) * t333) * MDP(21) + (t227 * t297 + t237 * t405 - t247 * t346 - t271 * t249 + t257 * t277 + t289 * t394 - t341 * t344 - t343 * t430) * MDP(22) + (-t217 * t341 - t220 * t430 + t221 * t275 + t222 * t297 + t228 * t347 + t231 * t289 + t233 * t250 + t303 * t418) * MDP(23) + (-t218 * t275 + t220 * t277 - t230 * t250 - t231 * t249 + t354 * t298 + (qJD(5) * t353 - t216 * t333 + t217 * t335) * t303) * MDP(24) + (t216 * t341 + t218 * t430 - t219 * t405 - t221 * t277 - t223 * t297 + t228 * t346 - t230 * t289 + t233 * t249) * MDP(25) + (t216 * t230 + t217 * t231 + t218 * t223 + t219 * t233 + t220 * t222 + t221 * t228) * MDP(26) + 0.2e1 * t428 * t329 + (-0.2e1 * MDP(7) * t370 - t338 * MDP(9) + (qJD(2) * t373 - t402) * MDP(12)) * t334; (-t298 * t299 + t303 * t392) * MDP(14) + (t413 - t416) * MDP(15) + (-t219 * t303 - t228 * t298) * MDP(26) - t428 * t339 + (MDP(14) * t296 - t254 * MDP(15) + (-t408 + t410) * MDP(24) + t353 * MDP(26)) * t297 + (-t306 * MDP(15) - t425) * qJD(1) + (t375 * (-qJD(1) * t335 + t297 * t333) + t427 * (qJD(1) * t333 + t297 * t335)) * t430 - (-t238 * MDP(15) + (t414 + t415) * MDP(24) + t355 * MDP(26) + (-t333 * t375 - t335 * t427 - MDP(14)) * t289 + ((-t333 * t427 + t335 * t375) * t430 + t425) * qJD(5)) * t341 + t375 * (-t303 * t250 - t298 * t275) + (MDP(12) * t334 + MDP(13) * t336) * (-t338 - t339) + t427 * (t249 * t303 - t277 * t298); t334 * MDP(7) * t403 - t339 * t429 + ((t254 - t260) * t299 - (-t261 + t253) * t296 + (t289 * t421 + t392 * t422) * pkin(3)) * MDP(14) + (t253 * t260 - t254 * t261 + (-t237 * t422 + t238 * t421 - t306 * t387) * pkin(3)) * MDP(15) + (t335 * t361 - t415) * MDP(16) + ((-t249 - t412) * t335 - t430 * t408 + t397) * MDP(17) + (-t277 * t299 + t335 * t406 + t393) * MDP(18) + (t275 * t299 - t382 * t430 + t395) * MDP(19) - t430 * t299 * MDP(20) + (-t226 * t299 + t322 * t250 - t260 * t275 + (-t237 + (-qJD(5) * t321 - t263) * t430) * t335 + (t261 * t430 + t342) * t333) * MDP(21) + (t227 * t299 + t237 * t333 - t322 * t249 - t260 * t277 + (t371 + t396) * t430 + t342 * t335) * MDP(22) + (-t219 * t335 + t222 * t299 + t250 * t301 + (-t321 * t381 + t225) * t430 - t399 * t275 + t426 * t333) * MDP(23) + (t224 * t275 - t225 * t277 + (-t250 * t321 + (t277 * t321 + t222) * qJD(5) + t363) * t335 + (-t249 * t321 + (t275 * t321 - t223) * qJD(5) + t362) * t333) * MDP(24) + (-t418 - t223 * t299 + t249 * t301 + (-t224 - t371) * t430 + t399 * t277 - t426 * t335) * MDP(25) + (t219 * t301 - t222 * t225 - t223 * t224 - t399 * t228 + (qJD(5) * t354 - t355) * t321) * MDP(26) + (MDP(13) * t334 * t339 - MDP(12) * t403) * qJ(2); -t296 ^ 2 * MDP(14) + (t254 * t296 + t391) * MDP(15) + t395 * MDP(21) + t397 * MDP(24) + t393 * MDP(25) + (-MDP(14) * t299 + MDP(15) * t253 - MDP(26) * t228 - t275 * t375 - t277 * t427) * t299 + (-t289 * MDP(23) + (t249 - t412) * MDP(24) + (qJD(5) * t223 - t362) * MDP(26) + (-MDP(22) * t430 + t296 * MDP(25)) * t430) * t335 + (t289 * MDP(22) + (qJD(5) * t222 + t363) * MDP(26) + MDP(24) * t361 + (-qJD(5) * MDP(21) - MDP(23) * t430) * t430) * t333; MDP(16) * t409 + (-t275 ^ 2 + t424) * MDP(17) + t229 * MDP(18) + (t364 + (-qJD(5) + t430) * t277) * MDP(19) - t289 * MDP(20) + (-t247 * t277 - t360 + t417) * MDP(21) + (t226 * t430 + t247 * t275 - t344) * MDP(22) + (-t242 * t275 - t345 + t417 - 0.2e1 * t423) * MDP(23) + (pkin(5) * t249 - qJ(6) * t250 + (t223 - t227) * t277 + (t222 - t377) * t275) * MDP(24) + (-0.2e1 * t419 - t228 * t275 + t242 * t277 + (0.2e1 * qJD(6) - t226) * t430 + t344) * MDP(25) + (-pkin(5) * t217 + qJ(6) * t216 - t222 * t227 + t223 * t377 - t228 * t242) * MDP(26); (t289 + t409) * MDP(23) + t229 * MDP(24) + (-t430 ^ 2 - t424) * MDP(25) + (-t223 * t430 + t345 + t423) * MDP(26);];
tauc  = t1;
