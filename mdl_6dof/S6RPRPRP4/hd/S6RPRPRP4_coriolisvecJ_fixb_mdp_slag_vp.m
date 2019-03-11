% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:06
% EndTime: 2019-03-09 03:13:14
% DurationCPUTime: 3.18s
% Computational Cost: add. (2313->363), mult. (4879->474), div. (0->0), fcn. (2662->6), ass. (0->160)
t309 = sin(pkin(9)) * pkin(1) + pkin(7);
t419 = pkin(4) + t309;
t299 = t309 * qJD(1);
t326 = sin(qJ(3));
t328 = cos(qJ(3));
t273 = -t328 * qJD(2) + t326 * t299;
t425 = qJD(4) + t273;
t398 = qJD(1) * t326;
t308 = qJD(5) + t398;
t424 = t308 ^ 2;
t321 = t326 ^ 2;
t322 = t328 ^ 2;
t423 = (t321 - t322) * MDP(6);
t329 = -pkin(3) - pkin(8);
t310 = -cos(pkin(9)) * pkin(1) - pkin(2);
t343 = -qJ(4) * t326 + t310;
t270 = t329 * t328 + t343;
t283 = t419 * t326;
t325 = sin(qJ(5));
t327 = cos(qJ(5));
t422 = t327 * t270 + t325 * t283;
t383 = qJD(1) * qJD(3);
t369 = t326 * t383;
t396 = qJD(3) * t325;
t397 = qJD(1) * t328;
t292 = t327 * t397 + t396;
t391 = qJD(5) * t292;
t256 = -t325 * t369 + t391;
t274 = t326 * qJD(2) + t328 * t299;
t320 = qJD(3) * qJ(4);
t263 = -t320 - t274;
t260 = -qJD(3) * pkin(3) + t425;
t385 = pkin(4) * t398 + t425;
t381 = MDP(21) + MDP(23);
t380 = MDP(22) - MDP(25);
t395 = qJD(3) * t326;
t314 = pkin(3) * t395;
t355 = pkin(8) * t326 - qJ(4) * t328;
t392 = qJD(4) * t326;
t333 = t355 * qJD(3) - t392;
t261 = t314 + t333;
t393 = qJD(3) * t328;
t280 = t419 * t393;
t420 = -qJD(5) * t422 - t261 * t325 + t280 * t327;
t319 = qJD(3) * qJD(4);
t382 = qJD(2) * qJD(3);
t400 = -t299 * t395 + t328 * t382;
t255 = -t319 - t400;
t237 = -pkin(4) * t369 - t255;
t368 = t327 * t383;
t370 = t325 * t397;
t389 = qJD(5) * t327;
t257 = qJD(3) * t389 - qJD(5) * t370 - t326 * t368;
t394 = qJD(3) * t327;
t294 = -t370 + t394;
t221 = pkin(5) * t257 + qJ(6) * t256 - qJD(6) * t294 + t237;
t418 = t221 * t325;
t417 = t221 * t327;
t416 = t237 * t325;
t415 = t237 * t327;
t414 = t256 * t327;
t413 = t257 * t325;
t412 = t292 * t308;
t411 = t292 * t327;
t410 = t294 * t325;
t409 = t308 * t329;
t248 = t326 * t257;
t408 = t326 * t327;
t330 = qJD(3) ^ 2;
t407 = t326 * t330;
t406 = t328 * t330;
t357 = pkin(5) * t327 + qJ(6) * t325;
t342 = -pkin(4) - t357;
t405 = -t357 * qJD(5) + qJD(6) * t327 + t342 * t398 - t425;
t259 = pkin(4) * t397 + t274;
t315 = pkin(3) * t398;
t275 = t355 * qJD(1) + t315;
t403 = t325 * t259 + t327 * t275;
t402 = -t326 * t256 + t294 * t393;
t266 = t299 * t393 + t326 * t382;
t390 = qJD(5) * t325;
t372 = t308 * t390;
t374 = t326 * t394;
t401 = t308 * t374 + t328 * t372;
t284 = t419 * t328;
t282 = -pkin(3) * t328 + t343;
t264 = qJD(1) * t282;
t300 = qJD(1) * t310;
t388 = qJD(6) * t308;
t250 = t320 + t259;
t387 = t250 * qJD(5);
t386 = t294 * MDP(24);
t243 = t329 * qJD(3) + t385;
t252 = t270 * qJD(1);
t225 = t243 * t327 - t252 * t325;
t384 = qJD(6) - t225;
t379 = t325 * t409;
t378 = t327 * t409;
t331 = qJD(1) ^ 2;
t377 = t326 * t328 * t331;
t307 = pkin(3) * t369;
t246 = t333 * qJD(1) + t307;
t367 = t328 * t383;
t247 = pkin(4) * t367 + t266;
t376 = -t243 * t389 - t327 * t246 - t325 * t247;
t375 = t325 * t395;
t373 = t329 * t393;
t371 = t328 * t389;
t366 = MDP(20) * t397;
t364 = pkin(5) * t367;
t363 = -t243 * t390 - t325 * t246 + t327 * t247 - t252 * t389;
t362 = qJ(6) * t367;
t361 = t322 * t325 * t383;
t303 = t327 * t367;
t219 = -t363 - t364;
t226 = t243 * t325 + t252 * t327;
t224 = qJ(6) * t308 + t226;
t360 = -t224 * t398 + t219;
t359 = qJD(3) * t273 + t400;
t358 = qJD(3) * t274 - t266;
t356 = -pkin(5) * t325 + qJ(6) * t327;
t337 = t252 * t390 + t376;
t218 = -t337 + t362 + t388;
t354 = t218 * t325 - t219 * t327;
t223 = -pkin(5) * t308 + t384;
t353 = -t223 * t327 + t224 * t325;
t352 = t223 * t325 + t224 * t327;
t351 = t259 * t327 - t275 * t325;
t348 = -t270 * t325 + t283 * t327;
t347 = -0.2e1 * qJD(3) * t264;
t346 = 0.2e1 * qJD(3) * t300;
t345 = t308 * t294;
t341 = t308 * t375 - t361;
t339 = -qJ(4) * t393 - t392;
t265 = t339 * qJD(1) + t307;
t281 = t314 + t339;
t340 = qJD(1) * t281 + t309 * t330 + t265;
t338 = t226 * t308 + t363;
t336 = t327 * t261 - t270 * t390 + t325 * t280 + t283 * t389;
t335 = t292 * t393 - t322 * t368 + t248;
t334 = t225 * t308 + t337;
t332 = -t255 * t328 + t266 * t326 + (t260 * t328 + t263 * t326) * qJD(3);
t298 = qJ(4) - t356;
t297 = t329 * t303;
t296 = -qJ(4) * t397 + t315;
t279 = t419 * t395;
t269 = t292 * t375;
t253 = t264 * t398;
t245 = t357 * t328 + t284;
t244 = pkin(5) * t294 + qJ(6) * t292;
t233 = t412 - t256;
t232 = -pkin(5) * t326 - t348;
t231 = qJ(6) * t326 + t422;
t230 = -pkin(5) * t397 - t351;
t229 = qJ(6) * t397 + t403;
t228 = pkin(5) * t292 - qJ(6) * t294 + t250;
t227 = (t356 * qJD(5) + qJD(6) * t325) * t328 + (-t309 + t342) * t395;
t222 = -pkin(5) * t393 - t420;
t220 = qJ(6) * t393 + qJD(6) * t326 + t336;
t1 = [0.2e1 * t326 * MDP(5) * t367 - 0.2e1 * t383 * t423 + MDP(7) * t406 - MDP(8) * t407 + (-t309 * t406 + t326 * t346) * MDP(10) + (t309 * t407 + t328 * t346) * MDP(11) + t332 * MDP(12) + (t326 * t347 + t340 * t328) * MDP(13) + (-t340 * t326 + t328 * t347) * MDP(14) + (t264 * t281 + t265 * t282 + t332 * t309) * MDP(15) + (t256 * t325 * t328 + (-t371 + t375) * t294) * MDP(16) + (t294 * t374 - t269 + (t414 + t413 + (t410 + t411) * qJD(5)) * t328) * MDP(17) + (-t308 * t371 + t341 + t402) * MDP(18) + (-t248 + (-qJD(1) * t322 * t327 - t292 * t328) * qJD(3) + t401) * MDP(19) + (t308 + t398) * MDP(20) * t393 + (t420 * t308 - t279 * t292 + t284 * t257 + (-t250 * t394 + t363) * t326 + (-t325 * t387 + t415 + (t348 * qJD(1) + t225) * qJD(3)) * t328) * MDP(21) + (-t336 * t308 - t279 * t294 - t284 * t256 + ((qJD(3) * t250 + qJD(5) * t252) * t325 + t376) * t326 + (-t327 * t387 - t416 + (-qJD(1) * t422 - t226) * qJD(3)) * t328) * MDP(22) + (-t222 * t308 + t227 * t292 + t245 * t257 + (-t228 * t394 - t219) * t326 + (-t228 * t390 + t417 + (-qJD(1) * t232 - t223) * qJD(3)) * t328) * MDP(23) + (-t220 * t292 + t222 * t294 - t231 * t257 - t232 * t256 + t352 * t395 + (t353 * qJD(5) - t218 * t327 - t219 * t325) * t328) * MDP(24) + (t220 * t308 - t227 * t294 + t245 * t256 + (-t228 * t396 + t218) * t326 + (t228 * t389 + t418 + (qJD(1) * t231 + t224) * qJD(3)) * t328) * MDP(25) + (t218 * t231 + t219 * t232 + t220 * t224 + t221 * t245 + t222 * t223 + t227 * t228) * MDP(26); (t335 + t401) * MDP(21) + (t361 + t402) * MDP(22) + t335 * MDP(23) - t269 * MDP(24) + t341 * MDP(25) + (-t255 * MDP(15) + t256 * MDP(25) + t221 * MDP(26) + (-MDP(10) + MDP(13)) * t330 + (t260 * MDP(15) - t327 * t386 + t353 * MDP(26) + (-t325 * MDP(22) + t327 * MDP(23)) * t308) * qJD(3)) * t326 + ((-qJD(3) * t263 - t266) * MDP(15) + (t413 - t414) * MDP(24) - qJD(3) * t294 * MDP(25) + (qJD(3) * t228 - t354) * MDP(26) + (-MDP(11) + MDP(14)) * t330 + ((-t410 + t411) * MDP(24) - t352 * MDP(26) + (t325 * MDP(23) + t380 * t327) * t308) * qJD(5)) * t328; -MDP(5) * t377 + t331 * t423 + (-t300 * t398 + t358) * MDP(10) + (-t300 * t397 - t359) * MDP(11) + (-t296 * t397 + t253 - t358) * MDP(13) + (0.2e1 * t319 + (t264 * t328 + t296 * t326) * qJD(1) + t359) * MDP(14) + (-pkin(3) * t266 - qJ(4) * t255 - t260 * t274 - t263 * t425 - t264 * t296) * MDP(15) + (-t325 * t345 - t414) * MDP(16) + ((-t257 - t345) * t327 + (t256 + t412) * t325) * MDP(17) + (-t372 + t303 + (-t308 * t325 * t326 - t294 * t328) * qJD(1)) * MDP(18) + (-t308 * t389 + (-t308 * t408 + (t292 - t396) * t328) * qJD(1)) * MDP(19) - t308 * t366 + (t297 + qJ(4) * t257 + t416 - t351 * t308 + t385 * t292 + (t250 * t327 - t379) * qJD(5) + (-t225 * t328 + t250 * t408) * qJD(1)) * MDP(21) + (-qJ(4) * t256 + t415 + t403 * t308 + t385 * t294 + (-t250 * t325 - t378) * qJD(5) + (t226 * t328 + (-t250 * t326 - t373) * t325) * qJD(1)) * MDP(22) + (t418 + t230 * t308 + t257 * t298 + t297 - t405 * t292 + (t228 * t327 - t379) * qJD(5) + (t223 * t328 + t228 * t408) * qJD(1)) * MDP(23) + (t229 * t292 - t230 * t294 + (t256 * t329 + (-t292 * t329 - t224) * qJD(5) + t360) * t327 + (-t223 * t398 - t257 * t329 - t218 + (t294 * t329 - t223) * qJD(5)) * t325) * MDP(24) + (-t417 - t229 * t308 + t256 * t298 + t405 * t294 + (t228 * t325 + t378) * qJD(5) + (-t224 * t328 + (t228 * t326 + t373) * t325) * qJD(1)) * MDP(25) + (t221 * t298 - t223 * t230 - t224 * t229 - t405 * t228 + (t352 * qJD(5) + t354) * t329) * MDP(26); MDP(13) * t377 + (-t321 * t331 - t330) * MDP(14) + (t253 + t266) * MDP(15) + t381 * t303 + (t263 * MDP(15) - t228 * MDP(26) - t381 * t292 - t380 * t294) * qJD(3) + ((-t292 * t398 + t256 - t391) * MDP(24) + (qJD(5) * t224 - t360) * MDP(26) - t380 * t424) * t327 + ((qJD(5) * t294 - t257) * MDP(24) + (qJD(5) * t223 + t218) * MDP(26) + ((t223 * MDP(26) + t386) * t326 - t380 * t393) * qJD(1) - t381 * t424) * t325; t233 * MDP(18) - t257 * MDP(19) + qJD(3) * t366 + t338 * MDP(21) + t334 * MDP(22) + (t338 + 0.2e1 * t364) * MDP(23) + (pkin(5) * t256 - qJ(6) * t257) * MDP(24) + (-t334 + 0.2e1 * t362 + 0.2e1 * t388) * MDP(25) + (-pkin(5) * t219 + qJ(6) * t218 - t223 * t226 + t384 * t224 - t228 * t244) * MDP(26) + (t308 * MDP(19) - t250 * MDP(21) - t228 * MDP(23) + (t224 - t226) * MDP(24) + t244 * MDP(25) + MDP(17) * t294) * t294 + (t294 * MDP(16) + t250 * MDP(22) - t244 * MDP(23) + (t223 - t384) * MDP(24) - t228 * MDP(25) - MDP(17) * t292) * t292; (t292 * t294 - t367) * MDP(23) + t233 * MDP(24) + (-t294 ^ 2 - t424) * MDP(25) + (-t224 * t308 + t228 * t294 + t219) * MDP(26);];
tauc  = t1;
