% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:07
% EndTime: 2019-03-09 08:32:12
% DurationCPUTime: 3.42s
% Computational Cost: add. (3319->363), mult. (8201->467), div. (0->0), fcn. (5667->6), ass. (0->168)
t341 = sin(pkin(9));
t342 = cos(pkin(9));
t346 = cos(qJ(2));
t404 = t342 * t346;
t378 = qJD(1) * t404;
t344 = sin(qJ(2));
t385 = t344 * qJD(1);
t307 = t341 * t385 - t378;
t343 = sin(qJ(5));
t345 = cos(qJ(5));
t285 = qJD(2) * t345 + t307 * t343;
t320 = t341 * t346 + t342 * t344;
t310 = t320 * qJD(1);
t429 = qJD(5) + t310;
t422 = t345 * t429;
t434 = t285 * t422;
t379 = -pkin(2) * t346 - pkin(1);
t364 = t379 * qJD(1);
t324 = qJD(3) + t364;
t351 = -qJ(4) * t310 + t324;
t416 = pkin(3) + pkin(8);
t241 = t307 * t416 + t351;
t413 = -qJ(3) - pkin(7);
t326 = t413 * t346;
t323 = qJD(1) * t326;
t313 = t341 * t323;
t325 = t413 * t344;
t322 = qJD(1) * t325;
t318 = qJD(2) * pkin(2) + t322;
t276 = t318 * t342 + t313;
t363 = qJD(4) - t276;
t414 = t310 * pkin(4);
t246 = -qJD(2) * t416 + t363 + t414;
t227 = t241 * t345 + t246 * t343;
t389 = qJD(5) * t343;
t309 = t320 * qJD(2);
t352 = qJD(1) * t309;
t388 = qJD(5) * t345;
t393 = t307 * t388 + t343 * t352;
t262 = qJD(2) * t389 - t393;
t381 = qJD(1) * qJD(2);
t377 = t344 * t381;
t327 = t341 * t377;
t376 = t346 * t381;
t300 = t342 * t376 - t327;
t331 = pkin(2) * t377;
t373 = -t300 * qJ(4) + t331;
t386 = t310 * qJD(4);
t231 = t416 * t381 * t320 + t373 - t386;
t375 = qJD(2) * t413;
t304 = qJD(3) * t346 + t344 * t375;
t295 = t304 * qJD(1);
t305 = -qJD(3) * t344 + t346 * t375;
t296 = t305 * qJD(1);
t258 = t295 * t341 - t342 * t296;
t238 = pkin(4) * t300 + t258;
t370 = -t231 * t343 + t345 * t238;
t216 = pkin(5) * t300 + qJ(6) * t262 - qJD(5) * t227 - qJD(6) * t285 + t370;
t283 = qJD(2) * t343 - t345 * t307;
t222 = -qJ(6) * t283 + t227;
t433 = t222 * t429 + t216;
t293 = t345 * t352;
t390 = qJD(5) * t285;
t263 = -t293 + t390;
t380 = -t345 * t231 - t343 * t238 - t246 * t388;
t356 = -t241 * t389 - t380;
t217 = -qJ(6) * t263 - qJD(6) * t283 + t356;
t226 = -t241 * t343 + t345 * t246;
t221 = -qJ(6) * t285 + t226;
t220 = pkin(5) * t429 + t221;
t432 = -t220 * t429 + t217;
t365 = t429 * t283;
t431 = t262 - t365;
t294 = t345 * t300;
t368 = t343 * t429;
t430 = t368 * t429 - t294;
t428 = -0.2e1 * t381;
t427 = MDP(4) * t344;
t426 = MDP(5) * (t344 ^ 2 - t346 ^ 2);
t319 = t341 * t344 - t404;
t361 = -qJ(4) * t320 + t379;
t257 = t319 * t416 + t361;
t280 = -t342 * t325 - t326 * t341;
t271 = pkin(4) * t320 + t280;
t396 = t345 * t257 + t343 * t271;
t353 = qJD(2) * t310;
t279 = t342 * t322 + t313;
t384 = -qJD(4) + t279;
t405 = t342 * t323;
t277 = t341 * t318 - t405;
t274 = -qJD(2) * qJ(4) - t277;
t415 = pkin(4) * t307;
t251 = -t274 - t415;
t334 = -pkin(2) * t342 - pkin(3);
t330 = -pkin(8) + t334;
t419 = t251 * t429 + t330 * t300;
t418 = t285 ^ 2;
t306 = t310 ^ 2;
t412 = qJ(6) * t345;
t411 = t258 * t280;
t410 = t262 * t345;
t409 = t283 * t307;
t408 = t285 * t307;
t407 = t319 * t343;
t403 = t343 * t300;
t347 = qJD(2) ^ 2;
t402 = t344 * t347;
t401 = t346 * t347;
t348 = qJD(1) ^ 2;
t400 = t346 * t348;
t399 = qJ(6) - t330;
t398 = t220 - t221;
t372 = pkin(2) * t385 + qJ(4) * t307;
t247 = t310 * t416 + t372;
t278 = t322 * t341 - t405;
t260 = t278 - t415;
t397 = t345 * t247 + t343 * t260;
t259 = t342 * t295 + t341 * t296;
t255 = t345 * t260;
t387 = qJD(6) * t345;
t395 = t399 * t389 - t387 + pkin(5) * t307 - t255 - (-qJ(6) * t310 - t247) * t343;
t317 = t399 * t345;
t394 = -qJD(5) * t317 - qJD(6) * t343 - t310 * t412 - t397;
t391 = qJD(2) * t344;
t383 = t414 - t384;
t336 = pkin(2) * t391;
t337 = qJD(2) * qJD(4);
t253 = -t337 - t259;
t332 = pkin(2) * t341 + qJ(4);
t374 = pkin(1) * t428;
t371 = -qJ(6) * t319 - t257;
t269 = t304 * t341 - t342 * t305;
t369 = t429 ^ 2;
t270 = t304 * t342 + t305 * t341;
t281 = t325 * t341 - t326 * t342;
t265 = pkin(3) * t307 + t351;
t360 = t265 * t310 + t258;
t359 = t309 * t343 + t319 * t388;
t358 = -t309 * t345 + t319 * t389;
t312 = qJD(2) * t404 - t341 * t391;
t357 = -qJ(4) * t312 - qJD(4) * t320 + t336;
t234 = t309 * t416 + t357;
t248 = pkin(4) * t312 + t269;
t355 = t345 * t234 + t343 * t248 - t257 * t389 + t271 * t388;
t249 = -pkin(4) * t309 + t270;
t354 = -t422 * t429 - t403;
t350 = t258 * t320 + t269 * t310 - t270 * t307 + t280 * t300;
t235 = -pkin(4) * t352 - t253;
t349 = pkin(3) * t353 + t373;
t228 = t263 * pkin(5) + t235;
t316 = t399 * t343;
t302 = qJD(2) * t307;
t282 = t283 ^ 2;
t275 = pkin(3) * t319 + t361;
t273 = -qJD(2) * pkin(3) + t363;
t272 = -pkin(4) * t319 + t281;
t268 = pkin(3) * t310 + t372;
t267 = t345 * t271;
t256 = t345 * t263;
t250 = pkin(3) * t309 + t357;
t245 = t345 * t248;
t239 = t349 - t386;
t232 = pkin(5) * t283 + qJD(6) + t251;
t229 = t319 * t412 + t396;
t225 = pkin(5) * t320 + t343 * t371 + t267;
t219 = -qJ(6) * t358 + t319 * t387 + t355;
t218 = pkin(5) * t312 + t245 + t371 * t388 + (-qJ(6) * t309 - qJD(5) * t271 - qJD(6) * t319 - t234) * t343;
t1 = [0.2e1 * t376 * t427 + t426 * t428 + MDP(6) * t401 - MDP(7) * t402 + (-pkin(7) * t401 + t344 * t374) * MDP(9) + (pkin(7) * t402 + t346 * t374) * MDP(10) + (-t259 * t319 - t276 * t312 - t277 * t309 - t281 * t353 + t350) * MDP(11) + (t411 + t259 * t281 - t276 * t269 + t277 * t270 + (t324 + t364) * t336) * MDP(12) + (t253 * t319 + t273 * t312 + t274 * t309 - t281 * t352 + t350) * MDP(13) + (-t239 * t319 - t250 * t307 - t265 * t309 + (-t275 * t310 + t269) * qJD(2)) * MDP(14) + (qJD(2) * t270 - t239 * t320 - t250 * t310 - t265 * t312 - t275 * t300) * MDP(15) + (t239 * t275 + t250 * t265 - t253 * t281 + t269 * t273 - t270 * t274 + t411) * MDP(16) + (-t262 * t407 + t285 * t359) * MDP(17) + ((-t283 * t343 + t285 * t345) * t309 + (-t410 - t263 * t343 + (-t283 * t345 - t285 * t343) * qJD(5)) * t319) * MDP(18) + (-t262 * t320 + t285 * t312 + t319 * t403 + t359 * t429) * MDP(19) + (-t263 * t320 - t283 * t312 + t294 * t319 - t358 * t429) * MDP(20) + (t300 * t320 + t312 * t429) * MDP(21) + ((-t234 * t343 + t245) * t429 + (-t257 * t343 + t267) * t300 + t370 * t320 + t226 * t312 + t249 * t283 + t272 * t263 + (-t235 * t319 - t251 * t309) * t345 + (-t227 * t320 + t251 * t407 - t396 * t429) * qJD(5)) * MDP(22) + (-t227 * t312 + t235 * t407 + t249 * t285 + t251 * t359 - t272 * t262 - t300 * t396 - t320 * t356 - t355 * t429) * MDP(23) + (-t218 * t285 - t219 * t283 + t225 * t262 - t229 * t263 + (-t220 * t343 + t222 * t345) * t309 + (-t216 * t343 + t217 * t345 + (-t220 * t345 - t222 * t343) * qJD(5)) * t319) * MDP(24) + (t217 * t229 + t222 * t219 + t216 * t225 + t220 * t218 + t228 * ((-pkin(5) * t345 - pkin(4)) * t319 + t281) + t232 * (pkin(5) * t358 + t249)) * MDP(25); -t400 * t427 + t348 * t426 + ((t277 - t278) * t310 + (t279 - t276) * t307 + (-t342 * t300 - t341 * t353) * pkin(2)) * MDP(11) + (t276 * t278 - t277 * t279 + (-t258 * t342 + t259 * t341 - t324 * t385) * pkin(2)) * MDP(12) + (t334 * t300 + (-t274 - t278) * t310 - t332 * t352 + (t273 + t384) * t307) * MDP(13) + (-qJD(2) * t278 + t268 * t307 + t360) * MDP(14) + (-qJD(2) * t279 - t265 * t307 + t268 * t310 + t259 + 0.2e1 * t337) * MDP(15) + (-t253 * t332 + t258 * t334 - t265 * t268 - t273 * t278 + t274 * t384) * MDP(16) + (-t285 * t368 - t410) * MDP(17) + (-t256 - t434 + (t262 + t365) * t343) * MDP(18) + (t408 - t430) * MDP(19) + (t354 - t409) * MDP(20) + t429 * t307 * MDP(21) + (t226 * t307 + t235 * t343 + t332 * t263 + (-t255 + (-qJD(5) * t330 + t247) * t343) * t429 + t383 * t283 + t419 * t345) * MDP(22) + (-t227 * t307 + t235 * t345 - t332 * t262 + (-t330 * t388 + t397) * t429 + t383 * t285 - t419 * t343) * MDP(23) + (-t262 * t317 + t263 * t316 - t394 * t283 - t395 * t285 - t343 * t432 - t345 * t433) * MDP(24) + (-t217 * t316 - t216 * t317 + t228 * (pkin(5) * t343 + t332) + (pkin(5) * t422 + t383) * t232 + t394 * t222 + t395 * t220) * MDP(25) + (MDP(9) * t344 * t348 + MDP(10) * t400) * pkin(1); (t276 * t310 + t277 * t307 + t331) * MDP(12) - 0.2e1 * MDP(14) * t353 + (-t300 + t302) * MDP(15) + (-t274 * t307 + (-qJD(4) - t273) * t310 + t349) * MDP(16) + (t354 + t409) * MDP(22) + (t408 + t430) * MDP(23) + (-t431 * t343 - t256 + t434) * MDP(24) + (t232 * t307 - t343 * t433 + t345 * t432) * MDP(25) + (MDP(11) + MDP(13)) * (-t307 ^ 2 - t306); (t302 - t327) * MDP(13) - t310 * t307 * MDP(14) + (-t306 - t347) * MDP(15) + t360 * MDP(16) + t294 * MDP(22) + (MDP(13) * t378 + t274 * MDP(16) - t283 * MDP(22) - t285 * MDP(23) - t232 * MDP(25)) * qJD(2) + (-MDP(23) * t369 + t431 * MDP(24) + t433 * MDP(25)) * t345 + (-t300 * MDP(23) + (t285 * t310 - t263 + t390) * MDP(24) + t432 * MDP(25) - MDP(22) * t369) * t343; (-t282 + t418) * MDP(18) + t393 * MDP(19) + (t285 * t429 + t293) * MDP(20) + t300 * MDP(21) + (t227 * t429 - t251 * t285 + t370) * MDP(22) + (t226 * t429 + t380) * MDP(23) + t398 * MDP(25) * t222 + (t262 * MDP(24) + (-t232 * t285 + t216) * MDP(25)) * pkin(5) + (t285 * MDP(17) + MDP(19) * t429 + t251 * MDP(23) - MDP(24) * t398) * t283 + ((-MDP(20) * qJD(2) - MDP(22) * t241) * t345 + (-MDP(19) * qJD(2) - t307 * MDP(20) - MDP(22) * t246 + MDP(23) * t241) * t343) * qJD(5); (-t282 - t418) * MDP(24) + (t220 * t285 + t222 * t283 + t228) * MDP(25);];
tauc  = t1;
