% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:35
% EndTime: 2019-03-08 21:06:43
% DurationCPUTime: 3.48s
% Computational Cost: add. (2223->337), mult. (5837->460), div. (0->0), fcn. (4383->10), ass. (0->157)
t316 = sin(pkin(11));
t321 = sin(qJ(3));
t376 = qJD(2) * t321;
t318 = cos(pkin(11));
t324 = cos(qJ(3));
t387 = t318 * t324;
t287 = -qJD(2) * t387 + t316 * t376;
t323 = cos(qJ(6));
t320 = sin(qJ(6));
t374 = qJD(3) * t320;
t268 = -t323 * t287 + t374;
t295 = t316 * t324 + t318 * t321;
t290 = t295 * qJD(2);
t415 = qJD(6) + t290;
t417 = t268 * t415;
t270 = qJD(3) * t323 + t287 * t320;
t350 = t415 * t270;
t322 = sin(qJ(2));
t317 = sin(pkin(6));
t378 = qJD(1) * t317;
t360 = t322 * t378;
t373 = qJD(3) * t321;
t416 = pkin(3) * t373 - t360;
t414 = MDP(12) + MDP(14);
t413 = MDP(5) * t321;
t412 = t320 * t415;
t411 = (t321 ^ 2 - t324 ^ 2) * MDP(6);
t400 = -qJ(4) - pkin(8);
t354 = qJD(3) * t400;
t283 = qJD(4) * t324 + t321 * t354;
t284 = -qJD(4) * t321 + t324 * t354;
t325 = cos(qJ(2));
t359 = t325 * t378;
t382 = t283 * t316 - t318 * t284 - t295 * t359;
t407 = -t316 * t321 + t387;
t381 = t283 * t318 + t284 * t316 - t407 * t359;
t372 = qJD(3) * t324;
t292 = -t316 * t373 + t318 * t372;
t409 = qJ(5) * t292 + qJD(5) * t295 - t416;
t408 = t321 * MDP(10) + t324 * MDP(11);
t298 = qJD(2) * pkin(8) + t360;
t349 = qJ(4) * qJD(2) + t298;
t319 = cos(pkin(6));
t377 = qJD(1) * t319;
t358 = t321 * t377;
t263 = t349 * t324 + t358;
t255 = t316 * t263;
t306 = t324 * t377;
t262 = -t349 * t321 + t306;
t232 = t262 * t318 - t255;
t365 = -qJD(5) + t232;
t259 = qJD(3) * pkin(3) + t262;
t388 = t318 * t263;
t228 = t316 * t259 + t388;
t225 = -qJD(3) * qJ(5) - t228;
t402 = pkin(5) * t287;
t219 = -t225 - t402;
t230 = t262 * t316 + t388;
t362 = qJD(2) * qJD(3);
t356 = t321 * t362;
t302 = t316 * t356;
t355 = t324 * t362;
t279 = t318 * t355 - t302;
t311 = -pkin(3) * t318 - pkin(4);
t307 = -pkin(9) + t311;
t406 = t307 * t279 + (t219 - t230 + t402) * t415;
t286 = t290 ^ 2;
t404 = pkin(4) + pkin(9);
t289 = t295 * qJD(3);
t278 = qJD(2) * t289;
t403 = pkin(4) * t278;
t401 = pkin(5) * t290;
t399 = qJD(2) * pkin(2);
t342 = qJD(4) + t359;
t241 = (-t298 * t321 + t306) * qJD(3) + (-qJ(4) * t373 + t342 * t324) * qJD(2);
t242 = (-t298 * t324 - t358) * qJD(3) + (-qJ(4) * t372 - t342 * t321) * qJD(2);
t215 = t241 * t316 - t318 * t242;
t389 = t317 * t322;
t293 = t319 * t321 + t324 * t389;
t339 = t319 * t324 - t321 * t389;
t251 = t293 * t316 - t318 * t339;
t398 = t215 * t251;
t300 = t400 * t321;
t301 = t400 * t324;
t266 = -t318 * t300 - t301 * t316;
t397 = t215 * t266;
t371 = qJD(6) * t320;
t370 = qJD(6) * t323;
t380 = t320 * t278 + t287 * t370;
t244 = -qJD(3) * t371 + t380;
t396 = t244 * t323;
t395 = t268 * t287;
t394 = t270 * t287;
t393 = t279 * t320;
t392 = t279 * t325;
t391 = t407 * t320;
t326 = qJD(3) ^ 2;
t386 = t321 * t326;
t273 = t323 * t279;
t385 = t324 * t326;
t384 = pkin(5) * t292 + t382;
t383 = -pkin(5) * t289 + t381;
t216 = t318 * t241 + t316 * t242;
t285 = pkin(3) * t356 + qJD(2) * t360;
t375 = qJD(2) * t322;
t369 = qJD(6) * t325;
t368 = t290 * MDP(15);
t364 = t401 - t365;
t361 = -pkin(3) * t324 - pkin(2);
t357 = qJD(2) * t317 * t325;
t353 = pkin(3) * t376 + qJ(5) * t287;
t213 = pkin(5) * t279 + t215;
t347 = -qJ(5) * t279 + t285;
t333 = -qJD(5) * t290 + t347;
t217 = t404 * t278 + t333;
t352 = t323 * t213 - t217 * t320;
t227 = t259 * t318 - t255;
t351 = t323 * t415;
t348 = MDP(23) * t415;
t346 = -t404 * t289 + t409;
t345 = -pkin(4) * t289 + t409;
t343 = qJD(5) - t227;
t214 = -qJD(3) * qJD(5) - t216;
t218 = -t404 * qJD(3) + t343 + t401;
t280 = t361 * qJD(2) + qJD(4) - t359;
t329 = -qJ(5) * t290 + t280;
t226 = t404 * t287 + t329;
t209 = t218 * t323 - t226 * t320;
t210 = t218 * t320 + t226 * t323;
t267 = t300 * t316 - t301 * t318;
t340 = -qJ(5) * t295 + t361;
t240 = pkin(4) * t287 + t329;
t338 = t240 * t290 + t215;
t337 = t289 * t320 - t370 * t407;
t211 = -pkin(5) * t278 - t214;
t334 = t211 + (-qJD(6) * t307 + t404 * t290 + t353) * t415;
t249 = pkin(5) * t295 + t266;
t332 = -t211 * t407 + t219 * t289 - t249 * t279;
t331 = -0.2e1 * qJD(3) * t399;
t328 = t215 * t295 + t266 * t279 - t267 * t278 - t381 * t287 + t382 * t290;
t327 = qJD(2) ^ 2;
t309 = pkin(3) * t316 + qJ(5);
t281 = qJD(3) * t287;
t272 = t323 * t278;
t261 = -t293 * qJD(3) - t321 * t357;
t260 = t339 * qJD(3) + t324 * t357;
t253 = -pkin(4) * t407 + t340;
t252 = t318 * t293 + t316 * t339;
t250 = pkin(5) * t407 + t267;
t246 = pkin(4) * t290 + t353;
t245 = t270 * qJD(6) - t272;
t243 = -t404 * t407 + t340;
t231 = t260 * t318 + t261 * t316;
t229 = t260 * t316 - t318 * t261;
t224 = -qJD(3) * pkin(4) + t343;
t222 = t333 + t403;
t1 = [(t216 * t252 - t227 * t229 + t228 * t231 + t398) * MDP(13) + (-t214 * t252 + t224 * t229 - t225 * t231 + t398) * MDP(17) + ((t229 * t323 - t251 * t371) * t415 + t251 * t273 + t231 * t268 + t252 * t245) * MDP(23) + (-(t229 * t320 + t251 * t370) * t415 - t251 * t393 + t231 * t270 + t252 * t244) * MDP(24) + (MDP(10) * t261 - MDP(11) * t260 + MDP(15) * t229 + MDP(16) * t231) * qJD(3) + ((t280 * t375 - t285 * t325) * MDP(13) + (t278 * t325 - t287 * t375) * MDP(15) + (-t290 * t375 + t392) * MDP(16) + (-t222 * t325 + t240 * t375) * MDP(17) + ((-t320 * t375 + t323 * t369) * t415 + t320 * t392) * MDP(23) + (-(t320 * t369 + t323 * t375) * t415 + t325 * t273) * MDP(24) - t408 * t325 * t362 + (-t325 * MDP(4) + (-MDP(10) * t324 + MDP(11) * t321 - MDP(3)) * t322) * t327) * t317 + t414 * (t229 * t290 - t231 * t287 + t251 * t279 - t252 * t278); 0.2e1 * t355 * t413 - 0.2e1 * t362 * t411 + MDP(7) * t385 - MDP(8) * t386 + (-pkin(8) * t385 + t321 * t331) * MDP(10) + (pkin(8) * t386 + t324 * t331) * MDP(11) + (t216 * t407 - t227 * t292 - t228 * t289 + t328) * MDP(12) + (t216 * t267 - t382 * t227 + t381 * t228 + t416 * t280 + t285 * t361 + t397) * MDP(13) + (-t214 * t407 + t224 * t292 + t225 * t289 + t328) * MDP(14) + (t382 * qJD(3) + t222 * t407 - t240 * t289 - t253 * t278 + t345 * t287) * MDP(15) + (t381 * qJD(3) - t222 * t295 - t240 * t292 - t253 * t279 + t345 * t290) * MDP(16) + (-t214 * t267 + t222 * t253 + t382 * t224 - t381 * t225 - t345 * t240 + t397) * MDP(17) + (-t244 * t391 + t337 * t270) * MDP(18) + ((-t268 * t320 + t270 * t323) * t289 - (t396 - t245 * t320 + (-t268 * t323 - t270 * t320) * qJD(6)) * t407) * MDP(19) + (t244 * t295 + t270 * t292 - t279 * t391 + t337 * t415) * MDP(20) + (-t407 * t273 - t245 * t295 - t268 * t292 + (t289 * t323 + t371 * t407) * t415) * MDP(21) + (t279 * t295 + t292 * t415) * MDP(22) + (-t243 * t393 + t352 * t295 + t209 * t292 + t250 * t245 - t332 * t323 + (t346 * t320 + t384 * t323) * t415 + t383 * t268 + ((-t243 * t323 - t249 * t320) * t415 - t210 * t295 - t219 * t391) * qJD(6)) * MDP(23) + (-t210 * t292 + t250 * t244 + t383 * t270 + (-t243 * t279 - (qJD(6) * t218 + t217) * t295 - t219 * qJD(6) * t407 + (-qJD(6) * t249 + t346) * t415) * t323 + (-(-qJD(6) * t226 + t213) * t295 + (qJD(6) * t243 - t384) * t415 + t332) * t320) * MDP(24); ((t228 - t230) * t290 + (-t227 + t232) * t287 + (-t278 * t316 - t279 * t318) * pkin(3)) * MDP(12) + (t227 * t230 - t228 * t232 + (-t215 * t318 + t216 * t316 - t280 * t376) * pkin(3)) * MDP(13) + (-t278 * t309 + t279 * t311 + (-t225 - t230) * t290 + (t224 + t365) * t287) * MDP(14) + (-qJD(3) * t230 + t246 * t287 + t338) * MDP(15) + (-t240 * t287 + t246 * t290 + (0.2e1 * qJD(5) - t232) * qJD(3) + t216) * MDP(16) + (-t214 * t309 + t215 * t311 - t224 * t230 + t365 * t225 - t240 * t246) * MDP(17) + (-t320 * t350 + t396) * MDP(18) + ((-t245 - t350) * t323 + (-t244 + t417) * t320) * MDP(19) + (-t412 * t415 + t273 + t394) * MDP(20) + (-t351 * t415 - t393 - t395) * MDP(21) + t415 * t287 * MDP(22) + (t209 * t287 + t309 * t245 + t364 * t268 + t334 * t320 + t323 * t406) * MDP(23) + (-t210 * t287 + t309 * t244 + t364 * t270 - t320 * t406 + t334 * t323) * MDP(24) + t408 * qJD(2) * t399 + (-t324 * t413 + t411) * t327; (t227 * t290 + t228 * t287 + t285) * MDP(13) + (t281 + t302) * MDP(16) + (t403 - t225 * t287 + (-qJD(5) - t224) * t290 + t347) * MDP(17) + (-t393 + t395) * MDP(23) + (-t273 + t394) * MDP(24) + (MDP(24) * t412 - t323 * t348) * t415 + (-t368 + (-t295 * MDP(15) - MDP(16) * t387) * qJD(2)) * qJD(3) + t414 * (-t287 ^ 2 - t286); (t279 + t281) * MDP(14) - t287 * t368 + (-t286 - t326) * MDP(16) + (qJD(3) * t225 + t338) * MDP(17) + (-qJD(3) * t268 + t273) * MDP(23) + (-qJD(3) * t270 - t393) * MDP(24) + (-MDP(24) * t351 - t320 * t348) * t415; t270 * t268 * MDP(18) + (-t268 ^ 2 + t270 ^ 2) * MDP(19) + (t380 + t417) * MDP(20) + (t272 + t350) * MDP(21) + t279 * MDP(22) + (t210 * t415 - t219 * t270 + t352) * MDP(23) + (t209 * t415 - t213 * t320 - t217 * t323 + t219 * t268) * MDP(24) + (-MDP(20) * t374 - t270 * MDP(21) - t210 * MDP(23) - t209 * MDP(24)) * qJD(6);];
tauc  = t1;
