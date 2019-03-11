% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:09
% EndTime: 2019-03-08 19:45:15
% DurationCPUTime: 3.14s
% Computational Cost: add. (1944->326), mult. (5108->432), div. (0->0), fcn. (4005->10), ass. (0->145)
t317 = cos(pkin(11));
t400 = cos(qJ(4));
t360 = t400 * t317;
t347 = qJD(2) * t360;
t320 = sin(qJ(4));
t315 = sin(pkin(11));
t371 = qJD(2) * t315;
t357 = t320 * t371;
t286 = -t347 + t357;
t322 = cos(qJ(6));
t319 = sin(qJ(6));
t369 = qJD(4) * t319;
t264 = -t322 * t286 + t369;
t297 = t400 * t315 + t320 * t317;
t328 = qJD(2) * t297;
t413 = qJD(6) + t328;
t415 = t264 * t413;
t266 = qJD(4) * t322 + t286 * t319;
t349 = t413 * t266;
t305 = qJD(4) * t347;
t368 = qJD(4) * t320;
t356 = t315 * t368;
t279 = qJD(2) * t356 - t305;
t270 = t322 * t279;
t351 = t319 * t413;
t414 = -t351 * t413 - t270;
t385 = t315 * t320;
t296 = -t360 + t385;
t321 = sin(qJ(2));
t316 = sin(pkin(6));
t373 = qJD(1) * t316;
t359 = t321 * t373;
t301 = qJD(2) * qJ(3) + t359;
t318 = cos(pkin(6));
t372 = qJD(1) * t318;
t307 = t317 * t372;
t396 = pkin(8) * qJD(2);
t258 = t307 + (-t301 - t396) * t315;
t268 = t317 * t301 + t315 * t372;
t259 = t317 * t396 + t268;
t377 = t400 * t258 - t320 * t259;
t407 = qJD(5) - t377;
t374 = t315 ^ 2 + t317 ^ 2;
t411 = t374 * MDP(7);
t323 = cos(qJ(2));
t383 = t316 * t323;
t327 = t297 * t383;
t397 = pkin(8) + qJ(3);
t302 = t397 * t315;
t303 = t397 * t317;
t337 = -t320 * t302 + t400 * t303;
t378 = -qJD(1) * t327 + t297 * qJD(3) + t337 * qJD(4);
t355 = qJD(4) * t400;
t290 = -t317 * t355 + t356;
t410 = -qJ(5) * t290 + qJD(5) * t297 + t359;
t409 = t378 * qJD(4);
t343 = t360 * t383;
t358 = t323 * t373;
t379 = qJD(1) * t343 + (qJD(3) * t315 + qJD(4) * t303) * t320 - t358 * t385 - qJD(3) * t360 + t302 * t355;
t408 = t379 * qJD(4);
t344 = qJD(3) - t358;
t406 = MDP(14) - MDP(17);
t405 = MDP(15) - MDP(18);
t232 = t320 * t258 + t400 * t259;
t226 = -qJD(4) * qJ(5) - t232;
t398 = pkin(5) * t286;
t221 = -t226 - t398;
t401 = pkin(4) + pkin(9);
t404 = t401 * t279 + (t221 - t232 + t398) * t413;
t403 = t286 ^ 2;
t402 = t328 ^ 2;
t291 = t297 * qJD(4);
t280 = qJD(2) * t291;
t399 = pkin(4) * t280;
t395 = qJD(2) * pkin(2);
t394 = qJ(5) * t286;
t312 = -pkin(3) * t317 - pkin(2);
t281 = t312 * qJD(2) + t344;
t326 = -qJ(5) * t328 + t281;
t235 = pkin(4) * t286 + t326;
t393 = t235 * t328;
t367 = qJD(6) * t319;
t366 = qJD(6) * t322;
t376 = t319 * t280 + t286 * t366;
t239 = -qJD(4) * t367 + t376;
t392 = t239 * t322;
t391 = t264 * t286;
t390 = t266 * t286;
t389 = t279 * t319;
t388 = t279 * t323;
t387 = t286 * t328;
t386 = t296 * t319;
t384 = t316 * t321;
t381 = -pkin(5) * t291 - t379;
t380 = -t290 * pkin(5) + t378;
t370 = qJD(2) * t321;
t365 = qJD(6) * t323;
t363 = pkin(5) * t328 + t407;
t295 = (qJD(3) + t358) * qJD(2);
t354 = t374 * t295;
t304 = qJD(2) * t359;
t353 = qJ(5) * t279 + t304;
t218 = t258 * t368 + t259 * t355 + t297 * t295;
t216 = -pkin(5) * t279 + t218;
t333 = -qJD(5) * t328 + t353;
t220 = t401 * t280 + t333;
t352 = t322 * t216 - t220 * t319;
t350 = t322 * t413;
t348 = -t258 * t355 + t259 * t368 + t296 * t295;
t260 = t400 * t302 + t320 * t303;
t346 = -t401 * t291 + t410;
t345 = -pkin(4) * t291 + t410;
t219 = -t401 * qJD(4) + t363;
t224 = t401 * t286 + t326;
t212 = t219 * t322 - t224 * t319;
t213 = t219 * t319 + t224 * t322;
t341 = (-t301 * t315 + t307) * t315 - t268 * t317;
t338 = -qJ(5) * t297 + t312;
t284 = -t315 * t384 + t317 * t318;
t285 = t315 * t318 + t317 * t384;
t246 = t320 * t284 + t400 * t285;
t334 = t291 * t319 + t296 * t366;
t332 = qJD(4) * t232 - t218;
t217 = -qJD(4) * qJD(5) + t348;
t214 = -pkin(5) * t280 - t217;
t331 = t214 + (t413 * t401 + t394) * t413;
t243 = t297 * pkin(5) + t260;
t330 = t214 * t296 + t221 * t291 + t243 * t279;
t329 = -t350 * t413 + t389;
t325 = qJD(2) ^ 2;
t300 = t344 - t395;
t282 = qJD(4) * t286;
t271 = t322 * t280;
t257 = t279 * t297;
t248 = pkin(4) * t296 + t338;
t247 = pkin(4) * t328 + t394;
t245 = -t400 * t284 + t320 * t285;
t244 = -t296 * pkin(5) + t337;
t240 = t266 * qJD(6) - t271;
t238 = t401 * t296 + t338;
t230 = t333 + t399;
t228 = qJD(2) * t327 + t246 * qJD(4);
t227 = -t284 * t355 - qJD(2) * t343 + (qJD(4) * t285 + t371 * t383) * t320;
t225 = -qJD(4) * pkin(4) + t407;
t1 = [(t227 * t286 + t228 * t328 - t245 * t279 - t246 * t280) * MDP(16) + (-t217 * t246 + t218 * t245 + t225 * t228 + t226 * t227) * MDP(19) + ((t228 * t322 - t245 * t367) * t413 - t245 * t270 - t227 * t264 + t246 * t240) * MDP(25) + (-(t228 * t319 + t245 * t366) * t413 + t245 * t389 - t227 * t266 + t246 * t239) * MDP(26) + (-t284 * t315 + t285 * t317) * MDP(8) * t295 + (t405 * t227 - t406 * t228) * qJD(4) + ((-t230 * t323 + t235 * t370) * MDP(19) + ((-t319 * t370 + t322 * t365) * t413 - t319 * t388) * MDP(25) + (-(t319 * t365 + t322 * t370) * t413 - t323 * t270) * MDP(26) + (t300 * t321 + (-t341 - t359) * t323) * qJD(2) * MDP(8) + ((-MDP(4) + t411) * t323 + (-MDP(5) * t317 + MDP(6) * t315 - MDP(3)) * t321) * t325 - t405 * (-t328 * t370 - t388) + t406 * (-t280 * t323 + t286 * t370)) * t316; (t344 * qJD(2) * t374 + t354) * MDP(7) + (-t341 * qJD(3) + qJ(3) * t354 + (t341 * t323 + (-t300 - t395) * t321) * t373) * MDP(8) + (-t290 * t328 - t257) * MDP(9) + (t279 * t296 - t280 * t297 + t286 * t290 - t291 * t328) * MDP(10) + (t280 * t312 + t281 * t291 - t409 + (qJD(2) * t296 - t286) * t359) * MDP(14) + (-t279 * t312 - t281 * t290 + t408) * MDP(15) + (t217 * t296 + t218 * t297 - t225 * t290 + t226 * t291 - t260 * t279 - t280 * t337 + t379 * t286 + t328 * t378) * MDP(16) + (-t230 * t296 - t235 * t291 - t248 * t280 + t345 * t286 + t409) * MDP(17) + (-t230 * t297 + t235 * t290 + t248 * t279 + t328 * t345 - t408) * MDP(18) + (-t217 * t337 + t218 * t260 + t378 * t225 + t379 * t226 + t230 * t248 - t345 * t235) * MDP(19) + (t239 * t386 + t334 * t266) * MDP(20) + ((-t264 * t319 + t266 * t322) * t291 + (t392 - t240 * t319 + (-t264 * t322 - t266 * t319) * qJD(6)) * t296) * MDP(21) + (t239 * t297 - t266 * t290 - t279 * t386 + t334 * t413) * MDP(22) + (-t296 * t270 - t240 * t297 + t264 * t290 + (t291 * t322 - t296 * t367) * t413) * MDP(23) + (-t290 * t413 - t257) * MDP(24) + (t238 * t389 + t352 * t297 - t212 * t290 + t244 * t240 - t330 * t322 + (t346 * t319 + t380 * t322) * t413 + t381 * t264 + ((-t238 * t322 - t243 * t319) * t413 - t213 * t297 + t221 * t386) * qJD(6)) * MDP(25) + (t213 * t290 + t244 * t239 + t381 * t266 + (t238 * t279 - (qJD(6) * t219 + t220) * t297 + t221 * qJD(6) * t296 + (-qJD(6) * t243 + t346) * t413) * t322 + (-(-qJD(6) * t224 + t216) * t297 + (qJD(6) * t238 - t380) * t413 + t330) * t319) * MDP(26) + (-t290 * MDP(11) - t291 * MDP(12)) * qJD(4); -t325 * t411 + (t341 * qJD(2) + t304) * MDP(8) + t305 * MDP(15) + (-t402 - t403) * MDP(16) + (t279 + t282) * MDP(18) + (t399 - t226 * t286 + (-qJD(5) - t225) * t328 + t353) * MDP(19) + (t329 + t391) * MDP(25) + (t390 - t414) * MDP(26) + ((-t286 - t357) * MDP(15) + (0.2e1 * MDP(14) - 0.2e1 * MDP(17)) * t328) * qJD(4); MDP(9) * t387 + (t402 - t403) * MDP(10) + (t305 + (t286 - t357) * qJD(4)) * MDP(11) + (-t281 * t328 + t332) * MDP(14) + (qJD(4) * t377 + t281 * t286 + t348) * MDP(15) + (pkin(4) * t279 - qJ(5) * t280 + (-t226 - t232) * t328 + (t225 - t407) * t286) * MDP(16) + (t247 * t286 - t332 + t393) * MDP(17) + (-t235 * t286 + t247 * t328 + (0.2e1 * qJD(5) - t377) * qJD(4) - t348) * MDP(18) + (-pkin(4) * t218 - qJ(5) * t217 - t225 * t232 - t226 * t407 - t235 * t247) * MDP(19) + (-t319 * t349 + t392) * MDP(20) + ((-t240 - t349) * t322 + (-t239 + t415) * t319) * MDP(21) + (t390 + t414) * MDP(22) + (t329 - t391) * MDP(23) + t413 * t286 * MDP(24) + (qJ(5) * t240 + t212 * t286 + t363 * t264 + t331 * t319 + t404 * t322) * MDP(25) + (qJ(5) * t239 - t213 * t286 + t363 * t266 - t404 * t319 + t331 * t322) * MDP(26); (-t279 + t282) * MDP(16) - MDP(17) * t387 + (-qJD(4) ^ 2 - t402) * MDP(18) + (qJD(4) * t226 + t218 + t393) * MDP(19) + (-qJD(4) * t264 - t270) * MDP(25) + (-qJD(4) * t266 + t389) * MDP(26) + (-MDP(25) * t351 - MDP(26) * t350) * t413; t266 * t264 * MDP(20) + (-t264 ^ 2 + t266 ^ 2) * MDP(21) + (t376 + t415) * MDP(22) + (t271 + t349) * MDP(23) - t279 * MDP(24) + (t213 * t413 - t221 * t266 + t352) * MDP(25) + (t212 * t413 - t216 * t319 - t220 * t322 + t221 * t264) * MDP(26) + (-MDP(22) * t369 - t266 * MDP(23) - t213 * MDP(25) - t212 * MDP(26)) * qJD(6);];
tauc  = t1;
