% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:32
% EndTime: 2019-12-31 21:20:38
% DurationCPUTime: 2.87s
% Computational Cost: add. (2215->308), mult. (5444->402), div. (0->0), fcn. (3684->6), ass. (0->144)
t323 = cos(qJ(2));
t393 = cos(qJ(3));
t358 = t393 * t323;
t344 = qJD(1) * t358;
t320 = sin(qJ(3));
t321 = sin(qJ(2));
t368 = qJD(1) * t321;
t357 = t320 * t368;
t278 = -t344 + t357;
t322 = cos(qJ(5));
t315 = qJD(2) + qJD(3);
t319 = sin(qJ(5));
t380 = t315 * t319;
t259 = -t322 * t278 + t380;
t291 = t320 * t323 + t321 * t393;
t369 = qJD(1) * t291;
t399 = qJD(5) + t369;
t401 = t259 * t399;
t261 = t278 * t319 + t315 * t322;
t349 = t399 * t261;
t400 = t399 * t322;
t361 = qJD(1) * qJD(2);
t398 = -0.2e1 * t361;
t397 = MDP(5) * (t321 ^ 2 - t323 ^ 2);
t394 = -pkin(7) - pkin(6);
t299 = t394 * t323;
t295 = qJD(1) * t299;
t285 = t393 * t295;
t298 = t394 * t321;
t293 = qJD(1) * t298;
t255 = t320 * t293 - t285;
t367 = qJD(3) * t320;
t342 = pkin(2) * t367 - t255;
t282 = t320 * t295;
t256 = t293 * t393 + t282;
t356 = qJD(3) * t393;
t373 = -pkin(2) * t356 - qJD(4) + t256;
t390 = qJD(2) * pkin(2);
t286 = t293 + t390;
t253 = -t393 * t286 - t282;
t363 = qJD(4) + t253;
t396 = t369 ^ 2;
t395 = pkin(3) + pkin(8);
t392 = t278 * pkin(4);
t391 = t369 * pkin(4);
t366 = qJD(5) * t319;
t258 = t315 * t291;
t250 = t258 * qJD(1);
t365 = qJD(5) * t322;
t374 = t319 * t250 + t278 * t365;
t221 = -t315 * t366 + t374;
t389 = t221 * t322;
t359 = qJD(2) * t394;
t294 = t321 * t359;
t296 = t323 * t359;
t226 = -t393 * t294 - t320 * t296 - t298 * t356 - t299 * t367;
t388 = t226 * t315;
t338 = t320 * t298 - t299 * t393;
t227 = qJD(3) * t338 + t320 * t294 - t296 * t393;
t387 = t227 * t315;
t379 = t320 * t321;
t290 = -t358 + t379;
t311 = -pkin(2) * t323 - pkin(1);
t340 = -qJ(4) * t291 + t311;
t237 = t290 * t395 + t340;
t341 = t315 * t379;
t372 = t315 * t344;
t249 = qJD(1) * t341 - t372;
t386 = t237 * t249;
t385 = t249 * t319;
t254 = t320 * t286 - t285;
t384 = t254 * t315;
t382 = t278 * t369;
t381 = t290 * t319;
t325 = qJD(2) ^ 2;
t378 = t321 * t325;
t247 = t322 * t249;
t377 = t323 * t325;
t326 = qJD(1) ^ 2;
t376 = t323 * t326;
t375 = t395 * t249;
t371 = t391 - t373;
t364 = t391 + t363;
t313 = t321 * t390;
t355 = t321 * t361;
t354 = pkin(1) * t398;
t223 = -t315 * t395 + t364;
t297 = t311 * qJD(1);
t329 = -qJ(4) * t369 + t297;
t225 = t278 * t395 + t329;
t207 = t223 * t319 + t225 * t322;
t345 = qJD(1) * t359;
t287 = t321 * t345;
t288 = t323 * t345;
t346 = t286 * t356 + t393 * t287 + t320 * t288 + t295 * t367;
t218 = -t315 * qJD(4) - t346;
t211 = -pkin(4) * t250 - t218;
t353 = -t207 * t278 + t211 * t322;
t352 = t319 * t399;
t245 = -qJ(4) * t315 - t254;
t231 = -t245 - t392;
t350 = t399 * t231;
t251 = pkin(3) * t369 + qJ(4) * t278;
t242 = pkin(2) * t368 + t251;
t274 = t369 * pkin(8);
t310 = -pkin(2) * t393 - pkin(3);
t307 = -pkin(8) + t310;
t348 = -qJD(5) * t307 + t242 + t274;
t347 = qJD(5) * t395 + t251 + t274;
t220 = t286 * t367 + t320 * t287 - t393 * t288 - t295 * t356;
t343 = t392 + t342;
t262 = -t298 * t393 - t320 * t299;
t206 = t223 * t322 - t225 * t319;
t339 = t206 * t278 + t211 * t319 + (t322 * t369 + t365) * t231;
t337 = t258 * t319 + t290 * t365;
t336 = pkin(2) * t355 + qJ(4) * t249 - qJD(4) * t369;
t257 = -qJD(2) * t358 - t323 * t356 + t341;
t335 = qJ(4) * t257 - qJD(4) * t291 + t313;
t240 = pkin(3) * t278 + t329;
t334 = t240 * t369 + t220;
t333 = -t297 * t369 - t220;
t332 = t297 * t278 - t346;
t331 = -t240 * t278 - t218;
t243 = t291 * pkin(4) + t262;
t330 = t211 * t290 + t231 * t258 + t243 * t249;
t328 = t372 + (t278 - t357) * t315;
t248 = t322 * t250;
t222 = qJD(5) * t261 - t248;
t327 = t399 * t278 * MDP(26) + MDP(11) * t382 + ((-t222 - t349) * t322 + (-t221 + t401) * t319) * MDP(23) + (-t319 * t349 + t389) * MDP(22) + (t261 * t278 - t352 * t399 - t247) * MDP(24) + (-t259 * t278 - t399 * t400 + t385) * MDP(25) + t328 * MDP(13) + (-t278 ^ 2 + t396) * MDP(12);
t308 = pkin(2) * t320 + qJ(4);
t252 = pkin(3) * t290 + t340;
t244 = -t290 * pkin(4) + t338;
t241 = -pkin(3) * t315 + t363;
t236 = t249 * t291;
t234 = t254 - t392;
t217 = pkin(3) * t258 + t335;
t216 = -t257 * pkin(4) + t227;
t215 = -pkin(4) * t258 - t226;
t214 = pkin(3) * t250 + t336;
t213 = -pkin(4) * t249 + t220;
t212 = t322 * t213;
t210 = t258 * t395 + t335;
t205 = t250 * t395 + t336;
t1 = [0.2e1 * t323 * MDP(4) * t355 + t397 * t398 + MDP(6) * t377 - MDP(7) * t378 + (-pkin(6) * t377 + t321 * t354) * MDP(9) + (pkin(6) * t378 + t323 * t354) * MDP(10) + (-t257 * t369 - t236) * MDP(11) + (t249 * t290 - t250 * t291 + t257 * t278 - t258 * t369) * MDP(12) + (-t387 + t250 * t311 + t258 * t297 + (qJD(1) * t290 + t278) * t313) * MDP(16) + (-t249 * t311 - t257 * t297 + 0.2e1 * t369 * t313 + t388) * MDP(17) + (t218 * t290 + t220 * t291 + t226 * t278 + t227 * t369 - t241 * t257 + t245 * t258 - t249 * t262 - t250 * t338) * MDP(18) + (-t214 * t290 - t217 * t278 - t240 * t258 - t250 * t252 + t387) * MDP(19) + (-t214 * t291 - t217 * t369 + t240 * t257 + t249 * t252 - t388) * MDP(20) + (t214 * t252 + t217 * t240 - t218 * t338 + t220 * t262 + t226 * t245 + t227 * t241) * MDP(21) + (t221 * t381 + t261 * t337) * MDP(22) + ((-t259 * t319 + t261 * t322) * t258 + (t389 - t222 * t319 + (-t259 * t322 - t261 * t319) * qJD(5)) * t290) * MDP(23) + (t221 * t291 - t249 * t381 - t257 * t261 + t337 * t399) * MDP(24) + (-t290 * t247 - t222 * t291 + t257 * t259 + (t258 * t322 - t290 * t366) * t399) * MDP(25) + (-t257 * t399 - t236) * MDP(26) + (-t206 * t257 + t212 * t291 + t215 * t259 + t244 * t222 + (-t205 * t291 - t210 * t399 + t386) * t319 + (t216 * t399 - t330) * t322 + ((-t237 * t322 - t243 * t319) * t399 - t207 * t291 + t231 * t381) * qJD(5)) * MDP(27) + (t207 * t257 + t215 * t261 + t244 * t221 + (-(qJD(5) * t243 + t210) * t399 + t386 - (qJD(5) * t223 + t205) * t291 + t231 * qJD(5) * t290) * t322 + (-(-qJD(5) * t237 + t216) * t399 - (-qJD(5) * t225 + t213) * t291 + t330) * t319) * MDP(28) + (-t257 * MDP(13) - t258 * MDP(14)) * t315; t327 + t326 * t397 + (-t218 * t308 + t220 * t310 - t240 * t242 + t241 * t342 + t245 * t373) * MDP(21) + (t255 * t315 + (-t278 * t368 - t315 * t367) * pkin(2) + t333) * MDP(16) + (t256 * t315 + (-t315 * t356 - t368 * t369) * pkin(2) + t332) * MDP(17) + (-t249 * t310 - t250 * t308 + (-t245 + t342) * t369 + (t241 + t373) * t278) * MDP(18) + (t242 * t278 + t315 * t342 + t334) * MDP(19) + (t242 * t369 - t315 * t373 + t331) * MDP(20) + (t308 * t221 + t348 * t400 + t371 * t261 + (t307 * t249 - t343 * t399 - t350) * t319 + t353) * MDP(28) + (-t307 * t247 + t308 * t222 + t371 * t259 + (t319 * t348 + t322 * t343) * t399 + t339) * MDP(27) - t321 * MDP(4) * t376 + (MDP(9) * t321 * t326 + MDP(10) * t376) * pkin(1); t327 + (-t253 * t315 + t332) * MDP(17) + (-pkin(3) * t220 - qJ(4) * t218 - t240 * t251 - t241 * t254 - t245 * t363) * MDP(21) + (t333 + t384) * MDP(16) + (pkin(3) * t249 - qJ(4) * t250 + (-t245 - t254) * t369 + (t241 - t363) * t278) * MDP(18) + (t251 * t278 + t334 - t384) * MDP(19) + (t251 * t369 + t315 * t363 + t331) * MDP(20) + (qJ(4) * t221 + t347 * t400 + t364 * t261 + (t234 * t399 - t350 - t375) * t319 + t353) * MDP(28) + (t322 * t375 + qJ(4) * t222 + (-t234 * t322 + t319 * t347) * t399 + t364 * t259 + t339) * MDP(27); t328 * MDP(18) - MDP(19) * t382 + (-t315 ^ 2 - t396) * MDP(20) + (t245 * t315 + t334) * MDP(21) + (-t259 * t315 - t247) * MDP(27) + (-t261 * t315 + t385) * MDP(28) + (-MDP(27) * t352 - MDP(28) * t400) * t399; t261 * t259 * MDP(22) + (-t259 ^ 2 + t261 ^ 2) * MDP(23) + (t374 + t401) * MDP(24) + (t248 + t349) * MDP(25) - t249 * MDP(26) + (-t205 * t319 + t207 * t399 - t231 * t261 + t212) * MDP(27) + (-t205 * t322 + t206 * t399 - t213 * t319 + t231 * t259) * MDP(28) + (-MDP(24) * t380 - MDP(25) * t261 - MDP(27) * t207 - MDP(28) * t206) * qJD(5);];
tauc = t1;
