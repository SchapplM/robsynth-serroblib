% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:31
% EndTime: 2019-12-31 21:09:38
% DurationCPUTime: 3.46s
% Computational Cost: add. (2103->399), mult. (5007->516), div. (0->0), fcn. (2860->4), ass. (0->157)
t318 = sin(qJ(3));
t319 = sin(qJ(2));
t320 = cos(qJ(3));
t375 = qJD(3) * t320;
t359 = t319 * t375;
t321 = cos(qJ(2));
t377 = qJD(2) * t321;
t423 = t318 * t377 + t359;
t283 = -pkin(2) * t321 - pkin(7) * t319 - pkin(1);
t265 = t283 * qJD(1);
t380 = qJD(1) * t321;
t309 = pkin(6) * t380;
t287 = qJD(2) * pkin(7) + t309;
t242 = -t320 * t265 + t287 * t318;
t379 = qJD(2) * t318;
t381 = qJD(1) * t319;
t275 = t320 * t381 + t379;
t335 = pkin(4) * t275 + t242;
t371 = qJD(4) + t335;
t369 = qJD(1) * qJD(2);
t422 = -0.2e1 * t369;
t315 = t319 ^ 2;
t421 = (-t321 ^ 2 + t315) * MDP(5);
t301 = -qJD(3) + t380;
t290 = qJD(4) * t301;
t306 = t319 * t369;
t299 = qJ(4) * t306;
t368 = qJD(2) * qJD(3);
t357 = t318 * t368;
t252 = t423 * qJD(1) + t357;
t343 = pkin(2) * t319 - pkin(7) * t321;
t281 = t343 * qJD(2);
t266 = qJD(1) * t281;
t348 = pkin(6) * t306;
t376 = qJD(3) * t318;
t347 = t265 * t375 + t318 * t266 - t287 * t376 - t320 * t348;
t329 = -pkin(4) * t252 + t347;
t217 = -t290 + t299 + t329;
t413 = pkin(3) + qJ(5);
t224 = t413 * t301 + t371;
t420 = -t224 * t301 + t217;
t346 = -t265 * t376 + t320 * t266 - t287 * t375 + t318 * t348;
t223 = -pkin(3) * t306 - t346;
t243 = t318 * t265 + t320 * t287;
t236 = qJ(4) * t301 - t243;
t419 = -t236 * t301 + t223;
t417 = qJD(4) * t318 + t309 + (t318 * t380 - t376) * pkin(3);
t416 = MDP(16) - MDP(19);
t298 = t301 ^ 2;
t415 = pkin(4) + pkin(7);
t373 = t320 * qJD(2);
t273 = t318 * t381 - t373;
t414 = pkin(4) * t273;
t412 = qJ(4) * t252;
t411 = qJ(4) * t273;
t410 = qJ(4) * t320;
t358 = t321 * t369;
t361 = t319 * t376;
t251 = qJD(1) * t361 + (-t358 - t368) * t320;
t331 = pkin(6) * t358 + qJ(4) * t251 - qJD(4) * t275;
t220 = pkin(3) * t252 + t331;
t409 = t220 * t318;
t408 = t220 * t320;
t405 = t251 * t318;
t404 = t273 * t275;
t403 = t273 * t301;
t402 = t275 * t301;
t278 = t343 * qJD(1);
t401 = t278 * t320;
t286 = -qJD(2) * pkin(2) + pkin(6) * t381;
t400 = t286 * t318;
t399 = t286 * t320;
t398 = t301 * t320;
t397 = t318 * t319;
t396 = t318 * t321;
t395 = t319 * t320;
t322 = qJD(2) ^ 2;
t394 = t319 * t322;
t393 = t320 * t321;
t392 = t321 * t322;
t323 = qJD(1) ^ 2;
t391 = t321 * t323;
t289 = t415 * t320;
t363 = -pkin(6) * t318 - pkin(3);
t325 = pkin(4) * t393 + (-qJ(5) + t363) * t319;
t390 = -t325 * qJD(1) + qJD(3) * t289 + t401;
t262 = t318 * t278;
t355 = pkin(6) * t320 - qJ(4);
t326 = -pkin(4) * t396 - t355 * t319;
t389 = t326 * qJD(1) + t415 * t376 + t262;
t338 = qJ(5) * t318 - t410;
t328 = t321 * t338;
t388 = qJD(1) * t328 - t338 * qJD(3) + qJD(5) * t320 + t417;
t387 = qJ(4) * t375 - t380 * t410 + t417;
t386 = t318 * t281 + t283 * t375;
t384 = pkin(3) * t397 + t319 * pkin(6);
t305 = pkin(6) * t393;
t383 = t318 * t283 + t305;
t378 = qJD(2) * t319;
t374 = t319 * MDP(15);
t372 = -qJD(4) - t242;
t231 = t243 - t414;
t370 = -qJD(5) - t231;
t367 = pkin(7) * t301 * t318;
t366 = pkin(7) * t398;
t304 = pkin(6) * t396;
t365 = pkin(7) * t378;
t364 = pkin(7) * t373;
t360 = t321 * t376;
t356 = -qJ(4) * t318 - pkin(2);
t354 = pkin(1) * t422;
t352 = t283 * t320 - t304;
t351 = t273 + t373;
t350 = -t275 + t379;
t345 = t423 * pkin(3) + pkin(6) * t377 + qJ(4) * t361;
t344 = t363 * t319;
t249 = qJ(4) * t321 - t383;
t342 = -qJD(3) * t305 + t281 * t320 - t283 * t376;
t341 = -qJD(4) * t321 + t386;
t340 = pkin(6) * (-t301 + t380);
t339 = -t290 + t347;
t235 = pkin(3) * t301 - t372;
t337 = t235 * t320 + t236 * t318;
t336 = qJD(1) * t315 - t301 * t321;
t334 = -qJ(4) * t275 + t286;
t215 = qJD(5) * t273 + t413 * t252 + t331;
t227 = t413 * t273 + t334;
t333 = t215 * t318 + t227 * t375;
t332 = -t215 * t320 + t227 * t376;
t330 = pkin(4) * t251 + t346;
t232 = -t251 - t403;
t324 = -t413 * t306 - t330;
t314 = t321 * pkin(3);
t297 = 0.2e1 * t299;
t288 = t415 * t318;
t282 = -pkin(3) * t320 + t356;
t264 = -t413 * t320 + t356;
t256 = -qJ(4) * t395 + t384;
t250 = t314 - t352;
t248 = t338 * t319 + t384;
t246 = pkin(3) * t275 + t411;
t245 = qJD(1) * t344 - t401;
t244 = t355 * t381 - t262;
t240 = -pkin(4) * t397 - t249;
t239 = qJ(5) * t321 + t304 + t314 + (pkin(4) * t319 - t283) * t320;
t238 = pkin(3) * t273 + t334;
t233 = t413 * t275 + t411;
t229 = (-qJ(4) * t377 - qJD(4) * t319) * t320 + t345;
t228 = qJD(2) * t344 - t342;
t226 = qJD(5) - t236 - t414;
t225 = -qJ(4) * t378 + (t319 * t373 + t360) * pkin(6) - t341;
t222 = qJD(2) * t328 + (qJD(5) * t318 + (qJ(5) * qJD(3) - qJD(4)) * t320) * t319 + t345;
t221 = (-pkin(4) * t395 - t304) * qJD(3) + t326 * qJD(2) + t341;
t219 = -t299 - t339;
t218 = -pkin(4) * t361 + t325 * qJD(2) + qJD(5) * t321 - t342;
t216 = qJD(5) * t301 + t324;
t1 = [0.2e1 * t321 * MDP(4) * t306 + t421 * t422 + MDP(6) * t392 - MDP(7) * t394 + (-pkin(6) * t392 + t319 * t354) * MDP(9) + (pkin(6) * t394 + t321 * t354) * MDP(10) + (-t251 * t395 + (t321 * t373 - t361) * t275) * MDP(11) + ((-t273 * t320 - t275 * t318) * t377 + (t405 - t252 * t320 + (t273 * t318 - t275 * t320) * qJD(3)) * t319) * MDP(12) + (t301 * t361 + t251 * t321 + (t275 * t319 + t336 * t320) * qJD(2)) * MDP(13) + (t301 * t359 + t252 * t321 + (-t273 * t319 - t336 * t318) * qJD(2)) * MDP(14) + (-t301 - t380) * qJD(2) * t374 + (-t342 * t301 - t346 * t321 + (pkin(6) * t252 + t286 * t375) * t319 + ((pkin(6) * t273 + t400) * t321 + (t352 * qJD(1) + t318 * t340 - t242) * t319) * qJD(2)) * MDP(16) + ((-pkin(6) * t360 + t386) * t301 + t347 * t321 + (-pkin(6) * t251 - t286 * t376) * t319 + ((pkin(6) * t275 + t399) * t321 + (-t383 * qJD(1) + t320 * t340 - t243) * t319) * qJD(2)) * MDP(17) + (t225 * t273 + t228 * t275 + t249 * t252 - t250 * t251 + t337 * t377 + (t219 * t318 + t223 * t320 + (-t235 * t318 + t236 * t320) * qJD(3)) * t319) * MDP(18) + (-t228 * t301 - t229 * t273 - t252 * t256 + (-t238 * t379 - t223) * t321 + (-t238 * t375 - t409 + (qJD(1) * t250 + t235) * qJD(2)) * t319) * MDP(19) + (t225 * t301 - t229 * t275 + t251 * t256 + (-t238 * t373 + t219) * t321 + (t238 * t376 - t408 + (-qJD(1) * t249 - t236) * qJD(2)) * t319) * MDP(20) + (t219 * t249 + t220 * t256 + t223 * t250 + t225 * t236 + t228 * t235 + t229 * t238) * MDP(21) + (t218 * t275 - t221 * t273 - t239 * t251 - t240 * t252 + (t224 * t320 - t226 * t318) * t377 + (t216 * t320 - t217 * t318 + (-t224 * t318 - t226 * t320) * qJD(3)) * t319) * MDP(22) + (-t221 * t301 - t222 * t275 + t248 * t251 + (-t227 * t373 - t217) * t321 + ((qJD(1) * t240 + t226) * qJD(2) + t332) * t319) * MDP(23) + (t218 * t301 + t222 * t273 + t248 * t252 + (t227 * t379 + t216) * t321 + ((-qJD(1) * t239 - t224) * qJD(2) + t333) * t319) * MDP(24) + (t215 * t248 + t216 * t239 + t217 * t240 + t218 * t224 + t221 * t226 + t222 * t227) * MDP(25); -t319 * MDP(4) * t391 + t323 * t421 + (-t275 * t398 - t405) * MDP(11) + ((-t251 + t403) * t320 + (-t252 + t402) * t318) * MDP(12) + (-t301 * t375 + (t301 * t393 + t350 * t319) * qJD(1)) * MDP(13) + (t301 * t376 + (-t301 * t396 + t351 * t319) * qJD(1)) * MDP(14) + t301 * qJD(1) * t374 + (t278 * t398 - pkin(2) * t252 + (t366 + t400) * qJD(3) + (t242 * t319 + (-t286 * t321 - t365) * t318 + (t301 * t397 - t351 * t321) * pkin(6)) * qJD(1)) * MDP(16) + (pkin(2) * t251 - t262 * t301 + (-t367 + t399) * qJD(3) + (-t286 * t393 + (t243 - t364) * t319 + (t301 * t395 + t350 * t321) * pkin(6)) * qJD(1)) * MDP(17) + (-t244 * t273 - t245 * t275 + (-t219 - t301 * t235 + (qJD(3) * t275 - t252) * pkin(7)) * t320 + ((qJD(3) * t273 - t251) * pkin(7) + t419) * t318) * MDP(18) + (t408 + t245 * t301 - t252 * t282 + t387 * t273 + (-t238 * t318 - t366) * qJD(3) + (-t235 * t319 + (t238 * t321 + t365) * t318) * qJD(1)) * MDP(19) + (-t409 - t244 * t301 + t251 * t282 + t387 * t275 + (-t238 * t320 + t367) * qJD(3) + (t238 * t393 + (t236 + t364) * t319) * qJD(1)) * MDP(20) + (t220 * t282 - t235 * t245 - t236 * t244 - t387 * t238 + (t337 * qJD(3) - t219 * t320 + t223 * t318) * pkin(7)) * MDP(21) + (-t251 * t288 - t252 * t289 + t390 * t275 + t389 * t273 + t420 * t320 + (t226 * t301 + t216) * t318) * MDP(22) + (t251 * t264 + t389 * t301 + t388 * t275 + (t227 * t393 + (qJD(2) * t289 - t226) * t319) * qJD(1) - t333) * MDP(23) + (t252 * t264 + t390 * t301 - t388 * t273 + (-t227 * t396 + (-qJD(2) * t288 + t224) * t319) * qJD(1) + t332) * MDP(24) + (t215 * t264 + t216 * t288 + t217 * t289 + t390 * t224 - t389 * t226 - t388 * t227) * MDP(25) + (t323 * t319 * MDP(9) + MDP(10) * t391) * pkin(1); t232 * MDP(13) - MDP(14) * t357 - t347 * MDP(17) + (pkin(3) * t251 - t412) * MDP(18) + (t297 + t339) * MDP(20) + (-pkin(3) * t223 - qJ(4) * t219 - t235 * t243 + t372 * t236 - t238 * t246) * MDP(21) + (t251 * t413 - t412) * MDP(22) + (-0.2e1 * t290 + t297 + t329) * MDP(23) + t330 * MDP(24) + (qJ(4) * t217 - t216 * t413 + t370 * t224 + t371 * t226 - t227 * t233) * MDP(25) + (t242 * MDP(17) + t372 * MDP(20) - t335 * MDP(23) + (-0.2e1 * qJD(5) - t231) * MDP(24) - t416 * t243) * t301 + (-t301 * MDP(14) - t286 * MDP(16) + (-t236 - t243) * MDP(18) + t238 * MDP(19) + t246 * MDP(20) + (t226 + t370) * MDP(22) + t233 * MDP(23) - t227 * MDP(24) + MDP(12) * t275) * t275 + (-MDP(14) * t359 + (-MDP(14) * t396 + (-0.2e1 * pkin(3) * MDP(19) + 0.2e1 * t413 * MDP(24) + MDP(15)) * t319) * qJD(2)) * qJD(1) + (t275 * MDP(11) + t286 * MDP(17) + (t235 + t372) * MDP(18) + t246 * MDP(19) - t238 * MDP(20) + (t224 - t371) * MDP(22) - t227 * MDP(23) - t233 * MDP(24) - MDP(12) * t273) * t273 + t416 * t346; (t238 * t275 + t419) * MDP(21) + (t227 * t275 + (qJD(5) + t226) * t301 + t324) * MDP(25) + (-MDP(19) + MDP(24)) * (-t306 + t404) + (MDP(18) + MDP(22)) * t232 + (MDP(20) + MDP(23)) * (-t275 ^ 2 - t298); (-t252 - t402) * MDP(22) + (t306 + t404) * MDP(23) + (-t273 ^ 2 - t298) * MDP(24) + (-t227 * t273 + t420) * MDP(25);];
tauc = t1;
