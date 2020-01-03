% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:38
% EndTime: 2019-12-31 19:24:44
% DurationCPUTime: 3.07s
% Computational Cost: add. (2781->374), mult. (8523->497), div. (0->0), fcn. (6005->6), ass. (0->155)
t433 = MDP(17) + MDP(20);
t352 = sin(pkin(5));
t355 = sin(qJ(2));
t356 = cos(qJ(2));
t331 = -qJ(3) * t352 * t355 - pkin(2) * t356 - pkin(1);
t354 = cos(pkin(5));
t385 = qJ(3) * t354 + pkin(7);
t332 = t385 * t355;
t353 = cos(pkin(8));
t432 = t353 * (t331 * t352 - t332 * t354);
t416 = t352 * t356;
t369 = pkin(2) * t355 - qJ(3) * t416;
t322 = t369 * qJD(1);
t375 = qJD(1) * t385;
t323 = t355 * t375;
t324 = t356 * t375;
t351 = sin(pkin(8));
t420 = t351 * t354;
t421 = t351 * t352;
t253 = t322 * t421 - t353 * t323 - t324 * t420;
t400 = qJD(3) * t352;
t405 = qJD(1) * t355;
t409 = qJ(4) * t352 * t405 - qJD(4) * t354 - t353 * t400 + t253;
t430 = t352 ^ 2 + t354 ^ 2;
t399 = qJD(3) * t355;
t290 = qJD(2) * t369 - t352 * t399;
t374 = qJD(2) * t385;
t413 = t354 * t356;
t291 = qJD(3) * t413 - t355 * t374;
t292 = -t354 * t399 - t356 * t374;
t243 = t290 * t421 + t353 * t291 + t292 * t420;
t401 = qJD(2) * t355;
t233 = -t352 * (qJ(4) * t401 - qJD(4) * t356) - t243;
t404 = qJD(1) * t356;
t390 = t354 * t404;
t403 = qJD(2) * t352;
t330 = t390 + t403;
t298 = pkin(7) * t404 + qJ(3) * t330;
t311 = qJD(2) * pkin(2) - t323;
t312 = t331 * qJD(1);
t246 = -t351 * t298 + t353 * (t311 * t354 + t312 * t352);
t429 = MDP(5) * (t355 ^ 2 - t356 ^ 2) + pkin(1) * (MDP(10) * t356 + MDP(9) * t355) - t355 * t356 * MDP(4);
t397 = t354 * qJD(2);
t328 = t352 * t404 - t397;
t325 = t328 ^ 2;
t426 = pkin(3) + qJ(5);
t277 = t290 * qJD(1);
t425 = t277 * t353;
t391 = t351 * t405;
t402 = qJD(2) * t353;
t279 = -t352 * t402 - t353 * t390 + t391;
t316 = t351 * t413 + t353 * t355;
t280 = qJD(1) * t316 + t351 * t403;
t424 = t279 * t280;
t423 = t290 * t353;
t422 = t322 * t353;
t419 = t351 * t355;
t418 = t351 * t356;
t417 = t352 * t353;
t415 = t353 * t354;
t414 = t354 * t355;
t380 = t354 * t391;
t304 = t353 * t404 - t380;
t308 = t351 * t323;
t376 = t324 * t415 - t308;
t388 = t426 * t355;
t412 = -pkin(4) * t304 - (-qJD(1) * t388 - t422) * t352 - t376 - qJD(5) * t354 + t351 * t400;
t368 = t353 * t414 + t418;
t303 = t368 * qJD(1);
t262 = t354 * t322 + t324 * t352;
t367 = -qJ(4) * t304 + t262;
t411 = t303 * t426 - (-qJD(4) * t351 - qJD(5) * t353) * t352 + t367;
t410 = -pkin(4) * t303 + t409;
t333 = t385 * t356;
t319 = t351 * t333;
t408 = pkin(3) * t416 + t319;
t318 = pkin(2) * t420 + qJ(3) * t417;
t278 = t292 * qJD(1);
t255 = t354 * t277 - t278 * t352;
t406 = MDP(14) * t255;
t247 = t353 * t298 + t311 * t420 + t312 * t421;
t398 = t247 * MDP(14);
t396 = qJD(1) * qJD(2);
t395 = MDP(11) - MDP(16);
t394 = MDP(13) + MDP(15);
t268 = qJD(1) * t291 + qJD(2) * t400;
t231 = t353 * t268 + t277 * t421 + t278 * t420;
t257 = t331 * t421 - t332 * t420 + t353 * t333;
t392 = -pkin(2) * t353 - pkin(3);
t389 = t356 * t402;
t387 = t355 * t396;
t386 = -qJ(4) * t351 - pkin(2);
t384 = MDP(21) + t395;
t383 = MDP(12) - t433;
t382 = MDP(19) + t394;
t258 = t354 * t290 - t292 * t352;
t265 = t354 * t331 + t332 * t352;
t342 = t352 * t387;
t379 = qJD(2) * t388;
t263 = t351 * t268;
t378 = -t278 * t415 + t263;
t282 = t351 * t291;
t377 = -t292 * t415 + t282;
t299 = -qJ(4) * t354 - t318;
t259 = -t311 * t352 + t354 * t312 + qJD(3);
t225 = -qJ(4) * t342 + t328 * qJD(4) - t231;
t360 = qJD(4) - t246;
t238 = pkin(3) * t328 + t360;
t370 = -t246 * MDP(14) + t238 * MDP(18);
t239 = qJ(4) * t328 - t247;
t334 = qJD(2) * t380;
t294 = qJD(1) * t389 - t334;
t366 = -qJ(4) * t316 + t265;
t251 = qJ(4) * t416 - t257;
t364 = -pkin(3) * t387 - t425;
t363 = -qJ(4) * t280 + t259;
t305 = t368 * qJD(2);
t293 = qJD(1) * t305;
t219 = -pkin(4) * t293 - t225;
t306 = -t397 * t419 + t389;
t361 = -qJ(4) * t306 - qJD(4) * t316 + t258;
t359 = pkin(4) * t294 + qJD(5) * t328 + t378;
t221 = pkin(3) * t293 - qJ(4) * t294 - qJD(4) * t280 + t255;
t217 = t293 * qJ(5) + t279 * qJD(5) + t221;
t344 = qJ(3) * t421;
t317 = pkin(2) * t415 - t344;
t315 = -t353 * t413 + t419;
t301 = (-pkin(3) * t353 + t386) * t352;
t300 = t354 * t392 + t344;
t275 = (-t353 * t426 + t386) * t352;
t274 = pkin(4) * t417 - t299;
t269 = pkin(4) * t421 + t344 + (-qJ(5) + t392) * t354;
t256 = -t319 + t432;
t254 = t408 - t432;
t252 = t308 + (t322 * t352 - t324 * t354) * t353;
t250 = pkin(3) * t315 + t366;
t249 = (-pkin(3) * t405 - t422) * t352 + t376;
t245 = pkin(3) * t303 + t367;
t244 = -pkin(4) * t315 - t251;
t242 = -t282 + (t290 * t352 + t292 * t354) * t353;
t241 = t315 * t426 + t366;
t240 = t332 * t415 + pkin(4) * t316 + (qJ(5) * t356 - t331 * t353) * t352 + t408;
t236 = (-pkin(3) * t401 - t423) * t352 + t377;
t232 = pkin(3) * t279 + t363;
t230 = -t263 + (t277 * t352 + t278 * t354) * t353;
t229 = pkin(3) * t305 + t361;
t228 = t352 * t364 + t378;
t227 = -pkin(4) * t279 + qJD(5) - t239;
t226 = -pkin(4) * t305 - t233;
t224 = pkin(4) * t306 + (qJD(5) * t356 - t379 - t423) * t352 + t377;
t223 = t279 * t426 + t363;
t222 = pkin(4) * t280 + t328 * t426 + t360;
t220 = qJD(5) * t315 + t305 * t426 + t361;
t218 = (-qJD(1) * t379 - t425) * t352 + t359;
t1 = [(-t242 * t328 + t255 * t315 + t258 * t279 + t259 * t305 + t265 * t293) * MDP(11) + (t243 * t328 + t255 * t316 + t258 * t280 + t259 * t306 + t265 * t294) * MDP(12) + (-t230 * t316 - t231 * t315 - t242 * t280 - t243 * t279 - t246 * t306 - t247 * t305 - t256 * t294 - t257 * t293) * MDP(13) + (t230 * t256 + t231 * t257 + t242 * t246 + t243 * t247 + t255 * t265 + t258 * t259) * MDP(14) + (t225 * t315 + t228 * t316 + t233 * t279 + t236 * t280 + t238 * t306 + t239 * t305 + t251 * t293 + t254 * t294) * MDP(15) + (-t221 * t315 - t229 * t279 - t232 * t305 - t236 * t328 - t250 * t293) * MDP(16) + (-t221 * t316 - t229 * t280 - t232 * t306 + t233 * t328 - t250 * t294) * MDP(17) + (t221 * t250 + t225 * t251 + t228 * t254 + t229 * t232 + t233 * t239 + t236 * t238) * MDP(18) + (t218 * t316 - t219 * t315 + t222 * t306 + t224 * t280 - t226 * t279 - t227 * t305 + t240 * t294 - t244 * t293) * MDP(19) + (-t217 * t316 - t220 * t280 - t223 * t306 - t226 * t328 - t241 * t294) * MDP(20) + (t217 * t315 + t220 * t279 + t223 * t305 + t224 * t328 + t241 * t293) * MDP(21) + (t217 * t241 + t218 * t240 + t219 * t244 + t220 * t223 + t222 * t224 + t226 * t227) * MDP(22) + (t356 * MDP(6) - t355 * MDP(7) + (MDP(10) * t355 - MDP(9) * t356) * pkin(7)) * qJD(2) ^ 2 - 0.2e1 * t429 * t396 + ((-MDP(11) * t230 + MDP(12) * t231 - MDP(16) * t228 + MDP(17) * t225 - MDP(20) * t219 + MDP(21) * t218) * t356 + ((qJD(1) * t256 + t246) * MDP(11) + (-qJD(1) * t257 - t247) * MDP(12) + (qJD(1) * t254 + t238) * MDP(16) + (-qJD(1) * t251 - t239) * MDP(17) + (qJD(1) * t244 + t227) * MDP(20) + (-qJD(1) * t240 - t222) * MDP(21)) * t401) * t352; (t230 * t354 - t259 * t303) * MDP(11) + (-t231 * t354 - t259 * t304) * MDP(12) + (t246 * t304 + t247 * t303 - t293 * t318 - t294 * t317) * MDP(13) + (t230 * t317 + t231 * t318 - t246 * t252 - t247 * t253 - t259 * t262) * MDP(14) + (-t238 * t304 - t239 * t303 + t293 * t299 + t294 * t300) * MDP(15) + (t228 * t354 + t232 * t303 - t293 * t301) * MDP(16) + (-t225 * t354 + t232 * t304 - t294 * t301) * MDP(17) + (t221 * t301 + t225 * t299 + t228 * t300 - t232 * t245 - t238 * t249 + t239 * t409) * MDP(18) + (-t222 * t304 + t227 * t303 + t269 * t294 - t274 * t293) * MDP(19) + (t219 * t354 + t223 * t304 - t275 * t294) * MDP(20) + (-t218 * t354 - t223 * t303 + t275 * t293) * MDP(21) + (t217 * t275 + t218 * t269 + t219 * t274 + t222 * t412 - t223 * t411 - t227 * t410) * MDP(22) + (-t262 * MDP(12) + t252 * MDP(13) - t249 * MDP(15) + t245 * MDP(17) + MDP(19) * t412 + MDP(20) * t411) * t280 + (t252 * MDP(11) - t253 * MDP(12) + t249 * MDP(16) + MDP(17) * t409 + MDP(20) * t410 + MDP(21) * t412) * t328 + (-t262 * MDP(11) + t253 * MDP(13) + MDP(15) * t409 + t245 * MDP(16) + MDP(19) * t410 - MDP(21) * t411) * t279 + t429 * qJD(1) ^ 2 + ((-MDP(11) * t293 - MDP(12) * t294 - t406) * pkin(2) + (-t255 * MDP(11) + t231 * MDP(13) - t225 * MDP(15) + t221 * MDP(16) + t219 * MDP(19) - t217 * MDP(21) + (MDP(12) * t328 - MDP(13) * t279 + t398) * qJD(3)) * t353 + (t255 * MDP(12) - t230 * MDP(13) + t228 * MDP(15) - t221 * MDP(17) + t218 * MDP(19) - t217 * MDP(20) + (MDP(16) * t279 + MDP(17) * t280 - MDP(18) * t232) * qJD(4) + (t280 * t394 + t328 * t395 + t370) * qJD(3)) * t351 + ((qJD(2) * t317 - t246) * MDP(11) + (-qJD(2) * t318 + t247) * MDP(12) + (qJD(2) * t300 - t238) * MDP(16) + (-qJD(2) * t299 + t239) * MDP(17) + (qJD(2) * t274 - t227) * MDP(20) + (-qJD(2) * t269 + t222) * MDP(21)) * t405) * t352; t406 + t221 * MDP(18) + t217 * MDP(22) - t383 * t334 + (t239 * MDP(18) - t227 * MDP(22) + t279 * t382 - t328 * t383 - t398) * (t330 * t353 - t391 * t430) + (-t222 * MDP(22) - t280 * t382 - t328 * t384 - t370) * (t353 * t405 * t430 + t330 * t351) + (t384 * t418 + (t356 * t383 + t384 * t414) * t353) * t396; (t232 * t280 - t239 * t328 + t378) * MDP(18) + (t223 * t280 + t227 * t328 + t359) * MDP(22) + (t364 * MDP(18) + (-qJ(5) * t387 + t364) * MDP(22)) * t352 + (-MDP(16) + MDP(21)) * (-t342 + t424) + (MDP(15) + MDP(19)) * (-t279 * t328 + t294) + t433 * (-t280 ^ 2 - t325); (-t280 * t328 - t368 * t396) * MDP(19) + (t342 + t424) * MDP(20) + (-t279 ^ 2 - t325) * MDP(21) + (-t222 * t328 - t223 * t279 + t219) * MDP(22);];
tauc = t1;
