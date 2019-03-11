% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:41
% EndTime: 2019-03-08 21:26:49
% DurationCPUTime: 3.70s
% Computational Cost: add. (3251->335), mult. (8333->481), div. (0->0), fcn. (6283->10), ass. (0->157)
t346 = sin(qJ(2));
t341 = sin(pkin(6));
t392 = qJD(1) * t341;
t380 = t346 * t392;
t345 = sin(qJ(3));
t388 = qJD(3) * t345;
t435 = pkin(3) * t388 - t380;
t348 = cos(qJ(3));
t422 = -qJ(4) - pkin(8);
t372 = qJD(3) * t422;
t310 = qJD(4) * t348 + t345 * t372;
t311 = -qJD(4) * t345 + t348 * t372;
t340 = sin(pkin(11));
t342 = cos(pkin(11));
t268 = t310 * t342 + t311 * t340;
t410 = t342 * t348;
t321 = t340 * t345 - t410;
t349 = cos(qJ(2));
t379 = t349 * t392;
t293 = t321 * t379;
t400 = t268 + t293;
t322 = t340 * t348 + t342 * t345;
t314 = t322 * qJD(3);
t317 = t321 * qJD(3);
t434 = pkin(4) * t314 + pkin(9) * t317 + t435;
t315 = t322 * qJD(2);
t347 = cos(qJ(5));
t384 = qJD(2) * qJD(3);
t373 = t348 * t384;
t374 = t345 * t384;
t364 = -t340 * t374 + t342 * t373;
t385 = t347 * qJD(3);
t344 = sin(qJ(5));
t387 = qJD(5) * t344;
t260 = -qJD(5) * t385 + t315 * t387 - t347 * t364;
t299 = qJD(3) * t344 + t315 * t347;
t307 = qJD(2) * t314;
t343 = cos(pkin(6));
t391 = qJD(1) * t343;
t331 = t348 * t391;
t324 = qJD(2) * pkin(8) + t380;
t369 = qJ(4) * qJD(2) + t324;
t290 = -t369 * t345 + t331;
t285 = qJD(3) * pkin(3) + t290;
t378 = t345 * t391;
t291 = t369 * t348 + t378;
t411 = t342 * t291;
t244 = t340 * t285 + t411;
t241 = qJD(3) * pkin(9) + t244;
t381 = -pkin(3) * t348 - pkin(2);
t308 = t381 * qJD(2) + qJD(4) - t379;
t389 = qJD(2) * t345;
t313 = qJD(2) * t410 - t340 * t389;
t255 = -pkin(4) * t313 - pkin(9) * t315 + t308;
t233 = t241 * t347 + t255 * t344;
t367 = qJD(4) + t379;
t259 = (-t324 * t345 + t331) * qJD(3) + (-qJ(4) * t388 + t367 * t348) * qJD(2);
t426 = (-t324 * t348 - t378) * qJD(3) + (-qJD(3) * t348 * qJ(4) - t367 * t345) * qJD(2);
t236 = t342 * t259 + t426 * t340;
t390 = qJD(2) * t341;
t377 = t346 * t390;
t312 = pkin(3) * t374 + qJD(1) * t377;
t252 = t307 * pkin(4) - t364 * pkin(9) + t312;
t250 = t347 * t252;
t353 = -t233 * qJD(5) - t236 * t344 + t250;
t220 = pkin(5) * t307 + qJ(6) * t260 - qJD(6) * t299 + t353;
t297 = t315 * t344 - t385;
t227 = -qJ(6) * t297 + t233;
t309 = qJD(5) - t313;
t433 = t309 * t227 + t220;
t261 = t299 * qJD(5) + t344 * t364;
t386 = qJD(5) * t347;
t358 = t347 * t236 - t241 * t387 + t344 * t252 + t255 * t386;
t221 = -qJ(6) * t261 - qJD(6) * t297 + t358;
t232 = -t241 * t344 + t347 * t255;
t226 = -qJ(6) * t299 + t232;
t224 = pkin(5) * t309 + t226;
t432 = -t309 * t224 + t221;
t431 = MDP(5) * t345;
t428 = (t345 ^ 2 - t348 ^ 2) * MDP(6);
t427 = -t293 * t344 + t347 * t434;
t280 = pkin(4) * t321 - pkin(9) * t322 + t381;
t425 = t280 * t386 + t434 * t344 + t347 * t400;
t424 = MDP(10) * t345 + MDP(11) * t348;
t423 = t299 ^ 2;
t421 = qJD(2) * pkin(2);
t420 = t260 * t344;
t419 = t297 * t313;
t418 = t299 * t309;
t417 = t309 * t344;
t416 = t313 * t344;
t415 = t322 * t344;
t414 = t322 * t347;
t281 = t340 * t291;
t413 = t341 * t346;
t412 = t341 * t349;
t350 = qJD(3) ^ 2;
t409 = t345 * t350;
t326 = t422 * t345;
t327 = t422 * t348;
t295 = t326 * t340 - t327 * t342;
t286 = t347 * t295;
t303 = t347 * t307;
t408 = t348 * t350;
t334 = pkin(3) * t340 + pkin(9);
t407 = qJ(6) + t334;
t365 = qJ(6) * t317 - qJD(6) * t322;
t406 = pkin(5) * t314 - t268 * t344 + t365 * t347 + (-t286 + (qJ(6) * t322 - t280) * t344) * qJD(5) + t427;
t375 = t322 * t386;
t405 = -qJ(6) * t375 + (-qJD(5) * t295 + t365) * t344 + t425;
t404 = t224 - t226;
t248 = t290 * t342 - t281;
t269 = pkin(3) * t389 + pkin(4) * t315 - pkin(9) * t313;
t403 = t347 * t248 + t344 * t269;
t402 = -t344 * t261 - t297 * t386;
t401 = t340 * t310 - t342 * t311 - t322 * t379;
t399 = t309 * t416 + t303;
t398 = t344 * t280 + t286;
t371 = qJD(5) * t407;
t397 = qJ(6) * t416 + qJD(6) * t347 - t344 * t371 - t403;
t265 = t347 * t269;
t396 = -pkin(5) * t315 - t265 + (qJ(6) * t313 - t371) * t347 + (-qJD(6) + t248) * t344;
t335 = -pkin(3) * t342 - pkin(4);
t376 = t349 * t390;
t235 = t259 * t340 - t342 * t426;
t243 = t285 * t342 - t281;
t246 = t340 * t290 + t411;
t294 = -t342 * t326 - t327 * t340;
t370 = t309 * t347;
t366 = t235 * t322 - t295 * t307;
t225 = pkin(5) * t261 + t235;
t240 = -qJD(3) * pkin(4) - t243;
t318 = t343 * t345 + t348 * t413;
t362 = t343 * t348 - t345 * t413;
t275 = t342 * t318 + t340 * t362;
t257 = -t275 * t344 - t347 * t412;
t363 = -t275 * t347 + t344 * t412;
t361 = -t317 * t344 + t375;
t360 = -t317 * t347 - t322 * t387;
t357 = t309 * t240 - t334 * t307;
t354 = -0.2e1 * qJD(3) * t421;
t351 = qJD(2) ^ 2;
t320 = t407 * t347;
t319 = t407 * t344;
t296 = t297 ^ 2;
t289 = -t318 * qJD(3) - t345 * t376;
t288 = t362 * qJD(3) + t348 * t376;
t277 = t347 * t280;
t274 = t318 * t340 - t342 * t362;
t247 = t288 * t342 + t289 * t340;
t245 = t288 * t340 - t342 * t289;
t239 = -qJ(6) * t415 + t398;
t238 = pkin(5) * t297 + qJD(6) + t240;
t237 = pkin(5) * t321 - qJ(6) * t414 - t295 * t344 + t277;
t231 = t363 * qJD(5) - t247 * t344 + t347 * t377;
t230 = t257 * qJD(5) + t247 * t347 + t344 * t377;
t1 = [t289 * qJD(3) * MDP(10) - t288 * qJD(3) * MDP(11) + (t245 * t315 + t247 * t313 + t274 * t364 - t275 * t307) * MDP(12) + (t235 * t274 + t236 * t275 - t243 * t245 + t244 * t247) * MDP(13) + (t231 * t309 + t245 * t297 + t257 * t307 + t261 * t274) * MDP(19) + (-t230 * t309 + t245 * t299 - t260 * t274 + t307 * t363) * MDP(20) + (-t230 * t297 - t231 * t299 + t257 * t260 + t261 * t363) * MDP(21) + (t220 * t257 - t221 * t363 + t224 * t231 + t225 * t274 + t227 * t230 + t238 * t245) * MDP(22) + (-t312 * t349 * MDP(13) + (t308 * t346 * MDP(13) - t424 * t349 * qJD(3)) * qJD(2) + (-t349 * MDP(4) + (-MDP(10) * t348 + MDP(11) * t345 - MDP(3)) * t346) * t351) * t341; 0.2e1 * t373 * t431 - 0.2e1 * t384 * t428 + MDP(7) * t408 - MDP(8) * t409 + (-pkin(8) * t408 + t345 * t354) * MDP(10) + (pkin(8) * t409 + t348 * t354) * MDP(11) + (-t236 * t321 + t243 * t317 - t244 * t314 + t294 * t364 + t400 * t313 + t401 * t315 + t366) * MDP(12) + (t235 * t294 + t236 * t295 - t401 * t243 + t400 * t244 + t435 * t308 + t312 * t381) * MDP(13) + (-t260 * t414 + t360 * t299) * MDP(14) + (-(-t297 * t347 - t299 * t344) * t317 + (t420 - t261 * t347 + (t297 * t344 - t299 * t347) * qJD(5)) * t322) * MDP(15) + (-t260 * t321 + t299 * t314 + t322 * t303 + t360 * t309) * MDP(16) + (-t261 * t321 - t297 * t314 - t307 * t415 - t361 * t309) * MDP(17) + (t307 * t321 + t309 * t314) * MDP(18) + (t277 * t307 + (-t241 * t386 + t250) * t321 + t232 * t314 + t294 * t261 + t240 * t375 + (-t295 * t386 + t427) * t309 + t401 * t297 + ((-qJD(5) * t280 - t268) * t309 + (-qJD(5) * t255 - t236) * t321 - t240 * t317 + t366) * t344) * MDP(19) + (-t398 * t307 - t358 * t321 - t233 * t314 - t294 * t260 + t235 * t414 + (t295 * t387 - t425) * t309 + t401 * t299 + t360 * t240) * MDP(20) + (t237 * t260 - t239 * t261 - (-t224 * t347 - t227 * t344) * t317 - t406 * t299 - t405 * t297 + (-t220 * t347 - t221 * t344 + (t224 * t344 - t227 * t347) * qJD(5)) * t322) * MDP(21) + (t221 * t239 + t220 * t237 + t225 * (pkin(5) * t415 + t294) + (t361 * pkin(5) + t401) * t238 + t405 * t227 + t406 * t224) * MDP(22); ((t244 - t246) * t315 + (-t248 + t243) * t313 + (-t340 * t307 - t342 * t364) * pkin(3)) * MDP(12) + (t243 * t246 - t244 * t248 + (-t235 * t342 + t236 * t340 - t308 * t389) * pkin(3)) * MDP(13) + (t299 * t370 - t420) * MDP(14) + ((-t260 + t419) * t347 - t299 * t417 + t402) * MDP(15) + (-t299 * t315 + t307 * t344 + t309 * t370) * MDP(16) + (t297 * t315 - t309 * t387 + t399) * MDP(17) - t309 * t315 * MDP(18) + (-t232 * t315 - t235 * t347 - t246 * t297 + t335 * t261 + (-t334 * t386 - t265) * t309 + (t248 * t309 + t357) * t344) * MDP(19) + (t233 * t315 + t235 * t344 - t246 * t299 - t335 * t260 + (t334 * t387 + t403) * t309 + t357 * t347) * MDP(20) + (-t260 * t319 - t261 * t320 - t397 * t297 - t396 * t299 - t433 * t344 + t432 * t347) * MDP(21) + (t221 * t320 - t220 * t319 + t225 * (-pkin(5) * t347 + t335) + (pkin(5) * t417 - t246) * t238 + t397 * t227 + t396 * t224) * MDP(22) + t424 * qJD(2) * t421 + (-t348 * t431 + t428) * t351; -t313 ^ 2 * MDP(12) + (-t244 * t313 + t312) * MDP(13) + t399 * MDP(19) + t402 * MDP(21) + (-MDP(12) * t315 + MDP(13) * t243 - MDP(19) * t297 - MDP(20) * t299 - MDP(22) * t238) * t315 + (-qJD(5) * t309 * MDP(19) - t307 * MDP(20) + MDP(21) * t418 + t432 * MDP(22)) * t344 + ((t260 + t419) * MDP(21) + t433 * MDP(22) - t309 ^ 2 * MDP(20)) * t347; t299 * t297 * MDP(14) + (-t296 + t423) * MDP(15) + (t297 * t309 - t260) * MDP(16) + (-t261 + t418) * MDP(17) + t307 * MDP(18) + (t233 * t309 - t240 * t299 + t353) * MDP(19) + (t232 * t309 + t240 * t297 - t358) * MDP(20) + (pkin(5) * t260 - t404 * t297) * MDP(21) + (t404 * t227 + (-t238 * t299 + t220) * pkin(5)) * MDP(22); (-t296 - t423) * MDP(21) + (t224 * t299 + t227 * t297 + t225) * MDP(22);];
tauc  = t1;
