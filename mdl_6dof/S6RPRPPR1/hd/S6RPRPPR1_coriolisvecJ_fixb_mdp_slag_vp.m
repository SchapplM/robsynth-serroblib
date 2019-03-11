% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:41
% EndTime: 2019-03-09 02:39:49
% DurationCPUTime: 4.22s
% Computational Cost: add. (3330->343), mult. (8128->484), div. (0->0), fcn. (5855->10), ass. (0->155)
t363 = cos(qJ(3));
t423 = cos(pkin(10));
t388 = t423 * t363;
t343 = qJD(1) * t388;
t356 = sin(pkin(10));
t361 = sin(qJ(3));
t395 = t361 * qJD(1);
t323 = t356 * t395 - t343;
t320 = qJD(6) + t323;
t337 = t356 * t363 + t361 * t423;
t326 = t337 * qJD(1);
t355 = sin(pkin(11));
t358 = cos(pkin(11));
t310 = t358 * qJD(3) - t326 * t355;
t362 = cos(qJ(6));
t430 = t362 * t310;
t348 = sin(pkin(9)) * pkin(1) + pkin(7);
t407 = qJ(4) + t348;
t360 = sin(qJ(6));
t338 = t355 * t362 + t358 * t360;
t402 = t320 * t338;
t309 = qJD(3) * t355 + t326 * t358;
t269 = t309 * t362 + t310 * t360;
t429 = MDP(6) * (t361 ^ 2 - t363 ^ 2);
t428 = t363 * MDP(5);
t386 = t407 * qJD(1);
t312 = t363 * qJD(2) - t386 * t361;
t325 = t337 * qJD(3);
t317 = qJD(1) * t325;
t409 = t362 * t358;
t411 = t355 * t360;
t336 = -t409 + t411;
t403 = t320 * t336;
t427 = -t317 * t338 + t320 * t403;
t313 = qJD(2) * t361 + t363 * t386;
t426 = pkin(8) * t358;
t345 = pkin(3) * t356 + qJ(5);
t425 = pkin(8) + t345;
t424 = qJD(3) * pkin(3);
t392 = qJD(1) * qJD(4);
t296 = t312 * qJD(3) + t363 * t392;
t367 = -t313 * qJD(3) - t361 * t392;
t251 = t296 * t356 - t423 * t367;
t333 = t407 * t361;
t334 = t407 * t363;
t294 = t423 * t333 + t334 * t356;
t422 = t251 * t294;
t371 = -t356 * t361 + t388;
t421 = t251 * t371;
t267 = t309 * t360 - t430;
t420 = t267 * t326;
t419 = t269 * t326;
t393 = qJD(1) * qJD(3);
t390 = t361 * t393;
t318 = qJD(3) * t343 - t356 * t390;
t417 = t318 * t355;
t416 = t318 * t358;
t415 = t323 * t355;
t328 = t371 * qJD(3);
t414 = t328 * t355;
t413 = t337 * t355;
t412 = t337 * t358;
t302 = t356 * t313;
t364 = qJD(3) ^ 2;
t410 = t361 * t364;
t408 = t363 * t364;
t376 = -qJD(6) * t309 - t417;
t396 = qJD(6) * t362;
t404 = t310 * t396 + t318 * t409;
t236 = t360 * t376 + t404;
t406 = -t236 * t371 + t269 * t325;
t397 = qJD(6) * t337;
t245 = t328 * t338 + t396 * t412 - t397 * t411;
t289 = t338 * t337;
t405 = -t245 * t320 - t289 * t317;
t252 = t423 * t296 + t356 * t367;
t249 = qJD(3) * qJD(5) + t252;
t344 = pkin(3) * t390;
t257 = pkin(4) * t317 - qJ(5) * t318 - qJD(5) * t326 + t344;
t225 = t358 * t249 + t355 * t257;
t305 = t312 + t424;
t389 = t423 * t313;
t265 = t356 * t305 + t389;
t260 = qJD(3) * qJ(5) + t265;
t350 = -cos(pkin(9)) * pkin(1) - pkin(2);
t375 = -pkin(3) * t363 + t350;
t370 = t375 * qJD(1);
t322 = qJD(4) + t370;
t279 = pkin(4) * t323 - qJ(5) * t326 + t322;
t232 = t358 * t260 + t355 * t279;
t272 = t312 * t423 - t302;
t286 = pkin(3) * t395 + pkin(4) * t326 + qJ(5) * t323;
t239 = t358 * t272 + t355 * t286;
t391 = t361 * t424;
t270 = pkin(4) * t325 - qJ(5) * t328 - qJD(5) * t337 + t391;
t387 = qJD(3) * t407;
t314 = qJD(4) * t363 - t361 * t387;
t315 = -qJD(4) * t361 - t363 * t387;
t278 = t314 * t423 + t356 * t315;
t235 = t355 * t270 + t358 * t278;
t288 = -pkin(4) * t371 - qJ(5) * t337 + t375;
t295 = -t356 * t333 + t334 * t423;
t247 = t355 * t288 + t358 * t295;
t400 = MDP(14) * t355;
t341 = qJD(1) * t350;
t226 = pkin(8) * t310 + t232;
t398 = qJD(6) * t226;
t394 = t363 * MDP(11);
t224 = -t249 * t355 + t358 * t257;
t231 = -t260 * t355 + t358 * t279;
t234 = t358 * t270 - t278 * t355;
t238 = -t272 * t355 + t358 * t286;
t246 = t358 * t288 - t295 * t355;
t271 = t312 * t356 + t389;
t277 = t314 * t356 - t423 * t315;
t221 = -pkin(8) * t417 + t225;
t222 = pkin(5) * t323 - pkin(8) * t309 + t231;
t385 = -qJD(6) * t222 - t221;
t349 = -pkin(3) * t423 - pkin(4);
t384 = -t336 * t317 - t320 * t402;
t264 = t305 * t423 - t302;
t217 = t222 * t362 - t226 * t360;
t218 = t222 * t360 + t226 * t362;
t383 = -t231 * t355 + t232 * t358;
t233 = -pkin(5) * t371 - pkin(8) * t412 + t246;
t241 = -pkin(8) * t413 + t247;
t382 = t233 * t362 - t241 * t360;
t381 = t233 * t360 + t241 * t362;
t237 = qJD(6) * t269 + t318 * t338;
t380 = t237 * t371 - t267 * t325;
t244 = -t328 * t336 - t338 * t397;
t290 = t336 * t337;
t379 = -t244 * t320 + t290 * t317;
t378 = t251 * t337 + t294 * t318;
t374 = 0.2e1 * qJD(3) * t341;
t332 = t425 * t358;
t373 = pkin(5) * t326 + qJD(5) * t355 + qJD(6) * t332 + t323 * t426 + t238;
t331 = t425 * t355;
t372 = pkin(8) * t415 - qJD(5) * t358 + qJD(6) * t331 + t239;
t259 = -qJD(3) * pkin(4) + qJD(5) - t264;
t369 = t259 * t328 + t378;
t368 = -t317 * t345 + t318 * t349 + (-qJD(5) + t259) * t323;
t366 = t265 * MDP(13) + (t309 * t355 + t310 * t358) * MDP(16) + t383 * MDP(17);
t339 = -t358 * pkin(5) + t349;
t321 = t323 ^ 2;
t276 = pkin(5) * t413 + t294;
t254 = pkin(5) * t414 + t277;
t250 = -pkin(5) * t415 + t271;
t243 = -pkin(5) * t310 + t259;
t242 = pkin(5) * t417 + t251;
t228 = -pkin(8) * t414 + t235;
t223 = pkin(5) * t325 - t328 * t426 + t234;
t220 = pkin(5) * t317 - pkin(8) * t416 + t224;
t219 = t362 * t220;
t1 = [0.2e1 * t390 * t428 - 0.2e1 * t393 * t429 + MDP(7) * t408 - MDP(8) * t410 + (-t348 * t408 + t361 * t374) * MDP(10) + (t348 * t410 + t363 * t374) * MDP(11) + (t252 * t371 - t264 * t328 - t265 * t325 + t277 * t326 - t278 * t323 - t295 * t317 + t378) * MDP(12) + (t422 + t252 * t295 - t264 * t277 + t265 * t278 + (t322 + t370) * t391) * MDP(13) + (-t224 * t371 + t231 * t325 + t234 * t323 + t246 * t317 - t277 * t310 + t355 * t369) * MDP(14) + (t225 * t371 - t232 * t325 - t235 * t323 - t247 * t317 + t277 * t309 + t358 * t369) * MDP(15) + (-t234 * t309 + t235 * t310 + (-t224 * t337 - t231 * t328 - t246 * t318) * t358 + (-t225 * t337 - t232 * t328 - t247 * t318) * t355) * MDP(16) + (t224 * t246 + t225 * t247 + t231 * t234 + t232 * t235 + t259 * t277 + t422) * MDP(17) + (-t236 * t290 + t244 * t269) * MDP(18) + (-t236 * t289 + t237 * t290 - t244 * t267 - t245 * t269) * MDP(19) + (-t379 + t406) * MDP(20) + (t380 + t405) * MDP(21) + (-t317 * t371 + t320 * t325) * MDP(22) + ((t223 * t362 - t228 * t360) * t320 + t382 * t317 - (-t221 * t360 + t219) * t371 + t217 * t325 + t254 * t267 + t276 * t237 + t242 * t289 + t243 * t245 + (t218 * t371 - t320 * t381) * qJD(6)) * MDP(23) + (-(t223 * t360 + t228 * t362) * t320 - t381 * t317 + (t220 * t360 + t221 * t362) * t371 - t218 * t325 + t254 * t269 + t276 * t236 - t242 * t290 + t243 * t244 + (t217 * t371 - t320 * t382) * qJD(6)) * MDP(24); (-t317 * t337 - t318 * t371 + t325 * t326) * MDP(12) + (t252 * t337 - t264 * t325 - t421) * MDP(13) + (-t310 * t325 - t317 * t413 - t371 * t417) * MDP(14) + (t309 * t325 - t317 * t412 - t371 * t416) * MDP(15) + (-t224 * t413 + t225 * t412 + t259 * t325 - t421) * MDP(17) + (-t380 + t405) * MDP(23) + (t379 + t406) * MDP(24) + (-t361 * MDP(10) - t394) * t364 + ((-MDP(15) * t358 - MDP(12) - t400) * t323 + t366) * t328; ((t265 - t271) * t326 + (-t264 + t272) * t323 + (-t317 * t356 - t318 * t423) * pkin(3)) * MDP(12) + (t264 * t271 - t265 * t272 + (-t251 * t423 + t252 * t356 - t322 * t395) * pkin(3)) * MDP(13) + (-t231 * t326 - t238 * t323 - t251 * t358 + t271 * t310 + t355 * t368) * MDP(14) + (t232 * t326 + t239 * t323 + t251 * t355 - t271 * t309 + t358 * t368) * MDP(15) + (t238 * t309 - t239 * t310 + (qJD(5) * t310 - t231 * t323 + t225) * t358 + (qJD(5) * t309 - t232 * t323 - t224) * t355) * MDP(16) + (-t231 * t238 - t232 * t239 + t251 * t349 - t259 * t271 + (-t224 * t355 + t225 * t358) * t345 + t383 * qJD(5)) * MDP(17) + (t236 * t338 - t269 * t403) * MDP(18) + (-t236 * t336 - t237 * t338 + t267 * t403 - t269 * t402) * MDP(19) + (-t419 - t427) * MDP(20) + (t384 + t420) * MDP(21) - t320 * t326 * MDP(22) + ((-t331 * t362 - t332 * t360) * t317 + t339 * t237 + t242 * t336 - t217 * t326 - t250 * t267 + (t360 * t372 - t362 * t373) * t320 + t402 * t243) * MDP(23) + (-(-t331 * t360 + t332 * t362) * t317 + t339 * t236 + t242 * t338 + t218 * t326 - t250 * t269 + (t360 * t373 + t362 * t372) * t320 - t403 * t243) * MDP(24) + (-t361 * t428 + t429) * qJD(1) ^ 2 + (-MDP(10) * t395 - qJD(1) * t394) * t341; (-t326 ^ 2 - t321) * MDP(12) + (t264 * t326 + t344) * MDP(13) + (t310 * t326 + t317 * t358) * MDP(14) + (-t309 * t326 - t317 * t355 - t321 * t358) * MDP(15) + (t224 * t358 + t225 * t355 - t259 * t326) * MDP(17) + (t384 - t420) * MDP(23) + (-t419 + t427) * MDP(24) + (-t355 ^ 2 - t358 ^ 2) * MDP(16) * t318 + (-t323 * t400 + t366) * t323; (t309 * t323 + t417) * MDP(14) + (t310 * t323 + t416) * MDP(15) + (-t309 ^ 2 - t310 ^ 2) * MDP(16) + (t231 * t309 - t232 * t310 + t251) * MDP(17) + (t269 * t320 + t237) * MDP(23) + (t320 * t430 + (-t309 * t320 + t376) * t360 + t404) * MDP(24); -t267 ^ 2 * MDP(19) + (t267 * t320 + t404) * MDP(20) + t317 * MDP(22) + (t218 * t320 + t219) * MDP(23) + (t217 * t320 + t243 * t267) * MDP(24) + (MDP(18) * t267 + MDP(19) * t269 + MDP(21) * t320 - MDP(23) * t243) * t269 + (MDP(21) * t376 - MDP(23) * t398 + MDP(24) * t385) * t362 + (t376 * MDP(20) + (-qJD(6) * t310 - t416) * MDP(21) + t385 * MDP(23) + (-t220 + t398) * MDP(24)) * t360;];
tauc  = t1;
