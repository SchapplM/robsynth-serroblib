% Calculate vector of inverse dynamics joint torques for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:55
% EndTime: 2019-12-31 19:42:00
% DurationCPUTime: 3.75s
% Computational Cost: add. (1148->378), mult. (2332->459), div. (0->0), fcn. (1258->6), ass. (0->169)
t323 = sin(qJ(2));
t397 = qJD(1) * t323;
t293 = pkin(6) * t397;
t257 = -qJ(4) * t397 + t293;
t372 = qJD(3) + t257;
t425 = pkin(2) + pkin(3);
t374 = t425 * qJD(2);
t240 = -t374 + t372;
t433 = t323 * t425;
t326 = cos(qJ(2));
t383 = qJD(1) * qJD(2);
t370 = t326 * t383;
t382 = qJDD(1) * t323;
t432 = t370 + t382;
t396 = qJD(1) * t326;
t431 = t425 * qJDD(2);
t401 = t326 * pkin(2) + t323 * qJ(3);
t430 = -pkin(1) - t401;
t422 = pkin(6) - qJ(4);
t294 = pkin(6) * t396;
t259 = -qJ(4) * t396 + t294;
t316 = qJD(2) * qJ(3);
t247 = -t259 - t316;
t324 = sin(qJ(1));
t327 = cos(qJ(1));
t429 = g(1) * t327 + g(2) * t324;
t371 = t323 * t383;
t270 = qJ(4) * t371;
t391 = qJD(2) * t323;
t347 = pkin(6) * t391 + qJD(4) * t326;
t381 = qJDD(1) * t326;
t290 = pkin(6) * t381;
t313 = qJDD(2) * qJ(3);
t314 = qJD(2) * qJD(3);
t376 = t290 + t313 + t314;
t228 = qJ(4) * t381 + t347 * qJD(1) - t270 - t376;
t226 = qJDD(2) * pkin(4) - t228;
t305 = t326 * pkin(3);
t375 = t305 + t401;
t253 = pkin(1) + t375;
t356 = pkin(4) * t323 + pkin(7) * t326;
t235 = t356 + t253;
t241 = qJD(2) * pkin(4) - t247;
t389 = qJD(4) * t323;
t390 = qJD(2) * t326;
t245 = t422 * t390 - t389;
t254 = qJDD(5) + t432;
t264 = t422 * t323;
t276 = qJD(5) + t397;
t315 = -pkin(7) - t425;
t275 = pkin(6) * t370;
t289 = pkin(6) * t382;
t369 = qJDD(3) + t275 + t289;
t331 = -qJ(4) * t432 - qJD(1) * t389 + t369;
t225 = t315 * qJDD(2) + t331;
t248 = -qJD(1) * pkin(1) - pkin(2) * t396 - qJ(3) * t397;
t238 = pkin(3) * t396 + qJD(4) - t248;
t233 = t356 * qJD(1) + t238;
t362 = -qJD(5) * t233 - t225;
t428 = -(qJD(5) * t235 + t245) * t276 + (qJD(2) * t241 + t362) * t323 - t226 * t326 - t264 * t254;
t421 = pkin(6) * qJDD(2);
t427 = (qJD(1) * t430 + t248) * qJD(2) - t421;
t286 = qJ(3) * t396;
t343 = pkin(4) * t326 + t315 * t323;
t423 = g(3) * t323;
t426 = (t343 * qJD(1) + qJD(5) * t315 + t286) * t276 + t429 * t326 - t226 + t423;
t311 = g(1) * t324;
t424 = g(2) * t327;
t309 = g(3) * t326;
t420 = qJ(3) * t326;
t319 = qJDD(1) * pkin(1);
t419 = qJDD(2) * pkin(2);
t322 = sin(qJ(5));
t325 = cos(qJ(5));
t349 = qJD(2) * qJD(5) + t381;
t387 = qJD(5) * t326;
t373 = t322 * t387;
t403 = qJD(1) * t373 + t325 * t371;
t229 = -qJDD(2) * t322 - t349 * t325 + t403;
t418 = t229 * t322;
t417 = t235 * t254;
t386 = t325 * qJD(2);
t255 = t322 * t396 - t386;
t416 = t255 * t276;
t415 = t255 * t326;
t392 = qJD(2) * t322;
t256 = t325 * t396 + t392;
t414 = t256 * t276;
t413 = t256 * t326;
t412 = t276 * t325;
t411 = t323 * t324;
t410 = t323 * t327;
t330 = qJD(1) ^ 2;
t409 = t323 * t330;
t408 = t324 * t325;
t407 = t324 * t326;
t406 = t325 * t327;
t405 = t326 * t327;
t404 = -t325 * qJDD(2) - t322 * t371;
t298 = t323 * qJD(3);
t402 = qJ(3) * t390 + t298;
t317 = t323 ^ 2;
t318 = t326 ^ 2;
t399 = t317 - t318;
t398 = t317 + t318;
t395 = qJD(2) * t247;
t394 = qJD(2) * t255;
t393 = qJD(2) * t256;
t388 = qJD(5) * t276;
t384 = -qJD(4) - t238;
t380 = t276 * t322 * t323;
t379 = t323 * t412;
t378 = t326 * t409;
t377 = t290 + 0.2e1 * t313 + 0.2e1 * t314;
t368 = t311 - t424;
t366 = -qJD(2) * pkin(2) + qJD(3);
t365 = qJD(1) * t253 + t238;
t338 = t343 * qJD(2);
t353 = pkin(2) * t381 + t432 * qJ(3) + qJD(1) * t298 + t319;
t340 = pkin(3) * t381 + qJDD(4) + t353;
t221 = qJD(1) * t338 + t356 * qJDD(1) + t340;
t236 = t315 * qJD(2) + t372;
t363 = qJD(5) * t236 - t221;
t359 = t327 * pkin(1) + pkin(2) * t405 + t324 * pkin(6) + qJ(3) * t410;
t358 = -g(1) * t410 - g(2) * t411 + t289 + t309;
t357 = t323 * t374;
t329 = qJD(2) ^ 2;
t355 = pkin(6) * t329 + t424;
t260 = t293 + t366;
t262 = t294 + t316;
t352 = t260 * t326 - t262 * t323;
t351 = -qJDD(3) - t358;
t348 = -0.2e1 * pkin(1) * t383 - t421;
t346 = t275 - t351;
t345 = -t254 * t322 - t325 * t388;
t344 = -t254 * t325 + t322 * t388;
t341 = -t355 + 0.2e1 * t319;
t224 = -qJD(1) * t357 + t340;
t237 = -t357 + t402;
t337 = -qJD(1) * t237 - qJDD(1) * t253 - t224 + t424;
t336 = t346 - t419;
t335 = -t315 * t254 + (-t241 + t259) * t276;
t231 = pkin(2) * t371 - t353;
t246 = pkin(2) * t391 - t402;
t334 = -qJD(1) * t246 - qJDD(1) * t430 - t231 - t355;
t239 = -pkin(6) * t371 + t376;
t242 = t369 - t419;
t333 = t352 * qJD(2) + t239 * t326 + t242 * t323;
t321 = qJ(3) + pkin(4);
t307 = t327 * pkin(6);
t280 = g(1) * t407;
t279 = g(1) * t411;
t274 = qJ(3) * t405;
t272 = qJ(3) * t407;
t265 = t422 * t326;
t258 = pkin(2) * t397 - t286;
t252 = -t322 * t324 + t323 * t406;
t251 = -t322 * t410 - t408;
t250 = -t322 * t327 - t323 * t408;
t249 = t322 * t411 - t406;
t244 = -qJ(4) * t391 + t347;
t243 = -t425 * t397 + t286;
t232 = t338 + t402;
t230 = t256 * qJD(5) + t322 * t381 + t404;
t227 = t331 - t431;
t223 = t233 * t322 + t236 * t325;
t222 = t233 * t325 - t236 * t322;
t220 = t325 * t221;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t317 + 0.2e1 * t323 * t370) * MDP(4) + 0.2e1 * (t323 * t381 - t399 * t383) * MDP(5) + (qJDD(2) * t323 + t326 * t329) * MDP(6) + (qJDD(2) * t326 - t323 * t329) * MDP(7) + (t348 * t323 + t341 * t326 + t280) * MDP(9) + (-t341 * t323 + t348 * t326 - t279) * MDP(10) + (t323 * t427 + t334 * t326 + t280) * MDP(11) + (t398 * qJDD(1) * pkin(6) + t333 - t429) * MDP(12) + (t334 * t323 - t326 * t427 + t279) * MDP(13) + (t333 * pkin(6) - g(1) * t307 - g(2) * t359 + t248 * t246 + (t231 - t311) * t430) * MDP(14) + (qJDD(2) * t265 + t279 + (t365 * t326 - t244) * qJD(2) - t337 * t323) * MDP(15) + (qJDD(2) * t264 - t280 + (t365 * t323 + t245) * qJD(2) + t337 * t326) * MDP(16) + ((-qJD(2) * t240 - qJDD(1) * t265 + t228 + (-qJD(2) * t264 + t244) * qJD(1)) * t326 + (-t395 - qJDD(1) * t264 - t227 + (qJD(2) * t265 - t245) * qJD(1)) * t323 + t429) * MDP(17) + (t227 * t264 + t240 * t245 - t228 * t265 + t247 * t244 + t224 * t253 + t238 * t237 - g(1) * (-qJ(4) * t327 + t307) - g(2) * (pkin(3) * t405 + t359) + (-g(1) * (t430 - t305) + g(2) * qJ(4)) * t324) * MDP(18) + (-t229 * t325 * t326 + (-t323 * t386 - t373) * t256) * MDP(19) + ((t255 * t325 + t256 * t322) * t391 + (t418 - t230 * t325 + (t255 * t322 - t256 * t325) * qJD(5)) * t326) * MDP(20) + ((t276 * t386 + t229) * t323 + (t344 - t393) * t326) * MDP(21) + ((-t276 * t392 + t230) * t323 + (-t345 + t394) * t326) * MDP(22) + (t254 * t323 + t276 * t390) * MDP(23) + (t222 * t390 - g(1) * t250 - g(2) * t252 + t220 * t323 - t265 * t230 + t244 * t255 + (t232 * t276 + t417 + (-t236 * t323 - t241 * t326 - t264 * t276) * qJD(5)) * t325 + t428 * t322) * MDP(24) + (-t223 * t390 - g(1) * t249 - g(2) * t251 + t265 * t229 + t244 * t256 + (-(-qJD(5) * t264 + t232) * t276 - t417 + t363 * t323 + t241 * t387) * t322 + t428 * t325) * MDP(25) + t368 * MDP(2) + t429 * MDP(3); -MDP(4) * t378 + t399 * t330 * MDP(5) + MDP(6) * t382 + MDP(7) * t381 + qJDD(2) * MDP(8) + (pkin(1) * t409 - t358) * MDP(9) + (t423 - t290 + (pkin(1) * t330 + t429) * t326) * MDP(10) + (0.2e1 * t419 + (-t248 * t323 + t258 * t326) * qJD(1) + t351) * MDP(11) + ((-pkin(2) * t323 + t420) * qJDD(1) + ((t262 - t316) * t323 + (-t260 + t366) * t326) * qJD(1)) * MDP(12) + ((qJD(1) * t258 - g(3)) * t323 + (qJD(1) * t248 - t429) * t326 + t377) * MDP(13) + (t239 * qJ(3) + t262 * qJD(3) - t242 * pkin(2) - t248 * t258 - g(1) * (-pkin(2) * t410 + t274) - g(2) * (-pkin(2) * t411 + t272) - g(3) * t401 - t352 * qJD(1) * pkin(6)) * MDP(14) + (qJD(2) * t257 + t270 + (-g(3) + (-pkin(6) * qJD(2) - t243) * qJD(1)) * t323 + (-qJ(4) * qJDD(1) + t384 * qJD(1) - t429) * t326 + t377) * MDP(15) + (-qJ(4) * t382 - qJD(2) * t259 - 0.2e1 * t431 + ((-qJ(4) * qJD(2) + t243) * t326 + t384 * t323) * qJD(1) + t346) * MDP(16) + (-t420 + t433) * qJDD(1) * MDP(17) + (-g(1) * t274 - g(2) * t272 - g(3) * t375 - t228 * qJ(3) - t227 * t425 - t238 * t243 - t240 * t259 - t247 * t372 + t429 * t433) * MDP(18) + (t256 * t412 - t418) * MDP(19) + ((-t229 - t416) * t325 + (-t230 - t414) * t322) * MDP(20) + ((-t379 + t413) * qJD(1) + t345) * MDP(21) + ((t380 - t415) * qJD(1) + t344) * MDP(22) - t276 * MDP(23) * t396 + (-t222 * t396 - t321 * t230 - t255 * t372 + t335 * t322 - t325 * t426) * MDP(24) + (t223 * t396 + t321 * t229 - t256 * t372 + t322 * t426 + t335 * t325) * MDP(25); (-qJD(2) * t262 + t336) * MDP(14) + (-qJDD(2) * pkin(3) - qJ(4) * t370 + t336 + t395) * MDP(18) + (t345 + t394) * MDP(24) + (t344 + t393) * MDP(25) + (MDP(13) + MDP(15)) * (-t317 * t330 - t329) + (-MDP(11) + MDP(16)) * (qJDD(2) + t378) + ((-MDP(18) * qJ(4) + MDP(12) - MDP(17)) * qJDD(1) + (t248 * MDP(14) + t384 * MDP(18) + (-MDP(24) * t325 + MDP(25) * t322) * t276) * qJD(1)) * t323; (t340 + t368) * MDP(18) - t344 * MDP(24) + t345 * MDP(25) - t398 * MDP(17) * t330 + (t323 * MDP(15) - MDP(16) * t326) * qJDD(1) + ((t240 * t323 - t247 * t326) * MDP(18) + (-t380 - t415) * MDP(24) + (-t379 - t413) * MDP(25) + (0.2e1 * t326 * MDP(15) + (-t425 * MDP(18) + 0.2e1 * MDP(16)) * t323) * qJD(2)) * qJD(1); t256 * t255 * MDP(19) + (-t255 ^ 2 + t256 ^ 2) * MDP(20) + (t403 - t416) * MDP(21) + (t404 - t414) * MDP(22) + t254 * MDP(23) + (-g(1) * t251 + g(2) * t249 + t223 * t276 + t241 * t256 + t220) * MDP(24) + (g(1) * t252 - g(2) * t250 + t222 * t276 - t241 * t255) * MDP(25) + (-MDP(21) * t381 + (-t225 - t309) * MDP(25) + (-qJD(2) * MDP(21) + MDP(22) * t396 - MDP(24) * t236 - MDP(25) * t233) * qJD(5)) * t325 + (-qJDD(2) * MDP(21) + t349 * MDP(22) + (t362 - t309) * MDP(24) + t363 * MDP(25)) * t322;];
tau = t1;
