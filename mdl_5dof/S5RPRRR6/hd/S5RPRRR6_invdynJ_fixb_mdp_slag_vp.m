% Calculate vector of inverse dynamics joint torques for
% S5RPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:48
% EndTime: 2019-12-31 19:01:53
% DurationCPUTime: 3.10s
% Computational Cost: add. (2235->301), mult. (4666->411), div. (0->0), fcn. (3286->14), ass. (0->151)
t320 = qJD(3) + qJD(4);
t328 = sin(qJ(4));
t332 = cos(qJ(3));
t412 = cos(qJ(4));
t377 = qJD(1) * t412;
t329 = sin(qJ(3));
t391 = qJD(1) * t329;
t418 = -t328 * t391 + t332 * t377;
t420 = t418 * t320;
t324 = qJ(3) + qJ(4);
t317 = sin(t324);
t321 = qJ(1) + pkin(9);
t313 = sin(t321);
t314 = cos(t321);
t416 = g(1) * t314 + g(2) * t313;
t419 = t416 * t317;
t325 = sin(pkin(9));
t306 = pkin(1) * t325 + pkin(6);
t408 = pkin(7) + t306;
t395 = t328 * t332;
t273 = -qJD(1) * t395 - t329 * t377;
t253 = -pkin(4) * t273 - pkin(8) * t418;
t270 = qJD(5) - t418;
t415 = (pkin(8) * qJD(5) + t253) * t270;
t371 = t408 * qJD(1);
t261 = t332 * qJD(2) - t371 * t329;
t262 = qJD(2) * t329 + t371 * t332;
t407 = qJD(3) * pkin(3);
t257 = t261 + t407;
t379 = t412 * t262;
t236 = t328 * t257 + t379;
t315 = t332 * qJDD(2);
t288 = t306 * qJDD(1);
t369 = pkin(7) * qJDD(1) + t288;
t238 = qJDD(3) * pkin(3) - t262 * qJD(3) - t369 * t329 + t315;
t243 = t261 * qJD(3) + t329 * qJDD(2) + t369 * t332;
t414 = t236 * qJD(4) - t412 * t238 + t328 * t243;
t319 = qJDD(3) + qJDD(4);
t214 = -t319 * pkin(4) + t414;
t373 = qJD(3) * t408;
t268 = t329 * t373;
t269 = t332 * t373;
t275 = t408 * t329;
t276 = t408 * t332;
t349 = -t412 * t275 - t328 * t276;
t225 = t349 * qJD(4) - t412 * t268 - t328 * t269;
t396 = t328 * t262;
t235 = t412 * t257 - t396;
t232 = -t320 * pkin(4) - t235;
t280 = t412 * t329 + t395;
t255 = t320 * t280;
t372 = qJDD(1) * t412;
t385 = qJDD(1) * t329;
t358 = t328 * t385 - t332 * t372;
t246 = t255 * qJD(1) + t358;
t242 = qJDD(5) + t246;
t326 = cos(pkin(9));
t307 = -pkin(1) * t326 - pkin(2);
t285 = -pkin(3) * t332 + t307;
t348 = -t328 * t329 + t412 * t332;
t250 = -pkin(4) * t348 - pkin(8) * t280 + t285;
t252 = -t328 * t275 + t412 * t276;
t254 = t320 * t348;
t376 = qJD(4) * t412;
t389 = qJD(4) * t328;
t340 = t328 * t238 + t412 * t243 + t257 * t376 - t262 * t389;
t213 = t319 * pkin(8) + t340;
t274 = t285 * qJD(1);
t248 = -pkin(4) * t418 + pkin(8) * t273 + t274;
t366 = qJD(5) * t248 + t213;
t413 = t214 * t280 + t232 * t254 - t252 * t242 - (qJD(5) * t250 + t225) * t270 + t366 * t348;
t308 = g(3) * t317;
t318 = cos(t324);
t409 = g(3) * t318;
t384 = qJDD(1) * t332;
t245 = t328 * t384 + t329 * t372 + t420;
t327 = sin(qJ(5));
t331 = cos(qJ(5));
t387 = qJD(5) * t331;
t380 = t331 * t245 + t327 * t319 + t320 * t387;
t388 = qJD(5) * t327;
t223 = t273 * t388 + t380;
t406 = t223 * t327;
t405 = t232 * t418;
t404 = t232 * t280;
t403 = t242 * t327;
t402 = t250 * t242;
t258 = -t273 * t327 - t331 * t320;
t401 = t258 * t270;
t260 = -t273 * t331 + t320 * t327;
t400 = t260 * t270;
t399 = t280 * t331;
t398 = t318 * t327;
t397 = t318 * t331;
t394 = qJDD(2) - g(3);
t393 = -t223 * t348 + t260 * t255;
t322 = t329 ^ 2;
t392 = -t332 ^ 2 + t322;
t291 = qJD(1) * t307;
t386 = qJD(1) * qJD(3);
t383 = t329 * t407;
t382 = t280 * t403;
t381 = t242 * t399;
t375 = t329 * t386;
t374 = -t214 - t409;
t367 = t270 * t331;
t263 = pkin(3) * t375 + t285 * qJDD(1);
t219 = pkin(4) * t246 - pkin(8) * t245 + t263;
t233 = t320 * pkin(8) + t236;
t365 = qJD(5) * t233 - t219;
t311 = pkin(3) * t328 + pkin(8);
t363 = pkin(3) * t391 + qJD(5) * t311 + t253;
t239 = t328 * t261 + t379;
t362 = pkin(3) * t389 - t239;
t360 = g(1) * t313 - g(2) * t314;
t330 = sin(qJ(1));
t333 = cos(qJ(1));
t359 = g(1) * t330 - g(2) * t333;
t303 = t331 * t319;
t224 = qJD(5) * t260 + t245 * t327 - t303;
t356 = t224 * t348 - t255 * t258;
t355 = -t242 * t311 - t405;
t354 = t254 * t320 + t280 * t319;
t240 = t412 * t261 - t396;
t353 = -pkin(3) * t376 + t240;
t221 = t233 * t331 + t248 * t327;
t352 = g(3) * t398 + t214 * t327 - t221 * t273 + t232 * t387;
t220 = -t233 * t327 + t248 * t331;
t351 = t220 * t273 + t232 * t388 + t331 * t419;
t347 = -t254 * t327 - t280 * t387;
t346 = -t254 * t331 + t280 * t388;
t345 = -pkin(8) * t242 + t235 * t270 - t405;
t343 = -qJD(1) * t291 - t288 + t416;
t342 = 0.2e1 * t291 * qJD(3) - qJDD(3) * t306;
t334 = qJD(3) ^ 2;
t341 = -0.2e1 * qJDD(1) * t307 - t306 * t334 + t360;
t339 = ((t223 - t401) * t331 + (-t224 - t400) * t327) * MDP(20) + (t260 * t367 + t406) * MDP(19) + (-t270 ^ 2 * t327 + t242 * t331 - t258 * t273) * MDP(22) + (t260 * t273 + t270 * t367 + t403) * MDP(21) + (t245 - t420) * MDP(14) + (-t358 + (-qJD(1) * t280 - t273) * t320) * MDP(15) + (t273 ^ 2 - t418 ^ 2) * MDP(13) + t319 * MDP(16) + (MDP(12) * t418 + t270 * MDP(23)) * t273;
t338 = -t274 * t418 + t416 * t318 + t308 - t340;
t337 = t274 * t273 - t409 - t414 + t419;
t312 = -t412 * pkin(3) - pkin(4);
t287 = qJDD(3) * t332 - t329 * t334;
t286 = qJDD(3) * t329 + t332 * t334;
t267 = t313 * t327 + t314 * t397;
t266 = t313 * t331 - t314 * t398;
t265 = -t313 * t397 + t314 * t327;
t264 = t313 * t398 + t314 * t331;
t244 = -t255 * t320 + t319 * t348;
t229 = pkin(4) * t255 - pkin(8) * t254 + t383;
t226 = t252 * qJD(4) - t328 * t268 + t412 * t269;
t218 = t331 * t219;
t1 = [qJDD(1) * MDP(1) + t359 * MDP(2) + (g(1) * t333 + g(2) * t330) * MDP(3) + (t359 + (t325 ^ 2 + t326 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t322 + 0.2e1 * t332 * t375) * MDP(5) + 0.2e1 * (t329 * t384 - t392 * t386) * MDP(6) + t286 * MDP(7) + t287 * MDP(8) + (t342 * t329 + t341 * t332) * MDP(10) + (-t341 * t329 + t342 * t332) * MDP(11) + (t245 * t280 - t254 * t273) * MDP(12) + (t245 * t348 - t246 * t280 + t254 * t418 + t255 * t273) * MDP(13) + t354 * MDP(14) + t244 * MDP(15) + (-t226 * t320 + t246 * t285 + t255 * t274 - t263 * t348 + t360 * t318 + t319 * t349 - t383 * t418) * MDP(17) + (-t225 * t320 + t245 * t285 - t252 * t319 + t254 * t274 + t263 * t280 - t273 * t383 - t360 * t317) * MDP(18) + (t223 * t399 - t346 * t260) * MDP(19) + ((-t258 * t331 - t260 * t327) * t254 + (-t406 - t224 * t331 + (t258 * t327 - t260 * t331) * qJD(5)) * t280) * MDP(20) + (-t346 * t270 + t381 + t393) * MDP(21) + (t347 * t270 + t356 - t382) * MDP(22) + (-t242 * t348 + t255 * t270) * MDP(23) + (-g(1) * t265 - g(2) * t267 - t218 * t348 + t220 * t255 - t349 * t224 + t226 * t258 + (t229 * t270 + t402 + (t233 * t348 - t252 * t270 + t404) * qJD(5)) * t331 + t413 * t327) * MDP(24) + (-g(1) * t264 - g(2) * t266 - t221 * t255 - t349 * t223 + t226 * t260 + (-(-qJD(5) * t252 + t229) * t270 - t402 - t365 * t348 - qJD(5) * t404) * t327 + t413 * t331) * MDP(25); t394 * MDP(4) + t287 * MDP(10) - t286 * MDP(11) + t244 * MDP(17) - t354 * MDP(18) + (-t356 - t382) * MDP(24) + (-t381 + t393) * MDP(25) + (t347 * MDP(24) + t346 * MDP(25)) * t270; (t239 * t320 + (t412 * t319 - t320 * t389 + t391 * t418) * pkin(3) + t337) * MDP(17) + (t240 * t320 + (t273 * t391 - t319 * t328 - t320 * t376) * pkin(3) + t338) * MDP(18) + (-g(3) * t332 + t343 * t329 + t315) * MDP(10) + (-t394 * t329 + t343 * t332) * MDP(11) + t339 + (t312 * t223 + t355 * t331 - t327 * t419 + t362 * t260 + (t363 * t327 + t353 * t331) * t270 + t352) * MDP(25) + qJDD(3) * MDP(9) + (t312 * t224 + t374 * t331 + t355 * t327 + t362 * t258 + (t353 * t327 - t363 * t331) * t270 + t351) * MDP(24) + MDP(7) * t385 + MDP(8) * t384 + (-t329 * t332 * MDP(5) + t392 * MDP(6)) * qJD(1) ^ 2; (-pkin(4) * t223 - t236 * t260 + t345 * t331 + (-t419 + t415) * t327 + t352) * MDP(25) + (t236 * t320 + t337) * MDP(17) + (t235 * t320 + t338) * MDP(18) + t339 + (-pkin(4) * t224 - t236 * t258 + t345 * t327 + (t374 - t415) * t331 + t351) * MDP(24); t260 * t258 * MDP(19) + (-t258 ^ 2 + t260 ^ 2) * MDP(20) + (t380 + t401) * MDP(21) + (t303 + t400) * MDP(22) + t242 * MDP(23) + (-g(1) * t266 + g(2) * t264 + t221 * t270 - t232 * t260 + t218) * MDP(24) + (g(1) * t267 - g(2) * t265 + t220 * t270 + t232 * t258) * MDP(25) + ((-t213 + t308) * MDP(25) + (MDP(22) * t273 - MDP(24) * t233 - MDP(25) * t248) * qJD(5)) * t331 + (qJD(5) * t273 * MDP(21) + (-qJD(5) * t320 - t245) * MDP(22) + (-t366 + t308) * MDP(24) + t365 * MDP(25)) * t327;];
tau = t1;
