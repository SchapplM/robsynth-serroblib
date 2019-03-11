% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:17
% EndTime: 2019-03-09 03:06:24
% DurationCPUTime: 2.97s
% Computational Cost: add. (3940->363), mult. (9210->475), div. (0->0), fcn. (6285->8), ass. (0->147)
t312 = sin(pkin(9)) * pkin(1) + pkin(7);
t384 = qJ(4) + t312;
t326 = cos(qJ(3));
t400 = cos(pkin(10));
t355 = t400 * t326;
t308 = qJD(1) * t355;
t320 = sin(pkin(10));
t324 = sin(qJ(3));
t371 = t324 * qJD(1);
t295 = -t320 * t371 + t308;
t292 = qJD(5) - t295;
t409 = MDP(6) * (t324 ^ 2 - t326 ^ 2);
t408 = t326 * MDP(5);
t304 = t320 * t326 + t324 * t400;
t297 = t304 * qJD(1);
t323 = sin(qJ(5));
t325 = cos(qJ(5));
t280 = qJD(3) * t323 + t297 * t325;
t366 = qJD(1) * qJD(3);
t358 = t324 * t366;
t331 = qJD(3) * t308 - t320 * t358;
t251 = qJD(5) * t280 + t323 * t331;
t369 = t325 * qJD(3);
t278 = t297 * t323 - t369;
t296 = t304 * qJD(3);
t337 = -t320 * t324 + t355;
t407 = -t251 * t337 + t296 * t278;
t373 = qJD(5) * t323;
t250 = -qJD(5) * t369 + t297 * t373 - t325 * t331;
t379 = t250 * t337 + t280 * t296;
t314 = -cos(pkin(9)) * pkin(1) - pkin(2);
t342 = -pkin(3) * t326 + t314;
t263 = -pkin(4) * t337 - pkin(8) * t304 + t342;
t301 = t384 * t324;
t302 = t384 * t326;
t266 = -t320 * t301 + t302 * t400;
t378 = t323 * t263 + t325 * t266;
t350 = t384 * qJD(1);
t282 = t326 * qJD(2) - t350 * t324;
t283 = t324 * qJD(2) + t326 * t350;
t365 = qJD(1) * qJD(4);
t406 = -t283 * qJD(3) - t324 * t365;
t405 = MDP(19) + MDP(21);
t273 = t320 * t283;
t401 = qJD(3) * pkin(3);
t276 = t282 + t401;
t236 = t276 * t400 - t273;
t233 = -qJD(3) * pkin(4) - t236;
t219 = t278 * pkin(5) - t280 * qJ(6) + t233;
t290 = qJD(1) * t296;
t311 = pkin(3) * t320 + pkin(8);
t387 = t311 * t290;
t404 = t219 * t292 - t387;
t403 = t280 ^ 2;
t402 = pkin(5) * t290;
t399 = qJ(6) * t290;
t270 = t282 * qJD(3) + t326 * t365;
t226 = t270 * t320 - t400 * t406;
t210 = pkin(5) * t251 + qJ(6) * t250 - qJD(6) * t280 + t226;
t398 = t210 * t323;
t356 = t400 * t283;
t237 = t320 * t276 + t356;
t234 = qJD(3) * pkin(8) + t237;
t332 = t342 * qJD(1);
t294 = qJD(4) + t332;
t247 = -pkin(4) * t295 - pkin(8) * t297 + t294;
t218 = t234 * t325 + t247 * t323;
t397 = t218 * t292;
t396 = t250 * t323;
t395 = t278 * t295;
t394 = t278 * t323;
t393 = t280 * t278;
t351 = t280 * t292;
t392 = t280 * t325;
t391 = t292 * t323;
t299 = t337 * qJD(3);
t390 = t299 * t323;
t389 = t299 * t325;
t388 = t304 * t325;
t287 = t323 * t290;
t327 = qJD(3) ^ 2;
t386 = t324 * t327;
t288 = t325 * t290;
t385 = t326 * t327;
t238 = t282 * t320 + t356;
t347 = pkin(5) * t323 - qJ(6) * t325;
t383 = qJD(6) * t323 - t292 * t347 + t238;
t382 = -t251 * t388 - t278 * t389;
t239 = t282 * t400 - t273;
t258 = pkin(3) * t371 + pkin(4) * t297 - pkin(8) * t295;
t381 = t325 * t239 + t323 * t258;
t372 = qJD(5) * t325;
t380 = -t323 * t251 - t278 * t372;
t377 = t295 * t391 + t288;
t268 = t292 * t389;
t376 = t304 * t288 + t268;
t375 = t292 * t372 + t287;
t306 = qJD(1) * t314;
t368 = t326 * MDP(11);
t217 = -t234 * t323 + t247 * t325;
t367 = qJD(6) - t217;
t363 = MDP(20) - MDP(23);
t362 = t324 * t401;
t361 = t280 * t390;
t360 = t311 * t373;
t359 = t292 * t373;
t227 = t400 * t270 + t406 * t320;
t309 = pkin(3) * t358;
t252 = t290 * pkin(4) - pkin(8) * t331 + t309;
t335 = t325 * t227 - t234 * t373 + t247 * t372 + t323 * t252;
t207 = qJD(6) * t292 + t335 + t399;
t213 = -pkin(5) * t292 + t367;
t354 = -t213 * t295 + t207;
t349 = t323 * t227 + t234 * t372 + t247 * t373 - t325 * t252;
t208 = t349 - t402;
t214 = qJ(6) * t292 + t218;
t353 = t214 * t295 + t208;
t352 = qJD(3) * t384;
t284 = qJD(4) * t326 - t324 * t352;
t285 = -qJD(4) * t324 - t326 * t352;
t245 = t284 * t320 - t400 * t285;
t265 = t400 * t301 + t302 * t320;
t313 = -pkin(3) * t400 - pkin(4);
t348 = t325 * pkin(5) + t323 * qJ(6);
t346 = t207 * t325 + t208 * t323;
t345 = t213 * t325 - t214 * t323;
t344 = t226 * t304 - t266 * t290;
t341 = 0.2e1 * qJD(3) * t306;
t339 = t304 * t372 + t390;
t338 = t304 * t373 - t389;
t336 = t219 * t280 + t349;
t246 = t284 * t400 + t320 * t285;
t259 = pkin(4) * t296 - pkin(8) * t299 + t362;
t334 = t325 * t246 + t323 * t259 + t263 * t372 - t266 * t373;
t333 = t233 * t292 - t387;
t300 = -t348 + t313;
t243 = pkin(5) * t280 + qJ(6) * t278;
t228 = t304 * t347 + t265;
t223 = t278 * t292 - t250;
t221 = pkin(5) * t337 - t263 * t325 + t266 * t323;
t220 = -qJ(6) * t337 + t378;
t216 = -pkin(5) * t297 + t239 * t323 - t258 * t325;
t215 = qJ(6) * t297 + t381;
t212 = t347 * t299 + (qJD(5) * t348 - qJD(6) * t325) * t304 + t245;
t211 = -pkin(5) * t296 + t378 * qJD(5) + t246 * t323 - t259 * t325;
t209 = qJ(6) * t296 - qJD(6) * t337 + t334;
t1 = [0.2e1 * t358 * t408 - 0.2e1 * t366 * t409 + MDP(7) * t385 - MDP(8) * t386 + (-t312 * t385 + t324 * t341) * MDP(10) + (t312 * t386 + t326 * t341) * MDP(11) + (t227 * t337 - t236 * t299 - t237 * t296 + t245 * t297 + t246 * t295 + t265 * t331 + t344) * MDP(12) + (t226 * t265 + t227 * t266 - t236 * t245 + t237 * t246 + (t294 + t332) * t362) * MDP(13) + (-t250 * t388 - t280 * t338) * MDP(14) + (-t361 + (t396 + (-t392 + t394) * qJD(5)) * t304 + t382) * MDP(15) + (-t304 * t359 + t376 + t379) * MDP(16) + (-t287 * t304 - t292 * t339 - t407) * MDP(17) + (-t290 * t337 + t292 * t296) * MDP(18) + (t349 * t337 + t217 * t296 + t245 * t278 + t265 * t251 + ((-qJD(5) * t266 + t259) * t292 + t263 * t290 + t233 * qJD(5) * t304) * t325 + ((-qJD(5) * t263 - t246) * t292 + t233 * t299 + t344) * t323) * MDP(19) + (-t218 * t296 + t226 * t388 - t233 * t338 + t245 * t280 - t265 * t250 - t290 * t378 - t292 * t334 + t335 * t337) * MDP(20) + (t208 * t337 - t211 * t292 + t212 * t278 - t213 * t296 + t219 * t339 - t221 * t290 + t228 * t251 + t304 * t398) * MDP(21) + (-t209 * t278 + t211 * t280 - t220 * t251 - t221 * t250 + t345 * t299 + (-t207 * t323 + t208 * t325 + (-t213 * t323 - t214 * t325) * qJD(5)) * t304) * MDP(22) + (-t207 * t337 + t209 * t292 - t210 * t388 - t212 * t280 + t214 * t296 + t219 * t338 + t220 * t290 + t228 * t250) * MDP(23) + (t207 * t220 + t208 * t221 + t209 * t214 + t210 * t228 + t211 * t213 + t212 * t219) * MDP(24); (t299 * t295 + t296 * t297 - t331 * t337) * MDP(12) + (-t226 * t337 - t236 * t296 + t237 * t299) * MDP(13) + (-t268 + t379) * MDP(20) + (t361 + t382) * MDP(22) + (t376 - t379) * MDP(23) + (-t210 * t337 + t213 * t390 + t214 * t389 + t219 * t296) * MDP(24) + (-t324 * MDP(10) - t368) * t327 + (t227 * MDP(13) - MDP(22) * t396 + t346 * MDP(24) + (-t325 * MDP(20) - t323 * t405 - MDP(12)) * t290 + ((t392 + t394) * MDP(22) + t345 * MDP(24) + (t323 * t363 - t325 * t405) * t292) * qJD(5)) * t304 + t405 * (-t292 * t390 + t407); ((t237 - t238) * t297 + (-t239 + t236) * t295 + (-t320 * t290 - t331 * t400) * pkin(3)) * MDP(12) + (t236 * t238 - t237 * t239 + (-t226 * t400 + t227 * t320 - t294 * t371) * pkin(3)) * MDP(13) + (t325 * t351 - t396) * MDP(14) + ((-t250 + t395) * t325 - t280 * t391 + t380) * MDP(15) + (-t292 * t295 * t325 - t280 * t297 + t375) * MDP(16) + (t278 * t297 - t359 + t377) * MDP(17) - t292 * t297 * MDP(18) + (-t217 * t297 - t238 * t278 + t313 * t251 + (-t226 + (-qJD(5) * t311 - t258) * t292) * t325 + (t239 * t292 + t333) * t323) * MDP(19) + (t218 * t297 + t226 * t323 - t238 * t280 - t313 * t250 + (t360 + t381) * t292 + t333 * t325) * MDP(20) + (-t210 * t325 + t213 * t297 + t251 * t300 + (-t311 * t372 + t216) * t292 - t383 * t278 + t404 * t323) * MDP(21) + (t215 * t278 - t216 * t280 + (-t251 * t311 + (t280 * t311 + t213) * qJD(5) + t354) * t325 + (-t250 * t311 + (t278 * t311 - t214) * qJD(5) + t353) * t323) * MDP(22) + (-t398 - t214 * t297 + t250 * t300 + (-t215 - t360) * t292 + t383 * t280 - t404 * t325) * MDP(23) + (t210 * t300 - t213 * t216 - t214 * t215 - t383 * t219 + (qJD(5) * t345 + t346) * t311) * MDP(24) + (-t324 * t408 + t409) * qJD(1) ^ 2 + (-MDP(10) * t371 - qJD(1) * t368) * t306; -t295 ^ 2 * MDP(12) + (-t237 * t295 + t309) * MDP(13) + t377 * MDP(19) + t380 * MDP(22) + t375 * MDP(23) + (-MDP(12) * t297 + t236 * MDP(13) - t219 * MDP(24) - t278 * t405 - t280 * t363) * t297 + (t290 * MDP(21) + (t250 + t395) * MDP(22) + (qJD(5) * t214 - t353) * MDP(24) + (-MDP(20) * t292 - t295 * MDP(23)) * t292) * t325 + (-t290 * MDP(20) + (qJD(5) * t213 + t354) * MDP(24) + MDP(22) * t351 + (-qJD(5) * MDP(19) - MDP(21) * t292) * t292) * t323; MDP(14) * t393 + (-t278 ^ 2 + t403) * MDP(15) + t223 * MDP(16) + (-t251 + t351) * MDP(17) + t290 * MDP(18) + (-t233 * t280 - t349 + t397) * MDP(19) + (t217 * t292 + t233 * t278 - t335) * MDP(20) + (-t243 * t278 - t336 + t397 + 0.2e1 * t402) * MDP(21) + (pkin(5) * t250 - qJ(6) * t251 + (t214 - t218) * t280 + (t213 - t367) * t278) * MDP(22) + (0.2e1 * t399 - t219 * t278 + t243 * t280 + (0.2e1 * qJD(6) - t217) * t292 + t335) * MDP(23) + (-pkin(5) * t208 + qJ(6) * t207 - t213 * t218 + t214 * t367 - t219 * t243) * MDP(24); (-qJD(3) * t297 + t393) * MDP(21) + t223 * MDP(22) + (-t292 ^ 2 - t403) * MDP(23) + (-t214 * t292 + t336 - t402) * MDP(24);];
tauc  = t1;
