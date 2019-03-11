% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:37
% EndTime: 2019-03-09 02:16:41
% DurationCPUTime: 3.14s
% Computational Cost: add. (3605->351), mult. (7830->443), div. (0->0), fcn. (5418->6), ass. (0->136)
t299 = sin(pkin(9));
t301 = -pkin(1) - qJ(3);
t387 = t301 * qJD(1);
t283 = qJD(2) + t387;
t335 = -pkin(7) * qJD(1) + t283;
t265 = t335 * t299;
t300 = cos(pkin(9));
t266 = t335 * t300;
t303 = sin(qJ(4));
t381 = cos(qJ(4));
t237 = t381 * t265 + t303 * t266;
t233 = qJD(4) * pkin(8) + t237;
t338 = t381 * t300;
t330 = qJD(1) * t338;
t353 = qJD(1) * t299;
t337 = t303 * t353;
t271 = t330 - t337;
t298 = qJD(1) * qJ(2);
t292 = qJD(3) + t298;
t279 = pkin(3) * t353 + t292;
t319 = -t381 * t299 - t303 * t300;
t386 = t319 * qJD(1);
t234 = -pkin(4) * t386 - pkin(8) * t271 + t279;
t302 = sin(qJ(5));
t304 = cos(qJ(5));
t207 = t233 * t304 + t234 * t302;
t391 = qJD(5) - t386;
t203 = qJ(6) * t391 + t207;
t394 = t203 * t391;
t253 = qJD(4) * t302 + t271 * t304;
t333 = t253 * t391;
t307 = t386 * qJD(4);
t393 = qJD(4) * qJD(5) + t307;
t354 = t299 ^ 2 + t300 ^ 2;
t392 = qJD(3) * t354;
t390 = qJ(2) * MDP(6) + t299 * MDP(7) + t300 * MDP(8) + MDP(5);
t288 = t299 * pkin(3) + qJ(2);
t318 = -t303 * t299 + t338;
t245 = -pkin(4) * t319 - pkin(8) * t318 + t288;
t378 = -pkin(7) + t301;
t277 = t378 * t299;
t278 = t378 * t300;
t249 = t381 * t277 + t303 * t278;
t355 = t302 * t245 + t304 * t249;
t389 = -t303 * t265 + t381 * t266;
t388 = -t303 * t277 + t381 * t278;
t385 = MDP(23) + MDP(25);
t343 = MDP(24) - MDP(27);
t309 = t319 * qJD(3);
t216 = qJD(1) * t309 + qJD(4) * t389;
t282 = qJD(4) * t337;
t264 = qJD(4) * t330 - t282;
t297 = qJD(1) * qJD(2);
t235 = t264 * pkin(4) - pkin(8) * t307 + t297;
t349 = qJD(5) * t304;
t350 = qJD(5) * t302;
t331 = -t302 * t216 - t233 * t349 - t234 * t350 + t304 * t235;
t380 = pkin(5) * t264;
t197 = -t331 - t380;
t384 = -t197 + t394;
t232 = -qJD(4) * pkin(4) - t389;
t251 = -t304 * qJD(4) + t271 * t302;
t208 = t251 * pkin(5) - t253 * qJ(6) + t232;
t379 = pkin(8) * t264;
t383 = t208 * t391 - t379;
t382 = t253 ^ 2;
t376 = qJ(6) * t264;
t225 = t271 * t350 - t304 * t393;
t226 = t271 * t349 + t302 * t393;
t336 = qJD(4) * t381;
t351 = qJD(4) * t303;
t311 = -qJD(3) * t271 - t265 * t336 - t266 * t351;
t199 = pkin(5) * t226 + qJ(6) * t225 - qJD(6) * t253 - t311;
t375 = t199 * t302;
t374 = t207 * t391;
t373 = t225 * t302;
t372 = t226 * t304;
t371 = t251 * t386;
t370 = t251 * t271;
t369 = t251 * t302;
t368 = t251 * t304;
t367 = t253 * t251;
t366 = t253 * t271;
t365 = t253 * t302;
t364 = t253 * t304;
t363 = t391 * t304;
t362 = t318 * t304;
t258 = t302 * t264;
t260 = t304 * t264;
t328 = pkin(5) * t302 - qJ(6) * t304;
t359 = qJD(6) * t302 - t328 * t391 + t237;
t357 = -t302 * t226 - t251 * t349;
t246 = pkin(4) * t271 - pkin(8) * t386;
t356 = t302 * t246 + t304 * t389;
t272 = -t299 * t336 - t300 * t351;
t352 = qJD(4) * t272;
t206 = -t233 * t302 + t234 * t304;
t347 = qJD(6) - t206;
t346 = qJD(4) * MDP(17);
t342 = pkin(8) * t350;
t334 = qJD(1) * t354;
t332 = t302 * t391;
t329 = -pkin(5) * t304 - qJ(6) * t302;
t313 = t304 * t216 - t233 * t350 + t234 * t349 + t302 * t235;
t196 = qJD(6) * t391 + t313 + t376;
t327 = t196 * t304 + t197 * t302;
t202 = -pkin(5) * t391 + t347;
t326 = t202 * t304 - t203 * t302;
t325 = t202 * t302 + t203 * t304;
t323 = t202 * t391 + t196;
t322 = t349 * t391 - t363 * t386 + t258;
t321 = t260 + (t302 * t386 - t350) * t391;
t317 = -t272 * t302 - t318 * t349;
t316 = t272 * t304 - t318 * t350;
t315 = t232 * t391 - t379;
t314 = t208 * t253 - t331;
t227 = qJD(4) * t388 + t309;
t273 = -t299 * t351 + t300 * t336;
t243 = pkin(4) * t273 - pkin(8) * t272 + qJD(2);
t312 = t304 * t227 + t302 * t243 + t245 * t349 - t249 * t350;
t308 = -t302 * t385 - t343 * t304;
t228 = t318 * qJD(3) + t249 * qJD(4);
t306 = (t364 + t369) * MDP(26) + t326 * MDP(28) + (t343 * t302 - t304 * t385) * t391;
t305 = qJD(1) ^ 2;
t280 = -pkin(4) + t329;
t219 = pkin(5) * t253 + qJ(6) * t251;
t215 = t318 * t328 - t388;
t211 = pkin(5) * t319 - t245 * t304 + t249 * t302;
t210 = -qJ(6) * t319 + t355;
t209 = t251 * t391 - t225;
t205 = -pkin(5) * t271 - t246 * t304 + t302 * t389;
t204 = qJ(6) * t271 + t356;
t201 = t328 * t272 - (t329 * qJD(5) + qJD(6) * t304) * t318 + t228;
t200 = -pkin(5) * t273 + qJD(5) * t355 + t227 * t302 - t243 * t304;
t198 = qJ(6) * t273 - qJD(6) * t319 + t312;
t1 = [0.2e1 * qJD(3) * MDP(9) * t334 + ((t292 + t298) * qJD(2) + (-t283 - t387) * t392) * MDP(10) + (t271 * t272 + t307 * t318) * MDP(11) + (-t264 * t318 - t271 * t273 + t272 * t386 + t307 * t319) * MDP(12) + MDP(13) * t352 - t273 * qJD(4) * MDP(14) + (-0.2e1 * qJD(2) * t386 - qJD(4) * t228 + t264 * t288 + t273 * t279) * MDP(16) + (qJD(2) * t271 - t227 * qJD(4) + t279 * t272 + t288 * t307 + t297 * t318) * MDP(17) + (-t225 * t362 + t316 * t253) * MDP(18) + ((-t365 - t368) * t272 - (-t373 + t372 + (t364 - t369) * qJD(5)) * t318) * MDP(19) + (t225 * t319 + t253 * t273 + t260 * t318 + t316 * t391) * MDP(20) + (t226 * t319 - t251 * t273 - t258 * t318 + t317 * t391) * MDP(21) + (-t264 * t319 + t273 * t391) * MDP(22) + (-t331 * t319 + t206 * t273 + t228 * t251 - t388 * t226 + ((-qJD(5) * t249 + t243) * t391 + t245 * t264 + t232 * qJD(5) * t318) * t304 + ((-qJD(5) * t245 - t227) * t391 - t249 * t264 - t311 * t318 + t232 * t272) * t302) * MDP(23) + (-t207 * t273 + t225 * t388 + t228 * t253 + t316 * t232 - t355 * t264 - t311 * t362 - t312 * t391 + t313 * t319) * MDP(24) + (t197 * t319 - t200 * t391 + t201 * t251 - t202 * t273 - t317 * t208 - t211 * t264 + t215 * t226 + t318 * t375) * MDP(25) + (-t198 * t251 + t200 * t253 - t210 * t226 - t211 * t225 + t326 * t272 - (t325 * qJD(5) + t196 * t302 - t197 * t304) * t318) * MDP(26) + (-t196 * t319 + t198 * t391 - t199 * t362 - t201 * t253 + t203 * t273 - t316 * t208 + t210 * t264 + t215 * t225) * MDP(27) + (t196 * t210 + t197 * t211 + t198 * t203 + t199 * t215 + t200 * t202 + t201 * t208) * MDP(28) + 0.2e1 * t390 * t297; MDP(16) * t352 + (-t199 * t318 - t208 * t272) * MDP(28) - t390 * t305 + (-t346 + (t365 - t368) * MDP(26) + t325 * MDP(28) + t308 * t391) * t273 + ((-t292 - t392) * MDP(10) + t386 * MDP(16) - t271 * MDP(17) + t306) * qJD(1) - ((-t372 - t373) * MDP(26) + t327 * MDP(28) + t308 * t264 + t306 * qJD(5)) * t319 + t385 * (-t226 * t318 - t272 * t251) - t343 * (-t225 * t318 + t253 * t272); -t354 * MDP(9) * t305 + (t283 * t334 + t297) * MDP(10) + (-t282 + (t330 + t271) * qJD(4)) * MDP(16) + 0.2e1 * t386 * t346 + (t321 - t370) * MDP(23) + (-t363 * t391 - t258 - t366) * MDP(24) + (-t332 * t391 + t260 - t370) * MDP(25) + ((t225 + t371) * t304 + t302 * t333 + t357) * MDP(26) + (t322 + t366) * MDP(27) + (-t208 * t271 + t323 * t302 + t304 * t384) * MDP(28); t271 ^ 2 * MDP(12) + (t282 + (-t330 + t271) * qJD(4)) * MDP(14) + (qJD(4) * t237 - t271 * t279 + t311) * MDP(16) + (t304 * t333 - t373) * MDP(18) + ((-t225 + t371) * t304 - t253 * t332 + t357) * MDP(19) + (t322 - t366) * MDP(20) + (t321 + t370) * MDP(21) - t391 * t271 * MDP(22) + (-pkin(4) * t226 - t206 * t271 - t237 * t251 + (t311 + (-pkin(8) * qJD(5) - t246) * t391) * t304 + (t389 * t391 + t315) * t302) * MDP(23) + (pkin(4) * t225 + t207 * t271 - t311 * t302 - t237 * t253 + (t342 + t356) * t391 + t315 * t304) * MDP(24) + (-t199 * t304 + t202 * t271 + t226 * t280 + (-pkin(8) * t349 + t205) * t391 - t359 * t251 + t383 * t302) * MDP(25) + (t204 * t251 - t205 * t253 + ((qJD(5) * t253 - t226) * pkin(8) + t323) * t304 + ((qJD(5) * t251 - t225) * pkin(8) - t384) * t302) * MDP(26) + (-t375 - t203 * t271 + t225 * t280 + (-t204 - t342) * t391 + t359 * t253 - t383 * t304) * MDP(27) + (t199 * t280 - t202 * t205 - t203 * t204 - t359 * t208 + (t326 * qJD(5) + t327) * pkin(8)) * MDP(28) + (-t271 * MDP(11) + (-qJD(3) - t279) * MDP(17) - MDP(12) * t386) * t386; MDP(18) * t367 + (-t251 ^ 2 + t382) * MDP(19) + t209 * MDP(20) + (-t226 + t333) * MDP(21) + t264 * MDP(22) + (-t232 * t253 + t331 + t374) * MDP(23) + (t206 * t391 + t232 * t251 - t313) * MDP(24) + (-t219 * t251 - t314 + t374 + 0.2e1 * t380) * MDP(25) + (pkin(5) * t225 - qJ(6) * t226 + (t203 - t207) * t253 + (t202 - t347) * t251) * MDP(26) + (0.2e1 * t376 - t208 * t251 + t219 * t253 + (0.2e1 * qJD(6) - t206) * t391 + t313) * MDP(27) + (-pkin(5) * t197 + qJ(6) * t196 - t202 * t207 + t347 * t203 - t208 * t219) * MDP(28); (-t264 + t367) * MDP(25) + t209 * MDP(26) + (-t391 ^ 2 - t382) * MDP(27) + (t314 - t380 - t394) * MDP(28);];
tauc  = t1;
