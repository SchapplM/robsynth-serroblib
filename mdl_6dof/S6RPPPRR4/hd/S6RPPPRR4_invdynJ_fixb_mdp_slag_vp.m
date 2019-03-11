% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:58
% EndTime: 2019-03-09 01:36:01
% DurationCPUTime: 2.29s
% Computational Cost: add. (1292->320), mult. (2158->403), div. (0->0), fcn. (1327->8), ass. (0->150)
t272 = sin(qJ(5));
t333 = qJD(1) * qJD(5);
t319 = t272 * t333;
t381 = qJD(5) * qJD(6) + t319;
t274 = cos(qJ(5));
t338 = qJD(6) * t274;
t299 = qJD(1) * t338 + qJDD(5);
t380 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4);
t379 = qJDD(1) * pkin(3) + qJDD(4);
t268 = sin(pkin(9));
t334 = qJD(1) * qJD(2);
t244 = t268 * t334;
t275 = -pkin(1) - pkin(2);
t239 = t275 * qJDD(1) + qJDD(2);
t269 = cos(pkin(9));
t331 = qJDD(1) * t268;
t354 = qJ(2) * t331 - t269 * t239;
t211 = -t244 - t354;
t209 = -t211 + t379;
t206 = qJDD(1) * pkin(7) + t209;
t240 = t275 * qJD(1) + qJD(2);
t349 = qJ(2) * qJD(1);
t216 = t240 * t269 - t268 * t349;
t213 = qJD(1) * pkin(3) + qJD(4) - t216;
t210 = qJD(1) * pkin(7) + t213;
t205 = qJD(3) * t274 + t210 * t272;
t344 = qJD(5) * t205;
t193 = -qJDD(5) * pkin(5) + qJDD(3) * t272 - t206 * t274 + t344;
t204 = -qJD(3) * t272 + t210 * t274;
t198 = -qJD(5) * pkin(5) - t204;
t306 = -pkin(5) * t272 + pkin(8) * t274;
t368 = qJ(2) * t269;
t290 = t306 + t368;
t361 = t268 * t275;
t215 = -qJ(4) + t290 + t361;
t230 = -t268 * qJ(2) + t269 * t275;
t226 = pkin(3) - t230;
t219 = pkin(7) + t226;
t318 = t274 * t333;
t329 = qJDD(1) * t272;
t220 = -qJDD(6) + t318 + t329;
t242 = qJD(1) * t272 - qJD(6);
t263 = qJD(1) * qJ(4);
t363 = t240 * t268;
t207 = t290 * qJD(1) - t263 + t363;
t345 = qJD(5) * t204;
t313 = -qJDD(5) * pkin(8) - qJD(6) * t207 - qJDD(3) * t274 - t206 * t272 - t345;
t340 = qJD(5) * t274;
t347 = qJD(2) * t268;
t373 = sin(qJ(1));
t374 = cos(qJ(1));
t222 = t374 * t268 - t373 * t269;
t371 = g(1) * t222;
t378 = (qJD(6) * t215 + t219 * t340) * t242 + (t198 * qJD(5) + t219 * t220 + t242 * t347 - t313) * t272 - t371 - t193 * t274;
t217 = t269 * t349 + t363;
t214 = t217 - t263;
t231 = t361 + t368;
t223 = -qJ(4) + t231;
t377 = (qJD(1) * t223 + t214 - t347) * qJD(5) - qJDD(5) * t219;
t271 = sin(qJ(6));
t273 = cos(qJ(6));
t328 = qJDD(1) * t274;
t355 = t381 * t271;
t201 = t271 * t328 + t299 * t273 - t355;
t337 = t271 * qJD(5);
t322 = t242 * t337;
t302 = t201 + t322;
t310 = t299 * t271 + t381 * t273;
t200 = -t273 * t328 + t310;
t336 = t273 * qJD(5);
t321 = t242 * t336;
t303 = -t200 + t321;
t376 = -t302 * MDP(25) - t303 * MDP(26);
t221 = -t373 * t268 - t374 * t269;
t305 = -g(2) * t221 + t371;
t307 = -pkin(5) * t274 - pkin(8) * t272;
t375 = (pkin(8) * qJD(6) + t307 * qJD(1)) * t242 - t305 * t274 - g(3) * t272 - t193;
t372 = g(1) * t221;
t370 = g(2) * t222;
t369 = pkin(1) * qJDD(1);
t367 = t200 * t271;
t366 = t221 * t272;
t348 = qJD(1) * t274;
t224 = t271 * t348 + t336;
t365 = t224 * t242;
t225 = t273 * t348 - t337;
t364 = t225 * t242;
t362 = t268 * t274;
t360 = t271 * t272;
t359 = t272 * t273;
t358 = t273 * t242;
t357 = qJDD(3) + g(3);
t330 = qJDD(1) * t269;
t356 = qJ(2) * t330 + t268 * t239;
t353 = t374 * pkin(1) + t373 * qJ(2);
t352 = g(1) * t373 - g(2) * t374;
t265 = t274 ^ 2;
t351 = t272 ^ 2 - t265;
t276 = qJD(5) ^ 2;
t277 = qJD(1) ^ 2;
t350 = t276 + t277;
t346 = qJD(2) * t269;
t343 = qJD(5) * t224;
t342 = qJD(5) * t225;
t341 = qJD(5) * t272;
t339 = qJD(6) * t242;
t335 = qJ(2) * qJDD(1);
t326 = qJDD(5) * t272;
t325 = qJDD(5) * t274;
t324 = 0.2e1 * t334;
t323 = t374 * pkin(2) + t353;
t320 = t269 * t334;
t199 = qJD(5) * pkin(8) + t205;
t314 = t219 * t242 + t199;
t312 = -0.2e1 * t318;
t311 = qJDD(2) - t369;
t309 = -t356 - t380;
t308 = -t373 * pkin(1) + t374 * qJ(2);
t304 = t370 + t372;
t301 = -qJD(5) * t210 - t357;
t300 = -qJD(6) * t199 + t370;
t283 = t307 * qJD(5) + t346;
t298 = -t215 * t220 - (-qJD(4) + t283) * t242;
t297 = t216 * t268 - t217 * t269;
t212 = t320 + t356;
t296 = t271 * t220 + t273 * t339;
t295 = -t273 * t220 + t271 * t339;
t292 = -t305 + t354;
t291 = g(1) * t374 + g(2) * t373;
t289 = -qJDD(1) * t223 - t304;
t288 = -t373 * pkin(2) + t308;
t287 = t292 + t379;
t286 = t296 - t343;
t285 = t295 + t342;
t284 = -g(2) * t366 - g(3) * t274 + t313;
t282 = pkin(8) * t220 + (-t198 - t204) * t242;
t281 = -qJD(1) * t214 + qJD(3) * qJD(5) - t206 + t305;
t208 = t212 + t380;
t241 = -qJD(4) + t346;
t280 = -qJD(1) * t241 - t219 * t276 - t208 + t289;
t278 = t286 * MDP(25) - t285 * MDP(26);
t235 = t272 * t276 - t325;
t234 = t274 * t276 + t326;
t203 = -t221 * t271 + t222 * t359;
t202 = -t221 * t273 - t222 * t360;
t197 = t283 * qJD(1) + t306 * qJDD(1) - t309;
t196 = t273 * t197;
t195 = t273 * t199 + t271 * t207;
t194 = -t271 * t199 + t273 * t207;
t1 = [(t377 * t272 + t280 * t274) * MDP(19) + (t280 * t272 - t377 * t274) * MDP(18) + (-g(2) * t203 + (-t219 * t343 - t196) * t272 + (-qJD(5) * t194 + t201 * t219 + t224 * t347) * t274 + (-g(1) * t366 + (-t198 * t274 + t314 * t272) * qJD(6) + t298) * t273 + t378 * t271) * MDP(25) + (-t219 * t225 * t341 - g(2) * t202 + (qJD(5) * t195 - t200 * t219 + t225 * t347) * t274 + (t198 * t338 + (-t314 * qJD(6) + t197 + t372) * t272 - t298) * t271 + t378 * t273) * MDP(26) + t291 * MDP(3) + (-qJDD(1) * t230 + 0.2e1 * t244 + t292) * MDP(7) + 0.2e1 * (-t272 * t328 + t351 * t333) * MDP(14) + t352 * MDP(2) + (-t311 * pkin(1) - g(1) * t308 - g(2) * t353 + (t324 + t335) * qJ(2)) * MDP(6) + (qJDD(1) * t231 + t304 + 0.2e1 * t320 + t356) * MDP(8) + (qJDD(1) * t265 + t272 * t312) * MDP(13) + ((t224 * t273 + t225 * t271) * t341 + (t367 - t201 * t273 + (t224 * t271 - t225 * t273) * qJD(6)) * t274) * MDP(21) + (-qJDD(2) + t352 + 0.2e1 * t369) * MDP(4) + qJDD(1) * MDP(1) + ((-t200 - t321) * t272 + (-t295 + t342) * t274) * MDP(22) + ((-t201 + t322) * t272 + (-t296 - t343) * t274) * MDP(23) + ((-t241 - t346) * qJD(1) + t289 + t309) * MDP(11) + (t208 * t223 + t214 * t241 + t209 * t226 + t213 * t347 - g(1) * (t222 * pkin(3) + t221 * qJ(4) + t288) - g(2) * (-pkin(3) * t221 + qJ(4) * t222 + t323)) * MDP(12) + (-t291 + t324 + 0.2e1 * t335) * MDP(5) + (-t200 * t273 * t274 + (-t271 * t338 - t272 * t336) * t225) * MDP(20) + (t220 * t272 + t242 * t340) * MDP(24) + (-g(1) * t288 - g(2) * t323 - t297 * qJD(2) + t211 * t230 + t212 * t231) * MDP(9) + (-qJDD(1) * t226 - 0.2e1 * t244 - t287) * MDP(10) + t234 * MDP(16) + t235 * MDP(15); -qJDD(1) * MDP(4) - t277 * MDP(5) + (-qJ(2) * t277 + t311 - t352) * MDP(6) + (t297 * qJD(1) + t211 * t269 + t212 * t268 - t352) * MDP(9) + (t208 * t268 - t209 * t269 + (-t213 * t268 - t214 * t269) * qJD(1) - t352) * MDP(12) + ((t312 - t329) * t268 + (t350 * t272 - t325) * t269) * MDP(18) + ((0.2e1 * t319 - t328) * t268 + (t350 * t274 + t326) * t269) * MDP(19) + (t295 * t268 + (-t286 * t272 - t302 * t274) * t269 + ((-t268 * t360 + t269 * t273) * t242 - t224 * t362) * qJD(1)) * MDP(25) + (t296 * t268 + (t285 * t272 - t303 * t274) * t269 + (-(t268 * t359 + t269 * t271) * t242 - t225 * t362) * qJD(1)) * MDP(26) + (-MDP(7) + MDP(10)) * (t268 * t277 + t330) + (MDP(8) - MDP(11)) * (-t269 * t277 + t331); -t234 * MDP(18) + t235 * MDP(19) + (MDP(9) + MDP(12)) * t357 + t376 * t272 + t278 * t274; -qJDD(1) * MDP(10) - t277 * MDP(11) + (t244 + t287) * MDP(12) + (t214 * MDP(12) + (-MDP(25) * t273 + MDP(26) * t271) * t242) * qJD(1) + (qJDD(5) * MDP(18) - t350 * MDP(19) - t376) * t274 + (-t350 * MDP(18) - qJDD(5) * MDP(19) + t278) * t272; -MDP(15) * t328 + MDP(16) * t329 + qJDD(5) * MDP(17) + (t301 * t272 - t281 * t274 + t344) * MDP(18) + (t281 * t272 + t301 * t274 + t345) * MDP(19) + (t225 * t358 + t367) * MDP(20) + ((t200 - t365) * t273 + (t201 - t364) * t271) * MDP(21) + ((-t225 * t274 + t272 * t358) * qJD(1) - t296) * MDP(22) + ((t224 * t274 - t242 * t360) * qJD(1) + t295) * MDP(23) - t242 * MDP(24) * t348 + (pkin(5) * t201 + t194 * t348 + t205 * t224 + t282 * t271 + t375 * t273) * MDP(25) + (-pkin(5) * t200 - t195 * t348 + t205 * t225 - t375 * t271 + t282 * t273) * MDP(26) + (t274 * t272 * MDP(13) - t351 * MDP(14)) * t277; t225 * t224 * MDP(20) + (-t224 ^ 2 + t225 ^ 2) * MDP(21) + (t310 + t365) * MDP(22) + (-t355 + t364) * MDP(23) - t220 * MDP(24) + (-g(1) * t202 - t195 * t242 + t198 * t225 + t196) * MDP(25) + (g(1) * t203 - t194 * t242 - t198 * t224) * MDP(26) + (-MDP(22) * t328 + t299 * MDP(23) + t300 * MDP(25) + t284 * MDP(26)) * t273 + (MDP(23) * t328 + t284 * MDP(25) + (-t197 - t300) * MDP(26)) * t271;];
tau  = t1;
