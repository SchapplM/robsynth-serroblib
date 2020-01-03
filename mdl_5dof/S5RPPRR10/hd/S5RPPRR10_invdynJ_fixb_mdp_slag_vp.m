% Calculate vector of inverse dynamics joint torques for
% S5RPPRR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:23
% EndTime: 2019-12-31 18:04:27
% DurationCPUTime: 2.37s
% Computational Cost: add. (1158->283), mult. (2635->375), div. (0->0), fcn. (2029->10), ass. (0->138)
t307 = qJD(4) + qJD(5);
t391 = t307 ^ 2;
t310 = sin(pkin(8));
t311 = cos(pkin(8));
t314 = sin(qJ(4));
t317 = cos(qJ(4));
t262 = t310 * t314 + t311 * t317;
t252 = t262 * qJD(1);
t316 = cos(qJ(5));
t369 = qJD(1) * t311;
t350 = t314 * t369;
t370 = qJD(1) * t310;
t352 = t317 * t370;
t254 = -t350 + t352;
t313 = sin(qJ(5));
t378 = t254 * t313;
t214 = t316 * t252 + t378;
t387 = t214 * t307;
t271 = qJ(2) * t370 + qJD(3);
t258 = -pkin(6) * t370 + t271;
t381 = -pkin(6) + qJ(2);
t269 = t381 * t311;
t264 = qJD(1) * t269;
t337 = -t258 * t314 - t264 * t317;
t209 = -pkin(7) * t252 - t337;
t246 = -qJD(1) * pkin(1) - pkin(2) * t369 - qJ(3) * t370 + qJD(2);
t230 = pkin(3) * t369 - t246;
t212 = pkin(4) * t252 + t230;
t315 = sin(qJ(1));
t308 = qJ(4) + qJ(5);
t296 = sin(t308);
t297 = cos(t308);
t336 = t296 * t310 + t297 * t311;
t232 = t336 * t315;
t318 = cos(qJ(1));
t234 = t336 * t318;
t376 = t297 * t310;
t248 = t296 * t311 - t376;
t363 = qJD(5) * t313;
t390 = g(1) * t234 + g(2) * t232 - g(3) * t248 + t209 * t363 + t212 * t214;
t304 = qJDD(4) + qJDD(5);
t338 = -t252 * t313 + t316 * t254;
t389 = t304 * MDP(23) + t214 * t338 * MDP(19) + (-t214 ^ 2 + t338 ^ 2) * MDP(20);
t386 = t338 * t307;
t385 = t317 * t258 - t264 * t314;
t268 = t381 * t310;
t373 = t314 * t268 + t317 * t269;
t303 = g(2) * t318;
t382 = g(1) * t315;
t348 = -t303 + t382;
t358 = qJD(1) * qJD(3);
t284 = t310 * t358;
t356 = qJDD(1) * t311;
t357 = qJDD(1) * t310;
t384 = -pkin(2) * t356 - qJ(3) * t357 - t284;
t250 = t262 * qJD(4);
t343 = -t314 * t356 + t317 * t357;
t221 = -qJD(1) * t250 + t343;
t359 = qJD(1) * qJD(2);
t256 = qJ(2) * t357 + t310 * t359 + qJDD(3);
t235 = -pkin(6) * t357 + t256;
t237 = (t381 * qJDD(1) + t359) * t311;
t346 = t317 * t235 - t314 * t237;
t197 = qJDD(4) * pkin(4) - pkin(7) * t221 + t337 * qJD(4) + t346;
t325 = qJD(4) * t350 - t262 * qJDD(1);
t364 = qJD(4) * t317;
t351 = t310 * t364;
t222 = qJD(1) * t351 - t325;
t339 = t314 * t235 + t317 * t237;
t198 = -pkin(7) * t222 + t385 * qJD(4) + t339;
t375 = t311 * t315;
t231 = t296 * t375 - t315 * t376;
t233 = t248 * t318;
t383 = g(1) * t233 + g(2) * t231 + g(3) * t336 + t316 * t197 - t313 * t198 - t212 * t338;
t347 = t313 * t221 + t316 * t222;
t200 = t338 * qJD(5) + t347;
t300 = t311 * pkin(2);
t309 = qJDD(1) * pkin(1);
t208 = -pkin(7) * t254 + t385;
t207 = qJD(4) * pkin(4) + t208;
t380 = t207 * t316;
t379 = t209 * t316;
t306 = t311 ^ 2;
t354 = 0.2e1 * t359;
t280 = t306 * t354;
t355 = qJDD(1) * qJ(2) ^ 2;
t372 = qJ(2) * t280 + t306 * t355;
t371 = t318 * pkin(1) + t315 * qJ(2);
t368 = qJD(2) * t314;
t367 = qJD(2) * t317;
t366 = qJD(3) * t310;
t365 = qJD(4) * t314;
t362 = qJD(5) * t316;
t360 = qJ(2) * qJDD(1);
t294 = qJDD(2) - t309;
t265 = -t310 * qJ(3) - pkin(1) - t300;
t353 = t316 * t221 - t313 * t222 - t252 * t362;
t349 = -pkin(1) * t315 + t318 * qJ(2);
t345 = t317 * t268 - t269 * t314;
t255 = t311 * pkin(3) - t265;
t344 = g(1) * t318 + g(2) * t315;
t342 = -t207 * t313 - t379;
t263 = t310 * t317 - t311 * t314;
t210 = -pkin(7) * t263 + t345;
t211 = -pkin(7) * t262 + t373;
t341 = t210 * t316 - t211 * t313;
t340 = t210 * t313 + t211 * t316;
t223 = t316 * t262 + t263 * t313;
t224 = -t262 * t313 + t263 * t316;
t335 = t313 * t317 + t314 * t316;
t334 = -t313 * t314 + t316 * t317;
t228 = t294 + t384;
t333 = -t348 + t294;
t332 = -t294 + t309 - t303;
t331 = -qJDD(1) * t265 - t228 - t303;
t225 = pkin(3) * t356 - t228;
t330 = 0.2e1 * t306 * t360 + t280 - t344;
t329 = t268 * t364 - t269 * t365 + t310 * t368 + t311 * t367;
t328 = t254 * t363 - t353;
t305 = t310 ^ 2;
t327 = (t359 + t360) * t305;
t322 = -t373 * qJD(4) + t310 * t367 - t311 * t368;
t320 = qJD(1) ^ 2;
t319 = qJD(4) ^ 2;
t281 = g(1) * t375;
t251 = -t311 * t365 + t351;
t244 = t262 * t318;
t243 = t263 * t318;
t242 = t262 * t315;
t241 = t263 * t315;
t229 = pkin(4) * t251 + t366;
t226 = pkin(4) * t262 + t255;
t205 = pkin(7) * t250 + t322;
t204 = -pkin(7) * t251 + t329;
t203 = pkin(4) * t222 + t225;
t202 = t224 * qJD(5) - t250 * t313 + t316 * t251;
t201 = -t223 * qJD(5) - t250 * t316 - t251 * t313;
t1 = [qJDD(1) * MDP(1) + t348 * MDP(2) + t344 * MDP(3) + (t332 * t311 + t281) * MDP(4) + (-t332 - t382) * t310 * MDP(5) + (0.2e1 * t327 + t330) * MDP(6) + (-t294 * pkin(1) - g(1) * t349 - g(2) * t371 + (qJ(2) * t354 + t355) * t305 + t372) * MDP(7) + (t281 + (t331 + t284) * t311) * MDP(8) + (t256 * t310 + t327 + t330) * MDP(9) + (t305 * t358 + (t331 + t382) * t310) * MDP(10) + (t228 * t265 - g(1) * (-pkin(2) * t375 + t349) - g(2) * (t318 * t300 + t371) + (t256 * qJ(2) + t348 * qJ(3) + t271 * qJD(2) - t246 * qJD(3)) * t310 + t372) * MDP(11) + (t221 * t263 - t250 * t254) * MDP(12) + (-t221 * t262 - t222 * t263 + t250 * t252 - t251 * t254) * MDP(13) + (-qJD(4) * t250 + qJDD(4) * t263) * MDP(14) + (-qJD(4) * t251 - qJDD(4) * t262) * MDP(15) + (g(1) * t242 - g(2) * t244 + t322 * qJD(4) + t345 * qJDD(4) + t255 * t222 + t225 * t262 + t230 * t251 + t252 * t366) * MDP(17) + (g(1) * t241 - g(2) * t243 - t329 * qJD(4) - t373 * qJDD(4) + t255 * t221 + t225 * t263 - t230 * t250 + t254 * t366) * MDP(18) + (t201 * t338 - t224 * t328) * MDP(19) + (-t200 * t224 - t201 * t214 - t202 * t338 + t223 * t328) * MDP(20) + (t201 * t307 + t224 * t304) * MDP(21) + (-t202 * t307 - t223 * t304) * MDP(22) + (t229 * t214 + t226 * t200 + t203 * t223 + t212 * t202 + (-t340 * qJD(5) - t204 * t313 + t205 * t316) * t307 + t341 * t304 + g(1) * t232 - g(2) * t234) * MDP(24) + (t229 * t338 - t226 * t328 + t203 * t224 + t212 * t201 - (t341 * qJD(5) + t204 * t316 + t205 * t313) * t307 - t340 * t304 - g(1) * t231 + g(2) * t233) * MDP(25); t333 * MDP(7) + (-qJ(2) * t306 * t320 - t271 * t370 + t333 + t384) * MDP(11) + ((-t254 - t352) * qJD(4) + t325) * MDP(17) + (0.2e1 * qJD(4) * t252 - t343) * MDP(18) + (-t200 - t386) * MDP(24) + (t328 + t387) * MDP(25) + (MDP(5) - MDP(10)) * t357 + (-MDP(4) - MDP(8)) * t356 + (qJ(2) * MDP(7) + MDP(6) + MDP(9)) * (-t305 - t306) * t320; -t305 * t320 * MDP(10) + (g(3) * t311 + t256) * MDP(11) + (qJDD(4) * t317 - t314 * t319) * MDP(17) + (-t314 * qJDD(4) - t317 * t319) * MDP(18) + (t334 * t304 - t391 * t335) * MDP(24) + (-t335 * t304 - t391 * t334) * MDP(25) + (-t320 * t311 * MDP(8) + qJDD(1) * MDP(9) - t344 * MDP(11) + (t246 * MDP(11) - t252 * MDP(17) - t254 * MDP(18) - t214 * MDP(24) - MDP(25) * t338) * qJD(1)) * t310; t254 * t252 * MDP(12) + (-t252 ^ 2 + t254 ^ 2) * MDP(13) + t343 * MDP(14) + ((t254 - t352) * qJD(4) + t325) * MDP(15) + qJDD(4) * MDP(16) + (-g(1) * t243 - g(2) * t241 + g(3) * t262 - t230 * t254 + t346) * MDP(17) + (g(1) * t244 + g(2) * t242 + g(3) * t263 + t230 * t252 - t339) * MDP(18) + (-t328 + t387) * MDP(21) + (-t200 + t386) * MDP(22) + (-(-t208 * t313 - t379) * t307 + t342 * qJD(5) + (-t214 * t254 + t316 * t304 - t307 * t363) * pkin(4) + t383) * MDP(24) + ((-t209 * t307 - t197) * t313 + (-qJD(5) * t207 + t208 * t307 - t198) * t316 + (-t254 * t338 - t313 * t304 - t307 * t362) * pkin(4) + t390) * MDP(25) + t389; (t353 + t387) * MDP(21) + (-t347 + t386) * MDP(22) + (-t342 * t307 + t383) * MDP(24) + (-t316 * t198 - t313 * t197 + (-t209 * t313 + t380) * t307 + t390) * MDP(25) + (-MDP(21) * t378 - t338 * MDP(22) + t342 * MDP(24) - MDP(25) * t380) * qJD(5) + t389;];
tau = t1;
