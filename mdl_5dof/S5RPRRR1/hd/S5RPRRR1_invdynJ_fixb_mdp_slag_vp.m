% Calculate vector of inverse dynamics joint torques for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:12
% EndTime: 2019-07-18 13:26:17
% DurationCPUTime: 3.29s
% Computational Cost: add. (1235->375), mult. (3033->524), div. (0->0), fcn. (2308->8), ass. (0->160)
t260 = sin(qJ(1));
t264 = cos(qJ(1));
t380 = g(1) * t260 - g(2) * t264;
t262 = cos(qJ(4));
t258 = sin(qJ(4));
t341 = qJD(3) * t258;
t259 = sin(qJ(3));
t343 = qJD(1) * t259;
t232 = t262 * t343 + t341;
t263 = cos(qJ(3));
t342 = qJD(1) * t263;
t245 = -qJD(4) + t342;
t257 = sin(qJ(5));
t261 = cos(qJ(5));
t209 = t232 * t257 + t261 * t245;
t252 = t262 * qJD(3);
t230 = t258 * t343 - t252;
t225 = qJD(5) + t230;
t379 = t209 * t225;
t291 = -qJDD(2) + t380;
t251 = t263 * qJDD(1);
t326 = qJD(1) * qJD(3);
t303 = t259 * t326;
t378 = -t303 + t251;
t374 = g(2) * t260;
t376 = g(1) * t264;
t377 = -t376 - t374;
t350 = t264 * t258;
t353 = t262 * t263;
t222 = t260 * t353 - t350;
t375 = g(2) * t222;
t372 = g(3) * t259;
t371 = qJ(2) * t259;
t302 = t263 * t326;
t320 = qJDD(1) * t262;
t323 = qJD(3) * qJD(4);
t296 = t258 * qJDD(3) + t259 * t320 + (t302 + t323) * t262;
t325 = qJD(1) * qJD(4);
t301 = t259 * t325;
t204 = -t258 * t301 + t296;
t370 = t204 * t258;
t369 = t225 * t257;
t229 = qJDD(4) - t378;
t368 = t229 * t263;
t367 = t230 * t245;
t366 = t230 * t259;
t365 = t232 * t259;
t364 = t232 * t262;
t249 = t262 * qJDD(3);
t298 = qJD(4) + t342;
t321 = qJDD(1) * t259;
t205 = t262 * t301 - t249 + (qJD(3) * t298 + t321) * t258;
t203 = qJDD(5) + t205;
t363 = t257 * t203;
t329 = qJD(5) * t261;
t317 = t261 * t204 + t257 * t229 - t245 * t329;
t330 = qJD(5) * t257;
t196 = -t232 * t330 + t317;
t362 = t258 * t196;
t361 = t258 * t259;
t360 = t258 * t263;
t359 = t259 * t260;
t358 = t259 * t261;
t357 = t259 * t262;
t356 = t259 * t264;
t355 = t261 * t203;
t354 = t261 * t263;
t218 = t262 * t229;
t352 = t262 * t264;
t266 = qJD(1) ^ 2;
t351 = t263 * t266;
t312 = t258 * t342;
t349 = -t245 * t312 + t218;
t348 = -t258 * qJ(2) * t303 - t262 * qJDD(2);
t255 = t259 ^ 2;
t347 = -t263 ^ 2 + t255;
t346 = MDP(20) * t262;
t211 = t232 * t261 - t245 * t257;
t345 = MDP(21) * t211;
t344 = qJ(2) * qJD(1);
t340 = qJD(3) * t259;
t339 = qJD(3) * t261;
t338 = qJD(3) * t263;
t337 = qJD(4) * t258;
t336 = qJD(4) * t261;
t335 = qJD(4) * t262;
t334 = qJD(4) * t263;
t333 = qJD(5) * t211;
t311 = t262 * t342;
t228 = qJ(2) * t311 + qJD(2) * t258;
t332 = qJD(5) * t228;
t331 = qJD(5) * t232;
t328 = t245 * MDP(19);
t327 = qJ(2) * qJDD(1);
t324 = qJD(2) * qJD(1);
t322 = qJDD(1) * t255;
t319 = qJDD(3) * t259;
t318 = qJDD(3) * t263;
t295 = qJ(2) * t302;
t304 = t259 * t324;
t316 = (qJ(2) * t321 + t295 + t304) * t261;
t315 = t225 * t344;
t314 = t211 * t342;
t313 = t230 * t343;
t310 = t245 * t337;
t309 = t259 * t337;
t308 = t258 * t336;
t307 = t261 * t343;
t306 = qJD(5) * t344;
t305 = qJ(2) * t251;
t299 = 0.2e1 * t302;
t297 = -qJD(5) + t252;
t294 = qJ(2) * t263 * t325;
t221 = t260 * t360 + t352;
t223 = t260 * t262 - t263 * t350;
t293 = -g(1) * t223 + g(2) * t221;
t214 = t257 * t311 - t307;
t289 = t257 * t335 - t214;
t283 = t257 * t259 + t261 * t353;
t215 = t283 * qJD(1);
t288 = t261 * t335 - t215;
t287 = -t324 + t374;
t286 = qJD(1) * t255 + t245 * t263;
t285 = t298 * qJD(2);
t198 = (qJDD(2) - t294) * t258 + (t378 * qJ(2) + t285) * t262;
t284 = t306 * t259 + t198;
t220 = -t257 * t263 + t261 * t357;
t282 = -t257 * t353 + t358;
t219 = t257 * t357 + t354;
t281 = t229 * t258 - t245 * t335;
t280 = -t225 * t329 - t363;
t279 = t225 * t330 - t355;
t278 = t286 * t262;
t277 = -t351 + (qJD(4) + t245) * qJD(1);
t213 = qJ(2) * t257 * t343 + t228 * t261;
t276 = t258 * t338 + t259 * t335;
t275 = t377 + 0.2e1 * t324;
t274 = t287 + t376;
t273 = -t330 * t258 + t288;
t265 = qJD(3) ^ 2;
t272 = -qJ(2) * t265 + t291;
t199 = t262 * t294 + (t285 + t305) * t258 + t348;
t271 = g(3) * t361 - t199 + t293;
t270 = -qJD(2) * qJD(4) - t305 + t372;
t269 = qJD(4) * t209 + t280;
t268 = -t275 - t327;
t267 = t322 - t368 + (-t245 + 0.2e1 * t342) * t340;
t227 = qJ(2) * t312 - t262 * qJD(2);
t224 = t258 * t260 + t263 * t352;
t217 = t261 * t229;
t212 = qJ(2) * t307 - t228 * t257;
t208 = t224 * t261 + t257 * t356;
t207 = -t224 * t257 + t261 * t356;
t206 = t211 * t337;
t201 = t297 * t354 + (-t308 + (-qJD(5) * t262 + qJD(3)) * t257) * t259;
t200 = -t257 * t309 - t263 * t330 - t259 * t339 + (t257 * t338 + t259 * t329) * t262;
t197 = t204 * t257 - t217 + t333;
t195 = -qJD(5) * t213 - t257 * t198 + t316;
t194 = t284 * t261 + (t304 - t332 + (t302 + t321) * qJ(2)) * t257;
t1 = [(-qJ(2) * t318 - t259 * t272) * MDP(13) + (-t259 * t265 + t318) * MDP(10) + (-qJ(2) * t319 + t263 * t272) * MDP(12) + (t263 * t265 + t319) * MDP(9) + (t195 * t361 + t199 * t219 + t227 * t200 - g(1) * (-t222 * t261 - t257 * t359) - g(2) * t208 + t276 * t212 + (t209 * t360 + t225 * t282) * qJD(2) + ((-t209 * t341 + t297 * t369 + t355) * t259 + ((t257 * t337 + t339) * t225 + t258 * t197 + t269 * t262) * t263) * qJ(2)) * MDP(26) + ((-t230 * t262 - t232 * t258) * t338 + (-t370 - t205 * t262 + (t230 * t258 - t364) * qJD(4)) * t259) * MDP(15) + (-t228 * t340 - g(1) * t221 - g(2) * t223 + t198 * t263 + (t278 + t365) * qJD(2) + (t204 * t259 + t232 * t338 + t262 * t267 - t286 * t337) * qJ(2)) * MDP(20) + (-t227 * t340 + g(1) * t222 - g(2) * t224 + t199 * t263 + (t258 * t286 + t366) * qJD(2) + (qJD(4) * t278 + t205 * t259 + t230 * t338 + t258 * t267) * qJ(2)) * MDP(19) + (-t245 * t340 - t368) * MDP(18) + (-t197 * t361 - t200 * t225 - t203 * t219 - t209 * t276) * MDP(24) + (t203 * t361 + t225 * t276) * MDP(25) + (t196 * t361 + t201 * t225 + t203 * t220 + t211 * t276) * MDP(23) + (-t194 * t361 + t199 * t220 + t227 * t201 - g(1) * (t222 * t257 - t260 * t358) - g(2) * t207 - t276 * t213 + (t211 * t360 - t225 * t283) * qJD(2) + ((t225 * t261 * t297 - t211 * t341 - t363) * t259 + (-(qJD(3) * t257 - t308) * t225 + t362 + (qJD(4) * t211 + t279) * t262) * t263) * qJ(2)) * MDP(27) + t380 * MDP(2) + t291 * MDP(4) + 0.2e1 * (t251 * t259 - t326 * t347) * MDP(8) + (t275 + 0.2e1 * t327) * MDP(5) + ((-t245 * t252 - t204) * t263 + (qJD(3) * t232 + t218 + t310) * t259) * MDP(16) - t268 * qJ(2) * MDP(6) + (t204 * t357 + (t252 * t263 - t309) * t232) * MDP(14) + (-t196 * t219 - t197 * t220 - t200 * t211 - t201 * t209) * MDP(22) + (t196 * t220 + t201 * t211) * MDP(21) + ((t245 * t341 + t205) * t263 + (-qJD(3) * t230 - t281) * t259) * MDP(17) + qJDD(1) * MDP(1) - t377 * MDP(3) + (t259 * t299 + t322) * MDP(7); -qJDD(1) * MDP(4) - t266 * MDP(5) + (-qJ(2) * t266 - t291) * MDP(6) + (-t251 + 0.2e1 * t303) * MDP(12) + (t299 + t321) * MDP(13) + (-t313 + t349) * MDP(19) - t232 * MDP(20) * t343 + t214 * t225 * MDP(26) + (t215 * t225 + t206) * MDP(27) + ((-qJD(4) * t369 - t197) * MDP(26) + (-t225 * t336 - t196) * MDP(27) - t245 ^ 2 * MDP(20)) * t262 + (qJD(4) * t328 - t229 * MDP(20) + (-t209 * t342 + t269) * MDP(26) + (t279 - t314) * MDP(27)) * t258; -t259 * MDP(7) * t351 + t347 * MDP(8) * t266 + MDP(9) * t321 + MDP(10) * t251 + qJDD(3) * MDP(11) + (-g(3) * t263 + t259 * t268) * MDP(12) + (t263 * t268 + t372) * MDP(13) + (-t245 * t364 + t370) * MDP(14) + ((t204 + t367) * t262 + (t232 * t245 - t205) * t258) * MDP(15) + ((t245 * t353 - t365) * qJD(1) + t281) * MDP(16) + (t310 + t313 + t349) * MDP(17) + t245 * MDP(18) * t343 + ((-g(3) * t262 + (-t230 - t252) * t344) * t263 + (qJD(1) * t227 + t274 * t262 + (t258 * t277 - t320) * qJ(2)) * t259) * MDP(19) + ((g(3) * t258 + (-t232 + t341) * t344) * t263 + (qJD(1) * t228 - t274 * t258 + (qJDD(1) * t258 + t262 * t277) * qJ(2)) * t259) * MDP(20) + (t211 * t273 + t261 * t362) * MDP(21) + (t209 * t215 + t211 * t214 + (-t209 * t261 - t211 * t257) * t335 + (-t196 * t257 - t197 * t261 + (t209 * t257 - t211 * t261) * qJD(5)) * t258) * MDP(22) + (-t196 * t262 + t206 + (-t314 + t355) * t258 + t273 * t225) * MDP(23) + (t197 * t262 - t289 * t225 + (t209 * t245 + t280) * t258) * MDP(24) + (-t225 * t245 * t258 - t262 * t203) * MDP(25) + (-t195 * t262 - g(3) * t283 + t289 * t227 - t219 * t315 + (t227 * t329 + t212 * qJD(4) + t199 * t257 + (t209 * t371 - t212 * t263) * qJD(1)) * t258 - t377 * t220) * MDP(26) + (t194 * t262 - g(3) * t282 + t288 * t227 - t220 * t315 + (-t227 * t330 - t213 * qJD(4) + t199 * t261 + (t211 * t371 + t213 * t263) * qJD(1)) * t258 + t377 * t219) * MDP(27); -t230 ^ 2 * MDP(15) + (t296 - t367) * MDP(16) + t249 * MDP(17) + t229 * MDP(18) + (t293 - t348) * MDP(19) + (g(1) * t224 + t227 * t245 + t375) * MDP(20) + t270 * t346 + (-t209 * MDP(26) - t211 * MDP(27) - t328) * t228 + ((-t321 - t323) * MDP(17) + t270 * MDP(19) - qJDD(2) * MDP(20)) * t258 + (MDP(14) * t230 + t232 * MDP(15) - t245 * MDP(17) - t211 * MDP(23) + t209 * MDP(24) - t225 * MDP(25) - t212 * MDP(26) + t213 * MDP(27)) * t232 + (-MDP(16) * t309 - t276 * MDP(17) + (-MDP(19) * t258 - t346) * t263 * qJD(2) + ((-t262 * t334 - t365) * MDP(19) + (t252 * t259 + t258 * t334 + t366) * MDP(20)) * qJ(2)) * qJD(1) + (t196 * MDP(21) + (-t211 * t230 - t197 - t333) * MDP(22) + t203 * MDP(23) - t271 * MDP(27) - t225 ^ 2 * MDP(24)) * t257 + ((t196 - t379) * MDP(22) + t203 * MDP(24) + t271 * MDP(26) + (MDP(23) * t225 + t345) * t225) * t261; t209 * t345 + (-t209 ^ 2 + t211 ^ 2) * MDP(22) + (t317 + t379) * MDP(23) + (t211 * t225 + t217) * MDP(24) + t203 * MDP(25) + (-g(1) * t207 + g(3) * t219 - t227 * t211 + t213 * t225 + t316) * MDP(26) + (g(1) * t208 + g(3) * t220 + t227 * t209 + t212 * t225) * MDP(27) + (-MDP(24) * t331 + (-g(2) * t359 - t332) * MDP(26) + (-t284 + t375) * MDP(27)) * t261 + (-MDP(23) * t331 + (qJD(5) * t245 - t204) * MDP(24) + (-t198 + t375) * MDP(26) + (-t295 + t332) * MDP(27) + (-MDP(26) * t306 + (t287 - t327) * MDP(27)) * t259) * t257;];
tau  = t1;
