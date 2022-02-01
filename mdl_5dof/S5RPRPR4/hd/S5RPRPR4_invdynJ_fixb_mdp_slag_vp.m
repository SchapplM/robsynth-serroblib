% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:13
% EndTime: 2022-01-23 09:23:17
% DurationCPUTime: 2.42s
% Computational Cost: add. (1688->287), mult. (3589->380), div. (0->0), fcn. (2526->16), ass. (0->140)
t310 = sin(pkin(9));
t312 = cos(pkin(9));
t319 = cos(qJ(3));
t354 = qJD(1) * t319;
t346 = t312 * t354;
t316 = sin(qJ(3));
t355 = qJD(1) * t316;
t263 = t310 * t355 - t346;
t318 = cos(qJ(5));
t255 = t318 * t263;
t272 = t310 * t319 + t312 * t316;
t266 = t272 * qJD(1);
t315 = sin(qJ(5));
t363 = t266 * t315;
t224 = -t255 - t363;
t305 = qJD(3) + qJD(5);
t364 = t224 * t305;
t333 = t263 * t315 - t266 * t318;
t365 = t333 * t305;
t311 = sin(pkin(8));
t292 = pkin(1) * t311 + pkin(6);
t359 = qJ(4) + t292;
t307 = qJ(1) + pkin(8);
t298 = sin(t307);
t300 = cos(t307);
t340 = g(1) * t300 + g(2) * t298;
t301 = t319 * qJDD(2);
t280 = t292 * qJDD(1);
t325 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t280;
t341 = t359 * qJD(1);
t332 = t341 * qJD(3);
t210 = qJDD(3) * pkin(3) - t316 * t325 - t319 * t332 + t301;
t214 = (qJDD(2) - t332) * t316 + t325 * t319;
t193 = t210 * t312 - t214 * t310;
t352 = qJD(1) * qJD(3);
t345 = t316 * t352;
t284 = t310 * t345;
t344 = t319 * t352;
t232 = t272 * qJDD(1) + t312 * t344 - t284;
t191 = qJDD(3) * pkin(4) - pkin(7) * t232 + t193;
t194 = t210 * t310 + t214 * t312;
t265 = t272 * qJD(3);
t350 = qJDD(1) * t319;
t285 = t312 * t350;
t351 = qJDD(1) * t316;
t231 = qJD(1) * t265 + t310 * t351 - t285;
t192 = -pkin(7) * t231 + t194;
t248 = qJD(2) * t319 - t316 * t341;
t366 = qJD(3) * pkin(3);
t242 = t248 + t366;
t249 = qJD(2) * t316 + t319 * t341;
t361 = t312 * t249;
t213 = t242 * t310 + t361;
t371 = pkin(7) * t263;
t200 = t213 - t371;
t296 = pkin(3) * t319 + pkin(2);
t313 = cos(pkin(8));
t374 = pkin(1) * t313;
t277 = -t296 - t374;
t261 = qJD(1) * t277 + qJD(4);
t233 = pkin(4) * t263 + t261;
t306 = qJ(3) + pkin(9);
t302 = qJ(5) + t306;
t290 = sin(t302);
t291 = cos(t302);
t353 = qJD(5) * t315;
t377 = g(3) * t290 - t315 * t191 - t318 * t192 + t200 * t353 - t233 * t224 + t291 * t340;
t376 = -g(3) * t291 + t318 * t191 - t315 * t192 + t233 * t333 + t290 * t340;
t304 = qJDD(3) + qJDD(5);
t375 = t304 * MDP(20) + t224 * MDP(16) * t333 + (-t224 ^ 2 + t333 ^ 2) * MDP(17);
t343 = t231 * t318 + t232 * t315;
t197 = -qJD(5) * t333 + t343;
t373 = pkin(3) * t310;
t372 = pkin(3) * t316;
t370 = pkin(7) * t266;
t367 = g(3) * t319;
t238 = t310 * t249;
t362 = t310 * t316;
t212 = t242 * t312 - t238;
t199 = qJD(3) * pkin(4) + t212 - t370;
t360 = t318 * t199;
t358 = qJDD(2) - g(3);
t216 = t248 * t312 - t238;
t342 = qJD(3) * t359;
t252 = qJD(4) * t319 - t316 * t342;
t253 = -qJD(4) * t316 - t319 * t342;
t218 = t252 * t312 + t253 * t310;
t269 = t359 * t316;
t270 = t359 * t319;
t230 = -t269 * t310 + t270 * t312;
t308 = t316 ^ 2;
t357 = -t319 ^ 2 + t308;
t294 = -pkin(2) - t374;
t283 = qJD(1) * t294;
t349 = pkin(3) * t345 + qJDD(4);
t348 = t316 * t366;
t347 = -qJD(5) * t255 - t231 * t315 + t232 * t318;
t215 = -t248 * t310 - t361;
t217 = -t252 * t310 + t253 * t312;
t229 = -t269 * t312 - t270 * t310;
t339 = g(1) * t298 - g(2) * t300;
t317 = sin(qJ(1));
t320 = cos(qJ(1));
t338 = g(1) * t317 - g(2) * t320;
t337 = -t315 * t199 - t318 * t200;
t271 = -t312 * t319 + t362;
t234 = t271 * t318 + t272 * t315;
t268 = t271 * qJD(3);
t201 = -qJD(5) * t234 - t265 * t315 - t268 * t318;
t235 = -t271 * t315 + t272 * t318;
t336 = t201 * t305 + t235 * t304;
t219 = -pkin(7) * t272 + t229;
t220 = -pkin(7) * t271 + t230;
t335 = t219 * t318 - t220 * t315;
t334 = t219 * t315 + t220 * t318;
t293 = pkin(3) * t312 + pkin(4);
t331 = t293 * t315 + t318 * t373;
t330 = t293 * t318 - t315 * t373;
t196 = -t266 * t353 + t347;
t328 = -qJD(1) * t283 - t280 + t340;
t327 = 0.2e1 * qJD(3) * t283 - qJDD(3) * t292;
t247 = qJDD(1) * t277 + t349;
t321 = qJD(3) ^ 2;
t324 = -0.2e1 * qJDD(1) * t294 - t292 * t321 + t339;
t314 = -qJ(4) - pkin(6);
t299 = cos(t306);
t297 = sin(t306);
t279 = qJDD(3) * t319 - t316 * t321;
t278 = qJDD(3) * t316 + t319 * t321;
t251 = pkin(4) * t265 + t348;
t250 = pkin(3) * t355 + pkin(4) * t266;
t246 = pkin(4) * t271 + t277;
t211 = pkin(4) * t231 + t247;
t206 = -pkin(7) * t265 + t218;
t205 = pkin(7) * t268 + t217;
t204 = t216 - t370;
t203 = t215 + t371;
t202 = qJD(5) * t235 + t265 * t318 - t268 * t315;
t195 = -t202 * t305 - t234 * t304;
t1 = [qJDD(1) * MDP(1) + t338 * MDP(2) + (g(1) * t320 + g(2) * t317) * MDP(3) + (t338 + (t311 ^ 2 + t313 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t308 + 0.2e1 * t316 * t344) * MDP(5) + 0.2e1 * (t316 * t350 - t352 * t357) * MDP(6) + t278 * MDP(7) + t279 * MDP(8) + (t316 * t327 + t319 * t324) * MDP(10) + (-t316 * t324 + t319 * t327) * MDP(11) + (qJDD(3) * t229 + t231 * t277 + t247 * t271 + t261 * t265 + t339 * t299 + (t263 * t372 + t217) * qJD(3)) * MDP(12) + (-qJDD(3) * t230 + t232 * t277 + t247 * t272 - t261 * t268 - t339 * t297 + (t266 * t372 - t218) * qJD(3)) * MDP(13) + (-t193 * t272 - t194 * t271 + t212 * t268 - t213 * t265 - t217 * t266 - t218 * t263 - t229 * t232 - t230 * t231 - t340) * MDP(14) + (t194 * t230 + t213 * t218 + t193 * t229 + t212 * t217 + t247 * t277 + t261 * t348 - g(1) * (-pkin(1) * t317 - t296 * t298 - t300 * t314) - g(2) * (pkin(1) * t320 + t296 * t300 - t298 * t314)) * MDP(15) + (t196 * t235 - t201 * t333) * MDP(16) + (-t196 * t234 - t197 * t235 + t201 * t224 + t202 * t333) * MDP(17) + t336 * MDP(18) + t195 * MDP(19) + (-t251 * t224 + t246 * t197 + t211 * t234 + t233 * t202 + (-qJD(5) * t334 + t205 * t318 - t206 * t315) * t305 + t335 * t304 + t339 * t291) * MDP(21) + (-t251 * t333 + t246 * t196 + t211 * t235 + t233 * t201 - (qJD(5) * t335 + t205 * t315 + t206 * t318) * t305 - t334 * t304 - t339 * t290) * MDP(22); t358 * MDP(4) + t279 * MDP(10) - t278 * MDP(11) + (-qJD(3) * t265 - qJDD(3) * t271) * MDP(12) + (qJD(3) * t268 - qJDD(3) * t272) * MDP(13) + (-t231 * t272 + t232 * t271 + t263 * t268 + t265 * t266) * MDP(14) + (-t193 * t271 + t194 * t272 - t212 * t265 - t213 * t268 - g(3)) * MDP(15) + t195 * MDP(21) - t336 * MDP(22); MDP(7) * t351 + MDP(8) * t350 + qJDD(3) * MDP(9) + (t316 * t328 + t301 - t367) * MDP(10) + (-t316 * t358 + t319 * t328) * MDP(11) + (-g(3) * t299 - qJD(3) * t215 - t261 * t266 + t340 * t297 + (qJDD(3) * t312 - t263 * t355) * pkin(3) + t193) * MDP(12) + (g(3) * t297 + qJD(3) * t216 + t261 * t263 + t340 * t299 + (-qJDD(3) * t310 - t266 * t355) * pkin(3) - t194) * MDP(13) + ((t213 + t215) * t266 + (-t212 + t216) * t263 + (-t231 * t310 - t232 * t312) * pkin(3)) * MDP(14) + (-t212 * t215 - t213 * t216 + (-t367 + t193 * t312 + t194 * t310 + (-qJD(1) * t261 + t340) * t316) * pkin(3)) * MDP(15) + (t196 - t364) * MDP(18) + (-t197 - t365) * MDP(19) + (t330 * t304 + t250 * t224 - (t203 * t318 - t204 * t315) * t305 + (-t305 * t331 + t337) * qJD(5) + t376) * MDP(21) + (-t331 * t304 + t250 * t333 + (t203 * t315 + t204 * t318) * t305 + (-t305 * t330 - t360) * qJD(5) + t377) * MDP(22) + (-MDP(5) * t316 * t319 + MDP(6) * t357) * qJD(1) ^ 2 + t375; -t285 * MDP(12) - t284 * MDP(13) + (-t263 ^ 2 - t266 ^ 2) * MDP(14) + (t212 * t266 + t213 * t263 - t339 + t349) * MDP(15) + (t197 - t365) * MDP(21) + (t196 + t364) * MDP(22) + (MDP(12) * t362 + MDP(13) * t272 + MDP(15) * t277) * qJDD(1) + ((t310 * t354 + t312 * t355 + t266) * MDP(12) + (-t263 + t346) * MDP(13)) * qJD(3); (t347 - t364) * MDP(18) + (-t343 - t365) * MDP(19) + (-t305 * t337 + t376) * MDP(21) + ((-t200 * t315 + t360) * t305 + t377) * MDP(22) + (-MDP(18) * t363 + t333 * MDP(19) + t337 * MDP(21) - MDP(22) * t360) * qJD(5) + t375;];
tau = t1;
