% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:46:12
% EndTime: 2021-01-15 12:46:19
% DurationCPUTime: 2.60s
% Computational Cost: add. (1855->279), mult. (3844->350), div. (0->0), fcn. (2524->12), ass. (0->146)
t301 = sin(qJ(3));
t303 = cos(qJ(3));
t298 = sin(pkin(8));
t276 = pkin(1) * t298 + pkin(6);
t378 = pkin(7) + t276;
t338 = t378 * qJD(1);
t229 = qJD(2) * t301 + t303 * t338;
t300 = sin(qJ(4));
t223 = t300 * t229;
t228 = t303 * qJD(2) - t338 * t301;
t377 = qJD(3) * pkin(3);
t226 = t228 + t377;
t382 = cos(qJ(4));
t336 = t226 * t382 - t223;
t251 = t300 * t303 + t301 * t382;
t245 = t251 * qJD(1);
t372 = t245 * qJ(5);
t199 = -t372 + t336;
t293 = qJD(3) + qJD(4);
t286 = t303 * qJDD(2);
t261 = t276 * qJDD(1);
t337 = pkin(7) * qJDD(1) + t261;
t207 = qJDD(3) * pkin(3) - qJD(3) * t229 - t301 * t337 + t286;
t208 = qJD(3) * t228 + t301 * qJDD(2) + t303 * t337;
t390 = t207 * t382 - t300 * t208;
t248 = t378 * t301;
t249 = t378 * t303;
t360 = -t248 * t300 + t249 * t382;
t291 = qJDD(3) + qJDD(4);
t343 = t382 * qJD(4);
t381 = pkin(3) * t293;
t389 = -pkin(3) * t291 * t300 - t343 * t381;
t299 = cos(pkin(8));
t277 = -pkin(1) * t299 - pkin(2);
t290 = t303 * pkin(3);
t388 = t277 - t290;
t283 = t291 * pkin(4);
t365 = t300 * t301;
t327 = t293 * t365;
t345 = t382 * t303;
t331 = qJD(1) * t345;
t340 = qJDD(1) * t382;
t349 = qJDD(1) * t303;
t332 = t293 * t331 + t300 * t349 + t301 * t340;
t211 = qJD(1) * t327 - t332;
t376 = t211 * qJ(5);
t387 = t283 + t376;
t386 = t382 * qJD(3) + t343;
t221 = t293 * t251;
t271 = t303 * t340;
t350 = qJDD(1) * t301;
t326 = t300 * t350 - t271;
t212 = qJD(1) * t221 + t326;
t385 = pkin(4) * t212 + qJDD(5);
t294 = qJ(1) + pkin(8);
t285 = cos(t294);
t297 = qJ(3) + qJ(4);
t288 = sin(t297);
t368 = t285 * t288;
t284 = sin(t294);
t370 = t284 * t288;
t289 = cos(t297);
t379 = g(1) * t289;
t384 = g(2) * t370 - g(3) * t368 - t379;
t383 = t245 ^ 2;
t375 = t212 * qJ(5);
t355 = qJD(1) * t301;
t243 = t300 * t355 - t331;
t374 = t243 * qJ(5);
t373 = t243 * t293;
t369 = t284 * t289;
t367 = t285 * t289;
t364 = qJDD(2) - g(1);
t197 = pkin(4) * t293 + t199;
t363 = t197 - t199;
t220 = -t303 * t386 + t327;
t362 = -t212 * t251 + t220 * t243;
t361 = t228 * t382 - t223;
t359 = g(2) * t368 + g(3) * t370;
t358 = pkin(4) * t289 + t290;
t295 = t301 ^ 2;
t357 = -t303 ^ 2 + t295;
t356 = MDP(19) * t301;
t264 = qJD(1) * t277;
t353 = qJD(4) * t300;
t247 = t388 * qJD(1);
t339 = pkin(4) * t243 + qJD(5);
t218 = t247 + t339;
t352 = qJD(5) + t218;
t351 = qJD(1) * qJD(3);
t347 = pkin(3) * t355;
t346 = t301 * t377;
t225 = t382 * t229;
t342 = t301 * t351;
t341 = qJD(3) * t378;
t335 = -t228 * t300 - t225;
t334 = -t248 * t382 - t249 * t300;
t333 = t293 * t301;
t330 = g(2) * t285 + g(3) * t284;
t329 = g(2) * t284 - g(3) * t285;
t302 = sin(qJ(1));
t304 = cos(qJ(1));
t328 = -g(2) * t304 - g(3) * t302;
t250 = -t345 + t365;
t324 = -t211 * t250 + t221 * t245;
t323 = t220 * t293 - t251 * t291;
t321 = -t226 * t300 - t225;
t240 = t301 * t341;
t241 = t303 * t341;
t320 = -t240 * t382 - t241 * t300 - t248 * t343 - t249 * t353;
t242 = t243 ^ 2;
t319 = t245 * t243 * MDP(12) + (-qJD(1) * t300 * t333 + t332 + t373) * MDP(14) - t326 * MDP(15) + (-t242 + t383) * MDP(13) + t291 * MDP(16);
t318 = -qJD(1) * t264 - t261 + t329;
t275 = pkin(3) * t342;
t231 = qJDD(1) * t388 + t275;
t317 = 0.2e1 * qJD(3) * t264 - qJDD(3) * t276;
t305 = qJD(3) ^ 2;
t316 = 0.2e1 * qJDD(1) * t277 + t276 * t305 + t330;
t315 = qJD(4) * t321 + t390;
t314 = -qJD(4) * t360 + t240 * t300 - t241 * t382;
t313 = t300 * t207 + t208 * t382 + t226 * t343 - t229 * t353;
t312 = g(1) * t288 + g(2) * t369 - g(3) * t367 - t313;
t311 = t315 + t384;
t310 = t247 * t243 + t312;
t309 = -t247 * t245 + t311;
t308 = t243 * t352 + t312 + t375;
t292 = qJ(5) + pkin(7) + pkin(6);
t282 = pkin(3) * t382 + pkin(4);
t260 = qJDD(3) * t303 - t301 * t305;
t259 = qJDD(3) * t301 + t303 * t305;
t252 = pkin(2) + t358;
t230 = pkin(4) * t245 + t347;
t227 = pkin(4) * t250 + t388;
t217 = pkin(4) * t221 + t346;
t214 = -qJ(5) * t250 + t360;
t213 = -qJ(5) * t251 + t334;
t210 = -t221 * t293 - t250 * t291;
t202 = -t372 + t361;
t201 = t335 + t374;
t200 = -t321 - t374;
t198 = t231 + t385;
t196 = qJ(5) * t220 - qJD(5) * t251 + t314;
t195 = -qJ(5) * t221 - qJD(5) * t250 + t320;
t194 = -qJD(5) * t243 + t313 - t375;
t193 = -qJD(5) * t245 + t315 + t387;
t1 = [qJDD(1) * MDP(1) + t328 * MDP(2) + (g(2) * t302 - g(3) * t304) * MDP(3) + (t328 + (t298 ^ 2 + t299 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t295 + 0.2e1 * t303 * t342) * MDP(5) + 0.2e1 * (t301 * t349 - t351 * t357) * MDP(6) + t259 * MDP(7) + t260 * MDP(8) + (t301 * t317 - t303 * t316) * MDP(10) + (t301 * t316 + t303 * t317) * MDP(11) + (-t211 * t251 - t220 * t245) * MDP(12) + (-t324 + t362) * MDP(13) - t323 * MDP(14) + t210 * MDP(15) + (-g(2) * t367 - g(3) * t369 + t212 * t388 + t247 * t221 + t231 * t250 + t243 * t346 + t291 * t334 + t293 * t314) * MDP(17) + (-t211 * t388 - t247 * t220 + t231 * t251 + t245 * t346 - t291 * t360 - t293 * t320 + t359) * MDP(18) + (t196 * t293 + t198 * t250 + t212 * t227 + t213 * t291 + t217 * t243 + t218 * t221 - t289 * t330) * MDP(19) + (-t195 * t293 + t198 * t251 - t211 * t227 - t214 * t291 + t217 * t245 - t218 * t220 + t359) * MDP(20) + (-t193 * t251 - t194 * t250 - t195 * t243 - t196 * t245 + t197 * t220 - t200 * t221 + t211 * t213 - t212 * t214 - t329) * MDP(21) + (t194 * t214 + t200 * t195 + t193 * t213 + t197 * t196 + t198 * t227 + t218 * t217 - g(2) * (pkin(1) * t304 + t252 * t285 + t284 * t292) - g(3) * (pkin(1) * t302 + t252 * t284 - t285 * t292)) * MDP(22); t364 * MDP(4) + t260 * MDP(10) - t259 * MDP(11) + (t324 + t362) * MDP(21) + (-t193 * t250 + t194 * t251 - t197 * t221 - t200 * t220 - g(1)) * MDP(22) + (MDP(17) + MDP(19)) * t210 + (MDP(18) + MDP(20)) * t323; MDP(7) * t350 + MDP(8) * t349 + qJDD(3) * MDP(9) + (-g(1) * t303 + t301 * t318 + t286) * MDP(10) + (-t301 * t364 + t303 * t318) * MDP(11) + (-t335 * t293 + (-t243 * t355 + t291 * t382 - t293 * t353) * pkin(3) + t309) * MDP(17) + (-t245 * t347 + t293 * t361 + t310 + t389) * MDP(18) + (-t201 * t293 - t230 * t243 + t282 * t291 - t352 * t245 + (-t225 + (-t226 - t381) * t300) * qJD(4) + t384 + t387 + t390) * MDP(19) + (t202 * t293 - t230 * t245 + t308 + t389) * MDP(20) + (t282 * t211 + (t200 + t201) * t245 + (-t197 + t202) * t243 + (-t212 * t300 + (-t243 * t382 + t245 * t300) * qJD(4)) * pkin(3)) * MDP(21) + (t193 * t282 - t200 * t202 - t197 * t201 - t218 * t230 - g(1) * t358 - t329 * (-pkin(3) * t301 - pkin(4) * t288) + (t194 * t300 + (-t197 * t300 + t200 * t382) * qJD(4)) * pkin(3)) * MDP(22) + t319 + (-MDP(5) * t301 * t303 + MDP(6) * t357) * qJD(1) ^ 2; (-t293 * t321 + t309) * MDP(17) + (t293 * t336 + t310) * MDP(18) + (t376 + t200 * t293 + 0.2e1 * t283 + (-t218 - t339) * t245 + t311) * MDP(19) + (-pkin(4) * t383 + t199 * t293 + t308) * MDP(20) + (pkin(4) * t211 - t243 * t363) * MDP(21) + (t363 * t200 + (-t218 * t245 + t288 * t329 + t193 - t379) * pkin(4)) * MDP(22) + t319; (t245 * t293 - t271) * MDP(19) + (t332 - t373) * MDP(20) + (-t242 - t383) * MDP(21) + (t197 * t245 + t200 * t243 + t275 + t330 + t385) * MDP(22) + (MDP(22) * t388 + t300 * t356) * qJDD(1) + (t386 * t356 + (MDP(19) * t293 * t303 - MDP(20) * t333) * t300) * qJD(1);];
tau = t1;
