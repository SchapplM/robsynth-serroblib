% Calculate vector of inverse dynamics joint torques for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:26
% EndTime: 2021-01-15 12:27:33
% DurationCPUTime: 2.79s
% Computational Cost: add. (1805->283), mult. (3389->335), div. (0->0), fcn. (2111->8), ass. (0->137)
t289 = qJ(3) + qJ(4);
t277 = sin(t289);
t278 = cos(t289);
t292 = sin(qJ(1));
t295 = cos(qJ(1));
t373 = g(1) * t292 - g(2) * t295;
t383 = g(3) * t277 - t278 * t373;
t381 = qJDD(2) - t373;
t294 = cos(qJ(3));
t341 = qJD(1) * qJD(3);
t333 = t294 * t341;
t291 = sin(qJ(3));
t340 = qJDD(1) * t291;
t380 = t333 + t340;
t339 = qJDD(1) * t294;
t379 = t291 * t341 - t339;
t284 = qJDD(1) * qJ(2);
t285 = (qJD(1) * qJD(2));
t322 = g(1) * t295 + g(2) * t292;
t310 = -t322 + (2 * t285);
t378 = 0.2e1 * t284 + t310;
t360 = qJDD(1) * pkin(1);
t377 = t360 - t381;
t296 = -pkin(1) - pkin(6);
t365 = pkin(7) - t296;
t257 = t296 * qJDD(1) + qJDD(2);
t247 = t294 * t257;
t258 = t296 * qJD(1) + qJD(2);
t346 = qJD(3) * t291;
t211 = qJDD(3) * pkin(3) + t379 * pkin(7) - t258 * t346 + t247;
t345 = qJD(3) * t294;
t213 = -t380 * pkin(7) + t257 * t291 + t258 * t345;
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t376 = t293 * t211 - t290 * t213;
t348 = qJD(1) * t291;
t335 = t290 * t348;
t347 = qJD(1) * t294;
t237 = t293 * t347 - t335;
t229 = t237 * qJ(5);
t227 = -pkin(7) * t348 + t258 * t291;
t221 = t290 * t227;
t228 = -pkin(7) * t347 + t294 * t258;
t224 = qJD(3) * pkin(3) + t228;
t329 = t293 * t224 - t221;
t201 = -t229 + t329;
t249 = t365 * t291;
t250 = t365 * t294;
t352 = -t293 * t249 - t290 * t250;
t282 = qJDD(3) + qJDD(4);
t343 = qJD(4) * t293;
t283 = qJD(3) + qJD(4);
t367 = pkin(3) * t283;
t375 = -t290 * pkin(3) * t282 - t343 * t367;
t275 = t282 * pkin(4);
t356 = t290 * t294;
t245 = t291 * t293 + t356;
t215 = t283 * t245;
t320 = t290 * t340 - t293 * t339;
t205 = t215 * qJD(1) + t320;
t363 = qJ(5) * t205;
t374 = t275 + t363;
t235 = t245 * qJD(1);
t344 = qJD(4) * t290;
t371 = (qJD(4) * t224 + t213) * t293 + t290 * t211 - t227 * t344;
t370 = t237 ^ 2;
t366 = g(3) * t291;
t298 = qJD(1) ^ 2;
t364 = qJ(2) * t298;
t316 = -qJD(4) * t335 - t379 * t290;
t325 = t283 * t294;
t206 = (qJD(1) * t325 + t340) * t293 + t316;
t362 = qJ(5) * t206;
t361 = qJ(5) * t235;
t358 = t237 * t283;
t222 = t293 * t227;
t200 = pkin(4) * t283 + t201;
t355 = t200 - t201;
t354 = t293 * t228 - t221;
t246 = -t290 * t291 + t293 * t294;
t353 = -t283 * t215 + t282 * t246;
t251 = pkin(3) * t348 + qJD(1) * qJ(2);
t288 = t294 ^ 2;
t350 = t291 ^ 2 - t288;
t297 = qJD(3) ^ 2;
t349 = -t297 - t298;
t331 = -pkin(4) * t235 - qJD(5);
t217 = -t331 + t251;
t342 = qJD(5) + t217;
t260 = pkin(3) * t345 + qJD(2);
t338 = qJDD(3) * t291;
t336 = pkin(3) * t347;
t266 = pkin(3) * t291 + qJ(2);
t332 = pkin(4) * t277 + t266;
t328 = -t228 * t290 - t222;
t327 = t249 * t290 - t293 * t250;
t323 = MDP(22) * t283;
t225 = t380 * pkin(3) + t284 + t285;
t319 = -t205 * t246 - t215 * t237;
t216 = -t290 * t346 - t291 * t344 + t293 * t325;
t318 = -t216 * t283 - t245 * t282;
t317 = -t224 * t290 - t222;
t315 = t322 * t277;
t314 = t322 * t278;
t313 = 0.2e1 * qJ(2) * t341 + qJDD(3) * t296;
t243 = t365 * t346;
t244 = qJD(3) * t250;
t312 = t290 * t243 - t293 * t244 + t249 * t344 - t250 * t343;
t311 = -t373 - t364;
t195 = pkin(4) * t206 + qJDD(5) + t225;
t234 = t235 ^ 2;
t309 = t235 * t237 * MDP(14) - t320 * MDP(16) + (t358 + (-t283 * t347 - t340) * t293 - t316) * MDP(17) + (-t234 + t370) * MDP(15) + t282 * MDP(18);
t307 = t317 * qJD(4) + t376;
t306 = -t352 * qJD(4) + t293 * t243 + t244 * t290;
t191 = -qJD(5) * t237 + t307 + t374;
t192 = -qJD(5) * t235 - t362 + t371;
t202 = -t317 - t361;
t305 = t191 * t246 + t192 * t245 - t200 * t215 + t202 * t216 - t373;
t304 = -t296 * t297 + t378;
t303 = g(3) * t278 + t373 * t277 - t371;
t302 = t307 + t383;
t301 = t251 * t235 + t303;
t300 = -t251 * t237 + t302;
t299 = t342 * t235 + t303 + t362;
t279 = qJ(5) + t365;
t276 = qJDD(3) * t294;
t270 = pkin(3) * t293 + pkin(4);
t226 = pkin(4) * t245 + t266;
t219 = pkin(4) * t237 + t336;
t212 = pkin(4) * t216 + t260;
t210 = -qJ(5) * t245 + t352;
t209 = -qJ(5) * t246 + t327;
t204 = -t229 + t354;
t203 = t328 + t361;
t194 = qJ(5) * t215 - qJD(5) * t246 + t306;
t193 = -qJ(5) * t216 - qJD(5) * t245 + t312;
t1 = [qJDD(1) * MDP(1) + t373 * MDP(2) + t322 * MDP(3) + (-0.2e1 * t360 + t381) * MDP(4) + t378 * MDP(5) + (t377 * pkin(1) + (t310 + t284) * qJ(2)) * MDP(6) + (qJDD(1) * t288 - 0.2e1 * t291 * t333) * MDP(7) + 0.2e1 * (-t291 * t339 + t350 * t341) * MDP(8) + (-t291 * t297 + t276) * MDP(9) + (-t294 * t297 - t338) * MDP(10) + (t304 * t291 + t313 * t294) * MDP(12) + (-t313 * t291 + t304 * t294) * MDP(13) + t319 * MDP(14) + (t205 * t245 - t206 * t246 + t215 * t235 - t216 * t237) * MDP(15) + t353 * MDP(16) + t318 * MDP(17) + (t266 * t206 + t251 * t216 + t225 * t245 + t260 * t235 + t327 * t282 + t306 * t283 - t315) * MDP(19) + (-t266 * t205 - t251 * t215 + t225 * t246 + t260 * t237 - t352 * t282 - t312 * t283 - t314) * MDP(20) + (t194 * t283 + t195 * t245 + t206 * t226 + t209 * t282 + t212 * t235 + t216 * t217 - t315) * MDP(21) + (-t193 * t283 + t195 * t246 - t205 * t226 - t210 * t282 + t212 * t237 - t215 * t217 - t314) * MDP(22) + (-t193 * t235 - t194 * t237 + t205 * t209 - t206 * t210 - t305) * MDP(23) + (t192 * t210 + t202 * t193 + t191 * t209 + t200 * t194 + t195 * t226 + t217 * t212 - g(1) * (-t279 * t292 + t332 * t295) - g(2) * (t279 * t295 + t332 * t292)) * MDP(24); qJDD(1) * MDP(4) - t298 * MDP(5) + (-t364 - t377) * MDP(6) + (t349 * t291 + t276) * MDP(12) + (t349 * t294 - t338) * MDP(13) + (-t206 * t245 - t216 * t235 - t319) * MDP(23) + (-qJD(1) * t217 + t305) * MDP(24) + (MDP(19) + MDP(21)) * (-qJD(1) * t235 + t353) + (MDP(20) + MDP(22)) * (-qJD(1) * t237 + t318); MDP(9) * t339 - MDP(10) * t340 + qJDD(3) * MDP(11) + (t311 * t294 + t247 + t366) * MDP(12) + (g(3) * t294 + (-t257 - t311) * t291) * MDP(13) + (-t328 * t283 + (-t235 * t347 + t282 * t293 - t283 * t344) * pkin(3) + t300) * MDP(19) + (-t237 * t336 + t354 * t283 + t301 + t375) * MDP(20) + (-t203 * t283 - t219 * t235 + t270 * t282 - t342 * t237 + (-t222 + (-t224 - t367) * t290) * qJD(4) + t383 + t374 + t376) * MDP(21) + (t204 * t283 - t219 * t237 + t299 + t375) * MDP(22) + (t205 * t270 + (t202 + t203) * t237 + (-t200 + t204) * t235 + (-t206 * t290 + (-t235 * t293 + t237 * t290) * qJD(4)) * pkin(3)) * MDP(23) + (t191 * t270 - t200 * t203 - t202 * t204 - t217 * t219 + t383 * pkin(4) + (t366 + t192 * t290 - t373 * t294 + (-t200 * t290 + t202 * t293) * qJD(4)) * pkin(3)) * MDP(24) + t309 + (t294 * t291 * MDP(7) - t350 * MDP(8)) * t298; (-t317 * t283 + t300) * MDP(19) + (t329 * t283 + t301) * MDP(20) + (t363 + t202 * t283 + 0.2e1 * t275 + (-t217 + t331) * t237 + t302) * MDP(21) + (-pkin(4) * t370 + t201 * t283 + t299) * MDP(22) + (pkin(4) * t205 - t355 * t235) * MDP(23) + (t355 * t202 + (-t217 * t237 + t191 + t383) * pkin(4)) * MDP(24) + t309; (t293 * t340 + t316 + t358) * MDP(21) + (-t235 * t283 - t320) * MDP(22) + (-t234 - t370) * MDP(23) + (t200 * t237 + t202 * t235 + t195 - t322) * MDP(24) + (-t323 * t356 + (MDP(21) * t325 - t291 * t323) * t293) * qJD(1);];
tau = t1;
