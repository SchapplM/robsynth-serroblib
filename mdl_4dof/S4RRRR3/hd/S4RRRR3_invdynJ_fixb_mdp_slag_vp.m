% Calculate vector of inverse dynamics joint torques for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:42
% EndTime: 2019-12-31 17:24:44
% DurationCPUTime: 1.91s
% Computational Cost: add. (1509->252), mult. (3506->343), div. (0->0), fcn. (2521->12), ass. (0->128)
t308 = cos(qJ(3));
t309 = cos(qJ(2));
t359 = qJD(1) * t309;
t345 = t308 * t359;
t304 = sin(qJ(3));
t305 = sin(qJ(2));
t360 = qJD(1) * t305;
t346 = t304 * t360;
t244 = -t345 + t346;
t307 = cos(qJ(4));
t246 = -t304 * t359 - t308 * t360;
t303 = sin(qJ(4));
t369 = t246 * t303;
t222 = -t307 * t244 + t369;
t350 = -qJD(3) - qJD(4);
t294 = qJD(2) - t350;
t380 = t222 * t294;
t328 = t244 * t303 + t307 * t246;
t379 = t294 * t328;
t353 = qJD(1) * qJD(2);
t343 = t309 * t353;
t352 = qJDD(1) * t305;
t378 = t343 + t352;
t306 = sin(qJ(1));
t310 = cos(qJ(1));
t334 = g(1) * t310 + g(2) * t306;
t375 = pkin(5) + pkin(6);
t270 = t375 * t309;
t262 = qJD(1) * t270;
t251 = t308 * t262;
t269 = t375 * t305;
t260 = qJD(1) * t269;
t370 = qJD(2) * pkin(2);
t253 = -t260 + t370;
t327 = -t304 * t253 - t251;
t373 = pkin(7) * t244;
t211 = -t327 - t373;
t292 = -pkin(2) * t309 - pkin(1);
t268 = t292 * qJD(1);
t233 = pkin(3) * t244 + t268;
t302 = qJ(2) + qJ(3);
t297 = qJ(4) + t302;
t289 = sin(t297);
t290 = cos(t297);
t356 = qJD(4) * t303;
t325 = g(3) * t289 + t211 * t356 - t233 * t222 + t290 * t334;
t299 = qJD(2) + qJD(3);
t351 = qJDD(1) * t309;
t212 = qJD(3) * t345 - t299 * t346 + t304 * t351 + t308 * t378;
t298 = qJDD(2) + qJDD(3);
t231 = qJDD(2) * pkin(2) - t375 * t378;
t344 = t305 * t353;
t232 = t375 * (-t344 + t351);
t317 = t327 * qJD(3) + t308 * t231 - t304 * t232;
t193 = pkin(3) * t298 - pkin(7) * t212 + t317;
t256 = t304 * t309 + t305 * t308;
t230 = t299 * t256;
t332 = t304 * t352 - t308 * t351;
t213 = qJD(1) * t230 + t332;
t358 = qJD(3) * t304;
t376 = t308 * (qJD(3) * t253 + t232) + t304 * t231 - t262 * t358;
t195 = -pkin(7) * t213 + t376;
t318 = -g(3) * t290 + t307 * t193 - t303 * t195 + t233 * t328 + t289 * t334;
t293 = qJDD(4) + t298;
t377 = t293 * MDP(22) + t222 * t328 * MDP(18) + (-t222 ^ 2 + t328 ^ 2) * MDP(19);
t240 = t246 * pkin(7);
t247 = t304 * t262;
t340 = t308 * t253 - t247;
t210 = t240 + t340;
t362 = -t304 * t269 + t308 * t270;
t342 = t212 * t303 + t307 * t213;
t197 = -qJD(4) * t328 + t342;
t374 = pkin(2) * t304;
t368 = t303 * t293;
t367 = t304 * t307;
t205 = pkin(3) * t299 + t210;
t366 = t307 * t205;
t365 = t307 * t211;
t364 = t307 * t293;
t363 = -t308 * t260 - t247;
t300 = t305 ^ 2;
t361 = -t309 ^ 2 + t300;
t357 = qJD(3) * t308;
t355 = qJD(4) * t307;
t349 = t305 * t370;
t348 = t307 * t212 - t303 * t213 - t244 * t355;
t347 = qJD(2) * t375;
t339 = t260 * t304 - t251;
t338 = -t308 * t269 - t270 * t304;
t336 = -qJD(4) * t205 - t195;
t333 = g(1) * t306 - g(2) * t310;
t331 = -t303 * t205 - t365;
t217 = -pkin(7) * t256 + t338;
t255 = t304 * t305 - t308 * t309;
t218 = -pkin(7) * t255 + t362;
t330 = t217 * t307 - t218 * t303;
t329 = t217 * t303 + t218 * t307;
t227 = t307 * t255 + t256 * t303;
t228 = -t255 * t303 + t256 * t307;
t326 = -0.2e1 * pkin(1) * t353 - pkin(5) * qJDD(2);
t241 = pkin(2) * t344 + qJDD(1) * t292;
t261 = t305 * t347;
t263 = t309 * t347;
t322 = -t308 * t261 - t304 * t263 - t269 * t357 - t270 * t358;
t196 = t246 * t356 + t348;
t311 = qJD(2) ^ 2;
t320 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t311 + t333;
t312 = qJD(1) ^ 2;
t319 = pkin(1) * t312 - pkin(5) * qJDD(1) + t334;
t316 = -qJD(3) * t362 + t261 * t304 - t308 * t263;
t295 = sin(t302);
t296 = cos(t302);
t315 = g(3) * t295 + t268 * t244 + t296 * t334 - t376;
t314 = -t246 * t244 * MDP(11) + (t196 - t380) * MDP(20) + (-t197 - t379) * MDP(21) + (t244 * t299 + t212) * MDP(13) + (-t332 + (-qJD(1) * t256 - t246) * t299) * MDP(14) + (-t244 ^ 2 + t246 ^ 2) * MDP(12) + t298 * MDP(15) + t377;
t313 = -g(3) * t296 + t268 * t246 + t295 * t334 + t317;
t291 = pkin(2) * t308 + pkin(3);
t236 = pkin(3) * t255 + t292;
t234 = pkin(2) * t360 - pkin(3) * t246;
t229 = t299 * t255;
t224 = pkin(3) * t230 + t349;
t216 = t240 + t363;
t215 = t339 + t373;
t204 = pkin(3) * t213 + t241;
t201 = pkin(7) * t229 + t316;
t200 = -pkin(7) * t230 + t322;
t199 = qJD(4) * t228 - t229 * t303 + t307 * t230;
t198 = -qJD(4) * t227 - t229 * t307 - t230 * t303;
t1 = [qJDD(1) * MDP(1) + t333 * MDP(2) + t334 * MDP(3) + (qJDD(1) * t300 + 0.2e1 * t305 * t343) * MDP(4) + 0.2e1 * (t305 * t351 - t353 * t361) * MDP(5) + (qJDD(2) * t305 + t309 * t311) * MDP(6) + (qJDD(2) * t309 - t305 * t311) * MDP(7) + (t305 * t326 + t309 * t320) * MDP(9) + (-t305 * t320 + t309 * t326) * MDP(10) + (t212 * t256 + t229 * t246) * MDP(11) + (-t212 * t255 - t213 * t256 + t229 * t244 + t230 * t246) * MDP(12) + (-t229 * t299 + t256 * t298) * MDP(13) + (-t230 * t299 - t255 * t298) * MDP(14) + (t292 * t213 + t268 * t230 + t241 * t255 + t244 * t349 + t296 * t333 + t298 * t338 + t299 * t316) * MDP(16) + (t292 * t212 - t268 * t229 + t241 * t256 - t246 * t349 - t295 * t333 - t298 * t362 - t299 * t322) * MDP(17) + (t196 * t228 - t198 * t328) * MDP(18) + (-t196 * t227 - t197 * t228 + t198 * t222 + t199 * t328) * MDP(19) + (t198 * t294 + t228 * t293) * MDP(20) + (-t199 * t294 - t227 * t293) * MDP(21) + (-t224 * t222 + t236 * t197 + t204 * t227 + t233 * t199 + (-qJD(4) * t329 - t200 * t303 + t201 * t307) * t294 + t330 * t293 + t333 * t290) * MDP(23) + (-t224 * t328 + t236 * t196 + t204 * t228 + t233 * t198 - (qJD(4) * t330 + t200 * t307 + t201 * t303) * t294 - t329 * t293 - t333 * t289) * MDP(24); (g(3) * t305 + t309 * t319) * MDP(10) + MDP(6) * t352 + (-g(3) * t309 + t305 * t319) * MDP(9) + (t363 * t299 + (t246 * t360 - t298 * t304 - t299 * t357) * pkin(2) + t315) * MDP(17) + (-t339 * t299 + (-t244 * t360 + t298 * t308 - t299 * t358) * pkin(2) + t313) * MDP(16) + qJDD(2) * MDP(8) + (t291 * t364 + t234 * t222 - (t215 * t307 - t216 * t303) * t294 + (-t304 * t368 + (-t303 * t308 - t367) * t294 * qJD(3)) * pkin(2) + ((-pkin(2) * t367 - t303 * t291) * t294 + t331) * qJD(4) + t318) * MDP(23) + t314 + (t234 * t328 + (-t291 * t293 - t193 + (-t350 * t374 + t215) * t294) * t303 + (-t293 * t374 + (-pkin(2) * t357 - qJD(4) * t291 + t216) * t294 + t336) * t307 + t325) * MDP(24) + MDP(7) * t351 + (-t305 * t309 * MDP(4) + MDP(5) * t361) * t312; (-(-t210 * t303 - t365) * t294 + t331 * qJD(4) + (-t222 * t246 - t294 * t356 + t364) * pkin(3) + t318) * MDP(23) + (-t299 * t327 + t313) * MDP(16) + t314 + (t299 * t340 + t315) * MDP(17) + ((-t211 * t294 - t193) * t303 + (t210 * t294 + t336) * t307 + (-t246 * t328 - t294 * t355 - t368) * pkin(3) + t325) * MDP(24); (t348 - t380) * MDP(20) + (-t342 - t379) * MDP(21) + (-t294 * t331 + t318) * MDP(23) + (-t307 * t195 - t303 * t193 + (-t211 * t303 + t366) * t294 + t325) * MDP(24) + (MDP(20) * t369 + MDP(21) * t328 + t331 * MDP(23) - MDP(24) * t366) * qJD(4) + t377;];
tau = t1;
