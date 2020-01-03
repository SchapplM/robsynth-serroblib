% Calculate vector of inverse dynamics joint torques for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:12
% EndTime: 2019-12-31 19:06:15
% DurationCPUTime: 2.27s
% Computational Cost: add. (1480->237), mult. (2158->311), div. (0->0), fcn. (1360->10), ass. (0->133)
t298 = cos(qJ(4));
t359 = qJD(4) * t298;
t404 = -qJD(5) * t298 - t359;
t387 = sin(qJ(1));
t388 = cos(qJ(1));
t403 = g(1) * t387 - g(2) * t388;
t296 = sin(qJ(3));
t301 = qJD(4) ^ 2;
t353 = qJD(1) - qJD(3);
t341 = t353 ^ 2;
t288 = qJDD(1) - qJDD(3);
t299 = cos(qJ(3));
t372 = t288 * t299;
t402 = t372 + (t301 + t341) * t296;
t294 = sin(qJ(5));
t295 = sin(qJ(4));
t297 = cos(qJ(5));
t245 = t294 * t298 + t295 * t297;
t289 = qJD(4) + qJD(5);
t394 = t289 * t245;
t401 = t353 * t394;
t300 = -pkin(1) - pkin(2);
t362 = qJ(2) * qJD(1);
t400 = -qJD(3) * t362 + qJDD(1) * t300 + qJDD(2);
t399 = -qJDD(2) + t403;
t318 = t404 * t297;
t366 = t297 * t298;
t370 = t294 * t295;
t397 = t366 - t370;
t243 = -t296 * t387 - t299 * t388;
t244 = t296 * t388 - t299 * t387;
t338 = g(1) * t243 + g(2) * t244;
t396 = qJD(4) * t353;
t264 = qJD(1) * t300 + qJD(2);
t354 = (qJD(1) * qJD(2));
t355 = (qJ(2) * qJDD(1));
t395 = qJD(3) * t264 + t354 + t355;
t208 = t397 * t288 - t401;
t342 = -qJ(2) * t296 + t299 * t300;
t364 = qJ(2) * t299 + t296 * t300;
t380 = pkin(1) * qJDD(1);
t393 = t380 + t399;
t323 = t296 * t395 - t299 * t400;
t386 = pkin(3) * t288;
t211 = t323 + t386;
t360 = qJD(4) * t295;
t373 = t288 * t298;
t324 = -t353 * t360 + t373;
t206 = pkin(4) * t324 + t211;
t237 = t264 * t299 - t296 * t362;
t384 = pkin(4) * t298;
t279 = -pkin(3) - t384;
t218 = -t279 * t353 - t237;
t219 = t289 * t370 + t318;
t293 = qJ(4) + qJ(5);
t281 = sin(t293);
t339 = g(1) * t244 - g(2) * t243;
t392 = -t206 * t245 + t218 * t219 + t339 * t281;
t282 = cos(t293);
t391 = t206 * t397 - t218 * t394 - t339 * t282;
t350 = t353 * t370;
t207 = -t245 * t288 + t289 * t350 + t318 * t353;
t227 = t353 * t366 - t350;
t229 = t245 * t353;
t291 = t295 ^ 2;
t348 = t353 * t359;
t363 = -t298 ^ 2 + t291;
t367 = t295 * t298;
t287 = qJDD(4) + qJDD(5);
t375 = t245 * t287;
t376 = t397 * t287;
t390 = (t207 * t245 + t219 * t229) * MDP(17) - (-t207 * t397 + t208 * t245 - t219 * t227 - t229 * t394) * MDP(18) - (t219 * t289 - t375) * MDP(19) - (t289 * t394 - t376) * MDP(20) + (-t288 * t291 - 0.2e1 * t295 * t348) * MDP(10) + 0.2e1 * (-t288 * t367 + t363 * t396) * MDP(11) + (qJDD(4) * t295 + t298 * t301) * MDP(12) + (qJDD(4) * t298 - t295 * t301) * MDP(13) - t288 * MDP(7);
t389 = pkin(7) + pkin(8);
t385 = pkin(3) * t353;
t247 = -pkin(7) + t364;
t381 = pkin(8) - t247;
t238 = t264 * t296 + t299 * t362;
t224 = -pkin(7) * t353 + t238;
t344 = -pkin(8) * t353 + t224;
t217 = t344 * t298;
t379 = t217 * t297;
t236 = qJD(2) * t296 + qJD(3) * t364;
t378 = t236 * t353;
t377 = t238 * t353;
t374 = t288 * t295;
t371 = t353 * t295;
t358 = qJD(5) * t294;
t357 = qJD(5) * t297;
t351 = pkin(4) * t360;
t349 = qJD(4) * t389;
t343 = qJD(4) * t381;
t246 = pkin(3) - t342;
t340 = -t238 + t351;
t216 = t344 * t295;
t215 = qJD(4) * pkin(4) - t216;
t332 = -t215 * t294 - t379;
t225 = t381 * t295;
t226 = t381 * t298;
t329 = t225 * t297 + t226 * t294;
t328 = t225 * t294 - t226 * t297;
t261 = t389 * t295;
t262 = t389 * t298;
t327 = -t261 * t297 - t262 * t294;
t326 = -t261 * t294 + t262 * t297;
t325 = -t211 + t339;
t321 = t289 * t397;
t320 = g(1) * t388 + g(2) * t387;
t317 = -t229 * t227 * MDP(17) + (t227 * t289 + t207) * MDP(19) + (-t229 * t289 - t208) * MDP(20) + (-t227 ^ 2 + t229 ^ 2) * MDP(18) + t287 * MDP(21);
t315 = t296 * t400 + t299 * t395;
t210 = -pkin(7) * t288 + t315;
t223 = -t237 + t385;
t316 = t223 * t353 - t210 - t338;
t314 = -pkin(7) * qJDD(4) + (t223 + t237 + t385) * qJD(4);
t313 = -qJDD(4) * t296 + 0.2e1 * t299 * t396;
t235 = qJD(2) * t299 + qJD(3) * t342;
t312 = -qJDD(4) * t247 + (-t246 * t353 - t223 - t235) * qJD(4);
t311 = -t320 + (2 * t354);
t309 = t323 - t339;
t308 = pkin(7) * t301 - t325 + t377 + t386;
t307 = -t246 * t288 + t247 * t301 + t325 - t378;
t305 = t315 + t338;
t202 = -t224 * t359 + qJDD(4) * pkin(4) - t210 * t295 + (t348 + t374) * pkin(8);
t304 = t217 * t358 + t218 * t227 + (-t217 * t289 - t202) * t294 - g(3) * t281 - t338 * t282;
t203 = -pkin(8) * t324 + t210 * t298 - t224 * t360;
t303 = g(3) * t282 + qJD(5) * t332 + t297 * t202 - t294 * t203 + t218 * t229 - t281 * t338;
t302 = qJD(1) ^ 2;
t249 = t298 * t349;
t248 = t295 * t349;
t239 = t246 + t384;
t222 = t236 - t351;
t213 = -t235 * t295 + t298 * t343;
t212 = t235 * t298 + t295 * t343;
t1 = [(t222 * t227 + t239 * t208 + (-qJD(5) * t328 - t212 * t294 + t213 * t297) * t289 + t329 * t287 + t391) * MDP(22) + (t295 * t307 + t298 * t312) * MDP(16) + qJDD(1) * MDP(1) + (t311 + (2 * t355)) * MDP(5) + (t235 * t353 + t288 * t364 + t305) * MDP(9) + (t393 * pkin(1) + (t311 + t355) * qJ(2)) * MDP(6) + (-t288 * t342 + t309 + t378) * MDP(8) + (0.2e1 * t380 + t399) * MDP(4) + t403 * MDP(2) + (t295 * t312 - t298 * t307) * MDP(15) + (-t222 * t229 + t239 * t207 - (qJD(5) * t329 + t212 * t297 + t213 * t294) * t289 - t328 * t287 + t392) * MDP(23) + t320 * MDP(3) - t390; -qJDD(1) * MDP(4) - t302 * MDP(5) + (-qJ(2) * t302 - t393) * MDP(6) + (-t296 * t341 - t372) * MDP(8) + (t288 * t296 - t299 * t341) * MDP(9) + (t313 * t295 - t298 * t402) * MDP(15) + (t295 * t402 + t313 * t298) * MDP(16) + ((-t208 + t401) * t299 + ((t294 * t360 + t295 * t358 + t318) * t289 - t375 - t353 * t227) * t296) * MDP(22) + ((t321 * t353 - t207) * t299 + (-(t294 * t404 - t295 * t357 - t297 * t360) * t289 - t376 + t353 * t229) * t296) * MDP(23); (-t309 - t377) * MDP(8) + (-t237 * t353 - t305) * MDP(9) + (t295 * t314 - t298 * t308) * MDP(15) + (t295 * t308 + t298 * t314) * MDP(16) + (t279 * t208 + (-qJD(5) * t326 + t248 * t294 - t249 * t297) * t289 + t327 * t287 + t237 * t394 + t340 * t227 - t391) * MDP(22) + (t279 * t207 - (qJD(5) * t327 - t248 * t297 - t249 * t294) * t289 - t326 * t287 + t237 * t321 - t340 * t229 - t392) * MDP(23) + t390; -MDP(12) * t374 - MDP(13) * t373 + qJDD(4) * MDP(14) + (g(3) * t298 + t295 * t316) * MDP(15) + (-g(3) * t295 + t298 * t316) * MDP(16) + (-(t216 * t294 - t379) * t289 + (t227 * t371 + t297 * t287 - t289 * t358) * pkin(4) + t303) * MDP(22) + ((-qJD(5) * t215 - t216 * t289 - t203) * t297 + (-t229 * t371 - t294 * t287 - t289 * t357) * pkin(4) + t304) * MDP(23) + t317 + (-MDP(10) * t367 + MDP(11) * t363) * t341; (-t289 * t332 + t303) * MDP(22) + ((-t203 + (-qJD(5) + t289) * t215) * t297 + t304) * MDP(23) + t317;];
tau = t1;
