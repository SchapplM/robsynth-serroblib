% Calculate vector of inverse dynamics joint torques for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:40:05
% EndTime: 2019-12-05 17:40:09
% DurationCPUTime: 2.10s
% Computational Cost: add. (1259->248), mult. (2497->310), div. (0->0), fcn. (1848->12), ass. (0->122)
t272 = sin(pkin(8));
t273 = cos(pkin(8));
t276 = sin(qJ(4));
t279 = cos(qJ(4));
t231 = t272 * t279 + t273 * t276;
t221 = t231 * qJD(1);
t278 = cos(qJ(5));
t318 = qJD(1) * t276;
t307 = t272 * t318;
t317 = qJD(1) * t279;
t309 = t273 * t317;
t223 = -t307 + t309;
t275 = sin(qJ(5));
t330 = t223 * t275;
t192 = t278 * t221 + t330;
t268 = qJD(4) + qJD(5);
t333 = t192 * t268;
t277 = sin(qJ(1));
t280 = cos(qJ(1));
t347 = g(1) * t277 - g(2) * t280;
t325 = t273 * t279;
t232 = -t272 * t276 + t325;
t274 = -pkin(1) - qJ(3);
t241 = t274 * qJD(1) + qJD(2);
t305 = -pkin(6) * qJD(1) + t241;
t215 = t305 * t272;
t216 = t305 * t273;
t290 = -t215 * t279 - t216 * t276;
t187 = -pkin(7) * t221 - t290;
t254 = qJD(1) * qJ(2) + qJD(3);
t258 = t272 * pkin(3);
t236 = qJD(1) * t258 + t254;
t203 = pkin(4) * t221 + t236;
t267 = pkin(8) + qJ(4);
t257 = qJ(5) + t267;
t249 = sin(t257);
t250 = cos(t257);
t314 = qJD(5) * t275;
t353 = g(3) * t250 + t187 * t314 + t203 * t192 + t249 * t347;
t264 = qJDD(4) + qJDD(5);
t289 = -t221 * t275 + t278 * t223;
t352 = t264 * MDP(22) + t192 * MDP(18) * t289 + (-t192 ^ 2 + t289 ^ 2) * MDP(19);
t320 = t272 ^ 2 + t273 ^ 2;
t350 = t241 * t320;
t334 = t289 * t268;
t349 = -t215 * t276 + t279 * t216;
t202 = -t231 * t275 + t278 * t232;
t339 = -pkin(6) + t274;
t234 = t339 * t272;
t235 = t339 * t273;
t323 = t279 * t234 + t276 * t235;
t269 = qJDD(1) * qJ(2);
t270 = qJD(1) * qJD(2);
t346 = t269 + t270;
t239 = qJDD(3) + t346;
t298 = g(1) * t280 + g(2) * t277;
t348 = t239 - t298;
t345 = t272 * MDP(7) + t273 * MDP(8);
t201 = t278 * t231 + t232 * t275;
t224 = t231 * qJD(4);
t315 = qJD(4) * t279;
t308 = t273 * t315;
t316 = qJD(4) * t276;
t225 = -t272 * t316 + t308;
t180 = -t201 * qJD(5) - t278 * t224 - t225 * t275;
t344 = t180 * t268 + t202 * t264;
t242 = qJDD(1) * t325;
t311 = qJDD(1) * t272;
t296 = -t276 * t311 + t242;
t199 = -qJD(1) * t224 + t296;
t341 = -qJD(1) * qJD(3) + qJDD(1) * t274;
t233 = qJDD(2) + t341;
t302 = -pkin(6) * qJDD(1) + t233;
t205 = t302 * t272;
t206 = t302 * t273;
t300 = -t276 * t205 + t279 * t206;
t176 = qJDD(4) * pkin(4) - pkin(7) * t199 + t290 * qJD(4) + t300;
t240 = qJD(4) * t307;
t342 = -t231 * qJDD(1) + t240;
t200 = qJD(1) * t308 - t342;
t291 = t279 * t205 + t276 * t206;
t177 = -pkin(7) * t200 + t349 * qJD(4) + t291;
t343 = g(3) * t249 + t278 * t176 - t275 * t177 - t203 * t289 - t250 * t347;
t301 = t199 * t275 + t278 * t200;
t179 = t289 * qJD(5) + t301;
t340 = 0.2e1 * t270;
t338 = pkin(1) * qJDD(1);
t186 = -pkin(7) * t223 + t349;
t185 = qJD(4) * pkin(4) + t186;
t336 = t185 * t278;
t335 = t187 * t278;
t248 = qJ(2) + t258;
t324 = -t224 * qJD(4) + t232 * qJDD(4);
t322 = t280 * pkin(1) + t277 * qJ(2);
t313 = qJD(5) * t278;
t310 = t278 * t199 - t275 * t200 - t221 * t313;
t306 = qJDD(2) - t347;
t304 = t320 * MDP(9);
t303 = t320 * t233;
t299 = -t234 * t276 + t279 * t235;
t227 = pkin(3) * t311 + t239;
t181 = t202 * qJD(5) - t224 * t275 + t278 * t225;
t295 = -t181 * t268 - t201 * t264;
t294 = -t185 * t275 - t335;
t189 = -pkin(7) * t232 + t299;
t190 = -pkin(7) * t231 + t323;
t293 = t189 * t278 - t190 * t275;
t292 = t189 * t275 + t190 * t278;
t288 = -qJD(4) * t225 - qJDD(4) * t231;
t178 = -t223 * t314 + t310;
t283 = -t231 * qJD(3) - t234 * t316 + t235 * t315;
t282 = -qJD(3) * t232 - t323 * qJD(4);
t281 = qJD(1) ^ 2;
t260 = t280 * qJ(2);
t256 = cos(t267);
t255 = sin(t267);
t210 = pkin(4) * t225 + qJD(2);
t208 = pkin(4) * t231 + t248;
t188 = pkin(4) * t200 + t227;
t184 = pkin(7) * t224 + t282;
t183 = -pkin(7) * t225 + t283;
t1 = [qJDD(1) * MDP(1) + t347 * MDP(2) + t298 * MDP(3) + (t306 - 0.2e1 * t338) * MDP(4) + (0.2e1 * t269 + t340 - t298) * MDP(5) + (-(qJDD(2) - t338) * pkin(1) - g(1) * (-pkin(1) * t277 + t260) - g(2) * t322 + (t269 + t340) * qJ(2)) * MDP(6) + (t347 + t320 * (-t233 - t341)) * MDP(9) + (t239 * qJ(2) + t254 * qJD(2) - g(1) * (t274 * t277 + t260) - g(2) * (qJ(3) * t280 + t322) + t274 * t303 - qJD(3) * t350) * MDP(10) + (t199 * t232 - t223 * t224) * MDP(11) + (-t199 * t231 - t200 * t232 + t221 * t224 - t223 * t225) * MDP(12) + t324 * MDP(13) + t288 * MDP(14) + (qJD(2) * t221 + t282 * qJD(4) + t299 * qJDD(4) + t248 * t200 + t236 * t225 + t227 * t231 - t298 * t255) * MDP(16) + (qJD(2) * t223 - t283 * qJD(4) - t323 * qJDD(4) + t248 * t199 - t236 * t224 + t227 * t232 - t298 * t256) * MDP(17) + (t178 * t202 + t180 * t289) * MDP(18) + (-t178 * t201 - t179 * t202 - t180 * t192 - t181 * t289) * MDP(19) + t344 * MDP(20) + t295 * MDP(21) + (t210 * t192 + t208 * t179 + t188 * t201 + t203 * t181 + (-t292 * qJD(5) - t183 * t275 + t184 * t278) * t268 + t293 * t264 - t298 * t249) * MDP(23) + (t210 * t289 + t208 * t178 + t188 * t202 + t203 * t180 - (t293 * qJD(5) + t183 * t278 + t184 * t275) * t268 - t292 * t264 - t298 * t250) * MDP(24) + t345 * (t346 + t348); t306 * MDP(6) + (-qJD(1) * t254 + t303 - t347) * MDP(10) + (-qJD(1) * t221 + t324) * MDP(16) + (-qJD(1) * t223 + t288) * MDP(17) + (-qJD(1) * t192 + t344) * MDP(23) + (-qJD(1) * t289 + t295) * MDP(24) + (-MDP(6) * qJ(2) - MDP(5) - t345) * t281 + (-pkin(1) * MDP(6) + MDP(4) - t304) * qJDD(1); (qJD(1) * t350 + t348) * MDP(10) - t240 * MDP(16) + t242 * MDP(17) + (t179 + t334) * MDP(23) + (t178 - t333) * MDP(24) - t281 * t304 + ((MDP(16) * t276 + MDP(8)) * t273 + (MDP(16) * t279 - MDP(17) * t276 + MDP(7)) * t272) * qJDD(1) + ((t223 + t309) * MDP(16) + (-t272 * t317 - t273 * t318 - t221) * MDP(17)) * qJD(4); t223 * t221 * MDP(11) + (-t221 ^ 2 + t223 ^ 2) * MDP(12) + t296 * MDP(13) + ((t223 - t309) * qJD(4) + t342) * MDP(14) + qJDD(4) * MDP(15) + (g(3) * t255 - t236 * t223 - t256 * t347 + t300) * MDP(16) + (g(3) * t256 + t236 * t221 + t255 * t347 - t291) * MDP(17) + (t178 + t333) * MDP(20) + (-t179 + t334) * MDP(21) + (-(-t186 * t275 - t335) * t268 + t294 * qJD(5) + (-t192 * t223 + t264 * t278 - t268 * t314) * pkin(4) + t343) * MDP(23) + ((-t187 * t268 - t176) * t275 + (-qJD(5) * t185 + t186 * t268 - t177) * t278 + (-t223 * t289 - t264 * t275 - t268 * t313) * pkin(4) + t353) * MDP(24) + t352; (t310 + t333) * MDP(20) + (-t301 + t334) * MDP(21) + (-t294 * t268 + t343) * MDP(23) + (-t278 * t177 - t275 * t176 + (-t187 * t275 + t336) * t268 + t353) * MDP(24) + (-MDP(20) * t330 - t289 * MDP(21) + t294 * MDP(23) - MDP(24) * t336) * qJD(5) + t352;];
tau = t1;
