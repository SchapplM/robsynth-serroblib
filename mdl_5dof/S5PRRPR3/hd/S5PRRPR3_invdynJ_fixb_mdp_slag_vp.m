% Calculate vector of inverse dynamics joint torques for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:47
% EndTime: 2019-12-05 16:19:51
% DurationCPUTime: 1.50s
% Computational Cost: add. (1295->230), mult. (2857->311), div. (0->0), fcn. (2079->10), ass. (0->116)
t267 = sin(pkin(9));
t268 = cos(pkin(9));
t271 = sin(qJ(3));
t273 = cos(qJ(3));
t235 = -t267 * t271 + t268 * t273;
t230 = t235 * qJD(2);
t272 = cos(qJ(5));
t222 = t272 * t230;
t236 = t267 * t273 + t268 * t271;
t232 = t236 * qJD(2);
t270 = sin(qJ(5));
t308 = t232 * t270;
t196 = t222 - t308;
t264 = qJD(3) + qJD(5);
t309 = t196 * t264;
t263 = pkin(8) + qJ(2);
t257 = sin(t263);
t258 = cos(t263);
t292 = g(1) * t258 + g(2) * t257;
t259 = t273 * qJDD(1);
t312 = qJ(4) + pkin(6);
t295 = t312 * qJD(3);
t293 = qJD(2) * t295;
t319 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t312 * qJDD(2);
t183 = qJDD(3) * pkin(3) - t319 * t271 - t273 * t293 + t259;
t186 = (qJDD(1) - t293) * t271 + t319 * t273;
t168 = t268 * t183 - t186 * t267;
t301 = qJD(2) * qJD(3);
t296 = t273 * t301;
t297 = t271 * t301;
t202 = qJDD(2) * t236 - t267 * t297 + t268 * t296;
t163 = qJDD(3) * pkin(4) - pkin(7) * t202 + t168;
t169 = t267 * t183 + t268 * t186;
t231 = t236 * qJD(3);
t201 = -qJD(2) * t231 + qJDD(2) * t235;
t164 = pkin(7) * t201 + t169;
t247 = t312 * t271;
t225 = t273 * qJD(1) - qJD(2) * t247;
t311 = qJD(3) * pkin(3);
t221 = t225 + t311;
t248 = t312 * t273;
t227 = qJD(1) * t271 + qJD(2) * t248;
t307 = t268 * t227;
t185 = t267 * t221 + t307;
t317 = pkin(7) * t230;
t174 = t185 + t317;
t256 = pkin(3) * t273 + pkin(2);
t244 = -t256 * qJD(2) + qJD(4);
t205 = -pkin(4) * t230 + t244;
t260 = qJ(3) + pkin(9) + qJ(5);
t252 = sin(t260);
t253 = cos(t260);
t303 = qJD(5) * t270;
t323 = g(3) * t252 - t270 * t163 - t272 * t164 + t174 * t303 - t205 * t196 + t253 * t292;
t262 = qJDD(3) + qJDD(5);
t286 = t230 * t270 + t272 * t232;
t322 = t262 * MDP(18) + (-t196 ^ 2 + t286 ^ 2) * MDP(15) - t196 * t286 * MDP(14);
t310 = t286 * t264;
t320 = -g(3) * t253 + t272 * t163 - t270 * t164 - t205 * t286 + t252 * t292;
t294 = -t272 * t201 + t202 * t270;
t167 = qJD(5) * t286 + t294;
t318 = pkin(3) * t267;
t316 = pkin(7) * t232;
t313 = g(3) * t273;
t214 = t267 * t227;
t184 = t268 * t221 - t214;
t173 = qJD(3) * pkin(4) + t184 - t316;
t306 = t272 * t173;
t305 = qJDD(1) - g(3);
t189 = t268 * t225 - t214;
t226 = qJD(4) * t273 - t271 * t295;
t228 = -qJD(4) * t271 - t273 * t295;
t190 = t268 * t226 + t267 * t228;
t207 = -t267 * t247 + t268 * t248;
t265 = t271 ^ 2;
t304 = -t273 ^ 2 + t265;
t300 = qJDD(2) * t273;
t299 = t271 * t311;
t298 = qJD(5) * t222 + t270 * t201 + t272 * t202;
t187 = -t225 * t267 - t307;
t188 = -t226 * t267 + t268 * t228;
t206 = -t268 * t247 - t248 * t267;
t291 = g(1) * t257 - g(2) * t258;
t234 = t235 * qJD(3);
t285 = t272 * t235 - t236 * t270;
t170 = qJD(5) * t285 - t231 * t270 + t234 * t272;
t204 = t235 * t270 + t236 * t272;
t290 = t170 * t264 + t204 * t262;
t289 = -t270 * t173 - t272 * t174;
t191 = -pkin(7) * t236 + t206;
t192 = pkin(7) * t235 + t207;
t288 = t191 * t272 - t192 * t270;
t287 = t191 * t270 + t192 * t272;
t254 = pkin(3) * t268 + pkin(4);
t283 = t254 * t270 + t272 * t318;
t282 = t254 * t272 - t270 * t318;
t281 = -0.2e1 * pkin(2) * t301 - pkin(6) * qJDD(3);
t166 = -t232 * t303 + t298;
t280 = pkin(3) * t297 - t256 * qJDD(2) + qJDD(4);
t274 = qJD(3) ^ 2;
t277 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t274 + t291;
t275 = qJD(2) ^ 2;
t276 = pkin(2) * t275 - pkin(6) * qJDD(2) + t292;
t246 = qJDD(3) * t273 - t271 * t274;
t245 = qJDD(3) * t271 + t273 * t274;
t213 = -pkin(4) * t235 - t256;
t209 = pkin(4) * t231 + t299;
t208 = pkin(3) * qJD(2) * t271 + pkin(4) * t232;
t179 = -pkin(4) * t201 + t280;
t178 = -pkin(7) * t231 + t190;
t177 = t189 - t316;
t176 = -pkin(7) * t234 + t188;
t175 = t187 - t317;
t171 = qJD(5) * t204 + t272 * t231 + t234 * t270;
t165 = -t171 * t264 + t262 * t285;
t1 = [t305 * MDP(1) + t246 * MDP(10) - t245 * MDP(11) + (t201 * t236 - t202 * t235 + t230 * t234 + t231 * t232) * MDP(12) + (t168 * t235 + t169 * t236 - t184 * t231 + t185 * t234 - g(3)) * MDP(13) + t165 * MDP(19) - t290 * MDP(20); qJDD(2) * MDP(2) + t291 * MDP(3) + t292 * MDP(4) + (qJDD(2) * t265 + 0.2e1 * t271 * t296) * MDP(5) + 0.2e1 * (t271 * t300 - t304 * t301) * MDP(6) + t245 * MDP(7) + t246 * MDP(8) + (t271 * t281 + t273 * t277) * MDP(10) + (-t271 * t277 + t273 * t281) * MDP(11) + (-t168 * t236 + t169 * t235 - t184 * t234 - t185 * t231 - t188 * t232 + t190 * t230 + t201 * t207 - t202 * t206 - t292) * MDP(12) + (t169 * t207 + t185 * t190 + t168 * t206 + t184 * t188 - t280 * t256 + t244 * t299 - g(1) * (-t256 * t257 + t258 * t312) - g(2) * (t256 * t258 + t257 * t312)) * MDP(13) + (t166 * t204 + t170 * t286) * MDP(14) + (t166 * t285 - t167 * t204 + t170 * t196 - t171 * t286) * MDP(15) + t290 * MDP(16) + t165 * MDP(17) + (-t209 * t196 + t213 * t167 - t179 * t285 + t205 * t171 + (-qJD(5) * t287 + t176 * t272 - t178 * t270) * t264 + t288 * t262 + t291 * t253) * MDP(19) + (t209 * t286 + t213 * t166 + t179 * t204 + t205 * t170 - (qJD(5) * t288 + t176 * t270 + t178 * t272) * t264 - t287 * t262 - t291 * t252) * MDP(20); t271 * qJDD(2) * MDP(7) + MDP(8) * t300 + qJDD(3) * MDP(9) + (t271 * t276 + t259 - t313) * MDP(10) + (-t305 * t271 + t276 * t273) * MDP(11) + ((t185 + t187) * t232 + (t184 - t189) * t230 + (t201 * t267 - t202 * t268) * pkin(3)) * MDP(12) + (-t184 * t187 - t185 * t189 + (-t313 + t168 * t268 + t169 * t267 + (-qJD(2) * t244 + t292) * t271) * pkin(3)) * MDP(13) + (t166 - t309) * MDP(16) + (-t167 + t310) * MDP(17) + (t282 * t262 + t208 * t196 - (t175 * t272 - t177 * t270) * t264 + (-t264 * t283 + t289) * qJD(5) + t320) * MDP(19) + (-t283 * t262 - t208 * t286 + (t175 * t270 + t177 * t272) * t264 + (-t264 * t282 - t306) * qJD(5) + t323) * MDP(20) + (-t271 * t273 * MDP(5) + t304 * MDP(6)) * t275 + t322; (-t230 ^ 2 - t232 ^ 2) * MDP(12) + (t184 * t232 - t185 * t230 + t280 - t291) * MDP(13) + (t167 + t310) * MDP(19) + (t166 + t309) * MDP(20); (t298 - t309) * MDP(16) + (-t294 + t310) * MDP(17) + (-t264 * t289 + t320) * MDP(19) + ((-t174 * t270 + t306) * t264 + t323) * MDP(20) + (-MDP(16) * t308 - MDP(17) * t286 + MDP(19) * t289 - MDP(20) * t306) * qJD(5) + t322;];
tau = t1;
