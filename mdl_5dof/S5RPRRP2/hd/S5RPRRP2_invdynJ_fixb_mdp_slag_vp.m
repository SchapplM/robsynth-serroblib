% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:15
% EndTime: 2022-01-23 09:28:17
% DurationCPUTime: 1.40s
% Computational Cost: add. (1348->247), mult. (2232->284), div. (0->0), fcn. (1229->12), ass. (0->133)
t333 = qJDD(2) - g(3);
t250 = cos(pkin(8));
t234 = pkin(1) * t250 + pkin(2);
t249 = sin(pkin(8));
t328 = pkin(1) * t249;
t298 = qJD(1) * t328;
t332 = -qJD(3) * t298 + t234 * qJDD(1);
t220 = t234 * qJD(1);
t331 = qJD(3) * t220 + qJDD(1) * t328;
t255 = cos(qJ(4));
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t196 = t220 * t253 + t256 * t298;
t245 = qJD(1) + qJD(3);
t188 = pkin(7) * t245 + t196;
t287 = qJ(5) * t245 + t188;
t274 = t287 * t255;
t246 = qJ(1) + pkin(8);
t238 = qJ(3) + t246;
t232 = sin(t238);
t227 = g(2) * t232;
t233 = cos(t238);
t229 = g(1) * t233;
t306 = -t227 - t229;
t284 = t234 * t256 - t253 * t328;
t305 = t253 * t234 + t256 * t328;
t252 = sin(qJ(4));
t300 = qJD(4) * t252;
t291 = t245 * t300;
t325 = pkin(4) * t255;
t235 = pkin(3) + t325;
t244 = qJDD(1) + qJDD(3);
t314 = t235 * t244;
t330 = -pkin(4) * t291 + t314;
t329 = -t253 * t332 - t256 * t331;
t327 = pkin(3) * t244;
t247 = t252 ^ 2;
t326 = pkin(4) * t247;
t228 = g(1) * t232;
t324 = g(2) * t233;
t323 = g(3) * t255;
t251 = -qJ(5) - pkin(7);
t322 = qJD(4) * pkin(4);
t321 = qJDD(4) * pkin(4);
t195 = t220 * t256 - t253 * t298;
t187 = -pkin(3) * t245 - t195;
t320 = t187 * t245;
t319 = t196 * t245;
t199 = t305 * qJD(3);
t318 = t199 * t245;
t317 = t232 * t255;
t316 = t233 * t252;
t313 = t244 * t252;
t312 = t244 * t255;
t311 = t245 * t252;
t310 = t245 * t255;
t309 = t252 * t255;
t201 = pkin(7) + t305;
t308 = -qJ(5) - t201;
t240 = t255 * qJD(2);
t177 = -t252 * t287 + t240;
t176 = t177 + t322;
t307 = t176 - t177;
t248 = t255 ^ 2;
t304 = -t247 - t248;
t303 = t247 - t248;
t301 = qJD(4) * t245;
t299 = qJD(4) * t255;
t297 = pkin(4) * t300;
t184 = t188 * t300;
t174 = pkin(7) * t244 - t329;
t281 = -qJD(4) * qJD(2) - t174;
t267 = -qJ(5) * t244 + t281;
t262 = qJD(5) * t245 - t267;
t168 = -t184 + (-qJ(5) * t301 + qJDD(2)) * t252 + t262 * t255;
t295 = t168 * t255 + t306;
t279 = -t253 * t331 + t256 * t332;
t170 = qJDD(5) - t279 - t330;
t181 = -t235 * t245 + qJD(5) - t195;
t213 = g(2) * t316;
t294 = t170 * t252 + t181 * t299 + t213;
t175 = -t279 - t327;
t293 = t175 * t252 + t187 * t299 + t213;
t214 = g(1) * t317;
t292 = t195 * t300 + t196 * t310 + t214;
t290 = -t170 - t324;
t289 = -t175 - t324;
t288 = qJD(4) * t251;
t286 = -t232 * t251 + t233 * t235;
t283 = qJD(4) * t308;
t282 = 0.2e1 * t245 * t299;
t200 = -pkin(3) - t284;
t258 = qJD(4) ^ 2;
t278 = -pkin(7) * t258 + t327;
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t277 = g(1) * t254 - g(2) * t257;
t276 = -t319 - t228;
t209 = qJDD(4) * t252 + t255 * t258;
t210 = qJDD(4) * t255 - t252 * t258;
t275 = 0.2e1 * (t244 * t309 - t301 * t303) * MDP(9) + (t244 * t247 + t252 * t282) * MDP(8) + t209 * MDP(10) + t210 * MDP(11) + t244 * MDP(5);
t178 = qJD(2) * t252 + t274;
t273 = t176 * t252 - t178 * t255;
t193 = t199 + t297;
t194 = t200 - t325;
t272 = t193 * t245 + t194 * t244;
t271 = -t232 * t235 - t233 * t251;
t237 = t255 * qJDD(2);
t270 = g(1) * t316 + t227 * t252 + t237 - t323;
t269 = -pkin(3) * t301 - pkin(7) * qJDD(4);
t266 = g(2) * t317 + t255 * t229 - t252 * t333 + t184;
t265 = -t228 - t279 + t324;
t264 = t200 * t244 + t201 * t258 + t318;
t198 = t284 * qJD(3);
t263 = -qJDD(4) * t201 + (t200 * t245 - t198) * qJD(4);
t261 = (-qJD(5) - t181) * t245 + t267;
t260 = -t306 + t329;
t243 = t245 ^ 2;
t241 = t255 * qJ(5);
t239 = t255 * qJD(5);
t222 = pkin(7) * t255 + t241;
t221 = t251 * t252;
t203 = -qJD(5) * t252 + t255 * t288;
t202 = t252 * t288 + t239;
t192 = t201 * t255 + t241;
t191 = t308 * t252;
t190 = t195 * t299;
t182 = t187 * t300;
t179 = t181 * t300;
t172 = (-qJD(5) - t198) * t252 + t255 * t283;
t171 = t198 * t255 + t252 * t283 + t239;
t167 = -qJD(4) * t274 - t252 * t262 + t237 + t321;
t1 = [qJDD(1) * MDP(1) + t277 * MDP(2) + (g(1) * t257 + g(2) * t254) * MDP(3) + (t277 + (t249 ^ 2 + t250 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t244 * t284 - t265 - t318) * MDP(6) + (-t198 * t245 - t244 * t305 + t260) * MDP(7) + (t182 + t214 + t263 * t252 + (-t264 + t289) * t255) * MDP(13) + (t263 * t255 + (t264 - t228) * t252 + t293) * MDP(14) + (qJDD(4) * t191 + t179 + t214 + (t194 * t311 + t172) * qJD(4) + (-t272 + t290) * t255) * MDP(15) + (-qJDD(4) * t192 + (t194 * t310 - t171) * qJD(4) + (t272 - t228) * t252 + t294) * MDP(16) + ((t171 * t245 + t192 * t244 + (-t191 * t245 - t176) * qJD(4)) * t255 + (-t172 * t245 - t191 * t244 - t167 + (-t192 * t245 - t178) * qJD(4)) * t252 + t295) * MDP(17) + (t168 * t192 + t178 * t171 + t167 * t191 + t176 * t172 + t170 * t194 + t181 * t193 - g(1) * (-pkin(2) * sin(t246) - t254 * pkin(1) + t271) - g(2) * (pkin(2) * cos(t246) + t257 * pkin(1) + t286)) * MDP(18) + t275; t333 * MDP(4) + (-qJD(4) * t273 + t167 * t255 + t168 * t252 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t210 + (-MDP(14) - MDP(16)) * t209; (-t265 + t319) * MDP(6) + (t195 * t245 + t260) * MDP(7) + (t182 + t269 * t252 + (t278 + t289) * t255 + t292) * MDP(13) + (t190 + t269 * t255 + (t276 - t278) * t252 + t293) * MDP(14) + (qJDD(4) * t221 + t179 + (-t235 * t311 + t203) * qJD(4) + (t290 + t330) * t255 + t292) * MDP(15) + (-qJDD(4) * t222 + t190 + (t276 - t314) * t252 + (-t202 + (-t235 * t255 + t326) * t245) * qJD(4) + t294) * MDP(16) + ((-qJD(4) * t176 + t222 * t244) * t255 + (-qJD(4) * t178 - t221 * t244 - t167) * t252 + (t202 * t255 - t203 * t252 + t304 * t195 + (-t221 * t255 - t222 * t252) * qJD(4)) * t245 + t295) * MDP(17) + (t168 * t222 + t167 * t221 - t170 * t235 - g(1) * t271 - g(2) * t286 + (-t196 + t297) * t181 + (-t195 * t255 + t202) * t178 + (t195 * t252 + t203) * t176) * MDP(18) + t275; -t243 * MDP(8) * t309 + t303 * t243 * MDP(9) + MDP(10) * t313 + MDP(11) * t312 + qJDD(4) * MDP(12) + ((-t174 - t320) * t252 + t270) * MDP(13) + ((-t188 * t252 + t240) * qJD(4) + (t281 - t320) * t255 + t266) * MDP(14) + (0.2e1 * t321 + (t178 - t274) * qJD(4) + (t243 * t325 + t261) * t252 + t270) * MDP(15) + (-t243 * t326 + (qJ(5) * t311 + t177) * qJD(4) + t261 * t255 + t266) * MDP(16) + (-pkin(4) * t313 + (t307 - t322) * t310) * MDP(17) + (t307 * t178 + (-t323 + t167 + (-t181 * t245 - t306) * t252) * pkin(4)) * MDP(18); (0.2e1 * t291 - t312) * MDP(15) + (t282 + t313) * MDP(16) + (t245 * t273 - t228 - t290) * MDP(18) + t304 * MDP(17) * t243;];
tau = t1;
