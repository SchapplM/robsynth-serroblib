% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:49
% EndTime: 2022-01-23 09:12:51
% DurationCPUTime: 2.28s
% Computational Cost: add. (1142->260), mult. (2242->343), div. (0->0), fcn. (1410->10), ass. (0->128)
t235 = sin(pkin(8));
t242 = cos(qJ(4));
t286 = qJD(1) * qJD(4);
t270 = t242 * t286;
t240 = sin(qJ(4));
t281 = qJDD(1) * t240;
t326 = t270 + t281;
t327 = t235 * t326;
t232 = qJ(1) + pkin(7);
t227 = sin(t232);
t228 = cos(t232);
t268 = -g(1) * t227 + g(2) * t228;
t236 = sin(pkin(7));
t220 = pkin(1) * t236 + qJ(3);
t287 = qJD(1) * qJD(3);
t203 = qJDD(1) * t220 + t287;
t237 = cos(pkin(8));
t325 = -t235 * (-qJ(5) - pkin(6)) + (pkin(4) * t242 + pkin(3)) * t237;
t306 = t237 * t240;
t190 = t227 * t306 + t228 * t242;
t192 = t227 * t242 - t228 * t306;
t324 = -g(1) * t192 + g(2) * t190;
t298 = qJD(1) * t237;
t215 = -qJD(4) + t298;
t323 = qJD(4) + t215;
t309 = t235 * t240;
t322 = g(3) * t309 + t324;
t282 = qJDD(1) * t237;
t214 = -qJDD(4) + t282;
t321 = pkin(4) * t214;
t319 = pkin(4) * t240;
t225 = t237 * qJDD(2);
t314 = t203 * t235;
t187 = -t225 + t314;
t315 = t187 * t235;
t312 = t228 * t240;
t230 = t235 ^ 2;
t244 = qJD(1) ^ 2;
t311 = t230 * t244;
t308 = t235 * t242;
t307 = t237 * MDP(5);
t305 = t237 * t242;
t238 = cos(pkin(7));
t223 = -pkin(1) * t238 - pkin(2);
t200 = -pkin(3) * t237 - pkin(6) * t235 + t223;
t189 = qJD(1) * t200 + qJD(3);
t208 = t220 * qJD(1);
t196 = qJD(2) * t235 + t208 * t237;
t267 = t242 * t189 - t196 * t240;
t296 = qJD(1) * t242;
t273 = t235 * t296;
t170 = -qJ(5) * t273 + t267;
t167 = -pkin(4) * t215 + t170;
t304 = -t170 + t167;
t290 = t237 * qJD(3);
t293 = qJD(4) * t242;
t303 = t200 * t293 + t242 * t290;
t204 = t220 * t305;
t302 = t240 * t200 + t204;
t301 = t237 ^ 2 + t230;
t233 = t240 ^ 2;
t234 = t242 ^ 2;
t300 = t233 - t234;
t299 = MDP(17) * t235;
t297 = qJD(1) * t240;
t295 = qJD(4) * t196;
t294 = qJD(4) * t240;
t292 = qJD(5) * t235;
t291 = t214 * MDP(12);
t226 = t237 * qJD(2);
t179 = qJD(5) - t226 + (pkin(4) * t297 + t208) * t235;
t289 = qJD(5) + t179;
t288 = qJ(5) * qJDD(1);
t285 = qJD(1) * qJD(5);
t283 = qJDD(1) * t223;
t280 = qJDD(1) * t242;
t279 = MDP(13) + MDP(15);
t278 = MDP(14) + MDP(16);
t277 = qJ(5) * t308;
t276 = t240 * t311;
t243 = cos(qJ(1));
t275 = t243 * pkin(1) + t228 * pkin(2) + t227 * qJ(3);
t274 = t235 * t297;
t272 = t215 * t294;
t271 = t220 * t294;
t241 = sin(qJ(1));
t269 = -pkin(1) * t241 + t228 * qJ(3);
t266 = MDP(17) * (-t233 - t234);
t199 = (t220 + t319) * t235;
t265 = qJD(1) * t199 + t179;
t188 = qJDD(2) * t235 + t203 * t237;
t264 = -qJD(4) * t189 - t188;
t263 = t214 - t282;
t262 = t214 + t282;
t186 = qJDD(1) * t200 + qJDD(3);
t261 = t240 * t186 + t242 * t188 + t189 * t293 - t196 * t294;
t260 = qJD(4) * t274;
t259 = pkin(4) * t327 + qJDD(5) - t225;
t258 = -g(1) * t190 - g(2) * t192;
t191 = -t227 * t305 + t312;
t193 = t227 * t240 + t228 * t305;
t257 = -g(1) * t191 - g(2) * t193;
t256 = -g(1) * t228 - g(2) * t227;
t254 = g(1) * t241 - g(2) * t243;
t253 = t188 * t237 + t315;
t252 = -t240 * t189 - t242 * t196;
t195 = t208 * t235 - t226;
t251 = t195 * t235 + t196 * t237;
t177 = t259 + t314;
t202 = (pkin(4) * t293 + qJD(3)) * t235;
t250 = qJD(1) * t202 + qJDD(1) * t199 + t177;
t182 = t242 * t186;
t248 = qJ(5) * t260 + t264 * t240 + t182;
t247 = -t215 ^ 2 - t311;
t246 = g(1) * t193 - g(2) * t191 + g(3) * t308 - t261;
t212 = t235 * t280;
t206 = qJDD(3) + t283;
t205 = t237 * t260;
t201 = t215 * t273;
t198 = t242 * t200;
t178 = -qJ(5) * t309 + t302;
t176 = -t277 + t198 + (-t220 * t240 - pkin(4)) * t237;
t171 = -qJ(5) * t274 - t252;
t169 = -t240 * t290 - t242 * t292 + (-t204 + (qJ(5) * t235 - t200) * t240) * qJD(4);
t168 = -t240 * t292 + (-t220 * t306 - t277) * qJD(4) + t303;
t166 = (-qJ(5) * t326 - t240 * t285) * t235 + t261;
t165 = -t321 + (-t295 + (-t285 - t288) * t235) * t242 + t248;
t1 = [qJDD(1) * MDP(1) + t254 * MDP(2) + (g(1) * t243 + g(2) * t241) * MDP(3) + (t254 + (t236 ^ 2 + t238 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t206 - t268 - t283) * t307 + (t203 * t301 + t253 + t256) * MDP(6) + (t206 * t223 - g(1) * (-pkin(2) * t227 + t269) - g(2) * t275 + t253 * t220 + t251 * qJD(3)) * MDP(7) + (qJDD(1) * t234 - 0.2e1 * t240 * t270) * t230 * MDP(8) + 0.2e1 * (-t240 * t280 + t300 * t286) * t230 * MDP(9) + (t205 + (-t262 * t242 + t272) * t235) * MDP(10) + (t262 * t240 + (t215 + t298) * t293) * t235 * MDP(11) + t237 * t291 + (-t182 * t237 - t198 * t214 + ((qJD(1) * t230 + t215 * t237) * t220 + t251) * t293 + (-(-qJD(4) * t200 - t290) * t215 - t264 * t237 + t230 * t287 + t315 + (qJDD(1) * t230 + t214 * t237) * t220) * t240 + t257) * MDP(13) + ((-t237 * t271 + t303) * t215 + t302 * t214 + t261 * t237 + (t187 * t242 - t195 * t294) * t235 + (t220 * t280 + (qJD(3) * t242 - t271) * qJD(1)) * t230 + t258) * MDP(14) + (-t165 * t237 - t169 * t215 - t176 * t214 + (t250 * t240 + t265 * t293) * t235 + t257) * MDP(15) + (t166 * t237 + t168 * t215 + t178 * t214 + (t250 * t242 - t265 * t294) * t235 + t258) * MDP(16) + ((-qJD(4) * t171 - qJDD(1) * t176 - t165 + (-qJD(4) * t178 - t169) * qJD(1)) * t242 + (qJD(4) * t167 - qJDD(1) * t178 - t166 + (qJD(4) * t176 - t168) * qJD(1)) * t240 - t268) * t299 + (t166 * t178 + t171 * t168 + t165 * t176 + t167 * t169 + t177 * t199 + t179 * t202 - g(1) * (pkin(4) * t312 + t269) - g(2) * (t325 * t228 + t275) + (-g(1) * (-pkin(2) - t325) - g(2) * t319) * t227) * MDP(18); (qJDD(2) - g(3)) * MDP(4) + (-t187 * t237 - g(3)) * MDP(7) + (-t177 * t237 - g(3)) * MDP(18) + t278 * t205 + (t279 * (t263 * t240 + (t215 - t298) * t293) + t278 * (t263 * t242 - t272) + t188 * MDP(7) + (-t165 * t240 + t166 * t242 - t167 * t293 - t171 * t294) * MDP(18)) * t235; (qJDD(3) + t268) * MDP(7) + (t165 * t242 + t166 * t240 - t167 * t294 + t171 * t293 + t268) * MDP(18) - t301 * MDP(6) * t244 + t278 * (t240 * t214 + t247 * t242) + t279 * (-t214 * t242 + t240 * t247) + (-t251 * MDP(7) + (t167 * t306 - t171 * t305 - t179 * t235) * MDP(18)) * qJD(1) + (t223 * MDP(7) + t235 * t266 - t307) * qJDD(1); t242 * MDP(8) * t276 - t300 * MDP(9) * t311 + (-t323 * t274 + t212) * MDP(10) + (-t201 - t327) * MDP(11) - t291 + (-t240 * t188 - t195 * t273 + t323 * t252 + t182 + t322) * MDP(13) + (t195 * t274 - t267 * t215 + t246) * MDP(14) + (-0.2e1 * t321 - t171 * t215 + (-pkin(4) * t276 - t295 + (-t289 * qJD(1) - t288) * t235) * t242 + t248 + t322) * MDP(15) + (-pkin(4) * t234 * t311 - t170 * t215 + (qJ(5) * t281 + (qJ(5) * t293 + t289 * t240) * qJD(1)) * t235 + t246) * MDP(16) + (-pkin(4) * t280 + (pkin(4) * qJD(4) - t304) * t297) * t299 + (t304 * t171 + (t165 + (g(3) * t240 - t179 * t296) * t235 + t324) * pkin(4)) * MDP(18); -t201 * MDP(15) + t212 * MDP(16) + (g(3) * t237 + t259) * MDP(18) + t266 * t311 + (MDP(15) * t281 + (t203 + t256) * MDP(18) + ((MDP(15) * qJD(4) + t167 * MDP(18)) * t242 + ((-qJD(4) + t215) * MDP(16) + t171 * MDP(18)) * t240) * qJD(1)) * t235;];
tau = t1;
