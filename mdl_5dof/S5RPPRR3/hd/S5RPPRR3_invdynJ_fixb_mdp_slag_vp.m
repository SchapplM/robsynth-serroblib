% Calculate vector of inverse dynamics joint torques for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:38
% EndTime: 2022-01-23 09:14:41
% DurationCPUTime: 1.48s
% Computational Cost: add. (1252->235), mult. (2642->309), div. (0->0), fcn. (2032->16), ass. (0->119)
t265 = cos(pkin(9));
t271 = cos(qJ(4));
t310 = t271 * t265;
t263 = sin(pkin(9));
t268 = sin(qJ(4));
t312 = t263 * t268;
t227 = -t310 + t312;
t219 = t227 * qJD(1);
t270 = cos(qJ(5));
t228 = t263 * t271 + t265 * t268;
t220 = t228 * qJD(1);
t267 = sin(qJ(5));
t313 = t220 * t267;
t182 = t270 * t219 + t313;
t261 = qJD(4) + qJD(5);
t315 = t182 * t261;
t262 = qJ(1) + pkin(8);
t253 = sin(t262);
t255 = cos(t262);
t289 = g(1) * t255 + g(2) * t253;
t264 = sin(pkin(8));
t242 = pkin(1) * t264 + qJ(3);
t235 = t242 * qJD(1);
t251 = t265 * qJD(2);
t201 = t251 + (-pkin(6) * qJD(1) - t235) * t263;
t208 = t263 * qJD(2) + t265 * t235;
t305 = qJD(1) * t265;
t202 = pkin(6) * t305 + t208;
t280 = -t201 * t268 - t202 * t271;
t173 = -pkin(7) * t219 - t280;
t266 = cos(pkin(8));
t246 = -pkin(1) * t266 - pkin(2);
t234 = -pkin(3) * t265 + t246;
t217 = t234 * qJD(1) + qJD(3);
t189 = pkin(4) * t219 + t217;
t260 = pkin(9) + qJ(4);
t256 = qJ(5) + t260;
t244 = sin(t256);
t245 = cos(t256);
t302 = qJD(5) * t267;
t328 = g(3) * t244 + t173 * t302 + t189 * t182 + t289 * t245;
t257 = qJDD(4) + qJDD(5);
t278 = -t219 * t267 + t270 * t220;
t327 = t257 * MDP(19) + t182 * MDP(15) * t278 + (-t182 ^ 2 + t278 ^ 2) * MDP(16);
t316 = t278 * t261;
t325 = t271 * t201 - t202 * t268;
t320 = pkin(6) + t242;
t223 = t320 * t263;
t224 = t320 * t265;
t309 = -t268 * t223 + t271 * t224;
t304 = qJD(1) * t268;
t294 = t263 * t304;
t297 = qJDD(1) * t271;
t298 = qJDD(1) * t268;
t303 = qJD(4) * t271;
t295 = t263 * t297 + t265 * t298 + t303 * t305;
t190 = -qJD(4) * t294 + t295;
t225 = qJD(1) * qJD(3) + t242 * qJDD(1);
t249 = t265 * qJDD(2);
t319 = pkin(6) * qJDD(1);
t197 = t249 + (-t225 - t319) * t263;
t204 = t263 * qJDD(2) + t265 * t225;
t198 = t265 * t319 + t204;
t292 = t271 * t197 - t268 * t198;
t165 = qJDD(4) * pkin(4) - pkin(7) * t190 + t280 * qJD(4) + t292;
t222 = t228 * qJD(4);
t239 = t265 * t297;
t286 = -t263 * t298 + t239;
t191 = qJD(1) * t222 - t286;
t281 = t268 * t197 + t271 * t198;
t166 = -pkin(7) * t191 + t325 * qJD(4) + t281;
t324 = -g(3) * t245 + t270 * t165 - t267 * t166 - t189 * t278 + t289 * t244;
t293 = t190 * t267 + t270 * t191;
t169 = t278 * qJD(5) + t293;
t323 = pkin(4) * t222;
t172 = -pkin(7) * t220 + t325;
t171 = qJD(4) * pkin(4) + t172;
t318 = t171 * t270;
t317 = t173 * t270;
t311 = t265 * MDP(5);
t308 = t263 ^ 2 + t265 ^ 2;
t306 = qJD(1) * t263;
t301 = qJD(5) * t270;
t299 = qJDD(1) * t246;
t296 = t270 * t190 - t267 * t191 - t219 * t301;
t291 = -t271 * t223 - t224 * t268;
t288 = g(1) * t253 - g(2) * t255;
t269 = sin(qJ(1));
t272 = cos(qJ(1));
t287 = g(1) * t269 - g(2) * t272;
t285 = -t171 * t267 - t317;
t194 = t270 * t227 + t228 * t267;
t221 = t227 * qJD(4);
t174 = -t194 * qJD(5) - t221 * t270 - t222 * t267;
t195 = -t227 * t267 + t228 * t270;
t284 = t174 * t261 + t195 * t257;
t179 = -pkin(7) * t228 + t291;
t180 = -pkin(7) * t227 + t309;
t283 = t179 * t270 - t180 * t267;
t282 = t179 * t267 + t180 * t270;
t203 = -t225 * t263 + t249;
t279 = -t203 * t263 + t204 * t265;
t168 = -t220 * t302 + t296;
t215 = t234 * qJDD(1) + qJDD(3);
t276 = -t223 * t303 + qJD(3) * t310 + (-qJD(3) * t263 - qJD(4) * t224) * t268;
t275 = -t228 * qJD(3) - t309 * qJD(4);
t254 = cos(t260);
t252 = sin(t260);
t233 = qJDD(3) + t299;
t207 = -t235 * t263 + t251;
t200 = pkin(4) * t227 + t234;
t193 = -qJD(4) * t222 - qJDD(4) * t227;
t192 = -qJD(4) * t221 + qJDD(4) * t228;
t178 = pkin(4) * t191 + t215;
t177 = pkin(7) * t221 + t275;
t176 = -pkin(7) * t222 + t276;
t175 = t195 * qJD(5) - t221 * t267 + t270 * t222;
t167 = -t175 * t261 - t194 * t257;
t1 = [qJDD(1) * MDP(1) + t287 * MDP(2) + (g(1) * t272 + g(2) * t269) * MDP(3) + (t287 + (t264 ^ 2 + t266 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t233 + t288 - t299) * t311 + (t225 * t308 + t279 - t289) * MDP(6) + (t233 * t246 - g(1) * (-pkin(1) * t269 - pkin(2) * t253 + qJ(3) * t255) - g(2) * (pkin(1) * t272 + pkin(2) * t255 + qJ(3) * t253) + t279 * t242 + (-t207 * t263 + t208 * t265) * qJD(3)) * MDP(7) + (t190 * t228 - t220 * t221) * MDP(8) + (-t190 * t227 - t191 * t228 + t219 * t221 - t220 * t222) * MDP(9) + t192 * MDP(10) + t193 * MDP(11) + (t275 * qJD(4) + t291 * qJDD(4) + t234 * t191 + t215 * t227 + t217 * t222 + t288 * t254) * MDP(13) + (-t276 * qJD(4) - t309 * qJDD(4) + t234 * t190 + t215 * t228 - t217 * t221 - t288 * t252) * MDP(14) + (t168 * t195 + t174 * t278) * MDP(15) + (-t168 * t194 - t169 * t195 - t174 * t182 - t175 * t278) * MDP(16) + t284 * MDP(17) + t167 * MDP(18) + (t182 * t323 + t200 * t169 + t178 * t194 + t189 * t175 + (-t282 * qJD(5) - t176 * t267 + t177 * t270) * t261 + t283 * t257 + t288 * t245) * MDP(20) + (t278 * t323 + t200 * t168 + t178 * t195 + t189 * t174 - (t283 * qJD(5) + t176 * t270 + t177 * t267) * t261 - t282 * t257 - t288 * t244) * MDP(21); (qJDD(2) - g(3)) * MDP(4) + (t203 * t265 + t204 * t263 - g(3)) * MDP(7) + t193 * MDP(13) - t192 * MDP(14) + t167 * MDP(20) - t284 * MDP(21); (t207 * t306 - t208 * t305 + qJDD(3) - t288) * MDP(7) - t239 * MDP(13) + t295 * MDP(14) + (t169 + t316) * MDP(20) + (t168 - t315) * MDP(21) - t308 * MDP(6) * qJD(1) ^ 2 + (MDP(13) * t312 + t246 * MDP(7) - t311) * qJDD(1) + ((t265 * t304 + t271 * t306 + t220) * MDP(13) + (-t219 - t294) * MDP(14)) * qJD(4); t220 * t219 * MDP(8) + (-t219 ^ 2 + t220 ^ 2) * MDP(9) + (t295 + (t219 - t294) * qJD(4)) * MDP(10) + t286 * MDP(11) + qJDD(4) * MDP(12) + (-g(3) * t254 - t217 * t220 + t289 * t252 + t292) * MDP(13) + (g(3) * t252 + t217 * t219 + t289 * t254 - t281) * MDP(14) + (t168 + t315) * MDP(17) + (-t169 + t316) * MDP(18) + (-(-t172 * t267 - t317) * t261 + t285 * qJD(5) + (-t182 * t220 + t257 * t270 - t261 * t302) * pkin(4) + t324) * MDP(20) + ((-t173 * t261 - t165) * t267 + (-qJD(5) * t171 + t172 * t261 - t166) * t270 + (-t220 * t278 - t257 * t267 - t261 * t301) * pkin(4) + t328) * MDP(21) + t327; (t296 + t315) * MDP(17) + (-t293 + t316) * MDP(18) + (-t285 * t261 + t324) * MDP(20) + (-t270 * t166 - t267 * t165 + (-t173 * t267 + t318) * t261 + t328) * MDP(21) + (-MDP(17) * t313 - t278 * MDP(18) + t285 * MDP(20) - MDP(21) * t318) * qJD(5) + t327;];
tau = t1;
