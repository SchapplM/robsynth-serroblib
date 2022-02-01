% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:57
% EndTime: 2022-01-23 09:20:59
% DurationCPUTime: 1.25s
% Computational Cost: add. (1161->207), mult. (1928->290), div. (0->0), fcn. (1143->14), ass. (0->112)
t224 = cos(pkin(8));
t208 = pkin(1) * t224 + pkin(2);
t222 = sin(pkin(8));
t295 = pkin(1) * t222;
t263 = qJD(3) * t295;
t299 = -qJD(1) * t263 + t208 * qJDD(1);
t193 = t208 * qJD(1);
t298 = qJD(3) * t193 + qJDD(1) * t295;
t217 = qJD(1) + qJD(3);
t214 = qJDD(1) + qJDD(3);
t292 = qJ(4) * t214;
t226 = sin(qJ(3));
t229 = cos(qJ(3));
t296 = -t299 * t226 - t298 * t229;
t162 = qJD(4) * t217 + t292 - t296;
t221 = sin(pkin(9));
t223 = cos(pkin(9));
t159 = qJDD(2) * t221 + t162 * t223;
t158 = -t223 * qJDD(2) + t162 * t221;
t291 = t158 * t221;
t297 = t159 * t223 + t291;
t274 = t226 * t208 + t229 * t295;
t225 = sin(qJ(5));
t228 = cos(qJ(5));
t242 = MDP(16) * t225 + MDP(17) * t228;
t264 = qJD(1) * t295;
t177 = t193 * t229 - t226 * t264;
t248 = qJD(4) - t177;
t294 = pkin(3) * t214;
t218 = qJ(1) + pkin(8);
t212 = qJ(3) + t218;
t206 = sin(t212);
t202 = g(1) * t206;
t207 = cos(t212);
t293 = g(2) * t207;
t178 = t193 * t226 + t229 * t264;
t170 = qJ(4) * t217 + t178;
t167 = -t223 * qJD(2) + t170 * t221;
t290 = t167 * t221;
t187 = -pkin(4) * t223 - pkin(7) * t221 - pkin(3);
t284 = t208 * t229;
t253 = -t226 * t295 + t284;
t171 = t187 - t253;
t289 = t171 * t228;
t288 = t178 * t217;
t241 = qJD(3) * t284 - t226 * t263;
t179 = qJD(4) + t241;
t287 = t179 * t217;
t286 = t179 * t225;
t180 = t274 * qJD(3);
t285 = t180 * t217;
t283 = t214 * t223;
t282 = t214 * t225;
t281 = t217 * t223;
t280 = t221 * t225;
t279 = t223 * t225;
t278 = t223 * t228;
t277 = t225 * t228;
t276 = t207 * pkin(3) + t206 * qJ(4);
t275 = -g(1) * t207 - g(2) * t206;
t215 = t221 ^ 2;
t216 = t223 ^ 2;
t273 = t215 + t216;
t220 = t228 ^ 2;
t272 = t225 ^ 2 - t220;
t268 = qJD(5) * t225;
t267 = qJD(5) * t228;
t192 = -qJDD(5) + t283;
t266 = t192 * MDP(15);
t195 = -qJD(5) + t281;
t265 = -qJD(5) - t195;
t261 = t217 * t267;
t260 = -t206 * pkin(3) + t207 * qJ(4);
t250 = -t298 * t226 + t299 * t229;
t243 = qJDD(4) - t250;
t163 = t243 - t294;
t259 = -t163 - t293;
t258 = t273 * t214;
t157 = t187 * t214 + t243;
t257 = t228 * t157 - t225 * t159;
t255 = t192 - t283;
t254 = t192 + t283;
t252 = t248 * t228;
t227 = sin(qJ(1));
t230 = cos(qJ(1));
t249 = g(1) * t227 - g(2) * t230;
t247 = t225 * t157 + t228 * t159;
t164 = t187 * t217 + t248;
t168 = qJD(2) * t221 + t170 * t223;
t246 = t164 * t228 - t168 * t225;
t245 = -t164 * t225 - t168 * t228;
t244 = t168 * t223 + t290;
t240 = t275 + t297;
t239 = -qJ(4) * t279 + t187 * t228;
t173 = t206 * t279 + t207 * t228;
t175 = t206 * t228 - t207 * t279;
t237 = -g(1) * t173 - g(2) * t175 + (t246 * qJD(5) + t247) * t223 + t228 * t291;
t174 = -t206 * t278 + t207 * t225;
t176 = t206 * t225 + t207 * t278;
t236 = -g(1) * t174 - g(2) * t176 + t158 * t280 + t267 * t290;
t235 = -t202 - t250 + t293;
t234 = qJ(4) * t267 + t248 * t225;
t184 = t221 * t268 * t281;
t233 = t223 * t266 + (t184 + (t195 * t268 - t254 * t228) * t221) * MDP(13) + (t254 * t225 + (t195 + t281) * t267) * t221 * MDP(14) + t214 * MDP(5) + (0.2e1 * (t272 * t217 * qJD(5) - t214 * t277) * MDP(12) + (t214 * t220 - 0.2e1 * t225 * t261) * MDP(11)) * t215;
t232 = -t275 + t296;
t213 = t217 ^ 2;
t188 = t223 * t202;
t182 = -pkin(3) - t253;
t181 = qJ(4) + t274;
t169 = -pkin(3) * t217 + t248;
t152 = t245 * qJD(5) + t257;
t1 = [qJDD(1) * MDP(1) + t249 * MDP(2) + (g(1) * t230 + g(2) * t227) * MDP(3) + (t249 + (t222 ^ 2 + t224 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t253 * t214 - t235 - t285) * MDP(6) + (-t274 * t214 - t241 * t217 + t232) * MDP(7) + (t188 + (-t182 * t214 + t259 - t285) * t223) * MDP(8) + (t181 * t258 + t273 * t287 + t240) * MDP(9) + (t163 * t182 + t169 * t180 - g(1) * (-pkin(2) * sin(t218) - t227 * pkin(1) + t260) - g(2) * (pkin(2) * cos(t218) + t230 * pkin(1) + t276) + t297 * t181 + t244 * t179) * MDP(10) + (-(-t171 * t268 + t180 * t228) * t195 - t192 * t289 + (-(-t181 * t267 - t286) * t195 + t181 * t225 * t192 - t152) * t223 + (t217 * t286 + (t261 + t282) * t181) * t215 + t236) * MDP(16) + ((t179 * t278 + t180 * t225) * t195 + (t171 * t225 + t181 * t278) * t192 + (t181 * t214 + t287) * t228 * t215 + (t195 * t289 + (-t290 + (-t195 * t223 - t215 * t217) * t181) * t225) * qJD(5) + t237) * MDP(17) + t233; (qJDD(2) - g(3)) * MDP(4) + (-t158 * t223 - g(3)) * MDP(10) + t184 * MDP(17) + (t159 * MDP(10) + (-qJD(5) * t195 * MDP(17) + t255 * MDP(16)) * t225 + (t255 * MDP(17) + (t195 - t281) * MDP(16) * qJD(5)) * t228) * t221; (-t235 + t288) * MDP(6) + (t177 * t217 + t232) * MDP(7) + (t188 + (t259 + t288 + t294) * t223) * MDP(8) + (t248 * t217 * t273 + qJ(4) * t258 + t240) * MDP(9) + (-t163 * pkin(3) - t169 * t178 - g(1) * t260 - g(2) * t276 + (t159 * qJ(4) + t248 * t168) * t223 + (t158 * qJ(4) + t248 * t167) * t221) * MDP(10) + (-t239 * t192 - t152 * t223 + (t228 * t178 + t187 * t268 + t234 * t223) * t195 + (qJ(4) * t282 + t234 * t217) * t215 + t236) * MDP(16) + ((qJ(4) * t278 + t187 * t225) * t192 + (-t225 * t178 + t223 * t252) * t195 + (-t167 * t280 + t239 * t195) * qJD(5) + (t228 * t292 + (-qJ(4) * t268 + t252) * t217) * t215 + t237) * MDP(17) + t233; -MDP(8) * t283 + (-t244 * t217 - t202 - t259) * MDP(10) + (-MDP(16) * t228 + MDP(17) * t225) * t192 + (-t216 * MDP(9) + (-MDP(9) - t242) * t215) * t213 - t242 * t195 ^ 2; -t266 + (-g(1) * t175 + g(2) * t173 + t245 * t195 + t257) * MDP(16) + (g(1) * t176 - g(2) * t174 - t246 * t195 - t247) * MDP(17) + (t245 * MDP(16) - t246 * MDP(17)) * qJD(5) + (MDP(11) * t277 - t272 * MDP(12)) * t215 * t213 + ((MDP(13) * t228 - MDP(14) * t225) * t214 + t242 * g(3) + ((t265 * MDP(14) - t167 * MDP(16)) * t228 + (t265 * MDP(13) + t167 * MDP(17)) * t225) * t217) * t221;];
tau = t1;
