% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:12
% EndTime: 2019-12-05 17:29:17
% DurationCPUTime: 1.84s
% Computational Cost: add. (904->230), mult. (1830->329), div. (0->0), fcn. (1301->14), ass. (0->128)
t231 = cos(pkin(7));
t300 = pkin(1) * t231;
t215 = -pkin(2) - t300;
t274 = qJDD(1) * t215;
t206 = qJDD(3) + t274;
t225 = qJ(1) + pkin(7);
t217 = sin(t225);
t219 = cos(t225);
t298 = g(2) * t219;
t258 = g(3) * t217 + t298;
t301 = -t206 + t258;
t228 = sin(pkin(7));
t214 = pkin(1) * t228 + qJ(3);
t202 = qJD(1) * qJD(3) + qJDD(1) * t214;
t227 = sin(pkin(8));
t299 = pkin(6) * t227;
t235 = cos(qJ(1));
t297 = g(2) * t235;
t229 = cos(pkin(9));
t234 = cos(qJ(5));
t287 = t229 * t234;
t269 = t227 * t287;
t207 = qJDD(1) * t269;
t232 = sin(qJ(5));
t288 = t229 * t232;
t226 = sin(pkin(9));
t289 = t226 * t234;
t247 = t288 + t289;
t240 = t247 * qJD(5);
t272 = qJDD(1) * t232;
t265 = t226 * t272;
t160 = t207 + (-qJD(1) * t240 - t265) * t227;
t230 = cos(pkin(8));
t296 = t160 * t230;
t283 = qJD(1) * t227;
t268 = t226 * t283;
t260 = t232 * t268;
t205 = qJD(5) * t260;
t276 = qJD(1) * qJD(5);
t267 = t234 * t276;
t259 = t229 * t267;
t161 = -t205 + (t247 * qJDD(1) + t259) * t227;
t295 = t161 * t230;
t185 = t247 * t283;
t212 = t230 * qJD(1) - qJD(5);
t294 = t185 * t212;
t187 = qJD(1) * t269 - t260;
t293 = t187 * t212;
t292 = t217 * t230;
t291 = t219 * t230;
t290 = t226 * t232;
t286 = t230 * t214;
t189 = (-t227 * t290 + t269) * qJD(5);
t193 = t247 * t227;
t271 = t230 * qJDD(1);
t211 = -qJDD(5) + t271;
t285 = t189 * t212 + t193 * t211;
t244 = pkin(3) * t230 + qJ(4) * t227 + pkin(2);
t201 = -t244 - t300;
t280 = qJD(4) * t227;
t169 = -qJD(1) * t280 + t201 * qJDD(1) + qJDD(3);
t182 = qJDD(2) * t227 + t202 * t230;
t156 = t226 * t169 + t229 * t182;
t183 = t201 * qJD(1) + qJD(3);
t210 = t214 * qJD(1);
t192 = qJD(2) * t227 + t210 * t230;
t159 = t226 * t183 + t229 * t192;
t168 = t226 * t201 + t229 * t286;
t221 = t227 ^ 2;
t284 = t230 ^ 2 + t221;
t282 = qJD(3) * t227;
t281 = qJD(3) * t230;
t279 = t211 * MDP(17);
t278 = t234 * MDP(18);
t273 = qJDD(1) * t227;
t270 = t229 * t299;
t266 = t226 * t273;
t233 = sin(qJ(1));
t264 = -pkin(1) * t233 + t219 * qJ(3);
t236 = qJD(1) ^ 2;
t263 = t284 * t236;
t155 = t229 * t169 - t182 * t226;
t245 = -pkin(4) * t230 - t270;
t152 = t245 * qJDD(1) + t155;
t153 = -pkin(6) * t266 + t156;
t262 = t234 * t152 - t232 * t153;
t158 = t229 * t183 - t192 * t226;
t261 = (-t226 ^ 2 - t229 ^ 2) * MDP(11);
t191 = qJD(2) * t230 - t227 * t210;
t181 = qJDD(2) * t230 - t227 * t202;
t257 = g(2) * t217 - g(3) * t219;
t256 = g(3) * t233 + t297;
t190 = qJD(4) - t191;
t255 = t232 * t152 + t234 * t153;
t154 = t245 * qJD(1) + t158;
t157 = -pkin(6) * t268 + t159;
t254 = t154 * t234 - t157 * t232;
t253 = -t154 * t232 - t157 * t234;
t196 = t229 * t201;
t162 = -t270 + t196 + (-t214 * t226 - pkin(4)) * t230;
t163 = -t226 * t299 + t168;
t252 = t162 * t234 - t163 * t232;
t251 = t162 * t232 + t163 * t234;
t250 = -t181 * t227 + t182 * t230;
t188 = t227 * t240;
t246 = -t287 + t290;
t194 = t246 * t227;
t249 = t188 * t212 + t194 * t211;
t248 = t191 * t227 - t192 * t230;
t176 = qJDD(4) - t181;
t167 = -t226 * t286 + t196;
t199 = -t226 * t281 - t229 * t280;
t243 = -t199 * qJD(1) - t167 * qJDD(1) - t155;
t200 = -t226 * t280 + t229 * t281;
t242 = t200 * qJD(1) + t168 * qJDD(1) + t156;
t241 = t212 * t246;
t239 = -t274 + t301;
t224 = pkin(9) + qJ(5);
t218 = cos(t224);
t216 = sin(t224);
t197 = (pkin(4) * t226 + t214) * t227;
t180 = -t216 * t217 - t218 * t291;
t179 = t216 * t291 - t217 * t218;
t178 = -t216 * t219 + t218 * t292;
t177 = t216 * t292 + t218 * t219;
t171 = pkin(4) * t268 + t190;
t166 = pkin(4) * t266 + t176;
t1 = [qJDD(1) * MDP(1) + t256 * MDP(2) + (-g(2) * t233 + g(3) * t235) * MDP(3) + (t256 + (t228 ^ 2 + t231 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t202 * t284 + t250 + t257) * MDP(7) + (t206 * t215 - g(2) * (-pkin(1) * t235 - pkin(2) * t219 - qJ(3) * t217) - g(3) * (-pkin(2) * t217 + t264) + t250 * t214 - t248 * qJD(3)) * MDP(8) + (t156 * t168 + t159 * t200 + t155 * t167 + t158 * t199 + pkin(1) * t297 - g(3) * t264 + t244 * t298 + (g(2) * qJ(3) + g(3) * t244) * t217) * MDP(12) + (-t160 * t194 - t187 * t188) * MDP(13) + (-t160 * t193 + t161 * t194 + t185 * t188 - t187 * t189) * MDP(14) + (t249 - t296) * MDP(15) + (t285 + t295) * MDP(16) + (-g(2) * t180 + g(3) * t178 + t197 * t161 + t166 * t193 + t171 * t189 + t185 * t282 - t252 * t211) * MDP(18) + (-g(2) * t179 - g(3) * t177 + t197 * t160 - t166 * t194 - t171 * t188 + t187 * t282 + t251 * t211) * MDP(19) + (t229 * MDP(10) + t226 * MDP(9)) * (t176 * t227 + t202 * t221 + t257) + ((t251 * qJD(5) - t199 * t234 + t200 * t232) * MDP(18) + (t252 * qJD(5) + t199 * t232 + t200 * t234) * MDP(19)) * t212 + (t239 * MDP(5) + (t258 * t229 + t243) * MDP(9) + (-t258 * t226 + t242) * MDP(10) + t279 + (-t253 * qJD(5) - t262) * MDP(18) + (t254 * qJD(5) + t255) * MDP(19)) * t230 + (-t239 * MDP(6) + (-t242 * t226 + t243 * t229 + t258) * MDP(11) + (qJD(3) * t190 + t176 * t214) * MDP(12)) * t227; (qJDD(2) - g(1)) * MDP(4) + (t181 * t230 + t182 * t227 - g(1)) * MDP(8) + (-t176 * t230 - g(1) + (-t155 * t226 + t156 * t229) * t227) * MDP(12) + (t285 - t295) * MDP(18) + (-t249 - t296) * MDP(19); -MDP(5) * t271 - MDP(7) * t263 + (t248 * qJD(1) - t301) * MDP(8) + (-t226 * t263 - t229 * t271) * MDP(9) + (t226 * t271 - t229 * t263) * MDP(10) + (t155 * t229 + t156 * t226 + (-t190 * t227 + (t158 * t226 - t159 * t229) * t230) * qJD(1) - t258) * MDP(12) + (t246 * t211 + t212 * t240 + (-t247 * t212 * t230 - t227 * t185) * qJD(1)) * MDP(18) + (t247 * t211 - qJD(5) * t241 + (-t227 * t187 + t230 * t241) * qJD(1)) * MDP(19) + (MDP(6) + t261) * t273; (g(1) * t230 + t176) * MDP(12) + (-t205 - t293) * MDP(18) + (t207 + t294) * MDP(19) + t236 * t221 * t261 + (t257 * MDP(12) + (MDP(10) * t226 - MDP(9) * t229) * t236 * t230 + ((t232 * MDP(18) + MDP(10)) * t229 + (-t232 * MDP(19) + MDP(9) + t278) * t226) * qJDD(1) + ((t158 * t229 + t159 * t226) * MDP(12) + (-t247 * MDP(19) + t229 * t278) * qJD(5)) * qJD(1)) * t227; t187 * t185 * MDP(13) + (-t185 ^ 2 + t187 ^ 2) * MDP(14) + (t207 - t294) * MDP(15) + (t205 - t293) * MDP(16) - t279 + (-g(2) * t177 + g(3) * t179 - t171 * t187 + t253 * t212 + t262) * MDP(18) + (-g(2) * t178 - g(3) * t180 + t171 * t185 - t254 * t212 - t255) * MDP(19) + (t253 * MDP(18) - t254 * MDP(19)) * qJD(5) + ((-t226 * t267 - t276 * t288 - t265) * MDP(15) + (-qJDD(1) * t289 - t229 * t272 - t259) * MDP(16) + (MDP(18) * t216 + MDP(19) * t218) * g(1)) * t227;];
tau = t1;
