% Calculate vector of inverse dynamics joint torques for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:15
% EndTime: 2019-12-05 15:22:19
% DurationCPUTime: 1.56s
% Computational Cost: add. (761->222), mult. (1606->314), div. (0->0), fcn. (1185->10), ass. (0->121)
t222 = pkin(7) + qJ(2);
t216 = cos(t222);
t210 = g(2) * t216;
t214 = sin(t222);
t291 = g(1) * t214;
t293 = -t210 + t291;
t289 = qJDD(2) * pkin(2);
t294 = t289 + t293;
t267 = qJD(2) * qJD(3);
t236 = qJ(3) * qJDD(2) + t267;
t224 = sin(pkin(8));
t292 = pkin(6) * t224;
t226 = cos(pkin(8));
t290 = qJ(3) * t226;
t225 = cos(pkin(9));
t228 = cos(qJ(5));
t279 = t225 * t228;
t261 = t224 * t279;
t199 = qJDD(2) * t261;
t227 = sin(qJ(5));
t280 = t225 * t227;
t223 = sin(pkin(9));
t281 = t223 * t228;
t239 = t280 + t281;
t232 = t239 * qJD(5);
t263 = qJDD(2) * t227;
t257 = t223 * t263;
t154 = t199 + (-qJD(2) * t232 - t257) * t224;
t288 = t154 * t226;
t275 = qJD(2) * t224;
t260 = t223 * t275;
t252 = t227 * t260;
t197 = qJD(5) * t252;
t266 = qJD(2) * qJD(5);
t259 = t228 * t266;
t251 = t225 * t259;
t155 = -t197 + (t239 * qJDD(2) + t251) * t224;
t287 = t155 * t226;
t175 = t239 * t275;
t274 = qJD(2) * t226;
t204 = -qJD(5) + t274;
t286 = t175 * t204;
t177 = qJD(2) * t261 - t252;
t285 = t177 * t204;
t284 = t214 * t226;
t283 = t216 * t226;
t282 = t223 * t227;
t179 = (-t224 * t282 + t261) * qJD(5);
t183 = t239 * t224;
t264 = qJDD(2) * t226;
t203 = -qJDD(5) + t264;
t278 = t179 * t204 + t183 * t203;
t198 = -pkin(3) * t226 - qJ(4) * t224 - pkin(2);
t271 = qJD(4) * t224;
t170 = -qJD(2) * t271 + t198 * qJDD(2) + qJDD(3);
t186 = qJDD(1) * t224 + t236 * t226;
t153 = t223 * t170 + t225 * t186;
t187 = t198 * qJD(2) + qJD(3);
t196 = qJ(3) * t274 + qJD(1) * t224;
t158 = t223 * t187 + t225 * t196;
t172 = t223 * t198 + t225 * t290;
t277 = t216 * pkin(2) + t214 * qJ(3);
t218 = t224 ^ 2;
t276 = t226 ^ 2 + t218;
t273 = qJD(3) * t224;
t272 = qJD(3) * t226;
t270 = t203 * MDP(17);
t269 = t228 * MDP(18);
t265 = qJDD(2) * t224;
t262 = t225 * t292;
t258 = t223 * t265;
t229 = qJD(2) ^ 2;
t255 = t276 * t229;
t152 = t225 * t170 - t186 * t223;
t237 = -pkin(4) * t226 - t262;
t149 = t237 * qJDD(2) + t152;
t150 = -pkin(6) * t258 + t153;
t254 = t228 * t149 - t227 * t150;
t157 = t225 * t187 - t196 * t223;
t253 = (-t223 ^ 2 - t225 ^ 2) * MDP(11);
t195 = -qJ(3) * t275 + qJD(1) * t226;
t250 = g(1) * t216 + g(2) * t214;
t193 = qJD(4) - t195;
t185 = -qJ(3) * t265 + qJDD(1) * t226 - t224 * t267;
t247 = t227 * t149 + t228 * t150;
t151 = t237 * qJD(2) + t157;
t156 = -pkin(6) * t260 + t158;
t246 = t151 * t228 - t156 * t227;
t245 = -t151 * t227 - t156 * t228;
t192 = t225 * t198;
t159 = -t262 + t192 + (-qJ(3) * t223 - pkin(4)) * t226;
t160 = -t223 * t292 + t172;
t244 = t159 * t228 - t160 * t227;
t243 = t159 * t227 + t160 * t228;
t178 = t224 * t232;
t238 = -t279 + t282;
t184 = t238 * t224;
t242 = t178 * t204 + t184 * t203;
t241 = -t185 * t224 + t186 * t226;
t240 = t195 * t224 - t196 * t226;
t182 = qJDD(4) - t185;
t171 = -t223 * t290 + t192;
t188 = -t223 * t272 - t225 * t271;
t235 = -t188 * qJD(2) - t171 * qJDD(2) - t152;
t189 = -t223 * t271 + t225 * t272;
t234 = t189 * qJD(2) + t172 * qJDD(2) + t153;
t233 = t204 * t238;
t212 = qJDD(3) - t289;
t231 = -t212 + t294;
t221 = pkin(9) + qJ(5);
t215 = cos(t221);
t213 = sin(t221);
t206 = t216 * qJ(3);
t194 = (pkin(4) * t223 + qJ(3)) * t224;
t174 = pkin(4) * t260 + t193;
t169 = t213 * t214 + t215 * t283;
t168 = -t213 * t283 + t214 * t215;
t167 = t213 * t216 - t215 * t284;
t166 = t213 * t284 + t215 * t216;
t162 = pkin(4) * t258 + t182;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t185 * t226 + t186 * t224 - g(3)) * MDP(8) + (-t182 * t226 - g(3) + (-t152 * t223 + t153 * t225) * t224) * MDP(12) + (t278 - t287) * MDP(18) + (-t242 - t288) * MDP(19); qJDD(2) * MDP(2) + t293 * MDP(3) + t250 * MDP(4) + (t236 * t276 + t241 - t250) * MDP(7) + (-t212 * pkin(2) - g(1) * (-pkin(2) * t214 + t206) - g(2) * t277 - t240 * qJD(3) + t241 * qJ(3)) * MDP(8) + (t153 * t172 + t158 * t189 + t152 * t171 + t157 * t188 - g(1) * t206 - g(2) * (pkin(3) * t283 + t277) - t198 * t291) * MDP(12) + (-t154 * t184 - t177 * t178) * MDP(13) + (-t154 * t183 + t155 * t184 + t175 * t178 - t177 * t179) * MDP(14) + (t242 - t288) * MDP(15) + (t278 + t287) * MDP(16) + (-g(1) * t167 - g(2) * t169 + t194 * t155 + t162 * t183 + t174 * t179 + t175 * t273 - t244 * t203) * MDP(18) + (-g(1) * t166 - g(2) * t168 + t194 * t154 - t162 * t184 - t174 * t178 + t177 * t273 + t243 * t203) * MDP(19) + (t225 * MDP(10) + t223 * MDP(9)) * (t182 * t224 + t236 * t218 - t250) + ((t243 * qJD(5) - t188 * t228 + t189 * t227) * MDP(18) + (t244 * qJD(5) + t188 * t227 + t189 * t228) * MDP(19)) * t204 + (t231 * MDP(5) + (t225 * t293 + t235) * MDP(9) + (-t223 * t293 + t234) * MDP(10) + t270 + (-t245 * qJD(5) - t254) * MDP(18) + (t246 * qJD(5) + t247) * MDP(19)) * t226 + (-t231 * MDP(6) + (-t234 * t223 + t235 * t225 + t293) * MDP(11) + (qJ(3) * t182 - qJ(4) * t210 + qJD(3) * t193) * MDP(12)) * t224; -MDP(5) * t264 - MDP(7) * t255 + (t240 * qJD(2) + qJDD(3) - t294) * MDP(8) + (-t223 * t255 - t225 * t264) * MDP(9) + (t223 * t264 - t225 * t255) * MDP(10) + (t152 * t225 + t153 * t223 + (-t193 * t224 + (t157 * t223 - t158 * t225) * t226) * qJD(2) - t293) * MDP(12) + (t238 * t203 + t204 * t232 + (-t239 * t204 * t226 - t224 * t175) * qJD(2)) * MDP(18) + (t239 * t203 - qJD(5) * t233 + (-t224 * t177 + t226 * t233) * qJD(2)) * MDP(19) + (MDP(6) + t253) * t265; (g(3) * t226 + t182) * MDP(12) + (-t197 - t285) * MDP(18) + (t199 + t286) * MDP(19) + t229 * t218 * t253 + (-t250 * MDP(12) + (MDP(10) * t223 - MDP(9) * t225) * t229 * t226 + ((t227 * MDP(18) + MDP(10)) * t225 + (-t227 * MDP(19) + MDP(9) + t269) * t223) * qJDD(2) + ((t157 * t225 + t158 * t223) * MDP(12) + (-t239 * MDP(19) + t225 * t269) * qJD(5)) * qJD(2)) * t224; t177 * t175 * MDP(13) + (-t175 ^ 2 + t177 ^ 2) * MDP(14) + (t199 - t286) * MDP(15) + (t197 - t285) * MDP(16) - t270 + (-g(1) * t168 + g(2) * t166 - t174 * t177 + t245 * t204 + t254) * MDP(18) + (g(1) * t169 - g(2) * t167 + t174 * t175 - t246 * t204 - t247) * MDP(19) + (t245 * MDP(18) - t246 * MDP(19)) * qJD(5) + ((-t223 * t259 - t266 * t280 - t257) * MDP(15) + (-qJDD(2) * t281 - t225 * t263 - t251) * MDP(16) + (MDP(18) * t213 + MDP(19) * t215) * g(3)) * t224;];
tau = t1;
