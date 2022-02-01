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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
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
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:48
% EndTime: 2022-01-20 09:12:52
% DurationCPUTime: 1.93s
% Computational Cost: add. (850->222), mult. (1721->317), div. (0->0), fcn. (1228->14), ass. (0->124)
t231 = cos(pkin(7));
t296 = pkin(1) * t231;
t216 = -pkin(2) - t296;
t272 = qJDD(1) * t216;
t204 = qJDD(3) + t272;
t225 = qJ(1) + pkin(7);
t220 = cos(t225);
t215 = g(2) * t220;
t218 = sin(t225);
t294 = g(1) * t218;
t260 = t215 - t294;
t297 = -t204 - t260;
t228 = sin(pkin(7));
t213 = pkin(1) * t228 + qJ(3);
t200 = qJD(1) * qJD(3) + qJDD(1) * t213;
t227 = sin(pkin(8));
t295 = pkin(6) * t227;
t229 = cos(pkin(9));
t234 = cos(qJ(5));
t284 = t229 * t234;
t267 = t227 * t284;
t205 = qJDD(1) * t267;
t226 = sin(pkin(9));
t232 = sin(qJ(5));
t285 = t229 * t232;
t244 = t226 * t234 + t285;
t239 = t244 * qJD(5);
t269 = qJDD(1) * t232;
t262 = t226 * t269;
t158 = t205 + (-qJD(1) * t239 - t262) * t227;
t230 = cos(pkin(8));
t293 = t158 * t230;
t281 = qJD(1) * t227;
t265 = t226 * t281;
t257 = t232 * t265;
t203 = qJD(5) * t257;
t274 = qJD(1) * qJD(5);
t264 = t234 * t274;
t256 = t229 * t264;
t159 = -t203 + (t244 * qJDD(1) + t256) * t227;
t292 = t159 * t230;
t183 = t244 * t281;
t210 = qJD(1) * t230 - qJD(5);
t291 = t183 * t210;
t185 = qJD(1) * t267 - t257;
t290 = t185 * t210;
t289 = t213 * t230;
t288 = t218 * t230;
t287 = t220 * t230;
t286 = t226 * t232;
t187 = (-t227 * t286 + t267) * qJD(5);
t191 = t244 * t227;
t270 = qJDD(1) * t230;
t209 = -qJDD(5) + t270;
t283 = t187 * t210 + t191 * t209;
t241 = -pkin(3) * t230 - qJ(4) * t227 - pkin(2);
t199 = t241 - t296;
t278 = qJD(4) * t227;
t167 = -qJD(1) * t278 + t199 * qJDD(1) + qJDD(3);
t180 = qJDD(2) * t227 + t200 * t230;
t154 = t226 * t167 + t229 * t180;
t181 = t199 * qJD(1) + qJD(3);
t208 = t213 * qJD(1);
t190 = qJD(2) * t227 + t208 * t230;
t157 = t226 * t181 + t229 * t190;
t166 = t226 * t199 + t229 * t289;
t222 = t227 ^ 2;
t282 = t230 ^ 2 + t222;
t280 = qJD(3) * t227;
t279 = qJD(3) * t230;
t277 = t209 * MDP(15);
t276 = t234 * MDP(16);
t271 = qJDD(1) * t226;
t268 = t229 * t295;
t235 = cos(qJ(1));
t266 = t235 * pkin(1) + t220 * pkin(2) + t218 * qJ(3);
t263 = t227 * t271;
t233 = sin(qJ(1));
t261 = -pkin(1) * t233 + t220 * qJ(3);
t236 = qJD(1) ^ 2;
t259 = t282 * t236;
t153 = t229 * t167 - t180 * t226;
t242 = -pkin(4) * t230 - t268;
t150 = t242 * qJDD(1) + t153;
t151 = -pkin(6) * t263 + t154;
t258 = t234 * t150 - t232 * t151;
t156 = t229 * t181 - t190 * t226;
t189 = qJD(2) * t230 - t227 * t208;
t179 = qJDD(2) * t230 - t227 * t200;
t255 = -g(1) * t220 - g(2) * t218;
t253 = g(1) * t233 - g(2) * t235;
t188 = qJD(4) - t189;
t252 = t232 * t150 + t234 * t151;
t152 = t242 * qJD(1) + t156;
t155 = -pkin(6) * t265 + t157;
t251 = t152 * t234 - t155 * t232;
t250 = -t152 * t232 - t155 * t234;
t194 = t229 * t199;
t160 = -t268 + t194 + (-t213 * t226 - pkin(4)) * t230;
t161 = -t226 * t295 + t166;
t249 = t160 * t234 - t161 * t232;
t248 = t160 * t232 + t161 * t234;
t247 = -t179 * t227 + t180 * t230;
t186 = t227 * t239;
t243 = -t284 + t286;
t192 = t243 * t227;
t246 = t186 * t210 + t192 * t209;
t245 = t189 * t227 - t190 * t230;
t174 = qJDD(4) - t179;
t240 = t210 * t243;
t224 = pkin(9) + qJ(5);
t219 = cos(t224);
t217 = sin(t224);
t198 = -t226 * t278 + t229 * t279;
t197 = -t226 * t279 - t229 * t278;
t195 = (pkin(4) * t226 + t213) * t227;
t178 = t217 * t218 + t219 * t287;
t177 = -t217 * t287 + t218 * t219;
t176 = t217 * t220 - t219 * t288;
t175 = t217 * t288 + t219 * t220;
t169 = pkin(4) * t265 + t188;
t165 = -t226 * t289 + t194;
t164 = pkin(4) * t263 + t174;
t1 = [qJDD(1) * MDP(1) + t253 * MDP(2) + (g(1) * t235 + g(2) * t233) * MDP(3) + (t253 + (t228 ^ 2 + t231 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t200 * t282 + t247 + t255) * MDP(6) + (t204 * t216 - g(1) * (-pkin(2) * t218 + t261) - g(2) * t266 + t247 * t213 - t245 * qJD(3)) * MDP(7) + (t154 * t166 + t157 * t198 + t153 * t165 + t156 * t197 - g(1) * t261 - g(2) * (pkin(3) * t287 + t266) + (-qJ(4) * t215 + t188 * qJD(3) + t174 * t213) * t227 - t241 * t294) * MDP(10) + (-t158 * t192 - t185 * t186) * MDP(11) + (-t158 * t191 + t159 * t192 + t183 * t186 - t185 * t187) * MDP(12) + (t246 - t293) * MDP(13) + (t283 + t292) * MDP(14) + (-g(1) * t176 - g(2) * t178 + t195 * t159 + t164 * t191 + t169 * t187 + t183 * t280 - t249 * t209) * MDP(16) + (-g(1) * t175 - g(2) * t177 + t195 * t158 - t164 * t192 - t169 * t186 + t185 * t280 + t248 * t209) * MDP(17) + (t226 * MDP(8) + t229 * MDP(9)) * (t174 * t227 + t200 * t222 + t255) + ((t248 * qJD(5) - t197 * t234 + t198 * t232) * MDP(16) + (t249 * qJD(5) + t197 * t232 + t198 * t234) * MDP(17)) * t210 + ((-t272 + t297) * MDP(5) + (-t197 * qJD(1) - t165 * qJDD(1) - t229 * t260 - t153) * MDP(8) + (t198 * qJD(1) + t166 * qJDD(1) + t226 * t260 + t154) * MDP(9) + t277 + (-qJD(5) * t250 - t258) * MDP(16) + (qJD(5) * t251 + t252) * MDP(17)) * t230; (qJDD(2) - g(3)) * MDP(4) + (t179 * t230 + t180 * t227 - g(3)) * MDP(7) + (-t174 * t230 - g(3) + (-t153 * t226 + t154 * t229) * t227) * MDP(10) + (t283 - t292) * MDP(16) + (-t246 - t293) * MDP(17); -MDP(5) * t270 - MDP(6) * t259 + (t245 * qJD(1) - t297) * MDP(7) + (-t226 * t259 - t229 * t270) * MDP(8) + (t226 * t270 - t229 * t259) * MDP(9) + (t153 * t229 + t154 * t226 + (-t188 * t227 + (t156 * t226 - t157 * t229) * t230) * qJD(1) + t260) * MDP(10) + (t243 * t209 + t210 * t239 + (-t244 * t210 * t230 - t227 * t183) * qJD(1)) * MDP(16) + (t244 * t209 - qJD(5) * t240 + (-t227 * t185 + t230 * t240) * qJD(1)) * MDP(17); (g(3) * t230 + t174) * MDP(10) + (-t203 - t290) * MDP(16) + (t205 + t291) * MDP(17) + (t255 * MDP(10) + (-MDP(8) * t229 + MDP(9) * t226) * t236 * t230 + ((t232 * MDP(16) + MDP(9)) * t229 + (-t232 * MDP(17) + MDP(8) + t276) * t226) * qJDD(1) + ((t156 * t229 + t157 * t226) * MDP(10) + (-t244 * MDP(17) + t229 * t276) * qJD(5)) * qJD(1)) * t227; t185 * t183 * MDP(11) + (-t183 ^ 2 + t185 ^ 2) * MDP(12) + (t205 - t291) * MDP(13) + (t203 - t290) * MDP(14) - t277 + (-g(1) * t177 + g(2) * t175 - t169 * t185 + t250 * t210 + t258) * MDP(16) + (g(1) * t178 - g(2) * t176 + t169 * t183 - t251 * t210 - t252) * MDP(17) + (t250 * MDP(16) - t251 * MDP(17)) * qJD(5) + ((-t226 * t264 - t274 * t285 - t262) * MDP(13) + (-t229 * t269 - t234 * t271 - t256) * MDP(14) + (MDP(16) * t217 + MDP(17) * t219) * g(3)) * t227;];
tau = t1;
