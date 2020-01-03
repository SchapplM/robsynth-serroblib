% Calculate vector of inverse dynamics joint torques for
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:21
% EndTime: 2019-12-31 17:55:23
% DurationCPUTime: 1.50s
% Computational Cost: add. (1134->231), mult. (2050->259), div. (0->0), fcn. (1300->8), ass. (0->107)
t228 = sin(qJ(1));
t230 = cos(qJ(1));
t297 = g(1) * t228 - g(2) * t230;
t223 = sin(pkin(7));
t224 = cos(pkin(7));
t227 = sin(qJ(4));
t229 = cos(qJ(4));
t189 = t223 * t229 + t224 * t227;
t294 = t189 * qJD(1);
t301 = t294 * qJD(4);
t300 = MDP(16) + MDP(18);
t272 = qJD(1) * t223;
t261 = t227 * t272;
t289 = qJD(4) * t261 - t189 * qJDD(1);
t226 = -pkin(1) - qJ(3);
t197 = t226 * qJD(1) + qJD(2);
t274 = t223 ^ 2 + t224 ^ 2;
t299 = t197 * t274;
t220 = qJDD(1) * qJ(2);
t221 = qJD(1) * qJD(2);
t296 = t220 + t221;
t194 = qJDD(3) + t296;
t252 = g(1) * t230 + g(2) * t228;
t298 = t194 - t252;
t246 = t223 * MDP(7) + t224 * MDP(8);
t292 = -MDP(17) + MDP(20);
t219 = pkin(7) + qJ(4);
t209 = sin(t219);
t210 = cos(t219);
t288 = -qJD(1) * qJD(3) + qJDD(1) * t226;
t190 = qJDD(2) + t288;
t255 = -pkin(6) * qJDD(1) + t190;
t174 = t255 * t223;
t175 = t255 * t224;
t258 = -pkin(6) * qJD(1) + t197;
t177 = t258 * t224;
t269 = qJD(4) * t229;
t264 = t229 * t174 + t227 * t175 + t177 * t269;
t291 = -g(3) * t210 - t297 * t209 + t264;
t285 = -pkin(6) + t226;
t191 = t285 * t223;
t192 = t285 * t224;
t165 = t191 * t227 - t192 * t229;
t155 = -qJD(3) * t189 - qJD(4) * t165;
t166 = t191 * t229 + t192 * t227;
t290 = -qJD(4) * t155 - qJDD(4) * t166 - t210 * t252;
t278 = t229 * t224;
t263 = qJD(1) * t278;
t182 = -t261 + t263;
t287 = t182 ^ 2;
t286 = 0.2e1 * t221;
t211 = t223 * pkin(3);
t284 = pkin(1) * qJDD(1);
t283 = qJDD(4) * pkin(4);
t176 = t258 * t223;
t282 = t176 * t227;
t281 = t294 * t182;
t203 = qJ(2) + t211;
t270 = qJD(4) * t227;
t184 = -t223 * t269 - t224 * t270;
t188 = t223 * t227 - t278;
t277 = t184 * qJD(4) - t188 * qJDD(4);
t276 = t230 * pkin(1) + t228 * qJ(2);
t160 = t176 * t229 + t177 * t227;
t271 = qJD(4) * t160;
t208 = qJD(1) * qJ(2) + qJD(3);
t159 = t177 * t229 - t282;
t268 = qJD(5) - t159;
t266 = qJDD(1) * t223;
t265 = qJDD(4) * qJ(5);
t262 = t224 * t269;
t193 = pkin(3) * t272 + t208;
t260 = g(2) * t276;
t259 = qJDD(2) - t297;
t257 = t274 * MDP(9);
t256 = t274 * t190;
t254 = t159 + t282;
t253 = t227 * t174 - t229 * t175 + t176 * t269 + t177 * t270;
t187 = pkin(3) * t266 + t194;
t248 = -pkin(4) * t209 + qJ(5) * t210;
t247 = -qJDD(1) * t278 + t227 * t266;
t163 = t247 + t301;
t245 = t163 * t188 + t182 * t184;
t185 = -t223 * t270 + t262;
t242 = -qJD(4) * t185 - qJDD(4) * t189;
t239 = -t248 + t211;
t164 = qJD(1) * t262 - t289;
t237 = pkin(4) * t164 + qJ(5) * t163 + t187;
t236 = g(3) * t209 - t297 * t210 - t253;
t156 = -qJD(3) * t188 + qJD(4) * t166;
t234 = -qJD(4) * t156 - qJDD(4) * t165 - t209 * t252;
t148 = t265 + (qJD(5) - t282) * qJD(4) + t264;
t150 = qJDD(5) + t253 - t283;
t154 = -qJD(4) * pkin(4) + t268;
t157 = qJD(4) * qJ(5) + t160;
t233 = t148 * t189 + t150 * t188 - t154 * t184 + t157 * t185 - t297;
t158 = pkin(4) * t294 - qJ(5) * t182 + t193;
t232 = t158 * t182 + qJDD(5) - t236;
t231 = qJD(1) ^ 2;
t225 = -pkin(6) - qJ(3);
t213 = t230 * qJ(2);
t179 = t294 ^ 2;
t162 = pkin(4) * t182 + qJ(5) * t294;
t161 = pkin(4) * t189 + qJ(5) * t188 + t203;
t151 = pkin(4) * t185 - qJ(5) * t184 + qJD(5) * t188 + qJD(2);
t149 = -qJD(5) * t182 + t237;
t1 = [qJDD(1) * MDP(1) + t297 * MDP(2) + t252 * MDP(3) + (t259 - 0.2e1 * t284) * MDP(4) + (0.2e1 * t220 + t286 - t252) * MDP(5) + (-(qJDD(2) - t284) * pkin(1) - g(1) * (-pkin(1) * t228 + t213) - t260 + (t220 + t286) * qJ(2)) * MDP(6) + (t297 + t274 * (-t190 - t288)) * MDP(9) + (t194 * qJ(2) + t208 * qJD(2) - g(1) * (t226 * t228 + t213) - g(2) * (qJ(3) * t230 + t276) + t226 * t256 - qJD(3) * t299) * MDP(10) + t245 * MDP(11) + (t163 * t189 + t164 * t188 - t182 * t185 - t184 * t294) * MDP(12) + t277 * MDP(13) + t242 * MDP(14) + (qJD(2) * t294 + t164 * t203 + t185 * t193 + t187 * t189 + t234) * MDP(16) + (qJD(2) * t182 - t163 * t203 + t184 * t193 - t187 * t188 + t290) * MDP(17) + (t149 * t189 + t151 * t294 + t158 * t185 + t161 * t164 + t234) * MDP(18) + (-t155 * t294 + t156 * t182 - t163 * t165 - t164 * t166 - t233) * MDP(19) + (t149 * t188 - t151 * t182 - t158 * t184 + t161 * t163 - t290) * MDP(20) + (t148 * t166 + t157 * t155 + t149 * t161 + t158 * t151 + t150 * t165 + t154 * t156 - g(1) * t213 - t260 + (-g(1) * t239 + g(2) * t225) * t230 + (-g(1) * (-pkin(1) + t225) - g(2) * t239) * t228) * MDP(21) + t246 * (t296 + t298); t259 * MDP(6) + (-qJD(1) * t208 + t256 - t297) * MDP(10) + (-t164 * t189 - t185 * t294 - t245) * MDP(19) + (-qJD(1) * t158 + t233) * MDP(21) + (-MDP(6) * qJ(2) - MDP(5) - t246) * t231 + (-pkin(1) * MDP(6) + MDP(4) - t257) * qJDD(1) + t292 * (qJD(1) * t182 - t242) + t300 * (-qJD(1) * t294 + t277); (qJD(1) * t299 + t298) * MDP(10) + (-t179 - t287) * MDP(19) + (t157 * t294 + (-qJD(5) - t154) * t182 + t237 - t252) * MDP(21) - t231 * t257 + t292 * (t247 + 0.2e1 * t301) + t246 * qJDD(1) + (qJD(4) * (t182 + t263) - t289) * t300; MDP(11) * t281 + (-t179 + t287) * MDP(12) - t247 * MDP(13) + ((t182 - t263) * qJD(4) + t289) * MDP(14) + qJDD(4) * MDP(15) + (-t182 * t193 + t236 + t271) * MDP(16) + (t254 * qJD(4) + t193 * t294 - t291) * MDP(17) + (-t162 * t294 - t232 + t271 + 0.2e1 * t283) * MDP(18) + (pkin(4) * t163 - qJ(5) * t164 + (t157 - t160) * t182 + (t154 - t268) * t294) * MDP(19) + (0.2e1 * t265 - t158 * t294 + t162 * t182 + (0.2e1 * qJD(5) - t254) * qJD(4) + t291) * MDP(20) + (-t150 * pkin(4) - g(3) * t248 + t148 * qJ(5) - t154 * t160 + t268 * t157 - t158 * t162 - t297 * (pkin(4) * t210 + qJ(5) * t209)) * MDP(21); (-qJDD(4) + t281) * MDP(18) - t247 * MDP(19) + (-qJD(4) ^ 2 - t287) * MDP(20) + (-qJD(4) * t157 + t232 - t283) * MDP(21);];
tau = t1;
