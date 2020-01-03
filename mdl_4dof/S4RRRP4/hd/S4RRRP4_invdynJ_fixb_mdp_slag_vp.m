% Calculate vector of inverse dynamics joint torques for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:45
% EndTime: 2019-12-31 17:15:47
% DurationCPUTime: 1.29s
% Computational Cost: add. (894->196), mult. (2104->261), div. (0->0), fcn. (1362->8), ass. (0->103)
t231 = cos(qJ(2));
t289 = pkin(5) + pkin(6);
t199 = t289 * t231;
t194 = qJD(1) * t199;
t228 = sin(qJ(3));
t178 = t228 * t194;
t229 = sin(qJ(2));
t198 = t289 * t229;
t192 = qJD(1) * t198;
t283 = qJD(2) * pkin(2);
t184 = -t192 + t283;
t288 = cos(qJ(3));
t257 = t184 * t288 - t178;
t187 = t228 * t231 + t229 * t288;
t176 = t187 * qJD(1);
t280 = t176 * qJ(4);
t293 = t280 - t257;
t227 = qJ(2) + qJ(3);
t219 = sin(t227);
t220 = cos(t227);
t230 = sin(qJ(1));
t232 = cos(qJ(1));
t252 = g(1) * t232 + g(2) * t230;
t292 = -g(3) * t220 + t252 * t219;
t224 = qJD(2) + qJD(3);
t221 = t231 * pkin(2);
t284 = pkin(1) + t221;
t274 = -t198 * t228 + t288 * t199;
t290 = t176 ^ 2;
t263 = t288 * t231;
t253 = qJD(1) * t263;
t271 = qJD(1) * t229;
t262 = t228 * t271;
t174 = -t253 + t262;
t197 = t284 * qJD(1);
t167 = pkin(3) * t174 + qJD(4) - t197;
t282 = t167 * t176;
t281 = t174 * qJ(4);
t279 = t228 * t229;
t278 = t230 * t220;
t277 = t232 * t220;
t153 = pkin(3) * t224 - t293;
t276 = t153 + t293;
t275 = -t192 * t288 - t178;
t273 = pkin(3) * t220 + t221;
t225 = t229 ^ 2;
t272 = -t231 ^ 2 + t225;
t270 = qJD(3) * t228;
t269 = qJD(1) * qJD(2);
t268 = qJDD(1) * t229;
t267 = qJDD(1) * t231;
t266 = t229 * t283;
t264 = qJD(2) * t289;
t182 = t288 * t194;
t261 = qJD(3) * t288;
t260 = t229 * t269;
t259 = t231 * t269;
t258 = qJDD(1) * t288;
t256 = t192 * t228 - t182;
t255 = -t198 * t288 - t199 * t228;
t254 = -t224 * t253 - t228 * t267 - t229 * t258;
t251 = g(1) * t230 - g(2) * t232;
t250 = t224 * t279;
t249 = t228 * t268 - t231 * t258;
t248 = -0.2e1 * pkin(1) * t269 - pkin(5) * qJDD(2);
t247 = -t184 * t228 - t182;
t171 = pkin(2) * t260 - qJDD(1) * t284;
t193 = t229 * t264;
t195 = t231 * t264;
t244 = -t193 * t288 - t195 * t228 - t198 * t261 - t199 * t270;
t173 = t174 ^ 2;
t222 = qJDD(2) + qJDD(3);
t243 = t176 * t174 * MDP(11) + (-t254 + (t174 - t262) * t224) * MDP(13) - t249 * MDP(14) + (-t173 + t290) * MDP(12) + t222 * MDP(15);
t233 = qJD(2) ^ 2;
t242 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t233 + t251;
t234 = qJD(1) ^ 2;
t241 = pkin(1) * t234 - pkin(5) * qJDD(1) + t252;
t165 = t224 * t187;
t157 = qJD(1) * t165 + t249;
t240 = t157 * pkin(3) + qJDD(4) + t171;
t166 = qJDD(2) * pkin(2) + t289 * (-t259 - t268);
t168 = t289 * (-t260 + t267);
t239 = qJD(3) * t247 + t166 * t288 - t228 * t168;
t238 = -qJD(3) * t274 + t228 * t193 - t195 * t288;
t237 = t228 * t166 + t168 * t288 + t184 * t261 - t194 * t270;
t236 = g(1) * t277 + g(2) * t278 + g(3) * t219 - t197 * t174 - t237;
t235 = t197 * t176 + t239 + t292;
t223 = -qJ(4) - t289;
t217 = pkin(2) * t288 + pkin(3);
t191 = pkin(1) + t273;
t186 = -t263 + t279;
t164 = -qJD(2) * t263 - t231 * t261 + t250;
t162 = -qJ(4) * t186 + t274;
t161 = -qJ(4) * t187 + t255;
t160 = -t280 + t275;
t159 = t256 + t281;
t156 = qJD(1) * t250 + t254;
t155 = -t247 - t281;
t150 = t164 * qJ(4) - t187 * qJD(4) + t238;
t149 = -qJ(4) * t165 - qJD(4) * t186 + t244;
t148 = -qJ(4) * t157 - qJD(4) * t174 + t237;
t147 = pkin(3) * t222 + qJ(4) * t156 - qJD(4) * t176 + t239;
t1 = [qJDD(1) * MDP(1) + t251 * MDP(2) + t252 * MDP(3) + (qJDD(1) * t225 + 0.2e1 * t229 * t259) * MDP(4) + 0.2e1 * (t229 * t267 - t269 * t272) * MDP(5) + (qJDD(2) * t229 + t231 * t233) * MDP(6) + (qJDD(2) * t231 - t229 * t233) * MDP(7) + (t229 * t248 + t231 * t242) * MDP(9) + (-t229 * t242 + t231 * t248) * MDP(10) + (-t156 * t187 - t164 * t176) * MDP(11) + (t156 * t186 - t157 * t187 + t164 * t174 - t165 * t176) * MDP(12) + (-t164 * t224 + t187 * t222) * MDP(13) + (-t165 * t224 - t186 * t222) * MDP(14) + (g(1) * t278 - g(2) * t277 - t157 * t284 - t197 * t165 + t171 * t186 + t174 * t266 + t222 * t255 + t224 * t238) * MDP(16) + (t156 * t284 + t197 * t164 + t171 * t187 + t176 * t266 - t219 * t251 - t222 * t274 - t224 * t244) * MDP(17) + (-t147 * t187 - t148 * t186 - t149 * t174 - t150 * t176 + t153 * t164 - t155 * t165 + t156 * t161 - t157 * t162 - t252) * MDP(18) + (t148 * t162 + t155 * t149 + t147 * t161 + t153 * t150 + t240 * (pkin(3) * t186 - t284) + t167 * (pkin(3) * t165 + t266) - g(1) * (-t191 * t230 - t223 * t232) - g(2) * (t191 * t232 - t223 * t230)) * MDP(19); MDP(6) * t268 + MDP(7) * t267 + qJDD(2) * MDP(8) + (-g(3) * t231 + t229 * t241) * MDP(9) + (g(3) * t229 + t231 * t241) * MDP(10) + (-t256 * t224 + (-t174 * t271 + t222 * t288 - t224 * t270) * pkin(2) + t235) * MDP(16) + (t275 * t224 + (-t176 * t271 - t222 * t228 - t224 * t261) * pkin(2) + t236) * MDP(17) + (t217 * t156 + (t155 + t159) * t176 + (-t153 + t160) * t174 + (-t157 * t228 + (-t174 * t288 + t176 * t228) * qJD(3)) * pkin(2)) * MDP(18) + (t147 * t217 - t155 * t160 - t153 * t159 - pkin(3) * t282 - g(3) * t273 - t252 * (-pkin(2) * t229 - pkin(3) * t219) + (-t167 * t271 + t148 * t228 + (-t153 * t228 + t155 * t288) * qJD(3)) * pkin(2)) * MDP(19) + t243 + (-MDP(4) * t229 * t231 + MDP(5) * t272) * t234; (-t224 * t247 + t235) * MDP(16) + (t224 * t257 + t236) * MDP(17) + (pkin(3) * t156 - t174 * t276) * MDP(18) + (t276 * t155 + (t147 - t282 + t292) * pkin(3)) * MDP(19) + t243; (-t173 - t290) * MDP(18) + (t153 * t176 + t155 * t174 + t240 - t251) * MDP(19);];
tau = t1;
