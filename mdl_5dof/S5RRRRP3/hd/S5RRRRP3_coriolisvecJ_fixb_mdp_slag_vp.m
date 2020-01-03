% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:28
% EndTime: 2019-12-31 21:49:31
% DurationCPUTime: 1.13s
% Computational Cost: add. (1397->202), mult. (2449->250), div. (0->0), fcn. (1241->6), ass. (0->107)
t197 = sin(qJ(4));
t195 = t197 ^ 2;
t200 = cos(qJ(4));
t196 = t200 ^ 2;
t276 = t195 + t196;
t198 = sin(qJ(3));
t199 = sin(qJ(2));
t262 = pkin(1) * qJD(1);
t227 = t199 * t262;
t185 = t198 * t227;
t201 = cos(qJ(3));
t202 = cos(qJ(2));
t226 = t202 * t262;
t173 = t201 * t226 - t185;
t238 = qJD(3) * t201;
t224 = pkin(2) * t238;
t272 = t173 - t224;
t261 = pkin(1) * qJD(2);
t223 = qJD(1) * t261;
t275 = qJD(3) * t227 + t199 * t223;
t246 = t197 * t200;
t274 = MDP(10) * t246 - (t195 - t196) * MDP(11);
t234 = qJD(5) * t197;
t235 = qJD(4) * t200;
t236 = qJD(4) * t197;
t174 = pkin(4) * t236 - qJ(5) * t235 - t234;
t239 = qJD(3) * t198;
t225 = pkin(2) * t239;
t167 = t174 + t225;
t244 = t199 * t201;
t211 = t198 * t202 + t244;
t172 = t211 * t262;
t273 = t167 - t172;
t190 = pkin(1) * t202 + pkin(2);
t245 = t198 * t199;
t271 = pkin(1) * t245 - t190 * t201;
t194 = qJD(1) + qJD(2);
t193 = qJD(3) + t194;
t269 = t276 * t193;
t178 = pkin(2) * t194 + t226;
t216 = t202 * t223;
t241 = t275 * t198;
t148 = (qJD(3) * t178 + t216) * t201 - t241;
t267 = 0.2e1 * qJD(4);
t266 = pkin(2) * t201;
t265 = pkin(3) * t193;
t203 = qJD(4) ^ 2;
t264 = pkin(8) * t203;
t263 = MDP(20) * pkin(8);
t260 = qJD(4) * pkin(4);
t259 = MDP(5) * t199;
t258 = MDP(9) * t201;
t154 = t190 * t239 + (t211 * qJD(2) + t199 * t238) * pkin(1);
t257 = t154 * t193;
t166 = t178 * t198 + t201 * t227;
t160 = pkin(8) * t193 + t166;
t256 = t160 * t197;
t255 = t160 * t200;
t254 = t166 * t193;
t171 = pkin(1) * t244 + t190 * t198 + pkin(8);
t253 = t171 * t203;
t179 = -pkin(4) * t200 - qJ(5) * t197 - pkin(3);
t252 = t179 * t193;
t188 = pkin(2) * t198 + pkin(8);
t251 = t188 * t203;
t247 = t193 * t200;
t144 = t197 * t148;
t140 = t160 * t235 + t144;
t165 = t178 * t201 - t185;
t243 = t165 * t236 + t166 * t247;
t242 = t172 * t247 + t173 * t236;
t237 = qJD(4) * qJ(5);
t233 = t171 * MDP(20);
t232 = t188 * MDP(20);
t231 = t203 * MDP(13);
t230 = -qJD(1) - t194;
t229 = -qJD(2) + t194;
t145 = t200 * t148;
t138 = t145 + (qJD(5) - t256) * qJD(4);
t218 = qJD(5) + t256;
t151 = t218 - t260;
t222 = t138 * t200 + t140 * t197 + t151 * t235;
t153 = t190 * t238 + (-t199 * t239 + (t201 * t202 - t245) * qJD(2)) * pkin(1);
t162 = t179 + t271;
t220 = t162 * t193 - t153;
t219 = (-pkin(3) + t271) * t193 - t153;
t149 = t178 * t239 + t198 * t216 + t275 * t201;
t214 = -t172 + t225;
t213 = pkin(4) * t197 - qJ(5) * t200;
t212 = t253 + t257;
t210 = -t178 * t238 + t241;
t141 = (t213 * qJD(4) - t234) * t193 + t149;
t209 = -t174 * t193 - t141 - t264;
t176 = t179 - t266;
t208 = t176 * t193 - t224;
t207 = (-pkin(3) - t266) * t193 - t224;
t143 = t154 + t174;
t206 = -t143 * t193 - t141 - t253;
t159 = -t165 - t265;
t205 = t203 * t200 * MDP(12) + (t149 * t197 + t159 * t235) * MDP(16) + t274 * t193 * t267;
t192 = t193 ^ 2;
t169 = t213 * t193;
t156 = t159 * t236;
t152 = t237 + t255;
t147 = -t165 + t252;
t142 = t147 * t236;
t1 = [(-t149 - t257) * MDP(8) + (-t153 * t193 + t210) * MDP(9) + t156 * MDP(15) + t142 * MDP(17) + (t269 * t153 + t222) * MDP(18) + (t141 * t162 + t143 * t147) * MDP(20) + (t230 * t259 + (t230 * MDP(6) - qJD(1) * t258) * t202) * t261 + (-t231 + t212 * MDP(16) + t206 * MDP(19) + (t140 * t171 + t151 * t153) * MDP(20) + (t219 * MDP(15) + t220 * MDP(17) + (-MDP(18) - t233) * t152) * qJD(4)) * t197 + ((-t149 - t212) * MDP(15) + t206 * MDP(17) + (t138 * t171 + t152 * t153) * MDP(20) + (t219 * MDP(16) + (-t147 - t220) * MDP(19) + t151 * t233) * qJD(4)) * t200 + t205; -t149 * MDP(8) + t210 * MDP(9) + (t156 + t242) * MDP(15) + (t142 + t242) * MDP(17) + t222 * MDP(18) + (t141 * t176 + t273 * t147) * MDP(20) + (t229 * t259 + (t229 * MDP(6) - qJD(2) * t258) * t202) * t262 + (-t214 * MDP(8) + (-t276 * MDP(18) + MDP(9)) * t272) * t193 + (-t141 * MDP(19) + (t140 * t188 - t272 * t151) * MDP(20) + (-MDP(13) + (MDP(16) - MDP(19)) * t188) * t203 + (t214 * MDP(16) - t273 * MDP(19)) * t193 + (t207 * MDP(15) + t208 * MDP(17) + (-MDP(18) - t232) * t152) * qJD(4)) * t197 + ((-t193 * t225 - t149 - t251) * MDP(15) + (-t167 * t193 - t141 - t251) * MDP(17) + (t138 * t188 - t272 * t152) * MDP(20) + ((t173 + t207) * MDP(16) + (-t147 - t173 - t208) * MDP(19) + t151 * t232) * qJD(4)) * t200 + t205; (-t149 + t254) * MDP(8) + (t165 * t193 - t148) * MDP(9) + (t156 + t243) * MDP(15) + (t142 + t243) * MDP(17) + (-t269 * t165 + t222) * MDP(18) + (t141 * t179 + (-t166 + t174) * t147) * MDP(20) + (-t231 + (-t254 + t264) * MDP(16) + (t209 + t254) * MDP(19) + (pkin(8) * t140 - t151 * t165) * MDP(20) + ((-MDP(15) * pkin(3) + MDP(17) * t179) * t193 + (-MDP(18) - t263) * t152) * qJD(4)) * t197 + ((-t149 - t264) * MDP(15) + t209 * MDP(17) + (pkin(8) * t138 - t152 * t165) * MDP(20) + ((t165 - t265) * MDP(16) + (-t147 - t165 - t252) * MDP(19) + t151 * t263) * qJD(4)) * t200 + t205; -t145 * MDP(16) + (qJD(5) * t267 + t145) * MDP(19) + (-pkin(4) * t140 + qJ(5) * t138 - t147 * t169 - t151 * t255 + t218 * t152) * MDP(20) + (-MDP(15) - MDP(17)) * t144 - t274 * t192 + ((-t159 * MDP(15) - t147 * MDP(17) + (t152 - t237) * MDP(18) + t169 * MDP(19)) * t197 + (-t159 * MDP(16) + t169 * MDP(17) + (qJD(5) - t151 - t260) * MDP(18) + t147 * MDP(19)) * t200) * t193; -t192 * MDP(17) * t246 + (-t192 * t195 - t203) * MDP(19) + (t147 * t193 * t197 - qJD(4) * t152 + t140) * MDP(20);];
tauc = t1;
