% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:40
% EndTime: 2019-12-31 18:24:44
% DurationCPUTime: 1.37s
% Computational Cost: add. (691->208), mult. (1569->299), div. (0->0), fcn. (850->6), ass. (0->112)
t187 = sin(pkin(8)) * pkin(1) + pkin(6);
t178 = t187 * qJD(1);
t201 = sin(qJ(3));
t203 = cos(qJ(3));
t251 = t203 * qJD(2) - t201 * t178;
t271 = -qJD(4) + t251;
t151 = -qJD(3) * pkin(3) - t271;
t247 = qJD(1) * t203;
t196 = t201 ^ 2;
t197 = t203 ^ 2;
t270 = (t196 - t197) * MDP(6);
t159 = t201 * qJD(2) + t203 * t178;
t153 = -qJD(3) * qJ(4) - t159;
t248 = qJD(1) * t201;
t233 = pkin(4) * t248;
t236 = t233 - t271;
t204 = -pkin(3) - pkin(7);
t269 = pkin(4) + t187;
t234 = qJD(2) * qJD(3);
t244 = qJD(3) * t201;
t254 = t178 * t244 - t203 * t234;
t137 = (qJD(4) - t233) * qJD(3) - t254;
t200 = sin(qJ(5));
t268 = t137 * t200;
t202 = cos(qJ(5));
t267 = t137 * t202;
t245 = qJD(3) * t200;
t174 = t202 * t247 + t245;
t235 = qJD(1) * qJD(3);
t226 = t201 * t235;
t180 = t200 * t226;
t147 = -t174 * qJD(5) + t180;
t266 = t147 * t202;
t240 = qJD(5) * t202;
t227 = t200 * t247;
t252 = qJD(5) * t227 + t202 * t226;
t148 = qJD(3) * t240 - t252;
t265 = t148 * t201;
t186 = qJD(5) + t248;
t264 = t174 * t186;
t263 = t174 * t203;
t243 = qJD(3) * t202;
t176 = -t227 + t243;
t262 = t176 * t186;
t261 = t186 * t200;
t260 = t186 * t201;
t259 = t186 * t202;
t258 = t186 * t204;
t205 = qJD(3) ^ 2;
t257 = t201 * t205;
t256 = t203 * t205;
t242 = qJD(3) * t203;
t255 = t147 * t201 + t176 * t242;
t156 = t178 * t242 + t201 * t234;
t229 = qJD(5) * t261;
t232 = t201 * t259;
t253 = qJD(3) * t232 + t203 * t229;
t188 = -cos(pkin(8)) * pkin(1) - pkin(2);
t212 = -qJ(4) * t201 + t188;
t165 = -pkin(3) * t203 + t212;
t154 = qJD(1) * t165;
t179 = qJD(1) * t188;
t249 = qJD(1) * t197;
t241 = qJD(4) * t201;
t191 = pkin(4) * t247;
t143 = -t153 + t191;
t239 = t143 * qJD(3);
t238 = t143 * qJD(5);
t237 = t203 * MDP(20);
t206 = qJD(1) ^ 2;
t231 = t201 * t203 * t206;
t230 = t202 * t249;
t228 = t203 * t240;
t167 = t269 * t203;
t225 = t203 * t235;
t224 = qJD(1) * t237;
t185 = pkin(3) * t226;
t219 = pkin(7) * t201 - qJ(4) * t203;
t208 = t219 * qJD(3) - t241;
t140 = t208 * qJD(1) + t185;
t141 = pkin(4) * t225 + t156;
t222 = -t140 * t200 + t202 * t141;
t221 = t186 * t228;
t220 = qJD(3) * t159 - t156;
t139 = t204 * qJD(3) + t236;
t157 = t204 * t203 + t212;
t144 = t157 * qJD(1);
t135 = t139 * t202 - t144 * t200;
t136 = t139 * t200 + t144 * t202;
t166 = t269 * t201;
t218 = t157 * t202 + t166 * t200;
t217 = t249 - t260;
t216 = -0.2e1 * qJD(3) * t154;
t215 = 0.2e1 * qJD(3) * t179;
t211 = t143 * t201 + t204 * t242;
t209 = -qJ(4) * t242 - t241;
t155 = t209 * qJD(1) + t185;
t192 = pkin(3) * t244;
t164 = t192 + t209;
t210 = qJD(1) * t164 + t187 * t205 + t155;
t146 = -qJD(3) * qJD(4) + t254;
t207 = -t146 * t203 + t156 * t201 + (t151 * t203 + t153 * t201) * qJD(3);
t193 = pkin(3) * t248;
t181 = t202 * t225;
t177 = -qJ(4) * t247 + t193;
t163 = qJD(3) * t167;
t162 = t269 * t244;
t160 = t219 * qJD(1) + t193;
t152 = t192 + t208;
t150 = t191 + t159;
t145 = t154 * t248;
t1 = [0.2e1 * t201 * MDP(5) * t225 - 0.2e1 * t235 * t270 + MDP(7) * t256 - MDP(8) * t257 + (-t187 * t256 + t201 * t215) * MDP(10) + (t187 * t257 + t203 * t215) * MDP(11) + t207 * MDP(12) + (t201 * t216 + t210 * t203) * MDP(13) + (-t210 * t201 + t203 * t216) * MDP(14) + (t154 * t164 + t155 * t165 + t207 * t187) * MDP(15) + (-t147 * t200 * t203 + (t200 * t244 - t228) * t176) * MDP(16) + ((-t174 * t200 + t176 * t202) * t244 + (-t266 + t148 * t200 + (t174 * t202 + t176 * t200) * qJD(5)) * t203) * MDP(17) + (-t217 * t245 - t221 + t255) * MDP(18) + (-t265 + (-t230 - t263) * qJD(3) + t253) * MDP(19) + (t186 + t248) * qJD(3) * t237 + ((-t152 * t200 + t163 * t202) * t186 - t162 * t174 + t167 * t148 + (-t202 * t239 + t222) * t201 + (-t136 * t201 - t218 * t186) * qJD(5) + (-t200 * t238 + t267 + ((-t157 * t200 + t166 * t202) * qJD(1) + t135) * qJD(3)) * t203) * MDP(21) + (t167 * t147 - t162 * t176 + (-(qJD(5) * t166 + t152) * t186 - (qJD(5) * t139 + t140) * t201) * t202 + (-(-qJD(5) * t157 + t163) * t186 + (qJD(5) * t144 - t141 + t239) * t201) * t200 + (-t202 * t238 - t268 + (-t218 * qJD(1) - t136) * qJD(3)) * t203) * MDP(22); (-t146 * t201 - t156 * t203) * MDP(15) + (t253 + t265) * MDP(21) + (t221 + t255) * MDP(22) + ((t151 * t201 - t153 * t203) * MDP(15) + (-t230 + t263) * MDP(21) + t217 * MDP(22) * t200) * qJD(3) + ((-MDP(11) + MDP(14)) * t203 + (-MDP(10) + MDP(13)) * t201) * t205; -MDP(5) * t231 + t206 * t270 + (-t179 * t248 + t220) * MDP(10) + (qJD(3) * t251 - t179 * t247 + t254) * MDP(11) + (-t177 * t247 + t145 - t220) * MDP(13) + ((0.2e1 * qJD(4) - t251) * qJD(3) + (t154 * t203 + t177 * t201) * qJD(1) - t254) * MDP(14) + (-pkin(3) * t156 - qJ(4) * t146 - t151 * t159 + t271 * t153 - t154 * t177) * MDP(15) + (-t176 * t261 + t266) * MDP(16) + ((-t148 - t262) * t202 + (-t147 + t264) * t200) * MDP(17) + (-t229 + t181 + (-t176 * t203 - t200 * t260) * qJD(1)) * MDP(18) + (-t186 * t240 + (-t232 + (t174 - t245) * t203) * qJD(1)) * MDP(19) - t186 * t224 + (qJ(4) * t148 + t268 - (t150 * t202 - t160 * t200) * t186 + t236 * t174 + (t143 * t202 - t200 * t258) * qJD(5) + (-t135 * t203 + t211 * t202) * qJD(1)) * MDP(21) + (qJ(4) * t147 + t267 + (t150 * t200 + t160 * t202) * t186 + t236 * t176 + (-t143 * t200 - t202 * t258) * qJD(5) + (t136 * t203 - t211 * t200) * qJD(1)) * MDP(22); MDP(13) * t231 + (-t196 * t206 - t205) * MDP(14) + (t145 + t156) * MDP(15) + t181 * MDP(21) + (t153 * MDP(15) - t174 * MDP(21) + (-t176 - t227) * MDP(22)) * qJD(3) + (-MDP(21) * t261 - MDP(22) * t259) * t186; t176 * t174 * MDP(16) + (-t174 ^ 2 + t176 ^ 2) * MDP(17) + (t180 + t264) * MDP(18) + (t252 + t262) * MDP(19) + qJD(3) * t224 + (t136 * t186 - t143 * t176 + t222) * MDP(21) + (t135 * t186 - t140 * t202 - t141 * t200 + t143 * t174) * MDP(22) + (-t174 * MDP(18) - MDP(19) * t243 - t136 * MDP(21) - t135 * MDP(22)) * qJD(5);];
tauc = t1;
