% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:55
% EndTime: 2019-12-31 18:12:59
% DurationCPUTime: 1.18s
% Computational Cost: add. (1031->215), mult. (2684->254), div. (0->0), fcn. (1817->4), ass. (0->92)
t209 = sin(pkin(7));
t210 = cos(pkin(7));
t259 = (qJ(2) * MDP(7) + MDP(6)) * (t209 ^ 2 + t210 ^ 2);
t251 = pkin(6) + qJ(2);
t193 = t251 * t209;
t190 = qJD(1) * t193;
t194 = t251 * t210;
t191 = qJD(1) * t194;
t212 = sin(qJ(3));
t213 = cos(qJ(3));
t244 = -t213 * t190 - t212 * t191;
t237 = -qJD(4) + t244;
t189 = t209 * t213 + t210 * t212;
t185 = t189 * qJD(3);
t166 = qJD(1) * t185;
t240 = qJD(1) * t213;
t228 = t210 * t240;
t241 = qJD(1) * t212;
t230 = t209 * t241;
t180 = -t228 + t230;
t257 = t166 * qJ(5) + t180 * qJD(5);
t182 = t189 * qJD(1);
t233 = MDP(13) - MDP(16);
t255 = -MDP(14) + MDP(20);
t238 = qJD(3) * t213;
t245 = t213 * t210;
t148 = (qJD(2) * t209 + qJD(3) * t194) * t212 - qJD(2) * t245 + t193 * t238;
t254 = t180 ^ 2;
t175 = t182 ^ 2;
t253 = pkin(4) * t180;
t252 = pkin(3) + qJ(5);
t249 = qJ(4) * t166;
t248 = qJ(4) * t180;
t247 = t180 * t182;
t246 = t209 * t212;
t160 = -t212 * t190 + t213 * t191;
t239 = qJD(3) * t212;
t145 = -pkin(4) * t182 + t244;
t236 = qJD(4) - t145;
t146 = t160 - t253;
t235 = -qJD(5) - t146;
t234 = qJD(1) * qJD(2);
t232 = MDP(14) - MDP(17);
t231 = (-MDP(16) + MDP(21));
t204 = -pkin(2) * t210 - pkin(1);
t229 = t209 * t239;
t227 = t212 * t234;
t226 = t213 * t234;
t161 = t193 * t213 + t212 * t194;
t225 = t190 * t238 + t191 * t239 + t209 * t227 - t210 * t226;
t142 = -t190 * t239 + t191 * t238 + t209 * t226 + t210 * t227;
t156 = -qJD(3) * qJ(4) - t160;
t199 = qJD(3) * t228;
t165 = qJD(1) * t229 - t199;
t224 = pkin(3) * t166 + qJ(4) * t165;
t223 = -t193 * t212 + t194 * t213;
t184 = -t210 * t238 + t229;
t222 = qJ(4) * t184 - qJD(4) * t189;
t220 = -qJ(4) * t189 + t204;
t219 = pkin(4) * t165 - t142;
t218 = -pkin(4) * t166 - t225;
t192 = t204 * qJD(1) + qJD(2);
t208 = qJD(3) * qJD(4);
t132 = t208 + t218;
t137 = -qJD(4) * t182 + t224;
t217 = -t209 * t240 - t210 * t241 + t182;
t216 = -qJ(4) * t182 + t192;
t149 = t189 * qJD(2) + t223 * qJD(3);
t214 = qJD(3) ^ 2;
t205 = 0.2e1 * t208;
t188 = -t245 + t246;
t168 = qJD(3) * t180;
t158 = pkin(3) * t188 + t220;
t157 = pkin(3) * t182 + t248;
t154 = t199 + (t180 - t230) * qJD(3);
t152 = -qJD(3) * pkin(3) - t237;
t151 = -pkin(4) * t188 + t223;
t150 = pkin(4) * t189 + t161;
t147 = pkin(3) * t180 + t216;
t144 = t252 * t188 + t220;
t143 = pkin(3) * t185 + t222;
t141 = t252 * t182 + t248;
t140 = qJD(5) - t156 - t253;
t139 = -t208 + t225;
t138 = -t252 * qJD(3) + t236;
t136 = -pkin(4) * t184 + t149;
t135 = -pkin(4) * t185 - t148;
t134 = t252 * t180 + t216;
t133 = -qJD(3) * qJD(5) - t219;
t131 = qJD(5) * t188 + t252 * t185 + t222;
t130 = t137 + t257;
t1 = [(-t165 * t189 - t182 * t184) * MDP(8) + (t165 * t188 - t166 * t189 + t180 * t184 - t182 * t185) * MDP(9) + (t166 * t204 + t185 * t192) * MDP(13) + (-t165 * t204 - t184 * t192) * MDP(14) + (t139 * t188 + t142 * t189 + t148 * t180 + t149 * t182 - t152 * t184 + t156 * t185 - t161 * t165 - t166 * t223) * MDP(15) + (-t137 * t188 - t143 * t180 - t147 * t185 - t158 * t166) * MDP(16) + (-t137 * t189 - t143 * t182 + t147 * t184 + t158 * t165) * MDP(17) + (t137 * t158 - t139 * t223 + t142 * t161 + t143 * t147 + t148 * t156 + t149 * t152) * MDP(18) + (-t132 * t188 + t133 * t189 - t135 * t180 + t136 * t182 - t138 * t184 - t140 * t185 - t150 * t165 - t151 * t166) * MDP(19) + (-t130 * t189 - t131 * t182 + t134 * t184 + t144 * t165) * MDP(20) + (t130 * t188 + t131 * t180 + t134 * t185 + t144 * t166) * MDP(21) + (t130 * t144 + t131 * t134 + t132 * t151 + t133 * t150 + t135 * t140 + t136 * t138) * MDP(22) + (-t184 * MDP(10) - t185 * MDP(11) + t135 * MDP(20) - t136 * MDP(21) + t232 * t148 - t233 * t149) * qJD(3) + 0.2e1 * t234 * t259; (t168 - t199) * MDP(17) + (-t156 * t180 + t224) * MDP(18) + (t140 * t180 + t224 + t257) * MDP(22) + ((-qJD(4) - t152) * MDP(18) + (-qJD(4) - t138) * MDP(22)) * t182 - t255 * t199 + (t255 * (t180 + t230) + (t189 * MDP(13) + MDP(17) * t246) * qJD(1) + ((2 * t231) + MDP(13)) * t182) * qJD(3) - qJD(1) ^ 2 * t259 + (MDP(15) + MDP(19)) * (-t254 - t175); -t254 * MDP(9) + t154 * MDP(10) + t225 * MDP(14) + (pkin(3) * t165 - t249) * MDP(15) + (t205 - t225) * MDP(17) + (-qJ(4) * t139 - t147 * t157 - t152 * t160 + t237 * t156) * MDP(18) + (t165 * t252 - t249) * MDP(19) + (t205 + t218) * MDP(20) + t219 * MDP(21) + (qJ(4) * t132 - t133 * t252 - t134 * t141 + t235 * t138 + t236 * t140) * MDP(22) + (-t192 * MDP(13) + (-t156 - t160) * MDP(15) + t147 * MDP(16) + t157 * MDP(17) + (t140 + t235) * MDP(19) + t141 * MDP(20) - t134 * MDP(21) + MDP(9) * t182) * t182 + (t217 * MDP(11) - t145 * MDP(20) + (0.2e1 * qJD(5) + t146) * MDP(21) + t233 * t160 + t232 * t244) * qJD(3) + (t182 * MDP(8) + t192 * MDP(14) + (t152 + t237) * MDP(15) + t157 * MDP(16) - t147 * MDP(17) + (t138 - t236) * MDP(19) - t134 * MDP(20) - t141 * MDP(21)) * t180 + (-pkin(3) * MDP(18) - t233) * t142; (-t165 + t168) * MDP(15) + (qJD(3) * t156 + t147 * t182 + t142) * MDP(18) + t154 * MDP(19) + (t134 * t182 + (-qJD(5) - t140) * qJD(3) - t219) * MDP(22) + t231 * t247 + (MDP(17) + MDP(20)) * (-t214 - t175); MDP(20) * t247 + (-t254 - t214) * MDP(21) + (-t134 * t180 + t132) * MDP(22) + (t217 * MDP(19) + t138 * MDP(22)) * qJD(3);];
tauc = t1;
