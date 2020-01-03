% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:07
% EndTime: 2019-12-31 17:21:11
% DurationCPUTime: 1.85s
% Computational Cost: add. (995->250), mult. (2491->360), div. (0->0), fcn. (1399->4), ass. (0->104)
t187 = cos(qJ(2));
t229 = qJD(1) * t187;
t175 = -qJD(3) + t229;
t218 = qJD(1) * qJD(2);
t255 = -0.2e1 * t218;
t185 = sin(qJ(2));
t254 = MDP(4) * t185;
t182 = t185 ^ 2;
t253 = MDP(5) * (-t187 ^ 2 + t182);
t170 = -pkin(2) * t187 - pkin(6) * t185 - pkin(1);
t158 = t170 * qJD(1);
t201 = pkin(2) * t185 - pkin(6) * t187;
t168 = t201 * qJD(2);
t159 = qJD(1) * t168;
t180 = pkin(5) * t229;
t173 = qJD(2) * pkin(6) + t180;
t186 = cos(qJ(3));
t208 = t185 * t218;
t224 = qJD(3) * t186;
t184 = sin(qJ(3));
t225 = qJD(3) * t184;
t251 = pkin(5) * t184;
t202 = -t158 * t225 + t186 * t159 - t173 * t224 + t208 * t251;
t132 = -pkin(3) * t208 - t202;
t141 = t158 * t184 + t173 * t186;
t138 = -qJ(4) * t175 + t141;
t252 = t138 * t175 + t132;
t250 = pkin(5) * t186;
t207 = t187 * t218;
t211 = t185 * t225;
t217 = qJD(2) * qJD(3);
t147 = qJD(1) * t211 + (-t207 - t217) * t186;
t206 = t184 * t217;
t209 = t185 * t224;
t226 = qJD(2) * t187;
t148 = t206 + (t184 * t226 + t209) * qJD(1);
t228 = qJD(2) * t184;
t230 = qJD(1) * t186;
t165 = t185 * t230 + t228;
t131 = pkin(3) * t148 + pkin(5) * t207 + qJ(4) * t147 - qJD(4) * t165;
t249 = t131 * t184;
t248 = t131 * t186;
t246 = t147 * t184;
t220 = t186 * qJD(2);
t231 = qJD(1) * t185;
t163 = t184 * t231 - t220;
t245 = t163 * t175;
t172 = -qJD(2) * pkin(2) + pkin(5) * t231;
t244 = t172 * t186;
t243 = t175 * t184;
t242 = t175 * t186;
t241 = t184 * t187;
t240 = t185 * t186;
t188 = qJD(2) ^ 2;
t239 = t185 * t188;
t238 = t186 * t187;
t237 = t187 * t188;
t189 = qJD(1) ^ 2;
t236 = t187 * t189;
t199 = pkin(3) * t184 - qJ(4) * t186;
t235 = qJD(4) * t184 + t175 * t199 + t180;
t234 = t184 * t168 + t170 * t224;
t233 = pkin(5) * t238 + t184 * t170;
t227 = qJD(2) * t185;
t223 = qJD(3) * t187;
t222 = qJD(4) * t175;
t221 = t185 * MDP(15);
t140 = t158 * t186 - t173 * t184;
t219 = qJD(4) - t140;
t216 = pkin(6) * t243;
t215 = pkin(6) * t242;
t214 = pkin(6) * t227;
t213 = pkin(6) * t220;
t212 = pkin(3) + t251;
t210 = t184 * t223;
t205 = pkin(1) * t255;
t204 = t163 + t220;
t203 = -t165 + t228;
t200 = pkin(3) * t186 + qJ(4) * t184;
t137 = pkin(3) * t175 + t219;
t198 = t137 * t186 - t138 * t184;
t197 = qJD(1) * t182 - t175 * t187;
t196 = pkin(5) + t199;
t195 = t158 * t224 + t184 * t159 - t173 * t225;
t194 = (qJ(4) - t250) * t231;
t193 = t168 * t186 - t170 * t225;
t191 = t184 * t227 - t186 * t223;
t190 = -t140 * t175 - t195;
t169 = -pkin(2) - t200;
t167 = t201 * qJD(1);
t156 = t184 * t167;
t152 = t196 * t185;
t146 = -t170 * t186 + t187 * t212;
t145 = -qJ(4) * t187 + t233;
t144 = pkin(3) * t165 + qJ(4) * t163;
t143 = -t167 * t186 - t212 * t231;
t142 = t156 + t194;
t139 = pkin(3) * t163 - qJ(4) * t165 + t172;
t136 = -t147 - t245;
t135 = (qJD(3) * t200 - qJD(4) * t186) * t185 + t196 * t226;
t134 = -pkin(3) * t227 - t191 * pkin(5) - t193;
t133 = qJ(4) * t227 - qJD(4) * t187 + (-t185 * t220 - t210) * pkin(5) + t234;
t130 = qJD(2) * t194 + t195 - t222;
t1 = [0.2e1 * t207 * t254 + t253 * t255 + MDP(6) * t237 - MDP(7) * t239 + (-pkin(5) * t237 + t185 * t205) * MDP(9) + (pkin(5) * t239 + t187 * t205) * MDP(10) + (-t147 * t240 + (t187 * t220 - t211) * t165) * MDP(11) + ((-t163 * t186 - t165 * t184) * t226 + (t246 - t148 * t186 + (t163 * t184 - t165 * t186) * qJD(3)) * t185) * MDP(12) + (t175 * t211 + t147 * t187 + (t165 * t185 + t186 * t197) * qJD(2)) * MDP(13) + (t175 * t209 + t148 * t187 + (-t163 * t185 - t184 * t197) * qJD(2)) * MDP(14) + (-t175 - t229) * qJD(2) * t221 + (-t193 * t175 - t202 * t187 + t172 * t209 + (t172 * t241 + (t170 * t230 + t140) * t185) * qJD(2) + (t185 * t148 + t163 * t226 - t175 * t191) * pkin(5)) * MDP(16) + ((-pkin(5) * t210 + t234) * t175 + t195 * t187 + (-pkin(5) * t147 - t172 * t225) * t185 + ((pkin(5) * t165 + t244) * t187 + (-pkin(5) * t242 - qJD(1) * t233 - t141) * t185) * qJD(2)) * MDP(17) + (t134 * t175 + t135 * t163 + t148 * t152 + (t139 * t228 + t132) * t187 + (t139 * t224 + t249 + (-qJD(1) * t146 - t137) * qJD(2)) * t185) * MDP(18) + (-t133 * t163 + t134 * t165 - t145 * t148 - t146 * t147 + t198 * t226 + (-t130 * t184 + t132 * t186 + (-t137 * t184 - t138 * t186) * qJD(3)) * t185) * MDP(19) + (-t133 * t175 - t135 * t165 + t147 * t152 + (-t139 * t220 - t130) * t187 + (t139 * t225 - t248 + (qJD(1) * t145 + t138) * qJD(2)) * t185) * MDP(20) + (t130 * t145 + t131 * t152 + t132 * t146 + t133 * t138 + t134 * t137 + t135 * t139) * MDP(21); -t236 * t254 + t189 * t253 + (-t165 * t242 - t246) * MDP(11) + ((-t147 + t245) * t186 + (t165 * t175 - t148) * t184) * MDP(12) + (-t175 * t224 + (t175 * t238 + t185 * t203) * qJD(1)) * MDP(13) + (t175 * t225 + (-t175 * t241 + t185 * t204) * qJD(1)) * MDP(14) + t175 * qJD(1) * t221 + (t167 * t242 - pkin(2) * t148 + (t172 * t184 + t215) * qJD(3) + (-t140 * t185 + (-t172 * t187 - t214) * t184 + (t185 * t243 - t187 * t204) * pkin(5)) * qJD(1)) * MDP(16) + (pkin(2) * t147 - t156 * t175 + (-t216 + t244) * qJD(3) + (-t172 * t238 + (t141 - t213) * t185 + (t175 * t240 + t187 * t203) * pkin(5)) * qJD(1)) * MDP(17) + (-t248 - t143 * t175 + t148 * t169 - t235 * t163 + (t139 * t184 + t215) * qJD(3) + (t137 * t185 + (-t139 * t187 - t214) * t184) * qJD(1)) * MDP(18) + (t142 * t163 - t143 * t165 + (t130 - t175 * t137 + (qJD(3) * t165 - t148) * pkin(6)) * t186 + ((qJD(3) * t163 - t147) * pkin(6) + t252) * t184) * MDP(19) + (-t249 + t142 * t175 + t147 * t169 + t235 * t165 + (-t139 * t186 + t216) * qJD(3) + (t139 * t238 + (-t138 + t213) * t185) * qJD(1)) * MDP(20) + (t131 * t169 - t137 * t143 - t138 * t142 - t235 * t139 + (qJD(3) * t198 + t130 * t186 + t132 * t184) * pkin(6)) * MDP(21) + (MDP(9) * t185 * t189 + MDP(10) * t236) * pkin(1); t136 * MDP(13) - MDP(14) * t206 + t190 * MDP(17) + (pkin(3) * t147 - qJ(4) * t148) * MDP(19) + (-t190 - 0.2e1 * t222) * MDP(20) + (-pkin(3) * t132 + qJ(4) * t130 - t137 * t141 + t138 * t219 - t139 * t144) * MDP(21) + (-t175 * MDP(14) - t172 * MDP(16) - t139 * MDP(18) + (t138 - t141) * MDP(19) + t144 * MDP(20) + MDP(12) * t165) * t165 + (t165 * MDP(11) + t172 * MDP(17) - t144 * MDP(18) + (t137 - t219) * MDP(19) - t139 * MDP(20) - MDP(12) * t163) * t163 + (-MDP(14) * t209 + (-MDP(14) * t241 + (MDP(15) + MDP(17) * t250 + 0.2e1 * pkin(3) * MDP(18) + (0.2e1 * qJ(4) - t250) * MDP(20)) * t185) * qJD(2)) * qJD(1) + (MDP(16) + MDP(18)) * (-t141 * t175 + t202); (t163 * t165 - t208) * MDP(18) + t136 * MDP(19) + (-t165 ^ 2 - t175 ^ 2) * MDP(20) + (t139 * t165 + t252) * MDP(21);];
tauc = t1;
