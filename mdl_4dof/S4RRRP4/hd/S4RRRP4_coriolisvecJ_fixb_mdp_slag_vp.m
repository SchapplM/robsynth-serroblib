% Calculate Coriolis joint torque vector for
% S4RRRP4
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
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:26
% EndTime: 2021-01-15 14:30:29
% DurationCPUTime: 1.08s
% Computational Cost: add. (888->173), mult. (2349->236), div. (0->0), fcn. (1496->4), ass. (0->100)
t193 = cos(qJ(2));
t248 = pkin(5) + pkin(6);
t176 = t248 * t193;
t172 = qJD(1) * t176;
t190 = sin(qJ(3));
t158 = t190 * t172;
t191 = sin(qJ(2));
t175 = t248 * t191;
t170 = qJD(1) * t175;
t246 = qJD(2) * pkin(2);
t164 = -t170 + t246;
t192 = cos(qJ(3));
t215 = t164 * t192 - t158;
t237 = t191 * t192;
t169 = t190 * t193 + t237;
t228 = qJD(1) * t169;
t241 = t228 * qJ(4);
t133 = -t241 + t215;
t187 = qJD(2) + qJD(3);
t209 = t187 * MDP(18);
t145 = t187 * t169;
t143 = t145 * qJD(1);
t220 = qJD(2) * t248;
t208 = qJD(1) * t220;
t165 = t191 * t208;
t166 = t193 * t208;
t226 = qJD(3) * t190;
t202 = -(qJD(3) * t164 - t165) * t192 + t190 * t166 + t172 * t226;
t253 = -qJ(4) * t143 - t202;
t223 = qJD(1) * qJD(2);
t252 = -0.2e1 * t223;
t251 = t191 * MDP(4);
t250 = (t191 ^ 2 - t193 ^ 2) * MDP(5);
t249 = t228 ^ 2;
t247 = pkin(2) * t187;
t238 = t190 * t191;
t207 = t187 * t238;
t218 = t193 * t223;
t235 = t192 * t193;
t219 = qJD(1) * t235;
t230 = qJD(3) * t219 + t192 * t218;
t142 = qJD(1) * t207 - t230;
t245 = qJ(4) * t142;
t227 = qJD(1) * t191;
t154 = t190 * t227 - t219;
t243 = qJ(4) * t154;
t242 = t154 * t187;
t186 = -pkin(2) * t193 - pkin(1);
t174 = t186 * qJD(1);
t239 = t174 * t228;
t194 = qJD(2) ^ 2;
t236 = t191 * t194;
t162 = t192 * t172;
t234 = t193 * t194;
t195 = qJD(1) ^ 2;
t233 = t193 * t195;
t132 = pkin(3) * t187 + t133;
t232 = t132 - t133;
t231 = -t170 * t192 - t158;
t225 = qJD(3) * t192;
t217 = pkin(3) * t154 + qJD(4);
t146 = t174 + t217;
t224 = qJD(4) + t146;
t222 = pkin(2) * t227;
t221 = t191 * t246;
t138 = pkin(3) * t143 + qJD(1) * t221;
t216 = pkin(1) * t252;
t214 = t190 * t165 - t166 * t192;
t212 = t170 * t190 - t162;
t210 = t187 * t191;
t206 = -t164 * t190 - t162;
t205 = t175 * t190 - t176 * t192;
t153 = t154 ^ 2;
t204 = t228 * t154 * MDP(11) + (-qJD(1) * t190 * t210 + t230 + t242) * MDP(13) + (-t153 + t249) * MDP(12);
t171 = t191 * t220;
t173 = t193 * t220;
t203 = -t171 * t192 - t173 * t190 - t175 * t225 - t176 * t226;
t201 = qJD(3) * t206 + t214;
t200 = qJD(3) * t205 + t171 * t190 - t173 * t192;
t199 = t174 * t154 + t202;
t198 = t201 + t245;
t197 = (-t162 + (-t164 - t247) * t190) * qJD(3) + t214;
t196 = t154 * t224 - t253;
t185 = pkin(2) * t192 + pkin(3);
t177 = t225 * t247;
t168 = -t235 + t238;
t148 = pkin(3) * t168 + t186;
t147 = pkin(3) * t228 + t222;
t144 = -qJD(2) * t235 - t193 * t225 + t207;
t141 = pkin(3) * t145 + t221;
t140 = -qJ(4) * t168 - t205;
t139 = -qJ(4) * t169 - t175 * t192 - t176 * t190;
t137 = -t241 + t231;
t136 = t212 + t243;
t134 = -t206 - t243;
t129 = qJ(4) * t144 - qJD(4) * t169 + t200;
t128 = -qJ(4) * t145 - qJD(4) * t168 + t203;
t127 = -qJD(4) * t228 + t198;
t126 = -qJD(4) * t154 + t253;
t1 = [0.2e1 * t218 * t251 + t250 * t252 + MDP(6) * t234 - MDP(7) * t236 + (-pkin(5) * t234 + t191 * t216) * MDP(9) + (pkin(5) * t236 + t193 * t216) * MDP(10) + (-t142 * t169 - t144 * t228) * MDP(11) + (t142 * t168 - t143 * t169 + t144 * t154 - t145 * t228) * MDP(12) + (t186 * t143 + t174 * t145 + (qJD(1) * t168 + t154) * t221) * MDP(16) + (-t186 * t142 - t174 * t144 + 0.2e1 * t221 * t228) * MDP(17) + (t138 * t168 + t141 * t154 + t143 * t148 + t145 * t146) * MDP(18) + (t138 * t169 + t141 * t228 - t142 * t148 - t144 * t146) * MDP(19) + (-t126 * t168 - t127 * t169 - t128 * t154 - t129 * t228 + t132 * t144 - t134 * t145 + t139 * t142 - t140 * t143) * MDP(20) + (t126 * t140 + t127 * t139 + t128 * t134 + t129 * t132 + t138 * t148 + t141 * t146) * MDP(21) + (-t144 * MDP(13) - t145 * MDP(14) + MDP(16) * t200 - MDP(17) * t203 + t129 * MDP(18) - t128 * MDP(19)) * t187; -t233 * t251 + t195 * t250 + (-t154 * t222 - t187 * t212 + t197 - t239) * MDP(16) + (t187 * t231 - t222 * t228 - t177 + t199) * MDP(17) + (-t136 * t187 - t147 * t154 - t224 * t228 + t197 + t245) * MDP(18) + (t137 * t187 - t147 * t228 - t177 + t196) * MDP(19) + (t142 * t185 + (t134 + t136) * t228 + (-t132 + t137) * t154 + (-t143 * t190 + (-t154 * t192 + t190 * t228) * qJD(3)) * pkin(2)) * MDP(20) + (t127 * t185 - t132 * t136 - t134 * t137 - t146 * t147 + (t126 * t190 + (-t132 * t190 + t134 * t192) * qJD(3)) * pkin(2)) * MDP(21) + t204 + (MDP(9) * t191 * t195 + MDP(10) * t233) * pkin(1); (-t187 * t206 + t201 - t239) * MDP(16) + (t187 * t215 + t199) * MDP(17) + (t134 * t187 + (-t146 - t217) * t228 + t198) * MDP(18) + (-pkin(3) * t249 + t133 * t187 + t196) * MDP(19) + (pkin(3) * t142 - t154 * t232) * MDP(20) + (t232 * t134 + (-t146 * t228 + t127) * pkin(3)) * MDP(21) + t204; t228 * t209 + (t230 - t242) * MDP(19) + (-t153 - t249) * MDP(20) + (t132 * t228 + t134 * t154 + t138) * MDP(21) + (t209 * t237 + (-MDP(19) * t210 + t193 * t209) * t190) * qJD(1);];
tauc = t1;
