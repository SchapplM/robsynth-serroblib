% Calculate Coriolis joint torque vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:42
% EndTime: 2021-01-15 10:56:47
% DurationCPUTime: 1.34s
% Computational Cost: add. (909->201), mult. (2442->289), div. (0->0), fcn. (1655->6), ass. (0->101)
t216 = (qJD(1) * qJD(2));
t245 = -2 * t216;
t188 = sin(qJ(2));
t190 = cos(qJ(2));
t244 = (t188 ^ 2 - t190 ^ 2) * MDP(5);
t240 = -qJ(3) - pkin(5);
t211 = qJD(2) * t240;
t158 = qJD(3) * t190 + t188 * t211;
t152 = t158 * qJD(1);
t186 = sin(pkin(7));
t197 = t188 * qJD(3) - t190 * t211;
t238 = cos(pkin(7));
t193 = t238 * t197;
t126 = qJD(1) * t193 + t152 * t186;
t208 = t238 * t190;
t176 = qJD(1) * t208;
t220 = qJD(1) * t188;
t159 = t186 * t220 - t176;
t157 = qJD(4) + t159;
t179 = pkin(2) * t186 + pkin(6);
t215 = pkin(2) * t220;
t209 = t238 * t188;
t169 = t186 * t190 + t209;
t222 = qJD(1) * t169;
t243 = (pkin(3) * t222 + pkin(6) * t159 + qJD(4) * t179 + t215) * t157 + t126;
t195 = t186 * t197;
t210 = t238 * t152;
t127 = -qJD(1) * t195 + t210;
t182 = -pkin(2) * t190 - pkin(1);
t221 = qJD(1) * t182;
t173 = qJD(3) + t221;
t131 = pkin(3) * t159 - pkin(6) * t222 + t173;
t133 = t238 * t158 - t195;
t213 = t240 * t188;
t171 = qJD(1) * t213;
t239 = qJD(2) * pkin(2);
t167 = t171 + t239;
t174 = t240 * t190;
t172 = qJD(1) * t174;
t229 = t186 * t172;
t140 = t238 * t167 + t229;
t136 = -qJD(2) * pkin(3) - t140;
t198 = -t186 * t188 + t208;
t139 = -pkin(3) * t198 - pkin(6) * t169 + t182;
t164 = t198 * qJD(2);
t145 = -t238 * t174 + t186 * t213;
t161 = t169 * qJD(2);
t154 = qJD(1) * t161;
t201 = t126 * t169 - t145 * t154;
t242 = t136 * t164 - (qJD(4) * t139 + t133) * t157 + (qJD(4) * t131 + t127) * t198 + t201;
t241 = pkin(2) * t188;
t187 = sin(qJ(4));
t218 = qJD(4) * t187;
t212 = t188 * t216;
t175 = t186 * t212;
t155 = qJD(2) * t176 - t175;
t189 = cos(qJ(4));
t217 = t189 * qJD(2);
t224 = qJD(4) * t217 + t189 * t155;
t128 = -t218 * t222 + t224;
t237 = t128 * t187;
t236 = t139 * t154;
t230 = t222 * t187;
t146 = -t217 + t230;
t235 = t146 * t157;
t234 = t146 * t222;
t148 = qJD(2) * t187 + t189 * t222;
t233 = t148 * t157;
t232 = t148 * t222;
t231 = t155 * t187;
t228 = t187 * t154;
t191 = qJD(2) ^ 2;
t227 = t188 * t191;
t150 = t189 * t154;
t226 = t190 * t191;
t192 = qJD(1) ^ 2;
t225 = t190 * t192;
t165 = t238 * t172;
t141 = t186 * t167 - t165;
t219 = qJD(4) * t169;
t214 = t188 * t239;
t207 = pkin(1) * t245;
t206 = t189 * t157;
t205 = 0.2e1 * t222;
t137 = qJD(2) * pkin(6) + t141;
t123 = t131 * t189 - t137 * t187;
t124 = t131 * t187 + t137 * t189;
t200 = t150 + (-t159 * t187 - t218) * t157;
t199 = t164 * t189 - t169 * t218;
t143 = t238 * t171 + t229;
t194 = -t179 * t154 + (t136 + t143) * t157;
t180 = -t238 * pkin(2) - pkin(3);
t177 = pkin(2) * t212;
t144 = -t174 * t186 - t240 * t209;
t142 = t171 * t186 - t165;
t135 = pkin(3) * t161 - pkin(6) * t164 + t214;
t132 = t158 * t186 + t193;
t130 = pkin(3) * t154 - pkin(6) * t155 + t177;
t129 = t148 * qJD(4) + t231;
t125 = t189 * t130;
t1 = [0.2e1 * t190 * MDP(4) * t212 + t244 * t245 + MDP(6) * t226 - MDP(7) * t227 + (-pkin(5) * t226 + t188 * t207) * MDP(9) + (pkin(5) * t227 + t190 * t207) * MDP(10) + (t154 * t182 + t161 * t173 + (-t132 + (-qJD(1) * t198 + t159) * t241) * qJD(2)) * MDP(11) + (t155 * t182 + t164 * t173 + (t205 * t241 - t133) * qJD(2)) * MDP(12) + (t127 * t198 + t132 * t222 - t133 * t159 - t140 * t164 - t141 * t161 + t144 * t155 + t201) * MDP(13) + (t126 * t144 + t127 * t145 - t132 * t140 + t133 * t141 + (t173 + t221) * t214) * MDP(14) + (t128 * t169 * t189 + t199 * t148) * MDP(15) + ((-t146 * t189 - t148 * t187) * t164 + (-t237 - t129 * t189 + (t146 * t187 - t148 * t189) * qJD(4)) * t169) * MDP(16) + (-t128 * t198 + t148 * t161 + t169 * t150 + t199 * t157) * MDP(17) + (-t169 * t228 + t129 * t198 - t146 * t161 + (-t164 * t187 - t189 * t219) * t157) * MDP(18) + (-t154 * t198 + t157 * t161) * MDP(19) + (t123 * t161 - t125 * t198 + t144 * t129 + t132 * t146 + (t135 * t157 + t236 + (t136 * t169 + t137 * t198 - t145 * t157) * qJD(4)) * t189 + t242 * t187) * MDP(20) + (-t124 * t161 + t144 * t128 + t132 * t148 + (-(-qJD(4) * t145 + t135) * t157 - t236 + (-qJD(4) * t137 + t130) * t198 - t136 * t219) * t187 + t242 * t189) * MDP(21); -t188 * MDP(4) * t225 + t192 * t244 + (qJD(2) * t142 - t159 * t215 - t173 * t222 - t126) * MDP(11) + (-t210 + t143 * qJD(2) + t173 * t159 + (-t222 * t241 + t195) * qJD(1)) * MDP(12) + ((t141 - t142) * t222 + (-t140 + t143) * t159 + (-t154 * t186 - t238 * t155) * pkin(2)) * MDP(13) + (t140 * t142 - t141 * t143 + (-t238 * t126 + t127 * t186 - t173 * t220) * pkin(2)) * MDP(14) + (t148 * t206 + t237) * MDP(15) + ((t128 - t235) * t189 + (-t129 - t233) * t187) * MDP(16) + (t157 * t206 + t228 - t232) * MDP(17) + (t200 + t234) * MDP(18) - t157 * t222 * MDP(19) + (-t123 * t222 + t180 * t129 - t142 * t146 + t194 * t187 - t243 * t189) * MDP(20) + (t124 * t222 + t180 * t128 - t142 * t148 + t243 * t187 + t194 * t189) * MDP(21) + (t192 * t188 * MDP(9) + MDP(10) * t225) * pkin(1); t205 * qJD(2) * MDP(11) + (-t175 + (t176 - t159) * qJD(2)) * MDP(12) + (-t159 ^ 2 - t222 ^ 2) * MDP(13) + (t140 * t222 + t141 * t159 + t177) * MDP(14) + (t200 - t234) * MDP(20) + (-t157 ^ 2 * t189 - t228 - t232) * MDP(21); t148 * t146 * MDP(15) + (-t146 ^ 2 + t148 ^ 2) * MDP(16) + (t224 + t235) * MDP(17) + (-t231 + t233) * MDP(18) + t154 * MDP(19) + (t124 * t157 - t127 * t187 - t136 * t148 + t125) * MDP(20) + (t123 * t157 - t127 * t189 - t130 * t187 + t136 * t146) * MDP(21) + (-MDP(17) * t230 - t148 * MDP(18) - t124 * MDP(20) - t123 * MDP(21)) * qJD(4);];
tauc = t1;
