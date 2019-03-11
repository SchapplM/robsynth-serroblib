% Calculate joint inertia matrix for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR10_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:28
% EndTime: 2019-03-09 11:09:31
% DurationCPUTime: 1.08s
% Computational Cost: add. (1581->231), mult. (3499->321), div. (0->0), fcn. (3887->10), ass. (0->106)
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t173 = cos(pkin(6));
t171 = sin(pkin(6));
t176 = sin(qJ(2));
t221 = t171 * t176;
t150 = t170 * t221 - t172 * t173;
t151 = t170 * t173 + t172 * t221;
t175 = sin(qJ(4));
t228 = cos(qJ(4));
t137 = -t175 * t150 + t228 * t151;
t235 = 0.2e1 * t137;
t178 = cos(qJ(2));
t220 = t171 * t178;
t234 = 2 * MDP(13);
t233 = 2 * MDP(22);
t232 = 2 * MDP(24);
t231 = 2 * MDP(31);
t230 = 2 * MDP(32);
t229 = pkin(4) + pkin(10);
t227 = pkin(1) * t176;
t226 = pkin(1) * t178;
t225 = pkin(9) + qJ(3);
t224 = MDP(25) * pkin(4);
t223 = pkin(2) * MDP(14);
t219 = t173 * MDP(8);
t218 = t176 * MDP(6);
t199 = pkin(8) * t220;
t146 = t199 + (qJ(3) + t227) * t173;
t147 = (-pkin(2) * t178 - qJ(3) * t176 - pkin(1)) * t171;
t131 = -t146 * t170 + t172 * t147;
t124 = -pkin(3) * t220 - pkin(9) * t151 + t131;
t132 = t172 * t146 + t170 * t147;
t127 = -pkin(9) * t150 + t132;
t119 = t175 * t124 + t228 * t127;
t216 = MDP(11) * t172;
t215 = MDP(12) * t170;
t161 = pkin(8) * t221;
t149 = t161 + (-pkin(2) - t226) * t173;
t214 = MDP(14) * t149;
t213 = MDP(19) * t178;
t212 = (MDP(25) * qJ(5));
t136 = t228 * t150 + t151 * t175;
t174 = sin(qJ(6));
t177 = cos(qJ(6));
t129 = t136 * t174 - t177 * t220;
t211 = MDP(26) * t129;
t210 = MDP(26) * t177;
t209 = MDP(31) * t174;
t208 = MDP(32) * t177;
t128 = t136 * t177 + t174 * t220;
t207 = t128 * MDP(29);
t206 = t136 * MDP(16);
t157 = t170 * t175 - t228 * t172;
t158 = t228 * t170 + t175 * t172;
t164 = -pkin(3) * t172 - pkin(2);
t188 = -qJ(5) * t158 + t164;
t139 = pkin(4) * t157 + t188;
t205 = t139 * MDP(23);
t204 = t164 * MDP(20);
t203 = MDP(15) + MDP(30);
t202 = MDP(20) - MDP(23);
t201 = MDP(24) - MDP(21);
t200 = 0.2e1 * t158;
t198 = qJ(5) * t220;
t197 = t177 * t174 * MDP(27);
t196 = MDP(23) - t224;
t159 = t225 * t170;
t160 = t225 * t172;
t140 = t228 * t159 + t160 * t175;
t195 = -t228 * t124 + t175 * t127;
t194 = -t131 * t170 + t132 * t172;
t193 = MDP(21) * t164 - MDP(24) * t139;
t192 = MDP(28) * t174 + MDP(29) * t177;
t130 = t229 * t157 + t188;
t133 = pkin(5) * t158 + t140;
t121 = -t130 * t174 + t133 * t177;
t122 = t130 * t177 + t133 * t174;
t191 = MDP(31) * t121 - MDP(32) * t122;
t190 = MDP(31) * t177 - t174 * MDP(32);
t189 = t208 + t209;
t116 = t198 - t119;
t162 = pkin(4) * t220;
t117 = t162 + t195;
t141 = -t175 * t159 + t228 * t160;
t187 = -MDP(16) + t192;
t186 = MDP(22) + t190;
t185 = -t189 - t201;
t138 = pkin(3) * t150 + t149;
t113 = t137 * pkin(5) + pkin(10) * t220 + t117;
t183 = -qJ(5) * t137 + t138;
t115 = t229 * t136 + t183;
t111 = t113 * t177 - t115 * t174;
t112 = t113 * t174 + t115 * t177;
t184 = t129 * MDP(28) + t111 * MDP(31) - t112 * MDP(32) + t207;
t182 = (-MDP(31) * t229 + MDP(28)) * t177 + (MDP(32) * t229 - MDP(29)) * t174;
t181 = -pkin(4) * MDP(22) + MDP(17) + t182;
t169 = t177 ^ 2;
t168 = t174 ^ 2;
t166 = t171 ^ 2;
t153 = t173 * t227 + t199;
t152 = t173 * t226 - t161;
t134 = -t157 * pkin(5) + t141;
t120 = pkin(4) * t136 + t183;
t114 = -pkin(5) * t136 - t116;
t1 = [(t116 ^ 2 + t117 ^ 2 + t120 ^ 2) * MDP(25) + t166 * t176 ^ 2 * MDP(4) + (t131 ^ 2 + t132 ^ 2 + t149 ^ 2) * MDP(14) + MDP(1) + (0.2e1 * t171 * t218 + t219) * t173 + t203 * t137 ^ 2 + (0.2e1 * t128 * MDP(27) + MDP(28) * t235 + t211) * t129 + (0.2e1 * MDP(5) * t176 + t213) * t166 * t178 + (t116 * t220 - t120 * t137) * t232 + 0.2e1 * (t119 * t220 + t137 * t138) * MDP(21) + 0.2e1 * (t136 * t138 + t195 * t220) * MDP(20) + 0.2e1 * (-t117 * t220 - t120 * t136) * MDP(23) + 0.2e1 * (t152 * t173 + t166 * t226) * MDP(9) + 0.2e1 * (-t131 * t220 + t149 * t150) * MDP(11) + 0.2e1 * (t132 * t220 + t149 * t151) * MDP(12) + 0.2e1 * (-t153 * t173 - t166 * t227) * MDP(10) + (-t131 * t151 - t132 * t150) * t234 + (t111 * t137 - t114 * t128) * t231 + (t116 * t136 + t117 * t137) * t233 + (-t112 * t137 + t114 * t129) * t230 + (-t206 + t207) * t235 + 0.2e1 * (-t137 * MDP(17) + t136 * MDP(18) + MDP(7) * t173) * t220; (-t116 * t141 + t117 * t140 + t120 * t139) * MDP(25) + t152 * MDP(9) - t153 * MDP(10) - pkin(2) * t214 + t194 * MDP(13) + t219 + (-pkin(2) * t150 - t149 * t172) * MDP(11) + (-pkin(2) * t151 + t149 * t170) * MDP(12) + (-t128 * MDP(31) + t129 * MDP(32)) * t134 + (MDP(22) * t140 + t191 + t193) * t137 + (-MDP(22) * t141 + t204 - t205) * t136 + (t218 + (t202 * t140 - t201 * t141 + MDP(7)) * t178) * t171 + (t194 * MDP(14) + (-t150 * t172 + t151 * t170) * MDP(13) + (MDP(11) * t170 + MDP(12) * t172) * t220) * qJ(3) + (-MDP(17) * t220 + t138 * MDP(21) + t117 * MDP(22) - t120 * MDP(24) + t203 * t137 + t184 - t206) * t158 + (t174 * t211 + MDP(18) * t220 + t116 * MDP(22) + (t128 * t174 + t129 * t177) * MDP(27) - t120 * MDP(23) + t138 * MDP(20) - t190 * t114 + t187 * t137) * t157; MDP(8) + (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) * MDP(25) + t193 * t200 + (-0.2e1 * t215 + 0.2e1 * t216 + t223) * pkin(2) + (t121 * t231 - t122 * t230 + t140 * t233 + t203 * t158) * t158 + (MDP(14) * qJ(3) + t234) * (t170 ^ 2 + t172 ^ 2) * qJ(3) + ((t174 * t230 - t177 * t231) * t134 - t141 * t233 + t187 * t200 + 0.2e1 * t204 - 0.2e1 * t205 + (MDP(26) * t168 + 0.2e1 * t197) * t157) * t157; t150 * MDP(11) + t151 * MDP(12) + MDP(25) * t120 + t202 * t136 + t185 * t137 + t214; MDP(25) * t139 + t202 * t157 + t185 * t158 + t215 - t216 - t223; MDP(14) + MDP(25); -t171 * t213 - t195 * MDP(20) - t119 * MDP(21) + (0.2e1 * t162 + t195) * MDP(23) + (-0.2e1 * t198 + t119) * MDP(24) + (-pkin(4) * t117 - qJ(5) * t116) * MDP(25) + t129 * t210 + (t128 * t177 - t129 * t174) * MDP(27) + (-qJ(5) * t128 + t114 * t174) * MDP(31) + (qJ(5) * t129 + t114 * t177) * MDP(32) + (-qJ(5) * MDP(22) - MDP(18)) * t136 + t181 * t137; (t201 + t212) * t141 + (-MDP(20) + t196) * t140 + t189 * t134 + t181 * t158 + (-MDP(18) + t174 * t210 + (-t168 + t169) * MDP(27) - t186 * qJ(5)) * t157; 0; -0.2e1 * t197 + t169 * MDP(26) + MDP(19) + (-0.2e1 * MDP(23) + t224) * pkin(4) + (t232 + 0.2e1 * t208 + 0.2e1 * t209 + t212) * qJ(5); -MDP(23) * t220 + MDP(25) * t117 + t186 * t137; MDP(25) * t140 + t186 * t158; 0; t196; MDP(25); MDP(30) * t137 + t184; MDP(30) * t158 + t192 * t157 + t191; -t189; t182; t190; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
