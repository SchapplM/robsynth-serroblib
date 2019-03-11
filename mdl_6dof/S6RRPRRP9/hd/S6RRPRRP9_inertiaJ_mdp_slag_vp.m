% Calculate joint inertia matrix for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP9_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:34:01
% EndTime: 2019-03-09 12:34:05
% DurationCPUTime: 1.08s
% Computational Cost: add. (1886->235), mult. (4177->341), div. (0->0), fcn. (4687->10), ass. (0->108)
t172 = sin(pkin(6));
t232 = 0.2e1 * t172;
t231 = 2 * MDP(13);
t230 = 2 * MDP(20);
t229 = 2 * MDP(21);
t228 = 2 * MDP(29);
t227 = cos(qJ(4));
t177 = sin(qJ(2));
t226 = pkin(1) * t177;
t179 = cos(qJ(2));
t225 = pkin(1) * t179;
t224 = pkin(9) + qJ(3);
t223 = -qJ(6) - pkin(10);
t222 = MDP(30) * pkin(5);
t221 = pkin(2) * MDP(14);
t171 = sin(pkin(11));
t173 = cos(pkin(11));
t176 = sin(qJ(4));
t154 = t227 * t171 + t176 * t173;
t220 = qJ(6) * t154;
t174 = cos(pkin(6));
t218 = t172 * t177;
t147 = t171 * t218 - t173 * t174;
t148 = t171 * t174 + t173 * t218;
t134 = -t176 * t147 + t227 * t148;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t217 = t172 * t179;
t129 = t134 * t178 - t175 * t217;
t219 = t129 * t178;
t216 = t174 * MDP(8);
t215 = t177 * MDP(6);
t155 = t224 * t171;
t156 = t224 * t173;
t139 = -t176 * t155 + t227 * t156;
t214 = t178 * t139;
t195 = pkin(8) * t217;
t143 = t195 + (qJ(3) + t226) * t174;
t144 = (-pkin(2) * t179 - qJ(3) * t177 - pkin(1)) * t172;
t131 = t173 * t143 + t171 * t144;
t169 = t175 ^ 2;
t170 = t178 ^ 2;
t212 = t169 + t170;
t211 = MDP(16) * t134;
t210 = MDP(19) * t179;
t209 = MDP(22) * t129;
t208 = MDP(22) * t178;
t207 = MDP(25) * t175;
t128 = t134 * t175 + t178 * t217;
t206 = t128 * MDP(25);
t205 = t129 * MDP(24);
t133 = t227 * t147 + t148 * t176;
t204 = t133 * MDP(26);
t203 = t134 * MDP(17);
t202 = t134 * MDP(21);
t201 = t139 * MDP(21);
t161 = pkin(8) * t218;
t146 = t161 + (-pkin(2) - t225) * t174;
t200 = t146 * MDP(14);
t153 = t171 * t176 - t227 * t173;
t199 = t153 * MDP(26);
t198 = t171 * MDP(12);
t197 = t173 * MDP(11);
t196 = t175 * MDP(28);
t164 = -pkin(3) * t173 - pkin(2);
t194 = MDP(18) * t217;
t193 = t175 * t178 * MDP(23);
t192 = -(MDP(29) * pkin(5)) + MDP(24);
t130 = -t143 * t171 + t173 * t144;
t122 = -pkin(3) * t217 - pkin(9) * t148 + t130;
t126 = -pkin(9) * t147 + t131;
t117 = t176 * t122 + t227 * t126;
t115 = -pkin(10) * t217 + t117;
t136 = t147 * pkin(3) + t146;
t119 = t133 * pkin(4) - t134 * pkin(10) + t136;
t111 = -t115 * t175 + t178 * t119;
t137 = pkin(4) * t153 - pkin(10) * t154 + t164;
t123 = t178 * t137 - t139 * t175;
t138 = t227 * t155 + t156 * t176;
t116 = t227 * t122 - t176 * t126;
t109 = pkin(5) * t133 - qJ(6) * t129 + t111;
t112 = t115 * t178 + t119 * t175;
t110 = -qJ(6) * t128 + t112;
t191 = t109 * t178 + t110 * t175;
t120 = pkin(5) * t153 - t178 * t220 + t123;
t121 = t214 + (t137 - t220) * t175;
t190 = t120 * t178 + t121 * t175;
t189 = -t130 * t171 + t131 * t173;
t157 = t223 * t175;
t158 = t223 * t178;
t188 = t157 * t178 - t158 * t175;
t124 = t137 * t175 + t214;
t187 = t123 * MDP(27) - t124 * MDP(28);
t186 = MDP(27) * t178 - t196;
t185 = MDP(27) * t175 + MDP(28) * t178;
t114 = pkin(4) * t217 - t116;
t184 = MDP(24) * t178 - MDP(16) - t207;
t183 = MDP(20) + t186;
t182 = t111 * MDP(27) - t112 * MDP(28) + t204 + t205 - t206;
t181 = t175 * MDP(24) + t178 * MDP(25) - t185 * pkin(10) - MDP(18);
t167 = t172 ^ 2;
t165 = -pkin(5) * t178 - pkin(4);
t150 = t174 * t226 + t195;
t149 = t174 * t225 - t161;
t132 = pkin(5) * t154 * t175 + t138;
t127 = t175 * t128;
t113 = t128 * pkin(5) + t114;
t1 = [(t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(30) + t167 * t177 ^ 2 * MDP(4) + t134 ^ 2 * MDP(15) + (t130 ^ 2 + t131 ^ 2 + t146 ^ 2) * MDP(14) + MDP(1) + (t215 * t232 + t216) * t174 + (-0.2e1 * MDP(23) * t128 + t209) * t129 + ((MDP(7) * t174 - t203) * t232 + (0.2e1 * MDP(5) * t177 + t210) * t167) * t179 + (0.2e1 * t194 + t204 + 0.2e1 * t205 - 0.2e1 * t206 - 0.2e1 * t211) * t133 + 0.2e1 * (-t130 * t217 + t146 * t147) * MDP(11) + 0.2e1 * (t131 * t217 + t146 * t148) * MDP(12) + (-t116 * t217 + t133 * t136) * t230 + (t117 * t217 + t134 * t136) * t229 + 0.2e1 * (t149 * t174 + t167 * t225) * MDP(9) + 0.2e1 * (-t150 * t174 - t167 * t226) * MDP(10) + (-t130 * t148 - t131 * t147) * t231 + 0.2e1 * (t111 * t133 + t114 * t128) * MDP(27) + 0.2e1 * (-t112 * t133 + t114 * t129) * MDP(28) + (-t109 * t129 - t110 * t128) * t228; (t123 * t133 + t128 * t138) * MDP(27) + (-t124 * t133 + t129 * t138) * MDP(28) + (-t120 * t129 - t121 * t128) * MDP(29) + (t109 * t120 + t110 * t121 + t113 * t132) * MDP(30) + t216 + t149 * MDP(9) - t150 * MDP(10) + (-pkin(2) * t147 - t146 * t173) * MDP(11) + (-pkin(2) * t148 + t146 * t171) * MDP(12) + t189 * MDP(13) - pkin(2) * t200 + (t133 * MDP(20) + t202) * t164 + (t215 + (t138 * MDP(20) + MDP(7) + t201) * t179) * t172 + ((-t147 * t173 + t148 * t171) * MDP(13) + t189 * MDP(14) + (MDP(11) * t171 + MDP(12) * t173) * t217) * qJ(3) + (MDP(20) * t136 + t182 + t194 - t211) * t153 + (-t191 * MDP(29) + t136 * MDP(21) + (-t128 * t178 - t129 * t175) * MDP(23) + t134 * MDP(15) - MDP(17) * t217 + t129 * t208 + t185 * t114 + t184 * t133) * t154; MDP(8) + (t120 ^ 2 + t121 ^ 2 + t132 ^ 2) * MDP(30) + (0.2e1 * t197 - 0.2e1 * t198 + t221) * pkin(2) + (-0.2e1 * t190 * MDP(29) + 0.2e1 * t185 * t138 + t164 * t229 + (MDP(22) * t170 + MDP(15) - 0.2e1 * t193) * t154) * t154 + (0.2e1 * t154 * t184 + t164 * t230 + 0.2e1 * t187 + t199) * t153 + (MDP(14) * qJ(3) + t231) * (t171 ^ 2 + t173 ^ 2) * qJ(3); t147 * MDP(11) + t148 * MDP(12) + t200 + t202 + (-t127 - t219) * MDP(29) + t191 * MDP(30) + t183 * t133; -t197 + t198 - t221 + t190 * MDP(30) + (-MDP(29) * t212 + MDP(21)) * t154 + t183 * t153; MDP(30) * t212 + MDP(14); t203 - t172 * t210 + t116 * MDP(20) - t117 * MDP(21) + t175 * t209 + (-t127 + t219) * MDP(23) + (-pkin(4) * t128 - t114 * t178) * MDP(27) + (-pkin(4) * t129 + t114 * t175) * MDP(28) + (-t109 * t175 + t110 * t178 + t128 * t158 - t129 * t157) * MDP(29) + (t109 * t157 - t110 * t158 + t113 * t165) * MDP(30) + t181 * t133; -t201 + (-t120 * t175 + t121 * t178) * MDP(29) + (t120 * t157 - t121 * t158 + t132 * t165) * MDP(30) - t183 * t138 + t181 * t153 + (MDP(17) + t175 * t208 + (-t169 + t170) * MDP(23) - t188 * MDP(29) - t185 * pkin(4)) * t154; t188 * MDP(30); MDP(19) + t169 * MDP(22) + 0.2e1 * t193 + (-t157 * t175 - t158 * t178) * t228 + (t157 ^ 2 + t158 ^ 2 + t165 ^ 2) * MDP(30) + 0.2e1 * t186 * pkin(4); (-MDP(29) * t129 + MDP(30) * t109) * pkin(5) + t182; t120 * t222 + t199 + (t178 * t192 - t207) * t154 + t187; -t196 + (MDP(27) + t222) * t178; t157 * t222 + (-MDP(28) * pkin(10) + MDP(25)) * t178 + (-MDP(27) * pkin(10) + t192) * t175; MDP(30) * (pkin(5) ^ 2) + MDP(26); t113 * MDP(30); t132 * MDP(30); 0; t165 * MDP(30); 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
