% Calculate joint inertia matrix for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRRPPP1_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:34
% EndTime: 2019-03-09 15:18:40
% DurationCPUTime: 1.96s
% Computational Cost: add. (2643->342), mult. (5935->492), div. (0->0), fcn. (6113->8), ass. (0->102)
t171 = sin(pkin(6));
t174 = sin(qJ(3));
t176 = cos(qJ(3));
t153 = -qJ(4) * t171 * t174 - pkin(3) * t176 - pkin(2);
t173 = cos(pkin(6));
t215 = qJ(4) * t173;
t191 = pkin(9) + t215;
t155 = t191 * t174;
t172 = cos(pkin(10));
t230 = t172 * (t153 * t171 - t155 * t173);
t177 = cos(qJ(2));
t175 = sin(qJ(2));
t207 = t174 * t175;
t181 = t171 * t177 + t173 * t207;
t229 = -(MDP(16) * t174 + MDP(17) * t176) * pkin(9) + t174 * MDP(13) + t176 * MDP(14);
t158 = -pkin(2) * t177 - t175 * pkin(9) - pkin(1);
t217 = pkin(8) * t177;
t140 = t174 * t158 + t176 * t217;
t127 = -qJ(4) * t181 + t140;
t154 = t176 * t158;
t205 = t175 * t176;
t219 = pkin(8) * t174;
t128 = -t205 * t215 + t154 + (-pkin(3) - t219) * t177;
t211 = t171 * t176;
t137 = (pkin(3) * t174 - qJ(4) * t211 + pkin(8)) * t175;
t170 = sin(pkin(10));
t113 = -t170 * t127 + (t128 * t173 + t137 * t171) * t172;
t228 = 2 * MDP(18);
t227 = 2 * MDP(19);
t226 = 2 * MDP(20);
t225 = 2 * MDP(22);
t224 = 2 * MDP(23);
t223 = 2 * MDP(24);
t222 = 2 * MDP(26);
t221 = 2 * MDP(27);
t220 = 2 * MDP(28);
t218 = pkin(8) * t176;
t216 = pkin(4) + qJ(6);
t214 = t170 * t171;
t213 = t170 * t173;
t212 = t171 * t172;
t209 = t172 * t173;
t208 = t173 * t176;
t206 = t174 * t176;
t156 = t191 * t176;
t150 = t170 * t156;
t204 = pkin(4) * t211 + t150;
t149 = pkin(3) * t213 + qJ(4) * t212;
t146 = t170 * t174 - t172 * t208;
t129 = t173 * t153 + t155 * t171;
t147 = t170 * t208 + t172 * t174;
t179 = -qJ(5) * t147 + t129;
t116 = t146 * t216 + t179;
t203 = t116 * MDP(29);
t118 = pkin(4) * t146 + t179;
t202 = t118 * MDP(25);
t119 = -t128 * t171 + t173 * t137;
t201 = t119 * MDP(21);
t200 = t129 * MDP(21);
t199 = t177 * MDP(15);
t198 = MDP(22) + MDP(26);
t197 = MDP(23) - MDP(28);
t196 = MDP(25) + MDP(29);
t114 = t172 * t127 + t128 * t213 + t137 * t214;
t123 = t153 * t214 - t155 * t213 + t172 * t156;
t194 = -pkin(3) * t172 - pkin(4);
t193 = MDP(12) * t206;
t192 = -qJ(5) * t170 - pkin(3);
t190 = MDP(18) - t197;
t189 = MDP(19) - MDP(24) - MDP(27);
t138 = -qJ(5) * t173 - t149;
t186 = MDP(13) * t176 - MDP(14) * t174;
t183 = MDP(23) * t172 - MDP(24) * t170;
t182 = -MDP(27) * t170 - MDP(28) * t172;
t152 = t171 * t207 - t177 * t173;
t110 = -qJ(5) * t152 - t114;
t132 = -t181 * t170 + t172 * t205;
t180 = -qJ(5) * t132 + t119;
t120 = qJ(5) * t211 - t123;
t169 = t176 ^ 2;
t168 = t175 ^ 2;
t167 = t174 ^ 2;
t166 = t171 ^ 2;
t161 = qJ(4) * t214;
t148 = pkin(3) * t209 - t161;
t142 = (-pkin(4) * t172 + t192) * t171;
t141 = t173 * t194 + t161;
t139 = -t174 * t217 + t154;
t134 = (-t172 * t216 + t192) * t171;
t133 = pkin(5) * t212 - t138;
t131 = t170 * t205 + t172 * t181;
t130 = pkin(5) * t214 + t161 + (-qJ(6) + t194) * t173;
t122 = -t150 + t230;
t121 = t204 - t230;
t117 = -pkin(5) * t146 - t120;
t115 = t155 * t209 + pkin(5) * t147 + (qJ(6) * t176 - t153 * t172) * t171 + t204;
t112 = pkin(4) * t131 + t180;
t111 = -pkin(4) * t152 - t113;
t109 = t131 * t216 + t180;
t108 = -pkin(5) * t131 - t110;
t107 = pkin(5) * t132 - t152 * t216 - t113;
t1 = [(t113 ^ 2 + t114 ^ 2 + t119 ^ 2) * MDP(21) + (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) * MDP(25) + (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) * MDP(29) - 0.2e1 * pkin(1) * t175 * MDP(10) + MDP(1) + (t169 * MDP(11) + MDP(4) - 0.2e1 * t193) * t168 + (t199 + 0.2e1 * pkin(1) * MDP(9) + 0.2e1 * (MDP(5) - t186) * t175) * t177 + 0.2e1 * (-t139 * t177 + t168 * t219) * MDP(16) + 0.2e1 * (t140 * t177 + t168 * t218) * MDP(17) + (t113 * t152 + t119 * t131) * t228 + (-t114 * t152 + t119 * t132) * t227 + (-t113 * t132 - t114 * t131) * t226 + (t110 * t131 + t111 * t132) * t225 + (t111 * t152 - t112 * t131) * t224 + (-t110 * t152 - t112 * t132) * t223 + (t107 * t132 - t108 * t131) * t222 + (t108 * t152 - t109 * t132) * t221 + (-t107 * t152 + t109 * t131) * t220; (-t113 * t211 + t119 * t146 + t122 * t152 + t129 * t131) * MDP(18) + (t114 * t211 + t119 * t147 - t123 * t152 + t129 * t132) * MDP(19) + (-t113 * t147 - t114 * t146 - t122 * t132 - t123 * t131) * MDP(20) + (t113 * t122 + t114 * t123 + t119 * t129) * MDP(21) + (t110 * t146 + t111 * t147 + t120 * t131 + t121 * t132) * MDP(22) + (-t111 * t211 - t112 * t146 - t118 * t131 + t121 * t152) * MDP(23) + (t110 * t211 - t112 * t147 - t118 * t132 - t120 * t152) * MDP(24) + (t110 * t120 + t111 * t121 + t112 * t118) * MDP(25) + (t107 * t147 - t108 * t146 + t115 * t132 - t117 * t131) * MDP(26) + (-t108 * t211 - t109 * t147 - t116 * t132 + t117 * t152) * MDP(27) + (t107 * t211 + t109 * t146 - t115 * t152 + t116 * t131) * MDP(28) + (t107 * t115 + t108 * t117 + t109 * t116) * MDP(29) + (-pkin(8) * MDP(10) + MDP(7) - t229) * t177 + (MDP(6) - pkin(8) * MDP(9) + MDP(11) * t206 + (-t167 + t169) * MDP(12) + (-pkin(2) * t174 - t218) * MDP(16) + (-pkin(2) * t176 + t219) * MDP(17)) * t175; MDP(8) + t167 * MDP(11) + 0.2e1 * t193 + (t122 ^ 2 + t123 ^ 2 + t129 ^ 2) * MDP(21) + (t118 ^ 2 + t120 ^ 2 + t121 ^ 2) * MDP(25) + (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) * MDP(29) + 0.2e1 * (t176 * MDP(16) - t174 * MDP(17)) * pkin(2) + (-t122 * t211 + t129 * t146) * t228 + (t123 * t211 + t129 * t147) * t227 + (-t122 * t147 - t123 * t146) * t226 + (t120 * t146 + t121 * t147) * t225 + (-t118 * t146 - t121 * t211) * t224 + (-t118 * t147 + t120 * t211) * t223 + (t115 * t147 - t117 * t146) * t222 + (-t116 * t147 - t117 * t211) * t221 + (t115 * t211 + t116 * t146) * t220; -t199 + t139 * MDP(16) - t140 * MDP(17) + (t113 * t173 + t148 * t152) * MDP(18) + (-t114 * t173 - t149 * t152) * MDP(19) + (-t131 * t149 - t132 * t148) * MDP(20) + (t113 * t148 + t114 * t149) * MDP(21) + (t131 * t138 + t132 * t141) * MDP(22) + (t111 * t173 - t131 * t142 + t141 * t152) * MDP(23) + (-t110 * t173 - t132 * t142 - t138 * t152) * MDP(24) + (t110 * t138 + t111 * t141 + t112 * t142) * MDP(25) + (t130 * t132 - t131 * t133) * MDP(26) + (t108 * t173 - t132 * t134 + t133 * t152) * MDP(27) + (-t107 * t173 - t130 * t152 + t131 * t134) * MDP(28) + (t107 * t130 + t108 * t133 + t109 * t134) * MDP(29) + t186 * t175 + ((-pkin(3) * t131 - t119 * t172) * MDP(18) + (-pkin(3) * t132 + t119 * t170) * MDP(19) + (-t113 * t170 + t114 * t172) * MDP(20) - pkin(3) * t201 + (-t110 * t172 + t111 * t170) * MDP(22) + (t107 * t170 + t108 * t172) * MDP(26) + t183 * t112 + t182 * t109) * t171; (-t146 * t149 - t147 * t148) * MDP(20) + (t122 * t148 + t123 * t149) * MDP(21) + (t138 * t146 + t141 * t147) * MDP(22) + (t120 * t138 + t121 * t141) * MDP(25) + (t130 * t147 - t133 * t146) * MDP(26) + (t115 * t130 + t117 * t133) * MDP(29) + (-t146 * MDP(23) - t147 * MDP(24) + t202) * t142 + (-t147 * MDP(27) + t146 * MDP(28) + t203) * t134 + (t122 * MDP(18) - t123 * MDP(19) + t121 * MDP(23) - t120 * MDP(24) + t117 * MDP(27) - t115 * MDP(28)) * t173 + ((-pkin(3) * t146 - t129 * t172 - t148 * t176) * MDP(18) + (-pkin(3) * t147 + t129 * t170 + t149 * t176) * MDP(19) + (-t122 * t170 + t123 * t172) * MDP(20) - pkin(3) * t200 + (-t120 * t172 + t121 * t170) * MDP(22) + (t118 * t172 - t141 * t176) * MDP(23) + (-t118 * t170 + t138 * t176) * MDP(24) + (t115 * t170 + t117 * t172) * MDP(26) + (-t116 * t170 - t133 * t176) * MDP(27) + (-t116 * t172 + t130 * t176) * MDP(28)) * t171 + t229; MDP(15) + (pkin(3) ^ 2 * t166 + t148 ^ 2 + t149 ^ 2) * MDP(21) + (t138 ^ 2 + t141 ^ 2 + t142 ^ 2) * MDP(25) + (t130 ^ 2 + t133 ^ 2 + t134 ^ 2) * MDP(29) + 0.2e1 * (MDP(18) * t172 - MDP(19) * t170) * t166 * pkin(3) + 0.2e1 * (MDP(18) * t148 - MDP(19) * t149 + MDP(23) * t141 - MDP(24) * t138 + MDP(27) * t133 - MDP(28) * t130) * t173 + 0.2e1 * ((-t148 * t170 + t149 * t172) * MDP(20) + (-t138 * t172 + t141 * t170) * MDP(22) + (t130 * t170 + t133 * t172) * MDP(26) + t183 * t142 + t182 * t134) * t171; t112 * MDP(25) + t109 * MDP(29) + t131 * t190 + t132 * t189 + t201; t146 * t190 + t147 * t189 + t200 + t202 + t203; MDP(25) * t142 + MDP(29) * t134 + (-MDP(21) * pkin(3) + t170 * t189 - t172 * t190) * t171; MDP(21) + t196; t111 * MDP(25) + t107 * MDP(29) + t132 * t198 + t152 * t197; t121 * MDP(25) + t115 * MDP(29) + t147 * t198 - t197 * t211; MDP(25) * t141 + MDP(29) * t130 + t173 * t197 + t198 * t214; 0; t196; -t131 * MDP(26) + t152 * MDP(27) + t108 * MDP(29); -t146 * MDP(26) - MDP(27) * t211 + t117 * MDP(29); MDP(26) * t212 + MDP(27) * t173 + MDP(29) * t133; 0; 0; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
