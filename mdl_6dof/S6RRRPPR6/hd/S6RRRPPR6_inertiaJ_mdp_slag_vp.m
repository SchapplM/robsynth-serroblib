% Calculate joint inertia matrix for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:54:07
% EndTime: 2019-03-09 15:54:10
% DurationCPUTime: 1.08s
% Computational Cost: add. (1637->229), mult. (3580->330), div. (0->0), fcn. (3946->10), ass. (0->102)
t165 = sin(pkin(6));
t225 = 0.2e1 * t165;
t172 = cos(qJ(3));
t212 = -qJ(4) - pkin(9);
t151 = t212 * t172;
t164 = sin(pkin(11));
t166 = cos(pkin(11));
t169 = sin(qJ(3));
t191 = t212 * t169;
t134 = -t151 * t164 - t166 * t191;
t136 = -t166 * t151 + t164 * t191;
t224 = t169 * MDP(13) + t172 * MDP(14) + t134 * MDP(21) + t136 * MDP(22) - (MDP(16) * t169 + MDP(17) * t172) * pkin(9);
t148 = t164 * t169 - t166 * t172;
t149 = t164 * t172 + t166 * t169;
t160 = -pkin(3) * t172 - pkin(2);
t183 = -qJ(5) * t149 + t160;
t133 = pkin(4) * t148 + t183;
t223 = -0.2e1 * t133;
t222 = 0.2e1 * MDP(16);
t221 = 2 * MDP(18);
t220 = 2 * MDP(20);
t219 = 0.2e1 * MDP(21);
t218 = 0.2e1 * MDP(22);
t217 = 2 * MDP(29);
t216 = 2 * MDP(30);
t215 = pkin(4) + pkin(10);
t170 = sin(qJ(2));
t214 = pkin(1) * t170;
t173 = cos(qJ(2));
t213 = pkin(1) * t173;
t129 = -t148 * pkin(5) + t136;
t211 = t129 * t148;
t210 = t165 * t170;
t209 = t165 * t173;
t167 = cos(pkin(6));
t208 = t167 * MDP(8);
t207 = t170 * MDP(6);
t194 = pkin(8) * t209;
t140 = t194 + (pkin(9) + t214) * t167;
t141 = (-pkin(2) * t173 - pkin(9) * t170 - pkin(1)) * t165;
t126 = -t140 * t169 + t172 * t141;
t144 = t167 * t169 + t172 * t210;
t119 = -pkin(3) * t209 - qJ(4) * t144 + t126;
t127 = t140 * t172 + t141 * t169;
t143 = -t167 * t172 + t169 * t210;
t122 = -qJ(4) * t143 + t127;
t114 = t164 * t119 + t166 * t122;
t159 = -pkin(3) * t166 - pkin(4);
t206 = MDP(23) * t159;
t205 = MDP(28) * t149;
t130 = t166 * t143 + t144 * t164;
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t123 = t130 * t171 + t168 * t209;
t204 = t123 * MDP(25);
t203 = t123 * MDP(27);
t124 = t130 * t168 - t171 * t209;
t202 = t124 * MDP(24);
t131 = -t143 * t164 + t144 * t166;
t201 = t131 * MDP(26);
t200 = t131 * MDP(28);
t199 = t144 * MDP(11);
t198 = t148 * MDP(21);
t197 = t168 * MDP(29);
t196 = t171 * MDP(24);
t195 = t171 * MDP(30);
t193 = t171 * t168 * MDP(25);
t192 = t134 ^ 2 + t136 ^ 2;
t113 = t119 * t166 - t164 * t122;
t112 = pkin(4) * t209 - t113;
t188 = t144 * MDP(13) - t143 * MDP(14);
t152 = pkin(8) * t210;
t139 = t152 + (-pkin(2) - t213) * t167;
t132 = pkin(3) * t143 + t139;
t178 = -qJ(5) * t131 + t132;
t115 = pkin(4) * t130 + t178;
t186 = -t130 * MDP(21) + t115 * MDP(23);
t185 = t171 * MDP(29) - t168 * MDP(30);
t184 = t195 + t197;
t111 = qJ(5) * t209 - t114;
t182 = MDP(20) + t185;
t181 = -MDP(22) - t184;
t180 = (MDP(26) * t168 + MDP(27) * t171) * t148;
t177 = t171 * MDP(26) - t168 * MDP(27) + t185 * (-pkin(10) + t159);
t108 = pkin(5) * t131 + pkin(10) * t209 + t112;
t110 = t130 * t215 + t178;
t106 = t108 * t171 - t110 * t168;
t107 = t108 * t168 + t110 * t171;
t176 = t124 * MDP(26) + t106 * MDP(29) - t107 * MDP(30) + t200 + t203;
t175 = t159 * MDP(20) + t177;
t163 = t171 ^ 2;
t162 = t168 ^ 2;
t161 = t165 ^ 2;
t155 = pkin(3) * t164 + qJ(5);
t146 = t167 * t214 + t194;
t145 = t167 * t213 - t152;
t128 = pkin(5) * t149 + t134;
t125 = t148 * t215 + t183;
t117 = t125 * t171 + t128 * t168;
t116 = -t125 * t168 + t128 * t171;
t109 = -pkin(5) * t130 - t111;
t1 = [t161 * t170 ^ 2 * MDP(4) + (t113 ^ 2 + t114 ^ 2 + t132 ^ 2) * MDP(19) + (t111 ^ 2 + t112 ^ 2 + t115 ^ 2) * MDP(23) + MDP(1) + (t207 * t225 + t208) * t167 + (-0.2e1 * t143 * MDP(12) + t199) * t144 + (t200 + 0.2e1 * t203) * t131 + (0.2e1 * t201 + t202 + 0.2e1 * t204) * t124 + ((t173 * MDP(15) + 0.2e1 * MDP(5) * t170) * t161 + (MDP(7) * t167 - t188) * t225) * t173 + (-t126 * t209 + t139 * t143) * t222 + 0.2e1 * (t127 * t209 + t139 * t144) * MDP(17) + (-t112 * t209 - t115 * t130) * t219 + (t111 * t209 - t115 * t131) * t218 + 0.2e1 * (t145 * t167 + t161 * t213) * MDP(9) + 0.2e1 * (-t146 * t167 - t161 * t214) * MDP(10) + (t111 * t130 + t112 * t131) * t220 + (t106 * t131 - t109 * t123) * t217 + (-t113 * t131 - t114 * t130) * t221 + (-t107 * t131 + t109 * t124) * t216; t208 + t145 * MDP(9) - t146 * MDP(10) + (-t143 * t169 + t144 * t172) * MDP(12) + (-pkin(2) * t143 - t139 * t172) * MDP(16) + (-pkin(2) * t144 + t139 * t169) * MDP(17) + (-t113 * t134 + t114 * t136 + t132 * t160) * MDP(19) + (-t111 * t136 + t112 * t134) * MDP(23) + (t116 * t131 - t123 * t129) * MDP(29) + (-t117 * t131 + t124 * t129) * MDP(30) + t169 * t199 + (-t131 * MDP(22) + t186) * t133 + (-t113 * MDP(18) + t112 * MDP(20) - t115 * MDP(22) + t176) * t149 + (-t114 * MDP(18) + t111 * MDP(20) - t115 * MDP(21) + (MDP(25) * t124 + MDP(27) * t131 - MDP(29) * t109) * t171 + (t109 * MDP(30) + t201 + t202 + t204) * t168) * t148 + (t207 + (MDP(7) - t224) * t173) * t165 + (MDP(18) + MDP(20)) * (-t130 * t136 + t131 * t134); MDP(8) + pkin(2) * t172 * t222 + (t160 ^ 2 + t192) * MDP(19) + t198 * t223 + (t133 ^ 2 + t192) * MDP(23) + (t162 * MDP(24) + 0.2e1 * t193) * t148 ^ 2 + (MDP(11) * t169 + 0.2e1 * t172 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t169 + (MDP(22) * t223 + 0.2e1 * t180 + t205) * t149 + (t116 * t149 - t171 * t211) * t217 + (-t117 * t149 + t168 * t211) * t216 + (t221 + t220) * (t134 * t149 - t136 * t148); t126 * MDP(16) - t127 * MDP(17) - t155 * t130 * MDP(20) + t112 * MDP(21) + t114 * MDP(22) + (-t111 * t155 + t112 * t159) * MDP(23) + t124 * t196 + (t123 * t171 - t124 * t168) * MDP(25) + (t109 * t168 - t123 * t155) * MDP(29) + (t109 * t171 + t124 * t155) * MDP(30) + (-MDP(15) - t159 * MDP(21) + (-qJ(5) - t155) * MDP(22)) * t209 + t175 * t131 + ((-t130 * t164 - t131 * t166) * MDP(18) + (t113 * t166 + t114 * t164) * MDP(19)) * pkin(3) + t188; (t134 * t159 + t136 * t155) * MDP(23) + t184 * t129 + t175 * t149 + (t168 * t196 + (-t162 + t163) * MDP(25) - t182 * t155) * t148 + ((-t148 * t164 - t149 * t166) * MDP(18) + (-t134 * t166 + t136 * t164) * MDP(19)) * pkin(3) + t224; -0.2e1 * t193 + t163 * MDP(24) + MDP(15) + (t164 ^ 2 + t166 ^ 2) * MDP(19) * pkin(3) ^ 2 + (t219 + t206) * t159 + (MDP(23) * t155 + 0.2e1 * t195 + 0.2e1 * t197 + t218) * t155; MDP(19) * t132 + t131 * t181 + t186; MDP(19) * t160 + MDP(23) * t133 + t149 * t181 - t198; 0; MDP(19) + MDP(23); -MDP(21) * t209 + t112 * MDP(23) + t131 * t182; t134 * MDP(23) + t149 * t182; MDP(21) + t206; 0; MDP(23); t176; t116 * MDP(29) - t117 * MDP(30) + t180 + t205; t177; -t184; t185; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
