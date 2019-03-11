% Calculate joint inertia matrix for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:45
% EndTime: 2019-03-09 05:17:47
% DurationCPUTime: 0.63s
% Computational Cost: add. (1339->162), mult. (2572->244), div. (0->0), fcn. (2996->10), ass. (0->83)
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t177 = -t168 * MDP(20) - t171 * MDP(21);
t198 = -qJ(5) - pkin(8);
t153 = t198 * t168;
t154 = t198 * t171;
t163 = sin(pkin(11));
t165 = cos(pkin(11));
t132 = t165 * t153 + t154 * t163;
t146 = t163 * t171 + t165 * t168;
t118 = -pkin(9) * t146 + t132;
t133 = t163 * t153 - t165 * t154;
t144 = -t163 * t168 + t165 * t171;
t119 = pkin(9) * t144 + t133;
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t128 = -t170 * t144 + t146 * t167;
t129 = t144 * t167 + t146 * t170;
t180 = t129 * MDP(26) - t128 * MDP(27) + (t118 * t170 - t119 * t167) * MDP(29) - (t118 * t167 + t119 * t170) * MDP(30);
t207 = t168 * MDP(17) + t171 * MDP(18) + pkin(8) * t177 + t180;
t166 = cos(pkin(10));
t157 = -pkin(2) * t166 - pkin(1);
t206 = 0.2e1 * t157;
t205 = 2 * MDP(22);
t204 = -2 * MDP(25);
t203 = 0.2e1 * MDP(30);
t202 = cos(qJ(3));
t201 = pkin(1) * MDP(7);
t200 = pkin(4) * t163;
t199 = pkin(7) + qJ(2);
t164 = sin(pkin(10));
t151 = t199 * t164;
t152 = t199 * t166;
t169 = sin(qJ(3));
t131 = -t169 * t151 + t152 * t202;
t197 = t131 * t171;
t147 = t164 * t202 + t169 * t166;
t196 = t147 * t168;
t195 = t147 * t171;
t193 = t164 * MDP(5);
t191 = t166 * MDP(4);
t190 = t168 * t171;
t145 = t164 * t169 - t166 * t202;
t127 = pkin(3) * t145 - pkin(8) * t147 + t157;
t114 = t171 * t127 - t131 * t168;
t102 = pkin(4) * t145 - qJ(5) * t195 + t114;
t110 = t197 + (-qJ(5) * t147 + t127) * t168;
t99 = t163 * t102 + t165 * t110;
t123 = t128 * MDP(29);
t189 = -t129 * MDP(30) - t123;
t120 = t146 * t147;
t121 = t144 * t147;
t111 = t170 * t120 + t121 * t167;
t108 = t111 * MDP(27);
t112 = -t120 * t167 + t121 * t170;
t109 = t112 * MDP(26);
t187 = t129 * MDP(24);
t156 = pkin(4) * t165 + pkin(5);
t137 = t156 * t170 - t167 * t200;
t186 = t137 * MDP(29);
t138 = t156 * t167 + t170 * t200;
t185 = t138 * MDP(30);
t184 = t147 * MDP(14);
t183 = MDP(19) + MDP(28);
t182 = t145 * MDP(28) - t108 + t109;
t158 = -pkin(4) * t171 - pkin(3);
t181 = MDP(16) * t190;
t98 = t165 * t102 - t110 * t163;
t96 = pkin(5) * t145 - pkin(9) * t121 + t98;
t97 = -pkin(9) * t120 + t99;
t93 = -t167 * t97 + t170 * t96;
t130 = t151 * t202 + t169 * t152;
t94 = t167 * t96 + t170 * t97;
t179 = t171 * MDP(17) - t168 * MDP(18);
t178 = t171 * MDP(20) - t168 * MDP(21);
t116 = pkin(4) * t196 + t130;
t175 = -MDP(13) - t178;
t162 = t171 ^ 2;
t161 = t168 ^ 2;
t134 = -pkin(5) * t144 + t158;
t115 = t127 * t168 + t197;
t113 = t120 * pkin(5) + t116;
t1 = [(t116 ^ 2 + t98 ^ 2 + t99 ^ 2) * MDP(23) + MDP(1) + t184 * t206 + t183 * t145 ^ 2 + (t112 * MDP(24) + t111 * t204) * t112 + (0.2e1 * t191 - 0.2e1 * t193 + t201) * pkin(1) + (t162 * MDP(15) + MDP(8) - 0.2e1 * t181) * t147 ^ 2 + (MDP(13) * t206 + 0.2e1 * t109 - 0.2e1 * t108 + 0.2e1 * (-MDP(9) + t179) * t147) * t145 + (-t120 * t99 - t121 * t98) * t205 + 0.2e1 * (t111 * t113 + t145 * t93) * MDP(29) + (t112 * t113 - t145 * t94) * t203 + 0.2e1 * (t114 * t145 + t130 * t196) * MDP(20) + 0.2e1 * (-t115 * t145 + t130 * t195) * MDP(21) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t164 ^ 2 + t166 ^ 2) * qJ(2); -t191 + t193 - t201 + t184 + (-t146 * t120 - t144 * t121) * MDP(22) + (t98 * t144 + t99 * t146) * MDP(23) + (-t175 + t189) * t145; MDP(7) + (t144 ^ 2 + t146 ^ 2) * MDP(23); -t131 * MDP(14) + (-t120 * t133 - t121 * t132 + t144 * t99 - t146 * t98) * MDP(22) + (t116 * t158 + t132 * t98 + t133 * t99) * MDP(23) + t112 * t187 + (-t111 * t129 - t112 * t128) * MDP(25) + (t111 * t134 + t113 * t128) * MDP(29) + (t112 * t134 + t113 * t129) * MDP(30) + t175 * t130 + (MDP(10) + MDP(15) * t190 + (-t161 + t162) * MDP(16) + t177 * pkin(3)) * t147 + (-MDP(11) + t207) * t145; (t132 * t144 + t133 * t146) * MDP(23); MDP(12) + t161 * MDP(15) + 0.2e1 * t181 + (-t132 * t146 + t133 * t144) * t205 + (t132 ^ 2 + t133 ^ 2 + t158 ^ 2) * MDP(23) + 0.2e1 * t134 * t123 + 0.2e1 * t178 * pkin(3) + (t128 * t204 + t134 * t203 + t187) * t129; t145 * MDP(19) + t114 * MDP(20) - t115 * MDP(21) + (t137 * t145 + t93) * MDP(29) + (-t138 * t145 - t94) * MDP(30) + t179 * t147 + ((-t120 * t163 - t121 * t165) * MDP(22) + (t163 * t99 + t165 * t98) * MDP(23)) * pkin(4) + t182; (t144 * t165 + t146 * t163) * MDP(23) * pkin(4) + t178 + t189; ((t144 * t163 - t146 * t165) * MDP(22) + (t132 * t165 + t133 * t163) * MDP(23)) * pkin(4) + t207; (t163 ^ 2 + t165 ^ 2) * MDP(23) * pkin(4) ^ 2 + 0.2e1 * t186 - 0.2e1 * t185 + t183; MDP(23) * t116 + t111 * MDP(29) + t112 * MDP(30); 0; MDP(23) * t158 - t189; 0; MDP(23); t93 * MDP(29) - t94 * MDP(30) + t182; t189; t180; MDP(28) - t185 + t186; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
