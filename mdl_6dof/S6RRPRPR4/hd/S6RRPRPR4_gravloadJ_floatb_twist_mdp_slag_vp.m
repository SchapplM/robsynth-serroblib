% Calculate Gravitation load on the joints for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:19
% EndTime: 2019-03-09 10:26:21
% DurationCPUTime: 0.75s
% Computational Cost: add. (395->119), mult. (877->200), div. (0->0), fcn. (1070->14), ass. (0->62)
t130 = sin(qJ(6));
t134 = cos(qJ(6));
t125 = qJ(4) + pkin(12);
t123 = sin(t125);
t124 = cos(t125);
t127 = sin(pkin(6));
t137 = cos(qJ(1));
t157 = t127 * t137;
t126 = sin(pkin(11));
t132 = sin(qJ(2));
t136 = cos(qJ(2));
t162 = cos(pkin(11));
t115 = -t136 * t126 - t132 * t162;
t128 = cos(pkin(6));
t108 = t115 * t128;
t133 = sin(qJ(1));
t142 = -t132 * t126 + t136 * t162;
t97 = -t108 * t137 + t133 * t142;
t90 = t123 * t157 - t97 * t124;
t139 = t128 * t142;
t96 = t133 * t115 + t137 * t139;
t170 = t130 * t90 - t134 * t96;
t169 = t130 * t96 + t134 * t90;
t152 = t133 * t136;
t155 = t132 * t137;
t112 = -t128 * t152 - t155;
t158 = t127 * t136;
t167 = -g(1) * t112 - g(3) * t158;
t161 = t124 * t130;
t160 = t124 * t134;
t159 = t127 * t133;
t122 = pkin(2) * t136 + pkin(1);
t154 = t133 * t122;
t153 = t133 * t132;
t151 = t136 * t137;
t149 = t128 * t151;
t131 = sin(qJ(4));
t135 = cos(qJ(4));
t147 = -t131 * t157 + t135 * t97;
t109 = pkin(2) * t128 * t132 + (-pkin(8) - qJ(3)) * t127;
t146 = pkin(4) * t127 * t131 - t109;
t144 = g(1) * t133 - g(2) * t137;
t98 = -t133 * t108 - t137 * t142;
t143 = t97 * t131 + t135 * t157;
t93 = t131 * t98 + t135 * t159;
t106 = t142 * t127;
t99 = t115 * t137 - t133 * t139;
t140 = g(1) * t99 + g(2) * t96 + g(3) * t106;
t129 = -qJ(5) - pkin(9);
t121 = pkin(4) * t135 + pkin(3);
t118 = t137 * t122;
t116 = pkin(2) * t149;
t113 = -t128 * t153 + t151;
t111 = -t128 * t155 - t152;
t110 = -t149 + t153;
t107 = t115 * t127;
t102 = -t107 * t124 + t123 * t128;
t94 = t131 * t159 - t135 * t98;
t92 = t123 * t159 - t124 * t98;
t87 = -t130 * t99 + t134 * t92;
t86 = -t130 * t92 - t134 * t99;
t1 = [t144 * MDP(2) + (-g(1) * t111 - g(2) * t113) * MDP(9) + (-g(1) * t110 - g(2) * t112) * MDP(10) + (-g(1) * (-t109 * t137 - t154) - g(2) * (-t109 * t133 + t118)) * MDP(12) + (g(1) * t147 - g(2) * t94) * MDP(18) + (-g(1) * t143 - g(2) * t93) * MDP(19) + (-g(1) * t96 + g(2) * t99) * MDP(20) + (-g(1) * (-t97 * t121 - t96 * t129 + t137 * t146 - t154) - g(2) * (-t121 * t98 + t129 * t99 + t133 * t146 + t118)) * MDP(21) + (-g(1) * t169 - g(2) * t87) * MDP(27) + (g(1) * t170 - g(2) * t86) * MDP(28) + (-t127 * MDP(11) + MDP(3)) * (g(1) * t137 + g(2) * t133); (g(2) * t110 + t167) * MDP(9) + (g(3) * t127 * t132 + g(1) * t113 - g(2) * t111) * MDP(10) + (-g(2) * t116 + (g(2) * t153 + t167) * pkin(2)) * MDP(12) + (g(1) * t98 - g(2) * t97 + g(3) * t107) * MDP(20) + (-g(1) * (pkin(2) * t112 + t99 * t121 + t98 * t129) - g(2) * (-pkin(2) * t153 + t121 * t96 - t129 * t97 + t116) - g(3) * (pkin(2) * t158 + t106 * t121 + t107 * t129)) * MDP(21) + (-g(1) * (-t130 * t98 + t160 * t99) - g(2) * (t130 * t97 + t160 * t96) - g(3) * (t106 * t160 - t107 * t130)) * MDP(27) + (-g(1) * (-t134 * t98 - t161 * t99) - g(2) * (t134 * t97 - t96 * t161) - g(3) * (-t106 * t161 - t107 * t134)) * MDP(28) + (-MDP(18) * t135 + MDP(19) * t131) * t140; (MDP(12) + MDP(21)) * (-g(3) * t128 - t127 * t144); (g(1) * t94 + g(2) * t147 - g(3) * (t107 * t135 - t128 * t131)) * MDP(19) + (-MDP(27) * t134 + MDP(28) * t130) * (g(1) * (t123 * t98 + t124 * t159) + g(2) * (-t97 * t123 - t124 * t157) + g(3) * (t107 * t123 + t124 * t128)) + (pkin(4) * MDP(21) + MDP(18)) * (g(2) * t143 - g(3) * (t107 * t131 + t128 * t135) - g(1) * t93); t140 * MDP(21); (-g(1) * t86 - g(2) * t170 - g(3) * (-t102 * t130 - t106 * t134)) * MDP(27) + (g(1) * t87 - g(2) * t169 - g(3) * (-t102 * t134 + t106 * t130)) * MDP(28);];
taug  = t1;
