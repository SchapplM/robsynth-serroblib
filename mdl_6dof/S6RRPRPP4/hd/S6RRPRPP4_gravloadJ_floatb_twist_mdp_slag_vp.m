% Calculate Gravitation load on the joints for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:51
% EndTime: 2019-03-09 10:01:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (265->111), mult. (457->144), div. (0->0), fcn. (421->8), ass. (0->65)
t179 = -MDP(10) + MDP(13);
t178 = MDP(9) - MDP(12) + MDP(22) + MDP(25);
t131 = -qJ(5) - pkin(8);
t136 = cos(qJ(2));
t132 = sin(qJ(4));
t134 = sin(qJ(1));
t164 = t134 * t132;
t156 = pkin(4) * t164;
t133 = sin(qJ(2));
t167 = t133 * t134;
t177 = t131 * t167 + t136 * t156;
t137 = cos(qJ(1));
t168 = t132 * t137;
t119 = pkin(4) * t168;
t166 = t133 * t137;
t176 = t136 * t119 + t131 * t166;
t109 = g(1) * t137 + g(2) * t134;
t96 = g(3) * t133 + t109 * t136;
t175 = pkin(4) * t132;
t135 = cos(qJ(4));
t174 = pkin(4) * t135;
t173 = g(1) * t134;
t169 = g(3) * t136;
t126 = t136 * pkin(2);
t124 = t133 * qJ(3);
t130 = qJ(4) + pkin(9);
t123 = cos(t130);
t165 = t134 * t123;
t163 = t134 * t135;
t162 = t135 * t137;
t161 = t136 * t131;
t160 = t136 * t137;
t152 = t133 * t163;
t159 = pkin(4) * t152 + t119;
t158 = t126 + t124;
t157 = -MDP(23) - MDP(27);
t155 = t135 * t169;
t154 = t132 * t166;
t153 = t133 * t162;
t121 = pkin(3) + t174;
t127 = t137 * pkin(7);
t151 = t137 * t121 + t134 * t161 + t127;
t150 = -pkin(1) - t126;
t149 = t137 * pkin(1) + pkin(2) * t160 + t134 * pkin(7) + qJ(3) * t166;
t115 = t134 * t136 * qJ(3);
t148 = -pkin(2) * t167 + t115;
t117 = qJ(3) * t160;
t147 = -pkin(2) * t166 + t117;
t146 = pkin(4) * t153 - t156;
t122 = sin(t130);
t144 = pkin(5) * t122 - qJ(6) * t123;
t143 = t133 * t175 + t158 - t161;
t91 = t122 * t134 - t123 * t166;
t93 = t122 * t137 + t133 * t165;
t142 = g(1) * t91 - g(2) * t93 + t123 * t169;
t139 = pkin(4) * t154 + t134 * t121 - t131 * t160 + t149;
t138 = ((-qJ(3) - t175) * t133 + t150) * t173;
t100 = -t133 * t164 + t162;
t99 = t152 + t168;
t98 = t154 + t163;
t97 = t153 - t164;
t95 = t109 * t133 - t169;
t94 = -t122 * t167 + t123 * t137;
t92 = t122 * t166 + t165;
t1 = [(-g(1) * t127 - g(2) * t149 - (t150 - t124) * t173) * MDP(14) + (-g(1) * t100 - g(2) * t98) * MDP(20) + (g(1) * t99 - g(2) * t97) * MDP(21) + (-g(1) * t151 - g(2) * t139 - t138) * MDP(23) + (-g(1) * t94 - g(2) * t92) * MDP(24) + (-g(1) * t93 - g(2) * t91) * MDP(26) + (-g(1) * (pkin(5) * t94 + qJ(6) * t93 + t151) - g(2) * (t92 * pkin(5) + t91 * qJ(6) + t139) - t138) * MDP(27) + (MDP(3) - MDP(11)) * t109 + (t179 * t133 + t178 * t136 + MDP(2)) * (-g(2) * t137 + t173); (-g(1) * t147 - g(2) * t148 - g(3) * t158) * MDP(14) + (-g(1) * (t147 + t176) - g(2) * (t148 + t177) - g(3) * t143) * MDP(23) + (-g(1) * (t117 + t176) - g(2) * (t115 + t177) - g(3) * (t144 * t133 + t143) + t109 * (pkin(2) * t133 - t144 * t136)) * MDP(27) + t178 * t95 + (-t132 * MDP(20) - t135 * MDP(21) - t122 * MDP(24) + t123 * MDP(26) - t179) * t96; (-MDP(14) + t157) * t95; (-g(1) * t97 - g(2) * t99 + t155) * MDP(20) + (g(1) * t98 - g(2) * t100 - t132 * t169) * MDP(21) + (pkin(4) * t155 - g(1) * t146 - g(2) * t159) * MDP(23) + t142 * MDP(24) + (-g(1) * t92 + g(2) * t94 + t122 * t169) * MDP(26) + (-g(1) * (-pkin(5) * t91 + qJ(6) * t92 + t146) - g(2) * (pkin(5) * t93 - qJ(6) * t94 + t159) - (-pkin(5) * t123 - qJ(6) * t122 - t174) * t169) * MDP(27); t157 * t96; -t142 * MDP(27);];
taug  = t1;
