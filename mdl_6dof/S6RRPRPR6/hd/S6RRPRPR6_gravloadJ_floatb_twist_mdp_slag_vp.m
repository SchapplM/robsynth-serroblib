% Calculate Gravitation load on the joints for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:56
% EndTime: 2019-03-09 10:42:59
% DurationCPUTime: 1.02s
% Computational Cost: add. (467->117), mult. (1170->192), div. (0->0), fcn. (1453->12), ass. (0->58)
t179 = MDP(18) - MDP(21);
t174 = MDP(19) - MDP(22);
t129 = sin(pkin(11));
t134 = sin(qJ(2));
t138 = cos(qJ(2));
t168 = cos(pkin(11));
t120 = -t138 * t129 - t134 * t168;
t135 = sin(qJ(1));
t139 = cos(qJ(1));
t131 = cos(pkin(6));
t144 = -t134 * t129 + t138 * t168;
t140 = t131 * t144;
t101 = t135 * t120 + t139 * t140;
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t113 = t120 * t131;
t102 = -t113 * t139 + t135 * t144;
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t130 = sin(pkin(6));
t163 = t130 * t139;
t94 = t102 * t133 + t137 * t163;
t177 = -t101 * t136 + t132 * t94;
t176 = t101 * t132 + t136 * t94;
t157 = t135 * t138;
t159 = t134 * t139;
t117 = -t131 * t157 - t159;
t164 = t130 * t138;
t175 = -g(1) * t117 - g(3) * t164;
t104 = t120 * t139 - t135 * t140;
t111 = t144 * t130;
t173 = -g(1) * t104 - g(2) * t101 - g(3) * t111;
t165 = t130 * t135;
t162 = t132 * t133;
t161 = t133 * t136;
t158 = t135 * t134;
t156 = t138 * t139;
t154 = t131 * t156;
t95 = t102 * t137 - t133 * t163;
t114 = pkin(2) * t131 * t134 + (-pkin(8) - qJ(3)) * t130;
t128 = pkin(2) * t138 + pkin(1);
t151 = -t114 * t135 + t139 * t128;
t147 = g(1) * t135 - g(2) * t139;
t103 = -t135 * t113 - t139 * t144;
t146 = -t114 * t139 - t135 * t128;
t112 = t120 * t130;
t106 = -t112 * t133 - t131 * t137;
t98 = -t103 * t133 - t137 * t165;
t143 = g(1) * t98 + g(2) * t94 + g(3) * t106;
t121 = pkin(2) * t154;
t118 = -t131 * t158 + t156;
t116 = -t131 * t159 - t157;
t115 = -t154 + t158;
t107 = -t112 * t137 + t131 * t133;
t99 = -t103 * t137 + t133 * t165;
t93 = -t104 * t136 + t132 * t98;
t92 = t104 * t132 + t136 * t98;
t1 = [t147 * MDP(2) + (-g(1) * t116 - g(2) * t118) * MDP(9) + (-g(1) * t115 - g(2) * t117) * MDP(10) + (-g(1) * t146 - g(2) * t151) * MDP(12) + (-g(1) * t101 + g(2) * t104) * MDP(20) + (-g(1) * (-t102 * pkin(3) - pkin(4) * t95 + t101 * pkin(9) - qJ(5) * t94 + t146) - g(2) * (-pkin(3) * t103 + pkin(4) * t99 - pkin(9) * t104 + qJ(5) * t98 + t151)) * MDP(23) + (g(1) * t177 - g(2) * t93) * MDP(29) + (g(1) * t176 - g(2) * t92) * MDP(30) + t174 * (-g(1) * t94 + g(2) * t98) - t179 * (-g(1) * t95 + g(2) * t99) + (-MDP(11) * t130 + MDP(3)) * (g(1) * t139 + g(2) * t135); (g(2) * t115 + t175) * MDP(9) + (g(3) * t130 * t134 + g(1) * t118 - g(2) * t116) * MDP(10) + (-g(2) * t121 + (g(2) * t158 + t175) * pkin(2)) * MDP(12) + (g(1) * t103 - g(2) * t102 + g(3) * t112) * MDP(20) + (-g(1) * (pkin(2) * t117 - t103 * pkin(9)) - g(2) * (-pkin(2) * t158 + pkin(9) * t102 + t121) - g(3) * (pkin(2) * t164 - pkin(9) * t112) + t173 * (pkin(4) * t137 + qJ(5) * t133 + pkin(3))) * MDP(23) + (-g(1) * (-t103 * t136 + t104 * t162) - g(2) * (t101 * t162 + t102 * t136) - g(3) * (t111 * t162 - t112 * t136)) * MDP(29) + (-g(1) * (t103 * t132 + t104 * t161) - g(2) * (t101 * t161 - t102 * t132) - g(3) * (t111 * t161 + t112 * t132)) * MDP(30) + (-t174 * t133 + t137 * t179) * t173; (MDP(12) + MDP(23)) * (-g(3) * t131 - t147 * t130); (-g(1) * (-pkin(4) * t98 + qJ(5) * t99) - g(2) * (-pkin(4) * t94 + qJ(5) * t95) - g(3) * (-pkin(4) * t106 + qJ(5) * t107)) * MDP(23) + (-MDP(29) * t132 - MDP(30) * t136 + t174) * (g(1) * t99 + g(2) * t95 + g(3) * t107) + t179 * t143; -t143 * MDP(23); (-g(1) * t92 - g(2) * t176 - g(3) * (t106 * t136 + t111 * t132)) * MDP(29) + (g(1) * t93 + g(2) * t177 - g(3) * (-t106 * t132 + t111 * t136)) * MDP(30);];
taug  = t1;
