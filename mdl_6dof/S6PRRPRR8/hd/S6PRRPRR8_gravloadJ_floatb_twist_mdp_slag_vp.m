% Calculate Gravitation load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:39
% EndTime: 2019-03-08 22:40:41
% DurationCPUTime: 0.86s
% Computational Cost: add. (482->126), mult. (1344->223), div. (0->0), fcn. (1698->14), ass. (0->64)
t160 = MDP(10) - MDP(13);
t159 = MDP(11) - MDP(14);
t158 = cos(qJ(3));
t118 = sin(pkin(7));
t157 = pkin(9) * t118;
t156 = cos(pkin(6));
t155 = sin(pkin(12));
t123 = sin(qJ(5));
t154 = t118 * t123;
t127 = cos(qJ(5));
t153 = t118 * t127;
t119 = sin(pkin(6));
t121 = cos(pkin(7));
t152 = t119 * t121;
t125 = sin(qJ(2));
t151 = t119 * t125;
t128 = cos(qJ(2));
t150 = t119 * t128;
t124 = sin(qJ(3));
t149 = t121 * t124;
t122 = sin(qJ(6));
t148 = t122 * t123;
t126 = cos(qJ(6));
t147 = t123 * t126;
t146 = t124 * t125;
t145 = t124 * t128;
t120 = cos(pkin(12));
t144 = t118 * t119 * t120;
t143 = t118 * t151;
t142 = t121 * t158;
t141 = t158 * t125;
t140 = t158 * t128;
t139 = t119 * t155;
t138 = t120 * t156;
t137 = t156 * t118;
t136 = t118 * t139;
t135 = t156 * t155;
t109 = -t155 * t125 + t128 * t138;
t110 = t125 * t138 + t155 * t128;
t87 = -t109 * t142 + t110 * t124 + t158 * t144;
t111 = -t120 * t125 - t128 * t135;
t112 = t120 * t128 - t125 * t135;
t89 = -t111 * t142 + t112 * t124 - t158 * t136;
t98 = t119 * t146 - t158 * t137 - t140 * t152;
t133 = g(1) * t89 + g(2) * t87 + g(3) * t98;
t108 = -t118 * t150 + t156 * t121;
t107 = (-t121 * t146 + t140) * t119;
t106 = (t121 * t141 + t145) * t119;
t101 = -t111 * t118 + t121 * t139;
t100 = -t109 * t118 - t120 * t152;
t99 = t124 * t137 + (t121 * t145 + t141) * t119;
t97 = t106 * t123 + t127 * t143;
t96 = t111 * t158 - t112 * t149;
t95 = t111 * t124 + t112 * t142;
t94 = t109 * t158 - t110 * t149;
t93 = t109 * t124 + t110 * t142;
t92 = t108 * t127 + t98 * t123;
t90 = t112 * t158 + (t111 * t121 + t136) * t124;
t88 = t109 * t149 + t110 * t158 - t124 * t144;
t86 = t112 * t153 + t95 * t123;
t85 = t110 * t153 + t93 * t123;
t84 = t101 * t127 + t89 * t123;
t82 = t100 * t127 + t87 * t123;
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * t111 - g(2) * t109 - g(3) * t150) * MDP(3) + (-g(1) * (t111 * pkin(2) + t96 * pkin(3) + t95 * qJ(4) + t112 * t157) - g(2) * (t109 * pkin(2) + t94 * pkin(3) + t93 * qJ(4) + t110 * t157) - g(3) * (t107 * pkin(3) + t106 * qJ(4) + (pkin(2) * t128 + t125 * t157) * t119)) * MDP(15) + (-g(1) * t86 - g(2) * t85 - g(3) * t97) * MDP(21) + (-g(1) * (-t112 * t154 + t95 * t127) - g(2) * (-t110 * t154 + t93 * t127) - g(3) * (t106 * t127 - t123 * t143)) * MDP(22) + (-g(1) * (t96 * t122 + t86 * t126) - g(2) * (t94 * t122 + t85 * t126) - g(3) * (t107 * t122 + t97 * t126)) * MDP(28) + (-g(1) * (-t86 * t122 + t96 * t126) - g(2) * (-t85 * t122 + t94 * t126) - g(3) * (t107 * t126 - t97 * t122)) * MDP(29) + t159 * (g(1) * t95 + g(2) * t93 + g(3) * t106) - t160 * (g(1) * t96 + g(2) * t94 + g(3) * t107) + (-t118 * MDP(12) + MDP(4)) * (g(1) * t112 + g(2) * t110 + g(3) * t151); (-g(1) * (-t89 * pkin(3) + t90 * qJ(4)) - g(2) * (-t87 * pkin(3) + t88 * qJ(4)) - g(3) * (-t98 * pkin(3) + t99 * qJ(4))) * MDP(15) + (-g(1) * (-t89 * t122 + t90 * t147) - g(2) * (-t87 * t122 + t88 * t147) - g(3) * (-t98 * t122 + t99 * t147)) * MDP(28) + (-g(1) * (-t89 * t126 - t90 * t148) - g(2) * (-t87 * t126 - t88 * t148) - g(3) * (-t98 * t126 - t99 * t148)) * MDP(29) + (-t123 * MDP(21) - MDP(22) * t127 + t159) * (g(1) * t90 + g(2) * t88 + g(3) * t99) + t160 * t133; -t133 * MDP(15); (g(1) * t84 + g(2) * t82 + g(3) * t92) * MDP(22) + (-MDP(28) * t126 + MDP(29) * t122 - MDP(21)) * (g(1) * (-t101 * t123 + t89 * t127) + g(2) * (-t100 * t123 + t87 * t127) + g(3) * (-t108 * t123 + t98 * t127)); (-g(1) * (-t84 * t122 + t90 * t126) - g(2) * (-t82 * t122 + t88 * t126) - g(3) * (-t92 * t122 + t99 * t126)) * MDP(28) + (-g(1) * (-t90 * t122 - t84 * t126) - g(2) * (-t88 * t122 - t82 * t126) - g(3) * (-t99 * t122 - t92 * t126)) * MDP(29);];
taug  = t1;
