% Calculate Gravitation load on the joints for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:42
% EndTime: 2019-03-08 19:11:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (955->119), mult. (2752->219), div. (0->0), fcn. (3625->18), ass. (0->70)
t128 = cos(pkin(13));
t130 = cos(pkin(7));
t123 = sin(pkin(14));
t124 = sin(pkin(13));
t127 = cos(pkin(14));
t131 = cos(pkin(6));
t158 = t128 * t131;
t150 = -t123 * t124 + t127 * t158;
t126 = sin(pkin(6));
t164 = sin(pkin(7));
t153 = t126 * t164;
t166 = -t128 * t153 + t150 * t130;
t165 = cos(qJ(4));
t163 = t124 * t131;
t125 = sin(pkin(8));
t133 = sin(qJ(5));
t162 = t125 * t133;
t137 = cos(qJ(5));
t161 = t125 * t137;
t160 = t126 * t130;
t159 = t127 * t130;
t129 = cos(pkin(8));
t134 = sin(qJ(4));
t157 = t129 * t134;
t132 = sin(qJ(6));
t156 = t132 * t137;
t136 = cos(qJ(6));
t155 = t136 * t137;
t154 = t129 * t165;
t152 = t131 * t164;
t149 = -t123 * t128 - t127 * t163;
t146 = -t127 * t153 + t131 * t130;
t144 = t146 * t125;
t143 = -t128 * t160 - t150 * t164;
t142 = t124 * t153 + t149 * t130;
t141 = t124 * t160 - t149 * t164;
t140 = t143 * t125;
t139 = t141 * t125;
t138 = cos(qJ(3));
t135 = sin(qJ(3));
t121 = -t123 * t163 + t127 * t128;
t120 = t123 * t158 + t124 * t127;
t118 = t135 * t152 + (t123 * t138 + t135 * t159) * t126;
t117 = t138 * t152 + (-t123 * t135 + t138 * t159) * t126;
t113 = -t117 * t125 + t146 * t129;
t112 = t121 * t138 + t142 * t135;
t111 = -t121 * t135 + t142 * t138;
t110 = t120 * t138 + t166 * t135;
t109 = -t120 * t135 + t166 * t138;
t106 = t117 * t165 - t118 * t157;
t105 = t117 * t134 + t118 * t154;
t104 = -t111 * t125 + t141 * t129;
t103 = -t109 * t125 + t143 * t129;
t102 = t118 * t165 + (t117 * t129 + t144) * t134;
t101 = -t117 * t154 + t118 * t134 - t165 * t144;
t100 = t106 * t137 + t118 * t162;
t99 = t111 * t165 - t112 * t157;
t98 = t111 * t134 + t112 * t154;
t97 = t109 * t165 - t110 * t157;
t96 = t109 * t134 + t110 * t154;
t95 = t102 * t137 + t113 * t133;
t93 = t112 * t165 + (t111 * t129 + t139) * t134;
t92 = -t111 * t154 + t112 * t134 - t165 * t139;
t91 = t110 * t165 + (t109 * t129 + t140) * t134;
t90 = -t109 * t154 + t110 * t134 - t165 * t140;
t89 = t112 * t162 + t137 * t99;
t88 = t110 * t162 + t137 * t97;
t87 = t104 * t133 + t137 * t93;
t85 = t103 * t133 + t137 * t91;
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(3) * t131 + (-g(1) * t124 + g(2) * t128) * t126) * MDP(2); (-g(1) * t111 - g(2) * t109 - g(3) * t117) * MDP(4) + (g(1) * t112 + g(2) * t110 + g(3) * t118) * MDP(5) + (-g(1) * t99 - g(2) * t97 - g(3) * t106) * MDP(11) + (g(1) * t98 + g(2) * t96 + g(3) * t105) * MDP(12) + (-g(1) * t89 - g(2) * t88 - g(3) * t100) * MDP(18) + (-g(1) * (t112 * t161 - t133 * t99) - g(2) * (t110 * t161 - t133 * t97) - g(3) * (-t106 * t133 + t118 * t161)) * MDP(19) + (-g(1) * (t132 * t98 + t136 * t89) - g(2) * (t132 * t96 + t136 * t88) - g(3) * (t100 * t136 + t105 * t132)) * MDP(25) + (-g(1) * (-t89 * t132 + t136 * t98) - g(2) * (-t132 * t88 + t136 * t96) - g(3) * (-t100 * t132 + t105 * t136)) * MDP(26); (g(1) * t93 + g(2) * t91 + g(3) * t102) * MDP(12) + (-g(1) * (t132 * t93 - t92 * t155) - g(2) * (t132 * t91 - t90 * t155) - g(3) * (-t101 * t155 + t102 * t132)) * MDP(25) + (-g(1) * (t136 * t93 + t92 * t156) - g(2) * (t136 * t91 + t90 * t156) - g(3) * (t101 * t156 + t102 * t136)) * MDP(26) + (t137 * MDP(18) - MDP(19) * t133 + MDP(11)) * (g(1) * t92 + g(2) * t90 + g(3) * t101); (g(1) * t87 + g(2) * t85 + g(3) * t95) * MDP(19) + (-MDP(25) * t136 + MDP(26) * t132 - MDP(18)) * (g(1) * (t104 * t137 - t133 * t93) + g(2) * (t103 * t137 - t133 * t91) + g(3) * (-t102 * t133 + t113 * t137)); (-g(1) * (-t132 * t87 + t136 * t92) - g(2) * (-t132 * t85 + t136 * t90) - g(3) * (t101 * t136 - t132 * t95)) * MDP(25) + (-g(1) * (-t132 * t92 - t136 * t87) - g(2) * (-t132 * t90 - t136 * t85) - g(3) * (-t101 * t132 - t136 * t95)) * MDP(26);];
taug  = t1;
