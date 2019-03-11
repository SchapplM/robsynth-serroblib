% Calculate Gravitation load on the joints for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:20
% EndTime: 2019-03-08 23:28:22
% DurationCPUTime: 0.90s
% Computational Cost: add. (534->138), mult. (1352->244), div. (0->0), fcn. (1703->16), ass. (0->65)
t167 = MDP(11) - MDP(19);
t132 = sin(qJ(2));
t135 = cos(qJ(2));
t127 = cos(pkin(12));
t162 = cos(pkin(6));
t145 = t127 * t162;
t160 = sin(pkin(12));
t113 = t132 * t145 + t135 * t160;
t141 = t162 * t160;
t115 = t127 * t135 - t132 * t141;
t126 = sin(pkin(6));
t153 = t126 * t132;
t166 = g(1) * t115 + g(2) * t113 + g(3) * t153;
t165 = cos(qJ(3));
t161 = cos(pkin(7));
t124 = qJ(4) + pkin(13);
t122 = sin(t124);
t125 = sin(pkin(7));
t159 = t122 * t125;
t123 = cos(t124);
t129 = sin(qJ(6));
t158 = t123 * t129;
t133 = cos(qJ(6));
t157 = t123 * t133;
t130 = sin(qJ(4));
t156 = t125 * t130;
t134 = cos(qJ(4));
t155 = t125 * t134;
t154 = t126 * t127;
t152 = t126 * t135;
t150 = t125 * t154;
t149 = t125 * t153;
t147 = t125 * t162;
t146 = t126 * t160;
t131 = sin(qJ(3));
t144 = t131 * t161;
t143 = t125 * t146;
t142 = t161 * t165;
t101 = t131 * t153 - t142 * t152 - t165 * t147;
t112 = -t132 * t160 + t135 * t145;
t92 = -t112 * t142 + t113 * t131 + t165 * t150;
t114 = -t127 * t132 - t135 * t141;
t94 = -t114 * t142 + t115 * t131 - t165 * t143;
t139 = g(1) * t94 + g(2) * t92 + g(3) * t101;
t128 = -qJ(5) - pkin(10);
t121 = pkin(4) * t134 + pkin(3);
t111 = -t125 * t152 + t161 * t162;
t110 = (-t132 * t144 + t165 * t135) * t126;
t109 = (t131 * t135 + t132 * t142) * t126;
t104 = -t114 * t125 + t146 * t161;
t103 = -t112 * t125 - t154 * t161;
t102 = t131 * t147 + (t165 * t132 + t135 * t144) * t126;
t100 = t110 * t123 + t122 * t149;
t99 = t114 * t165 - t115 * t144;
t98 = t114 * t131 + t115 * t142;
t97 = t112 * t165 - t113 * t144;
t96 = t112 * t131 + t113 * t142;
t95 = t115 * t165 + (t114 * t161 + t143) * t131;
t93 = t112 * t144 + t113 * t165 - t131 * t150;
t91 = t102 * t123 + t111 * t122;
t89 = t115 * t159 + t123 * t99;
t88 = t113 * t159 + t123 * t97;
t87 = t104 * t122 + t123 * t95;
t85 = t103 * t122 + t123 * t93;
t1 = [(-MDP(1) - MDP(20)) * g(3); (-g(1) * t114 - g(2) * t112 - g(3) * t152) * MDP(3) + t166 * MDP(4) + (-g(1) * t99 - g(2) * t97 - g(3) * t110) * MDP(10) + (-g(1) * (t115 * t156 + t134 * t99) - g(2) * (t113 * t156 + t134 * t97) - g(3) * (t110 * t134 + t130 * t149)) * MDP(17) + (-g(1) * (t115 * t155 - t130 * t99) - g(2) * (t113 * t155 - t130 * t97) - g(3) * (-t110 * t130 + t134 * t149)) * MDP(18) + (-g(1) * (pkin(2) * t114 + t121 * t99 - t128 * t98) - g(2) * (pkin(2) * t112 + t121 * t97 - t128 * t96) - g(3) * (pkin(2) * t152 - t109 * t128 + t110 * t121) - t166 * t125 * (pkin(4) * t130 + pkin(9))) * MDP(20) + (-g(1) * (t129 * t98 + t133 * t89) - g(2) * (t129 * t96 + t133 * t88) - g(3) * (t100 * t133 + t109 * t129)) * MDP(26) + (-g(1) * (-t129 * t89 + t133 * t98) - g(2) * (-t129 * t88 + t133 * t96) - g(3) * (-t100 * t129 + t109 * t133)) * MDP(27) + t167 * (g(1) * t98 + g(2) * t96 + g(3) * t109); (-g(1) * (-t121 * t94 - t128 * t95) - g(2) * (-t121 * t92 - t128 * t93) - g(3) * (-t101 * t121 - t102 * t128)) * MDP(20) + (-g(1) * (t129 * t95 - t157 * t94) - g(2) * (t129 * t93 - t157 * t92) - g(3) * (-t101 * t157 + t102 * t129)) * MDP(26) + (-g(1) * (t133 * t95 + t158 * t94) - g(2) * (t133 * t93 + t158 * t92) - g(3) * (t101 * t158 + t102 * t133)) * MDP(27) + t167 * (g(1) * t95 + g(2) * t93 + g(3) * t102) + (MDP(17) * t134 - MDP(18) * t130 + MDP(10)) * t139; (-g(1) * (-t104 * t130 - t134 * t95) - g(2) * (-t103 * t130 - t134 * t93) - g(3) * (-t102 * t134 - t111 * t130)) * MDP(18) + (-MDP(26) * t133 + MDP(27) * t129) * (g(1) * (t104 * t123 - t122 * t95) + g(2) * (t103 * t123 - t122 * t93) + g(3) * (-t102 * t122 + t111 * t123)) + (pkin(4) * MDP(20) + MDP(17)) * (-g(1) * (t104 * t134 - t130 * t95) - g(2) * (t103 * t134 - t130 * t93) - g(3) * (-t102 * t130 + t111 * t134)); -t139 * MDP(20); (-g(1) * (-t129 * t87 + t133 * t94) - g(2) * (-t129 * t85 + t133 * t92) - g(3) * (t101 * t133 - t129 * t91)) * MDP(26) + (-g(1) * (-t129 * t94 - t133 * t87) - g(2) * (-t129 * t92 - t133 * t85) - g(3) * (-t101 * t129 - t133 * t91)) * MDP(27);];
taug  = t1;
