% Calculate Gravitation load on the joints for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:45
% EndTime: 2019-03-08 22:07:47
% DurationCPUTime: 0.81s
% Computational Cost: add. (569->130), mult. (1562->242), div. (0->0), fcn. (2013->16), ass. (0->67)
t129 = sin(pkin(12));
t132 = cos(pkin(12));
t138 = sin(qJ(2));
t134 = cos(pkin(6));
t142 = cos(qJ(2));
t157 = t134 * t142;
t119 = -t129 * t138 + t132 * t157;
t121 = -t129 * t157 - t132 * t138;
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t131 = sin(pkin(6));
t161 = t131 * t142;
t164 = t131 * t132;
t167 = t129 * t131;
t170 = g(1) * (t121 * t133 + t130 * t167) - g(2) * (-t119 * t133 + t130 * t164) + g(3) * (t130 * t134 + t133 * t161);
t169 = g(3) * t131;
t168 = cos(pkin(13));
t136 = sin(qJ(5));
t166 = t130 * t136;
t140 = cos(qJ(5));
t165 = t130 * t140;
t163 = t131 * t133;
t162 = t131 * t138;
t137 = sin(qJ(3));
t160 = t133 * t137;
t141 = cos(qJ(3));
t159 = t133 * t141;
t158 = t134 * t138;
t135 = sin(qJ(6));
t156 = t135 * t140;
t139 = cos(qJ(6));
t155 = t139 * t140;
t154 = t130 * t162;
t128 = sin(pkin(13));
t124 = -t128 * t141 - t137 * t168;
t149 = -t128 * t137 + t141 * t168;
t120 = t129 * t142 + t132 * t158;
t122 = -t129 * t158 + t132 * t142;
t147 = g(1) * t122 + g(2) * t120 + g(3) * t162;
t114 = t124 * t130;
t116 = t124 * t133;
t146 = t114 * t164 - t116 * t119 + t120 * t149;
t145 = -t114 * t167 - t116 * t121 + t122 * t149;
t144 = -t134 * t114 + (-t116 * t142 + t138 * t149) * t131;
t127 = pkin(3) * t141 + pkin(2);
t118 = -t130 * t161 + t133 * t134;
t117 = pkin(3) * t160 + (-pkin(9) - qJ(4)) * t130;
t115 = t149 * t133;
t113 = t149 * t130;
t111 = -t121 * t130 + t129 * t163;
t110 = -t119 * t130 - t132 * t163;
t109 = (t116 * t138 + t142 * t149) * t131;
t108 = (t115 * t138 - t124 * t142) * t131;
t107 = t109 * t140 + t136 * t154;
t106 = t116 * t122 + t121 * t149;
t105 = t115 * t122 - t121 * t124;
t104 = t116 * t120 + t119 * t149;
t103 = t115 * t120 - t119 * t124;
t101 = t134 * t113 + (t115 * t142 + t124 * t138) * t131;
t99 = t118 * t136 + t140 * t144;
t96 = t113 * t167 + t115 * t121 + t122 * t124;
t93 = -t113 * t164 + t115 * t119 + t120 * t124;
t91 = t106 * t140 + t122 * t166;
t90 = t104 * t140 + t120 * t166;
t89 = t111 * t136 + t140 * t145;
t87 = t110 * t136 + t140 * t146;
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * t121 - g(2) * t119 - g(3) * t161) * MDP(3) + (-g(1) * (t121 * t141 - t122 * t160) - g(2) * (t119 * t141 - t120 * t160) - (-t138 * t160 + t141 * t142) * t169) * MDP(10) + (-g(1) * (-t121 * t137 - t122 * t159) - g(2) * (-t119 * t137 - t120 * t159) - (-t137 * t142 - t138 * t159) * t169) * MDP(11) + (-g(1) * (-t117 * t122 + t121 * t127) - g(2) * (-t117 * t120 + t119 * t127) - (-t117 * t138 + t127 * t142) * t169) * MDP(13) + (-g(1) * t91 - g(2) * t90 - g(3) * t107) * MDP(19) + (-g(1) * (-t106 * t136 + t122 * t165) - g(2) * (-t104 * t136 + t120 * t165) - g(3) * (-t109 * t136 + t140 * t154)) * MDP(20) + (-g(1) * (t105 * t135 + t139 * t91) - g(2) * (t103 * t135 + t139 * t90) - g(3) * (t107 * t139 + t108 * t135)) * MDP(26) + (-g(1) * (t105 * t139 - t135 * t91) - g(2) * (t103 * t139 - t135 * t90) - g(3) * (-t107 * t135 + t108 * t139)) * MDP(27) + (-MDP(12) * t130 + MDP(4)) * t147; (t137 * t170 + t147 * t141) * MDP(11) + (-g(1) * (t135 * t145 + t155 * t96) - g(2) * (t135 * t146 + t155 * t93) - g(3) * (t101 * t155 + t135 * t144)) * MDP(26) + (-g(1) * (t139 * t145 - t156 * t96) - g(2) * (t139 * t146 - t156 * t93) - g(3) * (-t101 * t156 + t139 * t144)) * MDP(27) + (-MDP(19) * t140 + MDP(20) * t136) * (g(1) * t96 + g(2) * t93 + g(3) * t101) + (MDP(13) * pkin(3) + MDP(10)) * (t137 * t147 - t141 * t170); (-g(1) * t111 - g(2) * t110 - g(3) * t118) * MDP(13); (g(1) * t89 + g(2) * t87 + g(3) * t99) * MDP(20) + (-MDP(26) * t139 + MDP(27) * t135 - MDP(19)) * (g(1) * (t111 * t140 - t136 * t145) + g(2) * (t110 * t140 - t136 * t146) + g(3) * (t118 * t140 - t136 * t144)); (-g(1) * (-t135 * t89 - t139 * t96) - g(2) * (-t135 * t87 - t139 * t93) - g(3) * (-t101 * t139 - t135 * t99)) * MDP(26) + (-g(1) * (t135 * t96 - t139 * t89) - g(2) * (t135 * t93 - t139 * t87) - g(3) * (t101 * t135 - t139 * t99)) * MDP(27);];
taug  = t1;
