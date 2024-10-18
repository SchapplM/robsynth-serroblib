% Calculate Gravitation load on the joints for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR15_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR15_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:26
% DurationCPUTime: 0.31s
% Computational Cost: add. (579->99), mult. (524->162), div. (0->0), fcn. (582->22), ass. (0->79)
t138 = sin(qJ(1));
t172 = -t138 / 0.2e1;
t133 = sin(pkin(5));
t171 = g(3) * t133;
t132 = sin(pkin(6));
t170 = t132 * t133;
t135 = cos(pkin(5));
t169 = t132 * t135;
t134 = cos(pkin(6));
t136 = sin(qJ(5));
t168 = t134 * t136;
t139 = cos(qJ(5));
t167 = t134 * t139;
t166 = t138 * t136;
t137 = sin(qJ(2));
t165 = t138 * t137;
t140 = cos(qJ(2));
t164 = t138 * t140;
t163 = t139 * t138;
t141 = cos(qJ(1));
t162 = t141 * t136;
t161 = t141 * t137;
t160 = t141 * t139;
t159 = t141 * t140;
t131 = qJ(2) + qJ(3);
t126 = pkin(5) + t131;
t122 = qJ(4) + t126;
t158 = sin(t122) / 0.2e1;
t157 = t138 * t170;
t156 = t141 * t170;
t155 = t135 * t162;
t154 = t135 * t166;
t153 = t135 * t163;
t152 = t135 * t160;
t127 = pkin(5) - t131;
t100 = t134 * t155 + t163;
t101 = t134 * t152 - t166;
t102 = t134 * t166 - t152;
t103 = t134 * t163 + t155;
t104 = t134 * t162 + t153;
t105 = -t134 * t160 + t154;
t116 = cos(t122) / 0.2e1;
t150 = -qJ(4) + t127;
t117 = cos(t150);
t130 = qJ(4) + t131;
t124 = sin(t130);
t125 = cos(t130);
t146 = sin(t150);
t106 = t158 - t146 / 0.2e1;
t147 = t138 * t106 - t141 * t125;
t148 = t141 * t106 + t138 * t125;
t107 = t117 / 0.2e1 + t116;
t94 = -t141 * t107 + t138 * t124;
t95 = -t138 * t107 - t141 * t124;
t98 = t134 * t154 - t160;
t99 = -t134 * t153 - t162;
t151 = (-g(1) * t95 + g(2) * t94 - g(3) * (t158 + t146 / 0.2e1)) * MDP(23) + (-g(1) * t147 + g(2) * t148 - g(3) * (t116 - t117 / 0.2e1)) * MDP(24) + (-g(1) * (-t104 * t125 + t98 * t124) - g(2) * (-t100 * t124 - t102 * t125) - (-t124 * t168 + t125 * t139) * t171) * MDP(30) + (-g(1) * (t105 * t125 - t99 * t124) - g(2) * (-t101 * t124 - t103 * t125) - (-t124 * t167 - t125 * t136) * t171) * MDP(31);
t118 = sin(t126);
t119 = sin(t127);
t120 = cos(t126);
t121 = cos(t127);
t112 = t118 - t119;
t129 = cos(t131);
t144 = t112 * t172 + t141 * t129;
t145 = t138 * t129 + t141 * t112 / 0.2e1;
t113 = t121 + t120;
t128 = sin(t131);
t96 = t138 * t128 - t141 * t113 / 0.2e1;
t97 = t113 * t172 - t141 * t128;
t149 = (-g(1) * t97 + g(2) * t96 - g(3) * (t118 / 0.2e1 + t119 / 0.2e1)) * MDP(16) + (g(1) * t144 + g(2) * t145 - g(3) * (t120 / 0.2e1 - t121 / 0.2e1)) * MDP(17) + t151;
t143 = -t104 * t124 - t98 * t125 + t136 * t157;
t142 = -t101 * t125 + t103 * t124 + t139 * t156;
t111 = -t135 * t165 + t159;
t110 = -t135 * t164 - t161;
t109 = -t135 * t161 - t164;
t108 = -t135 * t159 + t165;
t93 = -t100 * t125 + t102 * t124 + t136 * t156;
t92 = t105 * t124 + t99 * t125 + t139 * t157;
t1 = [(g(1) * t138 - g(2) * t141) * MDP(2) + (g(1) * t141 + g(2) * t138) * MDP(3) + (-g(1) * t109 - g(2) * t111) * MDP(9) + (-g(1) * t108 - g(2) * t110) * MDP(10) + (g(1) * t145 - g(2) * t144) * MDP(16) + (-g(1) * t96 - g(2) * t97) * MDP(17) + (g(1) * t148 + g(2) * t147) * MDP(23) + (-g(1) * t94 - g(2) * t95) * MDP(24) + (-g(1) * t93 - g(2) * t143) * MDP(30) + (-g(1) * t142 - g(2) * t92) * MDP(31); (-g(1) * t110 + g(2) * t108 - t140 * t171) * MDP(9) + (g(1) * t111 - g(2) * t109 + t137 * t171) * MDP(10) + t149; t149; t151; (-g(1) * t92 + g(2) * t142 - g(3) * (t139 * t169 + (-t124 * t136 + t125 * t167) * t133)) * MDP(30) + (g(1) * t143 - g(2) * t93 - g(3) * (-t136 * t169 + (-t124 * t139 - t125 * t168) * t133)) * MDP(31);];
taug = t1;
