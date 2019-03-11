% Calculate Gravitation load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:02
% EndTime: 2019-03-09 01:00:04
% DurationCPUTime: 0.65s
% Computational Cost: add. (640->129), mult. (1523->234), div. (0->0), fcn. (1950->16), ass. (0->62)
t119 = sin(pkin(7));
t125 = sin(qJ(2));
t128 = cos(qJ(2));
t121 = cos(pkin(13));
t156 = cos(pkin(6));
t140 = t121 * t156;
t154 = sin(pkin(13));
t131 = t154 * t125 - t128 * t140;
t120 = sin(pkin(6));
t147 = t120 * t121;
t155 = cos(pkin(7));
t159 = t119 * t147 + t131 * t155;
t135 = t156 * t154;
t132 = t121 * t125 + t128 * t135;
t141 = t120 * t154;
t158 = -t119 * t141 + t132 * t155;
t157 = cos(qJ(3));
t118 = qJ(4) + qJ(5);
t116 = sin(t118);
t153 = t116 * t119;
t117 = cos(t118);
t152 = t117 * t119;
t122 = sin(qJ(6));
t151 = t117 * t122;
t126 = cos(qJ(6));
t150 = t117 * t126;
t123 = sin(qJ(4));
t149 = t119 * t123;
t127 = cos(qJ(4));
t148 = t119 * t127;
t146 = t120 * t125;
t145 = t120 * t128;
t144 = t119 * t146;
t142 = t119 * t156;
t124 = sin(qJ(3));
t139 = t124 * t155;
t102 = t124 * t142 + (t157 * t125 + t128 * t139) * t120;
t103 = t131 * t119 - t155 * t147;
t104 = t132 * t119 + t155 * t141;
t109 = -t119 * t145 + t156 * t155;
t110 = t125 * t140 + t154 * t128;
t93 = t110 * t157 - t159 * t124;
t85 = t103 * t116 + t93 * t117;
t111 = t121 * t128 - t125 * t135;
t95 = t111 * t157 - t158 * t124;
t87 = t104 * t116 + t95 * t117;
t91 = t102 * t117 + t109 * t116;
t138 = (g(1) * t87 + g(2) * t85 + g(3) * t91) * MDP(25) + (-t126 * MDP(31) + t122 * MDP(32) - MDP(24)) * (g(1) * (t104 * t117 - t95 * t116) + g(2) * (t103 * t117 - t93 * t116) + g(3) * (-t102 * t116 + t109 * t117));
t136 = t155 * t157;
t108 = (-t125 * t139 + t157 * t128) * t120;
t107 = (t124 * t128 + t125 * t136) * t120;
t101 = t124 * t146 - t136 * t145 - t157 * t142;
t100 = t108 * t117 + t116 * t144;
t99 = -t111 * t139 - t132 * t157;
t98 = t111 * t136 - t132 * t124;
t97 = -t110 * t139 - t131 * t157;
t96 = t110 * t136 - t131 * t124;
t94 = t111 * t124 + t158 * t157;
t92 = t110 * t124 + t159 * t157;
t89 = t111 * t153 + t99 * t117;
t88 = t110 * t153 + t97 * t117;
t1 = [-g(3) * MDP(1); (g(1) * t132 + g(2) * t131 - g(3) * t145) * MDP(3) + (g(1) * t111 + g(2) * t110 + g(3) * t146) * MDP(4) + (-g(1) * t99 - g(2) * t97 - g(3) * t108) * MDP(10) + (g(1) * t98 + g(2) * t96 + g(3) * t107) * MDP(11) + (-g(1) * (t111 * t149 + t99 * t127) - g(2) * (t110 * t149 + t97 * t127) - g(3) * (t108 * t127 + t123 * t144)) * MDP(17) + (-g(1) * (t111 * t148 - t99 * t123) - g(2) * (t110 * t148 - t97 * t123) - g(3) * (-t108 * t123 + t127 * t144)) * MDP(18) + (-g(1) * t89 - g(2) * t88 - g(3) * t100) * MDP(24) + (-g(1) * (t111 * t152 - t99 * t116) - g(2) * (t110 * t152 - t97 * t116) - g(3) * (-t108 * t116 + t117 * t144)) * MDP(25) + (-g(1) * (t98 * t122 + t89 * t126) - g(2) * (t96 * t122 + t88 * t126) - g(3) * (t100 * t126 + t107 * t122)) * MDP(31) + (-g(1) * (-t89 * t122 + t98 * t126) - g(2) * (-t88 * t122 + t96 * t126) - g(3) * (-t100 * t122 + t107 * t126)) * MDP(32); (g(1) * t95 + g(2) * t93 + g(3) * t102) * MDP(11) + (-g(1) * (t95 * t122 - t94 * t150) - g(2) * (t93 * t122 - t92 * t150) - g(3) * (-t101 * t150 + t102 * t122)) * MDP(31) + (-g(1) * (t95 * t126 + t94 * t151) - g(2) * (t93 * t126 + t92 * t151) - g(3) * (t101 * t151 + t102 * t126)) * MDP(32) + (MDP(17) * t127 - MDP(18) * t123 + t117 * MDP(24) - MDP(25) * t116 + MDP(10)) * (g(1) * t94 + g(2) * t92 + g(3) * t101); (-g(1) * (t104 * t127 - t95 * t123) - g(2) * (t103 * t127 - t93 * t123) - g(3) * (-t102 * t123 + t109 * t127)) * MDP(17) + (-g(1) * (-t104 * t123 - t95 * t127) - g(2) * (-t103 * t123 - t93 * t127) - g(3) * (-t102 * t127 - t109 * t123)) * MDP(18) + t138; t138; (-g(1) * (-t87 * t122 + t94 * t126) - g(2) * (-t85 * t122 + t92 * t126) - g(3) * (t101 * t126 - t91 * t122)) * MDP(31) + (-g(1) * (-t94 * t122 - t87 * t126) - g(2) * (-t92 * t122 - t85 * t126) - g(3) * (-t101 * t122 - t91 * t126)) * MDP(32);];
taug  = t1;
