% Calculate Gravitation load on the joints for
% S6PRRPRR6
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
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:45
% EndTime: 2019-03-08 22:27:47
% DurationCPUTime: 0.74s
% Computational Cost: add. (582->137), mult. (1430->243), div. (0->0), fcn. (1808->16), ass. (0->64)
t121 = sin(pkin(7));
t127 = sin(qJ(2));
t129 = cos(qJ(2));
t124 = cos(pkin(12));
t158 = cos(pkin(6));
t142 = t124 * t158;
t156 = sin(pkin(12));
t132 = t156 * t127 - t129 * t142;
t122 = sin(pkin(6));
t149 = t122 * t124;
t157 = cos(pkin(7));
t163 = t121 * t149 + t132 * t157;
t138 = t158 * t156;
t133 = t124 * t127 + t129 * t138;
t143 = t122 * t156;
t162 = -t121 * t143 + t133 * t157;
t161 = MDP(11) - MDP(14);
t160 = cos(qJ(3));
t159 = pkin(9) * t121;
t119 = pkin(13) + qJ(5);
t117 = sin(t119);
t155 = t117 * t121;
t118 = cos(t119);
t154 = t118 * t121;
t125 = sin(qJ(6));
t153 = t118 * t125;
t128 = cos(qJ(6));
t152 = t118 * t128;
t120 = sin(pkin(13));
t151 = t120 * t121;
t123 = cos(pkin(13));
t150 = t121 * t123;
t148 = t122 * t127;
t147 = t122 * t129;
t145 = t121 * t148;
t144 = t121 * t158;
t126 = sin(qJ(3));
t141 = t126 * t157;
t139 = t157 * t160;
t110 = t127 * t142 + t156 * t129;
t90 = t110 * t126 + t163 * t160;
t111 = t124 * t129 - t127 * t138;
t92 = t111 * t126 + t162 * t160;
t99 = t126 * t148 - t139 * t147 - t160 * t144;
t136 = g(1) * t92 + g(2) * t90 + g(3) * t99;
t109 = -t121 * t147 + t158 * t157;
t108 = (-t127 * t141 + t160 * t129) * t122;
t107 = (t126 * t129 + t127 * t139) * t122;
t102 = t133 * t121 + t157 * t143;
t101 = t132 * t121 - t157 * t149;
t100 = t126 * t144 + (t160 * t127 + t129 * t141) * t122;
t98 = t108 * t118 + t117 * t145;
t97 = -t111 * t141 - t133 * t160;
t96 = t111 * t139 - t133 * t126;
t95 = -t110 * t141 - t132 * t160;
t94 = t110 * t139 - t132 * t126;
t93 = t111 * t160 - t162 * t126;
t91 = t110 * t160 - t163 * t126;
t89 = t100 * t118 + t109 * t117;
t87 = t111 * t155 + t118 * t97;
t86 = t110 * t155 + t118 * t95;
t85 = t102 * t117 + t118 * t93;
t83 = t101 * t117 + t118 * t91;
t1 = [(-MDP(1) - MDP(15)) * g(3); (g(1) * t133 + g(2) * t132 - g(3) * t147) * MDP(3) + (g(1) * t111 + g(2) * t110 + g(3) * t148) * MDP(4) + (-g(1) * t97 - g(2) * t95 - g(3) * t108) * MDP(10) + (-g(1) * (t111 * t151 + t123 * t97) - g(2) * (t110 * t151 + t123 * t95) - g(3) * (t108 * t123 + t120 * t145)) * MDP(12) + (-g(1) * (t111 * t150 - t120 * t97) - g(2) * (t110 * t150 - t120 * t95) - g(3) * (-t108 * t120 + t123 * t145)) * MDP(13) + (-g(1) * (-t133 * pkin(2) + t97 * pkin(3) + t96 * qJ(4) + t111 * t159) - g(2) * (-t132 * pkin(2) + t95 * pkin(3) + t94 * qJ(4) + t110 * t159) - g(3) * (t108 * pkin(3) + t107 * qJ(4) + (pkin(2) * t129 + t127 * t159) * t122)) * MDP(15) + (-g(1) * t87 - g(2) * t86 - g(3) * t98) * MDP(21) + (-g(1) * (t111 * t154 - t117 * t97) - g(2) * (t110 * t154 - t117 * t95) - g(3) * (-t108 * t117 + t118 * t145)) * MDP(22) + (-g(1) * (t125 * t96 + t128 * t87) - g(2) * (t125 * t94 + t128 * t86) - g(3) * (t107 * t125 + t128 * t98)) * MDP(28) + (-g(1) * (-t125 * t87 + t128 * t96) - g(2) * (-t125 * t86 + t128 * t94) - g(3) * (t107 * t128 - t125 * t98)) * MDP(29) + t161 * (g(1) * t96 + g(2) * t94 + g(3) * t107); (-g(1) * (-pkin(3) * t92 + qJ(4) * t93) - g(2) * (-pkin(3) * t90 + qJ(4) * t91) - g(3) * (-pkin(3) * t99 + qJ(4) * t100)) * MDP(15) + (-g(1) * (t125 * t93 - t92 * t152) - g(2) * (t125 * t91 - t90 * t152) - g(3) * (t100 * t125 - t99 * t152)) * MDP(28) + (-g(1) * (t128 * t93 + t92 * t153) - g(2) * (t128 * t91 + t90 * t153) - g(3) * (t100 * t128 + t99 * t153)) * MDP(29) + t161 * (g(1) * t93 + g(2) * t91 + g(3) * t100) + (MDP(12) * t123 - MDP(13) * t120 + t118 * MDP(21) - MDP(22) * t117 + MDP(10)) * t136; -t136 * MDP(15); (g(1) * t85 + g(2) * t83 + g(3) * t89) * MDP(22) + (-MDP(28) * t128 + MDP(29) * t125 - MDP(21)) * (g(1) * (t102 * t118 - t117 * t93) + g(2) * (t101 * t118 - t117 * t91) + g(3) * (-t100 * t117 + t109 * t118)); (-g(1) * (-t125 * t85 + t128 * t92) - g(2) * (-t125 * t83 + t128 * t90) - g(3) * (-t89 * t125 + t128 * t99)) * MDP(28) + (-g(1) * (-t125 * t92 - t128 * t85) - g(2) * (-t125 * t90 - t128 * t83) - g(3) * (-t125 * t99 - t128 * t89)) * MDP(29);];
taug  = t1;
