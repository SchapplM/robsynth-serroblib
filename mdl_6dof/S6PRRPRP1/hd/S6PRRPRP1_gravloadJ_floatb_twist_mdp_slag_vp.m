% Calculate Gravitation load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:35
% EndTime: 2019-03-08 21:26:37
% DurationCPUTime: 0.56s
% Computational Cost: add. (304->101), mult. (561->165), div. (0->0), fcn. (635->12), ass. (0->59)
t95 = sin(pkin(6));
t134 = g(3) * t95;
t104 = cos(qJ(2));
t101 = sin(qJ(2));
t119 = cos(pkin(6));
t113 = t101 * t119;
t94 = sin(pkin(10));
t96 = cos(pkin(10));
t79 = t94 * t104 + t96 * t113;
t99 = sin(qJ(5));
t133 = t79 * t99;
t81 = t104 * t96 - t94 * t113;
t132 = t81 * t99;
t93 = qJ(3) + pkin(11);
t92 = cos(t93);
t131 = t92 * t99;
t130 = t94 * t95;
t129 = t95 * t96;
t112 = t104 * t119;
t78 = t101 * t94 - t96 * t112;
t103 = cos(qJ(3));
t90 = pkin(3) * t103 + pkin(2);
t98 = -qJ(4) - pkin(8);
t128 = -t78 * t90 - t79 * t98;
t80 = t96 * t101 + t94 * t112;
t127 = -t80 * t90 - t81 * t98;
t100 = sin(qJ(3));
t126 = t100 * t81;
t125 = t100 * t95;
t124 = t101 * t95;
t102 = cos(qJ(5));
t123 = t102 * t92;
t122 = t103 * t95;
t121 = t104 * t95;
t120 = t104 * t99;
t118 = t102 * t104;
t117 = MDP(13) + MDP(22);
t116 = t94 * t122;
t115 = t96 * t122;
t114 = t100 * t124;
t111 = t119 * t103;
t89 = pkin(5) * t102 + pkin(4);
t91 = sin(t93);
t97 = -qJ(6) - pkin(9);
t110 = t89 * t92 - t91 * t97;
t109 = -t100 * t79 - t115;
t68 = t92 * t129 + t79 * t91;
t70 = -t92 * t130 + t81 * t91;
t74 = -t119 * t92 + t91 * t124;
t108 = g(1) * t70 + g(2) * t68 + g(3) * t74;
t107 = g(1) * t81 + g(2) * t79 + g(3) * t124;
t106 = -g(1) * t80 - g(2) * t78 + g(3) * t121;
t88 = pkin(3) * t111;
t84 = pkin(3) * t116;
t82 = t90 * t121;
t75 = t119 * t91 + t92 * t124;
t71 = t91 * t130 + t81 * t92;
t69 = -t91 * t129 + t79 * t92;
t1 = [(-MDP(1) - t117) * g(3); (-g(1) * t127 - g(2) * t128 - g(3) * (-t98 * t124 + t82)) * MDP(13) + (-g(1) * (-t80 * t123 + t132) - g(2) * (-t78 * t123 + t133) - (t101 * t99 + t92 * t118) * t134) * MDP(19) + (-g(1) * (t102 * t81 + t80 * t131) - g(2) * (t102 * t79 + t78 * t131) - (t101 * t102 - t92 * t120) * t134) * MDP(20) + (-g(1) * (pkin(5) * t132 - t110 * t80 + t127) - g(2) * (pkin(5) * t133 - t110 * t78 + t128) - g(3) * t82 - (t110 * t104 + (pkin(5) * t99 - t98) * t101) * t134) * MDP(22) + (MDP(4) - MDP(12)) * t107 + (-MDP(10) * t103 + MDP(11) * t100 - t91 * MDP(21) - MDP(3)) * t106; (-g(1) * (t116 - t126) - g(2) * t109 - g(3) * (t111 - t114)) * MDP(10) + (-g(1) * (-t103 * t81 - t94 * t125) - g(2) * (-t103 * t79 + t96 * t125) - g(3) * (-t119 * t100 - t101 * t122)) * MDP(11) + (-g(1) * t84 - g(3) * t88 + (g(2) * t115 + t107 * t100) * pkin(3)) * MDP(13) + (-g(1) * t71 - g(2) * t69 - g(3) * t75) * MDP(21) + (-g(1) * (-pkin(3) * t126 - t70 * t89 - t71 * t97 + t84) - g(2) * (t109 * pkin(3) - t68 * t89 - t69 * t97) - g(3) * (-pkin(3) * t114 - t74 * t89 - t75 * t97 + t88)) * MDP(22) + (t102 * MDP(19) - t99 * MDP(20)) * t108; t117 * t106; (-g(1) * (-t102 * t71 - t80 * t99) - g(2) * (-t102 * t69 - t78 * t99) - g(3) * (-t75 * t102 + t95 * t120)) * MDP(20) + (pkin(5) * MDP(22) + MDP(19)) * (-g(1) * (t102 * t80 - t71 * t99) - g(2) * (t102 * t78 - t69 * t99) - g(3) * (-t95 * t118 - t75 * t99)); -t108 * MDP(22);];
taug  = t1;
