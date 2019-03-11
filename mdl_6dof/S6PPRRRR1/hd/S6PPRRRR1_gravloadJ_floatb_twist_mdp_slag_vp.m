% Calculate Gravitation load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:30
% EndTime: 2019-03-08 19:01:31
% DurationCPUTime: 0.35s
% Computational Cost: add. (479->77), mult. (1138->138), div. (0->0), fcn. (1463->16), ass. (0->48)
t120 = sin(pkin(13));
t121 = sin(pkin(12));
t110 = t121 * t120;
t124 = cos(pkin(13));
t125 = cos(pkin(12));
t117 = t125 * t124;
t127 = cos(pkin(6));
t103 = -t127 * t117 + t110;
t122 = sin(pkin(7));
t123 = sin(pkin(6));
t114 = t123 * t122;
t126 = cos(pkin(7));
t133 = t103 * t126 + t125 * t114;
t111 = t121 * t124;
t115 = t125 * t120;
t104 = t127 * t111 + t115;
t113 = t123 * t121;
t132 = t104 * t126 - t122 * t113;
t131 = t124 * t126 * t123 + t127 * t122;
t130 = cos(qJ(3));
t95 = qJ(4) + qJ(5);
t94 = cos(t95);
t96 = sin(qJ(6));
t129 = t94 * t96;
t99 = cos(qJ(6));
t128 = t94 * t99;
t88 = t127 * t115 + t111;
t98 = sin(qJ(3));
t78 = t88 * t130 - t133 * t98;
t116 = t125 * t123;
t83 = t103 * t122 - t126 * t116;
t93 = sin(t95);
t72 = t78 * t94 + t83 * t93;
t89 = -t127 * t110 + t117;
t80 = t89 * t130 - t132 * t98;
t84 = t104 * t122 + t126 * t113;
t74 = t80 * t94 + t84 * t93;
t112 = t123 * t120;
t82 = t130 * t112 + t131 * t98;
t87 = -t124 * t114 + t127 * t126;
t76 = t82 * t94 + t87 * t93;
t119 = (g(1) * t74 + g(2) * t72 + g(3) * t76) * MDP(19) + (-t99 * MDP(25) + t96 * MDP(26) - MDP(18)) * (g(1) * (-t80 * t93 + t84 * t94) + g(2) * (-t78 * t93 + t83 * t94) + g(3) * (-t82 * t93 + t87 * t94));
t100 = cos(qJ(4));
t97 = sin(qJ(4));
t81 = t98 * t112 - t131 * t130;
t79 = t132 * t130 + t89 * t98;
t77 = t133 * t130 + t88 * t98;
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t113 + g(2) * t116 - g(3) * t127) * MDP(2); (g(1) * t80 + g(2) * t78 + g(3) * t82) * MDP(5) + (-g(1) * (-t79 * t128 + t80 * t96) - g(2) * (-t77 * t128 + t78 * t96) - g(3) * (-t81 * t128 + t82 * t96)) * MDP(25) + (-g(1) * (t79 * t129 + t80 * t99) - g(2) * (t77 * t129 + t78 * t99) - g(3) * (t81 * t129 + t82 * t99)) * MDP(26) + (MDP(11) * t100 - MDP(12) * t97 + t94 * MDP(18) - MDP(19) * t93 + MDP(4)) * (g(1) * t79 + g(2) * t77 + g(3) * t81); (-g(1) * (t84 * t100 - t80 * t97) - g(2) * (t83 * t100 - t78 * t97) - g(3) * (t87 * t100 - t82 * t97)) * MDP(11) + (-g(1) * (-t80 * t100 - t84 * t97) - g(2) * (-t78 * t100 - t83 * t97) - g(3) * (-t82 * t100 - t87 * t97)) * MDP(12) + t119; t119; (-g(1) * (-t74 * t96 + t79 * t99) - g(2) * (-t72 * t96 + t77 * t99) - g(3) * (-t76 * t96 + t81 * t99)) * MDP(25) + (-g(1) * (-t74 * t99 - t79 * t96) - g(2) * (-t72 * t99 - t77 * t96) - g(3) * (-t76 * t99 - t81 * t96)) * MDP(26);];
taug  = t1;
