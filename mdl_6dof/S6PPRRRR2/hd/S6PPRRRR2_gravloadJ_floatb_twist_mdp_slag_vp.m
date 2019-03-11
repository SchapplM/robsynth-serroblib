% Calculate Gravitation load on the joints for
% S6PPRRRR2
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
%   see S6PPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:24
% EndTime: 2019-03-08 19:05:26
% DurationCPUTime: 0.40s
% Computational Cost: add. (527->89), mult. (1360->160), div. (0->0), fcn. (1763->16), ass. (0->50)
t117 = sin(pkin(13));
t118 = sin(pkin(12));
t108 = t118 * t117;
t121 = cos(pkin(13));
t122 = cos(pkin(12));
t115 = t122 * t121;
t124 = cos(pkin(6));
t101 = -t124 * t115 + t108;
t119 = sin(pkin(7));
t120 = sin(pkin(6));
t112 = t120 * t119;
t123 = cos(pkin(7));
t133 = -t101 * t123 - t122 * t112;
t109 = t118 * t121;
t113 = t122 * t117;
t102 = t124 * t109 + t113;
t111 = t120 * t118;
t132 = t102 * t123 - t119 * t111;
t131 = t121 * t123 * t120 + t124 * t119;
t130 = cos(qJ(3));
t93 = qJ(5) + qJ(6);
t91 = sin(t93);
t98 = cos(qJ(4));
t129 = t91 * t98;
t92 = cos(t93);
t128 = t92 * t98;
t94 = sin(qJ(5));
t127 = t94 * t98;
t97 = cos(qJ(5));
t126 = t97 * t98;
t86 = t124 * t113 + t109;
t96 = sin(qJ(3));
t74 = t86 * t130 + t133 * t96;
t114 = t122 * t120;
t81 = t101 * t119 - t123 * t114;
t95 = sin(qJ(4));
t70 = t74 * t98 + t81 * t95;
t87 = -t124 * t108 + t115;
t76 = t87 * t130 - t132 * t96;
t82 = t102 * t119 + t123 * t111;
t72 = t76 * t98 + t82 * t95;
t73 = -t133 * t130 + t86 * t96;
t75 = t132 * t130 + t87 * t96;
t110 = t120 * t117;
t80 = t130 * t110 + t131 * t96;
t85 = -t121 * t112 + t124 * t123;
t78 = t80 * t98 + t85 * t95;
t79 = t96 * t110 - t131 * t130;
t125 = (-g(1) * (-t72 * t91 + t75 * t92) - g(2) * (-t70 * t91 + t73 * t92) - g(3) * (-t78 * t91 + t79 * t92)) * MDP(25) + (-g(1) * (-t72 * t92 - t75 * t91) - g(2) * (-t70 * t92 - t73 * t91) - g(3) * (-t78 * t92 - t79 * t91)) * MDP(26);
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t111 + g(2) * t114 - g(3) * t124) * MDP(2); (g(1) * t76 + g(2) * t74 + g(3) * t80) * MDP(5) + (-g(1) * (-t75 * t126 + t76 * t94) - g(2) * (-t73 * t126 + t74 * t94) - g(3) * (-t79 * t126 + t80 * t94)) * MDP(18) + (-g(1) * (t75 * t127 + t76 * t97) - g(2) * (t73 * t127 + t74 * t97) - g(3) * (t79 * t127 + t80 * t97)) * MDP(19) + (-g(1) * (-t75 * t128 + t76 * t91) - g(2) * (-t73 * t128 + t74 * t91) - g(3) * (-t79 * t128 + t80 * t91)) * MDP(25) + (-g(1) * (t75 * t129 + t76 * t92) - g(2) * (t73 * t129 + t74 * t92) - g(3) * (t79 * t129 + t80 * t92)) * MDP(26) + (t98 * MDP(11) - MDP(12) * t95 + MDP(4)) * (g(1) * t75 + g(2) * t73 + g(3) * t79); (g(1) * t72 + g(2) * t70 + g(3) * t78) * MDP(12) + (-MDP(18) * t97 + MDP(19) * t94 - MDP(25) * t92 + MDP(26) * t91 - MDP(11)) * (g(1) * (-t76 * t95 + t82 * t98) + g(2) * (-t74 * t95 + t81 * t98) + g(3) * (-t80 * t95 + t85 * t98)); (-g(1) * (-t72 * t94 + t75 * t97) - g(2) * (-t70 * t94 + t73 * t97) - g(3) * (-t78 * t94 + t79 * t97)) * MDP(18) + (-g(1) * (-t72 * t97 - t75 * t94) - g(2) * (-t70 * t97 - t73 * t94) - g(3) * (-t78 * t97 - t79 * t94)) * MDP(19) + t125; t125;];
taug  = t1;
