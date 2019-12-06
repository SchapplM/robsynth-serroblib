% Calculate Gravitation load on the joints for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:07
% EndTime: 2019-12-05 17:26:08
% DurationCPUTime: 0.44s
% Computational Cost: add. (363->104), mult. (1019->193), div. (0->0), fcn. (1304->14), ass. (0->60)
t140 = cos(qJ(3));
t139 = cos(pkin(5));
t138 = sin(pkin(11));
t105 = sin(pkin(6));
t110 = sin(qJ(4));
t137 = t105 * t110;
t114 = cos(qJ(4));
t136 = t105 * t114;
t106 = sin(pkin(5));
t108 = cos(pkin(6));
t135 = t106 * t108;
t112 = sin(qJ(2));
t134 = t106 * t112;
t115 = cos(qJ(2));
t133 = t106 * t115;
t111 = sin(qJ(3));
t132 = t108 * t111;
t109 = sin(qJ(5));
t131 = t109 * t114;
t130 = t111 * t112;
t129 = t111 * t115;
t113 = cos(qJ(5));
t128 = t113 * t114;
t107 = cos(pkin(11));
t127 = t105 * t106 * t107;
t126 = t105 * t134;
t125 = t108 * t140;
t124 = t140 * t112;
t123 = t140 * t115;
t122 = t106 * t138;
t121 = t107 * t139;
t120 = t139 * t105;
t119 = t105 * t122;
t118 = t139 * t138;
t100 = t107 * t115 - t112 * t118;
t99 = -t107 * t112 - t115 * t118;
t98 = t112 * t121 + t138 * t115;
t97 = -t138 * t112 + t115 * t121;
t96 = -t105 * t133 + t139 * t108;
t95 = (-t108 * t130 + t123) * t106;
t94 = (t108 * t124 + t129) * t106;
t91 = -t99 * t105 + t108 * t122;
t90 = -t97 * t105 - t107 * t135;
t89 = t111 * t120 + (t108 * t129 + t124) * t106;
t88 = t106 * t130 - t140 * t120 - t123 * t135;
t87 = t110 * t126 + t95 * t114;
t86 = -t100 * t132 + t99 * t140;
t85 = t100 * t125 + t99 * t111;
t84 = -t98 * t132 + t97 * t140;
t83 = t97 * t111 + t98 * t125;
t82 = t96 * t110 + t89 * t114;
t80 = t100 * t140 + (t108 * t99 + t119) * t111;
t79 = t100 * t111 - t140 * t119 - t99 * t125;
t78 = t98 * t140 + (t108 * t97 - t127) * t111;
t77 = t98 * t111 - t97 * t125 + t140 * t127;
t76 = t100 * t137 + t86 * t114;
t75 = t84 * t114 + t98 * t137;
t74 = t91 * t110 + t80 * t114;
t72 = t90 * t110 + t78 * t114;
t1 = [-g(3) * MDP(1); (-g(1) * t99 - g(2) * t97 - g(3) * t133) * MDP(3) + (g(1) * t100 + g(2) * t98 + g(3) * t134) * MDP(4) + (-g(1) * t86 - g(2) * t84 - g(3) * t95) * MDP(10) + (g(1) * t85 + g(2) * t83 + g(3) * t94) * MDP(11) + (-g(1) * t76 - g(2) * t75 - g(3) * t87) * MDP(17) + (-g(1) * (t100 * t136 - t86 * t110) - g(2) * (-t84 * t110 + t98 * t136) - g(3) * (-t95 * t110 + t114 * t126)) * MDP(18) + (-g(1) * (t85 * t109 + t76 * t113) - g(2) * (t83 * t109 + t75 * t113) - g(3) * (t94 * t109 + t87 * t113)) * MDP(24) + (-g(1) * (-t76 * t109 + t85 * t113) - g(2) * (-t75 * t109 + t83 * t113) - g(3) * (-t87 * t109 + t94 * t113)) * MDP(25); (g(1) * t80 + g(2) * t78 + g(3) * t89) * MDP(11) + (-g(1) * (t80 * t109 - t79 * t128) - g(2) * (t78 * t109 - t77 * t128) - g(3) * (t89 * t109 - t88 * t128)) * MDP(24) + (-g(1) * (t80 * t113 + t79 * t131) - g(2) * (t78 * t113 + t77 * t131) - g(3) * (t89 * t113 + t88 * t131)) * MDP(25) + (t114 * MDP(17) - MDP(18) * t110 + MDP(10)) * (g(1) * t79 + g(2) * t77 + g(3) * t88); (g(1) * t74 + g(2) * t72 + g(3) * t82) * MDP(18) + (-MDP(24) * t113 + MDP(25) * t109 - MDP(17)) * (g(1) * (-t80 * t110 + t91 * t114) + g(2) * (-t78 * t110 + t90 * t114) + g(3) * (-t89 * t110 + t96 * t114)); (-g(1) * (-t74 * t109 + t79 * t113) - g(2) * (-t72 * t109 + t77 * t113) - g(3) * (-t82 * t109 + t88 * t113)) * MDP(24) + (-g(1) * (-t79 * t109 - t74 * t113) - g(2) * (-t77 * t109 - t72 * t113) - g(3) * (-t88 * t109 - t82 * t113)) * MDP(25);];
taug = t1;
