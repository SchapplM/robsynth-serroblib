% Calculate Gravitation load on the joints for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:29
% EndTime: 2019-03-09 04:48:30
% DurationCPUTime: 0.52s
% Computational Cost: add. (234->98), mult. (375->134), div. (0->0), fcn. (346->8), ass. (0->49)
t108 = sin(qJ(3));
t110 = cos(qJ(4));
t112 = cos(qJ(1));
t127 = t110 * t112;
t107 = sin(qJ(4));
t109 = sin(qJ(1));
t130 = t109 * t107;
t78 = -t108 * t130 + t127;
t129 = t109 * t110;
t131 = t108 * t112;
t80 = t107 * t131 + t129;
t140 = MDP(13) - MDP(21) - MDP(24);
t111 = cos(qJ(3));
t137 = g(2) * t112;
t138 = g(1) * t109;
t87 = -t137 + t138;
t77 = -g(3) * t108 + t111 * t87;
t139 = pkin(4) * t107;
t135 = g(3) * t111;
t134 = t80 * pkin(4);
t133 = t107 * t112;
t132 = t108 * t109;
t128 = t109 * t111;
t126 = t111 * t112;
t125 = t112 * pkin(1) + t109 * qJ(2);
t124 = MDP(22) + MDP(26);
t123 = t107 * t135;
t101 = t112 * qJ(2);
t106 = -qJ(5) - pkin(8);
t97 = pkin(4) * t110 + pkin(3);
t122 = t106 * t126 + t97 * t131 + t101;
t105 = qJ(4) + pkin(9);
t98 = sin(t105);
t99 = cos(t105);
t119 = pkin(5) * t99 + qJ(6) * t98;
t118 = t78 * pkin(4);
t117 = t119 + t97;
t116 = pkin(4) * t133 + t112 * pkin(7) + t106 * t128 + t97 * t132 + t125;
t115 = (-pkin(1) - pkin(7) - t139) * t138;
t72 = -t112 * t99 + t132 * t98;
t74 = t109 * t99 + t131 * t98;
t114 = g(1) * t72 - g(2) * t74 + t135 * t98;
t90 = t106 * t131;
t84 = t97 * t128;
t81 = t108 * t127 - t130;
t79 = t108 * t129 + t133;
t75 = -t109 * t98 + t131 * t99;
t73 = t112 * t98 + t132 * t99;
t1 = [(-g(1) * (-t109 * pkin(1) + t101) - g(2) * t125) * MDP(6) + (-g(1) * t81 - g(2) * t79) * MDP(19) + (g(1) * t80 - g(2) * t78) * MDP(20) + (-g(1) * t122 - g(2) * t116 - t115) * MDP(22) + (-g(1) * t75 - g(2) * t73) * MDP(23) + (-g(1) * t74 - g(2) * t72) * MDP(25) + (-g(1) * (t75 * pkin(5) + t74 * qJ(6) + t122) - g(2) * (t73 * pkin(5) + t72 * qJ(6) + t116) - t115) * MDP(26) + (MDP(2) - MDP(4)) * t87 + (-MDP(12) * t108 - t140 * t111 + MDP(3) - MDP(5)) * (g(1) * t112 + g(2) * t109); (-MDP(6) - t124) * t87; (-g(1) * (-t106 * t132 + t84) - g(2) * (-t126 * t97 + t90) - g(3) * (-t106 * t111 - t108 * t97)) * MDP(22) + (-g(1) * t84 - g(2) * t90 + (g(3) * t117 + t106 * t138) * t108 + (g(3) * t106 + t117 * t137 - t119 * t138) * t111) * MDP(26) + t140 * (g(1) * t132 - g(2) * t131 + t135) + (-t110 * MDP(19) + t107 * MDP(20) - t99 * MDP(23) - t98 * MDP(25) - MDP(12)) * t77; (-g(1) * t78 - g(2) * t80 + t123) * MDP(19) + (g(1) * t79 - g(2) * t81 + t110 * t135) * MDP(20) + (pkin(4) * t123 - g(1) * t118 - g(2) * t134) * MDP(22) + t114 * MDP(23) + (-g(1) * t73 + g(2) * t75 - t135 * t99) * MDP(25) + (-g(1) * (-pkin(5) * t72 + qJ(6) * t73 + t118) - g(2) * (pkin(5) * t74 - qJ(6) * t75 + t134) - (-pkin(5) * t98 + qJ(6) * t99 - t139) * t135) * MDP(26); t124 * t77; -t114 * MDP(26);];
taug  = t1;
