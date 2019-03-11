% Calculate Gravitation load on the joints for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:21
% EndTime: 2019-03-10 01:03:22
% DurationCPUTime: 0.38s
% Computational Cost: add. (521->80), mult. (492->107), div. (0->0), fcn. (454->10), ass. (0->44)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t136 = pkin(5) * t105 + qJ(6) * t102;
t135 = MDP(24) - MDP(33);
t132 = MDP(30) + MDP(32);
t131 = MDP(31) - MDP(34);
t101 = qJ(2) + qJ(3);
t98 = qJ(4) + t101;
t94 = sin(t98);
t95 = cos(t98);
t134 = t95 * pkin(4) + t94 * pkin(10);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t116 = g(1) * t107 + g(2) * t104;
t133 = t116 * t94;
t96 = sin(t101);
t130 = pkin(3) * t96;
t129 = pkin(10) * t95;
t128 = g(3) * t94;
t123 = t104 * t102;
t122 = t104 * t105;
t121 = t105 * t107;
t120 = t107 * t102;
t119 = t136 * t95 + t134;
t97 = cos(t101);
t91 = pkin(3) * t97;
t118 = t91 + t119;
t106 = cos(qJ(2));
t99 = t106 * pkin(2);
t114 = t91 + t99 + pkin(1) + t134;
t112 = t135 * (t116 * t95 + t128) + (-t131 * t102 + t132 * t105 + MDP(23)) * (-g(3) * t95 + t133);
t79 = t95 * t123 + t121;
t81 = t95 * t120 - t122;
t64 = g(1) * t81 + g(2) * t79 + t102 * t128;
t110 = (-g(3) * t97 + t116 * t96) * MDP(16) + (g(3) * t96 + t116 * t97) * MDP(17) + t112;
t108 = (pkin(4) + t136) * t133;
t103 = sin(qJ(2));
t100 = -pkin(9) - pkin(8) - pkin(7);
t88 = t107 * t129;
t86 = t104 * t129;
t84 = -pkin(2) * t103 - t130;
t82 = t95 * t121 + t123;
t80 = t95 * t122 - t120;
t1 = [t116 * MDP(3) + (-g(1) * (-t80 * pkin(5) - t79 * qJ(6)) - g(2) * (t82 * pkin(5) + t81 * qJ(6)) + (g(1) * t100 - g(2) * t114) * t107 + (g(1) * t114 + g(2) * t100) * t104) * MDP(35) + t132 * (g(1) * t80 - g(2) * t82) - t131 * (g(1) * t79 - g(2) * t81) + (-t103 * MDP(10) + MDP(16) * t97 - MDP(17) * t96 + t95 * MDP(23) + t106 * MDP(9) - t135 * t94 + MDP(2)) * (g(1) * t104 - g(2) * t107); (-g(3) * t106 + t116 * t103) * MDP(9) + (g(3) * t103 + t116 * t106) * MDP(10) + (-g(1) * (t107 * t84 + t88) - g(2) * (t104 * t84 + t86) - g(3) * (t99 + t118) + t108) * MDP(35) + t110; (-g(1) * (-t107 * t130 + t88) - g(2) * (-t104 * t130 + t86) - g(3) * t118 + t108) * MDP(35) + t110; (-g(1) * t88 - g(2) * t86 - g(3) * t119 + t108) * MDP(35) + t112; (-g(1) * (-pkin(5) * t81 + qJ(6) * t82) - g(2) * (-pkin(5) * t79 + qJ(6) * t80) - (-pkin(5) * t102 + qJ(6) * t105) * t128) * MDP(35) + t132 * t64 + t131 * (g(1) * t82 + g(2) * t80 + t105 * t128); -t64 * MDP(35);];
taug  = t1;
