% Calculate Gravitation load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:13
% EndTime: 2019-03-08 21:02:16
% DurationCPUTime: 0.78s
% Computational Cost: add. (352->107), mult. (602->178), div. (0->0), fcn. (691->14), ass. (0->56)
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t101 = sin(pkin(6));
t105 = sin(qJ(2));
t126 = t101 * t105;
t129 = cos(pkin(6));
t140 = -t104 * t126 + t129 * t106;
t100 = sin(pkin(10));
t125 = t101 * t106;
t107 = cos(qJ(2));
t119 = t100 * t129;
t128 = cos(pkin(10));
t84 = -t105 * t119 + t107 * t128;
t139 = t100 * t125 - t104 * t84;
t138 = g(3) * t101;
t97 = pkin(12) + qJ(6);
t93 = sin(t97);
t98 = qJ(3) + pkin(11);
t96 = cos(t98);
t137 = t93 * t96;
t95 = cos(t97);
t136 = t95 * t96;
t99 = sin(pkin(12));
t135 = t96 * t99;
t103 = -qJ(4) - pkin(8);
t114 = t129 * t128;
t81 = t100 * t105 - t107 * t114;
t82 = t100 * t107 + t105 * t114;
t92 = pkin(3) * t106 + pkin(2);
t134 = -t82 * t103 - t81 * t92;
t83 = t105 * t128 + t107 * t119;
t133 = -t84 * t103 - t83 * t92;
t102 = cos(pkin(12));
t132 = t102 * t96;
t130 = t107 * t96;
t127 = t100 * t101;
t124 = t101 * t107;
t123 = t103 * t105;
t122 = MDP(13) + MDP(17);
t118 = t101 * t128;
t116 = t139 * pkin(3);
t94 = sin(t98);
t115 = pkin(4) * t96 + qJ(5) * t94;
t113 = t140 * pkin(3);
t71 = t118 * t96 + t82 * t94;
t73 = -t127 * t96 + t84 * t94;
t77 = t126 * t94 - t129 * t96;
t112 = g(1) * t73 + g(2) * t71 + g(3) * t77;
t111 = -t82 * t104 - t106 * t118;
t109 = -g(1) * t83 - g(2) * t81 + g(3) * t124;
t108 = t111 * pkin(3);
t85 = t92 * t124;
t78 = t126 * t96 + t129 * t94;
t74 = t127 * t94 + t84 * t96;
t72 = -t118 * t94 + t82 * t96;
t1 = [(-MDP(1) - t122) * g(3); (-g(1) * t133 - g(2) * t134 - g(3) * (-t101 * t123 + t85)) * MDP(13) + (-g(1) * (-t132 * t83 + t84 * t99) - g(2) * (-t132 * t81 + t82 * t99) - (t102 * t130 + t105 * t99) * t138) * MDP(14) + (-g(1) * (t102 * t84 + t135 * t83) - g(2) * (t102 * t82 + t135 * t81) - (t102 * t105 - t130 * t99) * t138) * MDP(15) + (-g(1) * (-t115 * t83 + t133) - g(2) * (-t115 * t81 + t134) - g(3) * t85 - (t107 * t115 - t123) * t138) * MDP(17) + (-g(1) * (-t136 * t83 + t84 * t93) - g(2) * (-t136 * t81 + t82 * t93) - (t105 * t93 + t130 * t95) * t138) * MDP(23) + (-g(1) * (t137 * t83 + t84 * t95) - g(2) * (t137 * t81 + t82 * t95) - (t105 * t95 - t130 * t93) * t138) * MDP(24) + (MDP(4) - MDP(12)) * (g(1) * t84 + g(2) * t82 + g(3) * t126) + (-t106 * MDP(10) + t104 * MDP(11) - t94 * MDP(16) - MDP(3)) * t109; (-g(1) * t139 - g(2) * t111 - g(3) * t140) * MDP(10) + (-g(1) * (-t104 * t127 - t106 * t84) - g(2) * (t104 * t118 - t82 * t106) - g(3) * (-t104 * t129 - t105 * t125)) * MDP(11) + (-g(1) * t116 - g(2) * t108 - g(3) * t113) * MDP(13) + (-g(1) * t74 - g(2) * t72 - g(3) * t78) * MDP(16) + (-g(1) * (-pkin(4) * t73 + qJ(5) * t74 + t116) - g(2) * (-t71 * pkin(4) + t72 * qJ(5) + t108) - g(3) * (-pkin(4) * t77 + qJ(5) * t78 + t113)) * MDP(17) + (MDP(14) * t102 - MDP(15) * t99 + MDP(23) * t95 - MDP(24) * t93) * t112; t122 * t109; -t112 * MDP(17); (-g(1) * (-t74 * t93 + t83 * t95) - g(2) * (-t72 * t93 + t81 * t95) - g(3) * (-t124 * t95 - t78 * t93)) * MDP(23) + (-g(1) * (-t74 * t95 - t83 * t93) - g(2) * (-t72 * t95 - t81 * t93) - g(3) * (t124 * t93 - t78 * t95)) * MDP(24);];
taug  = t1;
