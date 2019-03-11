% Calculate Gravitation load on the joints for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:52
% EndTime: 2019-03-09 09:47:54
% DurationCPUTime: 0.55s
% Computational Cost: add. (342->99), mult. (393->133), div. (0->0), fcn. (361->10), ass. (0->51)
t101 = -qJ(5) - pkin(8);
t100 = qJ(2) + pkin(9);
t94 = sin(t100);
t129 = t94 * t101;
t106 = cos(qJ(4));
t91 = pkin(4) * t106 + pkin(3);
t96 = cos(t100);
t138 = -t96 * t91 + t129;
t137 = MDP(20) + MDP(23);
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t118 = g(1) * t108 + g(2) * t105;
t113 = -g(3) * t96 + t118 * t94;
t135 = g(3) * t94;
t136 = t118 * t96 + t135;
t131 = t105 * t96;
t130 = t108 * t96;
t102 = -qJ(3) - pkin(7);
t128 = t102 * t108;
t103 = sin(qJ(4));
t127 = t103 * t108;
t126 = t105 * t103;
t125 = t105 * t106;
t124 = t106 * t108;
t123 = MDP(21) + MDP(25);
t122 = t96 * t127;
t107 = cos(qJ(2));
t97 = t107 * pkin(2);
t92 = t97 + pkin(1);
t121 = -t105 * t102 + t108 * t92;
t120 = t97 - t138;
t99 = qJ(4) + pkin(10);
t93 = sin(t99);
t95 = cos(t99);
t119 = pkin(5) * t95 + qJ(6) * t93;
t83 = g(1) * t105 - g(2) * t108;
t104 = sin(qJ(2));
t117 = -pkin(2) * t104 - t101 * t96;
t76 = t96 * t126 + t124;
t71 = t108 * t95 + t93 * t131;
t73 = -t105 * t95 + t93 * t130;
t115 = g(1) * t73 + g(2) * t71 + t93 * t135;
t114 = pkin(4) * t126 - t108 * t129 + t91 * t130 + t121;
t110 = pkin(4) * t127 - t128 + (-t92 + t138) * t105;
t89 = pkin(4) * t125;
t79 = t96 * t124 + t126;
t78 = -t122 + t125;
t77 = -t96 * t125 + t127;
t74 = t105 * t93 + t95 * t130;
t72 = -t108 * t93 + t95 * t131;
t1 = [(-g(1) * (-t105 * t92 - t128) - g(2) * t121) * MDP(12) + (-g(1) * t77 - g(2) * t79) * MDP(18) + (-g(1) * t76 - g(2) * t78) * MDP(19) + (-g(1) * t110 - g(2) * t114) * MDP(21) + (g(1) * t72 - g(2) * t74) * MDP(22) + (g(1) * t71 - g(2) * t73) * MDP(24) + (-g(1) * (-t72 * pkin(5) - t71 * qJ(6) + t110) - g(2) * (t74 * pkin(5) + t73 * qJ(6) + t114)) * MDP(25) + (MDP(3) - MDP(11)) * t118 + (-t104 * MDP(10) + t107 * MDP(9) + t137 * t94 + MDP(2)) * t83; (g(3) * t104 + t118 * t107) * MDP(10) + (-g(3) * t120 + t118 * (t91 * t94 - t117)) * MDP(21) + (-g(3) * (t119 * t96 + t120) + t118 * (-(-t119 - t91) * t94 - t117)) * MDP(25) - t137 * t136 + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t107 + t118 * t104) + (t106 * MDP(18) - t103 * MDP(19) + t95 * MDP(22) + t93 * MDP(24)) * t113; (-MDP(12) - t123) * t83; (-g(1) * t78 + g(2) * t76 + t103 * t135) * MDP(18) + (g(1) * t79 - g(2) * t77 + t106 * t135) * MDP(19) + (-g(1) * t89 + (g(2) * t124 + t136 * t103) * pkin(4)) * MDP(21) + t115 * MDP(22) + (-g(1) * t74 - g(2) * t72 - t95 * t135) * MDP(24) + (-g(1) * (-pkin(4) * t122 - t73 * pkin(5) + t74 * qJ(6) + t89) - g(2) * (-t76 * pkin(4) - t71 * pkin(5) + t72 * qJ(6)) - (-pkin(4) * t103 - pkin(5) * t93 + qJ(6) * t95) * t135) * MDP(25); -t123 * t113; -t115 * MDP(25);];
taug  = t1;
