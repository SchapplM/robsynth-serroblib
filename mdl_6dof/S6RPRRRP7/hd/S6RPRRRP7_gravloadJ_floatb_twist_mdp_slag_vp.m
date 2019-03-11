% Calculate Gravitation load on the joints for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:00
% EndTime: 2019-03-09 06:21:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (413->82), mult. (426->113), div. (0->0), fcn. (418->10), ass. (0->46)
t137 = MDP(14) - MDP(30);
t136 = MDP(27) + MDP(29);
t135 = MDP(28) - MDP(31);
t103 = sin(qJ(1));
t105 = cos(qJ(1));
t87 = g(1) * t105 + g(2) * t103;
t102 = sin(qJ(4));
t97 = pkin(10) + qJ(3);
t92 = sin(t97);
t131 = g(3) * t92;
t104 = cos(qJ(4));
t117 = t104 * t105;
t120 = t102 * t103;
t93 = cos(t97);
t81 = t93 * t120 + t117;
t118 = t103 * t104;
t119 = t102 * t105;
t83 = -t93 * t119 + t118;
t134 = -g(1) * t83 + g(2) * t81 + t102 * t131;
t98 = qJ(4) + qJ(5);
t94 = sin(t98);
t127 = t92 * t94;
t95 = cos(t98);
t126 = t92 * t95;
t125 = t103 * t94;
t124 = t103 * t95;
t123 = t105 * t94;
t122 = t105 * t95;
t106 = -pkin(9) - pkin(8);
t121 = t92 * t106;
t115 = pkin(4) * t102 + pkin(7) + qJ(2);
t76 = t93 * t125 + t122;
t78 = t93 * t123 - t124;
t69 = g(1) * t78 + g(2) * t76 + g(3) * t127;
t77 = t93 * t124 - t123;
t79 = t93 * t122 + t125;
t114 = t136 * t69 + t135 * (g(1) * t79 + g(2) * t77 + g(3) * t126);
t86 = g(1) * t103 - g(2) * t105;
t91 = pkin(4) * t104 + pkin(3);
t112 = pkin(5) * t95 + qJ(6) * t94 + t91;
t100 = cos(pkin(10));
t111 = -pkin(2) * t100 - t93 * t91 - pkin(1) + t121;
t107 = -g(1) * (-t78 * pkin(5) + qJ(6) * t79) - g(2) * (-t76 * pkin(5) + qJ(6) * t77) - g(3) * (-pkin(5) * t127 + qJ(6) * t126);
t84 = t93 * t117 + t120;
t82 = -t93 * t118 + t119;
t1 = [(-g(1) * (-pkin(1) * t103 + qJ(2) * t105) - g(2) * (pkin(1) * t105 + qJ(2) * t103)) * MDP(7) + (-g(1) * t82 - g(2) * t84) * MDP(20) + (-g(1) * t81 - g(2) * t83) * MDP(21) + (-g(1) * (-t77 * pkin(5) - t76 * qJ(6)) - g(2) * (t79 * pkin(5) + t78 * qJ(6)) + (-g(1) * t115 + g(2) * t111) * t105 + (-g(1) * t111 - g(2) * t115) * t103) * MDP(32) + (MDP(3) - MDP(6)) * t87 + t136 * (g(1) * t77 - g(2) * t79) - t135 * (g(1) * t76 - g(2) * t78) + (-t137 * t92 + t93 * MDP(13) + MDP(4) * t100 - MDP(5) * sin(pkin(10)) + MDP(2)) * t86; (-MDP(32) - MDP(7)) * t86; (-g(3) * (t112 * t93 - t121) + t87 * (t106 * t93 + t112 * t92)) * MDP(32) + t137 * (t87 * t93 + t131) + (t104 * MDP(20) - t102 * MDP(21) - t135 * t94 + t136 * t95 + MDP(13)) * (-g(3) * t93 + t87 * t92); t134 * MDP(20) + (g(1) * t84 - g(2) * t82 + t104 * t131) * MDP(21) + (t134 * pkin(4) + t107) * MDP(32) + t114; t107 * MDP(32) + t114; -t69 * MDP(32);];
taug  = t1;
