% Calculate Gravitation load on the joints for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:40
% EndTime: 2019-03-09 11:55:41
% DurationCPUTime: 0.46s
% Computational Cost: add. (409->87), mult. (432->120), div. (0->0), fcn. (422->10), ass. (0->50)
t140 = MDP(25) + MDP(27);
t139 = MDP(26) - MDP(29);
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t88 = g(1) * t108 + g(2) * t105;
t103 = sin(qJ(4));
t100 = qJ(2) + pkin(10);
t94 = sin(t100);
t135 = g(3) * t94;
t106 = cos(qJ(4));
t121 = t106 * t108;
t124 = t103 * t105;
t95 = cos(t100);
t82 = t95 * t124 + t121;
t122 = t105 * t106;
t123 = t103 * t108;
t84 = -t95 * t123 + t122;
t138 = -g(1) * t84 + g(2) * t82 + t103 * t135;
t101 = qJ(4) + qJ(5);
t96 = sin(t101);
t131 = t94 * t96;
t97 = cos(t101);
t130 = t94 * t97;
t129 = t105 * t96;
t128 = t105 * t97;
t127 = t108 * t96;
t126 = t108 * t97;
t109 = -pkin(9) - pkin(8);
t125 = t94 * t109;
t102 = -qJ(3) - pkin(7);
t119 = pkin(4) * t103 - t102;
t78 = t95 * t129 + t126;
t80 = t95 * t127 - t128;
t72 = g(1) * t80 + g(2) * t78 + g(3) * t131;
t79 = t95 * t128 - t127;
t81 = t95 * t126 + t129;
t118 = t140 * t72 + t139 * (g(1) * t81 + g(2) * t79 + g(3) * t130);
t87 = g(1) * t105 - g(2) * t108;
t92 = pkin(4) * t106 + pkin(3);
t116 = t95 * t92 - t125;
t115 = pkin(5) * t97 + qJ(6) * t96 + t92;
t110 = -g(1) * (-t80 * pkin(5) + qJ(6) * t81) - g(2) * (-t78 * pkin(5) + qJ(6) * t79) - g(3) * (-pkin(5) * t131 + qJ(6) * t130);
t107 = cos(qJ(2));
t104 = sin(qJ(2));
t98 = t107 * pkin(2);
t93 = t98 + pkin(1);
t89 = t108 * t93;
t85 = t95 * t121 + t124;
t83 = -t95 * t122 + t123;
t1 = [(-g(1) * (-t108 * t102 - t105 * t93) - g(2) * (-t105 * t102 + t89)) * MDP(12) + (-g(1) * t83 - g(2) * t85) * MDP(18) + (-g(1) * t82 - g(2) * t84) * MDP(19) + (-g(1) * (-t79 * pkin(5) - t78 * qJ(6)) - g(2) * (t81 * pkin(5) + t80 * qJ(6) + t89) + (-g(1) * t119 - g(2) * t116) * t108 + (-g(1) * (-t116 - t93) - g(2) * t119) * t105) * MDP(30) + (MDP(3) - MDP(11)) * t88 + t140 * (g(1) * t79 - g(2) * t81) - t139 * (g(1) * t78 - g(2) * t80) + (-t104 * MDP(10) + t94 * MDP(28) + t107 * MDP(9) + MDP(2)) * t87; (g(3) * t104 + t88 * t107) * MDP(10) + (-t88 * t95 - t135) * MDP(28) + (-g(3) * (t115 * t95 - t125 + t98) + t88 * (pkin(2) * t104 + t109 * t95 + t115 * t94)) * MDP(30) + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t107 + t88 * t104) + (t106 * MDP(18) - t103 * MDP(19) - t139 * t96 + t140 * t97) * (-g(3) * t95 + t88 * t94); (-MDP(12) - MDP(30)) * t87; t138 * MDP(18) + (g(1) * t85 - g(2) * t83 + t106 * t135) * MDP(19) + (pkin(4) * t138 + t110) * MDP(30) + t118; t110 * MDP(30) + t118; -t72 * MDP(30);];
taug  = t1;
