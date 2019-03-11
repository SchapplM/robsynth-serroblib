% Calculate Gravitation load on the joints for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:33:59
% EndTime: 2019-03-09 12:34:01
% DurationCPUTime: 0.78s
% Computational Cost: add. (451->123), mult. (793->189), div. (0->0), fcn. (914->12), ass. (0->55)
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t107 = sin(qJ(2));
t108 = sin(qJ(1));
t110 = cos(qJ(2));
t130 = cos(pkin(6));
t137 = cos(qJ(1));
t117 = t130 * t137;
t84 = t108 * t107 - t110 * t117;
t133 = t84 * t109;
t102 = sin(pkin(6));
t124 = t102 * t137;
t85 = t107 * t117 + t108 * t110;
t100 = pkin(11) + qJ(4);
t97 = sin(t100);
t98 = cos(t100);
t77 = -t97 * t124 + t85 * t98;
t145 = t77 * t106 - t133;
t134 = t84 * t106;
t144 = t77 * t109 + t134;
t142 = MDP(10) - MDP(13);
t141 = MDP(21) - MDP(29);
t139 = g(2) * t84;
t138 = g(2) * t85;
t136 = g(3) * t102;
t128 = t102 * t108;
t135 = t137 * pkin(1) + pkin(8) * t128;
t121 = t108 * t130;
t86 = t107 * t137 + t110 * t121;
t132 = t86 * t106;
t131 = t86 * t109;
t129 = t102 * t107;
t127 = t106 * t110;
t126 = t109 * t110;
t101 = sin(pkin(11));
t125 = t101 * t128;
t123 = -t108 * pkin(1) + pkin(8) * t124;
t122 = pkin(5) * t106 + pkin(9) + qJ(3);
t120 = t101 * t124;
t87 = -t107 * t121 + t110 * t137;
t81 = t128 * t97 + t87 * t98;
t73 = -t81 * t106 + t131;
t104 = -qJ(6) - pkin(10);
t103 = cos(pkin(11));
t95 = t103 * pkin(3) + pkin(2);
t96 = t109 * pkin(5) + pkin(4);
t116 = t104 * t97 - t96 * t98 - t95;
t76 = t124 * t98 + t85 * t97;
t80 = -t128 * t98 + t87 * t97;
t82 = t129 * t97 - t130 * t98;
t115 = g(1) * t80 + g(2) * t76 + g(3) * t82;
t112 = -g(1) * t86 + t110 * t136 - t139;
t83 = t129 * t98 + t130 * t97;
t74 = t81 * t109 + t132;
t1 = [(g(1) * t108 - g(2) * t137) * MDP(2) + (g(1) * t137 + g(2) * t108) * MDP(3) + (g(1) * t85 - g(2) * t87) * MDP(9) + (-g(1) * (-t85 * t103 + t120) - g(2) * (t87 * t103 + t125)) * MDP(11) + (-g(1) * (t85 * t101 + t103 * t124) - g(2) * (-t87 * t101 + t103 * t128)) * MDP(12) + (-g(1) * (-t85 * pkin(2) - t84 * qJ(3) + t123) - g(2) * (t87 * pkin(2) + t86 * qJ(3) + t135)) * MDP(14) + (g(1) * t77 - g(2) * t81) * MDP(20) + (g(1) * t144 - g(2) * t74) * MDP(27) + (-g(1) * t145 - g(2) * t73) * MDP(28) + (-g(1) * (pkin(3) * t120 + t104 * t76 - t122 * t84 - t77 * t96 - t85 * t95 + t123) - g(2) * (pkin(3) * t125 - t80 * t104 + t122 * t86 + t81 * t96 + t87 * t95 + t135)) * MDP(30) + t141 * (-g(1) * t76 + g(2) * t80) - t142 * (g(1) * t84 - g(2) * t86); (-g(1) * (-t86 * pkin(2) + t87 * qJ(3)) - g(2) * (-t84 * pkin(2) + t85 * qJ(3)) - (pkin(2) * t110 + qJ(3) * t107) * t136) * MDP(14) + (-g(1) * (t87 * t106 - t131 * t98) - g(2) * (t85 * t106 - t133 * t98) - (t106 * t107 + t126 * t98) * t136) * MDP(27) + (-g(1) * (t87 * t109 + t132 * t98) - g(2) * (t85 * t109 + t134 * t98) - (t107 * t109 - t127 * t98) * t136) * MDP(28) + (-g(1) * (t116 * t86 + t122 * t87) - t122 * t138 - t116 * t139 - (t107 * t122 - t110 * t116) * t136) * MDP(30) + t142 * (g(1) * t87 + g(3) * t129 + t138) + (-t103 * MDP(11) + t101 * MDP(12) - t98 * MDP(20) + t141 * t97 - MDP(9)) * t112; (MDP(14) + MDP(30)) * t112; (-g(1) * (-t81 * t104 - t80 * t96) - g(2) * (-t77 * t104 - t76 * t96) - g(3) * (-t83 * t104 - t82 * t96)) * MDP(30) + t141 * (g(1) * t81 + g(2) * t77 + g(3) * t83) + (MDP(27) * t109 - MDP(28) * t106 + MDP(20)) * t115; (g(1) * t74 + g(2) * t144 - g(3) * (t102 * t127 - t83 * t109)) * MDP(28) + (pkin(5) * MDP(30) + MDP(27)) * (g(2) * t145 - g(3) * (-t102 * t126 - t83 * t106) - g(1) * t73); -t115 * MDP(30);];
taug  = t1;
