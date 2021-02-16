% Calculate Gravitation load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:19
% EndTime: 2021-01-15 23:41:23
% DurationCPUTime: 0.82s
% Computational Cost: add. (303->106), mult. (577->177), div. (0->0), fcn. (655->12), ass. (0->59)
t105 = sin(qJ(5));
t158 = MDP(28) * t105 - MDP(18);
t107 = sin(qJ(2));
t108 = sin(qJ(1));
t140 = cos(pkin(5));
t144 = cos(qJ(2));
t124 = t140 * t144;
t145 = cos(qJ(1));
t112 = -t107 * t145 - t108 * t124;
t103 = sin(pkin(5));
t142 = g(3) * t103;
t91 = t107 * t108 - t124 * t145;
t114 = -g(1) * t112 + g(2) * t91 - t144 * t142;
t153 = MDP(10) - MDP(20);
t102 = qJ(3) + pkin(10);
t101 = cos(t102);
t130 = t107 * t140;
t131 = t108 * t144;
t90 = -t130 * t145 - t131;
t100 = sin(t102);
t133 = t103 * t145;
t94 = t100 * t133;
t118 = t90 * t101 + t94;
t136 = t103 * t108;
t126 = t145 * t144;
t89 = t108 * t130 - t126;
t121 = t100 * t136 - t101 * t89;
t137 = t103 * t107;
t84 = t100 * t140 + t101 * t137;
t152 = g(1) * t121 - g(2) * t118 + g(3) * t84;
t109 = cos(qJ(5));
t151 = -t105 * t112 + t109 * t121;
t106 = sin(qJ(3));
t143 = pkin(3) * t106;
t139 = t101 * t105;
t138 = t101 * t109;
t135 = t106 * t103;
t110 = cos(qJ(3));
t134 = t110 * t103;
t132 = t105 * t144;
t129 = t140 * t106;
t128 = t140 * t110;
t123 = -t105 * t91 + t109 * t94 + t138 * t90;
t122 = t106 * t126;
t79 = t100 * t89 + t101 * t136;
t120 = t108 * t135 - t110 * t89;
t104 = qJ(4) + pkin(8);
t99 = pkin(3) * t110 + pkin(2);
t119 = t104 * t107 + t144 * t99;
t117 = t90 * t100 - t101 * t133;
t116 = t107 * t129 + t134;
t113 = g(3) * (-t107 * t135 + t128);
t88 = t104 * t144 - t99 * t107;
t87 = pkin(1) + t119;
t86 = t116 * pkin(3);
t82 = t104 * t130 + t124 * t99;
t81 = -t90 * t106 + t110 * t133;
t80 = t104 * t124 - t99 * t130 + t103 * (pkin(7) + t143);
t1 = [(g(1) * t108 - g(2) * t145) * MDP(2) + (g(1) * t145 + g(2) * t108) * MDP(3) + (-g(1) * t90 + g(2) * t89) * MDP(9) + (-g(1) * ((-t107 * t128 + t135) * t145 - t110 * t131) - g(2) * t120) * MDP(16) + (-g(1) * t81 - g(2) * (t89 * t106 + t108 * t134)) * MDP(17) + (g(1) * t117 - g(2) * t79) * MDP(19) + (-g(1) * (-t87 * t108 + t145 * t80) - g(2) * (t80 * t108 + t145 * t87)) * MDP(21) + (-g(1) * t123 - g(2) * t151) * MDP(27) + (t109 * MDP(28) - t153) * (g(1) * t91 + g(2) * t112) + t158 * (g(1) * t118 + g(2) * t121); (-g(1) * (-t82 * t108 + t145 * t88) - g(2) * (t88 * t108 + t145 * t82) - t119 * t142) * MDP(21) + (-g(1) * (-t105 * t89 + t112 * t138) - g(2) * (-t105 * t90 - t138 * t91) - (t105 * t107 + t138 * t144) * t142) * MDP(27) + (-g(1) * (-t109 * t89 - t112 * t139) - g(2) * (-t109 * t90 + t139 * t91) - (-t101 * t132 + t107 * t109) * t142) * MDP(28) + t153 * (-g(1) * t89 - g(2) * t90 + g(3) * t137) + (MDP(16) * t110 - MDP(17) * t106 + MDP(18) * t101 - MDP(19) * t100 + MDP(9)) * t114; (-g(1) * (t108 * t116 - t122) + g(2) * t81 - t113) * MDP(16) + (g(1) * t120 - g(2) * (t106 * t133 + t90 * t110) - g(3) * (-t107 * t134 - t129)) * MDP(17) + t152 * MDP(19) + (-g(1) * (-pkin(3) * t122 + t86 * t108) - g(2) * (-t131 * t143 - t145 * t86) - pkin(3) * t113) * MDP(21) + (-MDP(27) * t109 + t158) * (g(2) * t117 + g(3) * (-t100 * t137 + t101 * t140) + g(1) * t79); -t114 * MDP(21); (t152 * t105 - t114 * t109) * MDP(27) + (g(1) * t151 - g(2) * t123 - g(3) * (t103 * t132 - t84 * t109)) * MDP(28);];
taug = t1;
