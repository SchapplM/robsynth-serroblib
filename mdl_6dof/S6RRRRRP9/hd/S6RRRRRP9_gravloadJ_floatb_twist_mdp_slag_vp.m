% Calculate Gravitation load on the joints for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:40
% EndTime: 2019-03-10 02:11:43
% DurationCPUTime: 0.86s
% Computational Cost: add. (484->133), mult. (966->218), div. (0->0), fcn. (1147->12), ass. (0->57)
t107 = qJ(4) + qJ(5);
t103 = sin(t107);
t104 = cos(t107);
t110 = sin(qJ(3));
t114 = cos(qJ(3));
t108 = sin(pkin(6));
t141 = cos(qJ(1));
t124 = t108 * t141;
t111 = sin(qJ(2));
t112 = sin(qJ(1));
t115 = cos(qJ(2));
t134 = cos(pkin(6));
t121 = t134 * t141;
t91 = t111 * t121 + t112 * t115;
t83 = -t110 * t124 + t114 * t91;
t90 = t111 * t112 - t115 * t121;
t151 = t103 * t83 - t104 * t90;
t150 = t103 * t90 + t104 * t83;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t149 = t109 * t83 - t113 * t90;
t148 = t109 * t90 + t113 * t83;
t128 = t108 * t115;
t130 = t108 * t112;
t123 = t112 * t134;
t93 = -t111 * t123 + t141 * t115;
t87 = t110 * t130 + t114 * t93;
t92 = t141 * t111 + t115 * t123;
t77 = -t103 * t87 + t104 * t92;
t129 = t108 * t114;
t89 = t134 * t110 + t111 * t129;
t147 = g(2) * t151 - g(3) * (-t89 * t103 - t104 * t128) - g(1) * t77;
t146 = MDP(17) - MDP(32);
t144 = g(2) * t90;
t143 = g(2) * t91;
t95 = pkin(4) * t109 + pkin(5) * t103;
t142 = pkin(9) + t95;
t140 = g(3) * t108;
t78 = t103 * t92 + t104 * t87;
t139 = t147 * MDP(30) + (g(1) * t78 + g(2) * t150 - g(3) * (t103 * t128 - t89 * t104)) * MDP(31);
t133 = t103 * t114;
t132 = t104 * t114;
t131 = t108 * t111;
t127 = t109 * t114;
t126 = t113 * t114;
t125 = t114 * t115;
t96 = t113 * pkin(4) + pkin(5) * t104;
t106 = -qJ(6) - pkin(11) - pkin(10);
t94 = pkin(3) + t96;
t120 = t106 * t110 - t114 * t94 - pkin(2);
t82 = t91 * t110 + t114 * t124;
t86 = t110 * t93 - t112 * t129;
t88 = t110 * t131 - t134 * t114;
t119 = g(1) * t86 + g(2) * t82 + g(3) * t88;
t80 = t109 * t92 + t113 * t87;
t79 = -t109 * t87 + t113 * t92;
t1 = [(g(1) * t112 - g(2) * t141) * MDP(2) + (g(1) * t141 + g(2) * t112) * MDP(3) + (g(1) * t91 - g(2) * t93) * MDP(9) + (-g(1) * t90 + g(2) * t92) * MDP(10) + (g(1) * t83 - g(2) * t87) * MDP(16) + (g(1) * t148 - g(2) * t80) * MDP(23) + (-g(1) * t149 - g(2) * t79) * MDP(24) + (g(1) * t150 - g(2) * t78) * MDP(30) + (-g(1) * t151 - g(2) * t77) * MDP(31) + (-g(1) * (-t112 * pkin(1) - t91 * pkin(2) + pkin(8) * t124 + t106 * t82 - t142 * t90 - t83 * t94) - g(2) * (t141 * pkin(1) + t93 * pkin(2) + pkin(8) * t130 - t86 * t106 + t142 * t92 + t87 * t94)) * MDP(33) + t146 * (-g(1) * t82 + g(2) * t86); (g(1) * t93 + g(3) * t131 + t143) * MDP(10) + (-g(1) * (t109 * t93 - t92 * t126) - g(2) * (t109 * t91 - t90 * t126) - (t109 * t111 + t113 * t125) * t140) * MDP(23) + (-g(1) * (t113 * t93 + t92 * t127) - g(2) * (t113 * t91 + t90 * t127) - (-t109 * t125 + t111 * t113) * t140) * MDP(24) + (-g(1) * (t103 * t93 - t92 * t132) - g(2) * (t103 * t91 - t90 * t132) - (t103 * t111 + t104 * t125) * t140) * MDP(30) + (-g(1) * (t104 * t93 + t92 * t133) - g(2) * (t104 * t91 + t90 * t133) - (-t103 * t125 + t104 * t111) * t140) * MDP(31) + (-g(1) * (t120 * t92 + t142 * t93) - t142 * t143 - t120 * t144 - (t142 * t111 - t120 * t115) * t140) * MDP(33) + (-t114 * MDP(16) + t146 * t110 - MDP(9)) * (-g(1) * t92 + g(3) * t128 - t144); (-g(1) * (-t106 * t87 - t86 * t94) - g(2) * (-t106 * t83 - t82 * t94) - g(3) * (-t106 * t89 - t88 * t94)) * MDP(33) + t146 * (g(1) * t87 + g(2) * t83 + g(3) * t89) + (MDP(23) * t113 - MDP(24) * t109 + MDP(30) * t104 - MDP(31) * t103 + MDP(16)) * t119; (-g(1) * t79 + g(2) * t149 - g(3) * (-t89 * t109 - t113 * t128)) * MDP(23) + (g(1) * t80 + g(2) * t148 - g(3) * (t109 * t128 - t89 * t113)) * MDP(24) + (-g(1) * (-t87 * t95 + t92 * t96) - g(2) * (-t83 * t95 + t90 * t96) - g(3) * (-t96 * t128 - t89 * t95)) * MDP(33) + t139; t147 * MDP(33) * pkin(5) + t139; -t119 * MDP(33);];
taug  = t1;
