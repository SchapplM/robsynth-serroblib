% Calculate Gravitation load on the joints for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:34:10
% EndTime: 2019-03-09 22:34:12
% DurationCPUTime: 0.68s
% Computational Cost: add. (445->116), mult. (686->197), div. (0->0), fcn. (788->14), ass. (0->56)
t157 = MDP(10) - MDP(25);
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t114 = qJ(3) + qJ(4);
t109 = pkin(12) + t114;
t106 = sin(t109);
t107 = cos(t109);
t115 = sin(pkin(6));
t123 = cos(qJ(1));
t139 = t115 * t123;
t118 = sin(qJ(2));
t119 = sin(qJ(1));
t122 = cos(qJ(2));
t145 = cos(pkin(6));
t135 = t123 * t145;
t96 = t118 * t135 + t119 * t122;
t86 = t106 * t139 - t96 * t107;
t95 = t118 * t119 - t122 * t135;
t156 = t116 * t86 + t120 * t95;
t155 = -t116 * t95 + t120 * t86;
t154 = g(1) * t123 + g(2) * t119;
t110 = sin(t114);
t111 = cos(t114);
t130 = t96 * t110 + t111 * t139;
t142 = t115 * t118;
t141 = t115 * t119;
t136 = t119 * t145;
t98 = -t118 * t136 + t122 * t123;
t89 = -t110 * t98 + t111 * t141;
t153 = -g(3) * (-t110 * t142 + t111 * t145) + g(2) * t130 - g(1) * t89;
t149 = g(3) * t115;
t144 = t107 * t116;
t143 = t107 * t120;
t121 = cos(qJ(3));
t140 = t115 * t121;
t138 = t116 * t122;
t137 = t120 * t122;
t101 = t121 * pkin(3) + pkin(4) * t111;
t133 = -t110 * t139 + t111 * t96;
t90 = t110 * t141 + t111 * t98;
t134 = t153 * MDP(23) + (g(1) * t90 + g(2) * t133 - g(3) * (-t110 * t145 - t111 * t142)) * MDP(24) + (-t120 * MDP(32) + t116 * MDP(33)) * (g(1) * (-t106 * t98 + t107 * t141) + g(2) * (-t96 * t106 - t107 * t139) + g(3) * (-t106 * t142 + t107 * t145));
t117 = sin(qJ(3));
t132 = -t117 * t139 + t121 * t96;
t129 = t96 * t117 + t121 * t139;
t97 = t123 * t118 + t122 * t136;
t125 = -g(1) * t97 - g(2) * t95 + t122 * t149;
t113 = -qJ(5) - pkin(10) - pkin(9);
t100 = pkin(3) * t117 + pkin(4) * t110;
t99 = pkin(2) + t101;
t94 = t106 * t145 + t107 * t142;
t92 = t117 * t141 + t121 * t98;
t91 = -t117 * t98 + t119 * t140;
t88 = t106 * t141 + t107 * t98;
t83 = t116 * t97 + t120 * t88;
t82 = -t116 * t88 + t120 * t97;
t1 = [(g(1) * t119 - g(2) * t123) * MDP(2) + t154 * MDP(3) + (g(1) * t96 - g(2) * t98) * MDP(9) + (g(1) * t132 - g(2) * t92) * MDP(16) + (-g(1) * t129 - g(2) * t91) * MDP(17) + (g(1) * t133 - g(2) * t90) * MDP(23) + (-g(1) * t130 - g(2) * t89) * MDP(24) + (-g(1) * (-t119 * pkin(1) + t95 * t113 - t96 * t99) - g(2) * (pkin(1) * t123 - t97 * t113 + t98 * t99) - t154 * t115 * (pkin(8) + t100)) * MDP(26) + (-g(1) * t155 - g(2) * t83) * MDP(32) + (g(1) * t156 - g(2) * t82) * MDP(33) - t157 * (g(1) * t95 - g(2) * t97); (-g(1) * (-t113 * t98 - t97 * t99) - g(2) * (-t113 * t96 - t95 * t99) - (-t113 * t118 + t122 * t99) * t149) * MDP(26) + (-g(1) * (t116 * t98 - t143 * t97) - g(2) * (t116 * t96 - t143 * t95) - (t107 * t137 + t116 * t118) * t149) * MDP(32) + (-g(1) * (t120 * t98 + t144 * t97) - g(2) * (t120 * t96 + t144 * t95) - (-t107 * t138 + t118 * t120) * t149) * MDP(33) + t157 * (g(1) * t98 + g(2) * t96 + g(3) * t142) + (-t121 * MDP(16) + t117 * MDP(17) - MDP(23) * t111 + MDP(24) * t110 - MDP(9)) * t125; (-g(1) * t91 + g(2) * t129 - g(3) * (-t117 * t142 + t121 * t145)) * MDP(16) + (g(1) * t92 + g(2) * t132 - g(3) * (-t117 * t145 - t118 * t140)) * MDP(17) + (-g(1) * (-t100 * t98 + t101 * t141) - g(2) * (-t96 * t100 - t101 * t139) - g(3) * (-t100 * t142 + t101 * t145)) * MDP(26) + t134; t153 * MDP(26) * pkin(4) + t134; t125 * MDP(26); (-g(1) * t82 - g(2) * t156 - g(3) * (-t115 * t137 - t116 * t94)) * MDP(32) + (g(1) * t83 - g(2) * t155 - g(3) * (t115 * t138 - t120 * t94)) * MDP(33);];
taug  = t1;
