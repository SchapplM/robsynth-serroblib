% Calculate Gravitation load on the joints for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:35
% EndTime: 2019-03-09 11:29:38
% DurationCPUTime: 0.99s
% Computational Cost: add. (380->131), mult. (850->208), div. (0->0), fcn. (995->12), ass. (0->52)
t157 = MDP(9) - MDP(12);
t156 = -MDP(10) + MDP(13);
t154 = MDP(21) - MDP(24);
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t120 = cos(qJ(1));
t117 = sin(qJ(1));
t146 = cos(pkin(6));
t135 = t117 * t146;
t100 = -t116 * t135 + t119 * t120;
t134 = t120 * t146;
t98 = t116 * t134 + t117 * t119;
t155 = -g(1) * t100 - g(2) * t98;
t111 = pkin(11) + qJ(6);
t108 = sin(t111);
t109 = cos(t111);
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t113 = sin(pkin(6));
t139 = t113 * t120;
t97 = t116 * t117 - t119 * t134;
t126 = -t115 * t97 + t118 * t139;
t153 = t108 * t126 + t109 * t98;
t152 = -t108 * t98 + t109 * t126;
t149 = g(3) * t113;
t145 = t108 * t115;
t144 = t109 * t115;
t112 = sin(pkin(11));
t143 = t112 * t115;
t142 = t113 * t116;
t141 = t113 * t117;
t140 = t113 * t119;
t114 = cos(pkin(11));
t138 = t114 * t115;
t137 = t115 * t116;
t136 = g(3) * (pkin(2) * t140 + qJ(3) * t142);
t129 = pkin(4) * t115 - qJ(5) * t118;
t99 = t116 * t120 + t119 * t135;
t128 = pkin(1) * t120 + pkin(2) * t100 + pkin(8) * t141 + qJ(3) * t99;
t125 = t115 * t139 + t118 * t97;
t82 = t115 * t141 - t118 * t99;
t95 = t115 * t146 + t118 * t140;
t124 = g(1) * t82 - g(2) * t125 + g(3) * t95;
t122 = -pkin(1) * t117 - pkin(2) * t98 + pkin(8) * t139 - qJ(3) * t97;
t81 = -g(1) * t99 - g(2) * t97 + g(3) * t140;
t96 = -t115 * t140 + t118 * t146;
t93 = t99 * pkin(2);
t91 = t97 * pkin(2);
t83 = t115 * t99 + t118 * t141;
t79 = t100 * t108 + t109 * t83;
t78 = t100 * t109 - t108 * t83;
t1 = [(g(1) * t117 - g(2) * t120) * MDP(2) + (-g(1) * t122 - g(2) * t128) * MDP(14) + (-g(1) * t126 - g(2) * t83) * MDP(20) + (-g(1) * (-t112 * t98 + t114 * t126) - g(2) * (t100 * t112 + t114 * t83)) * MDP(22) + (-g(1) * (-t112 * t126 - t114 * t98) - g(2) * (t100 * t114 - t112 * t83)) * MDP(23) + (-g(1) * (pkin(3) * t139 + pkin(4) * t126 - pkin(9) * t98 + qJ(5) * t125 + t122) - g(2) * (pkin(3) * t141 + pkin(4) * t83 + pkin(9) * t100 + qJ(5) * t82 + t128)) * MDP(25) + (-g(1) * t152 - g(2) * t79) * MDP(31) + (g(1) * t153 - g(2) * t78) * MDP(32) + t154 * (g(1) * t125 + g(2) * t82) + t156 * (g(1) * t97 - g(2) * t99) + t157 * (g(1) * t98 - g(2) * t100) + (-MDP(11) * t113 + MDP(3)) * (g(1) * t120 + g(2) * t117); (-g(1) * (qJ(3) * t100 - t93) - g(2) * (qJ(3) * t98 - t91) - t136) * MDP(14) + (-g(1) * (t100 * t138 - t112 * t99) - g(2) * (-t112 * t97 + t138 * t98) - (t112 * t119 + t114 * t137) * t149) * MDP(22) + (-g(1) * (-t100 * t143 - t114 * t99) - g(2) * (-t114 * t97 - t143 * t98) - (-t112 * t137 + t114 * t119) * t149) * MDP(23) + (-g(1) * (-pkin(9) * t99 - t93) - g(2) * (-pkin(9) * t97 - t91) - t136 - (pkin(9) * t119 + t116 * t129) * t149 + t155 * (qJ(3) + t129)) * MDP(25) + (-g(1) * (t100 * t144 - t108 * t99) - g(2) * (-t108 * t97 + t144 * t98) - (t108 * t119 + t109 * t137) * t149) * MDP(31) + (-g(1) * (-t100 * t145 - t109 * t99) - g(2) * (-t109 * t97 - t145 * t98) - (-t108 * t137 + t109 * t119) * t149) * MDP(32) - t157 * t81 + (-t115 * MDP(20) - t154 * t118 - t156) * (g(3) * t142 - t155); (MDP(14) + MDP(25)) * t81; (-g(1) * (-pkin(4) * t82 + qJ(5) * t83) - g(2) * (pkin(4) * t125 - qJ(5) * t126) - g(3) * (-pkin(4) * t95 + qJ(5) * t96)) * MDP(25) + t154 * (g(1) * t83 - g(2) * t126 + g(3) * t96) + (MDP(22) * t114 - MDP(23) * t112 + MDP(31) * t109 - MDP(32) * t108 + MDP(20)) * t124; -t124 * MDP(25); (-g(1) * t78 - g(2) * t153 - g(3) * (-t108 * t96 + t109 * t142)) * MDP(31) + (g(1) * t79 - g(2) * t152 - g(3) * (-t108 * t142 - t109 * t96)) * MDP(32);];
taug  = t1;
