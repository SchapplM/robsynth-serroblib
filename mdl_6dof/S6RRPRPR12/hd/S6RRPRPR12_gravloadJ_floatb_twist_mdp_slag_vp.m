% Calculate Gravitation load on the joints for
% S6RRPRPR12
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
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:48
% EndTime: 2019-03-09 11:21:50
% DurationCPUTime: 0.77s
% Computational Cost: add. (297->111), mult. (622->178), div. (0->0), fcn. (708->12), ass. (0->53)
t156 = -MDP(10) + MDP(13);
t155 = MDP(12) - MDP(9) - MDP(22);
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t111 = qJ(4) + pkin(11);
t108 = sin(t111);
t109 = cos(t111);
t112 = sin(pkin(6));
t121 = cos(qJ(1));
t139 = t112 * t121;
t116 = sin(qJ(2));
t117 = sin(qJ(1));
t120 = cos(qJ(2));
t145 = cos(pkin(6));
t132 = t121 * t145;
t96 = t117 * t116 - t120 * t132;
t86 = -t96 * t108 + t109 * t139;
t97 = t116 * t132 + t117 * t120;
t154 = t86 * t114 + t97 * t118;
t153 = -t97 * t114 + t86 * t118;
t115 = sin(qJ(4));
t149 = pkin(4) * t115;
t148 = g(3) * t112;
t144 = t108 * t114;
t143 = t108 * t118;
t142 = t112 * t116;
t141 = t112 * t117;
t140 = t112 * t120;
t138 = t114 * t116;
t137 = t116 * t118;
t136 = pkin(2) * t140 + qJ(3) * t142;
t133 = t117 * t145;
t99 = -t116 * t133 + t121 * t120;
t135 = t121 * pkin(1) + t99 * pkin(2) + pkin(8) * t141;
t134 = qJ(3) + t149;
t119 = cos(qJ(4));
t131 = t96 * t115 - t119 * t139;
t130 = -t117 * pkin(1) - t97 * pkin(2) + pkin(8) * t139;
t98 = t121 * t116 + t120 * t133;
t88 = -t115 * t141 + t98 * t119;
t126 = t115 * t139 + t96 * t119;
t123 = g(1) * t99 + g(2) * t97 + g(3) * t142;
t79 = -g(1) * t98 - g(2) * t96 + g(3) * t140;
t113 = -qJ(5) - pkin(9);
t107 = t119 * pkin(4) + pkin(3);
t94 = t98 * pkin(2);
t92 = t96 * pkin(2);
t91 = -t108 * t140 + t145 * t109;
t89 = t98 * t115 + t119 * t141;
t84 = t98 * t108 + t109 * t141;
t78 = t99 * t114 + t84 * t118;
t77 = -t84 * t114 + t99 * t118;
t1 = [(g(1) * t117 - g(2) * t121) * MDP(2) + (-g(1) * (-t96 * qJ(3) + t130) - g(2) * (t98 * qJ(3) + t135)) * MDP(14) + (g(1) * t131 - g(2) * t89) * MDP(20) + (g(1) * t126 - g(2) * t88) * MDP(21) + (-g(1) * (t107 * t139 + t97 * t113 - t134 * t96 + t130) - g(2) * (t107 * t141 - t99 * t113 + t134 * t98 + t135)) * MDP(23) + (-g(1) * t153 - g(2) * t78) * MDP(29) + (g(1) * t154 - g(2) * t77) * MDP(30) + t156 * (g(1) * t96 - g(2) * t98) - t155 * (g(1) * t97 - g(2) * t99) + (-t112 * MDP(11) + MDP(3)) * (g(1) * t121 + g(2) * t117); (-g(1) * (t99 * qJ(3) - t94) - g(2) * (t97 * qJ(3) - t92) - g(3) * t136) * MDP(14) + (-g(1) * (t98 * t113 + t134 * t99 - t94) - g(2) * (t96 * t113 + t134 * t97 - t92) - g(3) * ((-t113 * t120 + t116 * t149) * t112 + t136)) * MDP(23) + (-g(1) * (-t98 * t114 + t99 * t143) - g(2) * (-t96 * t114 + t97 * t143) - (t108 * t137 + t114 * t120) * t148) * MDP(29) + (-g(1) * (-t98 * t118 - t99 * t144) - g(2) * (-t96 * t118 - t97 * t144) - (-t108 * t138 + t118 * t120) * t148) * MDP(30) + t155 * t79 + (-t115 * MDP(20) - MDP(21) * t119 - t156) * t123; (MDP(14) + MDP(23)) * t79; (g(1) * t89 + g(2) * t131 - g(3) * (t115 * t140 - t145 * t119)) * MDP(21) + (-MDP(29) * t118 + MDP(30) * t114) * (g(1) * (-t108 * t141 + t98 * t109) + g(2) * (t108 * t139 + t96 * t109) + g(3) * (-t145 * t108 - t109 * t140)) + (MDP(23) * pkin(4) + MDP(20)) * (-g(3) * (-t145 * t115 - t119 * t140) - g(2) * t126 - g(1) * t88); -t123 * MDP(23); (-g(1) * t77 - g(2) * t154 - g(3) * (t112 * t137 - t91 * t114)) * MDP(29) + (g(1) * t78 - g(2) * t153 - g(3) * (-t112 * t138 - t91 * t118)) * MDP(30);];
taug  = t1;
