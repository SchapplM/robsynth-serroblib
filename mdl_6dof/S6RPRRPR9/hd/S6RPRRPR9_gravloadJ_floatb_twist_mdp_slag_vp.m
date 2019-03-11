% Calculate Gravitation load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:36
% EndTime: 2019-03-09 05:31:38
% DurationCPUTime: 0.90s
% Computational Cost: add. (549->120), mult. (1381->202), div. (0->0), fcn. (1742->16), ass. (0->63)
t138 = sin(qJ(1));
t141 = cos(qJ(1));
t166 = sin(pkin(12));
t170 = cos(pkin(6));
t150 = t170 * t166;
t168 = cos(pkin(12));
t119 = t138 * t168 + t141 * t150;
t137 = sin(qJ(3));
t133 = sin(pkin(6));
t167 = sin(pkin(7));
t156 = t133 * t167;
t171 = cos(qJ(3));
t148 = t171 * t156;
t152 = t170 * t168;
t118 = t138 * t166 - t141 * t152;
t169 = cos(pkin(7));
t158 = t118 * t169;
t102 = t119 * t137 + t141 * t148 + t171 * t158;
t135 = sin(qJ(6));
t139 = cos(qJ(6));
t103 = t119 * t171 + (-t141 * t156 - t158) * t137;
t157 = t133 * t169;
t111 = t118 * t167 - t141 * t157;
t132 = qJ(4) + pkin(13);
t129 = sin(t132);
t130 = cos(t132);
t94 = t103 * t130 + t111 * t129;
t181 = -t102 * t139 + t94 * t135;
t180 = t102 * t135 + t94 * t139;
t140 = cos(qJ(4));
t136 = sin(qJ(4));
t164 = t111 * t136;
t179 = t103 * t140 + t164;
t175 = t103 * t136 - t111 * t140;
t173 = MDP(14) - MDP(22);
t165 = qJ(2) * t133;
t144 = t138 * t152 + t141 * t166;
t112 = t138 * t157 + t144 * t167;
t163 = t112 * t136;
t162 = t130 * t135;
t161 = t130 * t139;
t160 = t141 * pkin(1) + t138 * t165;
t159 = -t138 * pkin(1) + t141 * t165;
t153 = g(1) * t138 - g(2) * t141;
t151 = t170 * t167;
t149 = t169 * t168;
t120 = -t138 * t150 + t141 * t168;
t143 = t144 * t169;
t107 = t120 * t171 + (t138 * t156 - t143) * t137;
t98 = -t107 * t136 + t112 * t140;
t106 = t120 * t137 - t138 * t148 + t171 * t143;
t108 = -t171 * t151 + (t137 * t166 - t149 * t171) * t133;
t146 = g(1) * t106 + g(2) * t102 + g(3) * t108;
t134 = -qJ(5) - pkin(10);
t128 = t140 * pkin(4) + pkin(3);
t117 = -t168 * t156 + t170 * t169;
t109 = t137 * t151 + (t137 * t149 + t171 * t166) * t133;
t101 = t109 * t130 + t117 * t129;
t99 = t107 * t140 + t163;
t97 = t107 * t130 + t112 * t129;
t92 = t106 * t135 + t97 * t139;
t91 = t106 * t139 - t97 * t135;
t1 = [t153 * MDP(2) + (g(1) * t119 - g(2) * t120) * MDP(4) + (-g(1) * t118 + g(2) * t144) * MDP(5) + (-g(1) * t159 - g(2) * t160) * MDP(7) + (g(1) * t103 - g(2) * t107) * MDP(13) + (g(1) * t179 - g(2) * t99) * MDP(20) + (-g(1) * t175 - g(2) * t98) * MDP(21) + (-g(1) * (-t119 * pkin(2) - pkin(4) * t164 + t102 * t134 - t103 * t128 + t159) - g(2) * (t120 * pkin(2) + pkin(4) * t163 - t106 * t134 + t107 * t128 + t160) + (g(1) * t111 - g(2) * t112) * pkin(9)) * MDP(23) + (g(1) * t180 - g(2) * t92) * MDP(29) + (-g(1) * t181 - g(2) * t91) * MDP(30) + t173 * (-g(1) * t102 + g(2) * t106) + (-t133 * MDP(6) + MDP(3)) * (g(1) * t141 + g(2) * t138); (MDP(23) + MDP(7)) * (-g(3) * t170 - t153 * t133); (-g(1) * (-t106 * t128 - t107 * t134) - g(2) * (-t102 * t128 - t103 * t134) - g(3) * (-t108 * t128 - t109 * t134)) * MDP(23) + (-g(1) * (-t106 * t161 + t107 * t135) - g(2) * (-t102 * t161 + t103 * t135) - g(3) * (-t108 * t161 + t109 * t135)) * MDP(29) + (-g(1) * (t106 * t162 + t107 * t139) - g(2) * (t102 * t162 + t103 * t139) - g(3) * (t108 * t162 + t109 * t139)) * MDP(30) + t173 * (g(1) * t107 + g(2) * t103 + g(3) * t109) + (MDP(20) * t140 - MDP(21) * t136 + MDP(13)) * t146; (g(1) * t99 + g(2) * t179 - g(3) * (-t109 * t140 - t117 * t136)) * MDP(21) + (-MDP(29) * t139 + MDP(30) * t135) * (g(1) * (-t107 * t129 + t112 * t130) + g(2) * (-t103 * t129 + t111 * t130) + g(3) * (-t109 * t129 + t117 * t130)) + (pkin(4) * MDP(23) + MDP(20)) * (-g(1) * t98 + g(2) * t175 - g(3) * (-t109 * t136 + t117 * t140)); -t146 * MDP(23); (-g(1) * t91 + g(2) * t181 - g(3) * (-t101 * t135 + t108 * t139)) * MDP(29) + (g(1) * t92 + g(2) * t180 - g(3) * (-t101 * t139 - t108 * t135)) * MDP(30);];
taug  = t1;
