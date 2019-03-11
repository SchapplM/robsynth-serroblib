% Calculate Gravitation load on the joints for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:50
% EndTime: 2019-03-09 13:34:52
% DurationCPUTime: 0.54s
% Computational Cost: add. (456->99), mult. (942->173), div. (0->0), fcn. (1177->14), ass. (0->56)
t127 = sin(qJ(6));
t131 = cos(qJ(6));
t123 = qJ(4) + qJ(5);
t121 = sin(t123);
t122 = cos(t123);
t124 = sin(pkin(12));
t129 = sin(qJ(2));
t133 = cos(qJ(2));
t159 = cos(pkin(12));
t115 = -t133 * t124 - t129 * t159;
t126 = cos(pkin(6));
t108 = t115 * t126;
t130 = sin(qJ(1));
t134 = cos(qJ(1));
t138 = -t129 * t124 + t133 * t159;
t143 = -t134 * t108 + t130 * t138;
t125 = sin(pkin(6));
t155 = t125 * t134;
t91 = -t121 * t155 + t122 * t143;
t135 = t138 * t126;
t98 = t130 * t115 + t134 * t135;
t166 = t91 * t127 + t98 * t131;
t165 = -t98 * t127 + t91 * t131;
t162 = g(3) * t125;
t158 = t122 * t127;
t157 = t122 * t131;
t156 = t125 * t130;
t153 = t130 * t129;
t152 = t130 * t133;
t151 = t134 * t129;
t150 = t134 * t133;
t107 = t115 * t125;
t104 = -t107 * t122 + t126 * t121;
t141 = t121 * t143 + t122 * t155;
t142 = t130 * t108 + t134 * t138;
t93 = -t121 * t142 + t122 * t156;
t94 = t121 * t156 + t122 * t142;
t147 = (g(1) * t94 + g(2) * t91 + g(3) * t104) * MDP(26) + (-t131 * MDP(32) + t127 * MDP(33) - MDP(25)) * (g(1) * t93 - g(2) * t141 + g(3) * (t107 * t121 + t126 * t122));
t128 = sin(qJ(4));
t132 = cos(qJ(4));
t146 = -t128 * t155 + t132 * t143;
t144 = g(1) * t130 - g(2) * t134;
t140 = t128 * t143 + t132 * t155;
t139 = t126 * t150 - t153;
t112 = -t126 * t152 - t151;
t120 = t133 * pkin(2) + pkin(1);
t113 = -t126 * t153 + t150;
t111 = -t126 * t151 - t152;
t109 = t126 * t129 * pkin(2) + (-pkin(8) - qJ(3)) * t125;
t106 = t138 * t125;
t101 = t134 * t115 - t130 * t135;
t96 = t128 * t156 + t132 * t142;
t95 = -t128 * t142 + t132 * t156;
t89 = -t101 * t127 + t94 * t131;
t88 = -t101 * t131 - t94 * t127;
t1 = [t144 * MDP(2) + (-g(1) * t111 - g(2) * t113) * MDP(9) + (g(1) * t139 - g(2) * t112) * MDP(10) + (-g(1) * (-t134 * t109 - t130 * t120) - g(2) * (-t130 * t109 + t134 * t120)) * MDP(12) + (g(1) * t146 - g(2) * t96) * MDP(18) + (-g(1) * t140 - g(2) * t95) * MDP(19) + (g(1) * t91 - g(2) * t94) * MDP(25) + (-g(1) * t141 - g(2) * t93) * MDP(26) + (g(1) * t165 - g(2) * t89) * MDP(32) + (-g(1) * t166 - g(2) * t88) * MDP(33) + (-t125 * MDP(11) + MDP(3)) * (g(1) * t134 + g(2) * t130); (g(1) * t113 - g(2) * t111 + t129 * t162) * MDP(10) + (-g(1) * (t101 * t157 + t127 * t142) - g(2) * (t127 * t143 + t98 * t157) - g(3) * (t106 * t157 - t107 * t127)) * MDP(32) + (-g(1) * (-t101 * t158 + t131 * t142) - g(2) * (t131 * t143 - t98 * t158) - g(3) * (-t106 * t158 - t107 * t131)) * MDP(33) + (-MDP(18) * t132 + MDP(19) * t128 - t122 * MDP(25) + MDP(26) * t121) * (g(1) * t101 + g(2) * t98 + g(3) * t106) + (pkin(2) * MDP(12) + MDP(9)) * (-g(1) * t112 - g(2) * t139 - t133 * t162); (-g(3) * t126 - t125 * t144) * MDP(12); (-g(1) * t95 + g(2) * t140 - g(3) * (t107 * t128 + t126 * t132)) * MDP(18) + (g(1) * t96 + g(2) * t146 - g(3) * (t107 * t132 - t126 * t128)) * MDP(19) + t147; t147; (-g(1) * t88 + g(2) * t166 - g(3) * (-t104 * t127 - t106 * t131)) * MDP(32) + (g(1) * t89 + g(2) * t165 - g(3) * (-t104 * t131 + t106 * t127)) * MDP(33);];
taug  = t1;
