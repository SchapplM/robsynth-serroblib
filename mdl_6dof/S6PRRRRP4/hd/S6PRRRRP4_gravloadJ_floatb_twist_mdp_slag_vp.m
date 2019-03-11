% Calculate Gravitation load on the joints for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:48
% EndTime: 2019-03-09 00:17:50
% DurationCPUTime: 0.77s
% Computational Cost: add. (582->109), mult. (1138->176), div. (0->0), fcn. (1387->12), ass. (0->58)
t177 = MDP(11) - MDP(27);
t176 = MDP(24) + MDP(26);
t175 = MDP(25) - MDP(28);
t133 = sin(pkin(11));
t137 = sin(qJ(2));
t140 = cos(qJ(2));
t166 = cos(pkin(11));
t167 = cos(pkin(6));
t151 = t167 * t166;
t118 = t133 * t137 - t140 * t151;
t154 = t133 * t167;
t120 = t166 * t137 + t140 * t154;
t174 = -g(1) * t120 - g(2) * t118;
t134 = sin(pkin(6));
t168 = g(3) * t134;
t119 = t133 * t140 + t137 * t151;
t135 = sin(qJ(4));
t165 = t119 * t135;
t121 = -t137 * t154 + t166 * t140;
t164 = t121 * t135;
t132 = qJ(4) + qJ(5);
t130 = sin(t132);
t139 = cos(qJ(3));
t163 = t130 * t139;
t131 = cos(t132);
t162 = t131 * t139;
t161 = t134 * t137;
t160 = t134 * t139;
t159 = t134 * t140;
t158 = t135 * t139;
t138 = cos(qJ(4));
t157 = t138 * t139;
t156 = t139 * t140;
t155 = t130 * t159;
t153 = t134 * t166;
t136 = sin(qJ(3));
t123 = t167 * t136 + t137 * t160;
t107 = t123 * t131 - t155;
t106 = t123 * t130 + t131 * t159;
t109 = t119 * t139 - t136 * t153;
t96 = t109 * t130 - t118 * t131;
t111 = t133 * t134 * t136 + t121 * t139;
t98 = t111 * t130 - t120 * t131;
t90 = g(1) * t98 + g(2) * t96 + g(3) * t106;
t97 = t109 * t131 + t118 * t130;
t99 = t111 * t131 + t120 * t130;
t152 = t176 * t90 + t175 * (g(1) * t99 + g(2) * t97 + g(3) * t107);
t143 = -g(1) * (-t98 * pkin(5) + t99 * qJ(6)) - g(2) * (-t96 * pkin(5) + t97 * qJ(6)) - g(3) * (-t106 * pkin(5) + t107 * qJ(6));
t142 = -g(1) * (-t111 * t135 + t120 * t138) - g(2) * (-t109 * t135 + t118 * t138) - g(3) * (-t123 * t135 - t138 * t159);
t141 = -pkin(10) - pkin(9);
t129 = pkin(4) * t138 + pkin(3);
t113 = (t130 * t137 + t131 * t156) * t134;
t112 = -t131 * t161 + t139 * t155;
t104 = -t120 * t162 + t121 * t130;
t103 = -t120 * t163 - t121 * t131;
t102 = -t118 * t162 + t119 * t130;
t101 = -t118 * t163 - t119 * t131;
t1 = [(-MDP(1) - MDP(29)) * g(3); (g(1) * t121 + g(2) * t119 + g(3) * t161) * MDP(4) + (-g(1) * (-t120 * t157 + t164) - g(2) * (-t118 * t157 + t165) - (t135 * t137 + t138 * t156) * t168) * MDP(17) + (-g(1) * (t120 * t158 + t121 * t138) - g(2) * (t118 * t158 + t119 * t138) - (-t135 * t156 + t137 * t138) * t168) * MDP(18) + (-g(1) * (pkin(4) * t164 + t104 * pkin(5) + t121 * pkin(8) + t103 * qJ(6)) - g(2) * (pkin(4) * t165 + t102 * pkin(5) + t119 * pkin(8) + t101 * qJ(6)) - g(3) * (t113 * pkin(5) + t112 * qJ(6)) - (pkin(4) * t135 + pkin(8)) * t137 * t168 + (-t140 * t168 - t174) * (t129 * t139 - t136 * t141 + pkin(2))) * MDP(29) + t176 * (-g(1) * t104 - g(2) * t102 - g(3) * t113) + t175 * (g(1) * t103 + g(2) * t101 + g(3) * t112) + (-t139 * MDP(10) + t136 * t177 - MDP(3)) * (g(3) * t159 + t174); (t141 * MDP(29) + t177) * (g(1) * t111 + g(2) * t109 + g(3) * t123) + (-t176 * t131 + t175 * t130 - MDP(29) * (pkin(5) * t131 + qJ(6) * t130 + t129) - MDP(17) * t138 + MDP(18) * t135 - MDP(10)) * (g(3) * (-t136 * t161 + t167 * t139) + g(2) * (-t119 * t136 - t139 * t153) + g(1) * (-t121 * t136 + t133 * t160)); t142 * MDP(17) + (-g(1) * (-t111 * t138 - t120 * t135) - g(2) * (-t109 * t138 - t118 * t135) - g(3) * (-t123 * t138 + t135 * t159)) * MDP(18) + (t142 * pkin(4) + t143) * MDP(29) + t152; t143 * MDP(29) + t152; -t90 * MDP(29);];
taug  = t1;
