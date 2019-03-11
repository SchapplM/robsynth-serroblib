% Calculate Gravitation load on the joints for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:11
% EndTime: 2019-03-09 07:34:14
% DurationCPUTime: 0.87s
% Computational Cost: add. (647->104), mult. (1542->184), div. (0->0), fcn. (1977->16), ass. (0->57)
t136 = cos(qJ(1));
t161 = sin(pkin(13));
t165 = cos(pkin(6));
t150 = t165 * t161;
t163 = cos(pkin(13));
t166 = sin(qJ(1));
t120 = t136 * t150 + t166 * t163;
t133 = sin(qJ(3));
t167 = cos(qJ(3));
t151 = t165 * t163;
t142 = -t136 * t151 + t166 * t161;
t130 = sin(pkin(6));
t162 = sin(pkin(7));
t156 = t130 * t162;
t164 = cos(pkin(7));
t177 = t136 * t156 + t142 * t164;
t106 = t120 * t133 + t167 * t177;
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t109 = -t120 * t167 + t133 * t177;
t157 = t130 * t164;
t115 = -t136 * t157 + t142 * t162;
t129 = qJ(4) + qJ(5);
t127 = sin(t129);
t128 = cos(t129);
t99 = t109 * t128 - t115 * t127;
t181 = t106 * t134 + t99 * t131;
t180 = -t106 * t131 + t99 * t134;
t176 = t109 * t127 + t115 * t128;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t175 = t109 * t132 + t115 * t135;
t174 = t109 * t135 - t115 * t132;
t138 = t136 * t161 + t166 * t151;
t168 = t138 * t164 - t166 * t156;
t160 = qJ(2) * t130;
t159 = t128 * t131;
t158 = t128 * t134;
t121 = t136 * t163 - t166 * t150;
t111 = t121 * t167 - t168 * t133;
t116 = t138 * t162 + t166 * t157;
t100 = -t111 * t127 + t116 * t128;
t101 = t111 * t128 + t116 * t127;
t148 = t162 * t165;
t149 = t164 * t163;
t113 = t133 * t148 + (t133 * t149 + t167 * t161) * t130;
t119 = -t163 * t156 + t165 * t164;
t105 = t113 * t128 + t119 * t127;
t155 = (g(1) * t101 - g(2) * t99 + g(3) * t105) * MDP(28) + (-t134 * MDP(34) + t131 * MDP(35) - MDP(27)) * (g(1) * t100 + g(2) * t176 + g(3) * (-t113 * t127 + t119 * t128));
t146 = g(1) * t166 - g(2) * t136;
t112 = -t167 * t148 + (t133 * t161 - t149 * t167) * t130;
t110 = t121 * t133 + t168 * t167;
t103 = t111 * t135 + t116 * t132;
t102 = -t111 * t132 + t116 * t135;
t96 = t101 * t134 + t110 * t131;
t95 = -t101 * t131 + t110 * t134;
t1 = [t146 * MDP(2) + (g(1) * t120 - g(2) * t121) * MDP(4) + (-g(1) * t142 + g(2) * t138) * MDP(5) + (-g(1) * (-t166 * pkin(1) + t136 * t160) - g(2) * (t136 * pkin(1) + t166 * t160)) * MDP(7) + (-g(1) * t109 - g(2) * t111) * MDP(13) + (-g(1) * t106 + g(2) * t110) * MDP(14) + (-g(1) * t174 - g(2) * t103) * MDP(20) + (g(1) * t175 - g(2) * t102) * MDP(21) + (-g(1) * t99 - g(2) * t101) * MDP(27) + (g(1) * t176 - g(2) * t100) * MDP(28) + (-g(1) * t180 - g(2) * t96) * MDP(34) + (g(1) * t181 - g(2) * t95) * MDP(35) + (t130 * MDP(6) - MDP(3)) * (-g(1) * t136 - g(2) * t166); (-g(3) * t165 - t146 * t130) * MDP(7); (g(1) * t111 - g(2) * t109 + g(3) * t113) * MDP(14) + (-g(1) * (-t110 * t158 + t111 * t131) - g(2) * (-t106 * t158 - t109 * t131) - g(3) * (-t112 * t158 + t113 * t131)) * MDP(34) + (-g(1) * (t110 * t159 + t111 * t134) - g(2) * (t106 * t159 - t109 * t134) - g(3) * (t112 * t159 + t113 * t134)) * MDP(35) + (MDP(20) * t135 - MDP(21) * t132 + t128 * MDP(27) - MDP(28) * t127 + MDP(13)) * (g(1) * t110 + g(2) * t106 + g(3) * t112); (-g(1) * t102 - g(2) * t175 - g(3) * (-t113 * t132 + t119 * t135)) * MDP(20) + (g(1) * t103 - g(2) * t174 - g(3) * (-t113 * t135 - t119 * t132)) * MDP(21) + t155; t155; (-g(1) * t95 - g(2) * t181 - g(3) * (-t105 * t131 + t112 * t134)) * MDP(34) + (g(1) * t96 - g(2) * t180 - g(3) * (-t105 * t134 - t112 * t131)) * MDP(35);];
taug  = t1;
