% Calculate Gravitation load on the joints for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:56
% EndTime: 2019-03-09 15:53:59
% DurationCPUTime: 0.87s
% Computational Cost: add. (428->130), mult. (791->209), div. (0->0), fcn. (907->12), ass. (0->66)
t178 = MDP(10) - MDP(18) - MDP(20);
t136 = sin(qJ(2));
t137 = sin(qJ(1));
t140 = cos(qJ(2));
t141 = cos(qJ(1));
t173 = cos(pkin(6));
t152 = t141 * t173;
t112 = t136 * t137 - t140 * t152;
t134 = sin(qJ(6));
t138 = cos(qJ(6));
t113 = t136 * t152 + t137 * t140;
t131 = qJ(3) + pkin(11);
t128 = sin(t131);
t129 = cos(t131);
t132 = sin(pkin(6));
t163 = t132 * t141;
t94 = t113 * t128 + t129 * t163;
t177 = t112 * t138 + t134 * t94;
t176 = -t112 * t134 + t138 * t94;
t174 = g(3) * t132;
t153 = t137 * t173;
t115 = -t136 * t153 + t140 * t141;
t135 = sin(qJ(3));
t170 = t115 * t135;
t169 = t128 * t134;
t168 = t128 * t138;
t167 = t132 * t136;
t166 = t132 * t137;
t139 = cos(qJ(3));
t165 = t132 * t139;
t164 = t132 * t140;
t133 = -qJ(4) - pkin(9);
t162 = t133 * t136;
t161 = t134 * t140;
t160 = t138 * t140;
t127 = pkin(3) * t139 + pkin(2);
t159 = -t112 * t127 - t113 * t133;
t114 = t141 * t136 + t140 * t153;
t158 = -t114 * t127 - t115 * t133;
t157 = t135 * t166;
t156 = t135 * t167;
t155 = t137 * t165;
t154 = t139 * t163;
t122 = t135 * t163;
t151 = t173 * t139;
t150 = t113 * t139 - t122;
t148 = pkin(4) * t129 + qJ(5) * t128;
t147 = t141 * pkin(1) + pkin(3) * t157 + pkin(8) * t166 - t114 * t133 + t115 * t127;
t97 = -t113 * t129 + t128 * t163;
t146 = t113 * t135 + t154;
t145 = -pkin(1) * t137 + pkin(3) * t122 + pkin(8) * t163 + t112 * t133 - t113 * t127;
t143 = -g(1) * t114 - g(2) * t112 + g(3) * t164;
t142 = g(1) * t115 + g(2) * t113 + g(3) * t167;
t126 = pkin(3) * t151;
t119 = pkin(3) * t155;
t116 = t127 * t164;
t107 = t173 * t128 + t129 * t167;
t106 = t128 * t167 - t173 * t129;
t101 = t115 * t139 + t157;
t100 = t155 - t170;
t99 = t115 * t129 + t128 * t166;
t98 = t115 * t128 - t129 * t166;
t90 = t114 * t138 + t134 * t98;
t89 = -t114 * t134 + t138 * t98;
t88 = -g(1) * t98 - g(2) * t94 - g(3) * t106;
t1 = [(g(1) * t137 - g(2) * t141) * MDP(2) + (g(1) * t141 + g(2) * t137) * MDP(3) + (g(1) * t113 - g(2) * t115) * MDP(9) + (g(1) * t150 - g(2) * t101) * MDP(16) + (-g(1) * t146 - g(2) * t100) * MDP(17) + (-g(1) * t145 - g(2) * t147) * MDP(19) + (g(1) * t97 + g(2) * t99) * MDP(21) + (g(1) * t94 - g(2) * t98) * MDP(22) + (-g(1) * (pkin(4) * t97 - qJ(5) * t94 + t145) - g(2) * (pkin(4) * t99 + qJ(5) * t98 + t147)) * MDP(23) + (g(1) * t177 - g(2) * t90) * MDP(29) + (g(1) * t176 - g(2) * t89) * MDP(30) - t178 * (g(1) * t112 - g(2) * t114); (-g(1) * t158 - g(2) * t159 - g(3) * (-t132 * t162 + t116)) * MDP(19) + (-g(1) * (-t114 * t148 + t158) - g(2) * (-t112 * t148 + t159) - g(3) * t116 - (t140 * t148 - t162) * t174) * MDP(23) + (-g(1) * (-t114 * t169 + t115 * t138) - g(2) * (-t112 * t169 + t113 * t138) - (t128 * t161 + t136 * t138) * t174) * MDP(29) + (-g(1) * (-t114 * t168 - t115 * t134) - g(2) * (-t112 * t168 - t113 * t134) - (t128 * t160 - t134 * t136) * t174) * MDP(30) + (-MDP(16) * t139 + MDP(17) * t135 + MDP(21) * t129 - MDP(22) * t128 - MDP(9)) * t143 + t178 * t142; (-g(1) * t100 + g(2) * t146 - g(3) * (t151 - t156)) * MDP(16) + (g(1) * t101 + g(2) * t150 - g(3) * (-t135 * t173 - t136 * t165)) * MDP(17) + (-g(1) * t119 - g(3) * t126 + (g(2) * t154 + t135 * t142) * pkin(3)) * MDP(19) + t88 * MDP(21) + (-g(1) * (-pkin(3) * t170 - pkin(4) * t98 + qJ(5) * t99 + t119) - g(2) * (-pkin(3) * t146 - t94 * pkin(4) - qJ(5) * t97) - g(3) * (-pkin(3) * t156 - pkin(4) * t106 + qJ(5) * t107 + t126)) * MDP(23) + (MDP(29) * t134 + MDP(30) * t138 + MDP(22)) * (-g(1) * t99 + g(2) * t97 - g(3) * t107); (MDP(19) + MDP(23)) * t143; t88 * MDP(23); (-g(1) * t89 - g(2) * t176 - g(3) * (t106 * t138 + t132 * t161)) * MDP(29) + (g(1) * t90 + g(2) * t177 - g(3) * (-t106 * t134 + t132 * t160)) * MDP(30);];
taug  = t1;
