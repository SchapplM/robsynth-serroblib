% Calculate Gravitation load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:27
% EndTime: 2019-03-08 21:16:29
% DurationCPUTime: 0.68s
% Computational Cost: add. (392->111), mult. (1033->183), div. (0->0), fcn. (1253->12), ass. (0->60)
t145 = sin(qJ(3));
t187 = qJ(4) * t145 + pkin(2);
t186 = MDP(11) - MDP(14) - MDP(17);
t185 = MDP(12) + MDP(16);
t184 = MDP(13) - MDP(18);
t146 = sin(qJ(2));
t148 = cos(qJ(2));
t177 = cos(pkin(10));
t178 = cos(pkin(6));
t159 = t178 * t177;
t176 = sin(pkin(10));
t126 = t146 * t159 + t176 * t148;
t142 = sin(pkin(6));
t166 = t142 * t177;
t182 = cos(qJ(3));
t110 = t126 * t145 + t182 * t166;
t158 = t178 * t176;
t128 = -t146 * t158 + t177 * t148;
t165 = t142 * t176;
t112 = t128 * t145 - t182 * t165;
t174 = t142 * t146;
t129 = t145 * t174 - t178 * t182;
t151 = g(1) * t112 + g(2) * t110 + g(3) * t129;
t173 = t142 * t148;
t172 = MDP(15) + MDP(19);
t125 = t176 * t146 - t148 * t159;
t171 = t125 * t182;
t127 = t177 * t146 + t148 * t158;
t170 = t127 * t182;
t141 = sin(pkin(11));
t169 = t141 * t182;
t143 = cos(pkin(11));
t168 = t143 * t182;
t167 = t148 * t182;
t111 = t126 * t182 - t145 * t166;
t164 = -t110 * pkin(3) + qJ(4) * t111;
t113 = t128 * t182 + t145 * t165;
t163 = -t112 * pkin(3) + qJ(4) * t113;
t130 = t178 * t145 + t182 * t174;
t162 = -t129 * pkin(3) + qJ(4) * t130;
t160 = t142 * t167;
t161 = pkin(3) * t160 + pkin(8) * t174 + t187 * t173;
t157 = -pkin(4) * t143 - qJ(5) * t141;
t154 = -pkin(3) * t171 + pkin(8) * t126 - t187 * t125;
t153 = -pkin(3) * t170 + pkin(8) * t128 - t187 * t127;
t147 = cos(qJ(6));
t144 = sin(qJ(6));
t115 = (t141 * t146 + t143 * t167) * t142;
t114 = t141 * t160 - t143 * t174;
t109 = t130 * t143 - t141 * t173;
t108 = t130 * t141 + t143 * t173;
t105 = -t127 * t168 + t128 * t141;
t104 = -t127 * t169 - t128 * t143;
t103 = -t125 * t168 + t126 * t141;
t102 = -t125 * t169 - t126 * t143;
t100 = t113 * t143 + t127 * t141;
t99 = t113 * t141 - t127 * t143;
t98 = t111 * t143 + t125 * t141;
t97 = t111 * t141 - t125 * t143;
t1 = [(-MDP(1) - t172) * g(3); (g(1) * t128 + g(2) * t126 + g(3) * t174) * MDP(4) + (g(1) * t170 + g(2) * t171 - g(3) * t160) * MDP(10) + (-g(1) * t153 - g(2) * t154 - g(3) * t161) * MDP(15) + (-g(1) * (pkin(4) * t105 + qJ(5) * t104 + t153) - g(2) * (pkin(4) * t103 + qJ(5) * t102 + t154) - g(3) * (pkin(4) * t115 + qJ(5) * t114 + t161)) * MDP(19) + (-g(1) * (t104 * t144 + t105 * t147) - g(2) * (t102 * t144 + t103 * t147) - g(3) * (t114 * t144 + t115 * t147)) * MDP(25) + (-g(1) * (t104 * t147 - t105 * t144) - g(2) * (t102 * t147 - t103 * t144) - g(3) * (t114 * t147 - t115 * t144)) * MDP(26) + t185 * (-g(1) * t105 - g(2) * t103 - g(3) * t115) + t184 * (g(1) * t104 + g(2) * t102 + g(3) * t114) + (t186 * t145 - MDP(3)) * (-g(1) * t127 - g(2) * t125 + g(3) * t173); (-g(1) * t163 - g(2) * t164 - g(3) * t162) * MDP(15) + (-g(1) * (t157 * t112 + t163) - g(2) * (t157 * t110 + t164) - g(3) * (t157 * t129 + t162)) * MDP(19) + t186 * (g(1) * t113 + g(2) * t111 + g(3) * t130) + (MDP(10) + t185 * t143 - t184 * t141 + MDP(26) * (t141 * t147 - t143 * t144) + MDP(25) * (t141 * t144 + t143 * t147)) * t151; -t172 * t151; (-g(1) * t99 - g(2) * t97 - g(3) * t108) * MDP(19); (-g(1) * (-t100 * t144 + t147 * t99) - g(2) * (-t144 * t98 + t147 * t97) - g(3) * (t108 * t147 - t109 * t144)) * MDP(25) + (-g(1) * (-t100 * t147 - t144 * t99) - g(2) * (-t144 * t97 - t147 * t98) - g(3) * (-t108 * t144 - t109 * t147)) * MDP(26);];
taug  = t1;
