% Calculate Gravitation load on the joints for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:27
% EndTime: 2019-03-09 12:05:29
% DurationCPUTime: 0.90s
% Computational Cost: add. (473->119), mult. (1161->192), div. (0->0), fcn. (1437->12), ass. (0->65)
t185 = MDP(19) - MDP(27);
t134 = sin(pkin(11));
t140 = sin(qJ(2));
t144 = cos(qJ(2));
t178 = cos(pkin(11));
t124 = -t144 * t134 - t140 * t178;
t141 = sin(qJ(1));
t145 = cos(qJ(1));
t136 = cos(pkin(6));
t158 = t144 * t178;
t151 = -t134 * t140 + t158;
t147 = t151 * t136;
t105 = t141 * t124 + t145 * t147;
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t117 = t124 * t136;
t106 = -t117 * t145 + t141 * t151;
t139 = sin(qJ(4));
t143 = cos(qJ(4));
t135 = sin(pkin(6));
t169 = t135 * t145;
t99 = t106 * t143 - t139 * t169;
t189 = t105 * t142 + t138 * t99;
t188 = -t105 * t138 + t142 * t99;
t165 = t141 * t144;
t167 = t140 * t145;
t121 = -t136 * t165 - t167;
t170 = t135 * t144;
t186 = -g(1) * t121 - g(3) * t170;
t108 = t124 * t145 - t141 * t147;
t172 = t135 * t140;
t115 = t134 * t172 - t135 * t158;
t184 = -g(1) * t108 - g(2) * t105 + g(3) * t115;
t177 = t106 * t138;
t107 = -t141 * t117 - t145 * t151;
t174 = t107 * t138;
t116 = t124 * t135;
t173 = t116 * t138;
t171 = t135 * t141;
t168 = t138 * t143;
t166 = t141 * t140;
t164 = t142 * t143;
t163 = t144 * t145;
t161 = t136 * t163;
t159 = pkin(5) * t138 + pkin(9);
t118 = pkin(2) * t136 * t140 + (-pkin(8) - qJ(3)) * t135;
t133 = pkin(2) * t144 + pkin(1);
t157 = -t118 * t141 + t145 * t133;
t154 = g(1) * t141 - g(2) * t145;
t103 = -t107 * t143 + t139 * t171;
t96 = -t103 * t138 - t108 * t142;
t153 = -t118 * t145 - t141 * t133;
t98 = t106 * t139 + t143 * t169;
t102 = -t107 * t139 - t143 * t171;
t110 = -t116 * t139 - t136 * t143;
t150 = g(1) * t102 + g(2) * t98 + g(3) * t110;
t137 = -qJ(6) - pkin(10);
t132 = pkin(5) * t142 + pkin(4);
t125 = pkin(2) * t161;
t122 = -t136 * t166 + t163;
t120 = -t136 * t167 - t165;
t119 = -t161 + t166;
t111 = -t116 * t143 + t136 * t139;
t97 = t103 * t142 - t108 * t138;
t1 = [t154 * MDP(2) + (-g(1) * t120 - g(2) * t122) * MDP(9) + (-g(1) * t119 - g(2) * t121) * MDP(10) + (-g(1) * t153 - g(2) * t157) * MDP(12) + (g(1) * t99 - g(2) * t103) * MDP(18) + (g(1) * t188 - g(2) * t97) * MDP(25) + (-g(1) * t189 - g(2) * t96) * MDP(26) + (-g(1) * (-t106 * pkin(3) + t159 * t105 - t132 * t99 + t137 * t98 + t153) - g(2) * (-pkin(3) * t107 - t102 * t137 + t103 * t132 - t159 * t108 + t157)) * MDP(28) + t185 * (-g(1) * t98 + g(2) * t102) + (-MDP(11) * t135 + MDP(3)) * (g(1) * t145 + g(2) * t141); (g(2) * t119 + t186) * MDP(9) + (g(1) * t122 - g(2) * t120 + g(3) * t172) * MDP(10) + (-g(2) * t125 + (g(2) * t166 + t186) * pkin(2)) * MDP(12) + (-g(1) * (t108 * t164 - t174) - g(2) * (t105 * t164 + t177) - g(3) * (-t115 * t164 - t173)) * MDP(25) + (-g(1) * (-t107 * t142 - t108 * t168) - g(2) * (-t105 * t168 + t106 * t142) - g(3) * (t115 * t168 - t116 * t142)) * MDP(26) + (-g(1) * (t121 * pkin(2) - pkin(5) * t174 - t107 * pkin(9)) - g(2) * (-pkin(2) * t166 + pkin(5) * t177 + pkin(9) * t106 + t125) - g(3) * (pkin(2) * t170 - pkin(5) * t173 - t116 * pkin(9)) + t184 * (t132 * t143 - t137 * t139 + pkin(3))) * MDP(28) + (t143 * MDP(18) - t185 * t139) * t184; (MDP(12) + MDP(28)) * (-g(3) * t136 - t154 * t135); (-g(1) * (-t102 * t132 - t103 * t137) - g(2) * (-t98 * t132 - t137 * t99) - g(3) * (-t110 * t132 - t111 * t137)) * MDP(28) + t185 * (g(1) * t103 + g(2) * t99 + g(3) * t111) + (MDP(25) * t142 - MDP(26) * t138 + MDP(18)) * t150; (g(1) * t97 + g(2) * t188 - g(3) * (-t111 * t142 - t115 * t138)) * MDP(26) + (pkin(5) * MDP(28) + MDP(25)) * (g(2) * t189 - g(3) * (-t111 * t138 + t115 * t142) - g(1) * t96); -t150 * MDP(28);];
taug  = t1;
