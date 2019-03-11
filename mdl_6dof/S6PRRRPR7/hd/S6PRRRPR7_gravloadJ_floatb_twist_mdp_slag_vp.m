% Calculate Gravitation load on the joints for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:39
% EndTime: 2019-03-08 23:44:42
% DurationCPUTime: 0.93s
% Computational Cost: add. (772->158), mult. (2048->270), div. (0->0), fcn. (2612->16), ass. (0->72)
t185 = cos(pkin(12));
t187 = cos(pkin(6));
t170 = t187 * t185;
t182 = sin(pkin(12));
t191 = sin(qJ(2));
t193 = cos(qJ(2));
t152 = -t170 * t193 + t182 * t191;
t183 = sin(pkin(7));
t184 = sin(pkin(6));
t166 = t184 * t183;
t186 = cos(pkin(7));
t197 = t152 * t186 + t166 * t185;
t168 = t187 * t182;
t153 = t168 * t193 + t185 * t191;
t165 = t184 * t182;
t196 = t153 * t186 - t165 * t183;
t167 = t186 * t184;
t195 = t167 * t193 + t183 * t187;
t194 = MDP(18) - MDP(21);
t192 = cos(qJ(3));
t141 = pkin(13) + qJ(6);
t139 = sin(t141);
t146 = cos(qJ(4));
t181 = t139 * t146;
t140 = cos(t141);
t180 = t140 * t146;
t142 = sin(pkin(13));
t179 = t142 * t146;
t143 = cos(pkin(13));
t178 = t143 * t146;
t177 = pkin(9) * t183;
t144 = sin(qJ(4));
t176 = t144 * t183;
t145 = sin(qJ(3));
t175 = t145 * t186;
t174 = t146 * t183;
t173 = t186 * t192;
t172 = t193 * t184;
t171 = t184 * t191;
t131 = t170 * t191 + t182 * t193;
t109 = t131 * t192 - t145 * t197;
t147 = t152 * t183 - t167 * t185;
t100 = t109 * t144 - t146 * t147;
t132 = -t168 * t191 + t185 * t193;
t111 = t132 * t192 - t145 * t196;
t148 = t153 * t183 + t165 * t186;
t102 = t111 * t144 - t146 * t148;
t123 = t145 * t195 + t171 * t192;
t151 = -t166 * t193 + t186 * t187;
t112 = t123 * t144 - t146 * t151;
t163 = g(1) * t102 + g(2) * t100 + g(3) * t112;
t157 = t191 * t167;
t156 = t191 * t166;
t129 = -t145 * t157 + t172 * t192;
t128 = t145 * t172 + t157 * t192;
t122 = t145 * t171 - t192 * t195;
t119 = t129 * t146 + t144 * t156;
t118 = t129 * t144 - t146 * t156;
t117 = -t132 * t175 - t153 * t192;
t116 = t132 * t173 - t145 * t153;
t115 = -t131 * t175 - t152 * t192;
t114 = t131 * t173 - t145 * t152;
t113 = t123 * t146 + t144 * t151;
t110 = t132 * t145 + t192 * t196;
t108 = t131 * t145 + t192 * t197;
t107 = t117 * t146 + t132 * t176;
t106 = t117 * t144 - t132 * t174;
t105 = t115 * t146 + t131 * t176;
t104 = t115 * t144 - t131 * t174;
t103 = t111 * t146 + t144 * t148;
t101 = t109 * t146 + t144 * t147;
t1 = [(-MDP(1) - MDP(22)) * g(3); (g(1) * t153 + g(2) * t152 - g(3) * t172) * MDP(3) + (g(1) * t132 + g(2) * t131 + g(3) * t171) * MDP(4) + (-g(1) * t117 - g(2) * t115 - g(3) * t129) * MDP(10) + (g(1) * t116 + g(2) * t114 + g(3) * t128) * MDP(11) + (-g(1) * t107 - g(2) * t105 - g(3) * t119) * MDP(17) + (-g(1) * (t107 * t143 + t116 * t142) - g(2) * (t105 * t143 + t114 * t142) - g(3) * (t119 * t143 + t128 * t142)) * MDP(19) + (-g(1) * (-t107 * t142 + t116 * t143) - g(2) * (-t105 * t142 + t114 * t143) - g(3) * (-t119 * t142 + t128 * t143)) * MDP(20) + (-g(1) * (-pkin(2) * t153 + t117 * pkin(3) + t107 * pkin(4) + t116 * pkin(10) + t106 * qJ(5) + t132 * t177) - g(2) * (-pkin(2) * t152 + t115 * pkin(3) + t105 * pkin(4) + t114 * pkin(10) + t104 * qJ(5) + t131 * t177) - g(3) * (pkin(2) * t172 + pkin(3) * t129 + pkin(4) * t119 + pkin(9) * t156 + pkin(10) * t128 + qJ(5) * t118)) * MDP(22) + (-g(1) * (t107 * t140 + t116 * t139) - g(2) * (t105 * t140 + t114 * t139) - g(3) * (t119 * t140 + t128 * t139)) * MDP(28) + (-g(1) * (-t107 * t139 + t116 * t140) - g(2) * (-t105 * t139 + t114 * t140) - g(3) * (-t119 * t139 + t128 * t140)) * MDP(29) + t194 * (g(1) * t106 + g(2) * t104 + g(3) * t118); (-g(1) * (-t110 * t178 + t111 * t142) - g(2) * (-t108 * t178 + t109 * t142) - g(3) * (-t122 * t178 + t123 * t142)) * MDP(19) + (-g(1) * (t110 * t179 + t111 * t143) - g(2) * (t108 * t179 + t109 * t143) - g(3) * (t122 * t179 + t123 * t143)) * MDP(20) + (-g(1) * (-t110 * t180 + t111 * t139) - g(2) * (-t108 * t180 + t109 * t139) - g(3) * (-t122 * t180 + t123 * t139)) * MDP(28) + (-g(1) * (t110 * t181 + t111 * t140) - g(2) * (t108 * t181 + t109 * t140) - g(3) * (t122 * t181 + t123 * t140)) * MDP(29) + (-MDP(22) * pkin(10) + MDP(11)) * (g(1) * t111 + g(2) * t109 + g(3) * t123) + (-t194 * t144 + MDP(22) * (pkin(4) * t146 + qJ(5) * t144 + pkin(3)) + MDP(17) * t146 + MDP(10)) * (g(1) * t110 + g(2) * t108 + g(3) * t122); (-g(1) * (-pkin(4) * t102 + qJ(5) * t103) - g(2) * (-pkin(4) * t100 + qJ(5) * t101) - g(3) * (-pkin(4) * t112 + qJ(5) * t113)) * MDP(22) + t194 * (g(1) * t103 + g(2) * t101 + g(3) * t113) + (MDP(19) * t143 - MDP(20) * t142 + MDP(28) * t140 - MDP(29) * t139 + MDP(17)) * t163; -t163 * MDP(22); (-g(1) * (-t103 * t139 + t110 * t140) - g(2) * (-t101 * t139 + t108 * t140) - g(3) * (-t113 * t139 + t122 * t140)) * MDP(28) + (-g(1) * (-t103 * t140 - t110 * t139) - g(2) * (-t101 * t140 - t108 * t139) - g(3) * (-t113 * t140 - t122 * t139)) * MDP(29);];
taug  = t1;
