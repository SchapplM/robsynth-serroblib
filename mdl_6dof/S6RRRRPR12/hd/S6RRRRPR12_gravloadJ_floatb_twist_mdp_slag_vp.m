% Calculate Gravitation load on the joints for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:43
% EndTime: 2019-03-09 23:44:48
% DurationCPUTime: 1.56s
% Computational Cost: add. (710->170), mult. (1786->292), div. (0->0), fcn. (2253->16), ass. (0->76)
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t168 = cos(qJ(1));
t198 = cos(pkin(6));
t178 = t168 * t198;
t202 = sin(qJ(1));
t144 = t164 * t202 - t167 * t178;
t145 = t164 * t178 + t167 * t202;
t163 = sin(qJ(3));
t197 = cos(pkin(7));
t179 = t163 * t197;
t158 = sin(pkin(7));
t159 = sin(pkin(6));
t187 = t159 * t168;
t184 = t158 * t187;
t203 = cos(qJ(3));
t121 = -t144 * t179 + t145 * t203 - t163 * t184;
t180 = t159 * t197;
t134 = t144 * t158 - t168 * t180;
t157 = qJ(4) + pkin(13);
t155 = sin(t157);
t156 = cos(t157);
t110 = t121 * t156 + t134 * t155;
t175 = t197 * t203;
t120 = t144 * t175 + t145 * t163 + t184 * t203;
t161 = sin(qJ(6));
t165 = cos(qJ(6));
t213 = t110 * t161 - t120 * t165;
t212 = t110 * t165 + t120 * t161;
t166 = cos(qJ(4));
t162 = sin(qJ(4));
t196 = t134 * t162;
t211 = t121 * t166 + t196;
t205 = MDP(17) - MDP(25);
t207 = t121 * t162 - t134 * t166;
t176 = t198 * t202;
t147 = -t164 * t176 + t168 * t167;
t189 = t159 * t164;
t204 = g(1) * t147 + g(2) * t145 + g(3) * t189;
t146 = -t168 * t164 - t167 * t176;
t135 = -t146 * t158 + t180 * t202;
t195 = t135 * t162;
t194 = t155 * t158;
t193 = t156 * t161;
t192 = t156 * t165;
t191 = t158 * t162;
t190 = t158 * t166;
t188 = t159 * t167;
t185 = t158 * t189;
t182 = t159 * t202;
t181 = t158 * t198;
t177 = t158 * t182;
t125 = t147 * t203 + (t146 * t197 + t177) * t163;
t114 = -t125 * t162 + t135 * t166;
t124 = -t146 * t175 + t147 * t163 - t177 * t203;
t131 = t163 * t189 - t175 * t188 - t181 * t203;
t172 = g(1) * t124 + g(2) * t120 + g(3) * t131;
t160 = -qJ(5) - pkin(11);
t154 = t166 * pkin(4) + pkin(3);
t143 = -t158 * t188 + t197 * t198;
t142 = (-t164 * t179 + t167 * t203) * t159;
t141 = (t163 * t167 + t164 * t175) * t159;
t132 = t163 * t181 + (t164 * t203 + t167 * t179) * t159;
t130 = t146 * t203 - t147 * t179;
t129 = t146 * t163 + t147 * t175;
t128 = -t144 * t203 - t145 * t179;
t127 = -t144 * t163 + t145 * t175;
t126 = t142 * t156 + t155 * t185;
t119 = t132 * t156 + t143 * t155;
t117 = t130 * t156 + t147 * t194;
t116 = t128 * t156 + t145 * t194;
t115 = t125 * t166 + t195;
t113 = t125 * t156 + t135 * t155;
t108 = t113 * t165 + t124 * t161;
t107 = -t113 * t161 + t124 * t165;
t1 = [(g(1) * t202 - g(2) * t168) * MDP(2) + (g(1) * t168 + g(2) * t202) * MDP(3) + (g(1) * t145 - g(2) * t147) * MDP(9) + (-g(1) * t144 - g(2) * t146) * MDP(10) + (g(1) * t121 - g(2) * t125) * MDP(16) + (g(1) * t211 - g(2) * t115) * MDP(23) + (-g(1) * t207 - g(2) * t114) * MDP(24) + (-g(1) * (-pkin(1) * t202 - t145 * pkin(2) - pkin(4) * t196 + pkin(9) * t187 + t120 * t160 - t121 * t154) - g(2) * (t168 * pkin(1) + t147 * pkin(2) + pkin(4) * t195 + pkin(9) * t182 - t124 * t160 + t125 * t154) + (g(1) * t134 - g(2) * t135) * pkin(10)) * MDP(26) + (g(1) * t212 - g(2) * t108) * MDP(32) + (-g(1) * t213 - g(2) * t107) * MDP(33) + t205 * (-g(1) * t120 + g(2) * t124); (-g(1) * t146 + g(2) * t144 - g(3) * t188) * MDP(9) + t204 * MDP(10) + (-g(1) * t130 - g(2) * t128 - g(3) * t142) * MDP(16) + (-g(1) * (t130 * t166 + t147 * t191) - g(2) * (t128 * t166 + t145 * t191) - g(3) * (t142 * t166 + t162 * t185)) * MDP(23) + (-g(1) * (-t130 * t162 + t147 * t190) - g(2) * (-t128 * t162 + t145 * t190) - g(3) * (-t142 * t162 + t166 * t185)) * MDP(24) + (-g(1) * (t146 * pkin(2) - t129 * t160 + t130 * t154) - g(2) * (-t144 * pkin(2) - t127 * t160 + t128 * t154) - g(3) * (pkin(2) * t188 - t141 * t160 + t142 * t154) - t204 * t158 * (pkin(4) * t162 + pkin(10))) * MDP(26) + (-g(1) * (t117 * t165 + t129 * t161) - g(2) * (t116 * t165 + t127 * t161) - g(3) * (t126 * t165 + t141 * t161)) * MDP(32) + (-g(1) * (-t117 * t161 + t129 * t165) - g(2) * (-t116 * t161 + t127 * t165) - g(3) * (-t126 * t161 + t141 * t165)) * MDP(33) + t205 * (g(1) * t129 + g(2) * t127 + g(3) * t141); (-g(1) * (-t124 * t154 - t125 * t160) - g(2) * (-t120 * t154 - t121 * t160) - g(3) * (-t131 * t154 - t132 * t160)) * MDP(26) + (-g(1) * (-t124 * t192 + t125 * t161) - g(2) * (-t120 * t192 + t121 * t161) - g(3) * (-t131 * t192 + t132 * t161)) * MDP(32) + (-g(1) * (t124 * t193 + t125 * t165) - g(2) * (t120 * t193 + t121 * t165) - g(3) * (t131 * t193 + t132 * t165)) * MDP(33) + t205 * (g(1) * t125 + g(2) * t121 + g(3) * t132) + (MDP(23) * t166 - MDP(24) * t162 + MDP(16)) * t172; (g(1) * t115 + g(2) * t211 - g(3) * (-t132 * t166 - t143 * t162)) * MDP(24) + (-MDP(32) * t165 + MDP(33) * t161) * (g(1) * (-t125 * t155 + t135 * t156) + g(2) * (-t121 * t155 + t134 * t156) + g(3) * (-t132 * t155 + t143 * t156)) + (pkin(4) * MDP(26) + MDP(23)) * (-g(1) * t114 + g(2) * t207 - g(3) * (-t132 * t162 + t143 * t166)); -t172 * MDP(26); (-g(1) * t107 + g(2) * t213 - g(3) * (-t119 * t161 + t131 * t165)) * MDP(32) + (g(1) * t108 + g(2) * t212 - g(3) * (-t119 * t165 - t131 * t161)) * MDP(33);];
taug  = t1;
