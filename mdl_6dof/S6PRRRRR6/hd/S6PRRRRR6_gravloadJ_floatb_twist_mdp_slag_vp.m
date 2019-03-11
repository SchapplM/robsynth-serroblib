% Calculate Gravitation load on the joints for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:05
% EndTime: 2019-03-09 01:25:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (1242->180), mult. (3605->327), div. (0->0), fcn. (4730->18), ass. (0->99)
t222 = cos(qJ(4));
t221 = cos(pkin(6));
t220 = cos(pkin(14));
t219 = sin(pkin(14));
t172 = sin(pkin(8));
t173 = sin(pkin(7));
t218 = t172 * t173;
t178 = sin(qJ(5));
t217 = t172 * t178;
t183 = cos(qJ(5));
t216 = t172 * t183;
t175 = cos(pkin(8));
t215 = t173 * t175;
t174 = sin(pkin(6));
t181 = sin(qJ(2));
t214 = t174 * t181;
t185 = cos(qJ(2));
t213 = t174 * t185;
t179 = sin(qJ(4));
t212 = t175 * t179;
t176 = cos(pkin(7));
t180 = sin(qJ(3));
t211 = t176 * t180;
t184 = cos(qJ(3));
t210 = t176 * t184;
t177 = sin(qJ(6));
t209 = t177 * t183;
t208 = t180 * t181;
t207 = t180 * t185;
t206 = t181 * t184;
t182 = cos(qJ(6));
t205 = t182 * t183;
t204 = t184 * t185;
t203 = t173 * t214;
t202 = t175 * t222;
t201 = t174 * t220;
t200 = t174 * t219;
t199 = t221 * t173;
t198 = t222 * t218;
t197 = t173 * t201;
t196 = t221 * t220;
t195 = t221 * t219;
t166 = -t219 * t181 + t185 * t196;
t192 = -t166 * t173 - t176 * t201;
t168 = -t220 * t181 - t185 * t195;
t191 = -t168 * t173 + t176 * t200;
t190 = -t173 * t213 + t221 * t176;
t189 = t168 * t176 + t173 * t200;
t188 = t192 * t172;
t187 = t191 * t172;
t186 = t190 * t172;
t169 = -t181 * t195 + t220 * t185;
t167 = t181 * t196 + t219 * t185;
t165 = (-t176 * t208 + t204) * t174;
t164 = (-t176 * t206 - t207) * t174;
t162 = t180 * t199 + (t176 * t207 + t206) * t174;
t161 = t184 * t199 + (t176 * t204 - t208) * t174;
t157 = -t164 * t172 + t175 * t203;
t156 = t168 * t184 - t169 * t211;
t155 = -t168 * t180 - t169 * t210;
t154 = t166 * t184 - t167 * t211;
t153 = -t166 * t180 - t167 * t210;
t152 = t169 * t184 + t180 * t189;
t151 = -t169 * t180 + t184 * t189;
t150 = t166 * t211 + t167 * t184 - t180 * t197;
t149 = -t167 * t180 + (t166 * t176 - t197) * t184;
t148 = -t161 * t172 + t175 * t190;
t145 = t165 * t222 + (t164 * t175 + t172 * t203) * t179;
t144 = -t164 * t202 + t165 * t179 - t198 * t214;
t143 = -t155 * t172 + t169 * t215;
t142 = -t153 * t172 + t167 * t215;
t141 = t161 * t222 - t162 * t212;
t140 = t161 * t179 + t162 * t202;
t139 = -t151 * t172 + t175 * t191;
t138 = -t149 * t172 + t192 * t175;
t137 = t162 * t222 + (t161 * t175 + t186) * t179;
t136 = -t161 * t202 + t162 * t179 - t222 * t186;
t135 = t141 * t183 + t162 * t217;
t134 = t145 * t183 + t157 * t178;
t133 = t151 * t222 - t152 * t212;
t132 = t151 * t179 + t152 * t202;
t131 = t149 * t222 - t150 * t212;
t130 = t149 * t179 + t150 * t202;
t129 = t156 * t222 + (t155 * t175 + t169 * t218) * t179;
t128 = -t155 * t202 + t156 * t179 - t169 * t198;
t127 = t154 * t222 + (t153 * t175 + t167 * t218) * t179;
t126 = -t153 * t202 + t154 * t179 - t167 * t198;
t125 = t152 * t222 + (t151 * t175 + t187) * t179;
t124 = -t151 * t202 + t152 * t179 - t222 * t187;
t123 = t150 * t222 + (t149 * t175 + t188) * t179;
t122 = -t149 * t202 + t150 * t179 - t222 * t188;
t121 = t137 * t183 + t148 * t178;
t119 = t133 * t183 + t152 * t217;
t118 = t131 * t183 + t150 * t217;
t117 = t129 * t183 + t143 * t178;
t116 = t127 * t183 + t142 * t178;
t115 = t125 * t183 + t139 * t178;
t113 = t123 * t183 + t138 * t178;
t1 = [-g(3) * MDP(1); (-g(1) * t168 - g(2) * t166 - g(3) * t213) * MDP(3) + (g(1) * t169 + g(2) * t167 + g(3) * t214) * MDP(4) + (-g(1) * t156 - g(2) * t154 - g(3) * t165) * MDP(10) + (-g(1) * t155 - g(2) * t153 - g(3) * t164) * MDP(11) + (-g(1) * t129 - g(2) * t127 - g(3) * t145) * MDP(17) + (g(1) * t128 + g(2) * t126 + g(3) * t144) * MDP(18) + (-g(1) * t117 - g(2) * t116 - g(3) * t134) * MDP(24) + (-g(1) * (-t129 * t178 + t143 * t183) - g(2) * (-t127 * t178 + t142 * t183) - g(3) * (-t145 * t178 + t157 * t183)) * MDP(25) + (-g(1) * (t117 * t182 + t128 * t177) - g(2) * (t116 * t182 + t126 * t177) - g(3) * (t134 * t182 + t144 * t177)) * MDP(31) + (-g(1) * (-t117 * t177 + t128 * t182) - g(2) * (-t116 * t177 + t126 * t182) - g(3) * (-t134 * t177 + t144 * t182)) * MDP(32); (-g(1) * t151 - g(2) * t149 - g(3) * t161) * MDP(10) + (g(1) * t152 + g(2) * t150 + g(3) * t162) * MDP(11) + (-g(1) * t133 - g(2) * t131 - g(3) * t141) * MDP(17) + (g(1) * t132 + g(2) * t130 + g(3) * t140) * MDP(18) + (-g(1) * t119 - g(2) * t118 - g(3) * t135) * MDP(24) + (-g(1) * (-t133 * t178 + t152 * t216) - g(2) * (-t131 * t178 + t150 * t216) - g(3) * (-t141 * t178 + t162 * t216)) * MDP(25) + (-g(1) * (t119 * t182 + t132 * t177) - g(2) * (t118 * t182 + t130 * t177) - g(3) * (t135 * t182 + t140 * t177)) * MDP(31) + (-g(1) * (-t119 * t177 + t132 * t182) - g(2) * (-t118 * t177 + t130 * t182) - g(3) * (-t135 * t177 + t140 * t182)) * MDP(32); (g(1) * t125 + g(2) * t123 + g(3) * t137) * MDP(18) + (-g(1) * (-t124 * t205 + t125 * t177) - g(2) * (-t122 * t205 + t123 * t177) - g(3) * (-t136 * t205 + t137 * t177)) * MDP(31) + (-g(1) * (t124 * t209 + t125 * t182) - g(2) * (t122 * t209 + t123 * t182) - g(3) * (t136 * t209 + t137 * t182)) * MDP(32) + (t183 * MDP(24) - MDP(25) * t178 + MDP(17)) * (g(1) * t124 + g(2) * t122 + g(3) * t136); (g(1) * t115 + g(2) * t113 + g(3) * t121) * MDP(25) + (-MDP(31) * t182 + MDP(32) * t177 - MDP(24)) * (g(1) * (-t125 * t178 + t139 * t183) + g(2) * (-t123 * t178 + t138 * t183) + g(3) * (-t137 * t178 + t148 * t183)); (-g(1) * (-t115 * t177 + t124 * t182) - g(2) * (-t113 * t177 + t122 * t182) - g(3) * (-t121 * t177 + t136 * t182)) * MDP(31) + (-g(1) * (-t115 * t182 - t124 * t177) - g(2) * (-t113 * t182 - t122 * t177) - g(3) * (-t121 * t182 - t136 * t177)) * MDP(32);];
taug  = t1;
