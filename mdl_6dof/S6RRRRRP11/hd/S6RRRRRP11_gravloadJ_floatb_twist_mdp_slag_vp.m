% Calculate Gravitation load on the joints for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:40
% EndTime: 2019-03-10 02:53:45
% DurationCPUTime: 1.50s
% Computational Cost: add. (871->172), mult. (2356->281), div. (0->0), fcn. (2995->14), ass. (0->76)
t173 = cos(qJ(1));
t212 = cos(pkin(6));
t221 = cos(qJ(2));
t198 = t212 * t221;
t218 = sin(qJ(2));
t219 = sin(qJ(1));
t156 = -t173 * t198 + t219 * t218;
t197 = t212 * t218;
t157 = t173 * t197 + t219 * t221;
t170 = sin(qJ(3));
t209 = sin(pkin(7));
t210 = sin(pkin(6));
t191 = t210 * t209;
t189 = t173 * t191;
t211 = cos(pkin(7));
t200 = t170 * t211;
t220 = cos(qJ(3));
t132 = -t156 * t200 + t157 * t220 - t170 * t189;
t169 = sin(qJ(4));
t172 = cos(qJ(4));
t192 = t211 * t210;
t179 = t156 * t209 - t173 * t192;
t120 = t132 * t172 + t169 * t179;
t196 = t211 * t220;
t131 = t156 * t196 + t157 * t170 + t220 * t189;
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t231 = t120 * t168 - t131 * t171;
t208 = t131 * t168;
t230 = t120 * t171 + t208;
t223 = MDP(24) - MDP(32);
t119 = t132 * t169 - t172 * t179;
t178 = t173 * t218 + t219 * t198;
t225 = t178 * t211 - t219 * t191;
t224 = t221 * t192 + t209 * t212;
t158 = t173 * t221 - t219 * t197;
t135 = t158 * t170 + t225 * t220;
t207 = t135 * t168;
t206 = t168 * t172;
t205 = t171 * t172;
t204 = pkin(5) * t168 + pkin(11);
t203 = pkin(9) * t210;
t202 = pkin(10) * t209;
t201 = t169 * t209;
t199 = t172 * t209;
t195 = t221 * t210;
t194 = t210 * t218;
t136 = t158 * t220 - t225 * t170;
t174 = -t178 * t209 - t219 * t192;
t124 = t136 * t172 - t169 * t174;
t117 = -t124 * t168 + t135 * t171;
t123 = t136 * t169 + t172 * t174;
t147 = t224 * t170 + t220 * t194;
t177 = -t221 * t191 + t212 * t211;
t129 = t147 * t169 - t172 * t177;
t187 = g(1) * t123 + g(2) * t119 + g(3) * t129;
t183 = t218 * t192;
t181 = t218 * t191;
t167 = -qJ(6) - pkin(12);
t166 = t171 * pkin(5) + pkin(4);
t154 = -t170 * t183 + t220 * t195;
t153 = t170 * t195 + t220 * t183;
t146 = t170 * t194 - t224 * t220;
t142 = t154 * t172 + t169 * t181;
t141 = t154 * t169 - t172 * t181;
t140 = -t158 * t200 - t178 * t220;
t139 = t158 * t196 - t178 * t170;
t138 = -t156 * t220 - t157 * t200;
t137 = -t156 * t170 + t157 * t196;
t130 = t147 * t172 + t169 * t177;
t128 = t140 * t172 + t158 * t201;
t127 = t140 * t169 - t158 * t199;
t126 = t138 * t172 + t157 * t201;
t125 = t138 * t169 - t157 * t199;
t118 = t124 * t171 + t207;
t1 = [(g(1) * t219 - g(2) * t173) * MDP(2) + (g(1) * t173 + g(2) * t219) * MDP(3) + (g(1) * t157 - g(2) * t158) * MDP(9) + (-g(1) * t156 + g(2) * t178) * MDP(10) + (g(1) * t132 - g(2) * t136) * MDP(16) + (-g(1) * t131 + g(2) * t135) * MDP(17) + (g(1) * t120 - g(2) * t124) * MDP(23) + (g(1) * t230 - g(2) * t118) * MDP(30) + (-g(1) * t231 - g(2) * t117) * MDP(31) + (-g(1) * (-t219 * pkin(1) - t157 * pkin(2) - pkin(3) * t132 - pkin(5) * t208 - pkin(11) * t131 + t119 * t167 - t120 * t166 + t173 * t203) - g(2) * (t173 * pkin(1) + t158 * pkin(2) + t136 * pkin(3) + pkin(5) * t207 + t135 * pkin(11) - t123 * t167 + t124 * t166 + t219 * t203) + (g(1) * t179 + g(2) * t174) * pkin(10)) * MDP(33) + t223 * (-g(1) * t119 + g(2) * t123); (g(1) * t178 + g(2) * t156 - g(3) * t195) * MDP(9) + (g(1) * t158 + g(2) * t157 + g(3) * t194) * MDP(10) + (-g(1) * t140 - g(2) * t138 - g(3) * t154) * MDP(16) + (g(1) * t139 + g(2) * t137 + g(3) * t153) * MDP(17) + (-g(1) * t128 - g(2) * t126 - g(3) * t142) * MDP(23) + (-g(1) * (t128 * t171 + t139 * t168) - g(2) * (t126 * t171 + t137 * t168) - g(3) * (t142 * t171 + t153 * t168)) * MDP(30) + (-g(1) * (-t128 * t168 + t139 * t171) - g(2) * (-t126 * t168 + t137 * t171) - g(3) * (-t142 * t168 + t153 * t171)) * MDP(31) + (-g(1) * (-t178 * pkin(2) + t140 * pkin(3) - t127 * t167 + t128 * t166 + t204 * t139 + t158 * t202) - g(2) * (-t156 * pkin(2) + t138 * pkin(3) - t125 * t167 + t126 * t166 + t204 * t137 + t157 * t202) - g(3) * (pkin(2) * t195 + t154 * pkin(3) + pkin(10) * t181 - t141 * t167 + t142 * t166 + t204 * t153)) * MDP(33) + t223 * (g(1) * t127 + g(2) * t125 + g(3) * t141); (-g(1) * (-t135 * t205 + t136 * t168) - g(2) * (-t131 * t205 + t132 * t168) - g(3) * (-t146 * t205 + t147 * t168)) * MDP(30) + (-g(1) * (t135 * t206 + t136 * t171) - g(2) * (t131 * t206 + t132 * t171) - g(3) * (t146 * t206 + t147 * t171)) * MDP(31) + (-t204 * MDP(33) + MDP(17)) * (g(1) * t136 + g(2) * t132 + g(3) * t147) + (MDP(16) + t172 * MDP(23) + (t166 * t172 - t167 * t169 + pkin(3)) * MDP(33) - t223 * t169) * (g(1) * t135 + g(2) * t131 + g(3) * t146); (-g(1) * (-t123 * t166 - t124 * t167) - g(2) * (-t119 * t166 - t120 * t167) - g(3) * (-t129 * t166 - t130 * t167)) * MDP(33) + t223 * (g(1) * t124 + g(2) * t120 + g(3) * t130) + (MDP(30) * t171 - MDP(31) * t168 + MDP(23)) * t187; (g(1) * t118 + g(2) * t230 - g(3) * (-t130 * t171 - t146 * t168)) * MDP(31) + (pkin(5) * MDP(33) + MDP(30)) * (-g(1) * t117 + g(2) * t231 - g(3) * (-t130 * t168 + t146 * t171)); -t187 * MDP(33);];
taug  = t1;
