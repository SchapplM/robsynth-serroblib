% Calculate Gravitation load on the joints for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:20
% EndTime: 2019-03-09 19:10:25
% DurationCPUTime: 1.59s
% Computational Cost: add. (749->169), mult. (2024->305), div. (0->0), fcn. (2603->16), ass. (0->84)
t163 = sin(pkin(13));
t169 = sin(qJ(3));
t174 = cos(qJ(3));
t206 = cos(pkin(13));
t157 = -t174 * t163 - t169 * t206;
t164 = sin(pkin(7));
t146 = t157 * t164;
t166 = cos(pkin(7));
t148 = t157 * t166;
t170 = sin(qJ(2));
t171 = sin(qJ(1));
t175 = cos(qJ(2));
t176 = cos(qJ(1));
t207 = cos(pkin(6));
t190 = t176 * t207;
t152 = t170 * t171 - t175 * t190;
t153 = t170 * t190 + t171 * t175;
t184 = -t169 * t163 + t174 * t206;
t165 = sin(pkin(6));
t200 = t165 * t176;
t125 = -t146 * t200 - t152 * t148 - t153 * t184;
t141 = -t152 * t164 + t166 * t200;
t168 = sin(qJ(5));
t173 = cos(qJ(5));
t114 = t125 * t173 + t141 * t168;
t167 = sin(qJ(6));
t172 = cos(qJ(6));
t145 = t184 * t164;
t147 = t184 * t166;
t185 = -t145 * t200 - t147 * t152 + t153 * t157;
t218 = t114 * t167 - t172 * t185;
t217 = t114 * t172 + t167 * t185;
t216 = t125 * t168 - t141 * t173;
t213 = g(1) * t176 + g(2) * t171;
t210 = g(3) * t165;
t208 = pkin(10) + qJ(4);
t205 = t164 * t168;
t204 = t164 * t173;
t203 = t165 * t170;
t202 = t165 * t171;
t201 = t165 * t175;
t199 = t166 * t169;
t198 = t166 * t174;
t197 = t167 * t173;
t196 = t169 * t175;
t195 = t170 * t174;
t194 = t172 * t173;
t193 = t164 * t203;
t192 = t164 * t200;
t191 = t171 * t207;
t189 = t207 * t164;
t188 = g(2) * t153 + g(3) * t203;
t187 = t152 * t166 + t192;
t154 = -t176 * t170 - t175 * t191;
t186 = t154 * t166 + t164 * t202;
t181 = t152 * t199 - t153 * t174 + t169 * t192;
t155 = -t170 * t191 + t175 * t176;
t180 = g(1) * t155 + t188;
t179 = -t146 * t202 - t148 * t154 + t155 * t184;
t178 = -t207 * t146 + (-t148 * t175 + t170 * t184) * t165;
t177 = g(2) * t187 - g(3) * (t166 * t201 + t189);
t162 = pkin(3) * t174 + pkin(2);
t151 = -t164 * t201 + t166 * t207;
t150 = pkin(3) * t199 - t164 * t208;
t143 = -t154 * t164 + t166 * t202;
t140 = (t148 * t170 + t175 * t184) * t165;
t139 = (t147 * t170 - t157 * t175) * t165;
t138 = t155 * t174 + t169 * t186;
t137 = -t155 * t169 + t174 * t186;
t136 = t140 * t173 + t168 * t193;
t135 = t148 * t155 + t154 * t184;
t134 = t147 * t155 - t154 * t157;
t133 = t148 * t153 - t152 * t184;
t132 = t147 * t153 + t152 * t157;
t130 = t207 * t145 + (t147 * t175 + t157 * t170) * t165;
t127 = t145 * t202 + t147 * t154 + t155 * t157;
t120 = t135 * t173 + t155 * t205;
t119 = t133 * t173 + t153 * t205;
t118 = t151 * t168 + t173 * t178;
t116 = t143 * t168 + t173 * t179;
t115 = t143 * t173 - t168 * t179;
t111 = t116 * t172 - t127 * t167;
t110 = -t116 * t167 - t127 * t172;
t1 = [(g(1) * t171 - g(2) * t176) * MDP(2) + t213 * MDP(3) + (g(1) * t153 - g(2) * t155) * MDP(9) + (-g(1) * t152 - g(2) * t154) * MDP(10) + (-g(1) * t181 - g(2) * t138) * MDP(16) + (-g(1) * (t153 * t169 + t174 * t187) - g(2) * t137) * MDP(17) + (-g(1) * t141 - g(2) * t143) * MDP(18) + (-g(1) * (-t171 * pkin(1) + t152 * t150 - t153 * t162) - g(2) * (pkin(1) * t176 + t154 * t150 + t155 * t162) - t213 * t165 * (pkin(3) * t164 * t169 + t166 * t208 + pkin(9))) * MDP(19) + (-g(1) * t114 - g(2) * t116) * MDP(25) + (g(1) * t216 - g(2) * t115) * MDP(26) + (-g(1) * t217 - g(2) * t111) * MDP(32) + (g(1) * t218 - g(2) * t110) * MDP(33); (-g(1) * t154 + g(2) * t152 - g(3) * t201) * MDP(9) + (-g(1) * (t154 * t174 - t155 * t199) - g(2) * (-t152 * t174 - t153 * t199) - (-t170 * t199 + t174 * t175) * t210) * MDP(16) + (-g(1) * (-t154 * t169 - t155 * t198) - g(2) * (t152 * t169 - t153 * t198) - (-t166 * t195 - t196) * t210) * MDP(17) + (-g(1) * (-t150 * t155 + t154 * t162) - g(2) * (-t150 * t153 - t152 * t162) - (-t150 * t170 + t162 * t175) * t210) * MDP(19) + (-g(1) * t120 - g(2) * t119 - g(3) * t136) * MDP(25) + (-g(1) * (-t135 * t168 + t155 * t204) - g(2) * (-t133 * t168 + t153 * t204) - g(3) * (-t140 * t168 + t173 * t193)) * MDP(26) + (-g(1) * (t120 * t172 + t134 * t167) - g(2) * (t119 * t172 + t132 * t167) - g(3) * (t136 * t172 + t139 * t167)) * MDP(32) + (-g(1) * (-t120 * t167 + t134 * t172) - g(2) * (-t119 * t167 + t132 * t172) - g(3) * (-t136 * t167 + t139 * t172)) * MDP(33) + (-MDP(18) * t164 + MDP(10)) * t180; (-g(1) * t137 + t169 * t188 + t174 * t177) * MDP(16) + (g(1) * t138 - g(2) * t181 - g(3) * (-t169 * t189 + (-t166 * t196 - t195) * t165)) * MDP(17) + (t180 * t169 + (-g(1) * t186 + t177) * t174) * pkin(3) * MDP(19) + (-g(1) * (t127 * t194 + t167 * t179) - g(2) * (-t125 * t167 + t185 * t194) - g(3) * (t130 * t194 + t167 * t178)) * MDP(32) + (-g(1) * (-t127 * t197 + t172 * t179) - g(2) * (-t125 * t172 - t185 * t197) - g(3) * (-t130 * t197 + t172 * t178)) * MDP(33) + (-MDP(25) * t173 + MDP(26) * t168) * (g(1) * t127 + g(2) * t185 + g(3) * t130); (-g(1) * t143 + g(2) * t141 - g(3) * t151) * MDP(19); (g(1) * t116 - g(2) * t114 + g(3) * t118) * MDP(26) + (-MDP(32) * t172 + MDP(33) * t167 - MDP(25)) * (g(1) * t115 + g(2) * t216 + g(3) * (t151 * t173 - t168 * t178)); (-g(1) * t110 - g(2) * t218 - g(3) * (-t118 * t167 - t130 * t172)) * MDP(32) + (g(1) * t111 - g(2) * t217 - g(3) * (-t118 * t172 + t130 * t167)) * MDP(33);];
taug  = t1;
