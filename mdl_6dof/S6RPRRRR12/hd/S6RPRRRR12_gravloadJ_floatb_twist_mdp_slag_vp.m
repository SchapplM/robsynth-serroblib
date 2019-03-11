% Calculate Gravitation load on the joints for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:45
% EndTime: 2019-03-09 08:00:50
% DurationCPUTime: 1.61s
% Computational Cost: add. (1263->146), mult. (3604->262), div. (0->0), fcn. (4743->18), ass. (0->79)
t163 = cos(pkin(14));
t160 = sin(pkin(14));
t171 = sin(qJ(1));
t195 = t171 * t160;
t166 = cos(pkin(6));
t175 = cos(qJ(1));
t197 = t166 * t175;
t156 = -t163 * t197 + t195;
t194 = t171 * t163;
t157 = t160 * t197 + t194;
t165 = cos(pkin(7));
t170 = sin(qJ(3));
t174 = cos(qJ(3));
t162 = sin(pkin(6));
t206 = sin(pkin(7));
t191 = t162 * t206;
t189 = t175 * t191;
t145 = (t156 * t165 + t189) * t174 + t157 * t170;
t198 = t165 * t170;
t146 = t156 * t198 - t157 * t174 + t170 * t189;
t164 = cos(pkin(8));
t169 = sin(qJ(4));
t207 = cos(qJ(4));
t161 = sin(pkin(8));
t200 = t162 * t175;
t183 = t156 * t206 - t165 * t200;
t209 = t183 * t161;
t125 = t146 * t207 + (t145 * t164 - t209) * t169;
t135 = t145 * t161 + t183 * t164;
t168 = sin(qJ(5));
t173 = cos(qJ(5));
t115 = t125 * t173 - t135 * t168;
t192 = t164 * t207;
t122 = t145 * t192 - t146 * t169 - t207 * t209;
t167 = sin(qJ(6));
t172 = cos(qJ(6));
t219 = t115 * t167 + t122 * t172;
t218 = t115 * t172 - t122 * t167;
t217 = t125 * t168 + t135 * t173;
t203 = t161 * t168;
t202 = t161 * t173;
t201 = t162 * t171;
t199 = t164 * t169;
t196 = t167 * t173;
t193 = t172 * t173;
t190 = t206 * t166;
t187 = g(1) * t171 - g(2) * t175;
t186 = t160 * t175 + t166 * t194;
t182 = -t163 * t191 + t166 * t165;
t179 = t182 * t161;
t178 = -t186 * t165 + t171 * t191;
t177 = t165 * t201 + t186 * t206;
t176 = t177 * t161;
t158 = t163 * t175 - t166 * t195;
t153 = t170 * t190 + (t160 * t174 + t163 * t198) * t162;
t152 = t174 * t190 + (t163 * t165 * t174 - t160 * t170) * t162;
t148 = t158 * t174 + t178 * t170;
t147 = -t158 * t170 + t178 * t174;
t142 = -t152 * t161 + t182 * t164;
t139 = t152 * t207 - t153 * t199;
t138 = t152 * t169 + t153 * t192;
t137 = -t147 * t161 + t177 * t164;
t134 = t153 * t207 + (t152 * t164 + t179) * t169;
t133 = -t152 * t192 + t153 * t169 - t207 * t179;
t132 = t139 * t173 + t153 * t203;
t131 = t147 * t207 - t148 * t199;
t130 = t147 * t169 + t148 * t192;
t129 = -t145 * t207 + t146 * t199;
t128 = -t145 * t169 - t146 * t192;
t127 = t148 * t207 + (t147 * t164 + t176) * t169;
t126 = -t147 * t192 + t148 * t169 - t207 * t176;
t121 = t134 * t173 + t142 * t168;
t119 = t131 * t173 + t148 * t203;
t118 = t129 * t173 - t146 * t203;
t117 = t127 * t173 + t137 * t168;
t116 = -t127 * t168 + t137 * t173;
t112 = t117 * t172 + t126 * t167;
t111 = -t117 * t167 + t126 * t172;
t1 = [t187 * MDP(2) + (g(1) * t157 - g(2) * t158) * MDP(4) + (-g(1) * t156 + g(2) * t186) * MDP(5) + (-g(1) * (-t171 * pkin(1) + qJ(2) * t200) - g(2) * (pkin(1) * t175 + qJ(2) * t201)) * MDP(7) + (-g(1) * t146 - g(2) * t148) * MDP(13) + (-g(1) * t145 - g(2) * t147) * MDP(14) + (-g(1) * t125 - g(2) * t127) * MDP(20) + (-g(1) * t122 + g(2) * t126) * MDP(21) + (-g(1) * t115 - g(2) * t117) * MDP(27) + (g(1) * t217 - g(2) * t116) * MDP(28) + (-g(1) * t218 - g(2) * t112) * MDP(34) + (g(1) * t219 - g(2) * t111) * MDP(35) + (-t162 * MDP(6) + MDP(3)) * (g(1) * t175 + g(2) * t171); (-g(3) * t166 - t187 * t162) * MDP(7); (-g(1) * t147 + g(2) * t145 - g(3) * t152) * MDP(13) + (g(1) * t148 - g(2) * t146 + g(3) * t153) * MDP(14) + (-g(1) * t131 - g(2) * t129 - g(3) * t139) * MDP(20) + (g(1) * t130 + g(2) * t128 + g(3) * t138) * MDP(21) + (-g(1) * t119 - g(2) * t118 - g(3) * t132) * MDP(27) + (-g(1) * (-t131 * t168 + t148 * t202) - g(2) * (-t129 * t168 - t146 * t202) - g(3) * (-t139 * t168 + t153 * t202)) * MDP(28) + (-g(1) * (t119 * t172 + t130 * t167) - g(2) * (t118 * t172 + t128 * t167) - g(3) * (t132 * t172 + t138 * t167)) * MDP(34) + (-g(1) * (-t119 * t167 + t130 * t172) - g(2) * (-t118 * t167 + t128 * t172) - g(3) * (-t132 * t167 + t138 * t172)) * MDP(35); (g(1) * t127 - g(2) * t125 + g(3) * t134) * MDP(21) + (-g(1) * (-t126 * t193 + t127 * t167) - g(2) * (-t122 * t193 - t125 * t167) - g(3) * (-t133 * t193 + t134 * t167)) * MDP(34) + (-g(1) * (t126 * t196 + t127 * t172) - g(2) * (t122 * t196 - t125 * t172) - g(3) * (t133 * t196 + t134 * t172)) * MDP(35) + (t173 * MDP(27) - MDP(28) * t168 + MDP(20)) * (g(1) * t126 + g(2) * t122 + g(3) * t133); (g(1) * t117 - g(2) * t115 + g(3) * t121) * MDP(28) + (-MDP(34) * t172 + MDP(35) * t167 - MDP(27)) * (g(1) * t116 + g(2) * t217 + g(3) * (-t134 * t168 + t142 * t173)); (-g(1) * t111 - g(2) * t219 - g(3) * (-t121 * t167 + t133 * t172)) * MDP(34) + (g(1) * t112 - g(2) * t218 - g(3) * (-t121 * t172 - t133 * t167)) * MDP(35);];
taug  = t1;
