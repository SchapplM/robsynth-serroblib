% Calculate Gravitation load on the joints for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% MDP [45x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [7x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S7RRRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1),zeros(45,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [45 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [45x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:42:59
% EndTime: 2019-03-10 07:43:08
% DurationCPUTime: 2.76s
% Computational Cost: add. (733->196), mult. (1998->342), div. (0->0), fcn. (2528->14), ass. (0->91)
t188 = sin(qJ(1));
t193 = cos(qJ(3));
t194 = cos(qJ(2));
t211 = t193 * t194;
t186 = sin(qJ(3));
t195 = cos(qJ(1));
t219 = t186 * t195;
t172 = t188 * t211 + t219;
t192 = cos(qJ(4));
t185 = sin(qJ(4));
t187 = sin(qJ(2));
t222 = t185 * t187;
t154 = t172 * t192 + t188 * t222;
t210 = t195 * t193;
t215 = t188 * t186;
t171 = t194 * t215 - t210;
t184 = sin(qJ(5));
t191 = cos(qJ(5));
t135 = t154 * t191 - t171 * t184;
t218 = t187 * t192;
t153 = t172 * t185 - t188 * t218;
t183 = sin(qJ(6));
t190 = cos(qJ(6));
t122 = t135 * t190 + t153 * t183;
t134 = t154 * t184 + t171 * t191;
t182 = sin(qJ(7));
t189 = cos(qJ(7));
t240 = t122 * t182 + t134 * t189;
t239 = t122 * t189 - t134 * t182;
t236 = t135 * t183 - t153 * t190;
t233 = -MDP(38) * t190 + MDP(24);
t232 = -MDP(38) * t183 + MDP(30);
t227 = t182 * t184;
t226 = t182 * t190;
t225 = t183 * t185;
t224 = t184 * t189;
t223 = t184 * t192;
t221 = t186 * t187;
t220 = t186 * t194;
t217 = t187 * t193;
t216 = t187 * t195;
t214 = t189 * t190;
t213 = t190 * t191;
t212 = t191 * t192;
t207 = t184 * t221;
t206 = t191 * t221;
t205 = t186 * t216;
t204 = g(1) * t195 + g(2) * t188;
t170 = -t185 * t194 + t192 * t217;
t169 = t185 * t217 + t192 * t194;
t176 = t194 * t210 - t215;
t175 = -t188 * t193 - t194 * t219;
t174 = t192 * t211 + t222;
t173 = t185 * t211 - t218;
t166 = t170 * t195;
t165 = t169 * t195;
t164 = t170 * t188;
t163 = t169 * t188;
t162 = (-t184 * t193 - t186 * t212) * t187;
t161 = -t191 * t217 + t192 * t207;
t160 = t176 * t192 + t185 * t216;
t159 = t176 * t185 - t192 * t216;
t158 = t174 * t191 - t184 * t220;
t157 = t174 * t184 + t191 * t220;
t152 = t170 * t191 - t207;
t151 = t170 * t184 + t206;
t150 = -t166 * t191 + t184 * t205;
t149 = -t166 * t184 - t191 * t205;
t148 = -t164 * t191 + t188 * t207;
t147 = -t164 * t184 - t188 * t206;
t146 = t162 * t190 - t221 * t225;
t145 = t175 * t212 - t176 * t184;
t144 = t175 * t223 + t176 * t191;
t143 = -t171 * t212 - t172 * t184;
t142 = -t171 * t223 + t172 * t191;
t141 = -t169 * t213 + t170 * t183;
t140 = t160 * t191 + t175 * t184;
t139 = t160 * t184 - t175 * t191;
t138 = t158 * t190 + t173 * t183;
t133 = t152 * t190 + t169 * t183;
t131 = t150 * t190 - t165 * t183;
t130 = t148 * t190 - t163 * t183;
t129 = t145 * t190 + t175 * t225;
t128 = t143 * t190 - t171 * t225;
t127 = -t159 * t213 + t160 * t183;
t126 = -t153 * t213 + t154 * t183;
t125 = t140 * t190 + t159 * t183;
t124 = -t140 * t183 + t159 * t190;
t120 = t125 * t189 - t139 * t182;
t119 = -t125 * t182 - t139 * t189;
t1 = [t204 * MDP(3) + (g(1) * t172 - g(2) * t176) * MDP(16) + (-g(1) * t171 - g(2) * t175) * MDP(17) + (g(1) * t154 - g(2) * t160) * MDP(23) + (-g(1) * t153 + g(2) * t159) * MDP(24) + (g(1) * t135 - g(2) * t140) * MDP(30) + (-g(1) * t134 + g(2) * t139) * MDP(31) + (g(1) * t122 - g(2) * t125) * MDP(37) + (-g(1) * t236 - g(2) * t124) * MDP(38) + (g(1) * t239 - g(2) * t120) * MDP(44) + (-g(1) * t240 - g(2) * t119) * MDP(45) + (-t187 * MDP(10) + t194 * MDP(9) + MDP(2)) * (g(1) * t188 - g(2) * t195); (g(3) * t187 + t204 * t194) * MDP(10) + (g(1) * t166 + g(2) * t164 - g(3) * t174) * MDP(23) + (-g(1) * t165 - g(2) * t163 + g(3) * t173) * MDP(24) + (-g(1) * t150 - g(2) * t148 - g(3) * t158) * MDP(30) + (g(1) * t149 + g(2) * t147 + g(3) * t157) * MDP(31) + (-g(1) * t131 - g(2) * t130 - g(3) * t138) * MDP(37) + (-g(1) * (-t150 * t183 - t165 * t190) - g(2) * (-t148 * t183 - t163 * t190) - g(3) * (-t158 * t183 + t173 * t190)) * MDP(38) + (-g(1) * (t131 * t189 - t149 * t182) - g(2) * (t130 * t189 - t147 * t182) - g(3) * (t138 * t189 - t157 * t182)) * MDP(44) + (-g(1) * (-t131 * t182 - t149 * t189) - g(2) * (-t130 * t182 - t147 * t189) - g(3) * (-t138 * t182 - t157 * t189)) * MDP(45) + (MDP(16) * t193 - MDP(17) * t186 + MDP(9)) * (-g(3) * t194 + t204 * t187); (g(1) * t176 + g(2) * t172 + g(3) * t217) * MDP(17) + (g(1) * t144 + g(2) * t142 - g(3) * t161) * MDP(31) + (-g(1) * t129 - g(2) * t128 - g(3) * t146) * MDP(37) + (-g(1) * (t129 * t189 - t144 * t182) - g(2) * (t128 * t189 - t142 * t182) - g(3) * (t146 * t189 + t161 * t182)) * MDP(44) + (-g(1) * (-t129 * t182 - t144 * t189) - g(2) * (-t128 * t182 - t142 * t189) - g(3) * (-t146 * t182 + t161 * t189)) * MDP(45) - t232 * (g(1) * t145 + g(2) * t143 + g(3) * t162) + (MDP(23) * t192 - t185 * t233 + MDP(16)) * (-g(1) * t175 + g(2) * t171 + g(3) * t221); (-g(1) * t127 - g(2) * t126 - g(3) * t141) * MDP(37) + (-g(1) * (t127 * t189 + t159 * t227) - g(2) * (t126 * t189 + t153 * t227) - g(3) * (t141 * t189 + t169 * t227)) * MDP(44) + (-g(1) * (-t127 * t182 + t159 * t224) - g(2) * (-t126 * t182 + t153 * t224) - g(3) * (-t141 * t182 + t169 * t224)) * MDP(45) + t233 * (g(1) * t160 + g(2) * t154 + g(3) * t170) + (-t184 * MDP(31) + t191 * t232 + MDP(23)) * (g(1) * t159 + g(2) * t153 + g(3) * t169); (g(1) * t140 + g(2) * t135 + g(3) * t152) * MDP(31) + (-g(1) * (-t139 * t214 - t140 * t182) - g(2) * (-t134 * t214 - t135 * t182) - g(3) * (-t151 * t214 - t152 * t182)) * MDP(44) + (-g(1) * (t139 * t226 - t140 * t189) - g(2) * (t134 * t226 - t135 * t189) - g(3) * (t151 * t226 - t152 * t189)) * MDP(45) + (t190 * MDP(37) + t232) * (g(1) * t139 + g(2) * t134 + g(3) * t151); (g(1) * t125 + g(2) * t122 + g(3) * t133) * MDP(38) + (-MDP(44) * t189 + MDP(45) * t182 - MDP(37)) * (g(1) * t124 - g(2) * t236 + g(3) * (-t152 * t183 + t169 * t190)); (-g(1) * t119 + g(2) * t240 - g(3) * (-t133 * t182 - t151 * t189)) * MDP(44) + (g(1) * t120 + g(2) * t239 - g(3) * (-t133 * t189 + t151 * t182)) * MDP(45);];
taug  = t1;
