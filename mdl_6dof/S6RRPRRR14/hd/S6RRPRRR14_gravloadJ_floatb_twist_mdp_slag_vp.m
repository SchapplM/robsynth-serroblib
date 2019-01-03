% Calculate Gravitation load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x6] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:20:49
% EndTime: 2019-01-03 10:20:58
% DurationCPUTime: 2.23s
% Computational Cost: add. (1301->179), mult. (3718->323), div. (0->0), fcn. (4867->18), ass. (0->96)
t197 = sin(qJ(2));
t198 = sin(qJ(1));
t201 = cos(qJ(2));
t202 = cos(qJ(1));
t244 = cos(pkin(6));
t222 = t202 * t244;
t182 = t197 * t222 + t198 * t201;
t187 = sin(pkin(14));
t191 = cos(pkin(14));
t181 = t197 * t198 - t201 * t222;
t189 = sin(pkin(7));
t193 = cos(pkin(7));
t190 = sin(pkin(6));
t232 = t190 * t202;
t219 = t181 * t193 + t189 * t232;
t166 = -t182 * t191 + t219 * t187;
t196 = sin(qJ(4));
t246 = cos(qJ(4));
t177 = -t181 * t189 + t193 * t232;
t188 = sin(pkin(8));
t192 = cos(pkin(8));
t247 = t182 * t187 + t219 * t191;
t249 = t177 * t188 + t192 * t247;
t144 = t166 * t246 + t249 * t196;
t154 = -t177 * t192 + t188 * t247;
t195 = sin(qJ(5));
t200 = cos(qJ(5));
t134 = t144 * t200 - t154 * t195;
t194 = sin(qJ(6));
t258 = t134 * t194;
t199 = cos(qJ(6));
t257 = t134 * t199;
t256 = t144 * t195 + t154 * t200;
t253 = t166 * t196;
t221 = t244 * t189;
t229 = t193 * t201;
t204 = t191 * t221 + (-t187 * t197 + t191 * t229) * t190;
t233 = t190 * t201;
t213 = t189 * t233 - t244 * t193;
t251 = -t213 * t188 + t204 * t192;
t223 = t198 * t244;
t184 = -t197 * t223 + t201 * t202;
t183 = -t202 * t197 - t201 * t223;
t234 = t190 * t198;
t217 = t183 * t193 + t189 * t234;
t207 = t184 * t187 - t217 * t191;
t218 = -t183 * t189 + t193 * t234;
t250 = -t218 * t188 + t207 * t192;
t248 = -g(1) * t177 - g(2) * t218;
t243 = qJ(3) * t189;
t238 = t187 * t193;
t237 = t188 * t189;
t236 = t189 * t192;
t235 = t190 * t197;
t231 = t191 * t193;
t230 = t193 * t197;
t228 = t194 * t200;
t227 = t199 * t200;
t226 = t189 * t235;
t225 = t188 * t246;
t224 = t192 * t246;
t220 = t189 * t225;
t180 = (-t187 * t230 + t191 * t201) * t190;
t179 = (-t187 * t201 - t191 * t230) * t190;
t176 = t191 * t235 + (t190 * t229 + t221) * t187;
t172 = -t179 * t188 + t192 * t226;
t171 = t183 * t191 - t184 * t238;
t170 = -t183 * t187 - t184 * t231;
t169 = -t181 * t191 - t182 * t238;
t168 = t181 * t187 - t182 * t231;
t167 = t184 * t191 + t217 * t187;
t163 = -t204 * t188 - t213 * t192;
t160 = -t170 * t188 + t184 * t236;
t159 = -t168 * t188 + t182 * t236;
t158 = t180 * t246 + (t179 * t192 + t188 * t226) * t196;
t157 = -t179 * t224 + t180 * t196 - t220 * t235;
t156 = t207 * t188 + t218 * t192;
t153 = t176 * t246 + t251 * t196;
t152 = t176 * t196 - t251 * t246;
t151 = t158 * t200 + t172 * t195;
t150 = t171 * t246 + (t170 * t192 + t184 * t237) * t196;
t149 = -t170 * t224 + t171 * t196 - t184 * t220;
t148 = t169 * t246 + (t168 * t192 + t182 * t237) * t196;
t147 = -t168 * t224 + t169 * t196 - t182 * t220;
t146 = t167 * t246 - t250 * t196;
t145 = t167 * t196 + t250 * t246;
t143 = -t177 * t225 - t224 * t247 + t253;
t141 = t249 * t246 - t253;
t140 = t153 * t200 + t163 * t195;
t138 = t150 * t200 + t160 * t195;
t137 = t148 * t200 + t159 * t195;
t136 = t146 * t200 + t156 * t195;
t135 = -t146 * t195 + t156 * t200;
t131 = t136 * t199 + t145 * t194;
t130 = -t136 * t194 + t145 * t199;
t1 = [(g(1) * t198 - g(2) * t202) * MDP(2) + (g(1) * t202 + g(2) * t198) * MDP(3) + (g(1) * t182 - g(2) * t184) * MDP(9) + (-g(1) * t181 - g(2) * t183) * MDP(10) + (-g(1) * t166 - g(2) * t167) * MDP(11) + (-g(1) * t247 + g(2) * t207) * MDP(12) + t248 * MDP(13) + (-g(1) * (-t198 * pkin(1) - t182 * pkin(2) + pkin(10) * t232) - g(2) * (pkin(1) * t202 + t184 * pkin(2) + pkin(10) * t234) + t248 * qJ(3)) * MDP(14) + (-g(1) * t144 - g(2) * t146) * MDP(20) + (g(1) * t143 + g(2) * t145) * MDP(21) + (-g(1) * t134 - g(2) * t136) * MDP(27) + (g(1) * t256 - g(2) * t135) * MDP(28) + (-g(1) * (t143 * t194 + t257) - g(2) * t131) * MDP(34) + (-g(1) * (t143 * t199 - t258) - g(2) * t130) * MDP(35); (-g(1) * t183 + g(2) * t181 - g(3) * t233) * MDP(9) + (-g(1) * t171 - g(2) * t169 - g(3) * t180) * MDP(11) + (-g(1) * t170 - g(2) * t168 - g(3) * t179) * MDP(12) + (-g(1) * (pkin(2) * t183 + t184 * t243) - g(2) * (-pkin(2) * t181 + t182 * t243) - g(3) * (pkin(2) * t201 + t197 * t243) * t190) * MDP(14) + (-g(1) * t150 - g(2) * t148 - g(3) * t158) * MDP(20) + (g(1) * t149 + g(2) * t147 + g(3) * t157) * MDP(21) + (-g(1) * t138 - g(2) * t137 - g(3) * t151) * MDP(27) + (-g(1) * (-t150 * t195 + t160 * t200) - g(2) * (-t148 * t195 + t159 * t200) - g(3) * (-t158 * t195 + t172 * t200)) * MDP(28) + (-g(1) * (t138 * t199 + t149 * t194) - g(2) * (t137 * t199 + t147 * t194) - g(3) * (t151 * t199 + t157 * t194)) * MDP(34) + (-g(1) * (-t138 * t194 + t149 * t199) - g(2) * (-t137 * t194 + t147 * t199) - g(3) * (-t151 * t194 + t157 * t199)) * MDP(35) + (-t189 * MDP(13) + MDP(10)) * (g(1) * t184 + g(2) * t182 + g(3) * t235); (-g(1) * t218 + g(2) * t177 + g(3) * t213) * MDP(14); (g(1) * t146 - g(2) * t144 + g(3) * t153) * MDP(21) + (-g(1) * (-t145 * t227 + t146 * t194) - g(2) * (-t141 * t227 - t144 * t194) - g(3) * (-t152 * t227 + t153 * t194)) * MDP(34) + (-g(1) * (t145 * t228 + t146 * t199) - g(2) * (t141 * t228 - t144 * t199) - g(3) * (t152 * t228 + t153 * t199)) * MDP(35) + (t200 * MDP(27) - MDP(28) * t195 + MDP(20)) * (g(1) * t145 + g(2) * t141 + g(3) * t152); (g(1) * t136 - g(2) * t134 + g(3) * t140) * MDP(28) + (-MDP(34) * t199 + MDP(35) * t194 - MDP(27)) * (g(1) * t135 + g(2) * t256 + g(3) * (-t153 * t195 + t163 * t200)); (-g(1) * t130 - g(2) * (t141 * t199 + t258) - g(3) * (-t140 * t194 + t152 * t199)) * MDP(34) + (g(1) * t131 - g(2) * (-t141 * t194 + t257) - g(3) * (-t140 * t199 - t152 * t194)) * MDP(35);];
taug  = t1;
