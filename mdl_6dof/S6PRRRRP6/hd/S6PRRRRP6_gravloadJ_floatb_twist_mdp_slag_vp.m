% Calculate Gravitation load on the joints for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:27
% EndTime: 2019-03-09 00:34:30
% DurationCPUTime: 1.00s
% Computational Cost: add. (1023->152), mult. (2865->251), div. (0->0), fcn. (3686->14), ass. (0->82)
t250 = MDP(18) - MDP(27);
t230 = MDP(24) + MDP(26);
t249 = MDP(25) - MDP(28);
t197 = sin(pkin(7));
t203 = sin(qJ(2));
t206 = cos(qJ(2));
t199 = cos(pkin(12));
t240 = cos(pkin(6));
t225 = t199 * t240;
t238 = sin(pkin(12));
t209 = t203 * t238 - t206 * t225;
t198 = sin(pkin(6));
t235 = t198 * t199;
t239 = cos(pkin(7));
t252 = t197 * t235 + t209 * t239;
t221 = t240 * t238;
t210 = t199 * t203 + t206 * t221;
t226 = t198 * t238;
t251 = -t197 * t226 + t210 * t239;
t248 = cos(qJ(3));
t247 = pkin(9) * t197;
t201 = sin(qJ(4));
t237 = t197 * t201;
t205 = cos(qJ(4));
t236 = t197 * t205;
t234 = t198 * t203;
t233 = t198 * t206;
t200 = sin(qJ(5));
t232 = t200 * t205;
t204 = cos(qJ(5));
t231 = t204 * t205;
t228 = t197 * t234;
t227 = t197 * t240;
t202 = sin(qJ(3));
t224 = t202 * t239;
t222 = t239 * t248;
t189 = t203 * t225 + t206 * t238;
t162 = t189 * t248 - t252 * t202;
t179 = t197 * t209 - t235 * t239;
t146 = t162 * t205 + t179 * t201;
t161 = t189 * t202 + t252 * t248;
t130 = t146 * t200 - t161 * t204;
t190 = t199 * t206 - t203 * t221;
t164 = t190 * t248 - t251 * t202;
t180 = t197 * t210 + t226 * t239;
t148 = t164 * t205 + t180 * t201;
t163 = t190 * t202 + t251 * t248;
t132 = t148 * t200 - t163 * t204;
t178 = t202 * t227 + (t203 * t248 + t206 * t224) * t198;
t188 = -t197 * t233 + t239 * t240;
t166 = t178 * t205 + t188 * t201;
t177 = t202 * t234 - t222 * t233 - t227 * t248;
t143 = t166 * t200 - t177 * t204;
t218 = g(1) * t132 + g(2) * t130 + g(3) * t143;
t187 = (-t203 * t224 + t206 * t248) * t198;
t186 = (t202 * t206 + t203 * t222) * t198;
t174 = t187 * t205 + t201 * t228;
t173 = t187 * t201 - t205 * t228;
t172 = -t190 * t224 - t210 * t248;
t171 = t190 * t222 - t202 * t210;
t170 = -t189 * t224 - t209 * t248;
t169 = t189 * t222 - t202 * t209;
t156 = t174 * t204 + t186 * t200;
t155 = t174 * t200 - t186 * t204;
t154 = t172 * t205 + t190 * t237;
t153 = t172 * t201 - t190 * t236;
t152 = t170 * t205 + t189 * t237;
t151 = t170 * t201 - t189 * t236;
t150 = -t177 * t231 + t178 * t200;
t149 = -t177 * t232 - t178 * t204;
t144 = t166 * t204 + t177 * t200;
t142 = t154 * t204 + t171 * t200;
t141 = t154 * t200 - t171 * t204;
t140 = t152 * t204 + t169 * t200;
t139 = t152 * t200 - t169 * t204;
t138 = -t163 * t231 + t164 * t200;
t137 = -t163 * t232 - t164 * t204;
t136 = -t161 * t231 + t162 * t200;
t135 = -t161 * t232 - t162 * t204;
t133 = t148 * t204 + t163 * t200;
t131 = t146 * t204 + t161 * t200;
t1 = [(-MDP(1) - MDP(29)) * g(3); (g(1) * t210 + g(2) * t209 - g(3) * t233) * MDP(3) + (g(1) * t190 + g(2) * t189 + g(3) * t234) * MDP(4) + (-g(1) * t172 - g(2) * t170 - g(3) * t187) * MDP(10) + (g(1) * t171 + g(2) * t169 + g(3) * t186) * MDP(11) + (-g(1) * t154 - g(2) * t152 - g(3) * t174) * MDP(17) + (-g(1) * (-pkin(2) * t210 + t172 * pkin(3) + t154 * pkin(4) + t142 * pkin(5) + t171 * pkin(10) + t153 * pkin(11) + t141 * qJ(6) + t190 * t247) - g(2) * (-pkin(2) * t209 + t170 * pkin(3) + t152 * pkin(4) + t140 * pkin(5) + t169 * pkin(10) + t151 * pkin(11) + t139 * qJ(6) + t189 * t247) - g(3) * (t187 * pkin(3) + t174 * pkin(4) + t156 * pkin(5) + t186 * pkin(10) + t173 * pkin(11) + t155 * qJ(6) + (pkin(2) * t206 + t203 * t247) * t198)) * MDP(29) + t230 * (-g(1) * t142 - g(2) * t140 - g(3) * t156) + t249 * (g(1) * t141 + g(2) * t139 + g(3) * t155) + t250 * (g(1) * t153 + g(2) * t151 + g(3) * t173); (g(1) * t164 + g(2) * t162 + g(3) * t178) * MDP(11) + (-g(1) * (t138 * pkin(5) + t164 * pkin(10) + t137 * qJ(6)) - g(2) * (t136 * pkin(5) + t162 * pkin(10) + t135 * qJ(6)) - g(3) * (t150 * pkin(5) + t178 * pkin(10) + t149 * qJ(6))) * MDP(29) + t249 * (g(1) * t137 + g(2) * t135 + g(3) * t149) + t230 * (-g(1) * t138 - g(2) * t136 - g(3) * t150) + (MDP(10) + t205 * MDP(17) + (pkin(4) * t205 + pkin(11) * t201 + pkin(3)) * MDP(29) - t250 * t201) * (g(1) * t163 + g(2) * t161 + g(3) * t177); (-MDP(29) * pkin(11) + t250) * (g(1) * t148 + g(2) * t146 + g(3) * t166) + (-MDP(17) - MDP(29) * (pkin(5) * t204 + qJ(6) * t200 + pkin(4)) - t230 * t204 + t249 * t200) * (g(3) * (-t178 * t201 + t188 * t205) + g(2) * (-t162 * t201 + t179 * t205) + g(1) * (-t164 * t201 + t180 * t205)); (-g(1) * (-pkin(5) * t132 + qJ(6) * t133) - g(2) * (-pkin(5) * t130 + qJ(6) * t131) - g(3) * (-pkin(5) * t143 + qJ(6) * t144)) * MDP(29) + t230 * t218 + t249 * (g(1) * t133 + g(2) * t131 + g(3) * t144); -t218 * MDP(29);];
taug  = t1;
