% Calculate Gravitation load on the joints for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:16:55
% EndTime: 2019-03-10 06:17:02
% DurationCPUTime: 2.69s
% Computational Cost: add. (1544->205), mult. (4444->366), div. (0->0), fcn. (5838->18), ass. (0->108)
t214 = sin(qJ(2));
t218 = cos(qJ(2));
t219 = cos(qJ(1));
t257 = cos(pkin(6));
t231 = t219 * t257;
t258 = sin(qJ(1));
t198 = t214 * t258 - t218 * t231;
t199 = t214 * t231 + t218 * t258;
t209 = cos(pkin(7));
t213 = sin(qJ(3));
t217 = cos(qJ(3));
t206 = sin(pkin(7));
t207 = sin(pkin(6));
t247 = t207 * t219;
t236 = t206 * t247;
t180 = (t198 * t209 + t236) * t217 + t199 * t213;
t245 = t209 * t213;
t181 = t198 * t245 - t199 * t217 + t213 * t236;
t208 = cos(pkin(8));
t212 = sin(qJ(4));
t194 = -t198 * t206 + t209 * t247;
t205 = sin(pkin(8));
t256 = t194 * t205;
t259 = cos(qJ(4));
t151 = t181 * t259 + (t180 * t208 + t256) * t212;
t166 = t180 * t205 - t194 * t208;
t211 = sin(qJ(5));
t216 = cos(qJ(5));
t139 = t151 * t216 - t166 * t211;
t210 = sin(qJ(6));
t268 = t139 * t210;
t215 = cos(qJ(6));
t267 = t139 * t215;
t266 = t151 * t211 + t166 * t216;
t233 = t208 * t259;
t262 = -t180 * t233 + t181 * t212;
t253 = t205 * t206;
t252 = t205 * t211;
t251 = t205 * t216;
t250 = t206 * t208;
t249 = t207 * t214;
t248 = t207 * t218;
t246 = t208 * t212;
t244 = t209 * t217;
t243 = t210 * t216;
t242 = t213 * t214;
t241 = t213 * t218;
t240 = t214 * t217;
t239 = t215 * t216;
t238 = t217 * t218;
t237 = t206 * t249;
t235 = t205 * t259;
t234 = t207 * t258;
t232 = t206 * t257;
t230 = t206 * t235;
t229 = t257 * t258;
t200 = -t214 * t219 - t218 * t229;
t225 = -t200 * t206 + t209 * t234;
t224 = t200 * t209 + t206 * t234;
t223 = -t206 * t248 + t209 * t257;
t221 = t225 * t205;
t220 = t223 * t205;
t201 = -t214 * t229 + t218 * t219;
t197 = (-t209 * t242 + t238) * t207;
t196 = (-t209 * t240 - t241) * t207;
t193 = t213 * t232 + (t209 * t241 + t240) * t207;
t192 = t217 * t232 + (t209 * t238 - t242) * t207;
t188 = -t196 * t205 + t208 * t237;
t187 = t200 * t217 - t201 * t245;
t186 = -t200 * t213 - t201 * t244;
t185 = -t198 * t217 - t199 * t245;
t184 = t198 * t213 - t199 * t244;
t183 = t201 * t217 + t213 * t224;
t182 = -t201 * t213 + t217 * t224;
t177 = -t192 * t205 + t208 * t223;
t174 = -t186 * t205 + t201 * t250;
t173 = -t184 * t205 + t199 * t250;
t172 = t197 * t259 + (t196 * t208 + t205 * t237) * t212;
t171 = -t196 * t233 + t197 * t212 - t230 * t249;
t170 = t192 * t259 - t193 * t246;
t169 = t192 * t212 + t193 * t233;
t168 = -t182 * t205 + t208 * t225;
t165 = t193 * t259 + (t192 * t208 + t220) * t212;
t164 = -t192 * t233 + t193 * t212 - t220 * t259;
t163 = t170 * t216 + t193 * t252;
t162 = t172 * t216 + t188 * t211;
t161 = t182 * t259 - t183 * t246;
t160 = t182 * t212 + t183 * t233;
t159 = -t180 * t259 + t181 * t246;
t158 = -t180 * t212 - t181 * t233;
t157 = t187 * t259 + (t186 * t208 + t201 * t253) * t212;
t156 = -t186 * t233 + t187 * t212 - t201 * t230;
t155 = t185 * t259 + (t184 * t208 + t199 * t253) * t212;
t154 = -t184 * t233 + t185 * t212 - t199 * t230;
t153 = t183 * t259 + (t182 * t208 + t221) * t212;
t152 = -t182 * t233 + t183 * t212 - t221 * t259;
t150 = -t194 * t235 + t262;
t148 = t256 * t259 - t262;
t147 = t165 * t216 + t177 * t211;
t145 = t161 * t216 + t183 * t252;
t144 = t159 * t216 - t181 * t252;
t143 = t157 * t216 + t174 * t211;
t142 = t155 * t216 + t173 * t211;
t141 = t153 * t216 + t168 * t211;
t140 = -t153 * t211 + t168 * t216;
t136 = t141 * t215 + t152 * t210;
t135 = -t141 * t210 + t152 * t215;
t1 = [(g(1) * t258 - g(2) * t219) * MDP(2) + (g(1) * t219 + g(2) * t258) * MDP(3) + (g(1) * t199 - g(2) * t201) * MDP(9) + (-g(1) * t198 - g(2) * t200) * MDP(10) + (-g(1) * t181 - g(2) * t183) * MDP(16) + (-g(1) * t180 - g(2) * t182) * MDP(17) + (-g(1) * t151 - g(2) * t153) * MDP(23) + (g(1) * t150 + g(2) * t152) * MDP(24) + (-g(1) * t139 - g(2) * t141) * MDP(30) + (g(1) * t266 - g(2) * t140) * MDP(31) + (-g(1) * (t150 * t210 + t267) - g(2) * t136) * MDP(37) + (-g(1) * (t150 * t215 - t268) - g(2) * t135) * MDP(38); (-g(1) * t200 + g(2) * t198 - g(3) * t248) * MDP(9) + (g(1) * t201 + g(2) * t199 + g(3) * t249) * MDP(10) + (-g(1) * t187 - g(2) * t185 - g(3) * t197) * MDP(16) + (-g(1) * t186 - g(2) * t184 - g(3) * t196) * MDP(17) + (-g(1) * t157 - g(2) * t155 - g(3) * t172) * MDP(23) + (g(1) * t156 + g(2) * t154 + g(3) * t171) * MDP(24) + (-g(1) * t143 - g(2) * t142 - g(3) * t162) * MDP(30) + (-g(1) * (-t157 * t211 + t174 * t216) - g(2) * (-t155 * t211 + t173 * t216) - g(3) * (-t172 * t211 + t188 * t216)) * MDP(31) + (-g(1) * (t143 * t215 + t156 * t210) - g(2) * (t142 * t215 + t154 * t210) - g(3) * (t162 * t215 + t171 * t210)) * MDP(37) + (-g(1) * (-t143 * t210 + t156 * t215) - g(2) * (-t142 * t210 + t154 * t215) - g(3) * (-t162 * t210 + t171 * t215)) * MDP(38); (-g(1) * t182 + g(2) * t180 - g(3) * t192) * MDP(16) + (g(1) * t183 - g(2) * t181 + g(3) * t193) * MDP(17) + (-g(1) * t161 - g(2) * t159 - g(3) * t170) * MDP(23) + (g(1) * t160 + g(2) * t158 + g(3) * t169) * MDP(24) + (-g(1) * t145 - g(2) * t144 - g(3) * t163) * MDP(30) + (-g(1) * (-t161 * t211 + t183 * t251) - g(2) * (-t159 * t211 - t181 * t251) - g(3) * (-t170 * t211 + t193 * t251)) * MDP(31) + (-g(1) * (t145 * t215 + t160 * t210) - g(2) * (t144 * t215 + t158 * t210) - g(3) * (t163 * t215 + t169 * t210)) * MDP(37) + (-g(1) * (-t145 * t210 + t160 * t215) - g(2) * (-t144 * t210 + t158 * t215) - g(3) * (-t163 * t210 + t169 * t215)) * MDP(38); (g(1) * t153 - g(2) * t151 + g(3) * t165) * MDP(24) + (-g(1) * (-t152 * t239 + t153 * t210) - g(2) * (-t148 * t239 - t151 * t210) - g(3) * (-t164 * t239 + t165 * t210)) * MDP(37) + (-g(1) * (t152 * t243 + t153 * t215) - g(2) * (t148 * t243 - t151 * t215) - g(3) * (t164 * t243 + t165 * t215)) * MDP(38) + (MDP(30) * t216 - MDP(31) * t211 + MDP(23)) * (g(1) * t152 + g(2) * t148 + g(3) * t164); (g(1) * t141 - g(2) * t139 + g(3) * t147) * MDP(31) + (-MDP(37) * t215 + MDP(38) * t210 - MDP(30)) * (g(1) * t140 + g(2) * t266 + g(3) * (-t165 * t211 + t177 * t216)); (-g(1) * t135 - g(2) * (t148 * t215 + t268) - g(3) * (-t147 * t210 + t164 * t215)) * MDP(37) + (g(1) * t136 - g(2) * (-t148 * t210 + t267) - g(3) * (-t147 * t215 - t164 * t210)) * MDP(38);];
taug  = t1;
