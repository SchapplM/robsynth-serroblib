% Calculate Gravitation load on the joints for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:44
% EndTime: 2019-03-09 21:26:47
% DurationCPUTime: 1.16s
% Computational Cost: add. (645->163), mult. (1360->254), div. (0->0), fcn. (1627->12), ass. (0->79)
t205 = sin(qJ(2));
t208 = cos(qJ(2));
t248 = cos(pkin(6));
t251 = cos(qJ(1));
t223 = t248 * t251;
t250 = sin(qJ(1));
t181 = t205 * t250 - t208 * t223;
t206 = cos(qJ(4));
t242 = t181 * t206;
t182 = t205 * t223 + t208 * t250;
t204 = sin(qJ(3));
t207 = cos(qJ(3));
t201 = sin(pkin(6));
t225 = t201 * t251;
t156 = t182 * t207 - t204 * t225;
t203 = sin(qJ(4));
t247 = t156 * t203;
t257 = t242 - t247;
t202 = -qJ(5) - pkin(10);
t256 = -t202 * t204 + pkin(2);
t196 = pkin(4) * t206 + pkin(3);
t255 = -t196 * t207 - t256;
t200 = qJ(4) + pkin(11);
t197 = sin(t200);
t198 = cos(t200);
t132 = t156 * t197 - t181 * t198;
t133 = t156 * t198 + t181 * t197;
t243 = t181 * t203;
t254 = t156 * t206 + t243;
t253 = MDP(17) - MDP(25) - MDP(28);
t252 = MDP(26) + MDP(30);
t249 = g(3) * t201;
t222 = t248 * t250;
t184 = -t205 * t222 + t208 * t251;
t224 = t201 * t250;
t160 = t184 * t207 + t204 * t224;
t246 = t160 * t203;
t241 = t182 * t203;
t183 = t205 * t251 + t208 * t222;
t240 = t183 * t203;
t239 = t183 * t206;
t238 = t184 * t203;
t236 = t197 * t207;
t235 = t198 * t207;
t234 = t201 * t205;
t233 = t201 * t208;
t231 = t203 * t205;
t230 = t203 * t207;
t229 = t206 * t207;
t228 = t207 * t208;
t226 = t206 * t233;
t180 = t204 * t248 + t207 * t234;
t220 = -t180 * t203 - t226;
t136 = t160 * t197 - t183 * t198;
t149 = t180 * t197 + t198 * t233;
t219 = g(1) * t136 + g(2) * t132 + g(3) * t149;
t155 = t182 * t204 + t207 * t225;
t159 = t184 * t204 - t207 * t224;
t179 = t204 * t234 - t207 * t248;
t218 = g(1) * t159 + g(2) * t155 + g(3) * t179;
t217 = g(1) * t160 + g(2) * t156 + g(3) * t180;
t215 = pkin(1) * t251 + t184 * pkin(2) + pkin(4) * t240 + pkin(8) * t224 + pkin(9) * t183 - t159 * t202 + t160 * t196;
t213 = pkin(9) * t234 + t256 * t233 + (pkin(4) * t231 + t196 * t228) * t201;
t212 = pkin(4) * t241 + pkin(9) * t182 + t255 * t181;
t211 = pkin(4) * t238 + pkin(9) * t184 + t255 * t183;
t210 = -pkin(1) * t250 - t182 * pkin(2) - pkin(4) * t243 + pkin(8) * t225 - t181 * pkin(9) + t155 * t202 - t156 * t196;
t172 = pkin(4) * t239;
t168 = pkin(4) * t242;
t162 = (t197 * t205 + t198 * t228) * t201;
t161 = (t197 * t228 - t198 * t205) * t201;
t150 = t180 * t198 - t197 * t233;
t144 = -t183 * t235 + t184 * t197;
t143 = -t183 * t236 - t184 * t198;
t142 = -t181 * t235 + t182 * t197;
t141 = -t181 * t236 - t182 * t198;
t139 = t160 * t206 + t240;
t138 = t239 - t246;
t137 = t160 * t198 + t183 * t197;
t1 = [(g(1) * t250 - g(2) * t251) * MDP(2) + (g(1) * t251 + g(2) * t250) * MDP(3) + (g(1) * t182 - g(2) * t184) * MDP(9) + (-g(1) * t181 + g(2) * t183) * MDP(10) + (g(1) * t156 - g(2) * t160) * MDP(16) + (g(1) * t254 - g(2) * t139) * MDP(23) + (g(1) * t257 - g(2) * t138) * MDP(24) + (-g(1) * t210 - g(2) * t215) * MDP(26) + (g(1) * t133 - g(2) * t137) * MDP(27) + (g(1) * t132 - g(2) * t136) * MDP(29) + (-g(1) * (-pkin(5) * t133 - qJ(6) * t132 + t210) - g(2) * (pkin(5) * t137 + qJ(6) * t136 + t215)) * MDP(30) + t253 * (-g(1) * t155 + g(2) * t159); (g(1) * t184 + g(2) * t182 + g(3) * t234) * MDP(10) + (-g(1) * (-t183 * t229 + t238) - g(2) * (-t181 * t229 + t241) - (t206 * t228 + t231) * t249) * MDP(23) + (-g(1) * (t183 * t230 + t184 * t206) - g(2) * (t181 * t230 + t182 * t206) - (-t203 * t228 + t205 * t206) * t249) * MDP(24) + (-g(1) * t211 - g(2) * t212 - g(3) * t213) * MDP(26) + (-g(1) * t144 - g(2) * t142 - g(3) * t162) * MDP(27) + (-g(1) * t143 - g(2) * t141 - g(3) * t161) * MDP(29) + (-g(1) * (pkin(5) * t144 + qJ(6) * t143 + t211) - g(2) * (pkin(5) * t142 + qJ(6) * t141 + t212) - g(3) * (t162 * pkin(5) + t161 * qJ(6) + t213)) * MDP(30) + (-MDP(16) * t207 + t253 * t204 - MDP(9)) * (-g(1) * t183 - g(2) * t181 + g(3) * t233); t253 * t217 + t252 * (-g(1) * (-t159 * t196 - t160 * t202) - g(2) * (-t155 * t196 - t156 * t202) - g(3) * (-t179 * t196 - t180 * t202)) + ((pkin(5) * t198 + qJ(6) * t197) * MDP(30) + MDP(23) * t206 - MDP(24) * t203 + MDP(27) * t198 + MDP(29) * t197 + MDP(16)) * t218; (-g(1) * t138 - g(2) * t257 - g(3) * t220) * MDP(23) + (g(1) * t139 + g(2) * t254 - g(3) * (-t180 * t206 + t203 * t233)) * MDP(24) + (-g(1) * t172 - g(2) * t168 + (g(3) * t226 + t203 * t217) * pkin(4)) * MDP(26) + t219 * MDP(27) + (-g(1) * t137 - g(2) * t133 - g(3) * t150) * MDP(29) + (-g(1) * (-pkin(4) * t246 - pkin(5) * t136 + qJ(6) * t137 + t172) - g(2) * (-pkin(4) * t247 - pkin(5) * t132 + qJ(6) * t133 + t168) - g(3) * (pkin(4) * t220 - t149 * pkin(5) + t150 * qJ(6))) * MDP(30); -t252 * t218; -t219 * MDP(30);];
taug  = t1;
