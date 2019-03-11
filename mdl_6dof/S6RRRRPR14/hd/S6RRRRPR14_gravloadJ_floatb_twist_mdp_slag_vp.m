% Calculate Gravitation load on the joints for
% S6RRRRPR14
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
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:34
% EndTime: 2019-03-10 00:25:40
% DurationCPUTime: 2.29s
% Computational Cost: add. (1028->200), mult. (2706->333), div. (0->0), fcn. (3462->16), ass. (0->76)
t223 = cos(pkin(6));
t227 = sin(qJ(2));
t209 = t223 * t227;
t228 = sin(qJ(1));
t230 = cos(qJ(2));
t231 = cos(qJ(1));
t162 = t231 * t209 + t228 * t230;
t177 = sin(qJ(3));
t229 = cos(qJ(3));
t210 = t223 * t230;
t188 = -t231 * t210 + t228 * t227;
t220 = sin(pkin(7));
t221 = sin(pkin(6));
t201 = t221 * t220;
t222 = cos(pkin(7));
t239 = t188 * t222 + t231 * t201;
t140 = -t162 * t229 + t239 * t177;
t176 = sin(qJ(4));
t178 = cos(qJ(4));
t202 = t222 * t221;
t232 = t188 * t220 - t231 * t202;
t128 = t140 * t178 - t176 * t232;
t137 = t162 * t177 + t239 * t229;
t173 = pkin(13) + qJ(6);
t171 = sin(t173);
t172 = cos(t173);
t243 = t128 * t171 + t137 * t172;
t242 = t128 * t172 - t137 * t171;
t127 = t140 * t176 + t178 * t232;
t185 = t228 * t210 + t231 * t227;
t235 = t185 * t222 - t228 * t201;
t234 = t230 * t202 + t223 * t220;
t233 = MDP(24) - MDP(27);
t219 = t171 * t178;
t218 = t172 * t178;
t174 = sin(pkin(13));
t217 = t174 * t178;
t175 = cos(pkin(13));
t216 = t175 * t178;
t215 = pkin(9) * t221;
t214 = pkin(10) * t220;
t213 = t176 * t220;
t212 = t177 * t222;
t211 = t178 * t220;
t208 = t222 * t229;
t207 = t230 * t221;
t206 = t221 * t227;
t163 = -t228 * t209 + t231 * t230;
t142 = t163 * t229 - t177 * t235;
t179 = -t185 * t220 - t228 * t202;
t129 = t142 * t176 + t179 * t178;
t153 = t177 * t234 + t229 * t206;
t183 = -t230 * t201 + t223 * t222;
t135 = t153 * t176 - t183 * t178;
t199 = g(1) * t129 - g(2) * t127 + g(3) * t135;
t194 = t227 * t202;
t190 = t227 * t201;
t160 = -t177 * t194 + t229 * t207;
t159 = t177 * t207 + t229 * t194;
t152 = t177 * t206 - t229 * t234;
t148 = t160 * t178 + t176 * t190;
t147 = t160 * t176 - t178 * t190;
t146 = -t163 * t212 - t185 * t229;
t145 = t163 * t208 - t185 * t177;
t144 = -t162 * t212 - t188 * t229;
t143 = t162 * t208 - t188 * t177;
t141 = t163 * t177 + t229 * t235;
t136 = t153 * t178 + t183 * t176;
t134 = t146 * t178 + t163 * t213;
t133 = t146 * t176 - t163 * t211;
t132 = t144 * t178 + t162 * t213;
t131 = t144 * t176 - t162 * t211;
t130 = t142 * t178 - t179 * t176;
t123 = t130 * t172 + t141 * t171;
t122 = -t130 * t171 + t141 * t172;
t1 = [(g(1) * t228 - g(2) * t231) * MDP(2) + (g(1) * t231 + g(2) * t228) * MDP(3) + (g(1) * t162 - g(2) * t163) * MDP(9) + (-g(1) * t188 + g(2) * t185) * MDP(10) + (-g(1) * t140 - g(2) * t142) * MDP(16) + (-g(1) * t137 + g(2) * t141) * MDP(17) + (-g(1) * t128 - g(2) * t130) * MDP(23) + (-g(1) * (t128 * t175 - t137 * t174) - g(2) * (t130 * t175 + t141 * t174)) * MDP(25) + (-g(1) * (-t128 * t174 - t137 * t175) - g(2) * (-t130 * t174 + t141 * t175)) * MDP(26) + (-g(1) * (-t228 * pkin(1) - t162 * pkin(2) + t140 * pkin(3) + t128 * pkin(4) - pkin(11) * t137 + t127 * qJ(5) + t231 * t215) - g(2) * (t231 * pkin(1) + t163 * pkin(2) + t142 * pkin(3) + t130 * pkin(4) + t141 * pkin(11) + t129 * qJ(5) + t228 * t215) + (g(1) * t232 + g(2) * t179) * pkin(10)) * MDP(28) + (-g(1) * t242 - g(2) * t123) * MDP(34) + (g(1) * t243 - g(2) * t122) * MDP(35) + t233 * (g(1) * t127 + g(2) * t129); (g(1) * t185 + g(2) * t188 - g(3) * t207) * MDP(9) + (g(1) * t163 + g(2) * t162 + g(3) * t206) * MDP(10) + (-g(1) * t146 - g(2) * t144 - g(3) * t160) * MDP(16) + (g(1) * t145 + g(2) * t143 + g(3) * t159) * MDP(17) + (-g(1) * t134 - g(2) * t132 - g(3) * t148) * MDP(23) + (-g(1) * (t134 * t175 + t145 * t174) - g(2) * (t132 * t175 + t143 * t174) - g(3) * (t148 * t175 + t159 * t174)) * MDP(25) + (-g(1) * (-t134 * t174 + t145 * t175) - g(2) * (-t132 * t174 + t143 * t175) - g(3) * (-t148 * t174 + t159 * t175)) * MDP(26) + (-g(1) * (-t185 * pkin(2) + t146 * pkin(3) + t134 * pkin(4) + t145 * pkin(11) + t133 * qJ(5) + t163 * t214) - g(2) * (-t188 * pkin(2) + t144 * pkin(3) + t132 * pkin(4) + t143 * pkin(11) + t131 * qJ(5) + t162 * t214) - g(3) * (pkin(2) * t207 + t160 * pkin(3) + t148 * pkin(4) + pkin(10) * t190 + t159 * pkin(11) + t147 * qJ(5))) * MDP(28) + (-g(1) * (t134 * t172 + t145 * t171) - g(2) * (t132 * t172 + t143 * t171) - g(3) * (t148 * t172 + t159 * t171)) * MDP(34) + (-g(1) * (-t134 * t171 + t145 * t172) - g(2) * (-t132 * t171 + t143 * t172) - g(3) * (-t148 * t171 + t159 * t172)) * MDP(35) + t233 * (g(1) * t133 + g(2) * t131 + g(3) * t147); (-g(1) * (-t141 * t216 + t142 * t174) - g(2) * (-t137 * t216 - t140 * t174) - g(3) * (-t152 * t216 + t153 * t174)) * MDP(25) + (-g(1) * (t141 * t217 + t142 * t175) - g(2) * (t137 * t217 - t140 * t175) - g(3) * (t152 * t217 + t153 * t175)) * MDP(26) + (-g(1) * (-t141 * t218 + t142 * t171) - g(2) * (-t137 * t218 - t140 * t171) - g(3) * (-t152 * t218 + t153 * t171)) * MDP(34) + (-g(1) * (t141 * t219 + t142 * t172) - g(2) * (t137 * t219 - t140 * t172) - g(3) * (t152 * t219 + t153 * t172)) * MDP(35) + (-pkin(11) * MDP(28) + MDP(17)) * (g(1) * t142 - g(2) * t140 + g(3) * t153) + (MDP(28) * (pkin(4) * t178 + qJ(5) * t176 + pkin(3)) + t178 * MDP(23) + MDP(16) - t233 * t176) * (g(1) * t141 + g(2) * t137 + g(3) * t152); (-g(1) * (-pkin(4) * t129 + qJ(5) * t130) - g(2) * (pkin(4) * t127 - qJ(5) * t128) - g(3) * (-pkin(4) * t135 + qJ(5) * t136)) * MDP(28) + t233 * (g(1) * t130 - g(2) * t128 + g(3) * t136) + (MDP(25) * t175 - MDP(26) * t174 + MDP(34) * t172 - MDP(35) * t171 + MDP(23)) * t199; -t199 * MDP(28); (-g(1) * t122 - g(2) * t243 - g(3) * (-t136 * t171 + t152 * t172)) * MDP(34) + (g(1) * t123 - g(2) * t242 - g(3) * (-t136 * t172 - t152 * t171)) * MDP(35);];
taug  = t1;
