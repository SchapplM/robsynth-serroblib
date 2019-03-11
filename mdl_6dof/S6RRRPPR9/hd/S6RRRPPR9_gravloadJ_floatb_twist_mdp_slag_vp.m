% Calculate Gravitation load on the joints for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:57
% EndTime: 2019-03-09 16:18:01
% DurationCPUTime: 1.19s
% Computational Cost: add. (573->147), mult. (1470->235), div. (0->0), fcn. (1801->12), ass. (0->66)
t189 = sin(qJ(2));
t191 = cos(qJ(2));
t227 = cos(pkin(6));
t233 = cos(qJ(1));
t206 = t227 * t233;
t231 = sin(qJ(1));
t168 = t189 * t206 + t191 * t231;
t188 = sin(qJ(3));
t185 = sin(pkin(6));
t216 = t185 * t233;
t232 = cos(qJ(3));
t145 = t168 * t232 - t188 * t216;
t167 = t189 * t231 - t191 * t206;
t184 = sin(pkin(11));
t186 = cos(pkin(11));
t127 = t145 * t184 - t167 * t186;
t128 = t145 * t186 + t167 * t184;
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t238 = t127 * t190 - t128 * t187;
t237 = t127 * t187 + t128 * t190;
t236 = qJ(4) * t188 + pkin(2);
t235 = MDP(17) - MDP(20) - MDP(23);
t221 = MDP(18) + MDP(22);
t234 = MDP(19) - MDP(24);
t215 = t185 * t232;
t144 = t168 * t188 + t215 * t233;
t205 = t227 * t231;
t170 = -t189 * t205 + t191 * t233;
t214 = t185 * t231;
t148 = t170 * t188 - t214 * t232;
t223 = t185 * t189;
t165 = t188 * t223 - t227 * t232;
t196 = g(1) * t148 + g(2) * t144 + g(3) * t165;
t222 = t185 * t191;
t219 = t167 * t232;
t169 = t189 * t233 + t191 * t205;
t218 = t169 * t232;
t217 = t184 * t232;
t213 = t186 * t232;
t212 = t191 * t232;
t211 = -t144 * pkin(3) + qJ(4) * t145;
t149 = t170 * t232 + t188 * t214;
t210 = -t148 * pkin(3) + qJ(4) * t149;
t166 = t188 * t227 + t189 * t215;
t209 = -t165 * pkin(3) + qJ(4) * t166;
t207 = t185 * t212;
t208 = pkin(3) * t207 + pkin(9) * t223 + t236 * t222;
t202 = -pkin(4) * t186 - qJ(5) * t184;
t199 = -pkin(3) * t219 + pkin(9) * t168 - t236 * t167;
t198 = -pkin(3) * t218 + pkin(9) * t170 - t236 * t169;
t193 = pkin(1) * t233 + t170 * pkin(2) + t149 * pkin(3) + pkin(8) * t214 + pkin(9) * t169 + qJ(4) * t148;
t192 = -pkin(1) * t231 - t168 * pkin(2) - pkin(3) * t145 + pkin(8) * t216 - t167 * pkin(9) - qJ(4) * t144;
t151 = (t184 * t189 + t186 * t212) * t185;
t150 = t184 * t207 - t186 * t223;
t143 = t166 * t186 - t184 * t222;
t142 = t166 * t184 + t186 * t222;
t137 = -t169 * t213 + t170 * t184;
t136 = -t169 * t217 - t170 * t186;
t135 = -t167 * t213 + t168 * t184;
t134 = -t167 * t217 - t168 * t186;
t132 = t149 * t186 + t169 * t184;
t131 = t149 * t184 - t169 * t186;
t119 = t131 * t187 + t132 * t190;
t118 = t131 * t190 - t132 * t187;
t1 = [(g(1) * t231 - g(2) * t233) * MDP(2) + (g(1) * t233 + g(2) * t231) * MDP(3) + (g(1) * t168 - g(2) * t170) * MDP(9) + (-g(1) * t167 + g(2) * t169) * MDP(10) + (g(1) * t145 - g(2) * t149) * MDP(16) + (-g(1) * t192 - g(2) * t193) * MDP(21) + (-g(1) * (-pkin(4) * t128 - qJ(5) * t127 + t192) - g(2) * (pkin(4) * t132 + qJ(5) * t131 + t193)) * MDP(25) + (g(1) * t237 - g(2) * t119) * MDP(31) + (g(1) * t238 - g(2) * t118) * MDP(32) + t221 * (g(1) * t128 - g(2) * t132) + t234 * (-g(1) * t127 + g(2) * t131) + t235 * (-g(1) * t144 + g(2) * t148); (g(1) * t170 + g(2) * t168 + g(3) * t223) * MDP(10) + (g(1) * t218 + g(2) * t219 - g(3) * t207) * MDP(16) + (-g(1) * t198 - g(2) * t199 - g(3) * t208) * MDP(21) + (-g(1) * (pkin(4) * t137 + qJ(5) * t136 + t198) - g(2) * (pkin(4) * t135 + qJ(5) * t134 + t199) - g(3) * (pkin(4) * t151 + qJ(5) * t150 + t208)) * MDP(25) + (-g(1) * (t136 * t187 + t137 * t190) - g(2) * (t134 * t187 + t135 * t190) - g(3) * (t150 * t187 + t151 * t190)) * MDP(31) + (-g(1) * (t136 * t190 - t137 * t187) - g(2) * (t134 * t190 - t135 * t187) - g(3) * (t150 * t190 - t151 * t187)) * MDP(32) + t221 * (-g(1) * t137 - g(2) * t135 - g(3) * t151) + t234 * (g(1) * t136 + g(2) * t134 + g(3) * t150) + (t235 * t188 - MDP(9)) * (-g(1) * t169 - g(2) * t167 + g(3) * t222); (-g(1) * t210 - g(2) * t211 - g(3) * t209) * MDP(21) + (-g(1) * (t148 * t202 + t210) - g(2) * (t144 * t202 + t211) - g(3) * (t165 * t202 + t209)) * MDP(25) + t235 * (g(1) * t149 + g(2) * t145 + g(3) * t166) + (MDP(16) + t221 * t186 - t234 * t184 + MDP(32) * (t184 * t190 - t186 * t187) + MDP(31) * (t184 * t187 + t186 * t190)) * t196; -(MDP(21) + MDP(25)) * t196; (-g(1) * t131 - g(2) * t127 - g(3) * t142) * MDP(25); (-g(1) * t118 - g(2) * t238 - g(3) * (t142 * t190 - t143 * t187)) * MDP(31) + (g(1) * t119 + g(2) * t237 - g(3) * (-t142 * t187 - t143 * t190)) * MDP(32);];
taug  = t1;
