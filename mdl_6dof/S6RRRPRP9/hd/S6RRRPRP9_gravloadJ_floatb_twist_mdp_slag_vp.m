% Calculate Gravitation load on the joints for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:16
% EndTime: 2019-03-09 17:27:19
% DurationCPUTime: 1.00s
% Computational Cost: add. (383->126), mult. (935->176), div. (0->0), fcn. (1053->8), ass. (0->63)
t171 = sin(qJ(3));
t175 = cos(qJ(3));
t177 = cos(qJ(1));
t173 = sin(qJ(1));
t176 = cos(qJ(2));
t207 = t173 * t176;
t148 = t171 * t207 + t175 * t177;
t204 = t177 * t171;
t206 = t175 * t176;
t149 = t173 * t206 - t204;
t170 = sin(qJ(5));
t174 = cos(qJ(5));
t122 = t148 * t170 + t149 * t174;
t150 = -t173 * t175 + t176 * t204;
t205 = t176 * t177;
t151 = t173 * t171 + t175 * t205;
t128 = t150 * t170 + t151 * t174;
t172 = sin(qJ(2));
t209 = t172 * t175;
t211 = t171 * t174;
t140 = t170 * t209 - t172 * t211;
t191 = -t150 * t174 + t151 * t170;
t224 = -t148 * t174 + t149 * t170;
t184 = g(1) * t191 + g(2) * t224 + g(3) * t140;
t203 = MDP(27) + MDP(29);
t221 = MDP(28) - MDP(31);
t190 = t170 * t171 + t174 * t175;
t225 = t190 * t172;
t233 = t221 * (g(1) * t128 + g(2) * t122 + g(3) * t225) + t203 * t184;
t223 = MDP(16) + MDP(18);
t222 = MDP(17) - MDP(20);
t231 = MDP(10) - MDP(19) + MDP(30);
t230 = -pkin(5) * t191 + qJ(6) * t128;
t229 = -pkin(5) * t224 + qJ(6) * t122;
t195 = g(1) * t177 + g(2) * t173;
t226 = t172 * t195;
t219 = -pkin(3) - pkin(4);
t218 = g(1) * t173;
t165 = t172 * pkin(8);
t167 = t176 * pkin(2);
t212 = t171 * t172;
t210 = t171 * t176;
t208 = t172 * t177;
t202 = -pkin(1) - t167;
t201 = -qJ(4) * t171 - pkin(2);
t200 = -t148 * pkin(3) + qJ(4) * t149;
t199 = -t150 * pkin(3) + qJ(4) * t151;
t198 = pkin(3) * t206 + qJ(4) * t210 + t165 + t167;
t193 = -t149 * pkin(3) + t177 * pkin(7) - qJ(4) * t148;
t192 = -pkin(5) * t140 + qJ(6) * t225;
t188 = t172 * (t170 * t175 - t211);
t181 = t177 * pkin(1) + pkin(2) * t205 + t151 * pkin(3) + t173 * pkin(7) + pkin(8) * t208 + t150 * qJ(4);
t180 = g(1) * t150 + g(2) * t148 + g(3) * t212;
t160 = pkin(8) * t205;
t157 = pkin(8) * t207;
t155 = qJ(4) * t209;
t143 = t190 * t176;
t142 = t170 * t206 - t174 * t210;
t135 = t177 * t225;
t134 = t177 * t188;
t133 = t173 * t225;
t132 = t173 * t188;
t1 = [t195 * MDP(3) + (-g(1) * t193 - g(2) * t181 - (t202 - t165) * t218) * MDP(21) + (-g(1) * (-pkin(4) * t149 - pkin(5) * t122 - qJ(6) * t224 + t193) - g(2) * (t151 * pkin(4) + t128 * pkin(5) - pkin(9) * t208 + qJ(6) * t191 + t181) - ((-pkin(8) + pkin(9)) * t172 + t202) * t218) * MDP(32) + t221 * (-g(1) * t224 + g(2) * t191) - t222 * (g(1) * t148 - g(2) * t150) + (t176 * MDP(9) + MDP(2)) * (-g(2) * t177 + t218) + t223 * (g(1) * t149 - g(2) * t151) + t203 * (g(1) * t122 - g(2) * t128) - t231 * (-g(2) * t208 + t172 * t218); (-g(1) * t160 - g(2) * t157 - g(3) * t198 + (pkin(3) * t175 - t201) * t226) * MDP(21) + (-g(1) * (-t135 * pkin(5) - pkin(9) * t205 - t134 * qJ(6) + t160) - g(2) * (-pkin(5) * t133 - pkin(9) * t207 - qJ(6) * t132 + t157) - g(3) * (pkin(4) * t206 + pkin(5) * t143 + qJ(6) * t142 + t198) + (g(3) * pkin(9) + t195 * (-t175 * t219 - t201)) * t172) * MDP(32) - t221 * (g(1) * t134 + g(2) * t132 - g(3) * t142) + t203 * (g(1) * t135 + g(2) * t133 - g(3) * t143) + t231 * (g(3) * t172 + t176 * t195) + (-t222 * t171 + t223 * t175 + MDP(9)) * (-g(3) * t176 + t226); (-g(1) * t199 - g(2) * t200 - g(3) * (-pkin(3) * t212 + t155)) * MDP(21) + (-g(1) * (-pkin(4) * t150 + t199 - t230) - g(2) * (-pkin(4) * t148 + t200 - t229) - g(3) * (t212 * t219 + t155 - t192)) * MDP(32) + t223 * t180 + t222 * (g(1) * t151 + g(2) * t149 + g(3) * t209) - t233; -(MDP(21) + MDP(32)) * t180; (-g(1) * t230 - g(2) * t229 - g(3) * t192) * MDP(32) + t233; -t184 * MDP(32);];
taug  = t1;
