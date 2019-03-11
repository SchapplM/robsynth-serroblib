% Calculate Gravitation load on the joints for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR15_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:15
% EndTime: 2019-03-10 00:50:20
% DurationCPUTime: 1.52s
% Computational Cost: add. (880->169), mult. (2421->281), div. (0->0), fcn. (3083->14), ass. (0->73)
t167 = sin(qJ(2));
t170 = cos(qJ(1));
t206 = cos(pkin(6));
t190 = t170 * t206;
t211 = sin(qJ(1));
t213 = cos(qJ(2));
t153 = t211 * t167 - t190 * t213;
t154 = t167 * t190 + t211 * t213;
t166 = sin(qJ(3));
t205 = cos(pkin(7));
t191 = t166 * t205;
t162 = sin(pkin(7));
t163 = sin(pkin(6));
t201 = t163 * t170;
t196 = t162 * t201;
t212 = cos(qJ(3));
t129 = -t153 * t191 + t154 * t212 - t166 * t196;
t165 = sin(qJ(4));
t169 = cos(qJ(4));
t192 = t163 * t205;
t175 = t153 * t162 - t170 * t192;
t116 = t129 * t165 - t169 * t175;
t187 = t205 * t212;
t128 = t153 * t187 + t154 * t166 + t196 * t212;
t164 = sin(qJ(6));
t168 = cos(qJ(6));
t224 = t116 * t164 + t128 * t168;
t223 = t116 * t168 - t128 * t164;
t222 = MDP(23) - MDP(26);
t117 = t129 * t169 + t165 * t175;
t188 = t206 * t211;
t174 = t170 * t167 + t188 * t213;
t194 = t163 * t211;
t217 = -t162 * t194 + t174 * t205;
t216 = MDP(17) - MDP(25);
t214 = MDP(24) - MDP(27);
t210 = pkin(10) * t162;
t204 = t162 * t165;
t203 = t162 * t169;
t202 = t163 * t167;
t200 = t164 * t165;
t199 = t165 * t168;
t198 = t167 * t162;
t197 = t163 * t198;
t195 = t163 * t213;
t193 = t162 * t206;
t155 = -t167 * t188 + t170 * t213;
t133 = t155 * t212 - t217 * t166;
t171 = -t162 * t174 - t192 * t211;
t120 = t133 * t165 + t169 * t171;
t144 = t166 * t193 + (t212 * t167 + t191 * t213) * t163;
t173 = -t162 * t195 + t205 * t206;
t126 = t144 * t165 - t169 * t173;
t182 = g(1) * t120 + g(2) * t116 + g(3) * t126;
t151 = (-t167 * t191 + t212 * t213) * t163;
t150 = (t166 * t213 + t167 * t187) * t163;
t143 = t166 * t202 - t187 * t195 - t193 * t212;
t139 = t151 * t169 + t165 * t197;
t138 = t151 * t165 - t169 * t197;
t137 = -t155 * t191 - t174 * t212;
t136 = t155 * t187 - t166 * t174;
t135 = -t153 * t212 - t154 * t191;
t134 = -t153 * t166 + t154 * t187;
t132 = t155 * t166 + t217 * t212;
t127 = t144 * t169 + t165 * t173;
t125 = t137 * t169 + t155 * t204;
t124 = t137 * t165 - t155 * t203;
t123 = t135 * t169 + t154 * t204;
t122 = t135 * t165 - t154 * t203;
t121 = t133 * t169 - t165 * t171;
t115 = t120 * t164 + t132 * t168;
t114 = t120 * t168 - t132 * t164;
t1 = [(g(1) * t211 - g(2) * t170) * MDP(2) + (g(1) * t170 + g(2) * t211) * MDP(3) + (g(1) * t154 - g(2) * t155) * MDP(9) + (-g(1) * t153 + g(2) * t174) * MDP(10) + (g(1) * t129 - g(2) * t133) * MDP(16) + (-g(1) * (-pkin(1) * t211 - t154 * pkin(2) - pkin(3) * t129 - pkin(4) * t117 + pkin(9) * t201 - pkin(11) * t128 - qJ(5) * t116) - g(2) * (t170 * pkin(1) + t155 * pkin(2) + t133 * pkin(3) + t121 * pkin(4) + pkin(9) * t194 + t132 * pkin(11) + t120 * qJ(5)) + (g(1) * t175 + g(2) * t171) * pkin(10)) * MDP(28) + (g(1) * t224 - g(2) * t115) * MDP(34) + (g(1) * t223 - g(2) * t114) * MDP(35) + t214 * (-g(1) * t116 + g(2) * t120) - t222 * (-g(1) * t117 + g(2) * t121) + t216 * (-g(1) * t128 + g(2) * t132); (g(1) * t174 + g(2) * t153 - g(3) * t195) * MDP(9) + (g(1) * t155 + g(2) * t154 + g(3) * t202) * MDP(10) + (-g(1) * t137 - g(2) * t135 - g(3) * t151) * MDP(16) + (-g(1) * (-pkin(2) * t174 + t137 * pkin(3) + t125 * pkin(4) + t136 * pkin(11) + t124 * qJ(5) + t155 * t210) - g(2) * (-pkin(2) * t153 + pkin(3) * t135 + pkin(4) * t123 + pkin(11) * t134 + qJ(5) * t122 + t154 * t210) - g(3) * (t151 * pkin(3) + t139 * pkin(4) + t150 * pkin(11) + t138 * qJ(5) + (pkin(2) * t213 + pkin(10) * t198) * t163)) * MDP(28) + (-g(1) * (t124 * t164 + t136 * t168) - g(2) * (t122 * t164 + t134 * t168) - g(3) * (t138 * t164 + t150 * t168)) * MDP(34) + (-g(1) * (t124 * t168 - t136 * t164) - g(2) * (t122 * t168 - t134 * t164) - g(3) * (t138 * t168 - t150 * t164)) * MDP(35) + t214 * (g(1) * t124 + g(2) * t122 + g(3) * t138) - t222 * (g(1) * t125 + g(2) * t123 + g(3) * t139) + t216 * (g(1) * t136 + g(2) * t134 + g(3) * t150); (-g(1) * (-t132 * t200 + t133 * t168) - g(2) * (-t128 * t200 + t129 * t168) - g(3) * (-t143 * t200 + t144 * t168)) * MDP(34) + (-g(1) * (-t132 * t199 - t133 * t164) - g(2) * (-t128 * t199 - t129 * t164) - g(3) * (-t143 * t199 - t144 * t164)) * MDP(35) + (-MDP(28) * pkin(11) + t216) * (g(1) * t133 + g(2) * t129 + g(3) * t144) + (MDP(16) + MDP(28) * (pkin(4) * t169 + qJ(5) * t165 + pkin(3)) + t222 * t169 - t214 * t165) * (g(1) * t132 + g(2) * t128 + g(3) * t143); (-g(1) * (-pkin(4) * t120 + qJ(5) * t121) - g(2) * (-pkin(4) * t116 + qJ(5) * t117) - g(3) * (-pkin(4) * t126 + qJ(5) * t127)) * MDP(28) + (-MDP(34) * t164 - MDP(35) * t168 + t214) * (g(1) * t121 + g(2) * t117 + g(3) * t127) + t222 * t182; -t182 * MDP(28); (-g(1) * t114 - g(2) * t223 - g(3) * (t126 * t168 - t143 * t164)) * MDP(34) + (g(1) * t115 + g(2) * t224 - g(3) * (-t126 * t164 - t143 * t168)) * MDP(35);];
taug  = t1;
