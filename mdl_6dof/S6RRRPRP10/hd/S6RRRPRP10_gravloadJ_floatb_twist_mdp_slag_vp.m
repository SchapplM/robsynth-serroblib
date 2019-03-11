% Calculate Gravitation load on the joints for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:37:08
% EndTime: 2019-03-09 17:37:12
% DurationCPUTime: 1.42s
% Computational Cost: add. (688->160), mult. (1405->243), div. (0->0), fcn. (1692->12), ass. (0->65)
t213 = MDP(17) - MDP(20) - MDP(30);
t170 = sin(qJ(2));
t172 = cos(qJ(2));
t203 = cos(pkin(6));
t209 = cos(qJ(1));
t186 = t203 * t209;
t208 = sin(qJ(1));
t146 = t170 * t186 + t172 * t208;
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t166 = sin(pkin(6));
t188 = t166 * t209;
t127 = t146 * t171 - t169 * t188;
t145 = t170 * t208 - t172 * t186;
t164 = pkin(11) + qJ(5);
t161 = sin(t164);
t162 = cos(t164);
t113 = t127 * t161 - t145 * t162;
t114 = t127 * t162 + t145 * t161;
t211 = MDP(27) + MDP(29);
t210 = MDP(28) - MDP(31);
t126 = t146 * t169 + t171 * t188;
t185 = t203 * t208;
t148 = -t170 * t185 + t172 * t209;
t187 = t166 * t208;
t130 = t148 * t169 - t171 * t187;
t196 = t166 * t170;
t143 = t169 * t196 - t171 * t203;
t175 = g(1) * t130 + g(2) * t126 + g(3) * t143;
t204 = g(3) * t166;
t200 = t161 * t171;
t199 = t162 * t171;
t165 = sin(pkin(11));
t198 = t165 * t170;
t197 = t165 * t171;
t195 = t166 * t172;
t167 = cos(pkin(11));
t194 = t167 * t171;
t193 = t171 * t172;
t192 = pkin(2) * t195 + pkin(9) * t196;
t191 = t161 * t195;
t190 = t209 * pkin(1) + t148 * pkin(2) + pkin(8) * t187;
t189 = pkin(4) * t165 + pkin(9);
t182 = pkin(3) * t171 + qJ(4) * t169;
t160 = pkin(4) * t167 + pkin(3);
t168 = -pkin(10) - qJ(4);
t181 = t160 * t171 - t168 * t169;
t180 = -pkin(1) * t208 - t146 * pkin(2) + pkin(8) * t188;
t131 = t148 * t171 + t169 * t187;
t147 = t170 * t209 + t172 * t185;
t117 = t131 * t161 - t147 * t162;
t144 = t169 * t203 + t171 * t196;
t124 = t144 * t161 + t162 * t195;
t178 = g(1) * t117 + g(2) * t113 + g(3) * t124;
t141 = t147 * pkin(2);
t139 = t145 * pkin(2);
t133 = (t161 * t170 + t162 * t193) * t166;
t132 = -t162 * t196 + t171 * t191;
t125 = t144 * t162 - t191;
t123 = -t147 * t199 + t148 * t161;
t122 = -t147 * t200 - t148 * t162;
t121 = -t145 * t199 + t146 * t161;
t120 = -t145 * t200 - t146 * t162;
t118 = t131 * t162 + t147 * t161;
t1 = [(g(1) * t208 - g(2) * t209) * MDP(2) + (g(1) * t209 + g(2) * t208) * MDP(3) + (g(1) * t146 - g(2) * t148) * MDP(9) + (-g(1) * t145 + g(2) * t147) * MDP(10) + (g(1) * t127 - g(2) * t131) * MDP(16) + (-g(1) * (-t127 * t167 - t145 * t165) - g(2) * (t131 * t167 + t147 * t165)) * MDP(18) + (-g(1) * (t127 * t165 - t145 * t167) - g(2) * (-t131 * t165 + t147 * t167)) * MDP(19) + (-g(1) * (-pkin(3) * t127 - t145 * pkin(9) - qJ(4) * t126 + t180) - g(2) * (pkin(3) * t131 + pkin(9) * t147 + qJ(4) * t130 + t190)) * MDP(21) + (-g(1) * (-pkin(5) * t114 - qJ(6) * t113 + t126 * t168 - t127 * t160 - t145 * t189 + t180) - g(2) * (pkin(5) * t118 + qJ(6) * t117 - t130 * t168 + t131 * t160 + t147 * t189 + t190)) * MDP(32) + t210 * (-g(1) * t113 + g(2) * t117) + t211 * (g(1) * t114 - g(2) * t118) + t213 * (-g(1) * t126 + g(2) * t130); (g(1) * t148 + g(2) * t146 + g(3) * t196) * MDP(10) + (-g(1) * (-t147 * t194 + t148 * t165) - g(2) * (-t145 * t194 + t146 * t165) - (t167 * t193 + t198) * t204) * MDP(18) + (-g(1) * (t147 * t197 + t148 * t167) - g(2) * (t145 * t197 + t146 * t167) - (-t165 * t193 + t167 * t170) * t204) * MDP(19) + (-g(1) * (pkin(9) * t148 - t147 * t182 - t141) - g(2) * (pkin(9) * t146 - t145 * t182 - t139) - g(3) * (t182 * t195 + t192)) * MDP(21) + (-g(1) * (pkin(5) * t123 + qJ(6) * t122 - t147 * t181 + t148 * t189 - t141) - g(2) * (pkin(5) * t121 + qJ(6) * t120 - t145 * t181 + t146 * t189 - t139) - g(3) * (t133 * pkin(5) + t132 * qJ(6) + t192) - (pkin(4) * t198 + t172 * t181) * t204) * MDP(32) + t210 * (g(1) * t122 + g(2) * t120 + g(3) * t132) + t211 * (-g(1) * t123 - g(2) * t121 - g(3) * t133) + (-t171 * MDP(16) + t213 * t169 - MDP(9)) * (-g(1) * t147 - g(2) * t145 + g(3) * t195); (-g(1) * (-pkin(3) * t130 + qJ(4) * t131) - g(2) * (-pkin(3) * t126 + qJ(4) * t127) - g(3) * (-pkin(3) * t143 + qJ(4) * t144)) * MDP(21) + (MDP(32) * t168 + t213) * (g(1) * t131 + g(2) * t127 + g(3) * t144) + (MDP(32) * (pkin(5) * t162 + qJ(6) * t161 + t160) + t211 * t162 - t210 * t161 + MDP(18) * t167 - MDP(19) * t165 + MDP(16)) * t175; -(MDP(21) + MDP(32)) * t175; (-g(1) * (-pkin(5) * t117 + qJ(6) * t118) - g(2) * (-pkin(5) * t113 + qJ(6) * t114) - g(3) * (-pkin(5) * t124 + qJ(6) * t125)) * MDP(32) + t211 * t178 + t210 * (g(1) * t118 + g(2) * t114 + g(3) * t125); -t178 * MDP(32);];
taug  = t1;
