% Calculate Gravitation load on the joints for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:33
% EndTime: 2019-03-10 02:27:37
% DurationCPUTime: 1.18s
% Computational Cost: add. (746->143), mult. (1468->223), div. (0->0), fcn. (1797->12), ass. (0->65)
t217 = MDP(17) - MDP(33);
t216 = MDP(30) + MDP(32);
t215 = MDP(31) - MDP(34);
t168 = sin(qJ(2));
t171 = cos(qJ(2));
t205 = cos(pkin(6));
t214 = cos(qJ(1));
t185 = t205 * t214;
t213 = sin(qJ(1));
t151 = t168 * t185 + t171 * t213;
t167 = sin(qJ(3));
t170 = cos(qJ(3));
t165 = sin(pkin(6));
t189 = t165 * t214;
t136 = t151 * t170 - t167 * t189;
t150 = t168 * t213 - t171 * t185;
t164 = qJ(4) + qJ(5);
t162 = sin(t164);
t163 = cos(t164);
t119 = t136 * t162 - t150 * t163;
t120 = t136 * t163 + t150 * t162;
t166 = sin(qJ(4));
t169 = cos(qJ(4));
t221 = t136 * t166 - t150 * t169;
t220 = t136 * t169 + t150 * t166;
t184 = t205 * t213;
t153 = -t168 * t184 + t171 * t214;
t188 = t165 * t213;
t140 = t153 * t170 + t167 * t188;
t152 = t168 * t214 + t171 * t184;
t125 = -t140 * t166 + t152 * t169;
t196 = t165 * t168;
t149 = t167 * t205 + t170 * t196;
t195 = t165 * t171;
t219 = g(2) * t221 - g(3) * (-t149 * t166 - t169 * t195) - g(1) * t125;
t218 = -g(1) * t152 - g(2) * t150;
t206 = g(3) * t165;
t200 = t151 * t166;
t199 = t153 * t166;
t198 = t162 * t170;
t197 = t163 * t170;
t194 = t166 * t170;
t193 = t169 * t170;
t192 = t170 * t171;
t191 = t162 * t195;
t190 = pkin(4) * t166 + pkin(9);
t187 = -t151 * t167 - t170 * t189;
t123 = t140 * t162 - t152 * t163;
t133 = t149 * t162 + t163 * t195;
t112 = g(1) * t123 + g(2) * t119 + g(3) * t133;
t124 = t140 * t163 + t152 * t162;
t134 = t149 * t163 - t191;
t186 = t215 * (g(1) * t124 + g(2) * t120 + g(3) * t134) + t216 * t112;
t173 = -g(1) * (-t123 * pkin(5) + qJ(6) * t124) - g(2) * (-t119 * pkin(5) + qJ(6) * t120) - g(3) * (-t133 * pkin(5) + qJ(6) * t134);
t172 = -pkin(11) - pkin(10);
t161 = pkin(4) * t169 + pkin(3);
t142 = (t162 * t168 + t163 * t192) * t165;
t141 = -t163 * t196 + t170 * t191;
t139 = t153 * t167 - t170 * t188;
t131 = -t152 * t197 + t153 * t162;
t130 = -t152 * t198 - t153 * t163;
t129 = -t150 * t197 + t151 * t162;
t128 = -t150 * t198 - t151 * t163;
t126 = t140 * t169 + t152 * t166;
t1 = [(g(1) * t213 - g(2) * t214) * MDP(2) + (g(1) * t214 + g(2) * t213) * MDP(3) + (g(1) * t151 - g(2) * t153) * MDP(9) + (-g(1) * t150 + g(2) * t152) * MDP(10) + (g(1) * t136 - g(2) * t140) * MDP(16) + (g(1) * t220 - g(2) * t126) * MDP(23) + (-g(1) * t221 - g(2) * t125) * MDP(24) + (-g(1) * (-pkin(1) * t213 - t151 * pkin(2) - pkin(5) * t120 + pkin(8) * t189 - qJ(6) * t119 - t136 * t161 - t150 * t190 - t172 * t187) - g(2) * (pkin(1) * t214 + t153 * pkin(2) + t124 * pkin(5) + pkin(8) * t188 + t123 * qJ(6) - t139 * t172 + t140 * t161 + t190 * t152)) * MDP(35) + t215 * (-g(1) * t119 + g(2) * t123) + t217 * (g(1) * t187 + g(2) * t139) + t216 * (g(1) * t120 - g(2) * t124); (g(1) * t153 + g(2) * t151 + g(3) * t196) * MDP(10) + (-g(1) * (-t152 * t193 + t199) - g(2) * (-t150 * t193 + t200) - (t166 * t168 + t169 * t192) * t206) * MDP(23) + (-g(1) * (t152 * t194 + t153 * t169) - g(2) * (t150 * t194 + t151 * t169) - (-t166 * t192 + t168 * t169) * t206) * MDP(24) + (-g(1) * (pkin(4) * t199 + t131 * pkin(5) + t153 * pkin(9) + t130 * qJ(6)) - g(2) * (pkin(4) * t200 + t129 * pkin(5) + t151 * pkin(9) + t128 * qJ(6)) - g(3) * (t142 * pkin(5) + t141 * qJ(6)) - t190 * t168 * t206 + (-t171 * t206 - t218) * (t161 * t170 - t167 * t172 + pkin(2))) * MDP(35) + t215 * (g(1) * t130 + g(2) * t128 + g(3) * t141) + t216 * (-g(1) * t131 - g(2) * t129 - g(3) * t142) + (-t170 * MDP(16) + t217 * t167 - MDP(9)) * (g(3) * t195 + t218); (MDP(35) * t172 + t217) * (g(1) * t140 + g(2) * t136 + g(3) * t149) + (MDP(35) * (pkin(5) * t163 + qJ(6) * t162 + t161) + t216 * t163 - t215 * t162 + MDP(23) * t169 - MDP(24) * t166 + MDP(16)) * (-g(3) * (-t167 * t196 + t170 * t205) - g(2) * t187 + g(1) * t139); t219 * MDP(23) + (g(1) * t126 + g(2) * t220 - g(3) * (-t149 * t169 + t166 * t195)) * MDP(24) + (t219 * pkin(4) + t173) * MDP(35) + t186; MDP(35) * t173 + t186; -t112 * MDP(35);];
taug  = t1;
