% Calculate Gravitation load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:55
% EndTime: 2019-03-08 22:46:57
% DurationCPUTime: 0.74s
% Computational Cost: add. (472->125), mult. (1001->201), div. (0->0), fcn. (1187->12), ass. (0->71)
t157 = -qJ(5) - pkin(9);
t159 = sin(qJ(3));
t202 = -t157 * t159 + pkin(2);
t161 = cos(qJ(4));
t152 = pkin(4) * t161 + pkin(3);
t162 = cos(qJ(3));
t201 = -t152 * t162 - t202;
t179 = MDP(20) + MDP(24);
t200 = MDP(11) - MDP(19) - MDP(22);
t156 = sin(pkin(6));
t199 = g(3) * t156;
t198 = cos(pkin(6));
t197 = cos(pkin(10));
t196 = sin(pkin(10));
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t175 = t198 * t197;
t138 = t160 * t175 + t196 * t163;
t177 = t156 * t197;
t120 = t138 * t162 - t159 * t177;
t158 = sin(qJ(4));
t195 = t120 * t158;
t174 = t198 * t196;
t140 = -t160 * t174 + t197 * t163;
t176 = t156 * t196;
t122 = t140 * t162 + t159 * t176;
t194 = t122 * t158;
t137 = t196 * t160 - t163 * t175;
t193 = t137 * t161;
t192 = t138 * t158;
t139 = t197 * t160 + t163 * t174;
t191 = t139 * t161;
t190 = t140 * t158;
t155 = qJ(4) + pkin(11);
t153 = sin(t155);
t188 = t153 * t162;
t154 = cos(t155);
t187 = t154 * t162;
t186 = t156 * t160;
t185 = t156 * t163;
t183 = t158 * t160;
t182 = t158 * t162;
t181 = t161 * t162;
t180 = t162 * t163;
t178 = t161 * t185;
t142 = t198 * t159 + t162 * t186;
t173 = -t142 * t158 - t178;
t104 = t120 * t153 - t137 * t154;
t106 = t122 * t153 - t139 * t154;
t115 = t142 * t153 + t154 * t185;
t172 = g(1) * t106 + g(2) * t104 + g(3) * t115;
t119 = t138 * t159 + t162 * t177;
t121 = t140 * t159 - t162 * t176;
t141 = t159 * t186 - t198 * t162;
t171 = g(1) * t121 + g(2) * t119 + g(3) * t141;
t170 = g(1) * t122 + g(2) * t120 + g(3) * t142;
t167 = pkin(8) * t186 + t202 * t185 + (pkin(4) * t183 + t152 * t180) * t156;
t166 = pkin(4) * t192 + pkin(8) * t138 + t201 * t137;
t165 = pkin(4) * t190 + pkin(8) * t140 + t201 * t139;
t132 = pkin(4) * t191;
t130 = pkin(4) * t193;
t124 = (t153 * t160 + t154 * t180) * t156;
t123 = (t153 * t180 - t154 * t160) * t156;
t116 = t142 * t154 - t153 * t185;
t112 = -t139 * t187 + t140 * t153;
t111 = -t139 * t188 - t140 * t154;
t110 = -t137 * t187 + t138 * t153;
t109 = -t137 * t188 - t138 * t154;
t107 = t122 * t154 + t139 * t153;
t105 = t120 * t154 + t137 * t153;
t1 = [(-MDP(1) - t179) * g(3); (g(1) * t140 + g(2) * t138 + g(3) * t186) * MDP(4) + (-g(1) * (-t139 * t181 + t190) - g(2) * (-t137 * t181 + t192) - (t161 * t180 + t183) * t199) * MDP(17) + (-g(1) * (t139 * t182 + t140 * t161) - g(2) * (t137 * t182 + t138 * t161) - (-t158 * t180 + t160 * t161) * t199) * MDP(18) + (-g(1) * t165 - g(2) * t166 - g(3) * t167) * MDP(20) + (-g(1) * t112 - g(2) * t110 - g(3) * t124) * MDP(21) + (-g(1) * t111 - g(2) * t109 - g(3) * t123) * MDP(23) + (-g(1) * (pkin(5) * t112 + qJ(6) * t111 + t165) - g(2) * (t110 * pkin(5) + t109 * qJ(6) + t166) - g(3) * (t124 * pkin(5) + t123 * qJ(6) + t167)) * MDP(24) + (-t162 * MDP(10) + t200 * t159 - MDP(3)) * (-g(1) * t139 - g(2) * t137 + g(3) * t185); t200 * t170 + t179 * (-g(1) * (-t121 * t152 - t122 * t157) - g(2) * (-t119 * t152 - t120 * t157) - g(3) * (-t141 * t152 - t142 * t157)) + ((pkin(5) * t154 + qJ(6) * t153) * MDP(24) + MDP(17) * t161 - MDP(18) * t158 + MDP(21) * t154 + MDP(23) * t153 + MDP(10)) * t171; (-g(1) * (t191 - t194) - g(2) * (t193 - t195) - g(3) * t173) * MDP(17) + (-g(1) * (-t122 * t161 - t139 * t158) - g(2) * (-t120 * t161 - t137 * t158) - g(3) * (-t142 * t161 + t158 * t185)) * MDP(18) + (-g(1) * t132 - g(2) * t130 + (g(3) * t178 + t170 * t158) * pkin(4)) * MDP(20) + t172 * MDP(21) + (-g(1) * t107 - g(2) * t105 - g(3) * t116) * MDP(23) + (-g(1) * (-pkin(4) * t194 - pkin(5) * t106 + qJ(6) * t107 + t132) - g(2) * (-pkin(4) * t195 - pkin(5) * t104 + qJ(6) * t105 + t130) - g(3) * (t173 * pkin(4) - t115 * pkin(5) + t116 * qJ(6))) * MDP(24); -t179 * t171; -t172 * MDP(24);];
taug  = t1;
