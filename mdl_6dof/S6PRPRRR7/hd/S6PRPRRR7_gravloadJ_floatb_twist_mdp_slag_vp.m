% Calculate Gravitation load on the joints for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:56
% EndTime: 2019-03-08 20:55:58
% DurationCPUTime: 0.81s
% Computational Cost: add. (979->144), mult. (2832->270), div. (0->0), fcn. (3707->18), ass. (0->86)
t153 = sin(pkin(8));
t158 = cos(pkin(8));
t152 = sin(pkin(14));
t155 = sin(pkin(6));
t156 = cos(pkin(14));
t163 = sin(qJ(2));
t154 = sin(pkin(7));
t204 = cos(pkin(6));
t186 = t204 * t154;
t159 = cos(pkin(7));
t166 = cos(qJ(2));
t193 = t159 * t166;
t169 = t156 * t186 + (-t152 * t163 + t156 * t193) * t155;
t196 = t155 * t166;
t178 = t154 * t196 - t204 * t159;
t208 = -t178 * t153 + t169 * t158;
t157 = cos(pkin(13));
t203 = sin(pkin(13));
t184 = t204 * t203;
t150 = t157 * t166 - t163 * t184;
t149 = -t157 * t163 - t166 * t184;
t188 = t155 * t203;
t177 = t149 * t159 + t154 * t188;
t171 = t150 * t152 - t177 * t156;
t179 = t149 * t154 - t159 * t188;
t207 = t179 * t153 + t171 * t158;
t187 = t157 * t204;
t148 = t163 * t187 + t203 * t166;
t147 = -t203 * t163 + t166 * t187;
t198 = t155 * t157;
t182 = t147 * t159 - t154 * t198;
t172 = t148 * t152 - t182 * t156;
t183 = t147 * t154 + t159 * t198;
t206 = t183 * t153 + t172 * t158;
t205 = cos(qJ(4));
t202 = qJ(3) * t154;
t201 = t152 * t159;
t200 = t154 * t153;
t199 = t154 * t158;
t197 = t155 * t163;
t195 = t156 * t159;
t194 = t159 * t163;
t160 = sin(qJ(6));
t165 = cos(qJ(5));
t192 = t160 * t165;
t164 = cos(qJ(6));
t191 = t164 * t165;
t190 = t154 * t197;
t189 = t158 * t205;
t185 = t205 * t200;
t162 = sin(qJ(4));
t161 = sin(qJ(5));
t146 = (-t152 * t194 + t156 * t166) * t155;
t145 = (-t152 * t166 - t156 * t194) * t155;
t143 = t156 * t197 + (t155 * t193 + t186) * t152;
t139 = -t145 * t153 + t158 * t190;
t138 = t149 * t156 - t150 * t201;
t137 = -t149 * t152 - t150 * t195;
t136 = t147 * t156 - t148 * t201;
t135 = -t147 * t152 - t148 * t195;
t134 = t150 * t156 + t177 * t152;
t133 = t148 * t156 + t182 * t152;
t132 = -t169 * t153 - t178 * t158;
t129 = t146 * t205 + (t145 * t158 + t153 * t190) * t162;
t128 = -t145 * t189 + t146 * t162 - t185 * t197;
t127 = -t137 * t153 + t150 * t199;
t126 = -t135 * t153 + t148 * t199;
t125 = t171 * t153 - t179 * t158;
t124 = t172 * t153 - t183 * t158;
t123 = t143 * t205 + t208 * t162;
t122 = t143 * t162 - t208 * t205;
t121 = t129 * t165 + t139 * t161;
t120 = t138 * t205 + (t137 * t158 + t150 * t200) * t162;
t119 = -t137 * t189 + t138 * t162 - t150 * t185;
t118 = t136 * t205 + (t135 * t158 + t148 * t200) * t162;
t117 = -t135 * t189 + t136 * t162 - t148 * t185;
t116 = t134 * t205 - t207 * t162;
t115 = t134 * t162 + t207 * t205;
t114 = t133 * t205 - t206 * t162;
t113 = t133 * t162 + t206 * t205;
t112 = t123 * t165 + t132 * t161;
t110 = t120 * t165 + t127 * t161;
t109 = t118 * t165 + t126 * t161;
t108 = t116 * t165 + t125 * t161;
t106 = t114 * t165 + t124 * t161;
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(1) * t149 - g(2) * t147 - g(3) * t196) * MDP(3) + (-g(1) * t138 - g(2) * t136 - g(3) * t146) * MDP(5) + (-g(1) * t137 - g(2) * t135 - g(3) * t145) * MDP(6) + (-g(1) * (pkin(2) * t149 + t150 * t202) - g(2) * (pkin(2) * t147 + t148 * t202) - g(3) * (pkin(2) * t166 + t163 * t202) * t155) * MDP(8) + (-g(1) * t120 - g(2) * t118 - g(3) * t129) * MDP(14) + (g(1) * t119 + g(2) * t117 + g(3) * t128) * MDP(15) + (-g(1) * t110 - g(2) * t109 - g(3) * t121) * MDP(21) + (-g(1) * (-t120 * t161 + t127 * t165) - g(2) * (-t118 * t161 + t126 * t165) - g(3) * (-t129 * t161 + t139 * t165)) * MDP(22) + (-g(1) * (t110 * t164 + t119 * t160) - g(2) * (t109 * t164 + t117 * t160) - g(3) * (t121 * t164 + t128 * t160)) * MDP(28) + (-g(1) * (-t110 * t160 + t119 * t164) - g(2) * (-t109 * t160 + t117 * t164) - g(3) * (-t121 * t160 + t128 * t164)) * MDP(29) + (-t154 * MDP(7) + MDP(4)) * (g(1) * t150 + g(2) * t148 + g(3) * t197); (g(1) * t179 + g(2) * t183 + g(3) * t178) * MDP(8); (g(1) * t116 + g(2) * t114 + g(3) * t123) * MDP(15) + (-g(1) * (-t115 * t191 + t116 * t160) - g(2) * (-t113 * t191 + t114 * t160) - g(3) * (-t122 * t191 + t123 * t160)) * MDP(28) + (-g(1) * (t115 * t192 + t116 * t164) - g(2) * (t113 * t192 + t114 * t164) - g(3) * (t122 * t192 + t123 * t164)) * MDP(29) + (t165 * MDP(21) - MDP(22) * t161 + MDP(14)) * (g(1) * t115 + g(2) * t113 + g(3) * t122); (g(1) * t108 + g(2) * t106 + g(3) * t112) * MDP(22) + (-MDP(28) * t164 + MDP(29) * t160 - MDP(21)) * (g(1) * (-t116 * t161 + t125 * t165) + g(2) * (-t114 * t161 + t124 * t165) + g(3) * (-t123 * t161 + t132 * t165)); (-g(1) * (-t108 * t160 + t115 * t164) - g(2) * (-t106 * t160 + t113 * t164) - g(3) * (-t112 * t160 + t122 * t164)) * MDP(28) + (-g(1) * (-t108 * t164 - t115 * t160) - g(2) * (-t106 * t164 - t113 * t160) - g(3) * (-t112 * t164 - t122 * t160)) * MDP(29);];
taug  = t1;
