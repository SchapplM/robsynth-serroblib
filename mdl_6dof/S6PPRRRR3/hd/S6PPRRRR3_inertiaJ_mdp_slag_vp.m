% Calculate joint inertia matrix for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:43
% EndTime: 2019-03-08 19:11:45
% DurationCPUTime: 0.71s
% Computational Cost: add. (825->192), mult. (2255->312), div. (0->0), fcn. (2694->16), ass. (0->102)
t148 = sin(pkin(8));
t208 = 0.2e1 * t148;
t155 = sin(qJ(6));
t159 = cos(qJ(6));
t164 = -(MDP(25) * t155 + MDP(26) * t159) * pkin(12) + t155 * MDP(22) + t159 * MDP(23);
t207 = 0.2e1 * MDP(25);
t206 = 0.2e1 * MDP(26);
t157 = sin(qJ(4));
t205 = pkin(3) * t157;
t161 = cos(qJ(4));
t204 = pkin(3) * t161;
t203 = pkin(11) * t155;
t202 = pkin(11) * t159;
t160 = cos(qJ(5));
t201 = pkin(11) * t160;
t200 = pkin(4) * MDP(18);
t199 = pkin(4) * MDP(19);
t152 = cos(pkin(8));
t195 = t148 * t161;
t176 = pkin(10) * t195;
t129 = t176 + (pkin(11) + t205) * t152;
t130 = (-pkin(4) * t161 - pkin(11) * t157 - pkin(3)) * t148;
t156 = sin(qJ(5));
t115 = -t156 * t129 + t160 * t130;
t113 = pkin(5) * t195 - t115;
t198 = t113 * t155;
t197 = t113 * t159;
t196 = t148 * t157;
t149 = sin(pkin(7));
t158 = sin(qJ(3));
t194 = t149 * t158;
t162 = cos(qJ(3));
t193 = t149 * t162;
t151 = cos(pkin(14));
t153 = cos(pkin(7));
t192 = t151 * t153;
t191 = t152 * t161;
t190 = t152 * t162;
t189 = t155 * t156;
t188 = t156 * t159;
t187 = t157 * MDP(8);
t186 = MDP(17) * t161;
t134 = t156 * t152 + t160 * t196;
t122 = t155 * t134 + t159 * t195;
t185 = t122 * MDP(23);
t123 = t159 * t134 - t155 * t195;
t184 = t123 * MDP(20);
t183 = t123 * MDP(22);
t133 = -t160 * t152 + t156 * t196;
t182 = t133 * MDP(24);
t181 = t134 * MDP(14);
t180 = t134 * MDP(15);
t179 = t152 * MDP(10);
t178 = t159 * MDP(20);
t177 = t160 * MDP(24);
t175 = t155 * t159 * MDP(21);
t174 = pkin(11) * MDP(18) - MDP(15);
t173 = pkin(11) * MDP(19) - MDP(16);
t116 = t160 * t129 + t156 * t130;
t172 = MDP(22) * t159 - MDP(23) * t155;
t138 = -t160 * pkin(5) - t156 * pkin(12) - pkin(4);
t125 = t159 * t138 - t155 * t201;
t126 = t155 * t138 + t159 * t201;
t170 = t125 * MDP(25) - t126 * MDP(26);
t169 = t159 * MDP(25) - t155 * MDP(26);
t140 = pkin(10) * t196;
t128 = t140 + (-pkin(4) - t204) * t152;
t167 = -t160 * MDP(18) + t156 * MDP(19) - MDP(11);
t166 = -MDP(14) + t172;
t165 = -MDP(18) - t169;
t112 = t133 * pkin(5) - t134 * pkin(12) + t128;
t114 = -pkin(12) * t195 + t116;
t102 = t159 * t112 - t155 * t114;
t103 = t155 * t112 + t159 * t114;
t163 = t102 * MDP(25) - t103 * MDP(26) + t182 + t183 - t185;
t154 = cos(pkin(6));
t150 = sin(pkin(6));
t147 = sin(pkin(14));
t146 = t159 ^ 2;
t145 = t156 ^ 2;
t144 = t155 ^ 2;
t142 = t148 ^ 2;
t136 = t152 * t205 + t176;
t135 = pkin(3) * t191 - t140;
t132 = -t148 * t193 + t152 * t153;
t131 = -t150 * t151 * t149 + t154 * t153;
t121 = t153 * t196 + (t157 * t190 + t158 * t161) * t149;
t120 = -t161 * t149 * t190 - t153 * t195 + t157 * t194;
t119 = t154 * t194 + (t147 * t162 + t158 * t192) * t150;
t118 = t154 * t193 + (-t147 * t158 + t162 * t192) * t150;
t111 = t160 * t121 + t156 * t132;
t110 = t156 * t121 - t160 * t132;
t109 = -t118 * t148 + t131 * t152;
t107 = t159 * t111 + t155 * t120;
t106 = -t155 * t111 + t159 * t120;
t105 = t119 * t161 + (t118 * t152 + t131 * t148) * t157;
t104 = -t118 * t191 + t119 * t157 - t131 * t195;
t101 = t105 * t160 + t109 * t156;
t100 = t105 * t156 - t109 * t160;
t99 = t101 * t159 + t104 * t155;
t98 = -t101 * t155 + t104 * t159;
t1 = [MDP(1) + (t154 ^ 2 + (t147 ^ 2 + t151 ^ 2) * t150 ^ 2) * MDP(2); t154 * MDP(2); MDP(2); t118 * MDP(4) - t119 * MDP(5) + (-t104 * t152 - t109 * t195) * MDP(11) + (-t105 * t152 + t109 * t196) * MDP(12) + (t100 * t195 + t104 * t133) * MDP(18) + (t101 * t195 + t104 * t134) * MDP(19) + (t100 * t122 + t98 * t133) * MDP(25) + (t100 * t123 - t99 * t133) * MDP(26); (-t120 * t152 - t132 * t195) * MDP(11) + (-t121 * t152 + t132 * t196) * MDP(12) + (t110 * t195 + t120 * t133) * MDP(18) + (t111 * t195 + t120 * t134) * MDP(19) + (t106 * t133 + t110 * t122) * MDP(25) + (-t107 * t133 + t110 * t123) * MDP(26) + (t162 * MDP(4) - t158 * MDP(5)) * t149; t142 * t157 ^ 2 * MDP(6) + t134 ^ 2 * MDP(13) + MDP(3) + (t187 * t208 + t179) * t152 + (-0.2e1 * t122 * MDP(21) + t184) * t123 + ((MDP(9) * t152 - t180) * t208 + (0.2e1 * MDP(7) * t157 + t186) * t142) * t161 + (0.2e1 * MDP(16) * t195 - 0.2e1 * t181 + t182 + 0.2e1 * t183 - 0.2e1 * t185) * t133 + 0.2e1 * (t135 * t152 + t142 * t204) * MDP(11) + 0.2e1 * (-t136 * t152 - t142 * t205) * MDP(12) + 0.2e1 * (-t115 * t195 + t128 * t133) * MDP(18) + 0.2e1 * (t116 * t195 + t128 * t134) * MDP(19) + (t102 * t133 + t113 * t122) * t207 + (-t103 * t133 + t113 * t123) * t206; -t105 * MDP(12) + (t100 * t189 - t98 * t160) * MDP(25) + (t100 * t188 + t99 * t160) * MDP(26) + t167 * t104; -t121 * MDP(12) + (-t106 * t160 + t110 * t189) * MDP(25) + (t107 * t160 + t110 * t188) * MDP(26) + t167 * t120; -t134 * t199 + t179 + t135 * MDP(11) - t136 * MDP(12) + (t161 * MDP(9) + t187) * t148 + (t170 - t200) * t133 + (-t128 * MDP(18) + t173 * t195 - t163 + t181) * t160 + (t134 * MDP(13) + t128 * MDP(19) + t123 * t178 + (-t122 * t159 - t123 * t155) * MDP(21) + (pkin(11) * t122 + t198) * MDP(25) + (pkin(11) * t123 + t197) * MDP(26) + t174 * t195 + t166 * t133) * t156; MDP(10) + (t177 + 0.2e1 * t200) * t160 + (t146 * MDP(20) + MDP(13) - 0.2e1 * t175) * t145 + (-t125 * t160 + t145 * t203) * t207 + (t126 * t160 + t145 * t202) * t206 + 0.2e1 * (-t160 * t166 - t199) * t156; -t101 * MDP(19) + t165 * t100; -t111 * MDP(19) + t165 * t110; t180 - t148 * t186 + t115 * MDP(18) - t116 * MDP(19) + t155 * t184 + (-t155 * t122 + t123 * t159) * MDP(21) + (-pkin(5) * t122 - t197) * MDP(25) + (-pkin(5) * t123 + t198) * MDP(26) + (-MDP(16) + t164) * t133; (-t164 - t173) * t160 + (t155 * t178 + (-t144 + t146) * MDP(21) + (-pkin(5) * t155 - t202) * MDP(25) + (-pkin(5) * t159 + t203) * MDP(26) - t174) * t156; t144 * MDP(20) + 0.2e1 * pkin(5) * t169 + MDP(17) + 0.2e1 * t175; t98 * MDP(25) - t99 * MDP(26); t106 * MDP(25) - t107 * MDP(26); t163; t172 * t156 + t170 - t177; t164; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
