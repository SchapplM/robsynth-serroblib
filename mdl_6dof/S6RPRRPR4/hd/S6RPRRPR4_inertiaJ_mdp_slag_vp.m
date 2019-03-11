% Calculate joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:18
% EndTime: 2019-03-09 05:10:20
% DurationCPUTime: 0.63s
% Computational Cost: add. (1444->159), mult. (2691->217), div. (0->0), fcn. (3192->10), ass. (0->89)
t176 = sin(qJ(4));
t159 = t176 * pkin(3) + qJ(5);
t171 = sin(pkin(11));
t173 = cos(pkin(11));
t208 = t171 ^ 2 + t173 ^ 2;
t210 = t208 * t159;
t175 = sin(qJ(6));
t178 = cos(qJ(6));
t146 = t175 * t171 - t178 * t173;
t147 = t178 * t171 + t175 * t173;
t212 = t147 * MDP(28) - t146 * MDP(29);
t186 = t146 * MDP(31) + t147 * MDP(32);
t227 = t173 * MDP(22) - t171 * MDP(23);
t185 = -t227 + t186;
t174 = cos(pkin(10));
t161 = -t174 * pkin(2) - pkin(1);
t172 = sin(pkin(10));
t177 = sin(qJ(3));
t179 = cos(qJ(3));
t198 = -t177 * t172 + t179 * t174;
t136 = -t198 * pkin(3) + t161;
t226 = 0.2e1 * t136;
t225 = 0.2e1 * t161;
t224 = 2 * MDP(24);
t223 = -2 * MDP(27);
t222 = cos(qJ(4));
t221 = pkin(1) * MDP(7);
t220 = t173 * pkin(5);
t166 = t173 * pkin(9);
t218 = pkin(7) + qJ(2);
t217 = pkin(4) * MDP(25);
t148 = t179 * t172 + t177 * t174;
t151 = t218 * t172;
t153 = t218 * t174;
t199 = -t179 * t151 - t177 * t153;
t125 = -t148 * pkin(8) + t199;
t192 = t177 * t151 - t179 * t153;
t126 = t198 * pkin(8) - t192;
t119 = -t222 * t125 + t176 * t126;
t216 = t119 * t173;
t132 = t222 * t148 + t176 * t198;
t215 = t132 * t171;
t214 = t172 * MDP(5);
t213 = t174 * MDP(4);
t131 = t176 * t148 - t222 * t198;
t118 = t131 * pkin(4) - t132 * qJ(5) + t136;
t120 = t176 * t125 + t222 * t126;
t107 = t171 * t118 + t173 * t120;
t209 = t208 * qJ(5);
t122 = t146 * t132;
t206 = MDP(26) * t122;
t121 = t147 * t132;
t205 = t121 * MDP(29);
t204 = t122 * MDP(28);
t203 = t131 * MDP(30);
t162 = -t222 * pkin(3) - pkin(4);
t202 = t162 * MDP(25);
t200 = MDP(19) + (MDP(26) * t147 + t146 * t223) * t147;
t106 = t173 * t118 - t171 * t120;
t197 = t208 * MDP(25);
t196 = -pkin(4) * t132 - qJ(5) * t131;
t195 = t106 * t173 + t107 * t171;
t194 = -t106 * t171 + t107 * t173;
t193 = -t131 * t159 + t132 * t162;
t191 = t198 * MDP(13);
t190 = 0.2e1 * t227;
t189 = t171 * MDP(22) + t173 * MDP(23);
t103 = t131 * pkin(5) - t132 * t166 + t106;
t104 = -pkin(9) * t215 + t107;
t188 = (t178 * t103 - t175 * t104) * MDP(31) - (t175 * t103 + t178 * t104) * MDP(32);
t187 = t121 * MDP(31) - t122 * MDP(32);
t184 = 0.2e1 * t186;
t183 = (t222 * MDP(20) - t176 * MDP(21)) * pkin(3);
t182 = t194 * MDP(24) + (-t147 * t121 + t122 * t146) * MDP(27) - t147 * t206 - t119 * MDP(20) - t120 * MDP(21) + t132 * MDP(17) + (-MDP(18) + t212) * t131;
t160 = -pkin(4) - t220;
t152 = t173 * qJ(5) + t166;
t150 = (-pkin(9) - qJ(5)) * t171;
t149 = t162 - t220;
t144 = t173 * t159 + t166;
t143 = (-pkin(9) - t159) * t171;
t135 = t175 * t150 + t178 * t152;
t134 = t178 * t150 - t175 * t152;
t128 = t175 * t143 + t178 * t144;
t127 = t178 * t143 - t175 * t144;
t115 = t119 * t171;
t110 = pkin(5) * t215 + t119;
t109 = t110 * t147;
t108 = t110 * t146;
t1 = [-t191 * t225 + (t106 ^ 2 + t107 ^ 2 + t119 ^ 2) * MDP(25) + MDP(1) + (MDP(15) * t132 + MDP(21) * t226) * t132 - (t121 * t223 - t206) * t122 + (0.2e1 * t213 - 0.2e1 * t214 + t221) * pkin(1) + (MDP(14) * t225 + MDP(8) * t148 + 0.2e1 * t198 * MDP(9)) * t148 + (-0.2e1 * t132 * MDP(16) + MDP(20) * t226 + t203 - 0.2e1 * t204 - 0.2e1 * t205) * t131 + 0.2e1 * t187 * t110 + 0.2e1 * (t106 * MDP(22) - t107 * MDP(23) + t188) * t131 + 0.2e1 * (-t195 * MDP(24) + t189 * t119) * t132 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t172 ^ 2 + t174 ^ 2) * qJ(2); -t213 + t214 - t221 - t191 + t148 * MDP(14) + t195 * MDP(25) + (-t208 * MDP(24) + MDP(21)) * t132 + (MDP(20) - t185) * t131; MDP(7) + t197; (t193 * t171 - t216) * MDP(22) + (t119 * t162 + t194 * t159) * MDP(25) + (t193 * t173 + t115) * MDP(23) + t182 + t199 * MDP(13) + t198 * MDP(11) + t192 * MDP(14) + t148 * MDP(10) + (t149 * t121 + t127 * t131 + t108) * MDP(31) + (-t149 * t122 - t128 * t131 + t109) * MDP(32); 0; MDP(12) + t210 * t224 + t159 ^ 2 * t197 + (-t190 + t202) * t162 + t149 * t184 + 0.2e1 * t183 + t200; (t196 * t171 - t216) * MDP(22) + (t196 * t173 + t115) * MDP(23) + (-t119 * pkin(4) + t194 * qJ(5)) * MDP(25) + (t160 * t121 + t134 * t131 + t108) * MDP(31) + (-t160 * t122 - t135 * t131 + t109) * MDP(32) + t182; 0; (t209 + t210) * MDP(24) + (-t162 * pkin(4) + qJ(5) * t210) * MDP(25) + t183 + t200 + t227 * (pkin(4) - t162) + t186 * (t149 + t160); t209 * t224 + qJ(5) ^ 2 * t197 + t160 * t184 + (t190 + t217) * pkin(4) + t200; t119 * MDP(25) + t189 * t132 + t187; 0; t185 + t202; t185 - t217; MDP(25); t188 + t203 - t204 - t205; -t186; t127 * MDP(31) - t128 * MDP(32) + t212; t134 * MDP(31) - t135 * MDP(32) + t212; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
