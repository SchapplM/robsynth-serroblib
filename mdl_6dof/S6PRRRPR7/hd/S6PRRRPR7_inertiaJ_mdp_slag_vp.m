% Calculate joint inertia matrix for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR7_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:40
% EndTime: 2019-03-08 23:44:43
% DurationCPUTime: 1.20s
% Computational Cost: add. (1553->269), mult. (3789->414), div. (0->0), fcn. (4329->14), ass. (0->115)
t187 = sin(pkin(7));
t253 = 0.2e1 * t187;
t227 = (MDP(22) * qJ(5));
t252 = MDP(21) + t227;
t186 = sin(pkin(13));
t189 = cos(pkin(13));
t192 = sin(qJ(6));
t196 = cos(qJ(6));
t167 = t186 * t192 - t196 * t189;
t168 = t186 * t196 + t189 * t192;
t238 = pkin(11) + qJ(5);
t172 = t238 * t186;
t173 = t238 * t189;
t203 = t168 * MDP(25) - t167 * MDP(26) + (-t172 * t196 - t173 * t192) * MDP(28) - (-t172 * t192 + t173 * t196) * MDP(29);
t209 = MDP(19) * t186 + MDP(20) * t189;
t251 = -t209 * qJ(5) - MDP(15) + t203;
t250 = 2 * MDP(17);
t249 = 0.2e1 * MDP(19);
t248 = 0.2e1 * MDP(20);
t247 = 2 * MDP(21);
t246 = -2 * MDP(24);
t245 = 0.2e1 * MDP(28);
t244 = 0.2e1 * MDP(29);
t194 = sin(qJ(3));
t243 = pkin(2) * t194;
t198 = cos(qJ(3));
t242 = pkin(2) * t198;
t241 = pkin(10) * t186;
t197 = cos(qJ(4));
t240 = pkin(10) * t197;
t193 = sin(qJ(4));
t239 = pkin(11) * t193;
t237 = pkin(4) * MDP(22);
t236 = pkin(10) * MDP(18);
t235 = t187 * t194;
t234 = t187 * t198;
t190 = cos(pkin(7));
t233 = t190 * MDP(9);
t199 = cos(qJ(2));
t232 = t190 * t199;
t231 = t194 * MDP(7);
t176 = pkin(9) * t235;
t155 = t176 + (-pkin(3) - t242) * t190;
t161 = -t197 * t190 + t193 * t235;
t162 = t190 * t193 + t197 * t235;
t135 = pkin(4) * t161 - qJ(5) * t162 + t155;
t213 = pkin(9) * t234;
t156 = t213 + (pkin(10) + t243) * t190;
t157 = (-pkin(3) * t198 - pkin(10) * t194 - pkin(2)) * t187;
t141 = t156 * t197 + t157 * t193;
t138 = -qJ(5) * t234 + t141;
t124 = t186 * t135 + t189 * t138;
t171 = -pkin(4) * t197 - qJ(5) * t193 - pkin(3);
t153 = t186 * t171 + t189 * t240;
t229 = MDP(16) * t198;
t228 = MDP(18) * t193;
t226 = MDP(27) * t197;
t146 = t162 * t186 + t189 * t234;
t147 = t162 * t189 - t186 * t234;
t130 = t196 * t146 + t147 * t192;
t225 = t130 * MDP(26);
t131 = -t146 * t192 + t147 * t196;
t224 = t131 * MDP(25);
t140 = -t193 * t156 + t157 * t197;
t139 = pkin(4) * t234 - t140;
t223 = t139 * MDP(22);
t158 = t168 * t193;
t222 = t158 * MDP(26);
t159 = t167 * t193;
t221 = t159 * MDP(23);
t220 = t159 * MDP(25);
t219 = t161 * MDP(27);
t218 = t162 * MDP(14);
t217 = t167 * MDP(28);
t216 = t168 * MDP(23);
t215 = t186 * MDP(20);
t214 = t189 * MDP(19);
t123 = t189 * t135 - t138 * t186;
t212 = -t123 * t186 + t124 * t189;
t188 = sin(pkin(6));
t191 = cos(pkin(6));
t195 = sin(qJ(2));
t145 = t191 * t235 + (t194 * t232 + t195 * t198) * t188;
t160 = -t187 * t188 * t199 + t190 * t191;
t137 = t145 * t197 + t160 * t193;
t144 = -t191 * t234 + (t194 * t195 - t198 * t232) * t188;
t125 = -t137 * t186 + t144 * t189;
t126 = t137 * t189 + t144 * t186;
t119 = t125 * t196 - t126 * t192;
t120 = t125 * t192 + t126 * t196;
t208 = t119 * MDP(28) - t120 * MDP(29);
t207 = t158 * MDP(28) - t159 * MDP(29);
t206 = -t214 + t215 - t237;
t205 = MDP(22) * pkin(10) + t209;
t204 = t146 * MDP(19) + t147 * MDP(20) + t223;
t201 = t168 * MDP(29) + t206 + t217;
t185 = t193 ^ 2;
t183 = t187 ^ 2;
t181 = -pkin(5) * t189 - pkin(4);
t169 = (pkin(5) * t186 + pkin(10)) * t193;
t166 = t189 * t171;
t164 = t190 * t243 + t213;
t163 = t190 * t242 - t176;
t152 = -t186 * t240 + t166;
t148 = -t186 * t239 + t153;
t142 = -t189 * t239 + t166 + (-pkin(5) - t241) * t197;
t136 = t145 * t193 - t160 * t197;
t129 = t142 * t192 + t148 * t196;
t128 = t142 * t196 - t148 * t192;
t127 = pkin(5) * t146 + t139;
t122 = -pkin(11) * t146 + t124;
t121 = pkin(5) * t161 - pkin(11) * t147 + t123;
t118 = t121 * t192 + t122 * t196;
t117 = t121 * t196 - t122 * t192;
t1 = [MDP(1) + (t125 ^ 2 + t126 ^ 2 + t136 ^ 2) * MDP(22); (-t144 * t190 - t160 * t234) * MDP(10) + (-t145 * t190 + t160 * t235) * MDP(11) + (t136 * t234 + t144 * t161) * MDP(17) + (t137 * t234 + t144 * t162) * MDP(18) + (t125 * t161 + t136 * t146) * MDP(19) + (-t126 * t161 + t136 * t147) * MDP(20) + (-t125 * t147 - t126 * t146) * MDP(21) + (t123 * t125 + t124 * t126 + t136 * t139) * MDP(22) + (t119 * t161 + t130 * t136) * MDP(28) + (-t120 * t161 + t131 * t136) * MDP(29) + (t199 * MDP(3) - t195 * MDP(4)) * t188; (t123 ^ 2 + t124 ^ 2 + t139 ^ 2) * MDP(22) + t162 ^ 2 * MDP(12) + t183 * t194 ^ 2 * MDP(5) + MDP(2) + (t231 * t253 + t233) * t190 + (t131 * MDP(23) + t130 * t246) * t131 + ((MDP(8) * t190 - t218) * t253 + (0.2e1 * MDP(6) * t194 + t229) * t183) * t198 + (-0.2e1 * t162 * MDP(13) + 0.2e1 * MDP(15) * t234 + t219 + 0.2e1 * t224 - 0.2e1 * t225) * t161 + 0.2e1 * (t163 * t190 + t183 * t242) * MDP(10) + (-t140 * t234 + t155 * t161) * t250 + 0.2e1 * (t141 * t234 + t155 * t162) * MDP(18) + 0.2e1 * (-t164 * t190 - t183 * t243) * MDP(11) + (t123 * t161 + t139 * t146) * t249 + (t117 * t161 + t127 * t130) * t245 + (-t124 * t161 + t139 * t147) * t248 + (-t118 * t161 + t127 * t131) * t244 + (-t123 * t147 - t124 * t146) * t247; -t145 * MDP(11) + (t125 * t152 + t126 * t153) * MDP(22) + t207 * t136 + (-t125 * MDP(19) + t126 * MDP(20) - t208) * t197 + ((-t125 * t189 - t126 * t186) * MDP(21) + t205 * t136) * t193 + (-MDP(17) * t197 - MDP(10) + t228) * t144; t233 + t163 * MDP(10) - t164 * MDP(11) + (-pkin(3) * t161 - t155 * t197) * MDP(17) + (-t123 * t197 + t152 * t161) * MDP(19) + (t124 * t197 - t153 * t161) * MDP(20) + (-t146 * t153 - t147 * t152) * MDP(21) + (t123 * t152 + t124 * t153) * MDP(22) - t131 * t221 + (t130 * t159 - t131 * t158) * MDP(24) + (-t131 * t197 - t159 * t161) * MDP(25) + (t130 * t197 - t158 * t161) * MDP(26) - t197 * t219 + (-t117 * t197 + t127 * t158 + t128 * t161 + t130 * t169) * MDP(28) + (t118 * t197 - t127 * t159 - t129 * t161 + t131 * t169) * MDP(29) + (t197 * MDP(13) - pkin(3) * MDP(18)) * t162 + (t231 + (MDP(8) + (-MDP(15) + t236) * t197) * t198) * t187 + (t162 * MDP(12) - t161 * MDP(13) - MDP(14) * t234 + t155 * MDP(18) + (-t123 * t189 - t124 * t186) * MDP(21) + t209 * t139 + (MDP(17) * t234 + t204) * pkin(10)) * t193; MDP(9) + t185 * MDP(12) - 0.2e1 * pkin(3) * t228 + (pkin(10) ^ 2 * t185 + t152 ^ 2 + t153 ^ 2) * MDP(22) - (t158 * t246 - t221) * t159 + (0.2e1 * t193 * MDP(13) + pkin(3) * t250 + 0.2e1 * t220 + 0.2e1 * t222 + t226) * t197 + (-t152 * t197 + t185 * t241) * t249 + (pkin(10) * t185 * t189 + t153 * t197) * t248 + (-t128 * t197 + t158 * t169) * t245 + (t129 * t197 - t159 * t169) * t244 + (-t152 * t189 - t153 * t186) * t193 * t247; -t137 * MDP(18) + (-MDP(17) + t201) * t136 + t252 * (-t125 * t186 + t126 * t189); t218 - t187 * t229 + t140 * MDP(17) - t141 * MDP(18) + (-pkin(4) * t146 - t139 * t189) * MDP(19) + (-pkin(4) * t147 + t139 * t186) * MDP(20) + t212 * MDP(21) - pkin(4) * t223 + t131 * t216 + (-t130 * t168 - t131 * t167) * MDP(24) + (t127 * t167 + t130 * t181) * MDP(28) + (t127 * t168 + t131 * t181) * MDP(29) + ((-t146 * t189 + t147 * t186) * MDP(21) + t212 * MDP(22)) * qJ(5) + t251 * t161; -t159 * t216 + (-t158 * t168 + t159 * t167) * MDP(24) + (t158 * t181 + t167 * t169) * MDP(28) + (-t159 * t181 + t168 * t169) * MDP(29) + (-t236 - t251) * t197 + (MDP(14) - t209 * pkin(4) + (-MDP(17) + t206) * pkin(10)) * t193 + t252 * (-t152 * t186 + t153 * t189); 0.2e1 * t181 * t217 + MDP(16) + (0.2e1 * t214 - 0.2e1 * t215 + t237) * pkin(4) + (t167 * t246 + t181 * t244 + t216) * t168 + (t247 + t227) * (t186 ^ 2 + t189 ^ 2) * qJ(5); t136 * MDP(22); t130 * MDP(28) + t131 * MDP(29) + t204; t205 * t193 + t207; t201; MDP(22); t208; t117 * MDP(28) - t118 * MDP(29) + t219 + t224 - t225; t128 * MDP(28) - t129 * MDP(29) - t220 - t222 - t226; t203; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
