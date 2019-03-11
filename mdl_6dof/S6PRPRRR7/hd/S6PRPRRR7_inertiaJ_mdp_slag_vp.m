% Calculate joint inertia matrix for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR7_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:58
% EndTime: 2019-03-08 20:56:00
% DurationCPUTime: 1.01s
% Computational Cost: add. (1619->218), mult. (4493->356), div. (0->0), fcn. (5313->16), ass. (0->108)
t184 = sin(qJ(6));
t188 = cos(qJ(6));
t193 = -(t184 * MDP(28) + t188 * MDP(29)) * pkin(12) + t184 * MDP(25) + t188 * MDP(26);
t176 = sin(pkin(14));
t178 = sin(pkin(7));
t180 = cos(pkin(14));
t222 = t178 * t180;
t182 = cos(pkin(7));
t234 = pkin(2) * t182;
t161 = qJ(3) * t222 + t176 * t234;
t177 = sin(pkin(8));
t181 = cos(pkin(8));
t220 = t180 * t181;
t207 = t178 * t220;
t143 = (t177 * t182 + t207) * pkin(10) + t161;
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t169 = t180 * t234;
t225 = t176 * t178;
t148 = pkin(3) * t182 + t169 + (-pkin(10) * t181 - qJ(3)) * t225;
t152 = (-pkin(10) * t176 * t177 - pkin(3) * t180 - pkin(2)) * t178;
t203 = t148 * t181 + t152 * t177;
t129 = -t186 * t143 + t203 * t190;
t237 = 0.2e1 * MDP(28);
t236 = 0.2e1 * MDP(29);
t172 = t178 ^ 2;
t235 = pkin(2) * t172;
t233 = pkin(11) * t184;
t232 = pkin(11) * t188;
t189 = cos(qJ(5));
t231 = pkin(11) * t189;
t230 = pkin(4) * MDP(21);
t229 = pkin(4) * MDP(22);
t135 = -t148 * t177 + t181 * t152;
t223 = t177 * t190;
t146 = -t182 * t223 + t186 * t225 - t190 * t207;
t224 = t177 * t186;
t147 = t182 * t224 + (t176 * t190 + t186 * t220) * t178;
t126 = pkin(4) * t146 - pkin(11) * t147 + t135;
t130 = t143 * t190 + t203 * t186;
t157 = t177 * t222 - t181 * t182;
t128 = -pkin(11) * t157 + t130;
t185 = sin(qJ(5));
t121 = t126 * t189 - t128 * t185;
t117 = -pkin(5) * t146 - t121;
t228 = t117 * t184;
t227 = t117 * t188;
t179 = sin(pkin(6));
t187 = sin(qJ(2));
t191 = cos(qJ(2));
t219 = t182 * t191;
t183 = cos(pkin(6));
t221 = t178 * t183;
t144 = t180 * t221 + (-t176 * t187 + t180 * t219) * t179;
t226 = t144 * t181;
t218 = t184 * t185;
t217 = t185 * t188;
t139 = t147 * t189 - t157 * t185;
t133 = t139 * t184 - t188 * t146;
t216 = t133 * MDP(26);
t134 = t139 * t188 + t146 * t184;
t215 = t134 * MDP(23);
t214 = t134 * MDP(25);
t138 = t147 * t185 + t189 * t157;
t213 = t138 * MDP(27);
t212 = t139 * MDP(17);
t211 = t139 * MDP(18);
t210 = t146 * MDP(20);
t209 = t188 * MDP(23);
t208 = t189 * MDP(27);
t206 = MDP(24) * t184 * t188;
t205 = -pkin(11) * MDP(21) + MDP(18);
t204 = -pkin(11) * MDP(22) + MDP(19);
t122 = t126 * t185 + t128 * t189;
t202 = t188 * MDP(25) - t184 * MDP(26);
t165 = -pkin(5) * t189 - pkin(12) * t185 - pkin(4);
t154 = t165 * t188 - t184 * t231;
t155 = t165 * t184 + t188 * t231;
t200 = t154 * MDP(28) - t155 * MDP(29);
t199 = MDP(28) * t188 - MDP(29) * t184;
t197 = t189 * MDP(21) - t185 * MDP(22) + MDP(14);
t196 = -MDP(17) + t202;
t195 = -MDP(21) - t199;
t194 = -t180 * MDP(5) + t176 * MDP(6) - pkin(2) * MDP(8);
t127 = pkin(4) * t157 - t129;
t118 = pkin(12) * t146 + t122;
t123 = pkin(5) * t138 - pkin(12) * t139 + t127;
t115 = -t118 * t184 + t123 * t188;
t116 = t118 * t188 + t123 * t184;
t192 = t115 * MDP(28) - t116 * MDP(29) + t213 + t214 - t216;
t175 = t188 ^ 2;
t174 = t185 ^ 2;
t173 = t184 ^ 2;
t163 = t181 * t185 + t189 * t224;
t162 = -t189 * t181 + t185 * t224;
t160 = -qJ(3) * t225 + t169;
t159 = -t178 * t179 * t191 + t183 * t182;
t150 = t163 * t188 - t184 * t223;
t149 = -t163 * t184 - t188 * t223;
t145 = t179 * t187 * t180 + (t179 * t219 + t221) * t176;
t137 = -t144 * t177 + t159 * t181;
t132 = t145 * t190 + (t159 * t177 + t226) * t186;
t131 = t145 * t186 - t159 * t223 - t190 * t226;
t125 = t132 * t189 + t137 * t185;
t124 = t132 * t185 - t137 * t189;
t120 = t125 * t188 + t131 * t184;
t119 = -t125 * t184 + t131 * t188;
t1 = [MDP(1) + (t144 ^ 2 + t145 ^ 2 + t159 ^ 2) * MDP(8); (t144 * t160 + t145 * t161) * MDP(8) + (t131 * t157 + t137 * t146) * MDP(14) + (t132 * t157 + t137 * t147) * MDP(15) + (-t124 * t146 + t131 * t138) * MDP(21) + (-t125 * t146 + t131 * t139) * MDP(22) + (t119 * t138 + t124 * t133) * MDP(28) + (-t120 * t138 + t124 * t134) * MDP(29) + (t144 * MDP(5) - t145 * MDP(6)) * t182 + (MDP(3) * t191 - MDP(4) * t187) * t179 + ((-t144 * t176 + t145 * t180) * MDP(7) + t194 * t159) * t178; (pkin(2) ^ 2 * t172 + t160 ^ 2 + t161 ^ 2) * MDP(8) + t139 ^ 2 * MDP(16) + t157 ^ 2 * MDP(13) + MDP(2) + (-0.2e1 * t157 * MDP(11) + MDP(9) * t147) * t147 + (-0.2e1 * t133 * MDP(24) + t215) * t134 + (-0.2e1 * t147 * MDP(10) + 0.2e1 * t157 * MDP(12) + t210 + 0.2e1 * t211) * t146 + (-0.2e1 * t146 * MDP(19) - 0.2e1 * t212 + t213 + 0.2e1 * t214 - 0.2e1 * t216) * t138 + 0.2e1 * (t160 * t182 + t180 * t235) * MDP(5) + 0.2e1 * (-t161 * t182 - t176 * t235) * MDP(6) + 0.2e1 * (t130 * t157 + t135 * t147) * MDP(15) + 0.2e1 * (-t129 * t157 + t135 * t146) * MDP(14) + 0.2e1 * (t121 * t146 + t127 * t138) * MDP(21) + 0.2e1 * (-t122 * t146 + t127 * t139) * MDP(22) + (t115 * t138 + t117 * t133) * t237 + (-t116 * t138 + t117 * t134) * t236 + 0.2e1 * (-t160 * t176 + t161 * t180) * MDP(7) * t178; t159 * MDP(8); (t146 * t181 - t157 * t223) * MDP(14) + (t147 * t181 + t157 * t224) * MDP(15) + (-t138 * t223 - t146 * t162) * MDP(21) + (-t139 * t223 - t146 * t163) * MDP(22) + (t133 * t162 + t138 * t149) * MDP(28) + (t134 * t162 - t138 * t150) * MDP(29) + t194 * t178; MDP(8); -t132 * MDP(15) + (-t119 * t189 + t124 * t218) * MDP(28) + (t120 * t189 + t124 * t217) * MDP(29) - t197 * t131; -t139 * t229 + t147 * MDP(11) - t146 * MDP(12) - t157 * MDP(13) + t129 * MDP(14) - t130 * MDP(15) + (t200 - t230) * t138 + (-t127 * MDP(21) + t204 * t146 - t192 + t212) * t189 + (t139 * MDP(16) + t127 * MDP(22) + t134 * t209 + (-t133 * t188 - t134 * t184) * MDP(24) + (pkin(11) * t133 + t228) * MDP(28) + (pkin(11) * t134 + t227) * MDP(29) + t205 * t146 + t196 * t138) * t185; (-t149 * t189 + t162 * t218) * MDP(28) + (t150 * t189 + t162 * t217) * MDP(29) + (-t186 * MDP(15) + t197 * t190) * t177; MDP(13) + (t208 + 0.2e1 * t230) * t189 + (MDP(23) * t175 + MDP(16) - 0.2e1 * t206) * t174 + (-t154 * t189 + t174 * t233) * t237 + (t155 * t189 + t174 * t232) * t236 + 0.2e1 * (-t189 * t196 - t229) * t185; -t125 * MDP(22) + t195 * t124; t211 + t210 + t121 * MDP(21) - t122 * MDP(22) + t184 * t215 + (-t133 * t184 + t134 * t188) * MDP(24) + (-pkin(5) * t133 - t227) * MDP(28) + (-pkin(5) * t134 + t228) * MDP(29) + (-MDP(19) + t193) * t138; -MDP(22) * t163 + t195 * t162; (-t193 + t204) * t189 + (t184 * t209 + (-t173 + t175) * MDP(24) + (-pkin(5) * t184 - t232) * MDP(28) + (-pkin(5) * t188 + t233) * MDP(29) + t205) * t185; MDP(23) * t173 + 0.2e1 * pkin(5) * t199 + MDP(20) + 0.2e1 * t206; t119 * MDP(28) - t120 * MDP(29); t192; MDP(28) * t149 - MDP(29) * t150; t202 * t185 + t200 - t208; t193; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
