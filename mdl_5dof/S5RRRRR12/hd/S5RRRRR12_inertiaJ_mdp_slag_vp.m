% Calculate joint inertia matrix for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR12_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:17
% EndTime: 2019-12-31 22:55:21
% DurationCPUTime: 1.34s
% Computational Cost: add. (1949->271), mult. (5185->397), div. (0->0), fcn. (5816->12), ass. (0->134)
t179 = sin(pkin(6));
t262 = 0.2e1 * t179;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t194 = -(MDP(30) * t183 + MDP(31) * t187) * pkin(11) + t183 * MDP(27) + t187 * MDP(28);
t182 = cos(pkin(5));
t186 = sin(qJ(2));
t180 = sin(pkin(5));
t190 = cos(qJ(2));
t238 = t180 * t190;
t166 = t182 * t186 * pkin(1) + pkin(8) * t238;
t181 = cos(pkin(6));
t237 = t181 * t190;
t210 = t180 * t237;
t145 = (t179 * t182 + t210) * pkin(9) + t166;
t185 = sin(qJ(3));
t189 = cos(qJ(3));
t254 = pkin(1) * t190;
t171 = t182 * t254;
t239 = t180 * t186;
t148 = pkin(2) * t182 + t171 + (-pkin(9) * t181 - pkin(8)) * t239;
t152 = (-pkin(9) * t179 * t186 - pkin(2) * t190 - pkin(1)) * t180;
t206 = t148 * t181 + t152 * t179;
t131 = -t185 * t145 + t189 * t206;
t261 = 2 * MDP(16);
t260 = 2 * MDP(17);
t259 = 2 * MDP(23);
t258 = 2 * MDP(24);
t257 = -2 * MDP(26);
t256 = 0.2e1 * MDP(30);
t255 = 0.2e1 * MDP(31);
t253 = pkin(2) * t185;
t252 = pkin(2) * t189;
t251 = pkin(10) * t183;
t250 = pkin(10) * t187;
t188 = cos(qJ(4));
t249 = pkin(10) * t188;
t248 = MDP(23) * pkin(3);
t247 = pkin(3) * MDP(24);
t135 = -t148 * t179 + t181 * t152;
t240 = t179 * t189;
t146 = -t182 * t240 + t185 * t239 - t189 * t210;
t241 = t179 * t185;
t147 = t182 * t241 + (t185 * t237 + t186 * t189) * t180;
t126 = pkin(3) * t146 - pkin(10) * t147 + t135;
t132 = t145 * t189 + t185 * t206;
t159 = t179 * t238 - t182 * t181;
t128 = -pkin(10) * t159 + t132;
t184 = sin(qJ(4));
t123 = t126 * t188 - t128 * t184;
t121 = -pkin(4) * t146 - t123;
t246 = t121 * t183;
t245 = t121 * t187;
t211 = pkin(9) * t240;
t157 = t211 + (pkin(10) + t253) * t181;
t158 = (-pkin(3) * t189 - pkin(10) * t185 - pkin(2)) * t179;
t141 = -t157 * t184 + t158 * t188;
t139 = pkin(4) * t240 - t141;
t244 = t139 * t183;
t243 = t139 * t187;
t175 = t180 ^ 2;
t242 = t175 * t186;
t236 = t182 * MDP(8);
t235 = MDP(14) * t181;
t137 = t147 * t188 - t159 * t184;
t234 = MDP(19) * t137;
t233 = MDP(22) * t189;
t134 = t137 * t187 + t146 * t183;
t232 = MDP(27) * t134;
t133 = t137 * t183 - t146 * t187;
t231 = MDP(28) * t133;
t136 = t147 * t184 + t159 * t188;
t230 = MDP(29) * t136;
t229 = t134 * MDP(25);
t228 = t137 * MDP(20);
t227 = t146 * MDP(22);
t226 = t147 * MDP(11);
t225 = t147 * MDP(12);
t162 = t181 * t184 + t188 * t241;
t149 = t162 * t183 + t187 * t240;
t224 = t149 * MDP(28);
t150 = t162 * t187 - t183 * t240;
t223 = t150 * MDP(25);
t222 = t150 * MDP(27);
t221 = t159 * MDP(13);
t220 = t159 * MDP(14);
t161 = -t188 * t181 + t184 * t241;
t219 = t161 * MDP(29);
t218 = t162 * MDP(18);
t217 = t162 * MDP(19);
t216 = t162 * MDP(20);
t215 = t181 * MDP(15);
t214 = t185 * MDP(13);
t213 = t187 * MDP(25);
t212 = t188 * MDP(29);
t209 = MDP(26) * t183 * t187;
t208 = pkin(10) * MDP(23) - MDP(20);
t207 = pkin(10) * MDP(24) - MDP(21);
t124 = t126 * t184 + t128 * t188;
t142 = t157 * t188 + t158 * t184;
t170 = pkin(9) * t241;
t163 = t181 * t252 - t170;
t165 = t181 * t253 + t211;
t205 = -t163 * MDP(16) + t165 * MDP(17);
t204 = MDP(27) * t187 - MDP(28) * t183;
t168 = -pkin(4) * t188 - pkin(11) * t184 - pkin(3);
t153 = t168 * t187 - t183 * t249;
t154 = t168 * t183 + t187 * t249;
t202 = MDP(30) * t153 - MDP(31) * t154;
t156 = t170 + (-pkin(3) - t252) * t181;
t200 = -MDP(19) + t204;
t199 = (MDP(6) * t186 + MDP(7) * t190) * t180;
t198 = t202 - t248;
t197 = t141 * MDP(23) - t142 * MDP(24) + t216;
t127 = pkin(3) * t159 - t131;
t196 = t147 * MDP(13) - t159 * MDP(15) + t131 * MDP(16) - t132 * MDP(17);
t195 = t123 * MDP(23) - t124 * MDP(24) + t227 + t228;
t122 = pkin(11) * t146 + t124;
t125 = pkin(4) * t136 - pkin(11) * t137 + t127;
t119 = -t122 * t183 + t125 * t187;
t120 = t122 * t187 + t125 * t183;
t193 = MDP(30) * t119 - MDP(31) * t120 + t230 - t231 + t232;
t138 = pkin(4) * t161 - pkin(11) * t162 + t156;
t140 = -pkin(11) * t240 + t142;
t129 = t138 * t187 - t140 * t183;
t130 = t138 * t183 + t140 * t187;
t192 = t129 * MDP(30) - t130 * MDP(31) + t219 + t222 - t224;
t191 = -MDP(21) + t194;
t178 = t187 ^ 2;
t177 = t184 ^ 2;
t176 = t183 ^ 2;
t174 = t179 ^ 2;
t164 = -pkin(8) * t239 + t171;
t1 = [t159 ^ 2 * MDP(15) + t137 ^ 2 * MDP(18) + MDP(1) + (MDP(4) * t186 + 0.2e1 * MDP(5) * t190) * t242 + (-0.2e1 * t221 + t226) * t147 + (t133 * t257 + t229) * t134 + (0.2e1 * t199 + t236) * t182 + (0.2e1 * t220 - 0.2e1 * t225 + t227 + 0.2e1 * t228) * t146 + (-0.2e1 * MDP(21) * t146 + t230 - 0.2e1 * t231 + 0.2e1 * t232 - 0.2e1 * t234) * t136 + 0.2e1 * (t164 * t182 + t175 * t254) * MDP(9) + 0.2e1 * (-pkin(1) * t242 - t166 * t182) * MDP(10) + (t132 * t159 + t135 * t147) * t260 + (-t131 * t159 + t135 * t146) * t261 + (t123 * t146 + t127 * t136) * t259 + (-t124 * t146 + t127 * t137) * t258 + (t119 * t136 + t121 * t133) * t256 + (-t120 * t136 + t121 * t134) * t255; (t127 * t162 + t137 * t156) * MDP(24) + (t127 * t161 + t136 * t156) * MDP(23) + t236 + (-t136 * t162 - t137 * t161) * MDP(19) + t164 * MDP(9) - t166 * MDP(10) + (t134 * t161 + t136 * t150) * MDP(27) + (t119 * t161 + t121 * t149 + t129 * t136 + t133 * t139) * MDP(30) + (-t133 * t161 - t136 * t149) * MDP(28) + (-t120 * t161 + t121 * t150 - t130 * t136 + t134 * t139) * MDP(31) + (-t133 * t150 - t134 * t149) * MDP(26) + t137 * t218 + t136 * t219 + t134 * t223 + t205 * t159 + t196 * t181 + t199 + (-t161 * MDP(21) + t197 - t235) * t146 + ((-MDP(16) * t146 - MDP(17) * t147) * pkin(2) + (-MDP(12) * t146 + MDP(17) * t135 - t221 + t226) * t185 + (-MDP(16) * t135 + t136 * MDP(21) - t195 - t220 + t225) * t189) * t179; t174 * t185 ^ 2 * MDP(11) + t162 ^ 2 * MDP(18) + MDP(8) + (t214 * t262 + t215) * t181 + (t149 * t257 + t223) * t150 + ((-t216 + t235) * t262 + (0.2e1 * MDP(12) * t185 + t233) * t174) * t189 + (0.2e1 * MDP(21) * t240 - 0.2e1 * t217 + t219 + 0.2e1 * t222 - 0.2e1 * t224) * t161 + (t163 * t181 + t174 * t252) * t261 + (-t165 * t181 - t174 * t253) * t260 + (-t141 * t240 + t156 * t161) * t259 + (t142 * t240 + t156 * t162) * t258 + (t129 * t161 + t139 * t149) * t256 + (-t130 * t161 + t139 * t150) * t255; -t137 * t247 - t146 * MDP(14) + t198 * t136 + (-MDP(23) * t127 - t146 * t207 - t193 + t234) * t188 + (t137 * MDP(18) + t127 * MDP(24) + t134 * t213 + (-t133 * t187 - t134 * t183) * MDP(26) + (pkin(10) * t133 + t246) * MDP(30) + (pkin(10) * t134 + t245) * MDP(31) - t208 * t146 + t200 * t136) * t184 + t196; -t162 * t247 + t215 + (MDP(14) * t189 + t214) * t179 + t198 * t161 + (-t156 * MDP(23) + t207 * t240 - t192 + t217) * t188 + (t218 + t156 * MDP(24) + t150 * t213 + (-t149 * t187 - t150 * t183) * MDP(26) + (pkin(10) * t149 + t244) * MDP(30) + (pkin(10) * t150 + t243) * MDP(31) + t208 * t240 + t200 * t161) * t184 - t205; MDP(15) + (t212 + (2 * t248)) * t188 + (MDP(25) * t178 + MDP(18) - 0.2e1 * t209) * t177 + (-t153 * t188 + t177 * t251) * t256 + (t154 * t188 + t177 * t250) * t255 + 0.2e1 * (-t188 * t200 - t247) * t184; t183 * t229 + (-t133 * t183 + t134 * t187) * MDP(26) + (-pkin(4) * t133 - t245) * MDP(30) + (-pkin(4) * t134 + t246) * MDP(31) + t191 * t136 + t195; -t179 * t233 + t183 * t223 + (-t149 * t183 + t150 * t187) * MDP(26) + (-pkin(4) * t149 - t243) * MDP(30) + (-pkin(4) * t150 + t244) * MDP(31) + t191 * t161 + t197; (-t194 - t207) * t188 + (t183 * t213 + (-t176 + t178) * MDP(26) + (-pkin(4) * t183 - t250) * MDP(30) + (-pkin(4) * t187 + t251) * MDP(31) - t208) * t184; 0.2e1 * t209 + MDP(25) * t176 + MDP(22) + 0.2e1 * (MDP(30) * t187 - MDP(31) * t183) * pkin(4); t193; t192; t184 * t204 + t202 - t212; t194; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
