% Calculate joint inertia matrix for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR15_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR15_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR15_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:26
% EndTime: 2024-09-27 22:26:28
% DurationCPUTime: 1.03s
% Computational Cost: add. (1813->202), mult. (4650->276), div. (0->0), fcn. (5137->12), ass. (0->104)
t199 = cos(pkin(6));
t189 = t199 * MDP(29);
t197 = sin(pkin(6));
t194 = t197 ^ 2;
t201 = sin(qJ(5));
t238 = t194 * t201;
t205 = cos(qJ(5));
t234 = t197 * t205;
t235 = t197 * t201;
t254 = MDP(27) * t235 + MDP(28) * t234;
t203 = sin(qJ(3));
t207 = cos(qJ(3));
t198 = sin(pkin(5));
t208 = cos(qJ(2));
t232 = t198 * t208;
t204 = sin(qJ(2));
t233 = t198 * t204;
t162 = t203 * t233 - t207 * t232;
t163 = (t203 * t208 + t204 * t207) * t198;
t202 = sin(qJ(4));
t206 = cos(qJ(4));
t147 = t206 * t162 + t163 * t202;
t148 = -t162 * t202 + t163 * t206;
t252 = t148 * MDP(20) - t147 * MDP(21);
t251 = t163 * MDP(13) - t162 * MDP(14);
t250 = 2 * MDP(30);
t249 = 2 * MDP(31);
t248 = pkin(8) + pkin(9);
t247 = pkin(1) * t208;
t200 = cos(pkin(5));
t246 = pkin(2) * t200;
t245 = pkin(2) * t203;
t244 = pkin(2) * t207;
t243 = pkin(3) * t200;
t192 = t197 * pkin(11);
t193 = t206 * pkin(3);
t184 = t200 * t247;
t152 = -t248 * t233 + t184 + t246;
t219 = t200 * t204 * pkin(1);
t156 = t248 * t232 + t219;
t137 = t207 * t152 - t156 * t203;
t134 = -pkin(10) * t163 + t137 + t243;
t239 = t156 * t207;
t138 = t152 * t203 + t239;
t135 = -pkin(10) * t162 + t138;
t127 = t206 * t134 - t135 * t202;
t124 = -pkin(11) * t148 * t199 + pkin(4) * t200 + t127;
t171 = (-pkin(2) * t208 - pkin(1)) * t198;
t149 = t162 * pkin(3) + t171;
t129 = t147 * pkin(4) - t148 * t192 + t149;
t121 = -t124 * t197 + t129 * t199;
t242 = t121 * t205;
t215 = -t147 * t199 + t197 * t200;
t132 = t148 * t205 + t215 * t201;
t241 = t132 * t197;
t240 = t135 * t206;
t237 = t194 * t205;
t195 = t198 ^ 2;
t236 = t195 * t204;
t231 = t199 * t201;
t230 = t199 * t205;
t187 = pkin(3) + t244;
t166 = t187 * t202 + t206 * t245;
t158 = t166 + t192;
t178 = t206 * t187;
t165 = -t202 * t245 + t178;
t164 = pkin(4) + t165;
t143 = -t158 * t201 + t164 * t230;
t229 = t143 * t199 + t164 * t237;
t176 = pkin(3) * t202 + t192;
t186 = t193 + pkin(4);
t153 = -t176 * t201 + t186 * t230;
t228 = t153 * t199 + t186 * t237;
t167 = pkin(4) * t230 - pkin(11) * t235;
t227 = pkin(4) * t237 + t167 * t199;
t226 = MDP(24) * t202;
t225 = MDP(25) * t132;
t131 = t147 * t230 + t148 * t201 - t200 * t234;
t224 = MDP(28) * t131;
t139 = -t147 * t197 - t199 * t200;
t223 = MDP(29) * t139;
t222 = t165 * MDP(23);
t221 = t166 * MDP(24);
t220 = t207 * MDP(16);
t218 = t189 + t254;
t128 = t134 * t202 + t240;
t123 = t215 * pkin(11) + t128;
t216 = t124 * t199 + t129 * t197;
t119 = t123 * t205 + t216 * t201;
t217 = -t119 * t199 + t121 * t235;
t214 = MDP(22) + (MDP(25) * t238 + 0.2e1 * MDP(26) * t237) * t201 + (0.2e1 * t254 + t189) * t199;
t213 = (MDP(23) * t206 - t226) * pkin(3);
t212 = (MDP(6) * t204 + MDP(7) * t208) * t198;
t211 = MDP(15) + t214;
t210 = (-t131 * t201 + t132 * t205) * t197 * MDP(26) + (t132 * t199 - t139 * t235) * MDP(27) + (-t131 * t199 - t139 * t234) * MDP(28) + t225 * t235 + t200 * MDP(22) - t139 * t189 + t252;
t209 = t200 * MDP(15) + t210 + t251;
t170 = pkin(8) * t232 + t219;
t169 = pkin(4) * t231 + pkin(11) * t234;
t168 = -pkin(8) * t233 + t184;
t154 = t176 * t205 + t186 * t231;
t144 = t158 * t205 + t164 * t231;
t118 = -t123 * t201 + t216 * t205;
t117 = t118 * t199;
t1 = [MDP(1) + (MDP(4) * t204 + 0.2e1 * MDP(5) * t208) * t236 + (MDP(11) * t163 - 0.2e1 * MDP(12) * t162) * t163 + (MDP(18) * t148 - 0.2e1 * MDP(19) * t147) * t148 + (t223 + 0.2e1 * t224) * t139 + (MDP(8) + MDP(15) + MDP(22)) * t200 ^ 2 + (-0.2e1 * MDP(26) * t131 - 0.2e1 * MDP(27) * t139 + t225) * t132 + 0.2e1 * (t212 + t251 + t252) * t200 + (t119 * t139 + t121 * t132) * t249 + (-t118 * t139 + t121 * t131) * t250 + 0.2e1 * (t127 * t200 + t147 * t149) * MDP(23) + 0.2e1 * (-t128 * t200 + t148 * t149) * MDP(24) + 0.2e1 * (t137 * t200 + t162 * t171) * MDP(16) + 0.2e1 * (-t138 * t200 + t163 * t171) * MDP(17) + 0.2e1 * (-pkin(1) * t236 - t170 * t200) * MDP(10) + 0.2e1 * (t168 * t200 + t195 * t247) * MDP(9); t212 + (t139 * t144 - t164 * t241 + t217) * MDP(31) + t200 * MDP(8) + t168 * MDP(9) - t170 * MDP(10) + t209 + (-t139 * t143 + t117 + (-t131 * t164 - t242) * t197) * MDP(30) + (-t239 + (-t152 - t246) * t203) * MDP(17) + (t200 * t244 + t137) * MDP(16) + (t165 * t200 + t127) * MDP(23) + (-t166 * t200 - t128) * MDP(24); MDP(8) + 0.2e1 * (-MDP(17) * t203 + t220) * pkin(2) + 0.2e1 * t222 - 0.2e1 * t221 + t229 * t250 + (-t144 * t199 - t164 * t238) * t249 + t211; (t139 * t154 - t186 * t241 + t217) * MDP(31) + t137 * MDP(16) - t138 * MDP(17) + t209 + (-t139 * t153 + t117 + (-t131 * t186 - t242) * t197) * MDP(30) + (-t240 + (-t134 - t243) * t202) * MDP(24) + (t200 * t193 + t127) * MDP(23); (t178 + t193) * MDP(23) + (t228 + t229) * MDP(30) + ((-t144 - t154) * t199 + (-t164 - t186) * t238) * MDP(31) + (-pkin(3) - t187) * t226 + (t220 + (-MDP(23) * t202 - MDP(24) * t206 - MDP(17)) * t203) * pkin(2) + t211; 0.2e1 * t213 + t228 * t250 + (-t154 * t199 - t186 * t238) * t249 + t211; t127 * MDP(23) - t128 * MDP(24) + (-t139 * t167 + t117 + (-pkin(4) * t131 - t242) * t197) * MDP(30) + (-pkin(4) * t241 + t139 * t169 + t217) * MDP(31) + t210; t222 - t221 + (t227 + t229) * MDP(30) + ((-t144 - t169) * t199 + (-pkin(4) - t164) * t238) * MDP(31) + t214; (t227 + t228) * MDP(30) + ((-t154 - t169) * t199 + (-pkin(4) - t186) * t238) * MDP(31) + t213 + t214; t227 * t250 + (-pkin(4) * t238 - t169 * t199) * t249 + t214; MDP(27) * t132 + t118 * MDP(30) - t119 * MDP(31) - t223 - t224; MDP(30) * t143 - MDP(31) * t144 + t218; MDP(30) * t153 - MDP(31) * t154 + t218; MDP(30) * t167 - MDP(31) * t169 + t218; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
