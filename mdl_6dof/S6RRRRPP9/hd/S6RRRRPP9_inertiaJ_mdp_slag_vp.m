% Calculate joint inertia matrix for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP9_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:52
% EndTime: 2019-03-09 21:49:57
% DurationCPUTime: 1.88s
% Computational Cost: add. (1933->330), mult. (4151->442), div. (0->0), fcn. (4318->8), ass. (0->107)
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t254 = pkin(5) + pkin(10);
t172 = t254 * t193;
t173 = t254 * t196;
t203 = t193 * MDP(20) + t196 * MDP(21) + t173 * MDP(30) - t172 * MDP(31);
t222 = MDP(24) - MDP(27);
t262 = -pkin(10) * ((MDP(23) - MDP(26)) * t193 + t222 * t196) - MDP(14) + t203;
t194 = sin(qJ(3));
t261 = 0.2e1 * t194;
t189 = sin(pkin(6));
t195 = sin(qJ(2));
t244 = t189 * t195;
t175 = pkin(8) * t244;
t190 = cos(pkin(6));
t198 = cos(qJ(2));
t252 = pkin(1) * t198;
t158 = t175 + (-pkin(2) - t252) * t190;
t197 = cos(qJ(3));
t163 = -t190 * t197 + t194 * t244;
t164 = t190 * t194 + t197 * t244;
t143 = pkin(3) * t163 - pkin(10) * t164 + t158;
t243 = t189 * t198;
t219 = pkin(8) * t243;
t253 = pkin(1) * t195;
t159 = t219 + (pkin(9) + t253) * t190;
t160 = (-pkin(2) * t198 - pkin(9) * t195 - pkin(1)) * t189;
t147 = t159 * t197 + t160 * t194;
t145 = -pkin(10) * t243 + t147;
t137 = t193 * t143 + t196 * t145;
t162 = t163 * qJ(5);
t134 = -t162 - t137;
t136 = t196 * t143 - t193 * t145;
t151 = t164 * t196 - t193 * t243;
t227 = t151 * MDP(20);
t150 = t164 * t193 + t196 * t243;
t229 = t150 * MDP(21);
t260 = -t136 * MDP(23) + t137 * MDP(24) + t134 * MDP(27) - t227 + t229;
t259 = MDP(27) + MDP(30);
t191 = pkin(4) + qJ(6);
t245 = qJ(5) * t193;
t258 = -t191 * t196 - t245;
t256 = 2 * MDP(25);
t255 = 0.2e1 * MDP(31);
t251 = pkin(9) * t197;
t250 = MDP(28) * pkin(10);
t249 = pkin(2) * MDP(16);
t248 = pkin(2) * MDP(17);
t247 = pkin(9) * MDP(17);
t246 = qJ(5) * t150;
t182 = qJ(5) * t196;
t242 = t190 * MDP(8);
t241 = t193 * t194;
t240 = t194 * t196;
t239 = t197 * qJ(5);
t238 = -qJ(6) - t191;
t237 = pkin(4) * t241 + t194 * pkin(9);
t171 = -pkin(3) * t197 - pkin(10) * t194 - pkin(2);
t236 = -t196 * t171 + t193 * t251;
t157 = t193 * t171 + t196 * t251;
t186 = t193 ^ 2;
t188 = t196 ^ 2;
t235 = t186 + t188;
t234 = MDP(15) * t198;
t233 = MDP(19) * t196;
t213 = -pkin(4) * t196 - t245;
t170 = -pkin(3) + t213;
t232 = MDP(28) * t170;
t230 = t150 * MDP(19);
t228 = t151 * MDP(18);
t226 = t163 * MDP(22);
t225 = t164 * MDP(12);
t224 = t164 * MDP(13);
t221 = MDP(25) + MDP(29);
t220 = MDP(26) - MDP(31);
t184 = t197 * pkin(4);
t154 = t184 + t236;
t217 = -0.2e1 * pkin(4) * MDP(26) + MDP(22);
t216 = -pkin(4) * MDP(28) + MDP(26);
t146 = -t194 * t159 + t160 * t197;
t215 = -pkin(5) * t150 + t137;
t214 = pkin(5) * t151 - t136;
t153 = -t157 + t239;
t144 = pkin(3) * t243 - t146;
t135 = -pkin(4) * t163 - t136;
t212 = -t134 * t196 + t135 * t193;
t210 = t196 * MDP(20) - t193 * MDP(21);
t209 = -t236 * MDP(23) - t157 * MDP(24);
t208 = -MDP(30) * t193 - MDP(31) * t196;
t207 = -qJ(5) * t151 + t144;
t167 = -pkin(3) + t258;
t206 = MDP(23) * pkin(3) + MDP(26) * t170 + MDP(29) * t173 - MDP(31) * t167;
t205 = -MDP(24) * pkin(3) - MDP(27) * t170 + MDP(29) * t172 - MDP(30) * t167;
t148 = pkin(5) * t240 + qJ(6) * t197 + t154;
t149 = -pkin(5) * t241 - t153;
t201 = -MDP(26) * t154 + MDP(27) * t153 - MDP(30) * t149 + MDP(31) * t148 - t209;
t199 = qJ(5) ^ 2;
t185 = t189 ^ 2;
t166 = t190 * t253 + t219;
t165 = t190 * t252 - t175;
t161 = -qJ(5) * t240 + t237;
t152 = (qJ(6) * t193 - t182) * t194 + t237;
t138 = pkin(4) * t150 + t207;
t133 = t191 * t150 + t207;
t132 = t162 + t215;
t131 = -t191 * t163 + t214;
t1 = [(t134 ^ 2 + t135 ^ 2 + t138 ^ 2) * MDP(28) + (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) * MDP(32) + t185 * t195 ^ 2 * MDP(4) + MDP(1) + t164 ^ 2 * MDP(11) + (0.2e1 * MDP(6) * t244 + t242) * t190 + (t228 - 0.2e1 * t230) * t151 + (0.2e1 * (MDP(7) * t190 - t224) * t189 + (0.2e1 * MDP(5) * t195 + t234) * t185) * t198 + (0.2e1 * MDP(14) * t243 - 0.2e1 * t225 + t226 + 0.2e1 * t227 - 0.2e1 * t229) * t163 + 0.2e1 * (-t146 * t243 + t158 * t163) * MDP(16) + 0.2e1 * (t147 * t243 + t158 * t164) * MDP(17) + 0.2e1 * (t165 * t190 + t185 * t252) * MDP(9) + 0.2e1 * (-t166 * t190 - t185 * t253) * MDP(10) + 0.2e1 * (t132 * t163 - t133 * t151) * MDP(30) + 0.2e1 * (t136 * t163 + t144 * t150) * MDP(23) + 0.2e1 * (t135 * t163 - t138 * t150) * MDP(26) + 0.2e1 * (-t134 * t163 - t138 * t151) * MDP(27) + (-t131 * t163 + t133 * t150) * t255 + 0.2e1 * (-t137 * t163 + t144 * t151) * MDP(24) + (t134 * t150 + t135 * t151) * t256 + 0.2e1 * (t131 * t151 - t132 * t150) * MDP(29); t242 + (t131 * t148 + t132 * t149 + t133 * t152) * MDP(32) - t164 * t248 + (t134 * t153 + t135 * t154 + t138 * t161) * MDP(28) + t165 * MDP(9) - t166 * MDP(10) + (MDP(6) * t195 + MDP(7) * t198) * t189 + (t154 * MDP(25) - t161 * MDP(27) + MDP(29) * t148 - t152 * MDP(30)) * t151 + (t153 * MDP(25) - t161 * MDP(26) - MDP(29) * t149 + t152 * MDP(31)) * t150 + (-t201 - t249) * t163 + (t225 - t158 * MDP(16) - t226 - t135 * MDP(26) - t132 * MDP(30) + t131 * MDP(31) + (-MDP(14) + t247) * t243 + t260) * t197 + (-MDP(13) * t243 + t164 * MDP(11) - t163 * MDP(12) + t158 * MDP(17) + (MDP(16) * t243 + t150 * MDP(23) + t151 * MDP(24)) * pkin(9) + (-t151 * MDP(19) - t163 * MDP(21) + t144 * MDP(23) + t134 * MDP(25) - t138 * MDP(26) - t132 * MDP(29) + t133 * MDP(31)) * t193 + (t163 * MDP(20) + t144 * MDP(24) + t135 * MDP(25) - t138 * MDP(27) + t131 * MDP(29) - t133 * MDP(30) + t228 - t230) * t196) * t194; MDP(8) + (t153 ^ 2 + t154 ^ 2 + t161 ^ 2) * MDP(28) + (t148 ^ 2 + t149 ^ 2 + t152 ^ 2) * MDP(32) + ((t153 * t193 + t154 * t196) * MDP(25) + (t148 * t196 - t149 * t193) * MDP(29) + (-t193 * MDP(26) - t196 * MDP(27)) * t161 + (-MDP(30) * t196 + MDP(31) * t193) * t152) * t261 + (0.2e1 * t249 + (MDP(12) - t210) * t261 + 0.2e1 * t201 + t197 * MDP(22)) * t197 + (-0.2e1 * t248 + (MDP(18) * t188 - 0.2e1 * t193 * t233 + MDP(11) + 0.2e1 * (t193 * MDP(23) + t196 * MDP(24)) * pkin(9)) * t194) * t194; t224 - t189 * t234 + t146 * MDP(16) - t147 * MDP(17) + t193 * t228 + (-t150 * t193 + t151 * t196) * MDP(19) + (-pkin(3) * t150 - t144 * t196) * MDP(23) + (-pkin(3) * t151 + t144 * t193) * MDP(24) + t212 * MDP(25) + (t138 * t196 - t150 * t170) * MDP(26) + (-t138 * t193 - t151 * t170) * MDP(27) + t138 * t232 + (t131 * t193 + t132 * t196 - t150 * t173 + t151 * t172) * MDP(29) + (-t133 * t193 - t151 * t167) * MDP(30) + (-t133 * t196 + t150 * t167) * MDP(31) + (t131 * t172 + t132 * t173 + t133 * t167) * MDP(32) + ((-t150 * t196 + t151 * t193) * MDP(25) + t212 * MDP(28)) * pkin(10) + t262 * t163; (t148 * t193 + t149 * t196) * MDP(29) + (t148 * t172 + t149 * t173) * MDP(32) + (MDP(26) * t196 - MDP(27) * t193 + t232) * t161 + (MDP(32) * t167 + t208) * t152 + (-t247 - t262) * t197 + (MDP(13) - pkin(9) * MDP(16) + (-t186 + t188) * MDP(19) + (-MDP(23) * pkin(9) + t205) * t196 + (t196 * MDP(18) + pkin(9) * MDP(24) - t206) * t193) * t194 + (MDP(25) + t250) * (-t153 * t196 + t154 * t193); MDP(15) + t186 * MDP(18) + (t235 * pkin(10) ^ 2 + t170 ^ 2) * MDP(28) + (t167 ^ 2 + t172 ^ 2 + t173 ^ 2) * MDP(32) + t235 * pkin(10) * t256 + 0.2e1 * t206 * t196 + 0.2e1 * (t205 + t233) * t193; (-pkin(4) * t151 - t246) * MDP(25) - t136 * MDP(26) + (-pkin(4) * t135 - qJ(5) * t134) * MDP(28) + (-t151 * t191 - t246) * MDP(29) + (0.2e1 * t162 + t215) * MDP(30) - t214 * MDP(31) + (qJ(5) * t132 - t131 * t191) * MDP(32) + (qJ(5) * MDP(27) + (pkin(4) - t238) * MDP(31) + t217) * t163 - t260; (0.2e1 * t184 + t236) * MDP(26) + (-pkin(4) * t154 - qJ(5) * t153) * MDP(28) - t154 * MDP(31) + (qJ(5) * t149 - t148 * t191) * MDP(32) + (t238 * MDP(31) - MDP(22)) * t197 + (t213 * MDP(25) + MDP(29) * t258 + t208 * pkin(5) + t210) * t194 + t209 + t259 * (-0.2e1 * t239 + t157); (-pkin(4) * t193 + t182) * MDP(25) + (-t191 * t193 + t182) * MDP(29) + (qJ(5) * t173 - t172 * t191) * MDP(32) + ((MDP(28) * qJ(5) - t222) * t196 + (-MDP(23) + t216) * t193) * pkin(10) + t203; (pkin(4) ^ 2 + t199) * MDP(28) + t191 * t255 + (t191 ^ 2 + t199) * MDP(32) + 0.2e1 * t259 * qJ(5) + t217; MDP(28) * t135 + MDP(32) * t131 + t221 * t151 + t220 * t163; MDP(28) * t154 + MDP(32) * t148 - t220 * t197 + t221 * t240; MDP(32) * t172 + (t221 + t250) * t193; -MDP(32) * t191 - MDP(31) + t216; MDP(28) + MDP(32); -t150 * MDP(29) + t163 * MDP(30) + MDP(32) * t132; -MDP(29) * t241 - t197 * MDP(30) + MDP(32) * t149; MDP(29) * t196 + MDP(32) * t173; MDP(32) * qJ(5) + MDP(30); 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
