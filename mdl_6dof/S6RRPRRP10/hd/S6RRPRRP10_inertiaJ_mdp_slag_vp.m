% Calculate joint inertia matrix for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP10_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:42:36
% EndTime: 2019-03-09 12:42:41
% DurationCPUTime: 1.54s
% Computational Cost: add. (2516->281), mult. (5514->394), div. (0->0), fcn. (6213->10), ass. (0->121)
t188 = sin(pkin(11));
t190 = cos(pkin(11));
t193 = sin(qJ(4));
t256 = cos(qJ(4));
t172 = t188 * t193 - t256 * t190;
t173 = t256 * t188 + t193 * t190;
t182 = -pkin(3) * t190 - pkin(2);
t156 = pkin(4) * t172 - pkin(10) * t173 + t182;
t251 = pkin(9) + qJ(3);
t175 = t251 * t188;
t176 = t251 * t190;
t159 = -t193 * t175 + t256 * t176;
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t142 = t192 * t156 + t195 * t159;
t247 = qJ(6) * t172;
t138 = t142 + t247;
t216 = -t195 * t156 + t159 * t192;
t252 = pkin(5) * t172;
t139 = t216 - t252;
t201 = -MDP(27) * t216 - t142 * MDP(28);
t271 = -t139 * MDP(29) + t138 * MDP(31) + t201;
t191 = cos(pkin(6));
t189 = sin(pkin(6));
t196 = cos(qJ(2));
t240 = t189 * t196;
t222 = pkin(8) * t240;
t194 = sin(qJ(2));
t255 = pkin(1) * t194;
t163 = t222 + (qJ(3) + t255) * t191;
t164 = (-pkin(2) * t196 - qJ(3) * t194 - pkin(1)) * t189;
t149 = -t163 * t188 + t190 * t164;
t241 = t189 * t194;
t167 = t188 * t191 + t190 * t241;
t140 = -pkin(3) * t240 - pkin(9) * t167 + t149;
t150 = t190 * t163 + t188 * t164;
t177 = t188 * t241;
t215 = t190 * t191 - t177;
t144 = t215 * pkin(9) + t150;
t134 = t193 * t140 + t256 * t144;
t132 = -pkin(10) * t240 + t134;
t151 = t167 * t193 - t256 * t215;
t152 = t256 * t167 + t193 * t215;
t179 = pkin(8) * t241;
t254 = pkin(1) * t196;
t155 = t177 * pkin(3) + t179 + (t182 - t254) * t191;
t137 = t151 * pkin(4) - t152 * pkin(10) + t155;
t128 = t195 * t132 + t192 * t137;
t217 = t132 * t192 - t195 * t137;
t148 = t152 * t195 - t192 * t240;
t232 = t148 * MDP(24);
t147 = t152 * t192 + t195 * t240;
t234 = t147 * MDP(25);
t270 = -MDP(27) * t217 - t128 * MDP(28) + t232 - t234;
t269 = 0.2e1 * t173;
t268 = t152 * MDP(16);
t267 = t152 * MDP(21);
t266 = t159 * MDP(21);
t265 = MDP(32) * pkin(10) + MDP(30);
t262 = 2 * MDP(13);
t261 = 2 * MDP(20);
t260 = 0.2e1 * MDP(21);
t259 = 0.2e1 * MDP(29);
t258 = 0.2e1 * MDP(30);
t257 = 0.2e1 * MDP(31);
t253 = pkin(5) * t151;
t249 = pkin(2) * MDP(14);
t248 = qJ(6) * t151;
t246 = t147 * t195;
t245 = t148 * t192;
t244 = t148 * t195;
t243 = t173 * t192;
t242 = t173 * t195;
t239 = t191 * MDP(8);
t186 = t192 ^ 2;
t187 = t195 ^ 2;
t237 = t186 + t187;
t236 = MDP(19) * t196;
t213 = -pkin(5) * t195 - qJ(6) * t192;
t174 = -pkin(4) + t213;
t235 = MDP(32) * t174;
t233 = t148 * MDP(22);
t231 = t151 * MDP(26);
t230 = t152 * MDP(17);
t229 = t172 * MDP(26);
t228 = t188 * MDP(12);
t227 = t190 * MDP(11);
t226 = t195 * MDP(22);
t225 = t195 * MDP(23);
t224 = MDP(27) + MDP(29);
t223 = -MDP(28) + MDP(31);
t221 = qJ(3) * t240;
t220 = MDP(6) * t241;
t219 = MDP(18) * t240;
t218 = -MDP(32) * pkin(5) - MDP(29);
t158 = t256 * t175 + t176 * t193;
t133 = t256 * t140 - t193 * t144;
t214 = MDP(27) - t218;
t212 = pkin(5) * t192 - qJ(6) * t195;
t211 = qJ(6) * MDP(32) + t223;
t210 = pkin(4) * MDP(27) - t174 * MDP(29);
t209 = -pkin(4) * MDP(28) - t174 * MDP(31);
t125 = t128 + t248;
t126 = t217 - t253;
t208 = t125 * t195 + t126 * t192;
t207 = t125 * t192 - t126 * t195;
t205 = t138 * t192 - t139 * t195;
t204 = -t149 * t188 + t150 * t190;
t203 = t195 * MDP(24) - t192 * MDP(25);
t202 = t192 * MDP(24) + t195 * MDP(25);
t131 = pkin(4) * t240 - t133;
t200 = t223 * t192 + t224 * t195 + MDP(20);
t199 = -MDP(18) + (-t224 * t192 + t223 * t195) * pkin(10) + t202;
t184 = t189 ^ 2;
t169 = t191 * t255 + t222;
t168 = t191 * t254 - t179;
t166 = t179 + (-pkin(2) - t254) * t191;
t146 = t192 * t147;
t145 = t212 * t173 + t158;
t129 = t147 * pkin(5) - t148 * qJ(6) + t131;
t1 = [(t125 ^ 2 + t126 ^ 2 + t129 ^ 2) * MDP(32) + t152 ^ 2 * MDP(15) + (t149 ^ 2 + t150 ^ 2 + t166 ^ 2) * MDP(14) + t184 * t194 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t220 + t239) * t191 + (-0.2e1 * t147 * MDP(23) + t233) * t148 + (0.2e1 * (MDP(7) * t191 - t230) * t189 + (0.2e1 * MDP(5) * t194 + t236) * t184) * t196 + (0.2e1 * t219 + t231 + 0.2e1 * t232 - 0.2e1 * t234 - 0.2e1 * t268) * t151 + 0.2e1 * (t168 * t191 + t184 * t254) * MDP(9) + 0.2e1 * (-t149 * t240 - t166 * t215) * MDP(11) + 0.2e1 * (t150 * t240 + t166 * t167) * MDP(12) + (t134 * t240 + t152 * t155) * t260 + (-t133 * t240 + t151 * t155) * t261 + 0.2e1 * (-t169 * t191 - t184 * t255) * MDP(10) + (-t149 * t167 + t150 * t215) * t262 + 0.2e1 * (-t128 * t151 + t131 * t148) * MDP(28) + (-t126 * t151 + t129 * t147) * t259 + (-t125 * t147 + t126 * t148) * t258 + (t125 * t151 - t129 * t148) * t257 + 0.2e1 * (t131 * t147 - t151 * t217) * MDP(27); t239 + t220 + (t125 * t138 + t126 * t139 + t129 * t145) * MDP(32) + t168 * MDP(9) - t169 * MDP(10) + t182 * t267 + (-t129 * t242 - t145 * t148) * MDP(31) + (t131 * t242 + t148 * t158) * MDP(28) + (t131 * t243 + t147 * t158) * MDP(27) + (t129 * t243 + t145 * t147) * MDP(29) + (pkin(2) * t215 - t166 * t190 + t188 * t221) * MDP(11) + (-pkin(2) * t167 + t166 * t188 + t190 * t221) * MDP(12) + ((t188 * t167 + t190 * t215) * qJ(3) + t204) * MDP(13) + (-t138 * t147 + t139 * t148) * MDP(30) + (-pkin(2) * t166 + t204 * qJ(3)) * MDP(14) + (t158 * MDP(20) + MDP(7) + t266) * t240 + (t182 * MDP(20) + t242 * MDP(24) - t243 * MDP(25) + t229 + t271) * t151 + (t155 * MDP(20) - t126 * MDP(29) + t125 * MDP(31) + t219 - t268 + t270) * t172 + (-t151 * MDP(16) + t148 * t226 + t152 * MDP(15) + t155 * MDP(21) - t207 * MDP(30) - MDP(17) * t240 + (-t245 - t246) * MDP(23)) * t173; MDP(8) + (t138 ^ 2 + t139 ^ 2 + t145 ^ 2) * MDP(32) + (0.2e1 * t227 - 0.2e1 * t228 + t249) * pkin(2) + (-t205 * MDP(30) + (t192 * MDP(27) + t195 * MDP(28)) * t158 + (t192 * MDP(29) - t195 * MDP(31)) * t145) * t269 + (MDP(14) * qJ(3) + t262) * (t188 ^ 2 + t190 ^ 2) * qJ(3) + (t182 * t260 + (t187 * MDP(22) - 0.2e1 * t192 * t225 + MDP(15)) * t173) * t173 + (t182 * t261 + t229 + (-MDP(16) + t203) * t269 + 0.2e1 * t271) * t172; -t215 * MDP(11) + t167 * MDP(12) + t166 * MDP(14) + t267 + (-t146 - t244) * MDP(30) + t207 * MDP(32) + t200 * t151; -t227 + t228 - t249 + t205 * MDP(32) + (-t237 * MDP(30) + MDP(21)) * t173 + t200 * t172; t237 * MDP(32) + MDP(14); t230 - t189 * t236 + t133 * MDP(20) - t134 * MDP(21) + t192 * t233 + (-t146 + t244) * MDP(23) + (-pkin(4) * t147 - t131 * t195) * MDP(27) + (-pkin(4) * t148 + t131 * t192) * MDP(28) + (-t129 * t195 + t147 * t174) * MDP(29) + t208 * MDP(30) + (-t129 * t192 - t148 * t174) * MDP(31) + t129 * t235 + ((t245 - t246) * MDP(30) + t208 * MDP(32)) * pkin(10) + t199 * t151; -t266 + (-t195 * MDP(29) - t192 * MDP(31) + t235) * t145 + (-t195 * MDP(27) + t192 * MDP(28) - MDP(20)) * t158 + t199 * t172 + (MDP(17) + (-t186 + t187) * MDP(23) + t209 * t195 + (-t210 + t226) * t192) * t173 + t265 * (t138 * t195 + t139 * t192); 0; MDP(19) + t186 * MDP(22) + (t237 * pkin(10) ^ 2 + t174 ^ 2) * MDP(32) + t237 * pkin(10) * t258 + 0.2e1 * t210 * t195 + 0.2e1 * (t209 + t225) * t192; t231 + (-t217 + 0.2e1 * t253) * MDP(29) + (-pkin(5) * t148 - qJ(6) * t147) * MDP(30) + (t128 + 0.2e1 * t248) * MDP(31) + (-pkin(5) * t126 + qJ(6) * t125) * MDP(32) + t270; t229 + (-t216 + 0.2e1 * t252) * MDP(29) + (t142 + 0.2e1 * t247) * MDP(31) + (-pkin(5) * t139 + qJ(6) * t138) * MDP(32) + (t213 * MDP(30) + t203) * t173 + t201; t211 * t192 + t214 * t195; -t212 * MDP(30) + (-t214 * t192 + t211 * t195) * pkin(10) + t202; MDP(26) + pkin(5) * t259 + qJ(6) * t257 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -t151 * MDP(29) + t148 * MDP(30) + t126 * MDP(32); -t172 * MDP(29) + MDP(30) * t242 + t139 * MDP(32); -t195 * MDP(32); t265 * t192; t218; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
