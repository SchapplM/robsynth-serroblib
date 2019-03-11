% Calculate joint inertia matrix for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR4_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:53
% EndTime: 2019-03-09 13:34:57
% DurationCPUTime: 1.25s
% Computational Cost: add. (1946->227), mult. (4628->332), div. (0->0), fcn. (5302->12), ass. (0->119)
t212 = sin(qJ(6));
t216 = cos(qJ(6));
t255 = t212 * MDP(29) + t216 * MDP(30);
t229 = t216 * MDP(32) - t212 * MDP(33);
t278 = -MDP(23) + t255;
t217 = cos(qJ(4));
t208 = sin(pkin(12));
t197 = t208 * pkin(2) + pkin(9);
t266 = pkin(10) + t197;
t188 = t266 * t217;
t213 = sin(qJ(5));
t214 = sin(qJ(4));
t234 = t266 * t214;
t271 = cos(qJ(5));
t161 = t213 * t188 + t234 * t271;
t162 = t188 * t271 - t213 * t234;
t277 = -t161 * MDP(25) - t162 * MDP(26);
t191 = t213 * t217 + t214 * t271;
t187 = t191 * MDP(26);
t190 = t213 * t214 - t217 * t271;
t276 = -t190 * MDP(25) - t187;
t275 = 2 * MDP(19);
t274 = 0.2e1 * MDP(25);
t273 = 0.2e1 * MDP(32);
t272 = 0.2e1 * MDP(33);
t218 = cos(qJ(2));
t270 = pkin(1) * t218;
t210 = cos(pkin(12));
t209 = sin(pkin(6));
t257 = t209 * t218;
t215 = sin(qJ(2));
t258 = t209 * t215;
t176 = t208 * t258 - t210 * t257;
t269 = t176 * pkin(4);
t267 = pkin(8) + qJ(3);
t211 = cos(pkin(6));
t196 = t211 * t270;
t168 = t211 * pkin(2) - t258 * t267 + t196;
t239 = t211 * t215 * pkin(1);
t171 = t257 * t267 + t239;
t154 = t208 * t168 + t210 * t171;
t152 = t211 * pkin(9) + t154;
t177 = (t208 * t218 + t210 * t215) * t209;
t192 = (-pkin(2) * t218 - pkin(1)) * t209;
t156 = t176 * pkin(3) - t177 * pkin(9) + t192;
t136 = -t214 * t152 + t217 * t156;
t167 = t177 * t217 + t211 * t214;
t134 = -t167 * pkin(10) + t136 + t269;
t137 = t217 * t152 + t214 * t156;
t166 = t177 * t214 - t211 * t217;
t135 = -t166 * pkin(10) + t137;
t129 = t271 * t134 - t213 * t135;
t127 = -t176 * pkin(5) - t129;
t265 = t127 * t216;
t149 = t166 * t271 + t213 * t167;
t264 = t149 * t212;
t263 = t149 * t216;
t262 = t161 * t216;
t261 = t191 * t212;
t260 = t191 * t216;
t205 = t209 ^ 2;
t259 = t205 * t215;
t256 = t211 * MDP(8);
t150 = -t213 * t166 + t167 * t271;
t142 = t216 * t150 + t176 * t212;
t254 = MDP(29) * t142;
t141 = t212 * t150 - t176 * t216;
t253 = MDP(30) * t141;
t252 = MDP(31) * t149;
t251 = t142 * MDP(27);
t250 = t150 * MDP(20);
t249 = t150 * MDP(21);
t248 = t166 * MDP(16);
t247 = t167 * MDP(13);
t246 = t176 * MDP(22);
t245 = t176 * MDP(23);
t244 = t190 * MDP(31);
t242 = t216 * MDP(27);
t240 = t217 * MDP(18);
t238 = t271 * pkin(4);
t198 = -t210 * pkin(2) - pkin(3);
t237 = t271 * t135;
t236 = t212 * t216 * MDP(28);
t206 = t212 ^ 2;
t235 = t206 * MDP(27) + MDP(24) + 0.2e1 * t236;
t153 = t210 * t168 - t208 * t171;
t233 = -pkin(5) * t191 - pkin(11) * t190;
t200 = t213 * pkin(4) + pkin(11);
t201 = -t238 - pkin(5);
t232 = -t190 * t200 + t191 * t201;
t193 = -t217 * pkin(4) + t198;
t231 = -t214 * MDP(19) + t240;
t230 = MDP(29) * t216 - MDP(30) * t212;
t228 = t212 * MDP(32) + t216 * MDP(33);
t151 = -t211 * pkin(3) - t153;
t130 = t213 * t134 + t237;
t227 = -MDP(21) + t230;
t226 = (t215 * MDP(6) + t218 * MDP(7)) * t209;
t225 = -t229 * t190 + t276;
t224 = (-t212 * t141 + t142 * t216) * MDP(28) + t212 * t251 + t150 * MDP(22) + t176 * MDP(24) + t278 * t149;
t223 = (MDP(25) * t271 - t213 * MDP(26)) * pkin(4);
t207 = t216 ^ 2;
t222 = t242 * t261 + t277 + ((-t206 + t207) * MDP(28) + MDP(22)) * t191 + t278 * t190;
t140 = t166 * pkin(4) + t151;
t128 = t176 * pkin(11) + t130;
t131 = t149 * pkin(5) - t150 * pkin(11) + t140;
t124 = -t212 * t128 + t216 * t131;
t125 = t216 * t128 + t212 * t131;
t221 = MDP(32) * t124 - MDP(33) * t125 + t252 - t253;
t220 = t214 * MDP(15) + t217 * MDP(16) + (-t214 * MDP(18) - t217 * MDP(19)) * t197;
t183 = pkin(8) * t257 + t239;
t182 = -pkin(8) * t258 + t196;
t158 = t190 * pkin(5) - t191 * pkin(11) + t193;
t157 = t161 * t212;
t144 = t212 * t158 + t216 * t162;
t143 = t216 * t158 - t212 * t162;
t138 = t142 * t190;
t126 = t127 * t212;
t1 = [(t153 ^ 2 + t154 ^ 2 + t192 ^ 2) * MDP(12) - 0.2e1 * t176 * t248 + MDP(1) + (MDP(4) * t215 + 0.2e1 * MDP(5) * t218) * t259 + (MDP(17) + MDP(24)) * t176 ^ 2 + (0.2e1 * t246 + t250) * t150 + (-0.2e1 * t141 * MDP(28) + t251) * t142 + (0.2e1 * t226 + t256) * t211 + (-0.2e1 * t166 * MDP(14) + 0.2e1 * t176 * MDP(15) + t247) * t167 + (-0.2e1 * t245 - 0.2e1 * t249 + t252 - 0.2e1 * t253 + 0.2e1 * t254) * t149 + 0.2e1 * (t182 * t211 + t205 * t270) * MDP(9) + 0.2e1 * (-pkin(1) * t259 - t183 * t211) * MDP(10) + (-t137 * t176 + t151 * t167) * t275 + 0.2e1 * (-t130 * t176 + t140 * t150) * MDP(26) + 0.2e1 * (-t153 * t177 - t154 * t176) * MDP(11) + 0.2e1 * (t136 * t176 + t151 * t166) * MDP(18) + (t129 * t176 + t140 * t149) * t274 + (t124 * t149 + t127 * t141) * t273 + (-t125 * t149 + t127 * t142) * t272; t182 * MDP(9) - t183 * MDP(10) + t256 + (t161 * t141 + t143 * t149) * MDP(32) + (t161 * t142 - t144 * t149) * MDP(33) + t138 * MDP(29) + (-t214 * t166 + t167 * t217) * MDP(14) + (-t151 * t217 + t198 * t166) * MDP(18) + (t151 * t214 + t198 * t167) * MDP(19) + t214 * t247 + t226 + (t149 * MDP(25) + t150 * MDP(26)) * t193 + (t220 + t277) * t176 + (t140 * MDP(25) + t221 - t245 - t249) * t190 + ((t153 * t210 + t154 * t208) * MDP(12) + (-t176 * t208 - t177 * t210) * MDP(11)) * pkin(2) + (t140 * MDP(26) + (-t141 * t216 - t142 * t212) * MDP(28) + t142 * t242 + t250 + t246 + t228 * t127 + t227 * t149) * t191; -0.2e1 * t198 * t240 + 0.2e1 * t193 * t187 + MDP(8) + (t208 ^ 2 + t210 ^ 2) * MDP(12) * pkin(2) ^ 2 + (0.2e1 * t227 * t191 + t193 * t274 + t244) * t190 + (t143 * t190 + t161 * t261) * t273 + (-t144 * t190 + t161 * t260) * t272 + (MDP(13) * t214 + 0.2e1 * t217 * MDP(14) + t198 * t275) * t214 + (t207 * MDP(27) + MDP(20) - 0.2e1 * t236) * t191 ^ 2; t192 * MDP(12) + (t190 * t141 - t149 * t261) * MDP(32) + (-t149 * t260 + t138) * MDP(33) + (t231 + t276) * t176; 0; MDP(12); t167 * MDP(15) - t248 + t176 * MDP(17) + t136 * MDP(18) - t137 * MDP(19) + (t176 * t238 + t129) * MDP(25) + (-t237 + (-t134 - t269) * t213) * MDP(26) + (t201 * t141 - t200 * t264 - t265) * MDP(32) + (t201 * t142 - t200 * t263 + t126) * MDP(33) + t224; (t232 * t212 - t262) * MDP(32) + (t232 * t216 + t157) * MDP(33) + t220 + t222; t225 + t231; -0.2e1 * t201 * t229 + MDP(17) + 0.2e1 * t223 + t235; t129 * MDP(25) - t130 * MDP(26) + (-pkin(5) * t141 - pkin(11) * t264 - t265) * MDP(32) + (-pkin(5) * t142 - pkin(11) * t263 + t126) * MDP(33) + t224; (t212 * t233 - t262) * MDP(32) + (t216 * t233 + t157) * MDP(33) + t222; t225; t223 + t235 + t229 * (pkin(5) - t201); 0.2e1 * pkin(5) * t229 + t235; t221 + t254; t143 * MDP(32) - t144 * MDP(33) + t191 * t230 + t244; -t228 * t191; -t200 * t228 + t255; -pkin(11) * t228 + t255; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
