% Calculate joint inertia matrix for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR11_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:18
% EndTime: 2019-03-09 19:34:21
% DurationCPUTime: 1.40s
% Computational Cost: add. (1428->256), mult. (3126->357), div. (0->0), fcn. (3383->10), ass. (0->117)
t216 = cos(qJ(5));
t211 = sin(qJ(6));
t215 = cos(qJ(6));
t246 = t215 * MDP(34);
t232 = -t211 * MDP(35) + t246;
t228 = MDP(27) + t232;
t290 = t216 * t228;
t231 = t211 * MDP(34) + t215 * MDP(35);
t210 = cos(pkin(6));
t213 = sin(qJ(3));
t217 = cos(qJ(3));
t209 = sin(pkin(6));
t214 = sin(qJ(2));
t270 = t209 * t214;
t175 = -t210 * t217 + t213 * t270;
t176 = t210 * t213 + t217 * t270;
t289 = t176 * MDP(13) - t175 * MDP(14);
t288 = -2 * qJ(4) * MDP(20) - MDP(15);
t244 = MDP(17) - MDP(20);
t279 = pkin(9) - pkin(10);
t190 = t279 * t217;
t212 = sin(qJ(5));
t242 = t279 * t213;
t163 = t190 * t212 - t216 * t242;
t182 = -t212 * t217 + t213 * t216;
t204 = t211 ^ 2;
t206 = t215 ^ 2;
t164 = t216 * t190 + t212 * t242;
t255 = t164 * MDP(28);
t287 = -t228 * t163 - (t204 - t206) * t182 * MDP(30) - t255;
t218 = cos(qJ(2));
t269 = t209 * t218;
t194 = pkin(3) * t269;
t243 = pkin(8) * t269;
t278 = pkin(1) * t214;
t170 = t243 + (pkin(9) + t278) * t210;
t171 = (-pkin(2) * t218 - pkin(9) * t214 - pkin(1)) * t209;
t267 = t213 * t170 - t217 * t171;
t154 = t194 + t267;
t145 = pkin(4) * t269 - pkin(10) * t176 + t154;
t158 = t217 * t170 + t213 * t171;
t153 = -qJ(4) * t269 + t158;
t151 = pkin(10) * t175 + t153;
t140 = t145 * t216 - t212 * t151;
t141 = t212 * t145 + t216 * t151;
t161 = t175 * t212 + t176 * t216;
t155 = t161 * t211 - t215 * t269;
t156 = t161 * t215 + t211 * t269;
t256 = t161 * MDP(24);
t260 = t156 * MDP(29);
t286 = -t141 * MDP(28) + (-t155 * t211 + t156 * t215) * MDP(30) + t211 * t260 + t256 + t140 * MDP(27);
t285 = 2 * MDP(18);
t284 = 2 * MDP(19);
t283 = 0.2e1 * MDP(27);
t282 = 0.2e1 * MDP(28);
t281 = 0.2e1 * MDP(34);
t280 = 0.2e1 * MDP(35);
t277 = pkin(1) * t218;
t219 = -pkin(3) - pkin(4);
t186 = qJ(4) * t212 - t216 * t219;
t184 = pkin(5) + t186;
t276 = pkin(5) + t184;
t275 = pkin(3) * MDP(21);
t138 = -pkin(5) * t269 - t140;
t274 = t138 * t211;
t273 = t138 * t215;
t160 = -t175 * t216 + t176 * t212;
t272 = t160 * t212;
t268 = t214 * MDP(6);
t177 = -pkin(8) * t270 + t210 * t277;
t205 = t213 ^ 2;
t266 = t217 ^ 2 + t205;
t265 = MDP(28) * t212;
t181 = t212 * t213 + t216 * t217;
t264 = MDP(33) * t181;
t261 = t155 * MDP(32);
t259 = t156 * MDP(31);
t258 = t160 * MDP(33);
t257 = t161 * MDP(23);
t253 = t176 * MDP(11);
t251 = t186 * MDP(27);
t187 = qJ(4) * t216 + t212 * t219;
t250 = t187 * MDP(28);
t248 = t213 * MDP(13);
t247 = t215 * MDP(30);
t188 = -t217 * pkin(3) - t213 * qJ(4) - pkin(2);
t169 = -t210 * pkin(2) - t177;
t241 = t211 * t247;
t198 = t204 * MDP(29);
t240 = MDP(26) + 0.2e1 * t241 + t198;
t179 = t217 * pkin(4) - t188;
t239 = -MDP(18) + t265;
t238 = -t215 * t211 * MDP(29) - MDP(24);
t152 = t175 * pkin(3) - t176 * qJ(4) + t169;
t237 = t153 * t217 + t154 * t213;
t234 = MDP(31) * t215 - MDP(32) * t211;
t233 = t211 * MDP(31) + t215 * MDP(32);
t230 = -MDP(23) + t234;
t229 = -MDP(26) - t250 - t251;
t150 = -pkin(4) * t175 - t152;
t227 = 0.2e1 * t232;
t225 = -pkin(11) * t231 + t233;
t224 = -(-pkin(11) + t187) * t231 - t233;
t139 = pkin(11) * t269 + t141;
t142 = pkin(5) * t160 - pkin(11) * t161 + t150;
t136 = -t139 * t211 + t142 * t215;
t137 = t139 * t215 + t142 * t211;
t223 = t136 * MDP(34) - t137 * MDP(35) + t258 + t259 - t261;
t222 = -MDP(25) + t225;
t221 = MDP(25) + t224;
t203 = t209 ^ 2;
t189 = t213 * pkin(9) * t269;
t178 = t210 * t278 + t243;
t159 = pkin(5) * t181 - pkin(11) * t182 + t179;
t147 = t159 * t211 + t164 * t215;
t146 = t159 * t215 - t164 * t211;
t1 = [t161 ^ 2 * MDP(22) + (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) * MDP(21) + t210 ^ 2 * MDP(8) + MDP(1) + (-0.2e1 * t175 * MDP(12) + t253) * t176 + (-0.2e1 * t155 * MDP(30) + t260) * t156 + ((MDP(4) * t214 + 0.2e1 * MDP(5) * t218) * t214 + (MDP(15) + MDP(26)) * t218 ^ 2) * t203 + (-0.2e1 * t257 + t258 + 0.2e1 * t259 - 0.2e1 * t261) * t160 + 0.2e1 * (t210 * t268 + (-t160 * MDP(25) + MDP(7) * t210 + t256 - t289) * t218) * t209 + 0.2e1 * (t169 * t175 + t267 * t269) * MDP(16) + (t152 * t175 + t154 * t269) * t285 + (t140 * t269 + t150 * t160) * t283 + (-t141 * t269 + t150 * t161) * t282 + 0.2e1 * (t177 * t210 + t203 * t277) * MDP(9) + 0.2e1 * (t158 * t269 + t169 * t176) * MDP(17) + 0.2e1 * (-t152 * t176 - t153 * t269) * MDP(20) + 0.2e1 * (-t178 * t210 - t203 * t278) * MDP(10) + (-t153 * t175 + t154 * t176) * t284 + (t136 * t160 + t138 * t155) * t281 + (-t137 * t160 + t138 * t156) * t280; t213 * t253 + t177 * MDP(9) - t178 * MDP(10) + t210 * MDP(8) + (t146 * t160 + t155 * t163) * MDP(34) + (-t147 * t160 + t156 * t163) * MDP(35) + (-pkin(2) * t175 - t169 * t217 + t189) * MDP(16) + (-t152 * t217 + t175 * t188 + t189) * MDP(18) + t152 * t188 * MDP(21) + t237 * MDP(19) + (-t175 * t213 + t176 * t217) * MDP(12) + (-pkin(2) * t176 + t169 * t213) * MDP(17) + (-t152 * t213 - t176 * t188) * MDP(20) + (t160 * MDP(27) + t161 * MDP(28)) * t179 + (t237 * MDP(21) + (-t175 * t217 + t176 * t213) * MDP(19)) * pkin(9) + (t150 * MDP(27) + t223 - t257) * t181 + (t161 * MDP(22) + t215 * t260 + (-t155 * t215 - t156 * t211) * MDP(30) + t150 * MDP(28) + t231 * t138 + t230 * t160) * t182 + (t268 + (-t248 + t182 * MDP(24) - t181 * MDP(25) - t163 * MDP(27) - t255 + MDP(7) + (t244 * pkin(9) - MDP(14)) * t217) * t218) * t209; MDP(8) + t205 * MDP(11) + (t266 * pkin(9) ^ 2 + t188 ^ 2) * MDP(21) + t266 * pkin(9) * t284 + 0.2e1 * (pkin(2) * MDP(16) - t188 * MDP(18)) * t217 + 0.2e1 * (t217 * MDP(12) - pkin(2) * MDP(17) - t188 * MDP(20)) * t213 + (t146 * t281 - t147 * t280 + t179 * t283 + t264) * t181 + (t179 * t282 + (t211 * t281 + t215 * t280) * t163 + 0.2e1 * t230 * t181 + (t206 * MDP(29) + MDP(22) - 0.2e1 * t241) * t182) * t182; -t267 * MDP(16) + (-0.2e1 * t194 - t267) * MDP(18) + (-pkin(3) * t176 - qJ(4) * t175) * MDP(19) + (-pkin(3) * t154 + qJ(4) * t153) * MDP(21) + (t155 * t184 + t273) * MDP(34) + (t156 * t184 - t274) * MDP(35) + (t229 + t288) * t269 + t221 * t160 - t244 * t158 - t286 + t289; t248 + t217 * MDP(14) + (-pkin(3) * t213 + qJ(4) * t217) * MDP(19) + (t184 * t231 + t238) * t182 + t221 * t181 + ((MDP(21) * qJ(4) - t244) * t217 + (-MDP(16) - MDP(18) - t275) * t213) * pkin(9) - t287; pkin(3) * t285 + (pkin(3) ^ 2 + (qJ(4) ^ 2)) * MDP(21) + t184 * t227 + 0.2e1 * t251 + 0.2e1 * t250 + t240 - t288; t176 * MDP(19) + t154 * MDP(21) + (-t155 * t216 - t211 * t272) * MDP(34) + (-t156 * t216 - t215 * t272) * MDP(35) + (MDP(27) * t216 - t239) * t269; (MDP(21) * pkin(9) + MDP(19)) * t213 + t231 * (-t181 * t212 - t182 * t216); t239 - t275 - t290; MDP(21); MDP(26) * t269 + (-pkin(5) * t155 - t273) * MDP(34) + (-pkin(5) * t156 + t274) * MDP(35) + t222 * t160 + t286; (-pkin(5) * t231 - t238) * t182 + t222 * t181 + t287; -t198 - t276 * t246 + (t276 * MDP(35) - 0.2e1 * t247) * t211 + t229; -t265 + t290; pkin(5) * t227 + t240; t223; t146 * MDP(34) - t147 * MDP(35) + t182 * t234 + t264; t224; -t231 * t212; t225; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
