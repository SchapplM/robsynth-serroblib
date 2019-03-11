% Calculate joint inertia matrix for
% S6RRPRRR5
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
%   see S6RRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR5_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:51
% EndTime: 2019-03-09 13:46:56
% DurationCPUTime: 1.61s
% Computational Cost: add. (2011->247), mult. (4850->357), div. (0->0), fcn. (5504->12), ass. (0->121)
t219 = sin(qJ(5));
t223 = cos(qJ(5));
t236 = MDP(25) * t219 + MDP(26) * t223;
t218 = sin(qJ(6));
t222 = cos(qJ(6));
t197 = t218 * t219 - t222 * t223;
t198 = t218 * t223 + t219 * t222;
t278 = pkin(10) + pkin(11);
t201 = t278 * t219;
t202 = t278 * t223;
t242 = t198 * MDP(29) - t197 * MDP(30) + (-t201 * t222 - t202 * t218) * MDP(32) - (-t201 * t218 + t202 * t222) * MDP(33);
t290 = t219 * MDP(22) + t223 * MDP(23) - t236 * pkin(10) + t242;
t216 = cos(pkin(12));
t206 = -pkin(2) * t216 - pkin(3);
t220 = sin(qJ(4));
t224 = cos(qJ(4));
t196 = -pkin(4) * t224 - pkin(10) * t220 + t206;
t190 = t223 * t196;
t264 = t220 * t223;
t214 = sin(pkin(12));
t205 = pkin(2) * t214 + pkin(9);
t271 = t205 * t219;
t164 = -pkin(11) * t264 + t190 + (-pkin(5) - t271) * t224;
t269 = t205 * t224;
t246 = t223 * t269;
t165 = t246 + (-pkin(11) * t220 + t196) * t219;
t151 = t222 * t164 - t165 * t218;
t152 = t164 * t218 + t165 * t222;
t189 = t197 * t220;
t184 = t189 * MDP(29);
t289 = MDP(32) * t151 - MDP(33) * t152 - t184;
t288 = -MDP(16) + t290;
t188 = t198 * t220;
t183 = t188 * MDP(30);
t287 = -t183 + t289;
t230 = (MDP(32) * t222 - MDP(33) * t218) * pkin(5);
t215 = sin(pkin(6));
t221 = sin(qJ(2));
t225 = cos(qJ(2));
t187 = (t214 * t225 + t216 * t221) * t215;
t217 = cos(pkin(6));
t172 = t187 * t224 + t217 * t220;
t266 = t215 * t225;
t267 = t215 * t221;
t186 = t214 * t267 - t216 * t266;
t159 = t172 * t219 - t186 * t223;
t160 = t172 * t223 + t186 * t219;
t149 = t222 * t159 + t160 * t218;
t150 = -t159 * t218 + t160 * t222;
t285 = t150 * MDP(29) - t149 * MDP(30);
t277 = pkin(1) * t225;
t203 = t217 * t277;
t275 = pkin(8) + qJ(3);
t175 = pkin(2) * t217 - t275 * t267 + t203;
t247 = pkin(1) * t217 * t221;
t180 = t275 * t266 + t247;
t158 = t214 * t175 + t216 * t180;
t155 = pkin(9) * t217 + t158;
t199 = (-pkin(2) * t225 - pkin(1)) * t215;
t161 = pkin(3) * t186 - pkin(9) * t187 + t199;
t146 = t155 * t224 + t161 * t220;
t142 = pkin(10) * t186 + t146;
t157 = t175 * t216 - t214 * t180;
t154 = -pkin(3) * t217 - t157;
t171 = t187 * t220 - t217 * t224;
t144 = pkin(4) * t171 - pkin(10) * t172 + t154;
t138 = -t142 * t219 + t223 * t144;
t139 = t142 * t223 + t144 * t219;
t284 = -t138 * MDP(25) + t139 * MDP(26);
t283 = 0.2e1 * MDP(25);
t282 = 0.2e1 * MDP(26);
t281 = -2 * MDP(28);
t280 = 0.2e1 * MDP(32);
t279 = 0.2e1 * MDP(33);
t276 = pkin(5) * t171;
t137 = -pkin(11) * t159 + t139;
t274 = t137 * t222;
t145 = -t220 * t155 + t161 * t224;
t141 = -pkin(4) * t186 - t145;
t273 = t141 * t219;
t272 = t141 * t223;
t270 = t205 * t223;
t209 = t215 ^ 2;
t268 = t209 * t221;
t265 = t217 * MDP(8);
t263 = -t188 * MDP(32) + t189 * MDP(33);
t261 = MDP(19) * t220;
t260 = MDP(20) * t160;
t259 = MDP(20) * t223;
t258 = MDP(27) * t189;
t257 = MDP(27) * t198;
t255 = MDP(32) * t197;
t251 = t172 * MDP(13);
t250 = t205 * MDP(19);
t249 = t206 * MDP(18);
t248 = MDP(24) + MDP(31);
t245 = t171 * MDP(31) + t285;
t244 = MDP(21) * t219 * t223;
t136 = -pkin(11) * t160 + t138 + t276;
t133 = t222 * t136 - t137 * t218;
t243 = -t205 * MDP(18) + MDP(15);
t241 = t160 * MDP(22) - t159 * MDP(23);
t240 = MDP(22) * t223 - MDP(23) * t219;
t168 = -t219 * t269 + t190;
t169 = t196 * t219 + t246;
t238 = t168 * MDP(25) - t169 * MDP(26);
t237 = MDP(25) * t223 - MDP(26) * t219;
t134 = t136 * t218 + t274;
t235 = -t133 * MDP(32) + t134 * MDP(33);
t232 = -MDP(14) + t240;
t229 = (MDP(6) * t221 + MDP(7) * t225) * t215;
t228 = t172 * MDP(14) - t241 - t285;
t212 = t223 ^ 2;
t210 = t219 ^ 2;
t208 = -pkin(5) * t223 - pkin(4);
t193 = pkin(8) * t266 + t247;
t192 = -pkin(8) * t267 + t203;
t191 = (pkin(5) * t219 + t205) * t220;
t162 = t188 * t171;
t140 = pkin(5) * t159 + t141;
t1 = [MDP(1) + t186 ^ 2 * MDP(17) + (t157 ^ 2 + t158 ^ 2 + t199 ^ 2) * MDP(12) + (MDP(4) * t221 + 0.2e1 * MDP(5) * t225) * t268 + (0.2e1 * t186 * MDP(15) + t251) * t172 + t248 * t171 ^ 2 + (-0.2e1 * t159 * MDP(21) + t260) * t160 + (MDP(27) * t150 + t149 * t281) * t150 + t265 * t217 + 0.2e1 * (t192 * t217 + t209 * t277) * MDP(9) + 0.2e1 * (-pkin(1) * t268 - t193 * t217) * MDP(10) + (-t134 * t171 + t140 * t150) * t279 + 0.2e1 * (t145 * t186 + t154 * t171) * MDP(18) + 0.2e1 * (-t146 * t186 + t154 * t172) * MDP(19) + 0.2e1 * (-t157 * t187 - t158 * t186) * MDP(11) + (t138 * t171 + t141 * t159) * t283 + (t133 * t171 + t140 * t149) * t280 + (-t139 * t171 + t141 * t160) * t282 + 0.2e1 * t229 * t217 + 0.2e1 * (-t186 * MDP(16) - t228) * t171; -t150 * t258 + (t149 * t189 - t150 * t188) * MDP(28) + t192 * MDP(9) - t193 * MDP(10) + t265 - t162 * MDP(30) + (t140 * t188 + t149 * t191) * MDP(32) + (-t140 * t189 + t150 * t191) * MDP(33) + t206 * t172 * MDP(19) + t229 + (t238 + t249 + t289) * t171 + ((t157 * t216 + t158 * t214) * MDP(12) + (-t186 * t214 - t187 * t216) * MDP(11)) * pkin(2) + (-t154 * MDP(18) + (MDP(16) - t250) * t186 - t248 * t171 + t228 + t235 + t284) * t224 + (t160 * t259 + t251 + (-t159 * t223 - t160 * t219) * MDP(21) + t154 * MDP(19) + (t159 * t205 + t273) * MDP(25) + (t160 * t205 + t272) * MDP(26) + t243 * t186 + t232 * t171) * t220; 0.2e1 * t206 * t261 + MDP(8) + (t214 ^ 2 + t216 ^ 2) * MDP(12) * pkin(2) ^ 2 + t188 * t191 * t280 - (t188 * t281 + t191 * t279 - t258) * t189 + (MDP(20) * t212 + t270 * t282 + t271 * t283 + MDP(13) - 0.2e1 * t244) * t220 ^ 2 + (-t151 * t280 + t152 * t279 - t168 * t283 + t169 * t282 - 0.2e1 * t232 * t220 + t248 * t224 + 0.2e1 * t183 + 0.2e1 * t184 - 0.2e1 * t249) * t224; t199 * MDP(12) + (-t171 * t219 * t220 - t159 * t224) * MDP(25) + (-t160 * t224 - t171 * t264) * MDP(26) + (-t149 * t224 - t162) * MDP(32) + (-t150 * t224 + t171 * t189) * MDP(33) + (MDP(18) * t224 - t261) * t186; 0; MDP(12); t172 * MDP(15) + t186 * MDP(17) + t145 * MDP(18) - t146 * MDP(19) + t219 * t260 + (-t159 * t219 + t160 * t223) * MDP(21) + (-pkin(4) * t159 - t272) * MDP(25) + (-pkin(4) * t160 + t273) * MDP(26) + t150 * t257 + (-t149 * t198 - t150 * t197) * MDP(28) + (t140 * t197 + t149 * t208) * MDP(32) + (t140 * t198 + t150 * t208) * MDP(33) + t288 * t171; -t189 * t257 + (-t188 * t198 + t189 * t197) * MDP(28) + (t188 * t208 + t191 * t197) * MDP(32) + (-t189 * t208 + t191 * t198) * MDP(33) + (-t250 - t288) * t224 + (t219 * t259 + (-t210 + t212) * MDP(21) + (-pkin(4) * t219 - t270) * MDP(25) + (-pkin(4) * t223 + t271) * MDP(26) + t243) * t220; -t261 + (-MDP(33) * t198 + MDP(18) + t237 - t255) * t224; 0.2e1 * t244 + 0.2e1 * t208 * t255 + MDP(20) * t210 + MDP(17) + 0.2e1 * t237 * pkin(4) + (t197 * t281 + t208 * t279 + t257) * t198; t171 * MDP(24) + (t222 * t276 + t133) * MDP(32) + (-t274 + (-t136 - t276) * t218) * MDP(33) + t241 + t245 - t284; (-t248 - t230) * t224 + t240 * t220 + t238 + t287; -t236 * t220 + t263; t290; 0.2e1 * t230 + t248; -t235 + t245; -t224 * MDP(31) + t287; t263; t242; MDP(31) + t230; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
