% Calculate joint inertia matrix for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR6_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:09
% EndTime: 2019-03-09 01:25:14
% DurationCPUTime: 1.67s
% Computational Cost: add. (2476->325), mult. (6739->489), div. (0->0), fcn. (7872->16), ass. (0->152)
t212 = sin(pkin(8));
t303 = 0.2e1 * t212;
t218 = sin(qJ(6));
t223 = cos(qJ(6));
t231 = -(t218 * MDP(31) + t223 * MDP(32)) * pkin(13) + t218 * MDP(28) + t223 * MDP(29);
t216 = cos(pkin(7));
t221 = sin(qJ(3));
t213 = sin(pkin(7));
t226 = cos(qJ(3));
t278 = t213 * t226;
t199 = pkin(2) * t216 * t221 + pkin(10) * t278;
t215 = cos(pkin(8));
t276 = t215 * t226;
t248 = t213 * t276;
t174 = (t212 * t216 + t248) * pkin(11) + t199;
t220 = sin(qJ(4));
t225 = cos(qJ(4));
t295 = pkin(2) * t226;
t204 = t216 * t295;
t279 = t213 * t221;
t179 = pkin(3) * t216 + t204 + (-pkin(11) * t215 - pkin(10)) * t279;
t184 = (-pkin(11) * t212 * t221 - pkin(3) * t226 - pkin(2)) * t213;
t244 = t179 * t215 + t184 * t212;
t155 = -t220 * t174 + t244 * t225;
t302 = 2 * MDP(17);
t301 = 2 * MDP(18);
t300 = 2 * MDP(24);
t299 = 2 * MDP(25);
t298 = -2 * MDP(27);
t297 = 0.2e1 * MDP(31);
t296 = 0.2e1 * MDP(32);
t294 = pkin(3) * t220;
t293 = pkin(3) * t225;
t292 = pkin(12) * t218;
t291 = pkin(12) * t223;
t224 = cos(qJ(5));
t290 = pkin(12) * t224;
t289 = pkin(4) * MDP(24);
t288 = pkin(4) * MDP(25);
t161 = -t179 * t212 + t215 * t184;
t280 = t212 * t225;
t175 = -t216 * t280 + t220 * t279 - t225 * t248;
t281 = t212 * t220;
t178 = t216 * t281 + (t220 * t276 + t221 * t225) * t213;
t150 = pkin(4) * t175 - pkin(12) * t178 + t161;
t156 = t174 * t225 + t244 * t220;
t191 = t212 * t278 - t215 * t216;
t152 = -pkin(12) * t191 + t156;
t219 = sin(qJ(5));
t145 = t150 * t224 - t152 * t219;
t141 = -pkin(5) * t175 - t145;
t287 = t141 * t218;
t286 = t141 * t223;
t214 = sin(pkin(6));
t217 = cos(pkin(6));
t222 = sin(qJ(2));
t227 = cos(qJ(2));
t274 = t216 * t227;
t176 = t217 * t278 + (-t221 * t222 + t226 * t274) * t214;
t177 = t217 * t279 + (t221 * t274 + t222 * t226) * t214;
t193 = -t213 * t214 * t227 + t217 * t216;
t158 = t177 * t225 + (t176 * t215 + t193 * t212) * t220;
t163 = -t176 * t212 + t193 * t215;
t148 = t158 * t219 - t163 * t224;
t285 = t148 * t219;
t249 = pkin(11) * t280;
t189 = t249 + (pkin(12) + t294) * t215;
t190 = (-pkin(4) * t225 - pkin(12) * t220 - pkin(3)) * t212;
t169 = -t189 * t219 + t190 * t224;
t167 = pkin(5) * t280 - t169;
t284 = t167 * t218;
t283 = t167 * t223;
t208 = t213 ^ 2;
t282 = t208 * t221;
t277 = t215 * t225;
t275 = t216 * MDP(9);
t273 = t218 * t223;
t272 = MDP(15) * t215;
t271 = MDP(23) * t225;
t165 = t178 * t224 - t191 * t219;
t159 = t165 * t218 - t223 * t175;
t270 = t159 * MDP(29);
t160 = t165 * t223 + t175 * t218;
t269 = t160 * MDP(26);
t268 = t160 * MDP(28);
t164 = t178 * t219 + t224 * t191;
t267 = t164 * MDP(30);
t266 = t165 * MDP(20);
t265 = t165 * MDP(21);
t264 = t175 * MDP(23);
t263 = t178 * MDP(12);
t262 = t178 * MDP(13);
t195 = t215 * t219 + t224 * t281;
t180 = t195 * t218 + t223 * t280;
t261 = t180 * MDP(29);
t181 = t195 * t223 - t218 * t280;
t260 = t181 * MDP(26);
t259 = t181 * MDP(28);
t258 = t191 * MDP(14);
t257 = t191 * MDP(15);
t194 = -t224 * t215 + t219 * t281;
t256 = t194 * MDP(30);
t255 = t195 * MDP(19);
t254 = t195 * MDP(20);
t253 = t195 * MDP(21);
t252 = t215 * MDP(16);
t251 = t220 * MDP(14);
t250 = t224 * MDP(30);
t247 = MDP(27) * t273;
t246 = pkin(12) * MDP(24) - MDP(21);
t245 = pkin(12) * MDP(25) - MDP(22);
t146 = t150 * t219 + t152 * t224;
t170 = t189 * t224 + t190 * t219;
t203 = pkin(11) * t281;
t196 = pkin(3) * t277 - t203;
t198 = t215 * t294 + t249;
t243 = -t196 * MDP(17) + t198 * MDP(18);
t242 = t223 * MDP(28) - t218 * MDP(29);
t201 = -pkin(5) * t224 - pkin(13) * t219 - pkin(4);
t185 = t201 * t223 - t218 * t290;
t186 = t201 * t218 + t223 * t290;
t240 = t185 * MDP(31) - t186 * MDP(32);
t239 = MDP(31) * t223 - MDP(32) * t218;
t188 = t203 + (-pkin(4) - t293) * t215;
t237 = -MDP(20) + t242;
t236 = (t221 * MDP(7) + t226 * MDP(8)) * t213;
t235 = t240 - t289;
t234 = t169 * MDP(24) - t170 * MDP(25) + t253;
t151 = pkin(4) * t191 - t155;
t233 = t178 * MDP(14) - t191 * MDP(16) + t155 * MDP(17) - t156 * MDP(18);
t232 = t145 * MDP(24) - t146 * MDP(25) + t264 + t265;
t142 = pkin(13) * t175 + t146;
t147 = pkin(5) * t164 - pkin(13) * t165 + t151;
t139 = -t142 * t218 + t147 * t223;
t140 = t142 * t223 + t147 * t218;
t230 = t139 * MDP(31) - t140 * MDP(32) + t267 + t268 - t270;
t166 = pkin(5) * t194 - pkin(13) * t195 + t188;
t168 = -pkin(13) * t280 + t170;
t153 = t166 * t223 - t168 * t218;
t154 = t166 * t218 + t168 * t223;
t229 = t153 * MDP(31) - t154 * MDP(32) + t256 + t259 - t261;
t228 = -MDP(22) + t231;
t211 = t223 ^ 2;
t210 = t219 ^ 2;
t209 = t218 ^ 2;
t207 = t212 ^ 2;
t197 = -pkin(10) * t279 + t204;
t157 = -t176 * t277 + t177 * t220 - t193 * t280;
t149 = t158 * t224 + t163 * t219;
t144 = t149 * t223 + t157 * t218;
t143 = -t149 * t218 + t157 * t223;
t1 = [MDP(1); (t176 * t216 - t193 * t278) * MDP(10) + (-t177 * t216 + t193 * t279) * MDP(11) + (t157 * t191 + t163 * t175) * MDP(17) + (t158 * t191 + t163 * t178) * MDP(18) + (-t148 * t175 + t157 * t164) * MDP(24) + (-t149 * t175 + t157 * t165) * MDP(25) + (t143 * t164 + t148 * t159) * MDP(31) + (-t144 * t164 + t148 * t160) * MDP(32) + (MDP(3) * t227 - MDP(4) * t222) * t214; t191 ^ 2 * MDP(16) + t165 ^ 2 * MDP(19) + MDP(2) + (MDP(5) * t221 + 0.2e1 * MDP(6) * t226) * t282 + (-0.2e1 * t258 + t263) * t178 + (t159 * t298 + t269) * t160 + (0.2e1 * t236 + t275) * t216 + (0.2e1 * t257 - 0.2e1 * t262 + t264 + 0.2e1 * t265) * t175 + (-0.2e1 * t175 * MDP(22) - 0.2e1 * t266 + t267 + 0.2e1 * t268 - 0.2e1 * t270) * t164 + 0.2e1 * (t197 * t216 + t208 * t295) * MDP(10) + 0.2e1 * (-pkin(2) * t282 - t199 * t216) * MDP(11) + (t156 * t191 + t161 * t178) * t301 + (-t155 * t191 + t161 * t175) * t302 + (-t146 * t175 + t151 * t165) * t299 + (t145 * t175 + t151 * t164) * t300 + (t139 * t164 + t141 * t159) * t297 + (-t140 * t164 + t141 * t160) * t296; t176 * MDP(10) - t177 * MDP(11) + (-t157 * t215 - t163 * t280) * MDP(17) + (-t158 * t215 + t163 * t281) * MDP(18) + (t148 * t280 + t157 * t194) * MDP(24) + (t149 * t280 + t157 * t195) * MDP(25) + (t143 * t194 + t148 * t180) * MDP(31) + (-t144 * t194 + t148 * t181) * MDP(32); (t151 * t195 + t165 * t188) * MDP(25) + (t151 * t194 + t164 * t188) * MDP(24) + t275 + (-t140 * t194 + t141 * t181 - t154 * t164 + t160 * t167) * MDP(32) + (-t164 * t195 - t165 * t194) * MDP(20) + t197 * MDP(10) - t199 * MDP(11) + (t160 * t194 + t164 * t181) * MDP(28) + (t139 * t194 + t141 * t180 + t153 * t164 + t159 * t167) * MDP(31) + (-t159 * t194 - t164 * t180) * MDP(29) + (-t159 * t181 - t160 * t180) * MDP(27) + t164 * t256 + t165 * t255 + t160 * t260 + t243 * t191 + t233 * t215 + t236 + (-t194 * MDP(22) + t234 - t272) * t175 + ((-t175 * MDP(17) - t178 * MDP(18)) * pkin(3) + (-t175 * MDP(13) + t161 * MDP(18) - t258 + t263) * t220 + (-t161 * MDP(17) + t164 * MDP(22) - t232 - t257 + t262) * t225) * t212; t207 * t220 ^ 2 * MDP(12) + t195 ^ 2 * MDP(19) + MDP(9) + (t251 * t303 + t252) * t215 + (t180 * t298 + t260) * t181 + ((-t253 + t272) * t303 + (0.2e1 * MDP(13) * t220 + t271) * t207) * t225 + (0.2e1 * MDP(22) * t280 - 0.2e1 * t254 + t256 + 0.2e1 * t259 - 0.2e1 * t261) * t194 + (t196 * t215 + t207 * t293) * t302 + (-t198 * t215 - t207 * t294) * t301 + (-t169 * t280 + t188 * t194) * t300 + (t170 * t280 + t188 * t195) * t299 + (t153 * t194 + t167 * t180) * t297 + (-t154 * t194 + t167 * t181) * t296; -t158 * MDP(18) + (-t143 * t224 + t218 * t285) * MDP(31) + (t144 * t224 + t223 * t285) * MDP(32) + (-MDP(24) * t224 + MDP(25) * t219 - MDP(17)) * t157; -t165 * t288 - t175 * MDP(15) + t235 * t164 + (-t151 * MDP(24) - t245 * t175 - t230 + t266) * t224 + (t165 * MDP(19) + t151 * MDP(25) + t223 * t269 + (-t159 * t223 - t160 * t218) * MDP(27) + (pkin(12) * t159 + t287) * MDP(31) + (pkin(12) * t160 + t286) * MDP(32) - t246 * t175 + t237 * t164) * t219 + t233; -t195 * t288 + t252 + (t225 * MDP(15) + t251) * t212 + t235 * t194 + (-t188 * MDP(24) + t245 * t280 - t229 + t254) * t224 + (t255 + t188 * MDP(25) + t223 * t260 + (-t180 * t223 - t181 * t218) * MDP(27) + (pkin(12) * t180 + t284) * MDP(31) + (pkin(12) * t181 + t283) * MDP(32) + t246 * t280 + t237 * t194) * t219 - t243; MDP(16) + (t250 + (2 * t289)) * t224 + (MDP(26) * t211 + MDP(19) - 0.2e1 * t247) * t210 + (-t185 * t224 + t210 * t292) * t297 + (t186 * t224 + t210 * t291) * t296 + 0.2e1 * (-t224 * t237 - t288) * t219; -MDP(25) * t149 + (-MDP(24) - t239) * t148; t218 * t269 + (-t159 * t218 + t160 * t223) * MDP(27) + (-pkin(5) * t159 - t286) * MDP(31) + (-pkin(5) * t160 + t287) * MDP(32) + t228 * t164 + t232; -t212 * t271 + t218 * t260 + (-t180 * t218 + t181 * t223) * MDP(27) + (-pkin(5) * t180 - t283) * MDP(31) + (-pkin(5) * t181 + t284) * MDP(32) + t228 * t194 + t234; (-t231 - t245) * t224 + (MDP(26) * t273 + (-t209 + t211) * MDP(27) + (-pkin(5) * t218 - t291) * MDP(31) + (-pkin(5) * t223 + t292) * MDP(32) - t246) * t219; MDP(26) * t209 + 0.2e1 * pkin(5) * t239 + MDP(23) + 0.2e1 * t247; MDP(31) * t143 - MDP(32) * t144; t230; t229; t242 * t219 + t240 - t250; t231; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
