% Calculate joint inertia matrix for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP11_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:44
% EndTime: 2019-03-10 02:53:50
% DurationCPUTime: 1.89s
% Computational Cost: add. (3768->361), mult. (9733->522), div. (0->0), fcn. (10953->12), ass. (0->153)
t221 = sin(pkin(7));
t310 = 0.2e1 * t221;
t225 = sin(qJ(5));
t229 = cos(qJ(5));
t233 = t225 * MDP(27) + t229 * MDP(28) - (t225 * MDP(30) + t229 * MDP(31)) * pkin(12) - MDP(21);
t224 = cos(pkin(6));
t228 = sin(qJ(2));
t222 = sin(pkin(6));
t232 = cos(qJ(2));
t281 = t222 * t232;
t201 = t224 * t228 * pkin(1) + pkin(9) * t281;
t223 = cos(pkin(7));
t280 = t223 * t232;
t251 = t222 * t280;
t177 = (t221 * t224 + t251) * pkin(10) + t201;
t227 = sin(qJ(3));
t231 = cos(qJ(3));
t301 = pkin(1) * t232;
t212 = t224 * t301;
t282 = t222 * t228;
t181 = pkin(2) * t224 + t212 + (-pkin(10) * t223 - pkin(9)) * t282;
t187 = (-pkin(10) * t221 * t228 - pkin(2) * t232 - pkin(1)) * t222;
t246 = t181 * t223 + t187 * t221;
t161 = -t227 * t177 + t246 * t231;
t309 = 2 * MDP(16);
t308 = 2 * MDP(17);
t307 = 2 * MDP(23);
t306 = 2 * MDP(24);
t305 = -2 * MDP(26);
t304 = 0.2e1 * MDP(30);
t303 = 0.2e1 * MDP(31);
t302 = 2 * MDP(32);
t300 = pkin(2) * t227;
t299 = pkin(2) * t231;
t298 = pkin(11) * t225;
t297 = pkin(11) * t229;
t230 = cos(qJ(4));
t296 = pkin(11) * t230;
t295 = -qJ(6) - pkin(12);
t294 = MDP(33) * pkin(5);
t293 = pkin(3) * MDP(23);
t292 = pkin(3) * MDP(24);
t291 = pkin(11) * MDP(24);
t226 = sin(qJ(4));
t290 = qJ(6) * t226;
t166 = -t181 * t221 + t223 * t187;
t283 = t221 * t231;
t179 = -t224 * t283 + t227 * t282 - t231 * t251;
t284 = t221 * t227;
t180 = t224 * t284 + (t227 * t280 + t228 * t231) * t222;
t155 = pkin(3) * t179 - pkin(11) * t180 + t166;
t162 = t177 * t231 + t246 * t227;
t194 = t221 * t281 - t224 * t223;
t158 = -pkin(11) * t194 + t162;
t149 = t155 * t230 - t226 * t158;
t147 = -pkin(4) * t179 - t149;
t289 = t147 * t225;
t288 = t147 * t229;
t253 = pkin(10) * t283;
t192 = t253 + (pkin(11) + t300) * t223;
t193 = (-pkin(3) * t231 - pkin(11) * t227 - pkin(2)) * t221;
t173 = -t226 * t192 + t193 * t230;
t171 = pkin(4) * t283 - t173;
t287 = t171 * t225;
t286 = t171 * t229;
t217 = t222 ^ 2;
t285 = t217 * t228;
t279 = t224 * MDP(8);
t278 = MDP(14) * t223;
t277 = MDP(22) * t231;
t276 = MDP(28) * t225;
t169 = t180 * t230 - t194 * t226;
t163 = t169 * t225 - t179 * t229;
t275 = t163 * MDP(28);
t164 = t169 * t229 + t179 * t225;
t274 = t164 * MDP(25);
t273 = t164 * MDP(27);
t168 = t180 * t226 + t194 * t230;
t272 = t168 * MDP(29);
t271 = t169 * MDP(19);
t270 = t169 * MDP(20);
t269 = t179 * MDP(22);
t268 = t180 * MDP(11);
t267 = t180 * MDP(12);
t197 = t223 * t226 + t230 * t284;
t182 = t197 * t225 + t229 * t283;
t266 = t182 * MDP(28);
t183 = t197 * t229 - t225 * t283;
t265 = t183 * MDP(25);
t264 = t183 * MDP(27);
t263 = t194 * MDP(13);
t262 = t194 * MDP(14);
t196 = -t230 * t223 + t226 * t284;
t261 = t196 * MDP(29);
t260 = t197 * MDP(18);
t259 = t197 * MDP(19);
t258 = t197 * MDP(20);
t257 = t223 * MDP(15);
t256 = t227 * MDP(13);
t255 = t229 * MDP(25);
t254 = t230 * MDP(29);
t252 = t229 * t296;
t250 = t225 * t229 * MDP(26);
t249 = pkin(11) * MDP(23) - MDP(20);
t248 = -MDP(21) + t291;
t247 = -(MDP(32) * pkin(5)) + MDP(27);
t150 = t155 * t226 + t158 * t230;
t148 = pkin(12) * t179 + t150;
t157 = pkin(3) * t194 - t161;
t152 = pkin(4) * t168 - pkin(12) * t169 + t157;
t144 = -t148 * t225 + t229 * t152;
t210 = pkin(10) * t284;
t191 = t210 + (-pkin(3) - t299) * t223;
t170 = pkin(4) * t196 - pkin(12) * t197 + t191;
t174 = t192 * t230 + t193 * t226;
t172 = -pkin(12) * t283 + t174;
t159 = t229 * t170 - t172 * t225;
t145 = t148 * t229 + t152 * t225;
t160 = t170 * t225 + t172 * t229;
t198 = t223 * t299 - t210;
t200 = t223 * t300 + t253;
t245 = -t198 * MDP(16) + t200 * MDP(17);
t205 = -pkin(4) * t230 - pkin(12) * t226 - pkin(3);
t202 = t229 * t205;
t188 = -t225 * t296 + t202;
t189 = t205 * t225 + t252;
t244 = t188 * MDP(30) - t189 * MDP(31);
t242 = MDP(27) * t229 - MDP(19) - t276;
t240 = (t228 * MDP(6) + t232 * MDP(7)) * t222;
t239 = t244 - t293;
t238 = t173 * MDP(23) - t174 * MDP(24) + t258;
t237 = t180 * MDP(13) - t194 * MDP(15) + t161 * MDP(16) - t162 * MDP(17);
t236 = t149 * MDP(23) - t150 * MDP(24) + t269 + t270;
t235 = t144 * MDP(30) - t145 * MDP(31) + t272 + t273 - t275;
t234 = t159 * MDP(30) - t160 * MDP(31) + t261 + t264 - t266;
t220 = t229 ^ 2;
t219 = t226 ^ 2;
t218 = t225 ^ 2;
t216 = t221 ^ 2;
t215 = -pkin(5) * t229 - pkin(4);
t207 = t295 * t229;
t206 = t295 * t225;
t204 = (pkin(5) * t225 + pkin(11)) * t226;
t199 = -pkin(9) * t282 + t212;
t184 = t252 + (t205 - t290) * t225;
t178 = -t229 * t290 + t202 + (-pkin(5) - t298) * t230;
t165 = pkin(5) * t182 + t171;
t154 = -qJ(6) * t182 + t160;
t153 = pkin(5) * t196 - qJ(6) * t183 + t159;
t146 = pkin(5) * t163 + t147;
t143 = -qJ(6) * t163 + t145;
t142 = pkin(5) * t168 - qJ(6) * t164 + t144;
t1 = [t169 ^ 2 * MDP(18) + (t142 ^ 2 + t143 ^ 2 + t146 ^ 2) * MDP(33) + t194 ^ 2 * MDP(15) + MDP(1) + (MDP(4) * t228 + 0.2e1 * MDP(5) * t232) * t285 + (-0.2e1 * t263 + t268) * t180 + (t163 * t305 + t274) * t164 + (0.2e1 * t240 + t279) * t224 + (0.2e1 * t262 - 0.2e1 * t267 + t269 + 0.2e1 * t270) * t179 + (-0.2e1 * t179 * MDP(21) - 0.2e1 * t271 + t272 + 0.2e1 * t273 - 0.2e1 * t275) * t168 + 0.2e1 * (t199 * t224 + t217 * t301) * MDP(9) + (-t145 * t168 + t147 * t164) * t303 + (-t142 * t164 - t143 * t163) * t302 + (t144 * t168 + t147 * t163) * t304 + (t162 * t194 + t166 * t180) * t308 + (-t161 * t194 + t166 * t179) * t309 + (t149 * t179 + t157 * t168) * t307 + (-t150 * t179 + t157 * t169) * t306 + 0.2e1 * (-pkin(1) * t285 - t201 * t224) * MDP(10); (t157 * t196 + t168 * t191) * MDP(23) + (t157 * t197 + t169 * t191) * MDP(24) + t279 - t201 * MDP(10) + (-t163 * t196 - t168 * t182) * MDP(28) + (t164 * t196 + t168 * t183) * MDP(27) + (t144 * t196 + t147 * t182 + t159 * t168 + t163 * t171) * MDP(30) + (-t145 * t196 + t147 * t183 - t160 * t168 + t164 * t171) * MDP(31) + (-t168 * t197 - t169 * t196) * MDP(19) + t199 * MDP(9) + (-t163 * t183 - t164 * t182) * MDP(26) + (-t142 * t183 - t143 * t182 - t153 * t164 - t154 * t163) * MDP(32) + (t142 * t153 + t143 * t154 + t146 * t165) * MDP(33) + t168 * t261 + t169 * t260 + t164 * t265 + t245 * t194 + t237 * t223 + t240 + (-t196 * MDP(21) + t238 - t278) * t179 + ((-t179 * MDP(16) - t180 * MDP(17)) * pkin(2) + (-t179 * MDP(12) + t166 * MDP(17) - t263 + t268) * t227 + (-t166 * MDP(16) + t168 * MDP(21) - t236 - t262 + t267) * t231) * t221; t197 ^ 2 * MDP(18) + (t153 ^ 2 + t154 ^ 2 + t165 ^ 2) * MDP(33) + t216 * t227 ^ 2 * MDP(11) + MDP(8) + (t256 * t310 + t257) * t223 + (t182 * t305 + t265) * t183 + ((-t258 + t278) * t310 + (0.2e1 * MDP(12) * t227 + t277) * t216) * t231 + (0.2e1 * MDP(21) * t283 - 0.2e1 * t259 + t261 + 0.2e1 * t264 - 0.2e1 * t266) * t196 + (t198 * t223 + t216 * t299) * t309 + (-t200 * t223 - t216 * t300) * t308 + (-t173 * t283 + t191 * t196) * t307 + (t174 * t283 + t191 * t197) * t306 + (t159 * t196 + t171 * t182) * t304 + (-t160 * t196 + t171 * t183) * t303 + (-t153 * t183 - t154 * t182) * t302; -t179 * MDP(14) - t169 * t292 + (-t163 * t184 - t164 * t178) * MDP(32) + (t142 * t178 + t143 * t184 + t146 * t204) * MDP(33) + t239 * t168 + (-t157 * MDP(23) - t248 * t179 - t235 + t271) * t230 + (t169 * MDP(18) + t157 * MDP(24) + t164 * t255 + (-t163 * t229 - t164 * t225) * MDP(26) + (pkin(11) * t163 + t289) * MDP(30) + (pkin(11) * t164 + t288) * MDP(31) + (-t142 * t229 - t143 * t225) * MDP(32) - t249 * t179 + t242 * t168) * t226 + t237; t257 - t197 * t292 + (-t178 * t183 - t182 * t184) * MDP(32) + (t153 * t178 + t154 * t184 + t165 * t204) * MDP(33) + (t231 * MDP(14) + t256) * t221 + t239 * t196 + (-t191 * MDP(23) + t248 * t283 - t234 + t259) * t230 + (t260 + t191 * MDP(24) + t183 * t255 + (-t182 * t229 - t183 * t225) * MDP(26) + (pkin(11) * t182 + t287) * MDP(30) + (pkin(11) * t183 + t286) * MDP(31) + (-t153 * t229 - t154 * t225) * MDP(32) + t249 * t283 + t242 * t196) * t226 - t245; MDP(15) + (t178 ^ 2 + t184 ^ 2 + t204 ^ 2) * MDP(33) + (t254 + (2 * t293)) * t230 + (MDP(25) * t220 + MDP(18) - 0.2e1 * t250) * t219 + (-t188 * t230 + t219 * t298) * t304 + (t189 * t230 + t219 * t297) * t303 + (-(2 * t292) + (-t178 * t229 - t184 * t225) * t302 - 0.2e1 * t242 * t230) * t226; t225 * t274 + (-t163 * t225 + t164 * t229) * MDP(26) + (-pkin(4) * t163 - t288) * MDP(30) + (-pkin(4) * t164 + t289) * MDP(31) + (-t142 * t225 + t143 * t229 + t163 * t207 - t164 * t206) * MDP(32) + (t142 * t206 - t143 * t207 + t146 * t215) * MDP(33) + t233 * t168 + t236; -t221 * t277 + t225 * t265 + (-t182 * t225 + t183 * t229) * MDP(26) + (-pkin(4) * t182 - t286) * MDP(30) + (-pkin(4) * t183 + t287) * MDP(31) + (-t153 * t225 + t154 * t229 + t182 * t207 - t183 * t206) * MDP(32) + (t153 * t206 - t154 * t207 + t165 * t215) * MDP(33) + t233 * t196 + t238; (-t178 * t225 + t184 * t229) * MDP(32) + (t178 * t206 - t184 * t207 + t204 * t215) * MDP(33) + (-t233 - t291) * t230 + (t225 * t255 + (-t218 + t220) * MDP(26) + (-pkin(4) * t225 - t297) * MDP(30) + (-pkin(4) * t229 + t298) * MDP(31) + (-t206 * t229 + t207 * t225) * MDP(32) - t249) * t226; MDP(22) + t218 * MDP(25) + 0.2e1 * t250 + (-t206 * t225 - t207 * t229) * t302 + (t206 ^ 2 + t207 ^ 2 + t215 ^ 2) * MDP(33) + 0.2e1 * (MDP(30) * t229 - MDP(31) * t225) * pkin(4); (-t164 * MDP(32) + t142 * MDP(33)) * pkin(5) + t235; (-t183 * MDP(32) + t153 * MDP(33)) * pkin(5) + t234; t178 * t294 - t254 + (t247 * t229 - t276) * t226 + t244; t206 * t294 + (-MDP(31) * pkin(12) + MDP(28)) * t229 + (-MDP(30) * pkin(12) + t247) * t225; MDP(33) * (pkin(5) ^ 2) + MDP(29); t146 * MDP(33); t165 * MDP(33); t204 * MDP(33); t215 * MDP(33); 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
