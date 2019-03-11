% Calculate joint inertia matrix for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR12_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:47
% EndTime: 2019-03-09 23:44:54
% DurationCPUTime: 2.27s
% Computational Cost: add. (4825->348), mult. (12411->527), div. (0->0), fcn. (14273->14), ass. (0->151)
t252 = sin(pkin(7));
t332 = 0.2e1 * t252;
t258 = sin(qJ(4));
t262 = cos(qJ(4));
t270 = -(MDP(23) * t258 + MDP(24) * t262) * pkin(11) + t258 * MDP(20) + t262 * MDP(21);
t331 = -MDP(14) + t270;
t256 = cos(pkin(6));
t260 = sin(qJ(2));
t253 = sin(pkin(6));
t264 = cos(qJ(2));
t313 = t253 * t264;
t231 = pkin(1) * t256 * t260 + pkin(9) * t313;
t255 = cos(pkin(7));
t312 = t255 * t264;
t284 = t253 * t312;
t209 = (t252 * t256 + t284) * pkin(10) + t231;
t259 = sin(qJ(3));
t263 = cos(qJ(3));
t322 = pkin(1) * t264;
t240 = t256 * t322;
t314 = t253 * t260;
t212 = pkin(2) * t256 + t240 + (-pkin(10) * t255 - pkin(9)) * t314;
t218 = (-pkin(10) * t252 * t260 - pkin(2) * t264 - pkin(1)) * t253;
t281 = t212 * t255 + t218 * t252;
t187 = -t259 * t209 + t281 * t263;
t330 = 2 * MDP(16);
t329 = 2 * MDP(17);
t328 = -2 * MDP(19);
t327 = 0.2e1 * MDP(23);
t326 = 0.2e1 * MDP(24);
t325 = 2 * MDP(25);
t324 = 2 * MDP(32);
t323 = 2 * MDP(33);
t321 = pkin(2) * t259;
t320 = pkin(2) * t263;
t319 = -qJ(5) - pkin(11);
t237 = t319 * t262;
t251 = sin(pkin(13));
t254 = cos(pkin(13));
t282 = t319 * t258;
t213 = -t237 * t251 - t254 * t282;
t234 = t251 * t262 + t254 * t258;
t318 = t213 * t234;
t248 = t253 ^ 2;
t317 = t248 * t260;
t316 = t252 * t259;
t315 = t252 * t263;
t311 = t256 * MDP(8);
t257 = sin(qJ(6));
t261 = cos(qJ(6));
t310 = t257 * t261;
t195 = -t212 * t252 + t255 * t218;
t210 = -t256 * t315 + t259 * t314 - t263 * t284;
t211 = t256 * t316 + (t259 * t312 + t260 * t263) * t253;
t181 = pkin(3) * t210 - pkin(11) * t211 + t195;
t188 = t209 * t263 + t281 * t259;
t224 = t252 * t313 - t256 * t255;
t184 = -pkin(11) * t224 + t188;
t171 = t262 * t181 - t184 * t258;
t198 = t211 * t262 - t224 * t258;
t166 = pkin(4) * t210 - qJ(5) * t198 + t171;
t172 = t181 * t258 + t184 * t262;
t197 = t211 * t258 + t224 * t262;
t168 = -qJ(5) * t197 + t172;
t163 = t251 * t166 + t254 * t168;
t285 = pkin(10) * t315;
t221 = t285 + (pkin(11) + t321) * t255;
t222 = (-pkin(3) * t263 - pkin(11) * t259 - pkin(2)) * t252;
t201 = -t221 * t258 + t262 * t222;
t227 = t255 * t258 + t262 * t316;
t190 = -pkin(4) * t315 - qJ(5) * t227 + t201;
t202 = t221 * t262 + t222 * t258;
t226 = -t262 * t255 + t258 * t316;
t194 = -qJ(5) * t226 + t202;
t179 = t251 * t190 + t254 * t194;
t186 = -t197 * t251 + t198 * t254;
t174 = t186 * t257 - t210 * t261;
t309 = t174 * MDP(28);
t308 = t174 * MDP(30);
t175 = t186 * t261 + t210 * t257;
t307 = t175 * MDP(27);
t185 = t254 * t197 + t198 * t251;
t306 = t185 * MDP(29);
t305 = t185 * MDP(31);
t304 = t197 * MDP(21);
t303 = t198 * MDP(18);
t302 = t198 * MDP(20);
t204 = -t226 * t251 + t227 * t254;
t199 = t204 * t257 + t261 * t315;
t301 = t199 * MDP(28);
t300 = t199 * MDP(30);
t200 = t204 * t261 - t257 * t315;
t299 = t200 * MDP(27);
t203 = t254 * t226 + t227 * t251;
t298 = t203 * MDP(29);
t297 = t203 * MDP(31);
t296 = t210 * MDP(22);
t295 = t211 * MDP(11);
t294 = t211 * MDP(12);
t293 = t224 * MDP(13);
t292 = t224 * MDP(14);
t291 = t227 * MDP(18);
t233 = t251 * t258 - t254 * t262;
t290 = t233 * MDP(31);
t289 = t255 * MDP(14);
t288 = t255 * MDP(15);
t287 = t259 * MDP(13);
t286 = t263 * MDP(22);
t246 = -pkin(4) * t262 - pkin(3);
t283 = MDP(28) * t310;
t162 = t166 * t254 - t168 * t251;
t178 = t190 * t254 - t194 * t251;
t239 = pkin(10) * t316;
t228 = t255 * t320 - t239;
t230 = t255 * t321 + t285;
t280 = -t228 * MDP(16) + t230 * MDP(17);
t279 = t227 * MDP(20) - t226 * MDP(21);
t276 = MDP(32) * t261 - t257 * MDP(33);
t275 = MDP(32) * t257 + MDP(33) * t261;
t220 = t239 + (-pkin(3) - t320) * t255;
t274 = (t260 * MDP(6) + t264 * MDP(7)) * t253;
t273 = (MDP(29) * t261 - MDP(30) * t257) * t234;
t205 = pkin(4) * t226 + t220;
t183 = pkin(3) * t224 - t187;
t272 = t211 * MDP(13) - t224 * MDP(15) + t187 * MDP(16) - t188 * MDP(17);
t271 = t201 * MDP(23) - t202 * MDP(24) + t279;
t173 = pkin(4) * t197 + t183;
t269 = t257 * MDP(29) + t261 * MDP(30) - t275 * (pkin(4) * t251 + pkin(12));
t268 = t171 * MDP(23) - t172 * MDP(24) + t296 + t302 - t304;
t161 = pkin(12) * t210 + t163;
t164 = pkin(5) * t185 - pkin(12) * t186 + t173;
t158 = -t161 * t257 + t164 * t261;
t159 = t161 * t261 + t164 * t257;
t267 = t175 * MDP(29) + t158 * MDP(32) - t159 * MDP(33) + t305 - t308;
t177 = -pkin(12) * t315 + t179;
t182 = pkin(5) * t203 - pkin(12) * t204 + t205;
t169 = -t177 * t257 + t182 * t261;
t170 = t177 * t261 + t182 * t257;
t266 = t200 * MDP(29) + t169 * MDP(32) - t170 * MDP(33) + t297 - t300;
t250 = t261 ^ 2;
t249 = t257 ^ 2;
t247 = t252 ^ 2;
t245 = -pkin(4) * t254 - pkin(5);
t229 = -pkin(9) * t314 + t240;
t215 = -t254 * t237 + t251 * t282;
t207 = pkin(5) * t233 - pkin(12) * t234 + t246;
t192 = t207 * t257 + t215 * t261;
t191 = t207 * t261 - t215 * t257;
t176 = pkin(5) * t315 - t178;
t160 = -pkin(5) * t210 - t162;
t1 = [t224 ^ 2 * MDP(15) + (t162 ^ 2 + t163 ^ 2 + t173 ^ 2) * MDP(26) + MDP(1) + (MDP(4) * t260 + 0.2e1 * MDP(5) * t264) * t317 + (-0.2e1 * t293 + t295) * t211 + (t197 * t328 + t303) * t198 + (t305 - 0.2e1 * t308) * t185 + (0.2e1 * t274 + t311) * t256 + (0.2e1 * t306 + t307 - 0.2e1 * t309) * t175 + (0.2e1 * t292 - 0.2e1 * t294 + t296 + 0.2e1 * t302 - 0.2e1 * t304) * t210 + (t171 * t210 + t183 * t197) * t327 + (-t172 * t210 + t183 * t198) * t326 + (t188 * t224 + t195 * t211) * t329 + (-t187 * t224 + t195 * t210) * t330 + 0.2e1 * (-pkin(1) * t317 - t231 * t256) * MDP(10) + (t158 * t185 + t160 * t174) * t324 + (-t159 * t185 + t160 * t175) * t323 + (-t162 * t186 - t163 * t185) * t325 + 0.2e1 * (t229 * t256 + t248 * t322) * MDP(9); (t183 * t227 + t198 * t220) * MDP(24) + (t183 * t226 + t197 * t220) * MDP(23) + (-t197 * t227 - t198 * t226) * MDP(19) + t229 * MDP(9) - t231 * MDP(10) + t175 * t299 + t185 * t297 + t198 * t291 + (-t174 * t200 - t175 * t199) * MDP(28) + (t158 * t203 + t160 * t199 + t169 * t185 + t174 * t176) * MDP(32) + (-t174 * t203 - t185 * t199) * MDP(30) + (t175 * t203 + t185 * t200) * MDP(29) + (-t159 * t203 + t160 * t200 - t170 * t185 + t175 * t176) * MDP(33) + (-t162 * t204 - t163 * t203 - t178 * t186 - t179 * t185) * MDP(25) + (t162 * t178 + t163 * t179 + t173 * t205) * MDP(26) + t311 + t280 * t224 + t272 * t255 + t274 + (t271 - t289) * t210 + ((-t210 * MDP(16) - t211 * MDP(17)) * pkin(2) + (-t210 * MDP(12) + t195 * MDP(17) - t293 + t295) * t259 + (-t195 * MDP(16) - t268 - t292 + t294) * t263) * t252; (t178 ^ 2 + t179 ^ 2 + t205 ^ 2) * MDP(26) + t247 * t259 ^ 2 * MDP(11) + MDP(8) + (t287 * t332 + t288) * t255 + (t226 * t328 + t291) * t227 + (t297 - 0.2e1 * t300) * t203 + (0.2e1 * t298 + t299 - 0.2e1 * t301) * t200 + ((0.2e1 * t259 * MDP(12) + t286) * t247 + (-t279 + t289) * t332) * t263 + (t228 * t255 + t247 * t320) * t330 + (-t230 * t255 - t247 * t321) * t329 + (-t201 * t315 + t220 * t226) * t327 + (t202 * t315 + t220 * t227) * t326 + (-t178 * t204 - t179 * t203) * t325 + (t169 * t203 + t176 * t199) * t324 + (-t170 * t203 + t176 * t200) * t323; t258 * t303 + (-t197 * t258 + t198 * t262) * MDP(19) + (-pkin(3) * t197 - t183 * t262) * MDP(23) + (-pkin(3) * t198 + t183 * t258) * MDP(24) + (-t185 * t215 + t186 * t213) * MDP(25) + (-t162 * t213 + t163 * t215 + t173 * t246) * MDP(26) + (t174 * t213 + t185 * t191) * MDP(32) + (t175 * t213 - t185 * t192) * MDP(33) + (-t163 * MDP(25) + t267) * t233 + (-t162 * MDP(25) + (-t175 * MDP(28) - t185 * MDP(30) + t160 * MDP(32)) * t257 + (t160 * MDP(33) + t306 + t307 - t309) * t261) * t234 + t331 * t210 + t272; t288 + t258 * t291 + (-t226 * t258 + t227 * t262) * MDP(19) + (-pkin(3) * t226 - t220 * t262) * MDP(23) + (-pkin(3) * t227 + t220 * t258) * MDP(24) + (-t203 * t215 + t204 * t213) * MDP(25) + (-t178 * t213 + t179 * t215 + t205 * t246) * MDP(26) + (t191 * t203 + t199 * t213) * MDP(32) + (-t192 * t203 + t200 * t213) * MDP(33) + (-t179 * MDP(25) + t266) * t233 + (-t178 * MDP(25) + (-MDP(28) * t200 - MDP(30) * t203 + MDP(32) * t176) * t257 + (t176 * MDP(33) + t298 + t299 - t301) * t261) * t234 + (-t263 * t331 + t287) * t252 - t280; MDP(15) + pkin(3) * t262 * t327 + (t213 ^ 2 + t215 ^ 2 + t246 ^ 2) * MDP(26) + (MDP(27) * t250 - 0.2e1 * t283) * t234 ^ 2 + (MDP(18) * t258 + 0.2e1 * MDP(19) * t262 - 0.2e1 * MDP(24) * pkin(3)) * t258 + (0.2e1 * t273 + t290) * t233 + (-t215 * t233 + t318) * t325 + (t191 * t233 + t257 * t318) * t324 + (-t192 * t233 + t261 * t318) * t323; t257 * t307 + (-t174 * t257 + t175 * t261) * MDP(28) + (-t160 * t261 + t174 * t245) * MDP(32) + (t160 * t257 + t175 * t245) * MDP(33) + t269 * t185 + ((-t185 * t251 - t186 * t254) * MDP(25) + (t162 * t254 + t163 * t251) * MDP(26)) * pkin(4) + t268; -t252 * t286 + t257 * t299 + (-t199 * t257 + t200 * t261) * MDP(28) + (-t176 * t261 + t199 * t245) * MDP(32) + (t176 * t257 + t200 * t245) * MDP(33) + t269 * t203 + ((-t203 * t251 - t204 * t254) * MDP(25) + (t178 * t254 + t179 * t251) * MDP(26)) * pkin(4) + t271; -t276 * t213 + t269 * t233 + (MDP(27) * t310 + (-t249 + t250) * MDP(28) + t275 * t245) * t234 + ((-t233 * t251 - t234 * t254) * MDP(25) + (-t213 * t254 + t215 * t251) * MDP(26)) * pkin(4) + t270; 0.2e1 * t283 + t249 * MDP(27) + MDP(22) + (t251 ^ 2 + t254 ^ 2) * MDP(26) * pkin(4) ^ 2 - 0.2e1 * t276 * t245; MDP(26) * t173 + t276 * t185; MDP(26) * t205 + t276 * t203; MDP(26) * t246 + t276 * t233; 0; MDP(26); t267; t266; MDP(32) * t191 - MDP(33) * t192 + t273 + t290; t269; t276; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
