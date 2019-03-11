% Calculate joint inertia matrix for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR15_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR15_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:19
% EndTime: 2019-03-10 00:50:26
% DurationCPUTime: 1.99s
% Computational Cost: add. (3368->359), mult. (8707->493), div. (0->0), fcn. (9729->12), ass. (0->158)
t240 = sin(pkin(7));
t335 = 0.2e1 * t240;
t241 = sin(pkin(6));
t243 = cos(pkin(6));
t246 = sin(qJ(3));
t247 = sin(qJ(2));
t250 = cos(qJ(3));
t242 = cos(pkin(7));
t318 = cos(qJ(2));
t277 = t242 * t318;
t311 = t240 * t246;
t201 = t243 * t311 + (t246 * t277 + t247 * t250) * t241;
t278 = t241 * t318;
t214 = t240 * t278 - t243 * t242;
t245 = sin(qJ(4));
t249 = cos(qJ(4));
t190 = t201 * t249 - t214 * t245;
t284 = MDP(18) + MDP(33);
t334 = t284 * t190;
t218 = t242 * t245 + t249 * t311;
t333 = t284 * t218;
t310 = t240 * t250;
t281 = pkin(10) * t310;
t317 = pkin(2) * t246;
t212 = t281 + (pkin(11) + t317) * t242;
t213 = (-pkin(3) * t250 - pkin(11) * t246 - pkin(2)) * t240;
t194 = t249 * t212 + t245 * t213;
t279 = qJ(5) * t310;
t191 = t279 - t194;
t332 = t191 * MDP(28);
t331 = t214 * MDP(15);
t222 = t243 * t247 * pkin(1) + pkin(9) * t278;
t272 = t241 * t277;
t196 = (t240 * t243 + t272) * pkin(10) + t222;
t280 = pkin(1) * t318;
t231 = t243 * t280;
t309 = t241 * t247;
t202 = pkin(2) * t243 + t231 + (-pkin(10) * t242 - pkin(9)) * t309;
t209 = (-pkin(10) * t240 * t247 - t318 * pkin(2) - pkin(1)) * t241;
t271 = t202 * t242 + t209 * t240;
t178 = -t246 * t196 + t271 * t250;
t330 = 2 * MDP(16);
t329 = 2 * MDP(17);
t328 = 2 * MDP(23);
t327 = 2 * MDP(24);
t326 = 2 * MDP(25);
t325 = 2 * MDP(26);
t324 = 2 * MDP(27);
t323 = 2 * MDP(31);
t322 = 2 * MDP(34);
t321 = 2 * MDP(35);
t320 = pkin(4) + pkin(12);
t319 = pkin(5) + pkin(11);
t316 = pkin(2) * t250;
t200 = -t243 * t310 + t246 * t309 - t250 * t272;
t315 = pkin(4) * t200;
t314 = MDP(28) * pkin(4);
t197 = t200 * qJ(5);
t228 = t319 * t249;
t313 = t228 * t249;
t235 = t241 ^ 2;
t312 = t235 * t247;
t308 = t243 * MDP(8);
t185 = -t202 * t240 + t242 * t209;
t173 = pkin(3) * t200 - pkin(11) * t201 + t185;
t179 = t196 * t250 + t271 * t246;
t177 = -pkin(11) * t214 + t179;
t168 = t245 * t173 + t249 * t177;
t307 = MDP(22) * t250;
t306 = MDP(28) * qJ(5);
t305 = MDP(28) * pkin(11) ^ 2;
t217 = -t249 * t242 + t245 * t311;
t244 = sin(qJ(6));
t248 = cos(qJ(6));
t204 = t217 * t244 - t248 * t310;
t304 = MDP(29) * t204;
t303 = MDP(29) * t248;
t302 = MDP(34) * t244;
t301 = MDP(35) * t248;
t167 = t173 * t249 - t245 * t177;
t166 = -t167 - t315;
t300 = t166 * MDP(28);
t189 = t201 * t245 + t214 * t249;
t180 = -t189 * t248 + t200 * t244;
t299 = t180 * MDP(32);
t181 = t189 * t244 + t200 * t248;
t298 = t181 * MDP(29);
t297 = t189 * MDP(19);
t296 = t189 * MDP(21);
t295 = t190 * MDP(20);
t193 = -t245 * t212 + t213 * t249;
t230 = pkin(4) * t310;
t192 = -t193 + t230;
t294 = t192 * MDP(28);
t293 = t200 * MDP(22);
t292 = t201 * MDP(11);
t291 = t201 * MDP(12);
t203 = t217 * t248 + t244 * t310;
t290 = t203 * MDP(32);
t289 = t214 * MDP(13);
t288 = t214 * MDP(14);
t287 = t217 * MDP(19);
t286 = t242 * MDP(15);
t285 = t246 * MDP(13);
t283 = MDP(23) - MDP(26);
t282 = MDP(24) - MDP(27);
t165 = -t197 - t168;
t276 = t248 * t244 * MDP(30);
t275 = -qJ(5) * t245 - pkin(3);
t274 = MDP(26) - t314;
t273 = -qJ(5) * MDP(25) - MDP(21);
t229 = pkin(10) * t311;
t219 = t242 * t316 - t229;
t221 = t242 * t317 + t281;
t270 = -t219 * MDP(16) + t221 * MDP(17);
t269 = t193 * MDP(23) - t194 * MDP(24);
t268 = -MDP(31) * t244 - MDP(32) * t248;
t267 = MDP(34) * t248 - t244 * MDP(35);
t211 = t229 + (-pkin(3) - t316) * t242;
t266 = MDP(19) + t268;
t265 = MDP(25) + t267;
t264 = (t247 * MDP(6) + t318 * MDP(7)) * t241;
t263 = MDP(14) * t242 - t218 * MDP(20) + t217 * MDP(21);
t262 = t167 * MDP(23) - t168 * MDP(24) + t293;
t261 = -qJ(5) * t218 + t211;
t176 = pkin(3) * t214 - t178;
t260 = t201 * MDP(13) + t178 * MDP(16) - t179 * MDP(17) - t331;
t162 = pkin(5) * t190 - t320 * t200 - t167;
t256 = -qJ(5) * t190 + t176;
t164 = t320 * t189 + t256;
t160 = t162 * t248 - t164 * t244;
t161 = t162 * t244 + t164 * t248;
t259 = t181 * MDP(31) + t160 * MDP(34) - t161 * MDP(35) - t299;
t182 = pkin(5) * t218 + pkin(12) * t310 + t192;
t183 = t320 * t217 + t261;
t170 = t182 * t248 - t183 * t244;
t171 = t182 * t244 + t183 * t248;
t258 = t204 * MDP(31) + t170 * MDP(34) - t171 * MDP(35) + t290;
t257 = (-MDP(34) * t320 + MDP(31)) * t248 + (MDP(35) * t320 - MDP(32)) * t244;
t255 = -pkin(4) * MDP(25) + MDP(20) + t257;
t169 = pkin(4) * t189 + t256;
t254 = t176 * MDP(24) + t166 * MDP(25) - t169 * MDP(27) + t259 - t297;
t188 = pkin(4) * t217 + t261;
t253 = t211 * MDP(24) + t192 * MDP(25) - t188 * MDP(27) + t258 - t287 + t333;
t239 = t249 ^ 2;
t238 = t248 ^ 2;
t237 = t245 ^ 2;
t236 = t244 ^ 2;
t234 = t240 ^ 2;
t227 = t319 * t245;
t225 = -pkin(4) * t249 + t275;
t223 = -t320 * t249 + t275;
t220 = -pkin(9) * t309 + t231;
t199 = t223 * t248 + t227 * t244;
t198 = -t223 * t244 + t227 * t248;
t184 = -pkin(5) * t217 - t191;
t163 = -pkin(5) * t189 - t165;
t1 = [MDP(1) + (t165 ^ 2 + t166 ^ 2 + t169 ^ 2) * MDP(28) + (MDP(4) * t247 + 0.2e1 * MDP(5) * t318) * t312 + (0.2e1 * t264 + t308) * t243 + (-0.2e1 * t180 * MDP(30) + t298) * t181 + 0.2e1 * (-pkin(1) * t312 - t222 * t243) * MDP(10) + 0.2e1 * (t220 * t243 + t235 * t280) * MDP(9) + (t185 * t329 - 0.2e1 * t289 + t292) * t201 + (-t178 * t330 + t179 * t329 + t331) * t214 + (t180 * t322 + t181 * t321) * t163 + (t165 * t326 - t169 * t325 + t176 * t328) * t189 + (-t165 * t324 + t166 * t325 + t167 * t328 - t168 * t327 + t185 * t330 + 0.2e1 * t288 - 0.2e1 * t291 + t293 + 0.2e1 * t295 - 0.2e1 * t296) * t200 + (t160 * t322 - t161 * t321 + t166 * t326 - t169 * t324 + t176 * t327 + t323 * t181 - 0.2e1 * t297 - 0.2e1 * t299 + t334) * t190; (t176 * t217 + t189 * t211) * MDP(23) + t204 * t298 + (-t169 * t217 - t188 * t189) * MDP(26) + (t165 * t191 + t166 * t192 + t169 * t188) * MDP(28) + (-t180 * t204 + t181 * t203) * MDP(30) + (-t163 * t203 + t180 * t184) * MDP(34) + (t165 * t217 + t189 * t191) * MDP(25) + (t163 * t204 + t181 * t184) * MDP(35) + t220 * MDP(9) - t222 * MDP(10) + t308 + t270 * t214 + t260 * t242 + t264 + t254 * t218 + t253 * t190 + (t192 * MDP(26) - t191 * MDP(27) - t263 + t269) * t200 + ((-t200 * MDP(16) - t201 * MDP(17)) * pkin(2) + (-t200 * MDP(12) + t185 * MDP(17) - t289 + t292) * t246 + (-t185 * MDP(16) - t166 * MDP(26) + t165 * MDP(27) - t262 - t288 + t291 - t295 + t296) * t250) * t240; t234 * t246 ^ 2 * MDP(11) + MDP(8) + (t188 ^ 2 + t192 ^ 2) * MDP(28) + (t285 * t335 + t286) * t242 + (0.2e1 * t203 * MDP(30) + t304) * t204 + ((0.2e1 * MDP(12) * t246 + t307) * t234 + t263 * t335) * t250 + (t219 * t242 + t234 * t316) * t330 + t194 * t310 * t327 + (-t193 * t310 + t211 * t217) * t328 + (-t188 * t217 - t192 * t310) * t325 + (-t221 * t242 - t234 * t317) * t329 + (t217 * t326 + t310 * t324 + t332) * t191 + (-t203 * t322 + t204 * t321) * t184 + (t170 * t322 - t171 * t321 - t188 * t324 + t192 * t326 + t323 * t204 + t211 * t327 - 0.2e1 * t287 + 0.2e1 * t290 + t333) * t218; -t200 * MDP(14) + (t180 * t228 + t190 * t198) * MDP(34) + (t181 * t228 - t190 * t199) * MDP(35) + (-t189 * MDP(23) - t190 * MDP(24)) * pkin(3) + (-t189 * MDP(26) - t190 * MDP(27) + t169 * MDP(28)) * t225 + (t200 * MDP(20) + t334 + (t190 * MDP(25) - t283 * t200 + t300) * pkin(11) + t254) * t245 + (t200 * MDP(21) - t176 * MDP(23) - t165 * MDP(25) + t169 * MDP(26) - t244 * t298 + (t180 * t244 - t181 * t248) * MDP(30) + t267 * t163 + t266 * t190 + (-t189 * MDP(25) - t165 * MDP(28) - t282 * t200) * pkin(11)) * t249 + t260; t286 + (t198 * t218 - t203 * t228) * MDP(34) + (-t199 * t218 + t204 * t228) * MDP(35) + (t250 * MDP(14) + t285) * t240 + (-t217 * MDP(23) - t218 * MDP(24)) * pkin(3) + (-t217 * MDP(26) - t218 * MDP(27) + t188 * MDP(28)) * t225 + (-MDP(20) * t310 + (t218 * MDP(25) + t283 * t310 + t294) * pkin(11) + t253) * t245 + (-MDP(21) * t310 - t211 * MDP(23) - t191 * MDP(25) + t188 * MDP(26) - t244 * t304 + (-t203 * t244 - t204 * t248) * MDP(30) + t267 * t184 + t266 * t218 + (-t217 * MDP(25) + t282 * t310 - t332) * pkin(11)) * t249 - t270; pkin(3) * t249 * t328 + MDP(15) + (MDP(28) * t225 + t249 * t325) * t225 + (MDP(29) * t236 + 0.2e1 * t276 + t305) * t239 + (t284 + t305) * t237 + 0.2e1 * (-pkin(3) * MDP(24) - t225 * MDP(27) + t266 * t249) * t245 + (t198 * t245 + t248 * t313) * t322 + (-t199 * t245 - t244 * t313) * t321 + (t237 + t239) * pkin(11) * t326; (-t167 - 0.2e1 * t315) * MDP(26) + (-t165 + t197) * MDP(27) + (-pkin(4) * t166 - qJ(5) * t165) * MDP(28) + t248 * t298 + (-t180 * t248 - t181 * t244) * MDP(30) + (qJ(5) * t180 + t163 * t244) * MDP(34) + (qJ(5) * t181 + t163 * t248) * MDP(35) + t273 * t189 + t255 * t190 + t262; -t240 * t307 + (-t193 + 0.2e1 * t230) * MDP(26) + (-0.2e1 * t279 + t194) * MDP(27) + (-pkin(4) * t192 - qJ(5) * t191) * MDP(28) + t204 * t303 + (t203 * t248 - t204 * t244) * MDP(30) + (-qJ(5) * t203 + t184 * t244) * MDP(34) + (qJ(5) * t204 + t184 * t248) * MDP(35) + t273 * t217 + t255 * t218 + t269; (t301 + t302) * t228 + t255 * t245 + (MDP(21) - t244 * t303 + (t236 - t238) * MDP(30) + t265 * qJ(5)) * t249 + ((-t282 + t306) * t249 + (-MDP(23) + t274) * t245) * pkin(11); -0.2e1 * t276 + t238 * MDP(29) + MDP(22) + (-(2 * MDP(26)) + t314) * pkin(4) + (t324 + 0.2e1 * t301 + 0.2e1 * t302 + t306) * qJ(5); t200 * MDP(26) + t190 * t265 + t300; -MDP(26) * t310 + t218 * t265 + t294; (MDP(28) * pkin(11) + t265) * t245; t274; MDP(28); t190 * MDP(33) + t259; t218 * MDP(33) + t258; MDP(33) * t245 + MDP(34) * t198 - MDP(35) * t199 + t268 * t249; t257; t267; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
