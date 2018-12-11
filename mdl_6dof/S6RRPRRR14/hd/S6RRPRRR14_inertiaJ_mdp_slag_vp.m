% Calculate joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR14_inertiaJ_mdp_slag_vp: MDP has to be [35x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:32:44
% EndTime: 2018-12-10 18:32:57
% DurationCPUTime: 2.96s
% Computational Cost: add. (6903->386), mult. (19234->583), div. (0->0), fcn. (22622->16), ass. (0->166)
t274 = sin(pkin(8));
t283 = sin(qJ(4));
t337 = t274 * t283;
t356 = cos(qJ(4));
t310 = t274 * t356;
t278 = cos(pkin(8));
t275 = sin(pkin(7));
t279 = cos(pkin(7));
t280 = cos(pkin(6));
t276 = sin(pkin(6));
t287 = cos(qJ(2));
t334 = t276 * t287;
t305 = -t275 * t334 + t279 * t280;
t273 = sin(pkin(14));
t284 = sin(qJ(2));
t335 = t276 * t284;
t259 = t273 * t335;
t277 = cos(pkin(14));
t311 = t279 * t334;
t298 = t275 * t280 + t311;
t371 = -t298 * t277 + t259;
t372 = -t274 * t305 + t278 * t371;
t281 = sin(qJ(6));
t285 = cos(qJ(6));
t292 = -(MDP(34) * t281 + MDP(35) * t285) * pkin(13) + t281 * MDP(31) + t285 * MDP(32);
t370 = 2 * MDP(11);
t369 = 2 * MDP(12);
t368 = 2 * MDP(13);
t367 = -2 * MDP(16);
t366 = -2 * MDP(17);
t365 = 2 * MDP(18);
t364 = 2 * MDP(20);
t363 = 2 * MDP(21);
t362 = -2 * MDP(25);
t361 = 2 * MDP(27);
t360 = 2 * MDP(28);
t359 = -2 * MDP(30);
t358 = 0.2e1 * MDP(34);
t357 = 0.2e1 * MDP(35);
t355 = pkin(1) * t287;
t268 = t275 ^ 2;
t354 = pkin(2) * t268;
t353 = pkin(2) * t275;
t352 = pkin(11) * t274;
t351 = pkin(11) * t278;
t350 = pkin(12) * t281;
t349 = pkin(12) * t285;
t286 = cos(qJ(5));
t348 = pkin(12) * t286;
t347 = MDP(27) * pkin(4);
t346 = pkin(4) * MDP(28);
t256 = t280 * t284 * pkin(1) + pkin(10) * t334;
t230 = t298 * qJ(3) + t256;
t265 = t280 * t355;
t235 = pkin(2) * t280 + t265 + (-qJ(3) * t279 - pkin(10)) * t335;
t243 = (-qJ(3) * t275 * t284 - pkin(2) * t287 - pkin(1)) * t276;
t338 = t273 * t279;
t339 = t273 * t275;
t210 = t277 * t230 + t235 * t338 + t243 * t339;
t198 = -pkin(11) * t372 + t210;
t332 = t277 * t279;
t336 = t275 * t277;
t209 = -t273 * t230 + t235 * t332 + t243 * t336;
t231 = t298 * t273 + t277 * t335;
t201 = t305 * pkin(3) - t231 * t351 + t209;
t217 = -t275 * t235 + t279 * t243;
t203 = pkin(3) * t371 - t231 * t352 + t217;
t188 = t356 * t198 + (t201 * t278 + t203 * t274) * t283;
t219 = -t274 * t371 - t278 * t305;
t185 = -t219 * pkin(12) + t188;
t189 = -t201 * t274 + t278 * t203;
t212 = t231 * t283 + t356 * t372;
t213 = t356 * t231 - t283 * t372;
t186 = pkin(4) * t212 - pkin(12) * t213 + t189;
t282 = sin(qJ(5));
t179 = -t185 * t282 + t186 * t286;
t177 = -pkin(5) * t212 - t179;
t345 = t177 * t281;
t344 = t177 * t285;
t264 = pkin(2) * t332;
t234 = pkin(3) * t279 + t264 + (-qJ(3) - t351) * t339;
t242 = (-pkin(3) * t277 - t273 * t352 - pkin(2)) * t275;
t216 = -t234 * t274 + t278 * t242;
t309 = t278 * t356;
t232 = -t279 * t310 + t283 * t339 - t309 * t336;
t333 = t277 * t278;
t233 = t279 * t337 + (t356 * t273 + t283 * t333) * t275;
t204 = pkin(4) * t232 - pkin(12) * t233 + t216;
t252 = pkin(2) * t338 + qJ(3) * t336;
t229 = (t274 * t279 + t275 * t333) * pkin(11) + t252;
t208 = t356 * t229 + (t234 * t278 + t242 * t274) * t283;
t249 = t274 * t336 - t278 * t279;
t206 = -t249 * pkin(12) + t208;
t194 = t204 * t286 - t206 * t282;
t190 = -pkin(5) * t232 - t194;
t343 = t190 * t281;
t342 = t190 * t285;
t253 = -t286 * t278 + t282 * t337;
t341 = t253 * t282;
t269 = t276 ^ 2;
t340 = t269 * t284;
t331 = t280 * MDP(8);
t330 = MDP(15) * t233;
t222 = t233 * t286 - t249 * t282;
t215 = t222 * t285 + t232 * t281;
t329 = MDP(29) * t215;
t328 = MDP(29) * t281;
t327 = MDP(29) * t285;
t326 = MDP(33) * t286;
t200 = t213 * t286 - t219 * t282;
t192 = t200 * t281 - t285 * t212;
t325 = t192 * MDP(32);
t193 = t200 * t285 + t212 * t281;
t324 = t193 * MDP(31);
t199 = t213 * t282 + t286 * t219;
t323 = t199 * MDP(33);
t322 = t200 * MDP(23);
t321 = t200 * MDP(24);
t320 = t212 * MDP(26);
t214 = t222 * t281 - t285 * t232;
t319 = t214 * MDP(32);
t318 = t215 * MDP(31);
t221 = t233 * t282 + t286 * t249;
t317 = t221 * MDP(33);
t316 = t222 * MDP(22);
t315 = t222 * MDP(23);
t314 = t222 * MDP(24);
t313 = t232 * MDP(26);
t312 = t249 * MDP(19);
t308 = MDP(30) * t281 * t285;
t307 = -pkin(12) * MDP(27) + MDP(24);
t306 = -pkin(12) * MDP(28) + MDP(25);
t180 = t185 * t286 + t186 * t282;
t195 = t204 * t282 + t206 * t286;
t303 = MDP(31) * t285 - MDP(32) * t281;
t258 = -pkin(5) * t286 - pkin(13) * t282 - pkin(4);
t245 = t258 * t285 - t281 * t348;
t246 = t258 * t281 + t285 * t348;
t301 = MDP(34) * t245 - MDP(35) * t246;
t300 = MDP(34) * t285 - MDP(35) * t281;
t297 = -MDP(23) + t303;
t296 = t301 - t347;
t187 = -t283 * t198 + t201 * t309 + t203 * t310;
t207 = -t283 * t229 + t234 * t309 + t242 * t310;
t184 = t219 * pkin(4) - t187;
t205 = t249 * pkin(4) - t207;
t178 = pkin(13) * t212 + t180;
t181 = t199 * pkin(5) - t200 * pkin(13) + t184;
t175 = -t178 * t281 + t181 * t285;
t176 = t178 * t285 + t181 * t281;
t291 = t175 * MDP(34) - t176 * MDP(35) + t323 + t324 - t325;
t191 = pkin(13) * t232 + t195;
t196 = t221 * pkin(5) - t222 * pkin(13) + t205;
t182 = -t191 * t281 + t196 * t285;
t183 = t191 * t285 + t196 * t281;
t290 = t182 * MDP(34) - t183 * MDP(35) + t317 + t318 - t319;
t289 = -MDP(25) + t292;
t272 = t285 ^ 2;
t271 = t282 ^ 2;
t270 = t281 ^ 2;
t255 = -pkin(10) * t335 + t265;
t254 = t278 * t282 + t286 * t337;
t251 = -qJ(3) * t339 + t264;
t239 = t285 * t254 - t281 * t310;
t238 = -t281 * t254 - t285 * t310;
t1 = [t219 ^ 2 * MDP(19) + t200 ^ 2 * MDP(22) + (t209 ^ 2 + t210 ^ 2 + t217 ^ 2) * MDP(14) + MDP(1) + (MDP(4) * t284 + 0.2e1 * MDP(5) * t287) * t340 + (MDP(15) * t213 + t219 * t366) * t213 + (MDP(29) * t193 + t192 * t359) * t193 + (t331 + 0.2e1 * (MDP(6) * t284 + MDP(7) * t287) * t276) * t280 + (t213 * t367 + t219 * t365 + t320 + 0.2e1 * t321) * t212 + (t212 * t362 - 0.2e1 * t322 + t323 + 0.2e1 * t324 - 0.2e1 * t325) * t199 + (t175 * t199 + t177 * t192) * t358 + (-t176 * t199 + t177 * t193) * t357 + (t179 * t212 + t184 * t199) * t361 + (-t180 * t212 + t184 * t200) * t360 + (t188 * t219 + t189 * t213) * t363 + (-t187 * t219 + t189 * t212) * t364 + (-t210 * t305 + t217 * t231) * t369 + 0.2e1 * (-pkin(1) * t340 - t256 * t280) * MDP(10) + 0.2e1 * (t255 * t280 + t269 * t355) * MDP(9) + (-t209 * t231 - t210 * t371) * t368 + (t209 * t305 + t217 * t371) * t370; (-t252 * t305 - t210 * t279 + (-pkin(2) * t231 + t217 * t273) * t275) * MDP(12) + (-t251 * t231 + t252 * (t277 * t311 - t259) + (-t209 * t273 + (t252 * t280 + t210) * t277) * t275) * MDP(13) + (t209 * t251 + t210 * t252 - t217 * t353) * MDP(14) + t331 + (-t192 * t215 - t193 * t214) * MDP(30) + (t175 * t221 + t177 * t214 + t182 * t199 + t190 * t192) * MDP(34) + (-t192 * t221 - t199 * t214) * MDP(32) + (t193 * t221 + t199 * t215) * MDP(31) + (-t176 * t221 + t177 * t215 - t183 * t199 + t190 * t193) * MDP(35) + (-t199 * t222 - t200 * t221) * MDP(23) + (t179 * t232 + t184 * t221 + t194 * t212 + t199 * t205) * MDP(27) + (-t199 * t232 - t212 * t221) * MDP(25) + (t200 * t232 + t212 * t222) * MDP(24) + (-t180 * t232 + t184 * t222 - t195 * t212 + t200 * t205) * MDP(28) + (-t212 * t233 - t213 * t232) * MDP(16) + (t188 * t249 + t189 * t233 + t208 * t219 + t213 * t216) * MDP(21) + (-t187 * t249 + t189 * t232 - t207 * t219 + t212 * t216) * MDP(20) + (t212 * t249 + t219 * t232) * MDP(18) + (-t213 * t249 - t219 * t233) * MDP(17) + t255 * MDP(9) - t256 * MDP(10) + (t209 * t279 - t217 * t336 + t251 * t305 - t353 * t371) * MDP(11) + t219 * t312 + t212 * t313 + t200 * t316 + t199 * t317 + t193 * t329 + t213 * t330 + MDP(7) * t334 + MDP(6) * t335; t249 ^ 2 * MDP(19) + t222 ^ 2 * MDP(22) + (pkin(2) ^ 2 * t268 + t251 ^ 2 + t252 ^ 2) * MDP(14) + MDP(8) + (t249 * t366 + t330) * t233 + (t214 * t359 + t329) * t215 + (t233 * t367 + t249 * t365 + t313 + 0.2e1 * t314) * t232 + (t232 * t362 - 0.2e1 * t315 + t317 + 0.2e1 * t318 - 0.2e1 * t319) * t221 + (t182 * t221 + t190 * t214) * t358 + (-t183 * t221 + t190 * t215) * t357 + (t194 * t232 + t205 * t221) * t361 + (-t195 * t232 + t205 * t222) * t360 + (t208 * t249 + t216 * t233) * t363 + (-t207 * t249 + t216 * t232) * t364 + (-t252 * t279 - t273 * t354) * t369 + (t251 * t279 + t277 * t354) * t370 + (-t251 * t273 + t252 * t277) * t275 * t368; t371 * MDP(11) + t231 * MDP(12) + t217 * MDP(14) + (t278 * t212 - t219 * t310) * MDP(20) + (t213 * t278 + t219 * t337) * MDP(21) + (-t199 * t310 - t253 * t212) * MDP(27) + (-t200 * t310 - t254 * t212) * MDP(28) + (t192 * t253 + t199 * t238) * MDP(34) + (t193 * t253 - t199 * t239) * MDP(35); (t278 * t232 - t249 * t310) * MDP(20) + (t233 * t278 + t249 * t337) * MDP(21) + (-t221 * t310 - t253 * t232) * MDP(27) + (-t222 * t310 - t254 * t232) * MDP(28) + (t214 * t253 + t221 * t238) * MDP(34) + (t215 * t253 - t221 * t239) * MDP(35) + (-MDP(11) * t277 + MDP(12) * t273 - MDP(14) * pkin(2)) * t275; MDP(14); -t200 * t346 + t213 * MDP(17) - t212 * MDP(18) - t219 * MDP(19) + t187 * MDP(20) - t188 * MDP(21) + t296 * t199 + (-t184 * MDP(27) + t306 * t212 - t291 + t322) * t286 + (t200 * MDP(22) + t184 * MDP(28) + t193 * t327 + (-t192 * t285 - t193 * t281) * MDP(30) + (pkin(12) * t192 + t345) * MDP(34) + (pkin(12) * t193 + t344) * MDP(35) + t307 * t212 + t297 * t199) * t282; -t222 * t346 + t233 * MDP(17) - t232 * MDP(18) - t312 + t207 * MDP(20) - t208 * MDP(21) + t296 * t221 + (-t205 * MDP(27) + t306 * t232 - t290 + t315) * t286 + (t316 + t205 * MDP(28) + t215 * t327 + (-t214 * t285 - t215 * t281) * MDP(30) + (pkin(12) * t214 + t343) * MDP(34) + (pkin(12) * t215 + t342) * MDP(35) + t307 * t232 + t297 * t221) * t282; (-t238 * t286 + t281 * t341) * MDP(34) + (t239 * t286 + t285 * t341) * MDP(35) - t337 * MDP(21) + (MDP(27) * t286 - MDP(28) * t282 + MDP(20)) * t310; MDP(19) + (t326 + (2 * t347)) * t286 + (MDP(29) * t272 + MDP(22) - 0.2e1 * t308) * t271 + (-t245 * t286 + t271 * t350) * t358 + (t246 * t286 + t271 * t349) * t357 + 0.2e1 * (-t286 * t297 - t346) * t282; t321 + t320 + t179 * MDP(27) - t180 * MDP(28) + t193 * t328 + (-t192 * t281 + t193 * t285) * MDP(30) + (-pkin(5) * t192 - t344) * MDP(34) + (-pkin(5) * t193 + t345) * MDP(35) + t289 * t199; t314 + t313 + t194 * MDP(27) - t195 * MDP(28) + t215 * t328 + (-t214 * t281 + t215 * t285) * MDP(30) + (-pkin(5) * t214 - t342) * MDP(34) + (-pkin(5) * t215 + t343) * MDP(35) + t289 * t221; -MDP(28) * t254 + (-MDP(27) - t300) * t253; (-t292 + t306) * t286 + (t281 * t327 + (-t270 + t272) * MDP(30) + (-pkin(5) * t281 - t349) * MDP(34) + (-pkin(5) * t285 + t350) * MDP(35) + t307) * t282; MDP(29) * t270 + 0.2e1 * pkin(5) * t300 + MDP(26) + 0.2e1 * t308; t291; t290; MDP(34) * t238 - MDP(35) * t239; t303 * t282 + t301 - t326; t292; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
