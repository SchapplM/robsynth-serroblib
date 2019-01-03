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
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
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
% StartTime: 2019-01-03 10:20:54
% EndTime: 2019-01-03 10:21:07
% DurationCPUTime: 3.01s
% Computational Cost: add. (6903->386), mult. (19234->583), div. (0->0), fcn. (22622->16), ass. (0->166)
t275 = sin(pkin(8));
t284 = sin(qJ(4));
t338 = t275 * t284;
t357 = cos(qJ(4));
t311 = t275 * t357;
t279 = cos(pkin(8));
t276 = sin(pkin(7));
t280 = cos(pkin(7));
t281 = cos(pkin(6));
t277 = sin(pkin(6));
t288 = cos(qJ(2));
t335 = t277 * t288;
t306 = -t276 * t335 + t280 * t281;
t274 = sin(pkin(14));
t285 = sin(qJ(2));
t336 = t277 * t285;
t261 = t274 * t336;
t278 = cos(pkin(14));
t312 = t280 * t335;
t299 = t276 * t281 + t312;
t372 = -t299 * t278 + t261;
t373 = -t306 * t275 + t279 * t372;
t282 = sin(qJ(6));
t286 = cos(qJ(6));
t293 = -(MDP(34) * t282 + MDP(35) * t286) * pkin(13) + t282 * MDP(31) + t286 * MDP(32);
t371 = 2 * MDP(11);
t370 = 2 * MDP(12);
t369 = 2 * MDP(13);
t368 = -2 * MDP(16);
t367 = -2 * MDP(17);
t366 = 2 * MDP(18);
t365 = 2 * MDP(20);
t364 = 2 * MDP(21);
t363 = -2 * MDP(25);
t362 = 2 * MDP(27);
t361 = 2 * MDP(28);
t360 = -2 * MDP(30);
t359 = 0.2e1 * MDP(34);
t358 = 0.2e1 * MDP(35);
t356 = pkin(1) * t288;
t269 = t276 ^ 2;
t355 = pkin(2) * t269;
t354 = pkin(2) * t276;
t353 = pkin(11) * t275;
t352 = pkin(11) * t279;
t351 = pkin(12) * t282;
t350 = pkin(12) * t286;
t287 = cos(qJ(5));
t349 = pkin(12) * t287;
t348 = MDP(27) * pkin(4);
t347 = MDP(28) * pkin(4);
t257 = t281 * t285 * pkin(1) + pkin(10) * t335;
t231 = t299 * qJ(3) + t257;
t266 = t281 * t356;
t236 = pkin(2) * t281 + t266 + (-qJ(3) * t280 - pkin(10)) * t336;
t244 = (-qJ(3) * t276 * t285 - pkin(2) * t288 - pkin(1)) * t277;
t339 = t274 * t280;
t340 = t274 * t276;
t211 = t278 * t231 + t236 * t339 + t244 * t340;
t199 = -pkin(11) * t373 + t211;
t333 = t278 * t280;
t337 = t276 * t278;
t210 = -t274 * t231 + t236 * t333 + t244 * t337;
t232 = t299 * t274 + t278 * t336;
t202 = t306 * pkin(3) - t232 * t352 + t210;
t218 = -t276 * t236 + t280 * t244;
t204 = pkin(3) * t372 - t232 * t353 + t218;
t189 = t357 * t199 + (t202 * t279 + t204 * t275) * t284;
t220 = -t275 * t372 - t306 * t279;
t186 = -t220 * pkin(12) + t189;
t190 = -t202 * t275 + t279 * t204;
t213 = t232 * t284 + t357 * t373;
t214 = t232 * t357 - t284 * t373;
t187 = pkin(4) * t213 - pkin(12) * t214 + t190;
t283 = sin(qJ(5));
t180 = -t186 * t283 + t187 * t287;
t178 = -pkin(5) * t213 - t180;
t346 = t178 * t282;
t345 = t178 * t286;
t265 = pkin(2) * t333;
t235 = pkin(3) * t280 + t265 + (-qJ(3) - t352) * t340;
t243 = (-pkin(3) * t278 - t274 * t353 - pkin(2)) * t276;
t217 = -t235 * t275 + t279 * t243;
t310 = t279 * t357;
t233 = -t280 * t311 + t284 * t340 - t310 * t337;
t334 = t278 * t279;
t234 = t280 * t338 + (t357 * t274 + t284 * t334) * t276;
t205 = pkin(4) * t233 - pkin(12) * t234 + t217;
t253 = pkin(2) * t339 + qJ(3) * t337;
t230 = (t275 * t280 + t276 * t334) * pkin(11) + t253;
t209 = t357 * t230 + (t235 * t279 + t243 * t275) * t284;
t250 = t275 * t337 - t279 * t280;
t207 = -t250 * pkin(12) + t209;
t195 = t205 * t287 - t207 * t283;
t191 = -pkin(5) * t233 - t195;
t344 = t191 * t282;
t343 = t191 * t286;
t254 = -t287 * t279 + t283 * t338;
t342 = t254 * t283;
t270 = t277 ^ 2;
t341 = t270 * t285;
t332 = t281 * MDP(8);
t331 = MDP(15) * t234;
t223 = t234 * t287 - t250 * t283;
t216 = t223 * t286 + t233 * t282;
t330 = MDP(29) * t216;
t329 = MDP(29) * t282;
t328 = MDP(29) * t286;
t327 = MDP(33) * t287;
t201 = t214 * t287 - t220 * t283;
t193 = t201 * t282 - t213 * t286;
t326 = t193 * MDP(32);
t194 = t201 * t286 + t213 * t282;
t325 = t194 * MDP(31);
t200 = t214 * t283 + t220 * t287;
t324 = t200 * MDP(33);
t323 = t201 * MDP(23);
t322 = t201 * MDP(24);
t321 = t213 * MDP(26);
t215 = t223 * t282 - t286 * t233;
t320 = t215 * MDP(32);
t319 = t216 * MDP(31);
t222 = t234 * t283 + t287 * t250;
t318 = t222 * MDP(33);
t317 = t223 * MDP(22);
t316 = t223 * MDP(23);
t315 = t223 * MDP(24);
t314 = t233 * MDP(26);
t313 = t250 * MDP(19);
t309 = MDP(30) * t282 * t286;
t308 = -pkin(12) * MDP(27) + MDP(24);
t307 = -pkin(12) * MDP(28) + MDP(25);
t181 = t186 * t287 + t187 * t283;
t196 = t205 * t283 + t207 * t287;
t304 = MDP(31) * t286 - MDP(32) * t282;
t259 = -pkin(5) * t287 - pkin(13) * t283 - pkin(4);
t246 = t259 * t286 - t282 * t349;
t247 = t259 * t282 + t286 * t349;
t302 = t246 * MDP(34) - t247 * MDP(35);
t301 = MDP(34) * t286 - MDP(35) * t282;
t298 = -MDP(23) + t304;
t297 = t302 - t348;
t188 = -t284 * t199 + t202 * t310 + t204 * t311;
t208 = -t284 * t230 + t235 * t310 + t243 * t311;
t185 = t220 * pkin(4) - t188;
t206 = t250 * pkin(4) - t208;
t179 = pkin(13) * t213 + t181;
t182 = t200 * pkin(5) - t201 * pkin(13) + t185;
t176 = -t179 * t282 + t182 * t286;
t177 = t179 * t286 + t182 * t282;
t292 = t176 * MDP(34) - t177 * MDP(35) + t324 + t325 - t326;
t192 = pkin(13) * t233 + t196;
t197 = t222 * pkin(5) - t223 * pkin(13) + t206;
t183 = -t192 * t282 + t197 * t286;
t184 = t192 * t286 + t197 * t282;
t291 = t183 * MDP(34) - t184 * MDP(35) + t318 + t319 - t320;
t290 = -MDP(25) + t293;
t273 = t286 ^ 2;
t272 = t283 ^ 2;
t271 = t282 ^ 2;
t256 = -pkin(10) * t336 + t266;
t255 = t279 * t283 + t287 * t338;
t252 = -qJ(3) * t340 + t265;
t240 = t286 * t255 - t282 * t311;
t239 = -t282 * t255 - t286 * t311;
t1 = [t220 ^ 2 * MDP(19) + (t210 ^ 2 + t211 ^ 2 + t218 ^ 2) * MDP(14) + t201 ^ 2 * MDP(22) + MDP(1) + (MDP(4) * t285 + 0.2e1 * MDP(5) * t288) * t341 + (MDP(15) * t214 + t220 * t367) * t214 + (MDP(29) * t194 + t193 * t360) * t194 + (t332 + 0.2e1 * (MDP(6) * t285 + MDP(7) * t288) * t277) * t281 + (t214 * t368 + t220 * t366 + t321 + 0.2e1 * t322) * t213 + (t213 * t363 - 0.2e1 * t323 + t324 + 0.2e1 * t325 - 0.2e1 * t326) * t200 + (-t210 * t232 - t211 * t372) * t369 + (t210 * t306 + t218 * t372) * t371 + 0.2e1 * (t256 * t281 + t270 * t356) * MDP(9) + 0.2e1 * (-pkin(1) * t341 - t257 * t281) * MDP(10) + (-t211 * t306 + t218 * t232) * t370 + (-t177 * t200 + t178 * t194) * t358 + (t176 * t200 + t178 * t193) * t359 + (t180 * t213 + t185 * t200) * t362 + (-t181 * t213 + t185 * t201) * t361 + (t189 * t220 + t190 * t214) * t364 + (-t188 * t220 + t190 * t213) * t365; (t210 * t280 - t218 * t337 + t252 * t306 - t354 * t372) * MDP(11) + (-t252 * t232 + t253 * (t278 * t312 - t261) + (-t210 * t274 + (t253 * t281 + t211) * t278) * t276) * MDP(13) + MDP(7) * t335 + MDP(6) * t336 + (-t253 * t306 - t211 * t280 + (-pkin(2) * t232 + t218 * t274) * t276) * MDP(12) + t220 * t313 + t213 * t314 + t201 * t317 + t200 * t318 + t194 * t330 + t214 * t331 + t256 * MDP(9) - t257 * MDP(10) + (t189 * t250 + t190 * t234 + t209 * t220 + t214 * t217) * MDP(21) + (-t214 * t250 - t220 * t234) * MDP(17) + (-t188 * t250 + t190 * t233 - t208 * t220 + t213 * t217) * MDP(20) + (t213 * t250 + t220 * t233) * MDP(18) + (-t213 * t234 - t214 * t233) * MDP(16) + (t201 * t233 + t213 * t223) * MDP(24) + (t180 * t233 + t185 * t222 + t195 * t213 + t200 * t206) * MDP(27) + (-t200 * t233 - t213 * t222) * MDP(25) + (-t181 * t233 + t185 * t223 - t196 * t213 + t201 * t206) * MDP(28) + (-t200 * t223 - t201 * t222) * MDP(23) + (t194 * t222 + t200 * t216) * MDP(31) + (t176 * t222 + t178 * t215 + t183 * t200 + t191 * t193) * MDP(34) + (-t193 * t222 - t200 * t215) * MDP(32) + (-t177 * t222 + t178 * t216 - t184 * t200 + t191 * t194) * MDP(35) + (-t193 * t216 - t194 * t215) * MDP(30) + (t210 * t252 + t211 * t253 - t218 * t354) * MDP(14) + t332; t250 ^ 2 * MDP(19) + (pkin(2) ^ 2 * t269 + t252 ^ 2 + t253 ^ 2) * MDP(14) + t223 ^ 2 * MDP(22) + MDP(8) + (t250 * t367 + t331) * t234 + (t215 * t360 + t330) * t216 + (t234 * t368 + t250 * t366 + t314 + 0.2e1 * t315) * t233 + (t233 * t363 - 0.2e1 * t316 + t318 + 0.2e1 * t319 - 0.2e1 * t320) * t222 + (t252 * t280 + t278 * t355) * t371 + (-t253 * t280 - t274 * t355) * t370 + (t183 * t222 + t191 * t215) * t359 + (-t196 * t233 + t206 * t223) * t361 + (t209 * t250 + t217 * t234) * t364 + (-t208 * t250 + t217 * t233) * t365 + (-t184 * t222 + t191 * t216) * t358 + (t195 * t233 + t206 * t222) * t362 + (-t252 * t274 + t253 * t278) * t276 * t369; t372 * MDP(11) + t232 * MDP(12) + t218 * MDP(14) + (t279 * t213 - t220 * t311) * MDP(20) + (t214 * t279 + t220 * t338) * MDP(21) + (-t200 * t311 - t254 * t213) * MDP(27) + (-t201 * t311 - t255 * t213) * MDP(28) + (t193 * t254 + t200 * t239) * MDP(34) + (t194 * t254 - t200 * t240) * MDP(35); (t279 * t233 - t250 * t311) * MDP(20) + (t234 * t279 + t250 * t338) * MDP(21) + (-t222 * t311 - t254 * t233) * MDP(27) + (-t223 * t311 - t255 * t233) * MDP(28) + (t215 * t254 + t222 * t239) * MDP(34) + (t216 * t254 - t222 * t240) * MDP(35) + (-MDP(11) * t278 + MDP(12) * t274 - MDP(14) * pkin(2)) * t276; MDP(14); -t201 * t347 + t214 * MDP(17) - t213 * MDP(18) - t220 * MDP(19) + t188 * MDP(20) - t189 * MDP(21) + t297 * t200 + (-t185 * MDP(27) + t307 * t213 - t292 + t323) * t287 + (t201 * MDP(22) + t185 * MDP(28) + t194 * t328 + (-t193 * t286 - t194 * t282) * MDP(30) + (pkin(12) * t193 + t346) * MDP(34) + (pkin(12) * t194 + t345) * MDP(35) + t308 * t213 + t298 * t200) * t283; -t223 * t347 + t234 * MDP(17) - t233 * MDP(18) - t313 + t208 * MDP(20) - t209 * MDP(21) + t297 * t222 + (-t206 * MDP(27) + t307 * t233 - t291 + t316) * t287 + (t317 + t206 * MDP(28) + t216 * t328 + (-t215 * t286 - t216 * t282) * MDP(30) + (pkin(12) * t215 + t344) * MDP(34) + (pkin(12) * t216 + t343) * MDP(35) + t308 * t233 + t298 * t222) * t283; (-t239 * t287 + t282 * t342) * MDP(34) + (t240 * t287 + t286 * t342) * MDP(35) - t338 * MDP(21) + (MDP(27) * t287 - MDP(28) * t283 + MDP(20)) * t311; MDP(19) + (t327 + (2 * t348)) * t287 + (MDP(29) * t273 + MDP(22) - 0.2e1 * t309) * t272 + (-t246 * t287 + t272 * t351) * t359 + (t247 * t287 + t272 * t350) * t358 + 0.2e1 * (-t287 * t298 - t347) * t283; t322 + t321 + t180 * MDP(27) - t181 * MDP(28) + t194 * t329 + (-t193 * t282 + t194 * t286) * MDP(30) + (-pkin(5) * t193 - t345) * MDP(34) + (-pkin(5) * t194 + t346) * MDP(35) + t290 * t200; t315 + t314 + t195 * MDP(27) - t196 * MDP(28) + t216 * t329 + (-t215 * t282 + t216 * t286) * MDP(30) + (-pkin(5) * t215 - t343) * MDP(34) + (-pkin(5) * t216 + t344) * MDP(35) + t290 * t222; -MDP(28) * t255 + (-MDP(27) - t301) * t254; (-t293 + t307) * t287 + (t282 * t328 + (-t271 + t273) * MDP(30) + (-pkin(5) * t282 - t350) * MDP(34) + (-pkin(5) * t286 + t351) * MDP(35) + t308) * t283; MDP(29) * t271 + 0.2e1 * pkin(5) * t301 + MDP(26) + 0.2e1 * t309; t292; t291; MDP(34) * t239 - MDP(35) * t240; t304 * t283 + t302 - t327; t293; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
