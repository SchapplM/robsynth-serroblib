% Calculate joint inertia matrix for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR10V2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR10V2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR10V2_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:55
% EndTime: 2019-04-11 14:52:01
% DurationCPUTime: 1.55s
% Computational Cost: add. (1246->261), mult. (2738->374), div. (0->0), fcn. (3076->10), ass. (0->137)
t262 = sin(qJ(3));
t263 = sin(qJ(2));
t267 = cos(qJ(2));
t346 = cos(qJ(3));
t223 = t262 * t263 - t346 * t267;
t224 = t262 * t267 + t346 * t263;
t260 = sin(qJ(5));
t265 = cos(qJ(5));
t266 = cos(qJ(4));
t331 = t265 * t266;
t188 = t260 * t223 + t224 * t331;
t259 = sin(qJ(6));
t261 = sin(qJ(4));
t264 = cos(qJ(6));
t334 = t261 * t264;
t172 = t259 * t188 - t224 * t334;
t320 = MDP(35) * t172;
t339 = t259 * t261;
t173 = t264 * t188 + t224 * t339;
t321 = MDP(34) * t173;
t359 = -t320 + t321;
t358 = -0.2e1 * t260 * t265;
t311 = t261 * MDP(28);
t323 = MDP(27) * t261;
t357 = t260 * t311 - t265 * t323;
t356 = t261 * MDP(20) + t266 * MDP(21);
t333 = t261 * t265;
t216 = t259 * t333 + t264 * t266;
t217 = -t259 * t266 + t264 * t333;
t337 = t260 * t261;
t299 = t217 * MDP(34) - t216 * MDP(35) + MDP(36) * t337;
t281 = t264 * MDP(37) - t259 * MDP(38);
t271 = t264 * MDP(34) - t259 * MDP(35) - pkin(6) * t281;
t355 = -MDP(26) + t271;
t280 = MDP(37) * t259 + MDP(38) * t264;
t272 = t259 * MDP(34) + t264 * MDP(35) - pkin(6) * t280;
t354 = -MDP(28) + t272;
t248 = -t267 * pkin(2) - pkin(1);
t353 = 0.2e1 * t248;
t352 = 0.2e1 * t267;
t351 = 2 * MDP(30);
t350 = 2 * MDP(31);
t349 = -2 * MDP(33);
t348 = 0.2e1 * MDP(37);
t347 = 0.2e1 * MDP(38);
t345 = pkin(6) * t265;
t344 = t260 * pkin(3);
t247 = -t346 * pkin(2) - pkin(3);
t343 = pkin(3) - t247;
t342 = t188 * t266;
t255 = t261 ^ 2;
t341 = t255 * t260;
t340 = t255 * t265;
t338 = t260 * t247;
t335 = t260 * t266;
t332 = t261 * t266;
t192 = t223 * pkin(3) - t224 * pkin(5) + t248;
t330 = t266 * t192;
t169 = -t188 * pkin(6) - t330;
t174 = (pkin(6) * t224 + t192 * t265) * t261;
t163 = t264 * t169 - t259 * t174;
t303 = t192 * t337;
t329 = t163 * t337 + t216 * t303;
t246 = t262 * pkin(2) + pkin(5);
t198 = t338 + (t246 * t265 - pkin(6)) * t266;
t213 = (t246 - t345) * t261;
t178 = -t259 * t198 + t264 * t213;
t203 = t246 * t335 - t265 * t247;
t328 = t178 * t337 + t203 * t216;
t214 = -t344 + (pkin(5) * t265 - pkin(6)) * t266;
t227 = (pkin(5) - t345) * t261;
t189 = -t259 * t214 + t264 * t227;
t225 = t265 * pkin(3) + pkin(5) * t335;
t327 = t189 * t337 + t225 * t216;
t204 = t246 * t331 + t338;
t326 = t204 * t266 + t246 * t340;
t226 = pkin(5) * t331 - t344;
t325 = pkin(5) * t340 + t226 * t266;
t324 = MDP(26) * t261;
t322 = MDP(30) * t266;
t187 = -t265 * t223 + t224 * t335;
t319 = MDP(36) * t187;
t318 = MDP(38) * t260;
t317 = t173 * MDP(32);
t316 = t188 * MDP(25);
t315 = t188 * MDP(26);
t314 = t223 * MDP(22);
t313 = t259 * MDP(32);
t310 = t264 * MDP(32);
t308 = t265 * MDP(36);
t307 = t266 * MDP(23);
t306 = t266 * MDP(24);
t305 = t266 * MDP(29);
t304 = pkin(6) * t337;
t254 = t260 ^ 2;
t302 = t254 * t339;
t301 = t254 * t334;
t300 = t224 * t340;
t298 = MDP(25) * t333;
t297 = MDP(34) * t337;
t296 = MDP(35) * t337;
t295 = t259 * t264 * MDP(33);
t293 = t224 * t323;
t292 = t224 * t311;
t291 = MDP(19) * t332;
t290 = t259 * t304;
t289 = t264 * t304;
t288 = pkin(6) * t301;
t287 = -pkin(3) * t224 - pkin(5) * t223;
t286 = -t223 * t246 + t224 * t247;
t285 = t266 * MDP(20) - t261 * MDP(21);
t284 = -t261 * MDP(24) + t307;
t164 = t259 * t169 + t264 * t174;
t282 = MDP(37) * t163 - MDP(38) * t164;
t279 = -t266 * MDP(27) - t261 * t308;
t278 = -t164 * t318 - t224 * t305;
t277 = 0.2e1 * t284;
t276 = (t346 * MDP(16) - t262 * MDP(17)) * pkin(2);
t275 = -t306 + (-MDP(30) * t265 + MDP(31) * t260 - MDP(23)) * t261;
t274 = (-t259 * t216 + t217 * t264) * MDP(33) + t217 * t313 + t259 * t297 + t264 * t296 - t305 - t357;
t257 = t265 ^ 2;
t273 = -MDP(28) * t331 + (t216 * t265 - t302) * MDP(35) + (-t217 * t265 + t301) * MDP(34) + (-t254 + t257) * t324 + ((-t216 * t264 - t217 * t259) * MDP(33) + t217 * t310 + t298) * t260 + t356;
t258 = t266 ^ 2;
t270 = t258 * MDP(29) - 0.2e1 * t216 * t296 + MDP(15) + 0.2e1 * t291 + 0.2e1 * t357 * t266 + (MDP(32) * t217 + t216 * t349 + 0.2e1 * t297) * t217 + (MDP(25) * t257 + MDP(26) * t358 + MDP(36) * t254 + MDP(18)) * t255;
t269 = t282 + t319 + t359;
t268 = (-t172 * t217 - t173 * t216) * MDP(33) + t217 * t317 + (t300 - t342) * MDP(27) + t359 * t337 + (-t260 * t324 + t298) * t188 + (t266 * MDP(28) - t265 * t324 + t299) * t187 + (-t341 * MDP(28) + (-t255 + t258) * MDP(19) + MDP(18) * t332 + MDP(13)) * t224 + (-MDP(14) + t356) * t223;
t256 = t264 ^ 2;
t253 = t259 ^ 2;
t243 = pkin(5) * t341;
t236 = pkin(6) * t302;
t228 = t246 * t341;
t197 = t225 * t217;
t190 = t264 * t214 + t259 * t227;
t186 = t203 * t217;
t179 = t264 * t198 + t259 * t213;
t177 = t217 * t303;
t1 = [t224 * MDP(17) * t353 + pkin(1) * MDP(9) * t352 + MDP(1) + (0.2e1 * t293 + t316) * t188 + (t172 * t349 + t317) * t173 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t263 + MDP(5) * t352) * t263 + (t258 * MDP(18) + t255 * MDP(29) + MDP(11) - 0.2e1 * t291) * t224 ^ 2 + (-0.2e1 * t292 - 0.2e1 * t315 + t319 - 0.2e1 * t320 + 0.2e1 * t321) * t187 + 0.2e1 * t282 * t187 + 0.2e1 * ((-t300 - t342) * MDP(31) - t187 * t322 + (-t255 * t224 * MDP(30) + (t172 * MDP(37) + t173 * MDP(38)) * t261) * t260) * t192 + (MDP(16) * t353 + t314 + t192 * t277 + 0.2e1 * (-MDP(12) + t285) * t224) * t223; (t203 * t173 - t179 * t187 + t177) * MDP(38) + (t286 * MDP(23) + (t187 * t246 - t203 * t224) * MDP(30) + (t188 * t246 - t204 * t224) * MDP(31) + t278) * t261 + (t172 * t203 + t178 * t187 + t329) * MDP(37) + t263 * MDP(6) + t267 * MDP(7) + t268 + t286 * t306; (t203 * t266 + t228) * t351 + t326 * t350 + t328 * t348 + (-t179 * t337 + t186) * t347 + t270 + MDP(8) - 0.2e1 * t284 * t247 + 0.2e1 * t276; (t172 * t225 + t189 * t187 + t329) * MDP(37) + t268 + (t225 * t173 - t190 * t187 + t177) * MDP(38) + (t287 * MDP(23) + (pkin(5) * t187 - t224 * t225) * MDP(30) + (pkin(5) * t188 - t224 * t226) * MDP(31) + t278) * t261 + t287 * t306; (-t343 * MDP(24) + (-t179 - t190) * t318) * t261 + (t327 + t328) * MDP(37) + (t325 + t326) * MDP(31) + (t228 + t243) * MDP(30) + (t343 * MDP(23) + (t203 + t225) * MDP(30)) * t266 + (t186 + t197) * MDP(38) + t270 + t276; pkin(3) * t277 + t270 + (t225 * t266 + t243) * t351 + t325 * t350 + t327 * t348 + (-t190 * t337 + t197) * t347; t314 + t285 * t224 + (t307 + (t254 * t280 - MDP(24)) * t261) * t192 + (t192 * t322 - t269 + t292 + t315) * t265 + (t316 + t293 - MDP(31) * t330 + t173 * t310 + (-t172 * t264 - t173 * t259) * MDP(33) + t355 * t187) * t260; (-t178 * t265 - t288) * MDP(37) + (t179 * t265 + t236) * MDP(38) + (t203 * t280 + t279) * t260 + t275 * t246 + t273; (-t189 * t265 - t288) * MDP(37) + (t190 * t265 + t236) * MDP(38) + (t225 * t280 + t279) * t260 + t275 * pkin(5) + t273; t257 * MDP(36) + MDP(22) + (MDP(32) * t256 + MDP(25) - 0.2e1 * t295) * t254 + t355 * t358; t188 * MDP(27) + t173 * t313 + (-t259 * t172 + t173 * t264) * MDP(33) + (t224 * MDP(29) + (-MDP(31) * t265 + (-MDP(30) - t281) * t260) * t192) * t261 + t354 * t187; -t203 * MDP(30) - t204 * MDP(31) + (-t203 * t264 - t290) * MDP(37) + (t203 * t259 - t289) * MDP(38) + t274; -t225 * MDP(30) - t226 * MDP(31) + (-t225 * t264 - t290) * MDP(37) + (t225 * t259 - t289) * MDP(38) + t274; (MDP(27) + t259 * t310 + (-t253 + t256) * MDP(33)) * t260 - t354 * t265; t253 * MDP(32) + MDP(29) + 0.2e1 * t295; t269; t178 * MDP(37) - t179 * MDP(38) + t299; t189 * MDP(37) - t190 * MDP(38) + t299; t260 * t271 - t308; t272; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
