% Calculate joint inertia matrix for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR6_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:40
% EndTime: 2019-03-10 04:21:44
% DurationCPUTime: 2.15s
% Computational Cost: add. (2768->316), mult. (6159->430), div. (0->0), fcn. (7170->12), ass. (0->159)
t305 = sin(qJ(5));
t310 = cos(qJ(5));
t304 = sin(qJ(6));
t309 = cos(qJ(6));
t274 = t304 * t305 - t309 * t310;
t276 = t304 * t310 + t305 * t309;
t357 = MDP(34) * t276 - MDP(35) * t274;
t398 = MDP(27) * t305 + MDP(28) * t310 + t357;
t302 = sin(pkin(6));
t312 = cos(qJ(2));
t360 = t302 * t312;
t303 = cos(pkin(6));
t307 = sin(qJ(3));
t311 = cos(qJ(3));
t308 = sin(qJ(2));
t361 = t302 * t308;
t261 = t303 * t307 + t311 * t361;
t306 = sin(qJ(4));
t284 = t307 * t361;
t334 = t303 * t311 - t284;
t379 = cos(qJ(4));
t230 = t261 * t379 + t306 * t334;
t216 = t230 * t305 + t310 * t360;
t217 = t230 * t310 - t305 * t360;
t197 = t216 * t309 + t217 * t304;
t198 = -t216 * t304 + t217 * t309;
t395 = MDP(34) * t198 - MDP(35) * t197;
t277 = t306 * t311 + t307 * t379;
t231 = t276 * t277;
t226 = t231 * MDP(35);
t232 = t274 * t277;
t227 = t232 * MDP(34);
t394 = -t227 - t226;
t290 = pkin(3) * t306 + pkin(11);
t270 = (-pkin(12) - t290) * t305;
t297 = t310 * pkin(12);
t362 = t290 * t310;
t271 = t297 + t362;
t237 = t270 * t309 - t271 * t304;
t238 = t270 * t304 + t271 * t309;
t393 = MDP(37) * t237 - MDP(38) * t238;
t280 = (-pkin(11) - pkin(12)) * t305;
t373 = pkin(11) * t310;
t282 = t297 + t373;
t249 = t280 * t309 - t282 * t304;
t251 = t280 * t304 + t282 * t309;
t392 = MDP(37) * t249 - MDP(38) * t251;
t340 = pkin(8) * t360;
t378 = pkin(1) * t308;
t256 = t340 + (pkin(9) + t378) * t303;
t257 = (-pkin(2) * t312 - pkin(9) * t308 - pkin(1)) * t302;
t218 = -t256 * t307 + t257 * t311;
t341 = pkin(3) * t360;
t204 = -pkin(10) * t261 + t218 - t341;
t219 = t311 * t256 + t307 * t257;
t208 = pkin(10) * t334 + t219;
t336 = t379 * t208;
t190 = t204 * t306 + t336;
t188 = -pkin(11) * t360 + t190;
t229 = t261 * t306 - t334 * t379;
t286 = pkin(8) * t361;
t293 = -pkin(3) * t311 - pkin(2);
t377 = pkin(1) * t312;
t236 = t284 * pkin(3) + t286 + (t293 - t377) * t303;
t193 = t229 * pkin(4) - t230 * pkin(11) + t236;
t178 = -t188 * t305 + t193 * t310;
t179 = t188 * t310 + t193 * t305;
t391 = t178 * MDP(30) - t179 * MDP(31);
t390 = MDP(37) * t274 + MDP(38) * t276;
t324 = MDP(30) * t310 - MDP(31) * t305;
t389 = -(MDP(16) * t307 + MDP(17) * t311) * pkin(9) + t307 * MDP(13) + t311 * MDP(14);
t380 = pkin(9) + pkin(10);
t281 = t380 * t307;
t283 = t380 * t311;
t250 = t281 * t379 + t283 * t306;
t252 = -t281 * t306 + t283 * t379;
t275 = t306 * t307 - t311 * t379;
t388 = t277 * MDP(20) - t275 * MDP(21) - t250 * MDP(23) - t252 * MDP(24);
t387 = 2 * MDP(12);
t386 = 0.2e1 * MDP(16);
t385 = 0.2e1 * MDP(30);
t384 = 0.2e1 * MDP(31);
t383 = -2 * MDP(33);
t382 = 0.2e1 * MDP(37);
t381 = 0.2e1 * MDP(38);
t376 = pkin(5) * t229;
t375 = pkin(5) * t275;
t374 = pkin(5) * t309;
t372 = t310 * pkin(5);
t177 = -pkin(12) * t216 + t179;
t370 = t177 * t309;
t358 = -t204 * t379 + t208 * t306;
t187 = pkin(4) * t360 + t358;
t369 = t187 * t310;
t240 = pkin(4) * t275 - pkin(11) * t277 + t293;
t365 = t252 * t310;
t202 = t365 + (-pkin(12) * t277 + t240) * t305;
t368 = t202 * t309;
t367 = t229 * t305;
t366 = t250 * t310;
t364 = t277 * t305;
t363 = t277 * t310;
t359 = t303 * t308;
t355 = MDP(24) * t293;
t348 = t198 * MDP(32);
t347 = t217 * MDP(25);
t225 = t230 * MDP(20);
t346 = t232 * MDP(32);
t345 = t293 * MDP(23);
t344 = t307 * MDP(11);
t343 = t310 * MDP(25);
t342 = MDP(29) + MDP(36);
t339 = t379 * pkin(3);
t338 = MDP(36) * t229 + t395;
t337 = MDP(36) * t275 + t394;
t335 = t305 * t310 * MDP(26);
t176 = -pkin(12) * t217 + t178 + t376;
t173 = t176 * t309 - t177 * t304;
t206 = t240 * t310 - t252 * t305;
t201 = -pkin(12) * t363 + t206 + t375;
t184 = t201 * t309 - t202 * t304;
t291 = -t339 - pkin(4);
t332 = -pkin(4) * t277 - pkin(11) * t275;
t299 = t305 ^ 2;
t331 = t299 * MDP(25) + MDP(22) + 0.2e1 * t335 + (MDP(32) * t276 + t274 * t383) * t276;
t330 = -t275 * t290 + t277 * t291;
t327 = t217 * MDP(27) - t216 * MDP(28);
t326 = MDP(27) * t310 - MDP(28) * t305;
t207 = t240 * t305 + t365;
t325 = t206 * MDP(30) - t207 * MDP(31);
t323 = MDP(30) * t305 + MDP(31) * t310;
t174 = t176 * t304 + t370;
t322 = t173 * MDP(37) - MDP(38) * t174;
t185 = t201 * t304 + t368;
t321 = MDP(37) * t184 - MDP(38) * t185;
t320 = -MDP(19) + t326;
t319 = (MDP(37) * t309 - MDP(38) * t304) * pkin(5);
t318 = 0.2e1 * t390;
t317 = (MDP(23) * t379 - MDP(24) * t306) * pkin(3);
t316 = -t261 * MDP(13) - MDP(14) * t334;
t300 = t310 ^ 2;
t315 = (-t231 * t276 + t232 * t274) * MDP(33) - t276 * t346 + (-t299 + t300) * t277 * MDP(26) + t343 * t364 + t388 + t398 * t275;
t314 = -MDP(19) * t230 + t327 + t395;
t313 = -MDP(22) * t360 + (-t197 * t276 - t198 * t274) * MDP(33) + t276 * t348 + (-t216 * t305 + t217 * t310) * MDP(26) + t305 * t347 + t225 + (-MDP(21) + t398) * t229;
t298 = t302 ^ 2;
t292 = -pkin(4) - t372;
t279 = t291 - t372;
t263 = pkin(1) * t359 + t340;
t262 = t303 * t377 - t286;
t255 = t286 + (-pkin(2) - t377) * t303;
t241 = t250 * t305;
t222 = pkin(5) * t364 + t250;
t211 = t222 * t276;
t210 = t222 * t274;
t186 = t187 * t305;
t183 = pkin(5) * t216 + t187;
t182 = t183 * t276;
t181 = t183 * t274;
t1 = [t230 ^ 2 * MDP(18) + t303 ^ 2 * MDP(8) + MDP(1) + (t261 * MDP(11) + t334 * t387) * t261 + t342 * t229 ^ 2 + (-0.2e1 * MDP(26) * t216 + t347) * t217 + (t197 * t383 + t348) * t198 + ((MDP(4) * t308 + 0.2e1 * MDP(5) * t312) * t308 + (MDP(15) + MDP(22)) * t312 ^ 2) * t298 + 0.2e1 * MDP(6) * t359 * t302 + (-t218 * t360 - t255 * t334) * t386 + 0.2e1 * (t219 * t360 + t255 * t261) * MDP(17) + 0.2e1 * (t229 * t236 + t358 * t360) * MDP(23) + 0.2e1 * (t190 * t360 + t230 * t236) * MDP(24) + 0.2e1 * (t262 * t303 + t298 * t377) * MDP(9) + 0.2e1 * (-t263 * t303 - t298 * t378) * MDP(10) + (t178 * t229 + t187 * t216) * t385 + (t173 * t229 + t183 * t197) * t382 + (-t179 * t229 + t187 * t217) * t384 + (-t174 * t229 + t183 * t198) * t381 + 0.2e1 * (t303 * MDP(7) - t225 + t316) * t360 + 0.2e1 * (MDP(21) * t360 + t314) * t229; (-pkin(2) * t261 + t255 * t307) * MDP(17) + (pkin(2) * t334 - t255 * t311) * MDP(16) + t230 * t355 + t261 * t344 - t198 * t346 + (t197 * t232 - t198 * t231) * MDP(33) + t262 * MDP(9) - t263 * MDP(10) + (t183 * t231 + t197 * t222) * MDP(37) + (-t183 * t232 + t198 * t222) * MDP(38) + t303 * MDP(8) + (t261 * t311 + t307 * t334) * MDP(12) + (MDP(30) * t216 + MDP(31) * t217) * t250 + (t321 + t325 + t345 + t394) * t229 + (t236 * MDP(23) + t342 * t229 + t314 + t322 + t391) * t275 + (t217 * t343 + t236 * MDP(24) + t230 * MDP(18) + (-t216 * t310 - t217 * t305) * MDP(26) + t323 * t187 + t320 * t229) * t277 + (MDP(6) * t308 + (MDP(7) - t388 - t389) * t312) * t302; pkin(2) * t311 * t386 + MDP(8) - (t231 * t383 - t346) * t232 + (-0.2e1 * MDP(17) * pkin(2) + t311 * t387 + t344) * t307 + (0.2e1 * t355 + (MDP(25) * t300 + MDP(18) - 0.2e1 * t335) * t277) * t277 + (t363 * t384 + t364 * t385) * t250 + (t231 * t382 - t232 * t381) * t222 + (t184 * t382 - t185 * t381 + t206 * t385 - t207 * t384 + t275 * t342 + 0.2e1 * t277 * t320 - 0.2e1 * t226 - 0.2e1 * t227 + 0.2e1 * t345) * t275; t313 - MDP(15) * t360 + t218 * MDP(16) - t219 * MDP(17) + (t197 * t279 + t229 * t237 + t181) * MDP(37) + (t198 * t279 - t229 * t238 + t182) * MDP(38) + (t216 * t291 - t290 * t367 - t369) * MDP(30) + (t217 * t291 - t229 * t362 + t186) * MDP(31) + (-t339 * t360 - t358) * MDP(23) + (-t336 + (-t204 + t341) * t306) * MDP(24) - t316; t315 + (t310 * t330 + t241) * MDP(31) + (t305 * t330 - t366) * MDP(30) + (t231 * t279 + t237 * t275 + t210) * MDP(37) + (-t232 * t279 - t238 * t275 + t211) * MDP(38) + t389; t279 * t318 - 0.2e1 * t291 * t324 + MDP(15) + 0.2e1 * t317 + t331; t313 - t358 * MDP(23) - t190 * MDP(24) + (t197 * t292 + t229 * t249 + t181) * MDP(37) + (t198 * t292 - t229 * t251 + t182) * MDP(38) + (-pkin(4) * t216 - pkin(11) * t367 - t369) * MDP(30) + (-pkin(4) * t217 - t229 * t373 + t186) * MDP(31); t315 + (t310 * t332 + t241) * MDP(31) + (t305 * t332 - t366) * MDP(30) + (t231 * t292 + t249 * t275 + t210) * MDP(37) + (-t232 * t292 - t251 * t275 + t211) * MDP(38); t317 + t331 + t324 * (pkin(4) - t291) + t390 * (t279 + t292); 0.2e1 * pkin(4) * t324 + t292 * t318 + t331; t229 * MDP(29) + (t229 * t374 + t173) * MDP(37) + (-t370 + (-t176 - t376) * t304) * MDP(38) + t327 + t338 + t391; t275 * MDP(29) + (t275 * t374 + t184) * MDP(37) + (-t368 + (-t201 - t375) * t304) * MDP(38) + t326 * t277 + t325 + t337; -t290 * t323 + t393 + t398; -pkin(11) * t323 + t392 + t398; 0.2e1 * t319 + t342; t322 + t338; t321 + t337; t357 + t393; t357 + t392; MDP(36) + t319; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
