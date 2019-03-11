% Calculate joint inertia matrix for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR8_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:52
% EndTime: 2019-03-10 05:07:01
% DurationCPUTime: 3.12s
% Computational Cost: add. (4972->386), mult. (12842->546), div. (0->0), fcn. (14965->14), ass. (0->181)
t298 = sin(qJ(4));
t302 = cos(qJ(4));
t407 = (MDP(23) * t298 + MDP(24) * t302) * pkin(11) - t298 * MDP(20) - t302 * MDP(21);
t296 = sin(qJ(6));
t301 = cos(qJ(6));
t367 = t296 * MDP(34) + t301 * MDP(35);
t388 = pkin(11) + pkin(12);
t275 = t388 * t302;
t297 = sin(qJ(5));
t386 = cos(qJ(5));
t333 = t386 * t298;
t248 = t275 * t297 + t388 * t333;
t368 = t297 * t298;
t249 = t386 * t275 - t388 * t368;
t270 = -t386 * t302 + t368;
t271 = t297 * t302 + t333;
t399 = -t271 * MDP(27) + t270 * MDP(28) + t248 * MDP(30) + t249 * MDP(31);
t406 = MDP(14) + t399 + t407;
t295 = cos(pkin(6));
t299 = sin(qJ(3));
t303 = cos(qJ(3));
t293 = sin(pkin(6));
t294 = cos(pkin(7));
t387 = cos(qJ(2));
t336 = t294 * t387;
t330 = t293 * t336;
t300 = sin(qJ(2));
t370 = t293 * t300;
t292 = sin(pkin(7));
t371 = t292 * t303;
t241 = -t295 * t371 + t299 * t370 - t303 * t330;
t342 = MDP(22) + MDP(29);
t405 = t342 * t241;
t372 = t292 * t299;
t242 = t295 * t372 + (t299 * t336 + t300 * t303) * t293;
t337 = t293 * t387;
t259 = t292 * t337 - t295 * t294;
t221 = t242 * t298 + t259 * t302;
t222 = t242 * t302 - t259 * t298;
t208 = t386 * t221 + t222 * t297;
t209 = -t297 * t221 + t386 * t222;
t404 = t209 * MDP(27) - t208 * MDP(28);
t261 = -t294 * t302 + t298 * t372;
t262 = t294 * t298 + t302 * t372;
t231 = t386 * t261 + t262 * t297;
t232 = -t297 * t261 + t386 * t262;
t403 = t232 * MDP(27) - t231 * MDP(28);
t340 = pkin(10) * t371;
t385 = pkin(2) * t299;
t254 = t340 + (pkin(11) + t385) * t294;
t255 = (-pkin(3) * t303 - pkin(11) * t299 - pkin(2)) * t292;
t225 = -t298 * t254 + t302 * t255;
t226 = t254 * t302 + t255 * t298;
t402 = t225 * MDP(23) - t226 * MDP(24);
t339 = pkin(1) * t387;
t278 = t295 * t339;
t243 = pkin(2) * t295 + t278 + (-pkin(10) * t294 - pkin(9)) * t370;
t251 = (-pkin(10) * t292 * t300 - t387 * pkin(2) - pkin(1)) * t293;
t218 = -t243 * t292 + t294 * t251;
t200 = pkin(3) * t241 - pkin(11) * t242 + t218;
t266 = t295 * t300 * pkin(1) + pkin(9) * t337;
t239 = (t292 * t295 + t330) * pkin(10) + t266;
t328 = t243 * t294 + t251 * t292;
t212 = t239 * t303 + t328 * t299;
t204 = -t259 * pkin(11) + t212;
t188 = t302 * t200 - t204 * t298;
t189 = t200 * t298 + t204 * t302;
t401 = t188 * MDP(23) - t189 * MDP(24);
t318 = MDP(37) * t301 - MDP(38) * t296;
t400 = t259 * MDP(15);
t211 = -t299 * t239 + t328 * t303;
t398 = 2 * MDP(16);
t397 = 2 * MDP(17);
t396 = -2 * MDP(19);
t395 = 0.2e1 * MDP(23);
t394 = 0.2e1 * MDP(24);
t393 = 0.2e1 * MDP(30);
t392 = 0.2e1 * MDP(31);
t391 = -2 * MDP(33);
t390 = 0.2e1 * MDP(37);
t389 = 0.2e1 * MDP(38);
t384 = pkin(2) * t303;
t383 = pkin(4) * t241;
t382 = pkin(13) * t296;
t381 = pkin(13) * t301;
t183 = -pkin(12) * t222 + t188 + t383;
t185 = -pkin(12) * t221 + t189;
t179 = t386 * t183 - t297 * t185;
t177 = -pkin(5) * t241 - t179;
t379 = t177 * t301;
t341 = pkin(4) * t371;
t214 = -t262 * pkin(12) + t225 - t341;
t217 = -pkin(12) * t261 + t226;
t197 = t386 * t214 - t297 * t217;
t195 = pkin(5) * t371 - t197;
t378 = t195 * t301;
t377 = t248 * t301;
t376 = t271 * t296;
t281 = pkin(4) * t297 + pkin(13);
t375 = t281 * t296;
t374 = t281 * t301;
t288 = t293 ^ 2;
t373 = t288 * t300;
t369 = t295 * MDP(8);
t366 = MDP(11) * t242;
t365 = MDP(13) * t259;
t364 = MDP(18) * t262;
t363 = MDP(18) * t298;
t362 = MDP(26) * t209;
t361 = MDP(26) * t232;
t224 = t301 * t232 - t296 * t371;
t360 = MDP(32) * t224;
t359 = MDP(32) * t296;
t358 = MDP(32) * t301;
t193 = t209 * t301 + t241 * t296;
t357 = MDP(34) * t193;
t356 = MDP(34) * t224;
t192 = t209 * t296 - t241 * t301;
t355 = MDP(35) * t192;
t223 = t296 * t232 + t301 * t371;
t354 = MDP(35) * t223;
t353 = MDP(36) * t208;
t352 = MDP(36) * t231;
t345 = t232 * MDP(25);
t344 = t270 * MDP(36);
t343 = t299 * MDP(13);
t338 = t386 * pkin(4);
t283 = -pkin(4) * t302 - pkin(3);
t335 = t386 * t185;
t334 = t386 * t217;
t332 = t296 * t301 * MDP(33);
t289 = t296 ^ 2;
t331 = t289 * MDP(32) + MDP(29) + 0.2e1 * t332;
t329 = -pkin(5) * t271 - pkin(13) * t270;
t282 = -t338 - pkin(5);
t327 = -t270 * t281 + t271 * t282;
t277 = pkin(10) * t372;
t263 = t294 * t384 - t277;
t265 = t294 * t385 + t340;
t326 = -t263 * MDP(16) + t265 * MDP(17);
t325 = t222 * MDP(20) - t221 * MDP(21);
t324 = t262 * MDP(20) - t261 * MDP(21);
t180 = t297 * t183 + t335;
t321 = t179 * MDP(30) - t180 * MDP(31);
t198 = t297 * t214 + t334;
t320 = t197 * MDP(30) - t198 * MDP(31);
t319 = MDP(34) * t301 - MDP(35) * t296;
t317 = MDP(37) * t296 + MDP(38) * t301;
t253 = t277 + (-pkin(3) - t384) * t294;
t316 = -MDP(26) + t319;
t314 = (-t192 * t296 + t193 * t301) * MDP(33) + t193 * t359 + t241 * MDP(29) + t404 + t367 * t208;
t313 = (t300 * MDP(6) + t387 * MDP(7)) * t293;
t312 = (t386 * MDP(30) - t297 * MDP(31)) * pkin(4);
t290 = t301 ^ 2;
t311 = (-t289 + t290) * t271 * MDP(33) + t358 * t376 - t399 + t367 * t270;
t233 = t261 * pkin(4) + t253;
t310 = -MDP(29) * t371 + (-t223 * t296 + t224 * t301) * MDP(33) + t224 * t359 + t403 + t367 * t231;
t203 = t259 * pkin(3) - t211;
t309 = t242 * MDP(13) + t211 * MDP(16) - t212 * MDP(17) - t400;
t190 = t221 * pkin(4) + t203;
t308 = MDP(14) * t294 - t324 - t403;
t178 = t241 * pkin(13) + t180;
t181 = t208 * pkin(5) - t209 * pkin(13) + t190;
t174 = -t178 * t296 + t181 * t301;
t175 = t178 * t301 + t181 * t296;
t307 = MDP(37) * t174 - MDP(38) * t175 + t353 - t355 + t357;
t196 = -pkin(13) * t371 + t198;
t207 = t231 * pkin(5) - t232 * pkin(13) + t233;
t186 = -t196 * t296 + t207 * t301;
t187 = t196 * t301 + t207 * t296;
t306 = MDP(37) * t186 - MDP(38) * t187 + t352 - t354 + t356;
t305 = t242 * MDP(12) - t259 * MDP(14) - t325 - t404;
t287 = t292 ^ 2;
t264 = -pkin(9) * t370 + t278;
t244 = t248 * t296;
t236 = pkin(5) * t270 - pkin(13) * t271 + t283;
t216 = t236 * t296 + t249 * t301;
t215 = t236 * t301 - t249 * t296;
t194 = t195 * t296;
t176 = t177 * t296;
t1 = [t209 ^ 2 * MDP(25) + MDP(1) + (MDP(4) * t300 + 0.2e1 * MDP(5) * t387) * t373 + (MDP(18) * t222 + t221 * t396) * t222 + (MDP(32) * t193 + t192 * t391) * t193 + (0.2e1 * t313 + t369) * t295 + (t353 - 0.2e1 * t355 + 0.2e1 * t357 - 0.2e1 * t362) * t208 + (t174 * t208 + t177 * t192) * t390 + (-t175 * t208 + t177 * t193) * t389 + 0.2e1 * (t264 * t295 + t288 * t339) * MDP(9) + 0.2e1 * (-pkin(1) * t373 - t266 * t295) * MDP(10) + (t218 * t397 - 0.2e1 * t365 + t366) * t242 + (-t211 * t398 + t212 * t397 + t400) * t259 + (t221 * t395 + t222 * t394) * t203 + (t208 * t393 + t209 * t392) * t190 + (t179 * t393 - t180 * t392 + t188 * t395 - t189 * t394 + t218 * t398 - 0.2e1 * t305 + t405) * t241; (-t221 * t262 - t222 * t261) * MDP(19) + t264 * MDP(9) + (-t192 * t224 - t193 * t223) * MDP(33) + (t174 * t231 + t177 * t223 + t186 * t208 + t192 * t195) * MDP(37) + (-t192 * t231 - t208 * t223) * MDP(35) + (t193 * t231 + t208 * t224) * MDP(34) + (-t175 * t231 + t177 * t224 - t187 * t208 + t193 * t195) * MDP(38) + (-t208 * t232 - t209 * t231) * MDP(26) + (t203 * t261 + t253 * t221) * MDP(23) + (t203 * t262 + t253 * t222) * MDP(24) + (t190 * t231 + t233 * t208) * MDP(30) + (t190 * t232 + t233 * t209) * MDP(31) + t193 * t360 + t208 * t352 + t209 * t345 - t266 * MDP(10) + t369 + t222 * t364 + t326 * t259 + t309 * t294 + t313 + (-t308 + t320 + t402) * t241 + ((-MDP(16) * t241 - MDP(17) * t242) * pkin(2) + (-MDP(12) * t241 + MDP(17) * t218 - t365 + t366) * t299 + (-t218 * MDP(16) + t305 - t321 - t401 - t405) * t303) * t292; t294 ^ 2 * MDP(15) + t232 ^ 2 * MDP(25) + MDP(8) + (t261 * t396 + t364) * t262 + (t223 * t391 + t360) * t224 + ((MDP(11) * t299 + 0.2e1 * MDP(12) * t303) * t299 + t342 * t303 ^ 2) * t287 + (t352 - 0.2e1 * t354 + 0.2e1 * t356 - 0.2e1 * t361) * t231 + 0.2e1 * (t294 * t343 + t303 * t308) * t292 + (t186 * t231 + t195 * t223) * t390 + (-t187 * t231 + t195 * t224) * t389 + (-t265 * t294 - t287 * t385) * t397 + (-t225 * t371 + t253 * t261) * t395 + (-t197 * t371 + t233 * t231) * t393 + (t198 * t371 + t233 * t232) * t392 + (t226 * t371 + t253 * t262) * t394 + (t263 * t294 + t287 * t384) * t398; (-t221 * t298 + t222 * t302) * MDP(19) + (-pkin(3) * t221 - t203 * t302) * MDP(23) + (-pkin(3) * t222 + t203 * t298) * MDP(24) + (t192 * t248 + t208 * t215) * MDP(37) + (t193 * t248 - t208 * t216) * MDP(38) + t222 * t363 + (MDP(30) * t208 + MDP(31) * t209) * t283 + (MDP(30) * t190 + t307 - t362) * t270 + (t190 * MDP(31) + (-t192 * t301 - t193 * t296) * MDP(33) + t193 * t358 + t209 * MDP(25) + t317 * t177 + t316 * t208) * t271 - t406 * t241 + t309; t294 * MDP(15) + (-t261 * t298 + t262 * t302) * MDP(19) + (-pkin(3) * t261 - t253 * t302) * MDP(23) + (-pkin(3) * t262 + t253 * t298) * MDP(24) + (t215 * t231 + t223 * t248) * MDP(37) + (-t216 * t231 + t224 * t248) * MDP(38) + t262 * t363 + (MDP(30) * t231 + MDP(31) * t232) * t283 + (MDP(30) * t233 + t306 - t361) * t270 + (t233 * MDP(31) + (-t223 * t301 - t224 * t296) * MDP(33) + t224 * t358 + t345 + t317 * t195 + t316 * t231) * t271 + (t406 * t303 + t343) * t292 - t326; pkin(3) * t302 * t395 + MDP(15) + t248 * t376 * t390 + (0.2e1 * MDP(19) * t302 - 0.2e1 * MDP(24) * pkin(3) + t363) * t298 + (t215 * t390 - t216 * t389 + t283 * t393 + t344) * t270 + (0.2e1 * t316 * t270 + t283 * t392 + t377 * t389 + (MDP(32) * t290 + MDP(25) - 0.2e1 * t332) * t271) * t271; t241 * MDP(22) + (t241 * t338 + t179) * MDP(30) + (-t335 + (-t183 - t383) * t297) * MDP(31) + (t192 * t282 - t208 * t375 - t379) * MDP(37) + (t193 * t282 - t208 * t374 + t176) * MDP(38) + t314 + t325 + t401; -MDP(22) * t371 + (-t338 * t371 + t197) * MDP(30) + (-t334 + (-t214 + t341) * t297) * MDP(31) + (t223 * t282 - t231 * t375 - t378) * MDP(37) + (t224 * t282 - t231 * t374 + t194) * MDP(38) + t310 + t324 + t402; (t327 * t296 - t377) * MDP(37) + (t327 * t301 + t244) * MDP(38) + t311 - t407; -0.2e1 * t282 * t318 + MDP(22) + 0.2e1 * t312 + t331; (-pkin(5) * t192 - t208 * t382 - t379) * MDP(37) + (-pkin(5) * t193 - t208 * t381 + t176) * MDP(38) + t314 + t321; (-pkin(5) * t223 - t231 * t382 - t378) * MDP(37) + (-pkin(5) * t224 - t231 * t381 + t194) * MDP(38) + t310 + t320; (t329 * t296 - t377) * MDP(37) + (t329 * t301 + t244) * MDP(38) + t311; t312 + t331 + t318 * (pkin(5) - t282); 0.2e1 * pkin(5) * t318 + t331; t307; t306; MDP(37) * t215 - MDP(38) * t216 + t271 * t319 + t344; -t281 * t317 + t367; -pkin(13) * t317 + t367; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
