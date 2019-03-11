% Calculate joint inertia matrix for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% MDP [45x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [7x7]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S7RRRRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1),zeros(45,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [45 1]), ...
  'S7RRRRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [45x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:43:11
% EndTime: 2019-03-10 07:43:26
% DurationCPUTime: 4.92s
% Computational Cost: add. (2563->482), mult. (7190->683), div. (0->0), fcn. (8265->12), ass. (0->195)
t291 = sin(qJ(7));
t297 = cos(qJ(7));
t337 = MDP(44) * t297 - MDP(45) * t291;
t315 = t297 * MDP(41) - t291 * MDP(42) - t337 * pkin(4);
t444 = MDP(33) + t315;
t299 = cos(qJ(5));
t293 = sin(qJ(5));
t298 = cos(qJ(6));
t422 = t293 * t298;
t264 = t291 * t422 - t297 * t299;
t270 = t291 * t299 + t297 * t422;
t364 = pkin(3) * t298 + pkin(4);
t274 = t364 * t299;
t276 = (pkin(4) * t298 + pkin(3)) * t293;
t317 = MDP(44) * (-t274 * t291 - t276 * t297) - MDP(45) * (t274 * t297 - t276 * t291) + MDP(41) * t270 - MDP(42) * t264;
t445 = -MDP(35) * t299 + t317;
t439 = 2 * pkin(2);
t369 = 0.2e1 * t293;
t295 = sin(qJ(3));
t443 = 0.2e1 * t295;
t336 = MDP(44) * t291 + MDP(45) * t297;
t316 = -t291 * MDP(41) - t297 * MDP(42) + t336 * pkin(4);
t310 = MDP(35) + t316;
t294 = sin(qJ(4));
t302 = cos(qJ(2));
t300 = cos(qJ(4));
t296 = sin(qJ(2));
t301 = cos(qJ(3));
t417 = t296 * t301;
t366 = t300 * t417;
t269 = -t294 * t302 + t366;
t419 = t295 * t296;
t258 = t269 * t299 - t293 * t419;
t268 = t294 * t417 + t300 * t302;
t292 = sin(qJ(6));
t239 = t258 * t292 - t268 * t298;
t241 = t258 * t298 + t268 * t292;
t257 = t269 * t293 + t299 * t419;
t442 = (t268 * MDP(30) - t239 * MDP(37) - t241 * MDP(38)) * pkin(3) - t258 * MDP(25) + t257 * MDP(26) - t268 * MDP(27);
t418 = t295 * t300;
t272 = -t293 * t301 - t299 * t418;
t421 = t294 * t295;
t253 = t272 * t292 + t298 * t421;
t256 = t272 * t298 - t292 * t421;
t266 = t293 * t418 - t299 * t301;
t377 = t272 * MDP(25);
t441 = -(MDP(37) * t253 + MDP(38) * t256) * pkin(3) - t266 * MDP(26) - t377;
t412 = MDP(27) * t299;
t440 = (MDP(30) * t299 - MDP(31) * t293) * pkin(3) + MDP(19) + MDP(28) * t293 - t412;
t438 = 0.2e1 * pkin(3);
t437 = -0.2e1 * t295;
t436 = 0.2e1 * MDP(37);
t435 = 0.2e1 * MDP(38);
t434 = -2 * MDP(40);
t433 = 0.2e1 * MDP(44);
t432 = 0.2e1 * MDP(45);
t431 = pkin(2) * t296;
t430 = pkin(2) * t301;
t429 = pkin(3) * t294;
t428 = MDP(38) * pkin(3);
t427 = pkin(2) * MDP(24);
t282 = t293 ^ 2;
t426 = t282 * t298;
t284 = t295 ^ 2;
t285 = t296 ^ 2;
t425 = t284 * t285;
t424 = t292 * t298;
t423 = t293 * t294;
t420 = t294 * t299;
t288 = t299 ^ 2;
t416 = t282 + t288;
t290 = t301 ^ 2;
t415 = -t284 - t290;
t414 = MDP(21) * t301;
t413 = MDP(24) * t300;
t410 = MDP(33) * t256;
t409 = MDP(35) * t266;
t407 = MDP(38) * t298;
t238 = t256 * t297 + t266 * t291;
t406 = MDP(39) * t238;
t367 = t298 * t420;
t271 = t292 * t300 - t367;
t255 = t271 * t297 + t291 * t423;
t405 = MDP(39) * t255;
t404 = MDP(39) * t270;
t403 = MDP(39) * t291;
t402 = MDP(39) * t297;
t401 = MDP(41) * t238;
t237 = t256 * t291 - t266 * t297;
t400 = MDP(42) * t237;
t399 = MDP(43) * t253;
t231 = t241 * t291 + t257 * t297;
t398 = t231 * MDP(42);
t232 = t241 * t297 - t257 * t291;
t397 = t232 * MDP(41);
t396 = t239 * MDP(43);
t395 = t241 * MDP(33);
t394 = t241 * MDP(34);
t252 = t271 * t291 - t297 * t423;
t393 = t252 * MDP(42);
t392 = t255 * MDP(41);
t391 = t256 * MDP(32);
t390 = t256 * MDP(34);
t389 = t257 * MDP(35);
t388 = t257 * MDP(36);
t387 = t258 * MDP(26);
t386 = t258 * MDP(27);
t265 = t292 * t420 + t298 * t300;
t385 = t265 * MDP(43);
t384 = t266 * MDP(36);
t383 = t268 * MDP(28);
t382 = t268 * MDP(29);
t381 = t269 * MDP(19);
t380 = t271 * MDP(32);
t379 = t271 * MDP(33);
t378 = t271 * MDP(34);
t376 = t272 * MDP(26);
t373 = t299 * MDP(34);
t372 = t300 * MDP(20);
t371 = t300 * MDP(21);
t370 = t301 * MDP(13);
t368 = t294 * t430;
t365 = -pkin(2) * t300 - pkin(3);
t363 = MDP(23) * t417;
t362 = MDP(37) * t420;
t361 = MDP(40) * t291 * t297;
t360 = t293 * t299 * MDP(26);
t263 = t365 * t419;
t335 = -pkin(2) * t417 - pkin(3) * t269;
t243 = t263 * t293 - t299 * t335;
t275 = t365 * t301;
t357 = (pkin(3) * t300 + pkin(2)) * t295;
t249 = t275 * t293 - t299 * t357;
t359 = pkin(2) * t294 * t419;
t358 = MDP(23) * t300 + MDP(16);
t356 = -t269 * MDP(18) + t268 * MDP(19);
t355 = -t269 * MDP(20) + t268 * MDP(21);
t354 = t294 * MDP(21) - t372;
t352 = t272 * MDP(27) + t266 * MDP(28);
t245 = t299 * t263 + t293 * t335;
t351 = -t243 * MDP(30) - t245 * MDP(31);
t251 = t299 * t275 + t293 * t357;
t350 = t249 * MDP(30) + t251 * MDP(31);
t349 = -MDP(30) * t257 - MDP(31) * t258;
t347 = -t265 * MDP(35) - t378;
t345 = -t392 + t393;
t260 = -pkin(3) * t420 + pkin(4) * t271;
t262 = t364 * t423;
t242 = -t260 * t297 - t262 * t291;
t244 = -t260 * t291 + t262 * t297;
t341 = MDP(44) * t242 - MDP(45) * t244;
t339 = MDP(44) * t252 + MDP(45) * t255;
t338 = MDP(44) * t264 + MDP(45) * t270;
t334 = MDP(37) * t298 - MDP(38) * t292 + MDP(30);
t333 = MDP(37) + t337;
t247 = t251 * t298 - t292 * t368;
t234 = t245 * t298 - t292 * t359;
t330 = -t257 * MDP(28) + t382 + t386;
t328 = t241 * MDP(32) + t257 * MDP(34) + t243 * MDP(38);
t327 = -t266 * MDP(34) + t249 * MDP(38) + t391;
t326 = -t234 * MDP(38) + t388 + t394;
t325 = -t247 * MDP(38) - t384 + t390;
t324 = -MDP(22) + (MDP(23) * t294 + t413) * pkin(2);
t283 = t294 ^ 2;
t323 = MDP(17) + (MDP(30) * t293 + MDP(31) * t299) * t283;
t281 = t292 ^ 2;
t322 = -t336 * t281 - MDP(31);
t321 = -t350 + t352;
t320 = -t239 * MDP(33) + t328;
t319 = -t253 * MDP(33) + t327;
t318 = -t341 + t345;
t314 = t330 + t351;
t229 = pkin(4) * t241 + t243;
t230 = -pkin(4) * t257 + t234;
t225 = -t229 * t297 - t230 * t291;
t226 = -t229 * t291 + t230 * t297;
t313 = t225 * MDP(44) - t226 * MDP(45) - t396 + t397 - t398;
t235 = pkin(4) * t266 + t247;
t236 = pkin(4) * t256 + t249;
t227 = -t235 * t291 - t236 * t297;
t228 = t235 * t297 - t236 * t291;
t312 = MDP(44) * t227 - MDP(45) * t228 - t399 - t400 + t401;
t311 = -t318 + t385;
t307 = t266 * t407 + (MDP(37) * t266 - MDP(44) * t237 - MDP(45) * t238) * t292;
t246 = -t251 * t292 - t298 * t368;
t306 = t253 * MDP(35) - t246 * MDP(37) - t325 + t376;
t305 = -t243 * MDP(37) + t313 + t389 + t395;
t304 = -MDP(37) * t249 + t312 - t409 + t410;
t233 = -t245 * t292 - t298 * t359;
t303 = t387 + t383 + t239 * MDP(35) - t233 * MDP(37) + (-t257 * t407 - t268 * MDP(31) + (-MDP(37) * t257 - MDP(44) * t231 - MDP(45) * t232) * t292) * pkin(3) - t326;
t289 = t300 ^ 2;
t287 = t298 ^ 2;
t286 = t297 ^ 2;
t280 = t291 ^ 2;
t1 = [t302 ^ 2 * MDP(15) + t269 ^ 2 * MDP(18) + t258 ^ 2 * MDP(25) + t241 ^ 2 * MDP(32) + MDP(1) + (MDP(39) * t232 + t231 * t434) * t232 + (t301 * MDP(12) * t437 + t290 * MDP(11) + t284 * MDP(22) + MDP(4)) * t285 + (-0.2e1 * t383 - 0.2e1 * t387 + t388 + 0.2e1 * t394) * t257 + (-0.2e1 * t389 - 0.2e1 * t395 + t396 - 0.2e1 * t397 + 0.2e1 * t398) * t239 + (t355 * t443 + (MDP(14) * t437 + (2 * MDP(5)) + 0.2e1 * t370 + (-MDP(16) * t301 + MDP(17) * t295) * t439) * t302) * t296 + (t233 * t257 + t239 * t243) * t436 + (-t234 * t257 + t241 * t243) * t435 + (-t225 * t239 + t231 * t233) * t433 + (t226 * t239 + t232 * t233) * t432 + ((-t269 * t417 - t300 * t425) * MDP(24) + (-MDP(23) * t425 + t349 * t419) * t294) * t439 + (-t363 * t439 + 0.2e1 * t351 - 0.2e1 * t381 + t382 + 0.2e1 * t386) * t268; (-t257 * t272 + t258 * t266) * MDP(26) + (t234 * t266 + t241 * t249 + t243 * t256 - t247 * t257) * MDP(38) + (-t241 * t266 + t256 * t257) * MDP(34) + (-t233 * t266 + t239 * t249 + t243 * t253 + t246 * t257) * MDP(37) + (t239 * t266 - t253 * t257) * MDP(35) + (-t239 * t256 - t241 * t253) * MDP(33) + (t226 * t253 + t228 * t239 + t232 * t246 + t233 * t238) * MDP(45) + (-t225 * t253 - t227 * t239 + t231 * t246 + t233 * t237) * MDP(44) + (t231 * t253 + t237 * t239) * MDP(42) + (-t232 * t253 - t238 * t239) * MDP(41) + t302 * MDP(7) + (-t231 * t238 - t232 * t237) * MDP(40) + t258 * t377 - t257 * t384 + t253 * t396 + t241 * t391 + t232 * t406 + t321 * t268 + (-MDP(14) * t302 + (MDP(17) * t302 + t349 * t294) * pkin(2) + t355) * t301 + (-t290 * MDP(12) + MDP(6) + (MDP(12) - t354) * t284) * t296 + (-t302 * MDP(13) + t356 * t300 + (-MDP(11) + MDP(22)) * t417 + (t302 * MDP(16) + (t269 - t366) * MDP(24) + t268 * MDP(23)) * pkin(2) + ((-t272 * t431 + t245) * MDP(31) + (t266 * t431 + t243) * MDP(30) - pkin(2) * t363 + t381 - t330) * t294) * t295; t272 ^ 2 * MDP(25) + t256 ^ 2 * MDP(32) + t290 * MDP(22) + MDP(8) + (t237 * t434 + t406) * t238 + (0.2e1 * t376 + t384 - 0.2e1 * t390) * t266 + (-0.2e1 * t294 * t300 * MDP(19) + MDP(18) * t289 + MDP(29) * t283 + MDP(11)) * t284 + ((MDP(12) + t372) * t301 + (-t352 - t414) * t294) * t443 + (t399 + 0.2e1 * t400 - 0.2e1 * t401 + 0.2e1 * t409 - 0.2e1 * t410) * t253 + (t247 * t266 + t249 * t256) * t435 + (-t246 * t266 + t249 * t253) * t436 + (-t227 * t253 + t237 * t246) * t433 + (t228 * t253 + t238 * t246) * t432 + t415 * t413 * t439 + 0.2e1 * (t350 * t295 + (t415 * MDP(23) + (MDP(30) * t266 - MDP(31) * t272) * t301) * pkin(2)) * t294; (-t231 * t255 - t232 * t252) * MDP(40) + t302 * MDP(15) + t232 * t405 + t339 * t233 + t318 * t239 + t320 * t271 + t305 * t265 + (t314 - t381) * t300 + (t370 + (-MDP(14) + t371) * t295 + (t323 * t295 - t358 * t301) * pkin(2)) * t296 + ((MDP(20) * t295 + t301 * t427) * t296 + t442 * t299 + t303 * t293 + t356) * t294; (-t237 * t255 - t238 * t252) * MDP(40) + t238 * t405 + t339 * t246 + t318 * t253 + (t323 * pkin(2) - MDP(14)) * t301 + (t321 + t414) * t300 + t319 * t271 + t304 * t265 + (t289 * MDP(19) + t358 * pkin(2) - t440 * t283 - MDP(13)) * t295 + (t301 * MDP(20) + (-t427 + (MDP(18) - MDP(29)) * t300) * t295 + t441 * t299 + (t307 * pkin(3) + t306) * t293) * t294; t271 ^ 2 * MDP(32) + t289 * MDP(29) + MDP(15) + (t252 * t434 + t405) * t255 + (MDP(25) * t288 + MDP(36) * t282 + MDP(18) - 0.2e1 * t360) * t283 + (0.2e1 * t300 * t440 + t347 * t369) * t294 + ((-t271 * t420 + t283 * t426) * MDP(38) + (t282 * t283 * MDP(37) - t339 * t423) * t292) * t438 + (t362 * t438 + 0.2e1 * t341 + 0.2e1 * t379 + t385 + 0.2e1 * t392 - 0.2e1 * t393) * t265; (-t231 * t270 - t232 * t264) * MDP(40) + t232 * t404 + t338 * t233 - t317 * t239 + t324 * t419 + (MDP(30) * t359 + t303) * t299 + (-MDP(31) * t359 - t305 * t292 + t320 * t298 - t442) * t293 - t355; (-t237 * t270 - t238 * t264) * MDP(40) + t238 * t404 + t354 * t295 + t338 * t246 - t317 * t253 + t324 * t301 + ((-MDP(28) * t295 + MDP(30) * t430) * t294 + (MDP(31) * t421 + t307) * pkin(3) + t306) * t299 + (t319 * t298 + (-MDP(31) * t430 + (pkin(3) * MDP(30) - MDP(27)) * t295) * t294 - t304 * t292 - t441) * t293; -t371 + t255 * t404 + (-t252 * t270 - t255 * t264) * MDP(40) + t317 * t265 + (t300 * MDP(28) + (-MDP(31) * t300 - t339 * t292) * pkin(3) + t347) * t299 + (-t288 * MDP(26) - MDP(20) + (-MDP(34) * t298 + MDP(35) * t292 + MDP(26)) * t282) * t294 + (t300 * MDP(27) + (t265 * MDP(33) + t380) * t298 + (-MDP(25) + MDP(36)) * t420 + (-t300 * MDP(30) - t265 * MDP(37) + (t271 + t367) * MDP(38)) * pkin(3) + (-t379 + pkin(3) * t362 - t385 + (-t264 * t429 - t242) * MDP(44) + (-t270 * t429 + t244) * MDP(45) + t345) * t292) * t293; 0.2e1 * t360 + t288 * MDP(36) + MDP(22) + (t264 * t434 + t404) * t270 + (MDP(32) * t287 + MDP(43) * t281 + MDP(25)) * t282 + 0.2e1 * (-t293 * t373 + t416 * t428) * t298 + (-0.2e1 * MDP(33) * t426 + (t416 * MDP(37) - t338 * t299) * t438 - t445 * t369) * t292; t305 * t298 + (t232 * t402 + (-t231 * t297 - t232 * t291) * MDP(40) + t336 * t233 - t444 * t239 + t328) * t292 + t314; -MDP(29) * t421 + t304 * t298 + (t238 * t402 + (-t237 * t297 - t238 * t291) * MDP(40) + t336 * t246 - t444 * t253 + t327) * t292 + t321; t300 * MDP(29) + (t311 + t379) * t298 + (t380 + t255 * t402 + (-t252 * t297 - t255 * t291) * MDP(40) + t444 * t265) * t292 + (-t412 + (-t292 * MDP(34) - t298 * MDP(35) + MDP(28)) * t293 + (t322 * t293 + t334 * t299) * pkin(3)) * t294; (t322 * pkin(3) + MDP(28)) * t299 + t445 * t298 + (-t373 + t270 * t402 + (-t264 * t297 - t270 * t291) * MDP(40)) * t292 + (t287 * MDP(33) + MDP(27) + (MDP(32) - MDP(43)) * t424 - t334 * pkin(3) - t444 * t281) * t293; MDP(43) * t287 + MDP(29) + (MDP(39) * t286 + MDP(32) - 0.2e1 * t361) * t281 + 0.2e1 * t444 * t424; -t232 * t403 + (t231 * t291 - t232 * t297) * MDP(40) + t333 * t233 - t310 * t239 + t326; -t238 * t403 + (t237 * t291 - t238 * t297) * MDP(40) + t333 * t246 - t310 * t253 + t325; t378 - t255 * t403 + (t252 * t291 - t255 * t297) * MDP(40) + (-MDP(36) + (-t292 * t333 - t407) * pkin(3)) * t423 + t310 * t265; -t299 * MDP(36) - t270 * t403 + (t264 * t291 - t270 * t297) * MDP(40) + (MDP(34) * t293 - t299 * t428) * t298 + (-pkin(3) * t299 * t333 - t293 * t310) * t292; (MDP(34) - t291 * t402 + (t280 - t286) * MDP(40)) * t292 + t310 * t298; MDP(39) * t280 + MDP(36) + 0.2e1 * t361; t313; t312; t311; -MDP(43) * t292 * t293 + t317; MDP(43) * t298 + t292 * t315; t316; MDP(43);];
%% Postprocessing: Reshape Output
% From vec2symmat_7_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16) t1(22); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17) t1(23); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18) t1(24); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19) t1(25); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20) t1(26); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21) t1(27); t1(22) t1(23) t1(24) t1(25) t1(26) t1(27) t1(28);];
Mq  = res;
