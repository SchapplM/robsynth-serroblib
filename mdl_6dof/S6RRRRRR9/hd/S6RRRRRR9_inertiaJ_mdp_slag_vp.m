% Calculate joint inertia matrix for
% S6RRRRRR9
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
%   see S6RRRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR9_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:22
% EndTime: 2019-03-10 05:39:31
% DurationCPUTime: 3.08s
% Computational Cost: add. (4977->394), mult. (12958->545), div. (0->0), fcn. (15021->14), ass. (0->191)
t297 = sin(qJ(5));
t302 = cos(qJ(5));
t296 = sin(qJ(6));
t301 = cos(qJ(6));
t271 = t296 * t297 - t301 * t302;
t272 = t296 * t302 + t297 * t301;
t393 = pkin(12) + pkin(13);
t278 = t393 * t297;
t279 = t393 * t302;
t332 = t272 * MDP(34) - t271 * MDP(35) + (-t278 * t301 - t279 * t296) * MDP(37) - (-t278 * t296 + t279 * t301) * MDP(38);
t419 = -(MDP(30) * t297 + MDP(31) * t302) * pkin(12) + t297 * MDP(27) + t302 * MDP(28) + t332;
t298 = sin(qJ(4));
t303 = cos(qJ(4));
t275 = -pkin(4) * t303 - pkin(12) * t298 - pkin(3);
t270 = t302 * t275;
t383 = pkin(13) * t298;
t386 = pkin(11) * t297;
t235 = -t302 * t383 + t270 + (-pkin(5) - t386) * t303;
t384 = pkin(11) * t303;
t339 = t302 * t384;
t241 = t339 + (t275 - t383) * t297;
t215 = t301 * t235 - t241 * t296;
t216 = t235 * t296 + t241 * t301;
t258 = t272 * t298;
t251 = t258 * MDP(35);
t259 = t271 * t298;
t253 = t259 * MDP(34);
t418 = MDP(37) * t215 - MDP(38) * t216 - t251 - t253;
t308 = -MDP(21) + t419;
t249 = -t297 * t384 + t270;
t250 = t275 * t297 + t339;
t417 = t249 * MDP(30) - t250 * MDP(31) + t418;
t292 = sin(pkin(7));
t415 = 0.2e1 * t292;
t414 = -0.2e1 * t298;
t317 = (MDP(37) * t301 - MDP(38) * t296) * pkin(5);
t293 = sin(pkin(6));
t295 = cos(pkin(6));
t299 = sin(qJ(3));
t300 = sin(qJ(2));
t304 = cos(qJ(3));
t294 = cos(pkin(7));
t305 = cos(qJ(2));
t368 = t294 * t305;
t372 = t292 * t299;
t237 = t295 * t372 + (t299 * t368 + t300 * t304) * t293;
t369 = t293 * t305;
t260 = -t292 * t369 + t294 * t295;
t223 = t237 * t298 - t260 * t303;
t341 = MDP(29) + MDP(36);
t413 = t341 * t223;
t262 = -t303 * t294 + t298 * t372;
t412 = t341 * t262;
t224 = t237 * t303 + t260 * t298;
t338 = t293 * t368;
t370 = t293 * t300;
t371 = t292 * t304;
t236 = -t295 * t371 + t299 * t370 - t304 * t338;
t210 = t224 * t297 - t236 * t302;
t211 = t224 * t302 + t236 * t297;
t196 = t301 * t210 + t211 * t296;
t197 = -t210 * t296 + t211 * t301;
t411 = t197 * MDP(34) - t196 * MDP(35);
t263 = t294 * t298 + t303 * t372;
t239 = t263 * t297 + t302 * t371;
t240 = t263 * t302 - t297 * t371;
t217 = t301 * t239 + t240 * t296;
t218 = -t239 * t296 + t240 * t301;
t410 = t218 * MDP(34) - t217 * MDP(35);
t280 = pkin(10) * t372;
t390 = pkin(2) * t304;
t254 = t280 + (-pkin(3) - t390) * t294;
t225 = pkin(4) * t262 - pkin(12) * t263 + t254;
t340 = pkin(10) * t371;
t391 = pkin(2) * t299;
t255 = t340 + (pkin(11) + t391) * t294;
t256 = (-pkin(3) * t304 - pkin(11) * t299 - pkin(2)) * t292;
t229 = t255 * t303 + t256 * t298;
t227 = -pkin(12) * t371 + t229;
t205 = t302 * t225 - t227 * t297;
t206 = t225 * t297 + t227 * t302;
t408 = t205 * MDP(30) - t206 * MDP(31);
t392 = pkin(1) * t305;
t282 = t295 * t392;
t238 = pkin(2) * t295 + t282 + (-pkin(10) * t294 - pkin(9)) * t370;
t248 = (-pkin(10) * t292 * t300 - pkin(2) * t305 - pkin(1)) * t293;
t219 = -t238 * t292 + t294 * t248;
t201 = pkin(3) * t236 - pkin(11) * t237 + t219;
t267 = t295 * t300 * pkin(1) + pkin(9) * t369;
t233 = (t292 * t295 + t338) * pkin(10) + t267;
t331 = t238 * t294 + t248 * t292;
t208 = t233 * t304 + t331 * t299;
t204 = pkin(11) * t260 + t208;
t191 = t201 * t298 + t204 * t303;
t189 = pkin(12) * t236 + t191;
t207 = -t299 * t233 + t331 * t304;
t203 = -pkin(3) * t260 - t207;
t193 = pkin(4) * t223 - pkin(12) * t224 + t203;
t183 = -t189 * t297 + t302 * t193;
t184 = t189 * t302 + t193 * t297;
t407 = t183 * MDP(30) - t184 * MDP(31);
t406 = t224 * MDP(18);
t405 = t260 * MDP(15);
t328 = t240 * MDP(27) - t239 * MDP(28);
t311 = t263 * MDP(19) - t328 - t410;
t388 = pkin(5) * t262;
t199 = -pkin(13) * t240 + t205 + t388;
t200 = -pkin(13) * t239 + t206;
t186 = t301 * t199 - t200 * t296;
t376 = t200 * t301;
t187 = t199 * t296 + t376;
t322 = t186 * MDP(37) - t187 * MDP(38);
t404 = t254 * MDP(23) - t311 + t322 + t408 + t412;
t403 = 2 * MDP(16);
t402 = 2 * MDP(17);
t401 = 0.2e1 * MDP(23);
t400 = 2 * MDP(24);
t399 = -2 * MDP(26);
t398 = 0.2e1 * MDP(30);
t397 = 0.2e1 * MDP(31);
t396 = -2 * MDP(33);
t395 = 0.2e1 * MDP(37);
t394 = 0.2e1 * MDP(38);
t389 = pkin(5) * t223;
t387 = pkin(5) * t301;
t385 = pkin(11) * t302;
t382 = MDP(24) * pkin(3);
t381 = pkin(3) * MDP(23);
t380 = pkin(11) * MDP(24);
t182 = -pkin(13) * t210 + t184;
t379 = t182 * t301;
t190 = t201 * t303 - t298 * t204;
t188 = -pkin(4) * t236 - t190;
t378 = t188 * t297;
t377 = t188 * t302;
t228 = -t298 * t255 + t256 * t303;
t226 = pkin(4) * t371 - t228;
t375 = t226 * t297;
t374 = t226 * t302;
t287 = t293 ^ 2;
t373 = t287 * t300;
t367 = t295 * MDP(8);
t366 = t297 * t302;
t364 = MDP(14) * t294;
t363 = MDP(32) * t272;
t358 = t197 * MDP(32);
t355 = t211 * MDP(25);
t354 = t218 * MDP(32);
t353 = t224 * MDP(20);
t352 = t236 * MDP(22);
t351 = t237 * MDP(11);
t350 = t237 * MDP(12);
t349 = t240 * MDP(25);
t348 = t260 * MDP(13);
t347 = t260 * MDP(14);
t346 = t263 * MDP(18);
t345 = t263 * MDP(20);
t344 = t294 * MDP(15);
t343 = t299 * MDP(13);
t342 = t304 * MDP(22);
t337 = t223 * MDP(36) + t411;
t336 = t262 * MDP(36) + t410;
t335 = MDP(26) * t366;
t334 = pkin(11) * MDP(23) - MDP(20);
t333 = -MDP(21) + t380;
t181 = -pkin(13) * t211 + t183 + t389;
t178 = t301 * t181 - t182 * t296;
t264 = t294 * t390 - t280;
t266 = t294 * t391 + t340;
t330 = t264 * MDP(16) - t266 * MDP(17);
t329 = t211 * MDP(27) - t210 * MDP(28);
t327 = MDP(27) * t302 - MDP(28) * t297;
t179 = t181 * t296 + t379;
t323 = t178 * MDP(37) - t179 * MDP(38);
t319 = -MDP(19) + t327;
t316 = (t300 * MDP(6) + t305 * MDP(7)) * t293;
t315 = t228 * MDP(23) - t229 * MDP(24) + t345;
t314 = t237 * MDP(13) + t207 * MDP(16) - t208 * MDP(17) + t405;
t313 = t190 * MDP(23) - t191 * MDP(24) + t352 + t353;
t312 = t224 * MDP(19) - t329 - t411;
t309 = -t381 + t417;
t307 = -t203 * MDP(23) + t312 - t323 - t407;
t290 = t302 ^ 2;
t288 = t297 ^ 2;
t286 = t292 ^ 2;
t285 = -pkin(5) * t302 - pkin(4);
t274 = (pkin(5) * t297 + pkin(11)) * t298;
t265 = -pkin(9) * t370 + t282;
t212 = pkin(5) * t239 + t226;
t185 = pkin(5) * t210 + t188;
t1 = [MDP(1) + (MDP(4) * t300 + 0.2e1 * MDP(5) * t305) * t373 + (t210 * t399 + t355) * t211 + (t196 * t396 + t358) * t197 + (0.2e1 * t316 + t367) * t295 + 0.2e1 * (-pkin(1) * t373 - t267 * t295) * MDP(10) + 0.2e1 * (t265 * t295 + t287 * t392) * MDP(9) + (t219 * t402 + 0.2e1 * t348 + t351) * t237 + (t207 * t403 - t208 * t402 + t405) * t260 + (t203 * t400 + t406) * t224 + (t190 * t401 - t191 * t400 + t219 * t403 - 0.2e1 * t347 - 0.2e1 * t350 + t352 + 0.2e1 * t353) * t236 + (t210 * t398 + t211 * t397) * t188 + (t196 * t395 + t197 * t394) * t185 + (-0.2e1 * t236 * MDP(21) + t178 * t395 - t179 * t394 + t183 * t398 - t184 * t397 + t203 * t401 - 0.2e1 * t312 + t413) * t223; t211 * t349 + (t188 * t240 + t211 * t226) * MDP(31) + (t185 * t218 + t197 * t212) * MDP(38) + (t188 * t239 + t210 * t226) * MDP(30) + (t185 * t217 + t196 * t212) * MDP(37) + t224 * t346 + t197 * t354 + t367 - t267 * MDP(10) + t265 * MDP(9) + (-t210 * t240 - t211 * t239) * MDP(26) + (t203 * t263 + t224 * t254) * MDP(24) + (-t196 * t218 - t197 * t217) * MDP(33) + t330 * t260 + t314 * t294 + t316 - t307 * t262 + t404 * t223 + (-t262 * MDP(21) + t315 - t364) * t236 + ((-t236 * MDP(16) - t237 * MDP(17)) * pkin(2) + (-t236 * MDP(12) + t219 * MDP(17) + t348 + t351) * t299 + (-t219 * MDP(16) + t223 * MDP(21) - t313 + t347 + t350) * t304) * t292; t286 * t299 ^ 2 * MDP(11) + t263 ^ 2 * MDP(18) + MDP(8) + (t343 * t415 + t344) * t294 + (t239 * t399 + t349) * t240 + (t217 * t396 + t354) * t218 + ((-t345 + t364) * t415 + (0.2e1 * t299 * MDP(12) + t342) * t286) * t304 + (-t266 * t294 - t286 * t391) * t402 + (t264 * t294 + t286 * t390) * t403 + (t229 * t371 + t254 * t263) * t400 - t228 * t371 * t401 + (t239 * t398 + t240 * t397) * t226 + (t217 * t395 + t218 * t394) * t212 + (0.2e1 * MDP(21) * t371 + t186 * t395 - t187 * t394 + t205 * t398 - t206 * t397 + t254 * t401 - 0.2e1 * t311 + t412) * t262; -t224 * t382 + (-t185 * t259 + t197 * t274) * MDP(38) + (t185 * t258 + t196 * t274) * MDP(37) + (t196 * t259 - t197 * t258) * MDP(33) - t236 * MDP(14) - t259 * t358 + t309 * t223 + (-t333 * t236 + t307 - t413) * t303 + (t406 + (pkin(11) * t211 + t377) * MDP(31) + (pkin(11) * t210 + t378) * MDP(30) + t203 * MDP(24) + (-t210 * t302 - t211 * t297) * MDP(26) + t302 * t355 - t334 * t236 + t319 * t223) * t298 + t314; (-t212 * t259 + t218 * t274) * MDP(38) + (t212 * t258 + t217 * t274) * MDP(37) - t263 * t382 + t344 + (t217 * t259 - t218 * t258) * MDP(33) - t259 * t354 + (t304 * MDP(14) + t343) * t292 + t309 * t262 + (t333 * t371 - t404) * t303 + (t346 + t254 * MDP(24) + (-t239 * t302 - t240 * t297) * MDP(26) + (pkin(11) * t240 + t374) * MDP(31) + (pkin(11) * t239 + t375) * MDP(30) + t302 * t349 + t334 * t371 + t319 * t262) * t298 + t330; t258 * t274 * t395 + t382 * t414 + MDP(15) - (-MDP(32) * t259 + t258 * t396 + t274 * t394) * t259 + (MDP(25) * t290 + t385 * t397 + t386 * t398 + MDP(18) - 0.2e1 * t335) * t298 ^ 2 + (-t215 * t395 + t216 * t394 - t249 * t398 + t250 * t397 + t341 * t303 + t319 * t414 + 0.2e1 * t251 + 0.2e1 * t253 + 0.2e1 * t381) * t303; t297 * t355 + (-t210 * t297 + t211 * t302) * MDP(26) + (-pkin(4) * t210 - t377) * MDP(30) + (-pkin(4) * t211 + t378) * MDP(31) + t272 * t358 + (-t196 * t272 - t197 * t271) * MDP(33) + (t185 * t271 + t196 * t285) * MDP(37) + (t185 * t272 + t197 * t285) * MDP(38) + t308 * t223 + t313; -t292 * t342 + t297 * t349 + (-t239 * t297 + t240 * t302) * MDP(26) + (-pkin(4) * t239 - t374) * MDP(30) + (-pkin(4) * t240 + t375) * MDP(31) + t272 * t354 + (-t217 * t272 - t218 * t271) * MDP(33) + (t212 * t271 + t217 * t285) * MDP(37) + (t212 * t272 + t218 * t285) * MDP(38) + t308 * t262 + t315; -t259 * t363 + (-t258 * t272 + t259 * t271) * MDP(33) + (t258 * t285 + t271 * t274) * MDP(37) + (-t259 * t285 + t272 * t274) * MDP(38) + (-t380 - t308) * t303 + (MDP(25) * t366 + (-t288 + t290) * MDP(26) + (-pkin(4) * t297 - t385) * MDP(30) + (-pkin(4) * t302 + t386) * MDP(31) - t334) * t298; 0.2e1 * t335 + t271 * t285 * t395 + MDP(25) * t288 + MDP(22) + 0.2e1 * (MDP(30) * t302 - MDP(31) * t297) * pkin(4) + (t271 * t396 + t285 * t394 + t363) * t272; t223 * MDP(29) + (t223 * t387 + t178) * MDP(37) + (-t379 + (-t181 - t389) * t296) * MDP(38) + t329 + t337 + t407; t262 * MDP(29) + (t262 * t387 + t186) * MDP(37) + (-t376 + (-t199 - t388) * t296) * MDP(38) + t328 + t336 + t408; (-t341 - t317) * t303 + t327 * t298 + t417; t419; 0.2e1 * t317 + t341; t323 + t337; t322 + t336; -t303 * MDP(36) + t418; t332; MDP(36) + t317; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
