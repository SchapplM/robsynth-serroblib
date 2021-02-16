% Calculate Coriolis joint torque vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:42
% EndTime: 2021-01-15 19:47:53
% DurationCPUTime: 3.86s
% Computational Cost: add. (2513->333), mult. (6625->473), div. (0->0), fcn. (4816->8), ass. (0->148)
t339 = cos(qJ(2));
t389 = cos(pkin(8));
t357 = t389 * t339;
t323 = qJD(1) * t357;
t334 = sin(pkin(8));
t337 = sin(qJ(2));
t367 = qJD(1) * t337;
t297 = t334 * t367 - t323;
t294 = qJD(5) + t297;
t333 = sin(pkin(9));
t335 = cos(pkin(9));
t313 = t334 * t339 + t389 * t337;
t369 = qJD(1) * t313;
t282 = t335 * qJD(2) - t333 * t369;
t338 = cos(qJ(5));
t398 = t282 * t338;
t336 = sin(qJ(5));
t314 = t333 * t338 + t335 * t336;
t371 = t294 * t314;
t281 = qJD(2) * t333 + t335 * t369;
t240 = t281 * t338 + t282 * t336;
t363 = qJD(1) * qJD(2);
t397 = -0.2e1 * t363;
t396 = (t337 ^ 2 - t339 ^ 2) * MDP(5);
t299 = t313 * qJD(2);
t291 = qJD(1) * t299;
t376 = t338 * t335;
t378 = t333 * t336;
t312 = -t376 + t378;
t372 = t294 * t312;
t395 = -t291 * t314 + t372 * t294;
t296 = t297 ^ 2;
t394 = pkin(2) * t337;
t393 = pkin(7) * t335;
t325 = pkin(2) * t334 + qJ(4);
t392 = pkin(7) + t325;
t391 = -qJ(3) - pkin(6);
t390 = qJD(2) * pkin(2);
t238 = t281 * t336 - t398;
t388 = t238 * t369;
t387 = t240 * t369;
t359 = qJD(2) * t391;
t295 = qJD(3) * t339 + t337 * t359;
t287 = t295 * qJD(1);
t344 = -qJD(3) * t337 + t339 * t359;
t288 = t344 * qJD(1);
t245 = t287 * t334 - t389 * t288;
t320 = t391 * t337;
t321 = t391 * t339;
t277 = -t389 * t320 - t321 * t334;
t386 = t245 * t277;
t360 = t337 * t363;
t322 = t334 * t360;
t292 = qJD(2) * t323 - t322;
t384 = t292 * t333;
t383 = t292 * t335;
t382 = t297 * t333;
t345 = -t334 * t337 + t357;
t302 = t345 * qJD(2);
t381 = t302 * t333;
t380 = t313 * t333;
t379 = t313 * t335;
t317 = qJD(1) * t321;
t305 = t334 * t317;
t340 = qJD(2) ^ 2;
t377 = t337 * t340;
t375 = t339 * t340;
t341 = qJD(1) ^ 2;
t374 = t339 * t341;
t324 = pkin(2) * t360;
t233 = pkin(3) * t291 - qJ(4) * t292 - qJD(4) * t369 + t324;
t246 = t389 * t287 + t334 * t288;
t243 = qJD(2) * qJD(4) + t246;
t211 = t333 * t233 + t335 * t243;
t361 = t337 * t390;
t242 = pkin(3) * t299 - qJ(4) * t302 - qJD(4) * t313 + t361;
t259 = t389 * t295 + t334 * t344;
t218 = t333 * t242 + t335 * t259;
t329 = -pkin(2) * t339 - pkin(1);
t368 = qJD(1) * t329;
t319 = qJD(3) + t368;
t249 = pkin(3) * t297 - qJ(4) * t369 + t319;
t316 = qJD(1) * t320;
t310 = t316 + t390;
t358 = t389 * t317;
t269 = t334 * t310 - t358;
t265 = qJD(2) * qJ(4) + t269;
t222 = t333 * t249 + t335 * t265;
t362 = pkin(2) * t367;
t257 = pkin(3) * t369 + qJ(4) * t297 + t362;
t273 = t389 * t316 + t305;
t227 = t333 * t257 + t335 * t273;
t267 = -pkin(3) * t345 - qJ(4) * t313 + t329;
t278 = t334 * t320 - t389 * t321;
t230 = t333 * t267 + t335 * t278;
t364 = qJD(5) * t338;
t373 = t282 * t364 + t292 * t376;
t214 = pkin(7) * t282 + t222;
t366 = qJD(5) * t214;
t365 = qJD(5) * t313;
t356 = pkin(1) * t397;
t210 = t335 * t233 - t243 * t333;
t217 = t335 * t242 - t259 * t333;
t221 = t335 * t249 - t265 * t333;
t226 = t335 * t257 - t273 * t333;
t229 = t335 * t267 - t278 * t333;
t258 = t295 * t334 - t389 * t344;
t272 = t316 * t334 - t358;
t355 = 0.2e1 * t369;
t207 = -pkin(7) * t384 + t211;
t209 = pkin(4) * t297 - pkin(7) * t281 + t221;
t354 = -qJD(5) * t209 - t207;
t328 = -t389 * pkin(2) - pkin(3);
t353 = -t312 * t291 - t294 * t371;
t268 = t389 * t310 + t305;
t203 = t209 * t338 - t214 * t336;
t204 = t209 * t336 + t214 * t338;
t219 = -pkin(4) * t345 - pkin(7) * t379 + t229;
t223 = -pkin(7) * t380 + t230;
t352 = t219 * t338 - t223 * t336;
t351 = t219 * t336 + t223 * t338;
t350 = -t221 * t333 + t222 * t335;
t349 = t245 * t313 + t277 * t292;
t348 = -qJD(5) * t281 - t384;
t309 = t392 * t335;
t347 = pkin(4) * t369 + qJD(4) * t333 + qJD(5) * t309 + t297 * t393 + t226;
t308 = t392 * t333;
t346 = pkin(7) * t382 - qJD(4) * t335 + qJD(5) * t308 + t227;
t260 = -qJD(2) * pkin(3) + qJD(4) - t268;
t343 = t260 * t302 + t349;
t342 = -t291 * t325 + t292 * t328 + (-qJD(4) + t260) * t297;
t216 = t240 * qJD(5) + t314 * t292;
t318 = -t335 * pkin(4) + t328;
t262 = t312 * t313;
t261 = t314 * t313;
t255 = pkin(4) * t380 + t277;
t244 = -pkin(4) * t382 + t272;
t235 = pkin(4) * t381 + t258;
t234 = -pkin(4) * t282 + t260;
t228 = pkin(4) * t384 + t245;
t225 = t314 * t302 + t364 * t379 - t365 * t378;
t224 = -t312 * t302 - t314 * t365;
t215 = t348 * t336 + t373;
t212 = -pkin(7) * t381 + t218;
t208 = pkin(4) * t299 - t302 * t393 + t217;
t206 = pkin(4) * t291 - pkin(7) * t383 + t210;
t205 = t338 * t206;
t1 = [0.2e1 * t339 * MDP(4) * t360 + t396 * t397 + MDP(6) * t375 - MDP(7) * t377 + (-pkin(6) * t375 + t337 * t356) * MDP(9) + (pkin(6) * t377 + t339 * t356) * MDP(10) + (t291 * t329 + t299 * t319 + (-t258 + (-qJD(1) * t345 + t297) * t394) * qJD(2)) * MDP(11) + (t292 * t329 + t302 * t319 + (t355 * t394 - t259) * qJD(2)) * MDP(12) + (t246 * t345 + t258 * t369 - t259 * t297 - t268 * t302 - t269 * t299 - t278 * t291 + t349) * MDP(13) + (t386 + t246 * t278 - t258 * t268 + t259 * t269 + (t319 + t368) * t361) * MDP(14) + (-t210 * t345 + t217 * t297 + t221 * t299 + t229 * t291 - t258 * t282 + t343 * t333) * MDP(15) + (t211 * t345 - t218 * t297 - t222 * t299 - t230 * t291 + t258 * t281 + t343 * t335) * MDP(16) + (-t217 * t281 + t218 * t282 + (-t210 * t313 - t221 * t302 - t229 * t292) * t335 + (-t211 * t313 - t222 * t302 - t230 * t292) * t333) * MDP(17) + (t210 * t229 + t211 * t230 + t217 * t221 + t218 * t222 + t258 * t260 + t386) * MDP(18) + (-t215 * t262 + t224 * t240) * MDP(19) + (-t215 * t261 + t216 * t262 - t224 * t238 - t225 * t240) * MDP(20) + (-t215 * t345 + t224 * t294 + t240 * t299 - t262 * t291) * MDP(21) + (t216 * t345 - t225 * t294 - t238 * t299 - t261 * t291) * MDP(22) + (-t291 * t345 + t294 * t299) * MDP(23) + ((t208 * t338 - t212 * t336) * t294 + t352 * t291 - (-t207 * t336 + t205) * t345 + t203 * t299 + t235 * t238 + t255 * t216 + t228 * t261 + t234 * t225 + (t204 * t345 - t351 * t294) * qJD(5)) * MDP(24) + (-(t208 * t336 + t212 * t338) * t294 - t351 * t291 + (t206 * t336 + t207 * t338) * t345 - t204 * t299 + t235 * t240 + t255 * t215 - t228 * t262 + t234 * t224 + (t203 * t345 - t352 * t294) * qJD(5)) * MDP(25); -t337 * MDP(4) * t374 + t341 * t396 + (qJD(2) * t272 - t297 * t362 - t319 * t369 - t245) * MDP(11) + (qJD(2) * t273 + t297 * t319 - t362 * t369 - t246) * MDP(12) + ((t269 - t272) * t369 + (-t268 + t273) * t297 + (-t291 * t334 - t389 * t292) * pkin(2)) * MDP(13) + (t268 * t272 - t269 * t273 + (-t389 * t245 + t246 * t334 - t319 * t367) * pkin(2)) * MDP(14) + (-t221 * t369 - t226 * t297 - t245 * t335 + t272 * t282 + t342 * t333) * MDP(15) + (t222 * t369 + t227 * t297 + t245 * t333 - t272 * t281 + t342 * t335) * MDP(16) + (t226 * t281 - t227 * t282 + (qJD(4) * t282 - t221 * t297 + t211) * t335 + (qJD(4) * t281 - t222 * t297 - t210) * t333) * MDP(17) + (-t221 * t226 - t222 * t227 + t245 * t328 - t260 * t272 + (-t210 * t333 + t211 * t335) * t325 + t350 * qJD(4)) * MDP(18) + (t215 * t314 - t372 * t240) * MDP(19) + (-t215 * t312 - t216 * t314 + t372 * t238 - t371 * t240) * MDP(20) + (-t387 - t395) * MDP(21) + (t353 + t388) * MDP(22) - t294 * t369 * MDP(23) + ((-t308 * t338 - t309 * t336) * t291 + t318 * t216 + t228 * t312 - t203 * t369 - t244 * t238 + (t346 * t336 - t347 * t338) * t294 + t371 * t234) * MDP(24) + (-(-t308 * t336 + t309 * t338) * t291 + t318 * t215 + t228 * t314 + t204 * t369 - t244 * t240 + (t347 * t336 + t346 * t338) * t294 - t372 * t234) * MDP(25) + (t341 * t337 * MDP(9) + MDP(10) * t374) * pkin(1); t355 * qJD(2) * MDP(11) + (-t322 + (t323 - t297) * qJD(2)) * MDP(12) + (-t369 ^ 2 - t296) * MDP(13) + (t268 * t369 + t269 * t297 + t324) * MDP(14) + (t282 * t369 + t291 * t335 - t296 * t333) * MDP(15) + (-t281 * t369 - t291 * t333 - t296 * t335) * MDP(16) + ((t281 * t333 + t282 * t335) * t297 + (-t333 ^ 2 - t335 ^ 2) * t292) * MDP(17) + (t210 * t335 + t211 * t333 - t260 * t369 + t350 * t297) * MDP(18) + (t353 - t388) * MDP(24) + (-t387 + t395) * MDP(25); (t281 * t297 + t384) * MDP(15) + (t282 * t297 + t383) * MDP(16) + (-t281 ^ 2 - t282 ^ 2) * MDP(17) + (t221 * t281 - t222 * t282 + t245) * MDP(18) + (t240 * t294 + t216) * MDP(24) + (t294 * t398 + (-t281 * t294 + t348) * t336 + t373) * MDP(25); -t238 ^ 2 * MDP(20) + (t238 * t294 + t373) * MDP(21) + t291 * MDP(23) + (t204 * t294 + t205) * MDP(24) + (t203 * t294 + t234 * t238) * MDP(25) + (MDP(19) * t238 + MDP(20) * t240 + MDP(22) * t294 - MDP(24) * t234) * t240 + (t348 * MDP(22) - MDP(24) * t366 + t354 * MDP(25)) * t338 + (t348 * MDP(21) + (-qJD(5) * t282 - t383) * MDP(22) + t354 * MDP(24) + (-t206 + t366) * MDP(25)) * t336;];
tauc = t1;
