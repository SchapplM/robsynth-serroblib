% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:17
% EndTime: 2019-12-31 21:51:21
% DurationCPUTime: 1.70s
% Computational Cost: add. (2381->239), mult. (4005->308), div. (0->0), fcn. (2460->6), ass. (0->128)
t377 = pkin(8) + pkin(7);
t302 = cos(qJ(4));
t303 = cos(qJ(3));
t339 = qJD(4) * t302;
t341 = qJD(3) * t303;
t376 = -t302 * t341 - t303 * t339;
t300 = sin(qJ(3));
t375 = t300 * MDP(12) + t303 * MDP(13);
t369 = -MDP(19) - MDP(21);
t374 = MDP(7) * t300;
t299 = sin(qJ(4));
t265 = t299 * t303 + t300 * t302;
t296 = qJD(1) + qJD(2);
t255 = t265 * t296;
t337 = t303 * MDP(12);
t373 = (-t300 * MDP(13) + MDP(5) + t337) * t296;
t372 = (t300 ^ 2 - t303 ^ 2) * MDP(8);
t354 = t302 * t303;
t356 = t299 * t300;
t264 = -t354 + t356;
t327 = qJD(3) * t377;
t266 = t300 * t327;
t267 = t303 * t327;
t276 = t377 * t300;
t293 = t303 * pkin(8);
t277 = pkin(7) * t303 + t293;
t313 = -t276 * t302 - t277 * t299;
t304 = cos(qJ(2));
t363 = pkin(1) * qJD(1);
t333 = t304 * t363;
t350 = t313 * qJD(4) + t264 * t333 - t266 * t302 - t267 * t299;
t245 = -t276 * t299 + t277 * t302;
t349 = t245 * qJD(4) - t265 * t333 - t266 * t299 + t267 * t302;
t301 = sin(qJ(2));
t322 = t377 * t296 + t301 * t363;
t248 = t322 * t300;
t249 = t322 * t303;
t358 = t249 * t299;
t371 = -pkin(3) * t339 - t248 * t302 - t358;
t370 = t375 * qJD(3);
t368 = -MDP(20) + MDP(23);
t295 = qJD(3) + qJD(4);
t367 = t255 ^ 2;
t365 = pkin(1) * t304;
t286 = pkin(1) * t301 + pkin(7);
t364 = -pkin(8) - t286;
t362 = pkin(1) * qJD(2);
t321 = qJD(3) * t364;
t332 = t304 * t362;
t246 = t300 * t321 + t303 * t332;
t247 = -t300 * t332 + t303 * t321;
t261 = t364 * t300;
t262 = t286 * t303 + t293;
t314 = t261 * t302 - t262 * t299;
t196 = t314 * qJD(4) + t246 * t302 + t247 * t299;
t361 = t196 * t295;
t233 = t261 * t299 + t262 * t302;
t197 = t233 * qJD(4) + t246 * t299 - t247 * t302;
t360 = t197 * t295;
t329 = t296 * t356;
t253 = -t296 * t354 + t329;
t289 = -pkin(3) * t303 - pkin(2);
t257 = t289 * t296 - t333;
t210 = pkin(4) * t253 - qJ(5) * t255 + t257;
t359 = t210 * t255;
t357 = t249 * t302;
t305 = qJD(3) ^ 2;
t355 = t300 * t305;
t353 = t303 * t305;
t316 = t295 * t356;
t344 = t376 * t296;
t220 = t296 * t316 + t344;
t238 = t295 * t265;
t221 = t238 * t296;
t330 = qJD(1) * t362;
t282 = t301 * t330;
t342 = qJD(3) * t300;
t326 = t296 * t342;
t258 = pkin(3) * t326 + t282;
t194 = pkin(4) * t221 + qJ(5) * t220 - qJD(5) * t255 + t258;
t352 = t194 * t264 + t210 * t238;
t237 = t316 + t376;
t351 = -t194 * t265 + t210 * t237;
t348 = qJD(5) - t371;
t347 = t257 * t238 + t258 * t264;
t346 = -t257 * t237 + t258 * t265;
t270 = -pkin(2) * t296 - t333;
t345 = t270 * t341 + t300 * t282;
t340 = qJD(4) * t299;
t243 = qJD(3) * pkin(3) - t248;
t213 = t243 * t302 - t358;
t335 = qJD(5) - t213;
t334 = pkin(3) * t296 * t300;
t290 = pkin(3) * t342;
t206 = -t344 + (t253 - t329) * t295;
t328 = t206 * MDP(16) + (-t253 ^ 2 + t367) * MDP(15);
t325 = t296 * t341;
t320 = t304 * t330;
t292 = t295 * qJD(5);
t315 = qJD(3) * t322;
t230 = -t300 * t315 + t303 * t320;
t231 = -t300 * t320 - t303 * t315;
t318 = t302 * t230 + t299 * t231 + t243 * t339 - t249 * t340;
t192 = t292 + t318;
t193 = t299 * t230 - t302 * t231 + t243 * t340 + t249 * t339;
t207 = -pkin(4) * t295 + t335;
t214 = t243 * t299 + t357;
t209 = qJ(5) * t295 + t214;
t319 = -t192 * t264 + t193 * t265 - t207 * t237 - t209 * t238;
t317 = pkin(3) * t340 + t248 * t299 - t357;
t222 = pkin(4) * t255 + qJ(5) * t253;
t310 = t213 * t295 - t318;
t309 = t214 * t295 - t193;
t234 = pkin(4) * t264 - qJ(5) * t265 + t289;
t307 = t255 * MDP(14) + t257 * MDP(20) - t210 * MDP(23);
t203 = pkin(4) * t238 + qJ(5) * t237 - qJD(5) * t265 + t290;
t306 = -MDP(10) * t355 - t282 * MDP(5) + (t220 * t264 - t221 * t265 + t237 * t253 - t238 * t255) * MDP(15) + (-t220 * t265 - t237 * t255) * MDP(14) - 0.2e1 * t296 * qJD(3) * t372 + 0.2e1 * t325 * t374 + MDP(9) * t353 + (-t237 * MDP(16) - t238 * MDP(17)) * t295;
t291 = t301 * t362;
t288 = -pkin(2) - t365;
t287 = -pkin(3) * t302 - pkin(4);
t283 = pkin(3) * t299 + qJ(5);
t273 = t289 - t365;
t268 = t291 + t290;
t259 = t270 * t342;
t227 = t234 - t365;
t217 = t222 + t334;
t202 = t203 + t291;
t1 = [(-t286 * t353 + t288 * t326 + t259) * MDP(12) + (t286 * t355 + t288 * t325 + t345) * MDP(13) + (t221 * t273 + t253 * t268 + t347 - t360) * MDP(19) + (t202 * t253 + t221 * t227 + t352 - t360) * MDP(21) + t306 + (-t220 * t273 + t268 * t255 + t346 - t361) * MDP(20) + (-t196 * t253 + t197 * t255 + t220 * t314 - t221 * t233 + t319) * MDP(22) + (-t202 * t255 + t220 * t227 + t351 + t361) * MDP(23) + (t192 * t233 - t193 * t314 + t194 * t227 + t196 * t209 + t197 * t207 + t202 * t210) * MDP(24) + (((-qJD(1) - t296) * MDP(6) - t370) * t304 + (-qJD(1) * t337 - t373) * t301) * t362; (-pkin(2) * t326 - pkin(7) * t353 + t259) * MDP(12) + (-pkin(2) * t325 + pkin(7) * t355 + t345) * MDP(13) + (t192 * t245 - t193 * t313 + t194 * t234 + t203 * t210 + t349 * t207 + t350 * t209) * MDP(24) + (((-qJD(2) + t296) * MDP(6) + t370) * t304 + (-t210 * MDP(24) - qJD(2) * t337 + t369 * t253 + t368 * t255 + t373) * t301) * t363 + (t369 * t349 + t368 * t350) * t295 + t306 + (-t220 * t289 + t255 * t290 + t346) * MDP(20) + (t220 * t313 - t221 * t245 - t350 * t253 + t349 * t255 + t319) * MDP(22) + (-t203 * t255 + t220 * t234 + t351) * MDP(23) + (t221 * t289 + t253 * t290 + t347) * MDP(19) + (t203 * t253 + t221 * t234 + t352) * MDP(21); -t318 * MDP(20) + (-t220 * t287 - t221 * t283) * MDP(22) + t192 * MDP(23) + (t192 * t283 + t317 * t207 + t348 * t209 - t210 * t217) * MDP(24) + (-t257 * MDP(19) - MDP(20) * t334 - t210 * MDP(21) + (t209 + t317) * MDP(22) + t217 * MDP(23)) * t255 + (-t303 * t374 + t372) * t296 ^ 2 + (-MDP(19) * t334 - t217 * MDP(21) + (t207 - t348) * MDP(22) + t307) * t253 + (t371 * MDP(20) + t348 * MDP(23) + t369 * t317) * t295 + t328 + t375 * (-t270 * t296 - t320) + (t287 * MDP(24) + t369) * t193; (-t255 * t257 + t309) * MDP(19) + t310 * MDP(20) + (t309 - t359) * MDP(21) + (pkin(4) * t220 - qJ(5) * t221 - (-t209 + t214) * t255) * MDP(22) + (t222 * t255 + 0.2e1 * t292 - t310) * MDP(23) + (-pkin(4) * t193 + qJ(5) * t192 - t207 * t214 + t335 * t209 - t210 * t222) * MDP(24) + (-t222 * MDP(21) + (t207 - t335) * MDP(22) + t307) * t253 + t328; t255 * t253 * MDP(21) + t206 * MDP(22) + (-t295 ^ 2 - t367) * MDP(23) + (-t209 * t295 + t193 + t359) * MDP(24);];
tauc = t1;
