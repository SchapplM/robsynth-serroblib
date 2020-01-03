% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:04
% EndTime: 2019-12-31 19:48:09
% DurationCPUTime: 2.52s
% Computational Cost: add. (1432->306), mult. (3417->431), div. (0->0), fcn. (2089->6), ass. (0->155)
t385 = pkin(3) + pkin(6);
t314 = sin(qJ(2));
t366 = qJD(1) * t314;
t297 = qJD(5) + t366;
t311 = cos(pkin(8));
t316 = cos(qJ(2));
t365 = qJD(1) * t316;
t351 = t311 * t365;
t310 = sin(pkin(8));
t364 = qJD(2) * t310;
t268 = t351 + t364;
t349 = t310 * t365;
t363 = qJD(2) * t311;
t270 = -t349 + t363;
t313 = sin(qJ(5));
t315 = cos(qJ(5));
t331 = t268 * t313 - t270 * t315;
t392 = t297 * t331;
t355 = qJD(1) * qJD(2);
t391 = -0.2e1 * t355;
t390 = t314 * MDP(4);
t308 = t314 ^ 2;
t389 = (-t316 ^ 2 + t308) * MDP(5);
t272 = t310 * t315 + t311 * t313;
t260 = t272 * qJD(5);
t299 = pkin(6) * t366;
t388 = qJD(3) + t299;
t307 = qJD(2) * qJ(3);
t387 = qJD(4) + t307;
t300 = pkin(6) * t365;
t301 = pkin(3) * t365;
t279 = t300 + t301;
t258 = t279 + t387;
t312 = -pkin(2) - qJ(4);
t361 = qJD(2) * t316;
t386 = t314 * (-t258 + t387) - t312 * t361;
t384 = -pkin(7) + t312;
t383 = qJD(2) * pkin(2);
t345 = -qJ(3) * t314 - pkin(1);
t267 = t312 * t316 + t345;
t242 = t267 * qJD(1);
t356 = pkin(3) * t366 + t388;
t247 = qJD(2) * t312 + t356;
t211 = -t242 * t310 + t311 * t247;
t382 = t211 * t316;
t212 = t311 * t242 + t310 * t247;
t381 = t212 * t316;
t378 = t270 * t313;
t224 = t315 * t268 + t378;
t380 = t224 * t297;
t362 = qJD(2) * t314;
t278 = t385 * t362;
t306 = qJD(2) * qJD(3);
t252 = -qJD(1) * t278 + t306;
t379 = t252 * t310;
t377 = t310 * t313;
t376 = t310 * t314;
t375 = t311 * t314;
t374 = t311 * t316;
t317 = qJD(2) ^ 2;
t373 = t314 * t317;
t372 = t316 * t317;
t318 = qJD(1) ^ 2;
t371 = t316 * t318;
t348 = t314 * t355;
t296 = pkin(2) * t348;
t337 = -qJ(3) * t316 + qJ(4) * t314;
t360 = qJD(3) * t314;
t320 = qJD(2) * t337 - qJD(4) * t316 - t360;
t223 = qJD(1) * t320 + t296;
t347 = t316 * t355;
t295 = pkin(6) * t347;
t253 = t295 + (-qJD(4) + t301) * qJD(2);
t206 = t311 * t223 + t310 * t253;
t302 = pkin(2) * t362;
t232 = t302 + t320;
t280 = t385 * t361;
t214 = t311 * t232 + t310 * t280;
t330 = -t311 * t315 + t377;
t370 = -t260 * t297 - t330 * t347;
t303 = pkin(2) * t366;
t254 = qJD(1) * t337 + t303;
t222 = t311 * t254 + t310 * t279;
t291 = t385 * t314;
t228 = t311 * t267 + t310 * t291;
t323 = t272 * t314;
t241 = qJD(1) * t323;
t369 = -t260 - t241;
t352 = t311 * t366;
t358 = qJD(5) * t315;
t359 = qJD(5) * t313;
t368 = -t310 * t359 + t311 * t358 + t315 * t352 - t366 * t377;
t292 = t385 * t316;
t288 = -pkin(2) * t316 + t345;
t266 = qJD(1) * t288;
t353 = -pkin(4) * t311 - pkin(3);
t357 = -t353 * t366 + t388;
t339 = t315 * t348;
t340 = t313 * t348;
t354 = -t268 * t358 + t310 * t339 + t311 * t340;
t346 = MDP(23) * t365;
t344 = pkin(1) * t391;
t343 = qJD(3) - t383;
t205 = -t223 * t310 + t311 * t253;
t328 = pkin(4) * t316 - pkin(7) * t376;
t321 = t328 * qJD(2);
t200 = qJD(1) * t321 + t205;
t341 = pkin(7) * t311 * t362;
t201 = qJD(1) * t341 + t206;
t342 = t315 * t200 - t201 * t313;
t213 = -t232 * t310 + t311 * t280;
t221 = -t254 * t310 + t311 * t279;
t336 = t200 * t313 + t201 * t315;
t202 = pkin(4) * t366 - pkin(7) * t270 + t211;
t203 = -pkin(7) * t268 + t212;
t197 = t202 * t315 - t203 * t313;
t198 = t202 * t313 + t203 * t315;
t335 = t205 * t311 + t206 * t310;
t334 = -t211 * t310 + t212 * t311;
t275 = t311 * t291;
t216 = pkin(4) * t314 + t275 + (pkin(7) * t316 - t267) * t310;
t219 = -pkin(7) * t374 + t228;
t333 = t216 * t315 - t219 * t313;
t332 = t216 * t313 + t219 * t315;
t329 = -0.2e1 * qJD(2) * t266;
t326 = -qJ(3) * t361 - t360;
t243 = qJD(1) * t326 + t296;
t256 = t302 + t326;
t327 = pkin(6) * t317 + qJD(1) * t256 + t243;
t245 = (-pkin(6) + t353) * t362;
t281 = t384 * t310;
t325 = qJD(1) * t328 + qJD(4) * t311 + qJD(5) * t281 + t221;
t282 = t384 * t311;
t324 = pkin(7) * t352 + qJD(4) * t310 - qJD(5) * t282 + t222;
t249 = t330 * t316;
t207 = -t270 * t359 + t354;
t322 = t310 * t340 - t311 * t339;
t286 = pkin(6) * t348 - t306;
t287 = t299 + t343;
t290 = -t300 - t307;
t319 = -t286 * t316 + (t287 * t316 + (t290 + t300) * t314) * qJD(2);
t208 = -qJD(5) * t331 + t322;
t298 = pkin(4) * t310 + qJ(3);
t276 = -qJ(3) * t365 + t303;
t259 = pkin(4) * t374 + t292;
t250 = t272 * t316;
t246 = t266 * t366;
t231 = qJD(1) * t245 + t306;
t229 = pkin(4) * t268 + t258;
t227 = -t267 * t310 + t275;
t218 = t316 * t260 - t330 * t362;
t217 = qJD(2) * t323 + qJD(5) * t249;
t209 = t341 + t214;
t204 = t213 + t321;
t1 = [0.2e1 * t347 * t390 + t389 * t391 + MDP(6) * t372 - MDP(7) * t373 + (-pkin(6) * t372 + t314 * t344) * MDP(9) + (pkin(6) * t373 + t316 * t344) * MDP(10) + t319 * MDP(11) + (t314 * t329 + t316 * t327) * MDP(12) + (-t314 * t327 + t316 * t329) * MDP(13) + (pkin(6) * t319 + t243 * t288 + t256 * t266) * MDP(14) + (t252 * t374 - t268 * t278 + (qJD(1) * t213 + t205) * t314 + (-t258 * t375 + t382 + (t227 * t316 - t292 * t375) * qJD(1)) * qJD(2)) * MDP(15) + (-t316 * t379 - t270 * t278 + (-qJD(1) * t214 - t206) * t314 + (t258 * t376 - t381 + (-t228 * t316 + t292 * t376) * qJD(1)) * qJD(2)) * MDP(16) + (-t213 * t270 - t214 * t268 + (t205 * t310 - t206 * t311) * t316 + ((-t227 * t310 + t228 * t311) * qJD(1) + t334) * t362) * MDP(17) + (t205 * t227 + t206 * t228 + t211 * t213 + t212 * t214 + t252 * t292 - t258 * t278) * MDP(18) + (-t207 * t250 - t217 * t331) * MDP(19) + (t207 * t249 + t208 * t250 - t217 * t224 - t218 * t331) * MDP(20) + (t207 * t314 + t217 * t297 + (-qJD(1) * t250 - t331) * t361) * MDP(21) + (-t208 * t314 + t218 * t297 + (qJD(1) * t249 - t224) * t361) * MDP(22) + (t297 + t366) * MDP(23) * t361 + ((t204 * t315 - t209 * t313) * t297 + t342 * t314 + t245 * t224 + t259 * t208 - t231 * t249 - t229 * t218 + (-t198 * t314 - t297 * t332) * qJD(5) + (qJD(1) * t333 + t197) * t361) * MDP(24) + (-(t204 * t313 + t209 * t315) * t297 - t336 * t314 - t245 * t331 + t259 * t207 - t231 * t250 + t229 * t217 + (-t197 * t314 - t297 * t333) * qJD(5) + (-qJD(1) * t332 - t198) * t361) * MDP(25); -t371 * t390 + t318 * t389 + ((-t290 - t307) * t314 + (-t287 + t343) * t316) * qJD(1) * MDP(11) + (-t276 * t365 + t246) * MDP(12) + (0.2e1 * t306 + (t266 * t316 + t276 * t314) * qJD(1)) * MDP(13) + (-qJ(3) * t286 - qJD(3) * t290 - t266 * t276 + (-t290 * t314 + (-t287 - t383) * t316) * qJD(1) * pkin(6)) * MDP(14) + (t379 + t356 * t268 + (-t221 * t314 - t386 * t311 - t382) * qJD(1)) * MDP(15) + (t252 * t311 + t356 * t270 + (t222 * t314 + t386 * t310 + t381) * qJD(1)) * MDP(16) + (t221 * t270 + t222 * t268 + (qJD(4) * t270 - t212 * t366 - t205) * t311 + (qJD(4) * t268 + t211 * t366 - t206) * t310) * MDP(17) + (qJ(3) * t252 - t211 * t221 - t212 * t222 + t335 * t312 + t356 * t258 + (-t211 * t311 - t212 * t310) * qJD(4)) * MDP(18) + (-t207 * t330 - t331 * t369) * MDP(19) + (-t207 * t272 + t208 * t330 - t224 * t369 + t331 * t368) * MDP(20) + (-t241 * t297 + t331 * t365 + t370) * MDP(21) + (-t368 * t297 + (-qJD(2) * t272 + t224) * t365) * MDP(22) - t297 * t346 + (t298 * t208 + t231 * t272 + (t313 * t324 - t315 * t325) * t297 + t368 * t229 + t357 * t224 + ((-t281 * t313 + t282 * t315) * qJD(2) - t197) * t365) * MDP(24) + (t298 * t207 - t231 * t330 + (t313 * t325 + t315 * t324) * t297 + t369 * t229 - t357 * t331 + (-(t281 * t315 + t282 * t313) * qJD(2) + t198) * t365) * MDP(25) + (MDP(9) * t314 * t318 + MDP(10) * t371) * pkin(1); -t317 * MDP(13) + (t246 + t295) * MDP(14) + t335 * MDP(18) + t370 * MDP(24) + (-MDP(15) * t310 - MDP(16) * t311 - MDP(13)) * t318 * t308 + (-t241 * MDP(24) - MDP(25) * t368) * t297 + (t290 * MDP(14) + (-t268 + t351) * MDP(15) + (-t270 - t349) * MDP(16) - t258 * MDP(18) - t224 * MDP(24) + (-t272 * t365 + t331) * MDP(25)) * qJD(2) + (MDP(12) * t371 + ((-t268 * t311 + t270 * t310) * MDP(17) + t334 * MDP(18)) * qJD(1)) * t314; (-t268 ^ 2 - t270 ^ 2) * MDP(17) + (t211 * t270 + t212 * t268 + t252) * MDP(18) + (t208 - t392) * MDP(24) + (t207 - t380) * MDP(25) + ((t270 - t363) * MDP(15) + (-t268 + t364) * MDP(16)) * t366; -t331 * t224 * MDP(19) + (-t224 ^ 2 + t331 ^ 2) * MDP(20) + (t354 + t380) * MDP(21) + (-t322 - t392) * MDP(22) + qJD(2) * t346 + (t198 * t297 + t229 * t331 + t342) * MDP(24) + (t197 * t297 + t224 * t229 - t336) * MDP(25) + (-MDP(21) * t378 + MDP(22) * t331 - MDP(24) * t198 - MDP(25) * t197) * qJD(5);];
tauc = t1;
