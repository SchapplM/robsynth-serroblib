% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:51
% EndTime: 2019-03-09 03:22:57
% DurationCPUTime: 3.17s
% Computational Cost: add. (3050->322), mult. (6528->442), div. (0->0), fcn. (4311->6), ass. (0->143)
t314 = sin(qJ(3));
t316 = cos(qJ(3));
t391 = sin(pkin(9));
t392 = cos(pkin(9));
t322 = t314 * t392 + t316 * t391;
t275 = t322 * qJD(1);
t399 = qJD(5) + t275;
t402 = pkin(5) * t399;
t283 = -t391 * t314 + t392 * t316;
t278 = t283 * qJD(1);
t315 = cos(qJ(5));
t352 = t315 * qJD(3);
t313 = sin(qJ(5));
t354 = qJD(5) * t313;
t339 = qJD(3) * t391;
t334 = qJD(1) * t339;
t340 = qJD(3) * t392;
t335 = qJD(1) * t340;
t364 = t314 * t335 + t316 * t334;
t225 = -qJD(5) * t352 + t278 * t354 + t315 * t364;
t257 = qJD(3) * t313 + t278 * t315;
t269 = t314 * t334 - t316 * t335;
t317 = -pkin(1) - pkin(7);
t292 = qJD(1) * t317 + qJD(2);
t359 = qJD(1) * t316;
t273 = -qJ(4) * t359 + t316 * t292;
t268 = qJD(3) * pkin(3) + t273;
t360 = qJD(1) * t314;
t272 = -qJ(4) * t360 + t292 * t314;
t343 = t392 * t272;
t231 = t391 * t268 + t343;
t224 = qJD(3) * pkin(8) + t231;
t286 = pkin(3) * t360 + qJD(1) * qJ(2) + qJD(4);
t232 = pkin(4) * t275 - pkin(8) * t278 + t286;
t212 = t224 * t315 + t232 * t313;
t356 = qJD(4) * t314;
t357 = qJD(3) * t316;
t253 = t292 * t357 + (-qJ(4) * t357 - t356) * qJD(1);
t355 = qJD(4) * t316;
t358 = qJD(3) * t314;
t320 = -t292 * t358 + (qJ(4) * t358 - t355) * qJD(1);
t217 = t253 * t392 + t320 * t391;
t309 = qJD(1) * qJD(2);
t349 = qJD(1) * qJD(3);
t344 = t316 * t349;
t363 = pkin(3) * t344 + t309;
t221 = -t269 * pkin(4) + pkin(8) * t364 + t363;
t219 = t315 * t221;
t321 = -qJD(5) * t212 - t217 * t313 + t219;
t200 = -pkin(5) * t269 + qJ(6) * t225 - qJD(6) * t257 + t321;
t255 = t278 * t313 - t352;
t207 = -qJ(6) * t255 + t212;
t401 = t207 * t399 + t200;
t338 = t313 * t364;
t226 = qJD(5) * t257 - t338;
t353 = qJD(5) * t315;
t324 = t315 * t217 + t313 * t221 - t224 * t354 + t232 * t353;
t201 = -qJ(6) * t226 - qJD(6) * t255 + t324;
t211 = -t224 * t313 + t315 * t232;
t206 = -qJ(6) * t257 + t211;
t204 = t206 + t402;
t400 = -t204 * t399 + t201;
t398 = MDP(8) * (t314 ^ 2 - t316 ^ 2);
t395 = MDP(6) * qJ(2) + MDP(5);
t333 = t204 * t315 + t207 * t313;
t381 = t257 * t315;
t384 = t255 * t313;
t394 = -MDP(23) * (t381 + t384) + MDP(24) * t333 + t399 * (MDP(21) * t315 - MDP(22) * t313);
t393 = t257 ^ 2;
t216 = t253 * t391 - t392 * t320;
t389 = t216 * t283;
t388 = t225 * t313;
t387 = t226 * t315;
t265 = t391 * t272;
t230 = t268 * t392 - t265;
t277 = -t314 * t340 - t316 * t339;
t386 = t230 * t277;
t385 = t255 * t275;
t383 = t255 * t315;
t382 = t257 * t313;
t380 = t275 * t313;
t379 = t283 * t313;
t378 = t283 * t315;
t377 = t313 * t269;
t373 = qJ(4) - t317;
t287 = t373 * t314;
t288 = t373 * t316;
t251 = -t287 * t392 - t288 * t391;
t248 = t315 * t251;
t263 = t315 * t269;
t319 = qJD(1) ^ 2;
t376 = t316 * t319;
t318 = qJD(3) ^ 2;
t375 = t317 * t318;
t374 = t314 * pkin(3) + qJ(2);
t301 = pkin(3) * t391 + pkin(8);
t372 = qJ(6) + t301;
t371 = t204 - t206;
t370 = -t313 * t226 - t255 * t353;
t240 = t273 * t392 - t265;
t242 = pkin(3) * t359 + pkin(4) * t278 + pkin(8) * t275;
t369 = t315 * t240 + t313 * t242;
t368 = -t380 * t399 - t263;
t247 = pkin(4) * t322 - pkin(8) * t283 + t374;
t367 = t313 * t247 + t248;
t337 = qJD(5) * t372;
t366 = -qJ(6) * t380 + qJD(6) * t315 - t313 * t337 - t369;
t238 = t315 * t242;
t365 = -pkin(5) * t278 - t238 + (-qJ(6) * t275 - t337) * t315 + (-qJD(6) + t240) * t313;
t351 = pkin(3) * t357 + qJD(2);
t348 = 0.2e1 * qJD(1);
t270 = t358 * t373 - t355;
t271 = -qJD(3) * t288 - t356;
t236 = t270 * t391 + t271 * t392;
t276 = t314 * t339 - t316 * t340;
t241 = -pkin(4) * t276 - pkin(8) * t277 + t351;
t346 = t315 * t236 + t313 * t241 + t247 * t353;
t345 = t283 * t353;
t336 = t399 * t315;
t302 = -pkin(3) * t392 - pkin(4);
t235 = -t392 * t270 + t271 * t391;
t239 = t273 * t391 + t343;
t250 = -t287 * t391 + t392 * t288;
t332 = t204 * t313 - t207 * t315;
t331 = t251 * t269 + t389;
t329 = -qJ(6) * t277 - qJD(6) * t283;
t327 = MDP(21) * t313 + MDP(22) * t315;
t326 = t277 * t313 + t345;
t325 = t277 * t315 - t283 * t354;
t223 = -qJD(3) * pkin(4) - t230;
t323 = t223 * t399 + t301 * t269;
t209 = t226 * pkin(5) + t216;
t281 = t372 * t315;
t280 = t372 * t313;
t254 = t255 ^ 2;
t245 = t315 * t247;
t234 = t315 * t241;
t214 = t255 * pkin(5) + qJD(6) + t223;
t213 = -qJ(6) * t379 + t367;
t210 = pkin(5) * t322 - qJ(6) * t378 - t251 * t313 + t245;
t203 = -qJ(6) * t345 + (-qJD(5) * t251 + t329) * t313 + t346;
t202 = -pkin(5) * t276 - t236 * t313 + t234 + t329 * t315 + (-t248 + (qJ(6) * t283 - t247) * t313) * qJD(5);
t1 = [0.2e1 * t349 * t398 - t318 * t316 * MDP(10) + qJ(2) * t357 * t348 * MDP(12) + (-t316 * t375 + (-qJ(2) * t358 + qJD(2) * t316) * t348) * MDP(13) + (-t217 * t322 + t231 * t276 + t235 * t278 - t236 * t275 - t250 * t364 + t331 - t386) * MDP(14) + (t216 * t250 + t217 * t251 - t230 * t235 + t231 * t236 + t286 * t351 + t363 * t374) * MDP(15) + (-t225 * t378 + t257 * t325) * MDP(16) + ((-t382 - t383) * t277 + (t388 - t387 + (-t381 + t384) * qJD(5)) * t283) * MDP(17) + (-t225 * t322 - t257 * t276 - t263 * t283 + t325 * t399) * MDP(18) + (-t226 * t322 + t255 * t276 + t283 * t377 - t326 * t399) * MDP(19) + (-t269 * t322 - t276 * t399) * MDP(20) + ((-t251 * t353 + t234) * t399 - t245 * t269 + (-t224 * t353 + t219) * t322 - t211 * t276 + t235 * t255 + t250 * t226 + t223 * t345 + ((-qJD(5) * t247 - t236) * t399 + (-qJD(5) * t232 - t217) * t322 + t223 * t277 + t331) * t313) * MDP(21) + (-(-t251 * t354 + t346) * t399 + t367 * t269 - t324 * t322 + t212 * t276 + t235 * t257 - t250 * t225 + t216 * t378 + t325 * t223) * MDP(22) + (-t202 * t257 - t203 * t255 + t210 * t225 - t213 * t226 - t333 * t277 + (qJD(5) * t332 - t200 * t315 - t201 * t313) * t283) * MDP(23) + (t201 * t213 + t207 * t203 + t200 * t210 + t204 * t202 + t209 * (pkin(5) * t379 + t250) + t214 * (pkin(5) * t326 + t235)) * MDP(24) + 0.2e1 * t395 * t309 + (-0.2e1 * MDP(7) * t344 - t318 * MDP(9) + (qJD(2) * t348 - t375) * MDP(12)) * t314; (-t277 * t278 + t283 * t364) * MDP(14) + (t386 - t389) * MDP(15) + (-t226 * t283 - t255 * t277) * MDP(21) + (t225 * t283 - t257 * t277) * MDP(22) + (-t209 * t283 - t214 * t277) * MDP(24) - t395 * t319 + (MDP(14) * t275 - t231 * MDP(15) + (-t382 + t383) * MDP(23) + t332 * MDP(24) + t327 * t399) * t276 + (-t286 * MDP(15) - t394) * qJD(1) - (-t217 * MDP(15) + (t387 + t388) * MDP(23) + (t200 * t313 - t201 * t315) * MDP(24) + (-MDP(14) - t327) * t269 + t394 * qJD(5)) * t322 + (MDP(12) * t314 + MDP(13) * t316) * (-t318 - t319); t314 * MDP(7) * t376 - t319 * t398 + ((t231 - t239) * t278 - (-t240 + t230) * t275 + (t269 * t391 + t364 * t392) * pkin(3)) * MDP(14) + (t230 * t239 - t231 * t240 + (-t216 * t392 + t217 * t391 - t286 * t359) * pkin(3)) * MDP(15) + (t257 * t336 - t388) * MDP(16) + ((-t225 - t385) * t315 - t399 * t382 + t370) * MDP(17) + (-t257 * t278 + t336 * t399 - t377) * MDP(18) + (t255 * t278 - t354 * t399 + t368) * MDP(19) - t399 * t278 * MDP(20) + (-t211 * t278 - t216 * t315 + t302 * t226 - t239 * t255 + (-t301 * t353 - t238) * t399 + (t240 * t399 + t323) * t313) * MDP(21) + (t212 * t278 + t216 * t313 - t302 * t225 - t239 * t257 + (t301 * t354 + t369) * t399 + t323 * t315) * MDP(22) + (-t225 * t280 - t226 * t281 - t366 * t255 - t365 * t257 - t401 * t313 + t400 * t315) * MDP(23) + (t201 * t281 - t200 * t280 + t209 * (-t315 * pkin(5) + t302) + (t313 * t402 - t239) * t214 + t366 * t207 + t365 * t204) * MDP(24) + (MDP(13) * t314 * t319 - MDP(12) * t376) * qJ(2); -t275 ^ 2 * MDP(14) + (t231 * t275 + t363) * MDP(15) + t368 * MDP(21) + t370 * MDP(23) + (-MDP(14) * t278 + MDP(15) * t230 - MDP(21) * t255 - MDP(22) * t257 - MDP(24) * t214) * t278 + (t269 * MDP(22) + t400 * MDP(24) + (-MDP(21) * qJD(5) + MDP(23) * t257) * t399) * t313 + ((t225 - t385) * MDP(23) + t401 * MDP(24) - t399 ^ 2 * MDP(22)) * t315; t257 * t255 * MDP(16) + (-t254 + t393) * MDP(17) + (t255 * t399 - t225) * MDP(18) + (t338 + (-qJD(5) + t399) * t257) * MDP(19) - t269 * MDP(20) + (t212 * t399 - t223 * t257 + t321) * MDP(21) + (t211 * t399 + t223 * t255 - t324) * MDP(22) + (pkin(5) * t225 - t255 * t371) * MDP(23) + (t371 * t207 + (-t214 * t257 + t200) * pkin(5)) * MDP(24); (-t254 - t393) * MDP(23) + (t204 * t257 + t207 * t255 + t209) * MDP(24);];
tauc  = t1;
