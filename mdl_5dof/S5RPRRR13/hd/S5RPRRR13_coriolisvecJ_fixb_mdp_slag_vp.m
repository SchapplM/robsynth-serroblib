% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:32
% EndTime: 2019-12-31 19:15:39
% DurationCPUTime: 3.19s
% Computational Cost: add. (1616->319), mult. (3605->464), div. (0->0), fcn. (2338->6), ass. (0->145)
t282 = sin(qJ(4));
t285 = cos(qJ(4));
t331 = t285 * qJD(3);
t286 = cos(qJ(3));
t342 = qJD(1) * t286;
t251 = t282 * t342 - t331;
t284 = cos(qJ(5));
t322 = t285 * t342;
t340 = qJD(3) * t282;
t253 = t322 + t340;
t281 = sin(qJ(5));
t361 = t253 * t281;
t210 = t284 * t251 + t361;
t283 = sin(qJ(3));
t343 = qJD(1) * t283;
t274 = qJD(4) + t343;
t271 = qJD(5) + t274;
t378 = t210 * t271;
t301 = t251 * t281 - t284 * t253;
t377 = t271 * t301;
t254 = t281 * t282 - t284 * t285;
t371 = qJD(4) + qJD(5);
t376 = t371 * t254;
t255 = t281 * t285 + t282 * t284;
t290 = t371 * t255;
t260 = pkin(3) * t283 - pkin(7) * t286 + qJ(2);
t240 = t260 * qJD(1);
t287 = -pkin(1) - pkin(6);
t272 = t287 * qJD(1) + qJD(2);
t259 = t283 * t272;
t243 = qJD(3) * pkin(7) + t259;
t207 = t240 * t282 + t243 * t285;
t202 = -pkin(8) * t251 + t207;
t334 = qJD(5) * t281;
t200 = t202 * t334;
t360 = t272 * t286;
t244 = -qJD(3) * pkin(3) - t360;
t217 = pkin(4) * t251 + t244;
t375 = t210 * t217 + t200;
t335 = qJD(4) * t286;
t318 = t282 * t335;
t294 = -t283 * t331 - t318;
t221 = t294 * qJD(1) + qJD(4) * t331;
t305 = pkin(3) * t286 + pkin(7) * t283;
t249 = t305 * qJD(3) + qJD(2);
t231 = t249 * qJD(1);
t225 = t285 * t231;
t291 = -t207 * qJD(4) + t225;
t338 = qJD(3) * t286;
t191 = -pkin(8) * t221 + (pkin(4) * qJD(1) - t272 * t282) * t338 + t291;
t316 = t282 * t343;
t222 = -qJD(3) * t316 + qJD(4) * t253;
t319 = t286 * t331;
t336 = qJD(4) * t285;
t337 = qJD(4) * t282;
t297 = t282 * t231 + t240 * t336 - t243 * t337 + t272 * t319;
t192 = -pkin(8) * t222 + t297;
t312 = t284 * t191 - t281 * t192;
t374 = t217 * t301 + t312;
t314 = MDP(25) * t342;
t373 = qJD(3) * t314 + (-t210 ^ 2 + t301 ^ 2) * MDP(22) - t210 * t301 * MDP(21);
t329 = 0.2e1 * qJD(1);
t280 = t286 ^ 2;
t372 = MDP(8) * (t283 ^ 2 - t280);
t353 = t283 * t287;
t270 = t285 * t353;
t346 = t282 * t260 + t270;
t311 = t221 * t281 + t284 * t222;
t194 = -t301 * qJD(5) + t311;
t370 = pkin(7) + pkin(8);
t289 = qJD(1) ^ 2;
t369 = qJ(2) * t289;
t206 = t285 * t240 - t243 * t282;
t201 = -pkin(8) * t253 + t206;
t199 = pkin(4) * t274 + t201;
t368 = t199 * t284;
t367 = t202 * t284;
t366 = t221 * t282;
t365 = t222 * t286;
t364 = t244 * t282;
t363 = t251 * t274;
t362 = t253 * t274;
t359 = t274 * t282;
t358 = t274 * t283;
t357 = t274 * t285;
t356 = t282 * t286;
t355 = t282 * t287;
t354 = t283 * t285;
t352 = t285 * t286;
t288 = qJD(3) ^ 2;
t350 = t287 * t288;
t298 = t254 * t283;
t349 = -qJD(1) * t298 - t376;
t296 = t255 * qJD(1);
t348 = t283 * t296 + t290;
t256 = t305 * qJD(1);
t347 = t282 * t256 + t272 * t352;
t341 = qJD(3) * t271;
t339 = qJD(3) * t283;
t333 = qJD(5) * t284;
t330 = qJD(1) * qJD(3);
t328 = pkin(8) * t354;
t326 = t272 * t356;
t325 = t282 * t353;
t324 = t284 * t221 - t281 * t222 - t251 * t333;
t323 = qJD(4) * t370;
t321 = t282 * t339;
t317 = t285 * t335;
t315 = t286 * t330;
t313 = pkin(4) - t355;
t310 = t251 + t331;
t309 = -t253 + t340;
t308 = qJD(5) * t199 + t192;
t307 = qJD(4) * t283 + qJD(1);
t269 = t283 * t315;
t306 = -t259 + (t316 + t337) * pkin(4);
t242 = t285 * t256;
t267 = t370 * t285;
t304 = qJD(5) * t267 - t326 + t242 + (pkin(4) * t286 + t328) * qJD(1) + t285 * t323;
t266 = t370 * t282;
t303 = pkin(8) * t316 + qJD(5) * t266 + t282 * t323 + t347;
t189 = t199 * t281 + t367;
t248 = t285 * t260;
t208 = -pkin(8) * t352 + t313 * t283 + t248;
t216 = -pkin(8) * t356 + t346;
t302 = t208 * t281 + t216 * t284;
t300 = qJD(1) * t280 - t358;
t299 = -pkin(7) * t338 + t244 * t283;
t193 = -t253 * t334 + t324;
t295 = t254 * qJD(1);
t293 = -t317 + t321;
t292 = -qJD(4) * t325 + t282 * t249 + t260 * t336 + t287 * t319;
t277 = -pkin(4) * t285 - pkin(3);
t250 = (pkin(4) * t282 - t287) * t286;
t236 = t285 * t249;
t233 = t254 * t286;
t232 = t255 * t286;
t223 = -t293 * pkin(4) + t287 * t339;
t204 = pkin(4) * t222 + t272 * t339;
t198 = -t334 * t356 + (t371 * t352 - t321) * t284 + t294 * t281;
t197 = qJD(3) * t298 - t290 * t286;
t196 = t293 * pkin(8) + t292;
t195 = t236 + (-t270 + (pkin(8) * t286 - t260) * t282) * qJD(4) + (t313 * t286 + t328) * qJD(3);
t188 = -t202 * t281 + t368;
t1 = [-0.2e1 * MDP(7) * t269 + 0.2e1 * t330 * t372 + (-t283 * t350 + (qJ(2) * t338 + qJD(2) * t283) * t329) * MDP(12) + (-t286 * t350 + (-qJ(2) * t339 + qJD(2) * t286) * t329) * MDP(13) + (t221 * t352 + t294 * t253) * MDP(14) + ((t251 * t285 + t253 * t282) * t339 + (-t366 - t222 * t285 + (t251 * t282 - t253 * t285) * qJD(4)) * t286) * MDP(15) + (-t274 * t318 + t221 * t283 + (t253 * t286 + t300 * t285) * qJD(3)) * MDP(16) + (-t274 * t317 - t222 * t283 + (-t251 * t286 - t300 * t282) * qJD(3)) * MDP(17) + (t274 * t338 + t269) * MDP(18) + (-t287 * t365 + t225 * t283 + t236 * t274 + (-t207 * t283 + t244 * t352 - t346 * t274) * qJD(4) + ((t251 * t287 - t364) * t283 + (-t274 * t355 + (t248 - t325) * qJD(1) + t206) * t286) * qJD(3)) * MDP(19) + (-t292 * t274 - t297 * t283 + (-t287 * t221 - t244 * t337) * t286 + ((-t346 * qJD(1) - t207) * t286 + (t287 * t253 + (-t244 + t360) * t285) * t283) * qJD(3)) * MDP(20) + (-t193 * t233 - t197 * t301) * MDP(21) + (-t193 * t232 + t194 * t233 - t197 * t210 + t198 * t301) * MDP(22) + (t193 * t283 + t197 * t271 + (-qJD(1) * t233 - t301) * t338) * MDP(23) + (-t194 * t283 - t198 * t271 + (-qJD(1) * t232 - t210) * t338) * MDP(24) + (t271 * t338 + t269) * MDP(25) + ((t195 * t284 - t196 * t281) * t271 + t312 * t283 + t223 * t210 + t250 * t194 + t204 * t232 + t217 * t198 + (-t189 * t283 - t302 * t271) * qJD(5) + ((t208 * t284 - t216 * t281) * qJD(1) + t188) * t338) * MDP(26) + (t250 * t193 + t217 * t197 + t200 * t283 - t204 * t233 - t223 * t301 + (-(-qJD(5) * t216 + t195) * t271 - t191 * t283) * t281 + (-(qJD(5) * t208 + t196) * t271 - t308 * t283) * t284 + (-t302 * qJD(1) - t189) * t338) * MDP(27) + (qJ(2) * MDP(6) + MDP(5)) * qJD(2) * t329 + (-t286 * MDP(10) - t283 * MDP(9)) * t288; -t289 * MDP(5) - MDP(6) * t369 + (-t365 - t307 * t357 + (t251 * t283 + (-t274 - t343) * t356) * qJD(3)) * MDP(19) + (-t221 * t286 + t307 * t359 + (-t274 * t352 + (t253 - t322) * t283) * qJD(3)) * MDP(20) + (t271 * t295 + (-t255 * t341 - t194) * t286 + ((-t255 * t342 + t210) * qJD(3) + t271 * t376) * t283) * MDP(26) + (t271 * t296 + (t254 * t341 - t193) * t286 + (t290 * t271 + (t286 * t295 - t301) * qJD(3)) * t283) * MDP(27) + (t283 * MDP(12) + t286 * MDP(13)) * (-t288 - t289); t283 * MDP(13) * t369 + (t253 * t357 + t366) * MDP(14) + ((t221 - t363) * t285 + (-t222 - t362) * t282) * MDP(15) + (t274 * t336 + (t274 * t354 + t309 * t286) * qJD(1)) * MDP(16) + (-t274 * t337 + (-t282 * t358 + t310 * t286) * qJD(1)) * MDP(17) - t274 * MDP(18) * t342 + (-pkin(3) * t222 - t242 * t274 + (t274 * t356 - t310 * t283) * t272 + (-pkin(7) * t357 + t364) * qJD(4) + (-t206 * t286 + t299 * t282) * qJD(1)) * MDP(19) + (-pkin(3) * t221 + t347 * t274 + t309 * t259 + (pkin(7) * t359 + t244 * t285) * qJD(4) + (t207 * t286 + t299 * t285) * qJD(1)) * MDP(20) + (t193 * t255 - t301 * t349) * MDP(21) + (-t193 * t254 - t194 * t255 - t349 * t210 + t301 * t348) * MDP(22) + (t349 * t271 + (qJD(3) * t255 + t301) * t342) * MDP(23) + (-t348 * t271 + (-qJD(3) * t254 + t210) * t342) * MDP(24) - t271 * t314 + (t277 * t194 + t204 * t254 + (t303 * t281 - t304 * t284) * t271 + t348 * t217 + t306 * t210 + ((-t266 * t284 - t267 * t281) * qJD(3) - t188) * t342) * MDP(26) + (t277 * t193 + t204 * t255 + (t304 * t281 + t303 * t284) * t271 + t349 * t217 - t306 * t301 + (-(-t266 * t281 + t267 * t284) * qJD(3) + t189) * t342) * MDP(27) + (-t372 + (-qJ(2) * MDP(12) + t283 * MDP(7)) * t286) * t289; t253 * t251 * MDP(14) + (-t251 ^ 2 + t253 ^ 2) * MDP(15) + (t221 + t363) * MDP(16) + (-t222 + t362) * MDP(17) + MDP(18) * t315 + (-qJD(3) * t326 + t207 * t274 - t244 * t253 + t291) * MDP(19) + (t206 * t274 + t244 * t251 - t297) * MDP(20) + (t193 + t378) * MDP(23) + (-t194 - t377) * MDP(24) + (-(-t201 * t281 - t367) * t271 - t189 * qJD(5) + (-t210 * t253 - t271 * t334 + t284 * t315) * pkin(4) + t374) * MDP(26) + ((-t202 * t271 - t191) * t281 + (t201 * t271 - t308) * t284 + (t253 * t301 - t271 * t333 - t281 * t315) * pkin(4) + t375) * MDP(27) + t373; (t324 + t378) * MDP(23) + (-t311 - t377) * MDP(24) + (t189 * t271 + t374) * MDP(26) + (t188 * t271 - t281 * t191 - t284 * t192 + t375) * MDP(27) + (-MDP(23) * t361 + t301 * MDP(24) - t189 * MDP(26) - MDP(27) * t368) * qJD(5) + t373;];
tauc = t1;
