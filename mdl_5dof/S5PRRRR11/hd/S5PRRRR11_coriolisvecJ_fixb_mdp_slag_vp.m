% Calculate Coriolis joint torque vector for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:59
% EndTime: 2024-09-27 21:46:01
% DurationCPUTime: 0.95s
% Computational Cost: add. (1811->240), mult. (5940->343), div. (0->0), fcn. (4564->8), ass. (0->140)
t302 = sin(qJ(3));
t298 = sin(pkin(5));
t383 = pkin(7) + pkin(8);
t346 = t298 * t383;
t330 = t302 * t346;
t305 = cos(qJ(3));
t360 = qJD(1) * t298;
t291 = t305 * t360;
t299 = cos(pkin(5));
t357 = qJD(2) * t305;
t347 = pkin(2) * t357;
t362 = t299 * t347 + t291;
t254 = -qJD(2) * t330 + t362;
t358 = qJD(2) * t299;
t292 = qJD(3) + t358;
t245 = pkin(3) * t292 + t254;
t348 = pkin(2) * t299 * t302;
t369 = t298 * t305;
t269 = t383 * t369 + t348;
t345 = t302 * t360;
t255 = qJD(2) * t269 + t345;
t304 = cos(qJ(4));
t253 = t304 * t255;
t301 = sin(qJ(4));
t321 = -t301 * t245 - t253;
t342 = t298 * t357;
t359 = qJD(2) * t298;
t343 = t302 * t359;
t391 = -t301 * t343 + t304 * t342;
t379 = pkin(9) * t391;
t214 = -t321 + t379;
t300 = sin(qJ(5));
t352 = qJD(5) * t300;
t212 = t214 * t352;
t303 = cos(qJ(5));
t266 = t303 * t391;
t272 = -t301 * t342 - t304 * t343;
t232 = t272 * t300 + t266;
t276 = (-pkin(3) * t305 - pkin(2)) * t298;
t294 = t299 * qJD(1);
t273 = qJD(2) * t276 + t294;
t241 = -pkin(4) * t391 + t273;
t393 = -t241 * t232 + t212;
t386 = qJD(3) + qJD(4);
t234 = t391 * t386;
t323 = qJD(3) * t330;
t355 = qJD(3) * t305;
t344 = t299 * t355;
t325 = qJD(2) * t344;
t363 = pkin(2) * t325 + qJD(3) * t291;
t248 = -qJD(2) * t323 + t363;
t310 = -t305 * t346 - t348;
t249 = (qJD(2) * t310 - t345) * qJD(3);
t337 = -t301 * t248 + t304 * t249;
t309 = t321 * qJD(4) + t337;
t199 = -pkin(9) * t234 + t309;
t288 = qJD(4) + t292;
t284 = qJD(5) + t288;
t392 = (-t214 * t284 - t199) * t300 + t393;
t295 = t298 ^ 2;
t390 = (-t302 * t305 * MDP(5) + (t302 ^ 2 - t305 ^ 2) * MDP(6)) * t295 * qJD(2);
t319 = t303 * t272 - t300 * t391;
t318 = t301 * t305 + t302 * t304;
t384 = t298 * t386;
t247 = t318 * t384;
t235 = qJD(2) * t247;
t353 = qJD(4) * t304;
t354 = qJD(4) * t301;
t328 = -t245 * t353 - t304 * t248 - t301 * t249 + t255 * t354;
t198 = -pkin(9) * t235 - t328;
t340 = -t300 * t198 + t303 * t199;
t315 = t241 * t319 + t340;
t206 = qJD(5) * t266 + t303 * t234 - t300 * t235 + t272 * t352;
t308 = qJD(5) * t319 - t234 * t300 - t303 * t235;
t389 = t232 * t319 * MDP(19) + (-t232 ^ 2 + t319 ^ 2) * MDP(20) + (-t232 * t284 + t206) * MDP(21) + (-t284 * t319 + t308) * MDP(22);
t268 = t272 * pkin(9);
t251 = t301 * t255;
t338 = t304 * t245 - t251;
t213 = t268 + t338;
t388 = -t301 * t302 + t304 * t305;
t350 = qJD(3) - t292;
t387 = qJD(5) - t284;
t382 = pkin(3) * t284;
t381 = pkin(4) * t272;
t380 = pkin(7) * t302;
t377 = t308 * t299;
t274 = t388 * t298;
t275 = t318 * t298;
t239 = -t274 * t303 + t275 * t300;
t246 = t388 * t384;
t209 = -qJD(5) * t239 + t246 * t303 - t247 * t300;
t376 = t209 * t284;
t375 = t235 * t299;
t373 = t246 * t288;
t372 = t273 * t272;
t278 = -pkin(2) * t359 + t294;
t371 = t278 * t305;
t370 = t298 * t302;
t367 = t301 * t303;
t366 = t303 * t214;
t364 = t304 * t254 - t251;
t356 = qJD(3) * t302;
t211 = pkin(4) * t288 + t213;
t341 = -pkin(4) * t284 - t211;
t336 = -t254 * t301 - t253;
t290 = pkin(2) * t344;
t263 = t290 - t323;
t264 = t310 * qJD(3);
t335 = -t263 * t301 + t304 * t264;
t334 = t292 + t358;
t333 = qJD(5) * t211 + t198;
t332 = pkin(3) * t343;
t331 = t298 * pkin(3) * t356;
t324 = t350 * t370;
t322 = -t300 * t211 - t366;
t262 = (pkin(2) * t305 + pkin(3)) * t299 - t330;
t320 = -t262 * t301 - t269 * t304;
t240 = t274 * t300 + t275 * t303;
t316 = t298 * t334;
t313 = -t273 * t391 + t328;
t311 = t262 * t353 + t304 * t263 + t301 * t264 - t269 * t354;
t307 = t272 * t391 * MDP(12) + (-t288 * t391 + t234) * MDP(14) + (-t386 * t359 * t318 - t272 * t288) * MDP(15) + (t272 ^ 2 - t391 ^ 2) * MDP(13) + t389;
t293 = pkin(3) * t304 + pkin(4);
t277 = t298 * t325;
t257 = t332 - t381;
t256 = -pkin(4) * t274 + t276;
t225 = t247 * t288;
t224 = t234 * t299;
t223 = pkin(4) * t247 + t331;
t222 = pkin(4) * t235 + qJD(2) * t331;
t220 = pkin(9) * t274 - t320;
t217 = pkin(4) * t299 - pkin(9) * t275 + t262 * t304 - t269 * t301;
t216 = t268 + t364;
t215 = t336 - t379;
t210 = qJD(5) * t240 + t246 * t300 + t303 * t247;
t208 = t210 * t284;
t205 = t206 * t299;
t204 = -pkin(9) * t246 + t320 * qJD(4) + t335;
t203 = -pkin(9) * t247 + t311;
t1 = [t277 * MDP(11) + (-t225 + t375) * MDP(17) + (t224 - t373) * MDP(18) + (-t208 - t377) * MDP(24) + (t205 - t376) * MDP(25) + (-t305 * t292 * MDP(11) + (-t292 + t358) * MDP(10) * t302) * t298 * qJD(3); (t292 * t298 * t355 + t277) * MDP(7) - MDP(8) * t316 * t356 + (-t290 * t292 - t363 * t299) * MDP(11) + (t234 * t275 - t246 * t272) * MDP(12) + (t234 * t274 - t235 * t275 + t246 * t391 + t247 * t272) * MDP(13) + (t224 + t373) * MDP(14) + (-t225 - t375) * MDP(15) + (t335 * t288 + t337 * t299 + t276 * t235 + t273 * t247 + (t288 * t320 + t299 * t321) * qJD(4) + (-qJD(2) * t274 - t391) * t331) * MDP(17) + (-t311 * t288 + t328 * t299 + t276 * t234 + t273 * t246 + (qJD(2) * t275 - t272) * t331) * MDP(18) + (t206 * t240 - t209 * t319) * MDP(19) + (-t206 * t239 + t209 * t232 + t210 * t319 + t240 * t308) * MDP(20) + (t205 + t376) * MDP(21) + (-t208 + t377) * MDP(22) + ((-t203 * t300 + t204 * t303) * t284 + t340 * t299 - t223 * t232 - t256 * t308 + t222 * t239 + t241 * t210 + ((-t217 * t300 - t220 * t303) * t284 + t322 * t299) * qJD(5)) * MDP(24) + (t256 * t206 + t241 * t209 + t212 * t299 + t222 * t240 - t223 * t319 + (-(-qJD(5) * t220 + t204) * t284 - t199 * t299) * t300 + (-(qJD(5) * t217 + t203) * t284 - t333 * t299) * t303) * MDP(25) + ((-t305 * pkin(7) * t316 + ((t278 - t294) * t298 + (-t299 * t292 + (-t299 ^ 2 - t295) * qJD(2)) * pkin(2)) * t302) * MDP(10) + (-t295 * t347 + (t334 * t380 + t371) * t298) * MDP(11) - 0.2e1 * t390) * qJD(3); -qJD(1) * t324 * MDP(10) + (t362 * t292 + (t350 * t380 - t371) * t359 - t363) * MDP(11) + t307 + (t257 * t319 + (t215 * t284 - t199 - (-qJD(4) - qJD(5)) * t301 * t382) * t300 + ((-pkin(3) * t353 - qJD(5) * t293 + t216) * t284 - t333) * t303 + t393) * MDP(25) + (-(t215 * t303 - t216 * t300) * t284 + t257 * t232 + (-t300 * t304 - t367) * qJD(4) * t382 + ((-pkin(3) * t367 - t293 * t300) * t284 + t322) * qJD(5) + t315) * MDP(24) + (t364 * t288 + (t272 * t343 - t288 * t353) * pkin(3) + t313) * MDP(18) + (-t336 * t288 + t391 * t332 + t372 + (-t253 + (-pkin(3) * t288 - t245) * t301) * qJD(4) + t337) * MDP(17) + t350 * MDP(7) * t342 + ((-t278 * t370 + t350 * (-pkin(7) * t369 - t348)) * MDP(10) - MDP(8) * t324 + t390) * qJD(2); (-t288 * t321 + t309 + t372) * MDP(17) + (t288 * t338 + t313) * MDP(18) + (-(-t213 * t300 - t366) * t284 - t232 * t381 + (t300 * t341 - t366) * qJD(5) + t315) * MDP(24) + (-t319 * t381 + (qJD(5) * t341 + t213 * t284 - t198) * t303 + t392) * MDP(25) + t307; (t387 * t322 + t315) * MDP(24) + ((-t387 * t211 - t198) * t303 + t392) * MDP(25) + t389;];
tauc = t1;
