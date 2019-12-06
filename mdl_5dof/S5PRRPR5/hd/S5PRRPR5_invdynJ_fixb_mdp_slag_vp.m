% Calculate vector of inverse dynamics joint torques for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:28:08
% EndTime: 2019-12-05 16:28:14
% DurationCPUTime: 3.20s
% Computational Cost: add. (1829->337), mult. (4433->475), div. (0->0), fcn. (3520->14), ass. (0->155)
t403 = qJ(4) + pkin(7);
t313 = sin(pkin(10));
t316 = cos(pkin(10));
t321 = sin(qJ(3));
t324 = cos(qJ(3));
t286 = t313 * t324 + t316 * t321;
t279 = t286 * qJD(3);
t368 = qJDD(2) * t324;
t369 = qJDD(2) * t321;
t344 = -t313 * t369 + t316 * t368;
t240 = qJD(2) * t279 + qJDD(5) - t344;
t323 = cos(qJ(5));
t237 = t323 * t240;
t375 = qJD(2) * t321;
t384 = t316 * t324;
t278 = qJD(2) * t384 - t313 * t375;
t272 = qJD(5) - t278;
t320 = sin(qJ(5));
t374 = qJD(5) * t320;
t409 = t272 * t374 - t237;
t318 = cos(pkin(5));
t367 = t318 * qJDD(1);
t296 = t324 * t367;
t315 = sin(pkin(5));
t322 = sin(qJ(2));
t325 = cos(qJ(2));
t371 = qJD(1) * qJD(2);
t266 = qJDD(2) * pkin(7) + (qJDD(1) * t322 + t325 * t371) * t315;
t377 = qJD(1) * t318;
t330 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t377 + t266;
t376 = qJD(1) * t322;
t364 = t315 * t376;
t351 = t403 * qJD(2) + t364;
t342 = t351 * qJD(3);
t218 = qJDD(3) * pkin(3) - t321 * t330 - t324 * t342 + t296;
t219 = (-t342 + t367) * t321 + t330 * t324;
t210 = t218 * t316 - t219 * t313;
t208 = -qJDD(3) * pkin(4) - t210;
t314 = sin(pkin(9));
t317 = cos(pkin(9));
t383 = t318 * t322;
t275 = t314 * t325 + t317 * t383;
t277 = -t314 * t383 + t317 * t325;
t280 = t286 * qJD(2);
t301 = pkin(3) * t313 + pkin(8);
t310 = qJ(3) + pkin(10);
t305 = sin(t310);
t306 = cos(t310);
t387 = t315 * t322;
t388 = t315 * t317;
t389 = t314 * t315;
t408 = t272 * (pkin(3) * t375 + pkin(4) * t280 - pkin(8) * t278 + qJD(5) * t301) + g(1) * (-t277 * t305 + t306 * t389) + g(2) * (-t275 * t305 - t306 * t388) + g(3) * (-t305 * t387 + t306 * t318) + t208;
t211 = t313 * t218 + t316 * t219;
t209 = qJDD(3) * pkin(8) + t211;
t256 = -t351 * t321 + t324 * t377;
t401 = qJD(3) * pkin(3);
t253 = t256 + t401;
t257 = t321 * t377 + t324 * t351;
t395 = t257 * t313;
t224 = t253 * t316 - t395;
t220 = -qJD(3) * pkin(4) - t224;
t304 = pkin(3) * t324 + pkin(2);
t385 = t315 * t325;
t363 = qJD(1) * t385;
t271 = -t304 * qJD(2) + qJD(4) - t363;
t230 = -pkin(4) * t278 - pkin(8) * t280 + t271;
t285 = t313 * t321 - t384;
t242 = pkin(4) * t285 - pkin(8) * t286 - t304;
t290 = t403 * t324;
t360 = t403 * t321;
t261 = t316 * t290 - t313 * t360;
t282 = t285 * qJD(3);
t347 = g(1) * t277 + g(2) * t275;
t336 = -g(3) * t387 - t347;
t355 = qJD(3) * t403;
t273 = qJD(4) * t324 - t321 * t355;
t338 = -qJD(4) * t321 - t324 * t355;
t379 = t316 * t273 + t285 * t363 + t313 * t338;
t407 = -(qJD(5) * t230 + t209) * t285 + t208 * t286 - t220 * t282 + (-qJD(5) * t242 - t379) * t272 - t261 * t240 + t336;
t382 = t318 * t325;
t274 = t314 * t322 - t317 * t382;
t276 = t314 * t382 + t317 * t322;
t348 = g(1) * t276 + g(2) * t274;
t406 = -(g(3) * t385 - t348) * t306 + t242 * t240;
t359 = t322 * t371;
t294 = t315 * t359;
t326 = qJD(3) ^ 2;
t356 = qJDD(1) * t385;
t405 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t326 + t315 * (-g(3) * t325 + t359) - t294 + t348 + t356;
t283 = t318 * t324 - t321 * t387;
t404 = g(3) * t283;
t402 = qJD(2) * pkin(2);
t399 = t220 * t286;
t370 = qJD(2) * qJD(3);
t357 = t324 * t370;
t358 = t321 * t370;
t245 = qJDD(2) * t286 - t313 * t358 + t316 * t357;
t372 = t323 * qJD(3);
t365 = qJD(5) * t372 + t320 * qJDD(3) + t323 * t245;
t222 = -t280 * t374 + t365;
t398 = t222 * t320;
t397 = t240 * t320;
t390 = t280 * t320;
t262 = -t372 + t390;
t394 = t262 * t272;
t393 = t262 * t280;
t264 = qJD(3) * t320 + t280 * t323;
t392 = t264 * t272;
t391 = t264 * t280;
t386 = t315 * t324;
t251 = t316 * t257;
t381 = qJDD(1) - g(3);
t380 = t273 * t313 - t286 * t363 - t316 * t338;
t225 = t313 * t253 + t251;
t311 = t321 ^ 2;
t378 = -t324 ^ 2 + t311;
t373 = qJD(5) * t323;
t366 = t321 * t401;
t362 = qJD(2) * t385;
t354 = t315 * t381;
t353 = -t323 * qJDD(3) + t245 * t320;
t352 = t272 * t323;
t346 = g(1) * t314 - g(2) * t317;
t345 = pkin(4) * t279 + pkin(8) * t282 - t364 + t366;
t221 = qJD(3) * pkin(8) + t225;
t213 = t221 * t323 + t230 * t320;
t343 = t221 * t320 - t230 * t323;
t341 = t320 * t278 * t272 - t409;
t284 = t318 * t321 + t322 * t386;
t339 = -t282 * t323 - t286 * t374;
t289 = -t363 - t402;
t334 = -qJD(2) * t289 - t266 + t347;
t333 = pkin(3) * t358 - t304 * qJDD(2) + qJDD(4) + t294;
t229 = t256 * t316 - t395;
t332 = -t301 * t240 + (t220 + t229) * t272;
t329 = -pkin(7) * qJDD(3) + (t289 + t363 - t402) * qJD(3);
t243 = t333 - t356;
t327 = qJD(2) ^ 2;
t302 = -pkin(3) * t316 - pkin(4);
t269 = t305 * t318 + t306 * t387;
t260 = t290 * t313 + t316 * t360;
t255 = -qJD(3) * t284 - t321 * t362;
t254 = qJD(3) * t283 + t324 * t362;
t249 = t277 * t306 + t305 * t389;
t247 = t275 * t306 - t305 * t388;
t244 = -qJD(3) * t280 + t344;
t236 = t283 * t313 + t284 * t316;
t235 = -t316 * t283 + t284 * t313;
t228 = t254 * t316 + t255 * t313;
t227 = t256 * t313 + t251;
t226 = t254 * t313 - t316 * t255;
t223 = qJD(5) * t264 + t353;
t215 = -pkin(4) * t244 - pkin(8) * t245 + t243;
t214 = t323 * t215;
t1 = [t381 * MDP(1) + (qJD(3) * t255 + qJDD(3) * t283) * MDP(10) + (-qJD(3) * t254 - qJDD(3) * t284) * MDP(11) + (t226 * t280 + t228 * t278 + t235 * t245 + t236 * t244) * MDP(12) + (-t210 * t235 + t211 * t236 - t224 * t226 + t225 * t228 - g(3)) * MDP(13) + ((-t228 * t320 - t236 * t373) * t272 - t236 * t397 + t226 * t262 + t235 * t223) * MDP(19) + (-(t228 * t323 - t236 * t374) * t272 - t236 * t237 + t226 * t264 + t235 * t222) * MDP(20) + ((-qJDD(2) * MDP(4) + (-t324 * MDP(10) + t321 * MDP(11) - MDP(3)) * t327 + (MDP(13) * t271 + (MDP(19) * t323 - MDP(20) * t320) * t272) * qJD(2)) * t322 + (qJDD(2) * MDP(3) - t327 * MDP(4) + (-t358 + t368) * MDP(10) + (-t357 - t369) * MDP(11) - t243 * MDP(13) + t409 * MDP(19) + (t272 * t373 + t397) * MDP(20)) * t325) * t315; qJDD(2) * MDP(2) + (t381 * t385 + t348) * MDP(3) + (-t322 * t354 + t347) * MDP(4) + (qJDD(2) * t311 + 0.2e1 * t321 * t357) * MDP(5) + 0.2e1 * (t321 * t368 - t378 * t370) * MDP(6) + (qJDD(3) * t321 + t324 * t326) * MDP(7) + (qJDD(3) * t324 - t321 * t326) * MDP(8) + (t329 * t321 + t324 * t405) * MDP(10) + (-t321 * t405 + t329 * t324) * MDP(11) + (-t210 * t286 - t211 * t285 + t224 * t282 - t225 * t279 + t244 * t261 + t245 * t260 + t379 * t278 + t380 * t280 + t336) * MDP(12) + (t211 * t261 - t210 * t260 - t243 * t304 + t271 * t366 - g(1) * (-t276 * t304 + t277 * t403) - g(2) * (-t274 * t304 + t275 * t403) + t379 * t225 - t380 * t224 + (-t271 * t376 - g(3) * (t304 * t325 + t322 * t403)) * t315) * MDP(13) + (t222 * t286 * t323 + t264 * t339) * MDP(14) + (-(-t262 * t323 - t264 * t320) * t282 + (-t398 - t223 * t323 + (t262 * t320 - t264 * t323) * qJD(5)) * t286) * MDP(15) + (t222 * t285 + t286 * t237 + t264 * t279 + t272 * t339) * MDP(16) + (-t286 * t397 - t223 * t285 - t262 * t279 + (t282 * t320 - t286 * t373) * t272) * MDP(17) + (t240 * t285 + t272 * t279) * MDP(18) + (-t343 * t279 + t214 * t285 + t260 * t223 + t380 * t262 + (t345 * t272 + (-t221 * t285 - t261 * t272 + t399) * qJD(5) + t406) * t323 + t407 * t320) * MDP(19) + (-t213 * t279 + t260 * t222 + t380 * t264 + (-(-qJD(5) * t221 + t215) * t285 - qJD(5) * t399 + (qJD(5) * t261 - t345) * t272 - t406) * t320 + t407 * t323) * MDP(20); MDP(7) * t369 + MDP(8) * t368 + qJDD(3) * MDP(9) + (t321 * t334 - t346 * t386 + t296 - t404) * MDP(10) + (g(3) * t284 + (t315 * t346 - t367) * t321 + t334 * t324) * MDP(11) + ((t225 - t227) * t280 + (t224 - t229) * t278 + (t244 * t313 - t245 * t316) * pkin(3)) * MDP(12) + (t224 * t227 - t225 * t229 + (t211 * t313 + t210 * t316 - t271 * t375 - g(1) * (-t277 * t321 + t314 * t386) - g(2) * (-t275 * t321 - t317 * t386) - t404) * pkin(3)) * MDP(13) + (t264 * t352 + t398) * MDP(14) + ((t222 - t394) * t323 + (-t223 - t392) * t320) * MDP(15) + (t272 * t352 - t391 + t397) * MDP(16) + (t341 + t393) * MDP(17) - t272 * t280 * MDP(18) + (t302 * t223 - t227 * t262 + t280 * t343 + t332 * t320 - t323 * t408) * MDP(19) + (t213 * t280 + t302 * t222 - t227 * t264 + t320 * t408 + t332 * t323) * MDP(20) + (-t321 * t324 * MDP(5) + t378 * MDP(6)) * t327; (-t278 ^ 2 - t280 ^ 2) * MDP(12) + (t224 * t280 - t225 * t278 - t325 * t354 + t333 - t348) * MDP(13) + (t341 - t393) * MDP(19) + (-t272 ^ 2 * t323 - t391 - t397) * MDP(20); t264 * t262 * MDP(14) + (-t262 ^ 2 + t264 ^ 2) * MDP(15) + (t365 + t394) * MDP(16) + (-t353 + t392) * MDP(17) + t240 * MDP(18) + (-t320 * t209 + t214 + t213 * t272 - t220 * t264 - g(1) * (-t249 * t320 + t276 * t323) - g(2) * (-t247 * t320 + t274 * t323) - g(3) * (-t269 * t320 - t323 * t385)) * MDP(19) + (-t323 * t209 - t320 * t215 - t343 * t272 + t220 * t262 - g(1) * (-t249 * t323 - t276 * t320) - g(2) * (-t247 * t323 - t274 * t320) - g(3) * (-t269 * t323 + t320 * t385)) * MDP(20) + (-MDP(16) * t390 - MDP(17) * t264 - MDP(19) * t213 + MDP(20) * t343) * qJD(5);];
tau = t1;
