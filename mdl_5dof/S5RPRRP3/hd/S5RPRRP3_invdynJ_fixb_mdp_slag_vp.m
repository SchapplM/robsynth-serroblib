% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:23
% EndTime: 2022-01-23 09:30:26
% DurationCPUTime: 2.39s
% Computational Cost: add. (1855->279), mult. (3844->349), div. (0->0), fcn. (2524->12), ass. (0->147)
t306 = sin(qJ(3));
t308 = cos(qJ(3));
t303 = sin(pkin(8));
t281 = pkin(1) * t303 + pkin(6);
t383 = pkin(7) + t281;
t345 = t383 * qJD(1);
t232 = qJD(2) * t306 + t308 * t345;
t305 = sin(qJ(4));
t226 = t305 * t232;
t231 = t308 * qJD(2) - t345 * t306;
t382 = qJD(3) * pkin(3);
t229 = t231 + t382;
t387 = cos(qJ(4));
t343 = t229 * t387 - t226;
t254 = t305 * t308 + t306 * t387;
t248 = t254 * qJD(1);
t377 = t248 * qJ(5);
t202 = -t377 + t343;
t298 = qJD(3) + qJD(4);
t291 = t308 * qJDD(2);
t266 = t281 * qJDD(1);
t344 = pkin(7) * qJDD(1) + t266;
t210 = qJDD(3) * pkin(3) - qJD(3) * t232 - t306 * t344 + t291;
t211 = qJD(3) * t231 + t306 * qJDD(2) + t308 * t344;
t395 = t210 * t387 - t305 * t211;
t251 = t383 * t306;
t252 = t383 * t308;
t365 = -t251 * t305 + t252 * t387;
t296 = qJDD(3) + qJDD(4);
t350 = t387 * qJD(4);
t386 = pkin(3) * t298;
t394 = -pkin(3) * t296 * t305 - t350 * t386;
t304 = cos(pkin(8));
t282 = -pkin(1) * t304 - pkin(2);
t295 = t308 * pkin(3);
t393 = t282 - t295;
t288 = t296 * pkin(4);
t370 = t305 * t306;
t332 = t298 * t370;
t352 = t387 * t308;
t338 = qJD(1) * t352;
t347 = qJDD(1) * t387;
t355 = qJDD(1) * t308;
t339 = t298 * t338 + t305 * t355 + t306 * t347;
t214 = qJD(1) * t332 - t339;
t381 = t214 * qJ(5);
t392 = t288 + t381;
t391 = t387 * qJD(3) + t350;
t224 = t298 * t254;
t276 = t308 * t347;
t356 = qJDD(1) * t306;
t331 = t305 * t356 - t276;
t215 = qJD(1) * t224 + t331;
t390 = pkin(4) * t215 + qJDD(5);
t299 = qJ(1) + pkin(8);
t290 = cos(t299);
t302 = qJ(3) + qJ(4);
t293 = sin(t302);
t373 = t290 * t293;
t289 = sin(t299);
t375 = t289 * t293;
t294 = cos(t302);
t384 = g(3) * t294;
t389 = g(1) * t373 + g(2) * t375 - t384;
t388 = t248 ^ 2;
t380 = t215 * qJ(5);
t361 = qJD(1) * t306;
t246 = t305 * t361 - t338;
t379 = t246 * qJ(5);
t378 = t246 * t298;
t374 = t289 * t294;
t372 = t290 * t294;
t369 = qJDD(2) - g(3);
t200 = pkin(4) * t298 + t202;
t368 = t200 - t202;
t223 = -t308 * t391 + t332;
t367 = -t215 * t254 + t223 * t246;
t366 = t231 * t387 - t226;
t364 = pkin(4) * t294 + t295;
t300 = t306 ^ 2;
t363 = -t308 ^ 2 + t300;
t362 = MDP(19) * t306;
t269 = qJD(1) * t282;
t359 = qJD(4) * t305;
t250 = t393 * qJD(1);
t346 = pkin(4) * t246 + qJD(5);
t221 = t250 + t346;
t358 = qJD(5) + t221;
t357 = qJD(1) * qJD(3);
t354 = pkin(3) * t361;
t353 = t306 * t382;
t228 = t387 * t232;
t349 = t306 * t357;
t348 = qJD(3) * t383;
t342 = -t231 * t305 - t228;
t341 = -t251 * t387 - t252 * t305;
t340 = t298 * t306;
t337 = -g(1) * t375 + g(2) * t373;
t336 = g(1) * t374 - g(2) * t372;
t335 = g(1) * t290 + g(2) * t289;
t334 = g(1) * t289 - g(2) * t290;
t307 = sin(qJ(1));
t309 = cos(qJ(1));
t333 = g(1) * t307 - g(2) * t309;
t253 = -t352 + t370;
t329 = -t214 * t253 + t224 * t248;
t328 = t223 * t298 - t254 * t296;
t326 = -t229 * t305 - t228;
t243 = t306 * t348;
t244 = t308 * t348;
t325 = -t243 * t387 - t244 * t305 - t251 * t350 - t252 * t359;
t245 = t246 ^ 2;
t324 = t248 * t246 * MDP(12) + (-qJD(1) * t305 * t340 + t339 + t378) * MDP(14) - t331 * MDP(15) + (-t245 + t388) * MDP(13) + t296 * MDP(16);
t323 = -qJD(1) * t269 - t266 + t335;
t280 = pkin(3) * t349;
t234 = qJDD(1) * t393 + t280;
t322 = 0.2e1 * qJD(3) * t269 - qJDD(3) * t281;
t310 = qJD(3) ^ 2;
t321 = -0.2e1 * qJDD(1) * t282 - t281 * t310 + t334;
t320 = qJD(4) * t326 + t395;
t319 = -qJD(4) * t365 + t243 * t305 - t244 * t387;
t318 = t305 * t210 + t211 * t387 + t229 * t350 - t232 * t359;
t317 = g(1) * t372 + g(2) * t374 + g(3) * t293 - t318;
t316 = t320 + t389;
t315 = t250 * t246 + t317;
t314 = -t250 * t248 + t316;
t313 = t246 * t358 + t317 + t380;
t297 = -qJ(5) - pkin(7) - pkin(6);
t287 = pkin(3) * t387 + pkin(4);
t265 = qJDD(3) * t308 - t306 * t310;
t264 = qJDD(3) * t306 + t308 * t310;
t255 = pkin(2) + t364;
t233 = pkin(4) * t248 + t354;
t230 = pkin(4) * t253 + t393;
t220 = pkin(4) * t224 + t353;
t217 = -qJ(5) * t253 + t365;
t216 = -qJ(5) * t254 + t341;
t213 = -t224 * t298 - t253 * t296;
t205 = -t377 + t366;
t204 = t342 + t379;
t203 = -t326 - t379;
t201 = t234 + t390;
t199 = t223 * qJ(5) - t254 * qJD(5) + t319;
t198 = -qJ(5) * t224 - qJD(5) * t253 + t325;
t197 = -qJD(5) * t246 + t318 - t380;
t196 = -t248 * qJD(5) + t320 + t392;
t1 = [qJDD(1) * MDP(1) + t333 * MDP(2) + (g(1) * t309 + g(2) * t307) * MDP(3) + (t333 + (t303 ^ 2 + t304 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t300 + 0.2e1 * t308 * t349) * MDP(5) + 0.2e1 * (t306 * t355 - t357 * t363) * MDP(6) + t264 * MDP(7) + t265 * MDP(8) + (t306 * t322 + t308 * t321) * MDP(10) + (-t306 * t321 + t308 * t322) * MDP(11) + (-t214 * t254 - t223 * t248) * MDP(12) + (-t329 + t367) * MDP(13) - t328 * MDP(14) + t213 * MDP(15) + (t215 * t393 + t250 * t224 + t234 * t253 + t246 * t353 + t296 * t341 + t298 * t319 + t336) * MDP(17) + (-t214 * t393 - t250 * t223 + t234 * t254 + t248 * t353 - t296 * t365 - t298 * t325 + t337) * MDP(18) + (t199 * t298 + t201 * t253 + t215 * t230 + t216 * t296 + t220 * t246 + t221 * t224 + t336) * MDP(19) + (-t198 * t298 + t201 * t254 - t214 * t230 - t217 * t296 + t220 * t248 - t221 * t223 + t337) * MDP(20) + (-t196 * t254 - t197 * t253 - t198 * t246 - t199 * t248 + t200 * t223 - t203 * t224 + t214 * t216 - t215 * t217 - t335) * MDP(21) + (t197 * t217 + t203 * t198 + t196 * t216 + t200 * t199 + t201 * t230 + t221 * t220 - g(1) * (-pkin(1) * t307 - t255 * t289 - t290 * t297) - g(2) * (pkin(1) * t309 + t255 * t290 - t289 * t297)) * MDP(22); t369 * MDP(4) + t265 * MDP(10) - t264 * MDP(11) + (t329 + t367) * MDP(21) + (-t196 * t253 + t197 * t254 - t200 * t224 - t203 * t223 - g(3)) * MDP(22) + (MDP(17) + MDP(19)) * t213 + (MDP(18) + MDP(20)) * t328; MDP(7) * t356 + MDP(8) * t355 + qJDD(3) * MDP(9) + (-g(3) * t308 + t306 * t323 + t291) * MDP(10) + (-t306 * t369 + t308 * t323) * MDP(11) + (-t342 * t298 + (-t246 * t361 + t296 * t387 - t298 * t359) * pkin(3) + t314) * MDP(17) + (-t248 * t354 + t298 * t366 + t315 + t394) * MDP(18) + (-t204 * t298 - t233 * t246 + t287 * t296 - t358 * t248 + (-t228 + (-t229 - t386) * t305) * qJD(4) + t389 + t392 + t395) * MDP(19) + (t205 * t298 - t233 * t248 + t313 + t394) * MDP(20) + (t287 * t214 + (t203 + t204) * t248 + (-t200 + t205) * t246 + (-t215 * t305 + (-t246 * t387 + t248 * t305) * qJD(4)) * pkin(3)) * MDP(21) + (t196 * t287 - t203 * t205 - t200 * t204 - t221 * t233 - g(3) * t364 - t335 * (-pkin(3) * t306 - pkin(4) * t293) + (t197 * t305 + (-t200 * t305 + t203 * t387) * qJD(4)) * pkin(3)) * MDP(22) + t324 + (-t306 * t308 * MDP(5) + MDP(6) * t363) * qJD(1) ^ 2; (-t298 * t326 + t314) * MDP(17) + (t298 * t343 + t315) * MDP(18) + (t381 + t203 * t298 + 0.2e1 * t288 + (-t221 - t346) * t248 + t316) * MDP(19) + (-pkin(4) * t388 + t202 * t298 + t313) * MDP(20) + (pkin(4) * t214 - t246 * t368) * MDP(21) + (t368 * t203 + (-t221 * t248 + t293 * t335 + t196 - t384) * pkin(4)) * MDP(22) + t324; (t248 * t298 - t276) * MDP(19) + (t339 - t378) * MDP(20) + (-t245 - t388) * MDP(21) + (t200 * t248 + t203 * t246 + t280 - t334 + t390) * MDP(22) + (MDP(22) * t393 + t305 * t362) * qJDD(1) + (t391 * t362 + (MDP(19) * t298 * t308 - MDP(20) * t340) * t305) * qJD(1);];
tau = t1;
