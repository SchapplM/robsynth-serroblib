% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:54:14
% EndTime: 2020-01-03 11:54:17
% DurationCPUTime: 1.44s
% Computational Cost: add. (1347->208), mult. (2211->279), div. (0->0), fcn. (1414->14), ass. (0->129)
t385 = pkin(7) + pkin(8);
t299 = sin(qJ(5));
t300 = sin(qJ(4));
t303 = cos(qJ(5));
t304 = cos(qJ(4));
t242 = t299 * t304 + t300 * t303;
t293 = qJD(1) + qJD(3);
t236 = t242 * t293;
t301 = sin(qJ(3));
t305 = cos(qJ(3));
t298 = cos(pkin(9));
t277 = pkin(1) * t298 + pkin(2);
t263 = t277 * qJD(1);
t297 = sin(pkin(9));
t379 = pkin(1) * t297;
t383 = qJD(3) * t263 + qJDD(1) * t379;
t349 = qJD(1) * t379;
t384 = -qJD(3) * t349 + t277 * qJDD(1);
t334 = -t383 * t301 + t384 * t305;
t291 = qJDD(1) + qJDD(3);
t378 = pkin(3) * t291;
t211 = -t334 - t378;
t284 = qJ(1) + pkin(9) + qJ(3);
t275 = sin(t284);
t276 = cos(t284);
t332 = -g(2) * t276 - g(3) * t275;
t322 = -t211 + t332;
t382 = g(2) * t275 - g(3) * t276;
t358 = t303 * t304;
t362 = t299 * t300;
t241 = -t358 + t362;
t338 = t277 * t305 - t301 * t379;
t356 = t301 * t277 + t305 * t379;
t230 = t263 * t301 + t305 * t349;
t341 = t385 * t293 + t230;
t215 = t304 * qJD(2) - t341 * t300;
t292 = qJD(4) + qJD(5);
t216 = qJD(2) * t300 + t341 * t304;
t381 = t384 * t301 + t383 * t305;
t377 = pkin(3) * t293;
t376 = pkin(4) * t304;
t238 = pkin(7) + t356;
t372 = -pkin(8) - t238;
t371 = t216 * t303;
t229 = t263 * t305 - t301 * t349;
t370 = t229 * t292;
t369 = t230 * t293;
t233 = t356 * qJD(3);
t368 = t233 * t293;
t296 = qJ(4) + qJ(5);
t286 = sin(t296);
t367 = t275 * t286;
t366 = t276 * t286;
t364 = t291 * t304;
t363 = t293 * t300;
t359 = t300 * t304;
t357 = qJDD(2) - g(1);
t294 = t300 ^ 2;
t355 = -t304 ^ 2 + t294;
t352 = qJD(4) * t300;
t351 = qJD(4) * t304;
t350 = qJD(5) * t299;
t348 = pkin(4) * t352;
t346 = t293 * t362;
t345 = t293 * t358;
t281 = -pkin(3) - t376;
t344 = qJD(4) * t385;
t343 = t293 * t351;
t210 = pkin(7) * t291 + t381;
t342 = pkin(8) * t291 + t210;
t340 = qJD(4) * t372;
t199 = (t293 * t352 - t364) * pkin(4) + t211;
t217 = t281 * t293 - t229;
t221 = t292 * t241;
t336 = g(2) * t366 + g(3) * t367 + t199 * t242 - t217 * t221;
t223 = -t229 - t377;
t335 = t223 * t351 - t322 * t300;
t237 = -pkin(3) - t338;
t333 = -t230 + t348;
t302 = sin(qJ(1));
t306 = cos(qJ(1));
t331 = -g(2) * t306 - g(3) * t302;
t330 = t241 * t291;
t214 = qJD(4) * pkin(4) + t215;
t328 = -t214 * t299 - t371;
t290 = qJDD(4) + qJDD(5);
t327 = -t221 * t292 + t242 * t290;
t225 = t372 * t300;
t288 = t304 * pkin(8);
t226 = t238 * t304 + t288;
t326 = t225 * t303 - t226 * t299;
t325 = t225 * t299 + t226 * t303;
t268 = t385 * t300;
t269 = pkin(7) * t304 + t288;
t324 = -t268 * t303 - t269 * t299;
t323 = -t268 * t299 + t269 * t303;
t307 = qJD(4) ^ 2;
t320 = pkin(7) * t307 - t369 - t378;
t319 = t237 * t291 + t238 * t307 + t368;
t202 = qJD(5) * t345 + t242 * t291 - t292 * t346 + t303 * t343;
t234 = -t345 + t346;
t318 = t236 * t234 * MDP(15) + (t234 * t292 + t202) * MDP(17) - t330 * MDP(18) + (-t234 ^ 2 + t236 ^ 2) * MDP(16) + t290 * MDP(19);
t317 = -pkin(7) * qJDD(4) + (t229 - t377) * qJD(4);
t316 = -t223 * t293 - t210 + t382;
t315 = t332 + t334;
t232 = t338 * qJD(3);
t314 = -qJDD(4) * t238 + (t237 * t293 - t232) * qJD(4);
t222 = t292 * t242;
t287 = cos(t296);
t313 = t199 * t241 + t217 * t222 + t332 * t287;
t203 = t222 * t293 + t330;
t212 = -t222 * t292 - t241 * t290;
t253 = qJDD(4) * t300 + t304 * t307;
t254 = qJDD(4) * t304 - t300 * t307;
t312 = (-t202 * t241 - t203 * t242 + t221 * t234 - t222 * t236) * MDP(16) + (t202 * t242 - t221 * t236) * MDP(15) + t327 * MDP(17) + t212 * MDP(18) + 0.2e1 * (-t355 * t293 * qJD(4) + t291 * t359) * MDP(9) + (t291 * t294 + 0.2e1 * t300 * t343) * MDP(8) + t253 * MDP(10) + t254 * MDP(11) + t291 * MDP(5);
t283 = t304 * qJDD(2);
t193 = qJDD(4) * pkin(4) - t216 * qJD(4) - t342 * t300 + t283;
t311 = t216 * t350 + g(1) * t286 + t217 * t234 + (-t216 * t292 - t193) * t299 + t382 * t287;
t310 = -t381 + t382;
t194 = t215 * qJD(4) + t300 * qJDD(2) + t342 * t304;
t309 = -g(1) * t287 + g(2) * t367 - g(3) * t366 + t328 * qJD(5) + t303 * t193 - t299 * t194 - t217 * t236;
t250 = t304 * t344;
t249 = t300 * t344;
t228 = t237 - t376;
t227 = t233 + t348;
t218 = t223 * t352;
t209 = -t232 * t300 + t304 * t340;
t208 = t232 * t304 + t300 * t340;
t1 = [qJDD(1) * MDP(1) + (t218 + t314 * t300 + (-t319 + t322) * t304) * MDP(13) + (t331 + (t297 ^ 2 + t298 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t232 * t293 - t356 * t291 + t310) * MDP(7) + (t227 * t236 + t228 * t202 - (t326 * qJD(5) + t208 * t303 + t209 * t299) * t292 - t325 * t290 + t336) * MDP(21) + (t227 * t234 + t228 * t203 + (-t325 * qJD(5) - t208 * t299 + t209 * t303) * t292 + t326 * t290 + t313) * MDP(20) + (t319 * t300 + t314 * t304 + t335) * MDP(14) + t312 + (t338 * t291 + t315 - t368) * MDP(6) + t331 * MDP(2) + (g(2) * t302 - g(3) * t306) * MDP(3); t254 * MDP(13) - t253 * MDP(14) + t212 * MDP(20) - t327 * MDP(21) + t357 * MDP(4); (t229 * t293 + t310) * MDP(7) + (t218 + t317 * t300 + (-t320 + t322) * t304) * MDP(13) + (t320 * t300 + t317 * t304 + t335) * MDP(14) + (t281 * t202 - (t324 * qJD(5) - t249 * t303 - t250 * t299) * t292 - t323 * t290 + t333 * t236 - t241 * t370 + t336) * MDP(21) + (t281 * t203 + (-t323 * qJD(5) + t249 * t299 - t250 * t303) * t292 + t324 * t290 + t333 * t234 + t242 * t370 + t313) * MDP(20) + t312 + (t315 + t369) * MDP(6); t300 * t291 * MDP(10) + MDP(11) * t364 + qJDD(4) * MDP(12) + (-g(1) * t304 + t316 * t300 + t283) * MDP(13) + (-t357 * t300 + t316 * t304) * MDP(14) + (-(-t215 * t299 - t371) * t292 + (-t234 * t363 + t303 * t290 - t292 * t350) * pkin(4) + t309) * MDP(20) + ((-qJD(5) * t214 + t215 * t292 - t194) * t303 + (-qJD(5) * t303 * t292 - t236 * t363 - t299 * t290) * pkin(4) + t311) * MDP(21) + t318 + (-MDP(8) * t359 + t355 * MDP(9)) * t293 ^ 2; (-t328 * t292 + t309) * MDP(20) + ((-t194 + (-qJD(5) + t292) * t214) * t303 + t311) * MDP(21) + t318;];
tau = t1;
