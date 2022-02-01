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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:49:06
% EndTime: 2022-01-20 09:49:10
% DurationCPUTime: 1.58s
% Computational Cost: add. (1347->212), mult. (2211->284), div. (0->0), fcn. (1414->14), ass. (0->132)
t394 = pkin(7) + pkin(8);
t307 = sin(qJ(5));
t308 = sin(qJ(4));
t311 = cos(qJ(5));
t312 = cos(qJ(4));
t246 = t307 * t312 + t308 * t311;
t301 = qJD(1) + qJD(3);
t240 = t246 * t301;
t309 = sin(qJ(3));
t313 = cos(qJ(3));
t306 = cos(pkin(9));
t285 = pkin(1) * t306 + pkin(2);
t269 = t285 * qJD(1);
t305 = sin(pkin(9));
t387 = pkin(1) * t305;
t391 = qJD(3) * t269 + qJDD(1) * t387;
t356 = qJD(1) * t387;
t392 = -qJD(3) * t356 + t285 * qJDD(1);
t341 = -t391 * t309 + t392 * t313;
t299 = qJDD(1) + qJDD(3);
t386 = pkin(3) * t299;
t215 = -t341 - t386;
t292 = qJ(1) + pkin(9) + qJ(3);
t284 = cos(t292);
t383 = g(2) * t284;
t393 = t215 + t383;
t366 = t311 * t312;
t370 = t307 * t308;
t245 = -t366 + t370;
t283 = sin(t292);
t390 = g(1) * t284 + g(2) * t283;
t343 = t285 * t313 - t309 * t387;
t363 = t309 * t285 + t313 * t387;
t234 = t269 * t309 + t313 * t356;
t347 = t394 * t301 + t234;
t219 = t312 * qJD(2) - t347 * t308;
t300 = qJD(4) + qJD(5);
t220 = qJD(2) * t308 + t347 * t312;
t389 = t392 * t309 + t313 * t391;
t385 = pkin(3) * t301;
t384 = pkin(4) * t312;
t279 = g(1) * t283;
t242 = pkin(7) + t363;
t382 = -pkin(8) - t242;
t381 = t220 * t311;
t233 = t269 * t313 - t309 * t356;
t380 = t233 * t300;
t379 = t234 * t301;
t237 = t363 * qJD(3);
t378 = t237 * t301;
t304 = qJ(4) + qJ(5);
t294 = sin(t304);
t377 = t283 * t294;
t295 = cos(t304);
t376 = t283 * t295;
t375 = t284 * t294;
t374 = t284 * t295;
t372 = t299 * t312;
t371 = t301 * t308;
t367 = t308 * t312;
t365 = qJDD(2) - g(3);
t227 = -t233 - t385;
t359 = qJD(4) * t308;
t364 = t227 * t359 + t312 * t279;
t302 = t308 ^ 2;
t362 = -t312 ^ 2 + t302;
t358 = qJD(4) * t312;
t357 = qJD(5) * t307;
t355 = pkin(4) * t359;
t354 = t301 * t370;
t352 = t301 * t366;
t351 = t227 * t358 + t393 * t308;
t289 = -pkin(3) - t384;
t350 = qJD(4) * t394;
t349 = t301 * t358;
t214 = pkin(7) * t299 + t389;
t348 = pkin(8) * t299 + t214;
t345 = qJD(4) * t382;
t241 = -pkin(3) - t343;
t340 = -t234 + t355;
t310 = sin(qJ(1));
t314 = cos(qJ(1));
t339 = g(1) * t310 - g(2) * t314;
t338 = t245 * t299;
t218 = qJD(4) * pkin(4) + t219;
t336 = -t218 * t307 - t381;
t225 = t300 * t245;
t298 = qJDD(4) + qJDD(5);
t335 = -t225 * t300 + t246 * t298;
t229 = t382 * t308;
t296 = t312 * pkin(8);
t230 = t242 * t312 + t296;
t334 = t229 * t311 - t230 * t307;
t333 = t229 * t307 + t230 * t311;
t274 = t394 * t308;
t275 = pkin(7) * t312 + t296;
t332 = -t274 * t311 - t275 * t307;
t331 = -t274 * t307 + t275 * t311;
t203 = (t301 * t359 - t372) * pkin(4) + t215;
t221 = t289 * t301 - t233;
t329 = -g(1) * t377 + g(2) * t375 + t203 * t246 - t221 * t225;
t226 = t300 * t246;
t328 = g(1) * t376 - g(2) * t374 + t203 * t245 + t221 * t226;
t315 = qJD(4) ^ 2;
t327 = pkin(7) * t315 - t379 - t386;
t326 = t279 + t341 - t383;
t325 = t241 * t299 + t242 * t315 + t378;
t206 = qJD(5) * t352 + t246 * t299 - t300 * t354 + t311 * t349;
t238 = -t352 + t354;
t324 = t240 * t238 * MDP(15) + (t238 * t300 + t206) * MDP(17) - t338 * MDP(18) + (-t238 ^ 2 + t240 ^ 2) * MDP(16) + t298 * MDP(19);
t323 = -pkin(7) * qJDD(4) + (t233 - t385) * qJD(4);
t322 = -t227 * t301 - t214 + t390;
t236 = t343 * qJD(3);
t321 = -qJDD(4) * t242 + (t241 * t301 - t236) * qJD(4);
t207 = t226 * t301 + t338;
t216 = -t226 * t300 - t245 * t298;
t259 = qJDD(4) * t308 + t312 * t315;
t260 = qJDD(4) * t312 - t308 * t315;
t320 = (-t206 * t245 - t207 * t246 + t225 * t238 - t226 * t240) * MDP(16) + (t206 * t246 - t225 * t240) * MDP(15) + t335 * MDP(17) + t216 * MDP(18) + 0.2e1 * (-t362 * t301 * qJD(4) + t299 * t367) * MDP(9) + (t299 * t302 + 0.2e1 * t308 * t349) * MDP(8) + t259 * MDP(10) + t260 * MDP(11) + t299 * MDP(5);
t291 = t312 * qJDD(2);
t197 = qJDD(4) * pkin(4) - qJD(4) * t220 - t348 * t308 + t291;
t319 = t221 * t238 + t220 * t357 + g(2) * t376 + g(1) * t374 + g(3) * t294 + (-t220 * t300 - t197) * t307;
t318 = -t389 + t390;
t198 = qJD(4) * t219 + t308 * qJDD(2) + t348 * t312;
t317 = g(1) * t375 + g(2) * t377 - g(3) * t295 + t336 * qJD(5) + t311 * t197 - t307 * t198 - t221 * t240;
t256 = t312 * t350;
t255 = t308 * t350;
t232 = t241 - t384;
t231 = t237 + t355;
t213 = -t236 * t308 + t312 * t345;
t212 = t236 * t312 + t308 * t345;
t1 = [(t343 * t299 + t326 - t378) * MDP(6) + t339 * MDP(2) + (g(1) * t314 + g(2) * t310) * MDP(3) + (-t236 * t301 - t363 * t299 + t318) * MDP(7) + t320 + (t321 * t312 + (t325 - t279) * t308 + t351) * MDP(14) + (t231 * t238 + t232 * t207 + (-t333 * qJD(5) - t212 * t307 + t213 * t311) * t300 + t334 * t298 + t328) * MDP(20) + (t231 * t240 + t232 * t206 - (t334 * qJD(5) + t212 * t311 + t213 * t307) * t300 - t333 * t298 + t329) * MDP(21) + (t321 * t308 + (-t325 - t393) * t312 + t364) * MDP(13) + (t339 + (t305 ^ 2 + t306 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + qJDD(1) * MDP(1); t260 * MDP(13) - t259 * MDP(14) + t216 * MDP(20) - t335 * MDP(21) + t365 * MDP(4); (t326 + t379) * MDP(6) + (t323 * t308 + (-t327 - t393) * t312 + t364) * MDP(13) + (t233 * t301 + t318) * MDP(7) + (t323 * t312 + (t327 - t279) * t308 + t351) * MDP(14) + t320 + (t289 * t207 + (-qJD(5) * t331 + t255 * t307 - t256 * t311) * t300 + t332 * t298 + t340 * t238 + t246 * t380 + t328) * MDP(20) + (t289 * t206 - (qJD(5) * t332 - t255 * t311 - t256 * t307) * t300 - t331 * t298 + t340 * t240 - t245 * t380 + t329) * MDP(21); t308 * t299 * MDP(10) + MDP(11) * t372 + qJDD(4) * MDP(12) + (-g(3) * t312 + t308 * t322 + t291) * MDP(13) + (-t365 * t308 + t322 * t312) * MDP(14) + (-(-t219 * t307 - t381) * t300 + (-t238 * t371 + t311 * t298 - t300 * t357) * pkin(4) + t317) * MDP(20) + ((-qJD(5) * t218 + t219 * t300 - t198) * t311 + (-qJD(5) * t311 * t300 - t240 * t371 - t307 * t298) * pkin(4) + t319) * MDP(21) + t324 + (-MDP(8) * t367 + t362 * MDP(9)) * t301 ^ 2; (-t336 * t300 + t317) * MDP(20) + ((-t198 + (-qJD(5) + t300) * t218) * t311 + t319) * MDP(21) + t324;];
tau = t1;
