% Calculate vector of inverse dynamics joint torques for
% S5RPRRP9
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:18
% EndTime: 2019-12-31 18:49:24
% DurationCPUTime: 3.64s
% Computational Cost: add. (3065->322), mult. (7253->385), div. (0->0), fcn. (5517->12), ass. (0->152)
t344 = sin(pkin(8));
t347 = sin(qJ(3));
t397 = qJD(1) * qJD(3);
t388 = t347 * t397;
t345 = cos(pkin(8));
t349 = cos(qJ(3));
t395 = qJDD(1) * t349;
t443 = qJDD(1) * t347 + t349 * t397;
t393 = t344 * t395 + t443 * t345;
t272 = -t344 * t388 + t393;
t411 = t345 * t349;
t413 = t344 * t347;
t296 = -t411 + t413;
t289 = t296 * qJD(1);
t297 = t344 * t349 + t345 * t347;
t290 = t297 * qJD(1);
t346 = sin(qJ(4));
t394 = -t443 * t344 - t345 * t388;
t366 = t345 * t395 + t394;
t429 = cos(qJ(4));
t389 = qJD(4) * t429;
t400 = qJD(4) * t346;
t236 = -t429 * t272 + t289 * t389 + t290 * t400 - t346 * t366;
t264 = t429 * t289 + t290 * t346;
t343 = qJD(3) + qJD(4);
t420 = t264 * t343;
t228 = -t236 + t420;
t368 = -t346 * t289 + t290 * t429;
t237 = qJD(4) * t368 + t346 * t272 - t429 * t366;
t338 = qJDD(3) + qJDD(4);
t421 = t264 ^ 2;
t422 = t368 * t343;
t441 = t368 ^ 2;
t448 = t264 * t368;
t449 = t228 * MDP(17) + t338 * MDP(19) + MDP(15) * t448 + (-t237 + t422) * MDP(18) + (-t421 + t441) * MDP(16);
t392 = pkin(2) * t345 + pkin(1);
t302 = -qJD(1) * t392 + qJD(2);
t275 = pkin(3) * t289 + t302;
t241 = pkin(4) * t264 - qJ(5) * t368 + t275;
t446 = t241 * t264;
t445 = t264 * t275;
t342 = pkin(8) + qJ(3);
t335 = cos(t342);
t336 = qJ(4) + t342;
t324 = sin(t336);
t325 = cos(t336);
t405 = t325 * pkin(4) + t324 * qJ(5);
t444 = pkin(3) * t335 + t405;
t425 = qJDD(1) * pkin(1);
t348 = sin(qJ(1));
t350 = cos(qJ(1));
t437 = g(1) * t348 - g(2) * t350;
t371 = -qJDD(2) + t425 + t437;
t248 = pkin(4) * t368 + qJ(5) * t264;
t426 = pkin(6) + qJ(2);
t310 = t426 * t344;
t298 = qJD(1) * t310;
t311 = t426 * t345;
t299 = qJD(1) * t311;
t439 = -t349 * t298 - t299 * t347;
t407 = -t347 * t310 + t349 * t311;
t330 = t338 * qJ(5);
t331 = t343 * qJD(5);
t438 = t330 + t331;
t436 = t345 * MDP(4) - t344 * MDP(5);
t435 = qJ(2) * qJDD(1);
t333 = t338 * pkin(4);
t434 = t333 - qJDD(5);
t372 = t298 * t347 - t299 * t349;
t398 = qJD(1) * qJD(2);
t431 = qJDD(1) * t426 + t398;
t279 = t431 * t344;
t280 = t431 * t345;
t383 = -t349 * t279 - t347 * t280;
t233 = qJDD(3) * pkin(3) - pkin(7) * t272 + qJD(3) * t372 + t383;
t373 = -t347 * t279 + t349 * t280;
t235 = t366 * pkin(7) + t439 * qJD(3) + t373;
t258 = -pkin(7) * t290 + t439;
t257 = qJD(3) * pkin(3) + t258;
t259 = -pkin(7) * t289 - t372;
t380 = -t429 * t233 + t346 * t235 + t257 * t400 + t259 * t389;
t417 = t324 * t350;
t418 = t324 * t348;
t362 = g(1) * t417 + g(2) * t418 - g(3) * t325 - t380;
t355 = t241 * t368 - t362 - t434;
t433 = -t275 * t368 + t362;
t402 = qJD(2) * t344;
t294 = t349 * t310;
t408 = qJD(2) * t411 - qJD(3) * t294;
t409 = t347 * t311;
t428 = pkin(7) * t297;
t251 = -t347 * t402 + (-t409 - t428) * qJD(3) + t408;
t291 = t296 * qJD(3);
t353 = -t297 * qJD(2) - t407 * qJD(3);
t252 = pkin(7) * t291 + t353;
t382 = -t294 - t409;
t261 = t382 - t428;
t262 = -pkin(7) * t296 + t407;
t369 = t261 * t429 - t346 * t262;
t225 = qJD(4) * t369 + t251 * t429 + t346 * t252;
t247 = t346 * t261 + t262 * t429;
t432 = t225 * t343 + t247 * t338 + t324 * t437;
t391 = t429 * t259;
t240 = t346 * t257 + t391;
t424 = t240 * t343;
t416 = t325 * t348;
t415 = t325 * t350;
t410 = t346 * t259;
t243 = t258 * t429 - t410;
t406 = pkin(3) * t389 + qJD(5) - t243;
t404 = t344 ^ 2 + t345 ^ 2;
t401 = qJD(3) * t290;
t239 = t257 * t429 - t410;
t399 = qJD(5) - t239;
t390 = qJD(1) * t413;
t385 = t404 * qJD(1) ^ 2;
t381 = t346 * t233 + t429 * t235 + t257 * t389 - t259 * t400;
t379 = 0.2e1 * t404;
t242 = t346 * t258 + t391;
t378 = pkin(3) * t400 - t242;
t334 = sin(t342);
t377 = -pkin(3) * t334 - pkin(4) * t324;
t376 = g(1) * t350 + g(2) * t348;
t277 = pkin(3) * t296 - t392;
t370 = t392 + t444;
t274 = -t346 * t296 + t297 * t429;
t367 = t297 * qJD(3);
t365 = -g(1) * t415 - g(2) * t416 - g(3) * t324 + t381;
t364 = pkin(3) * t367;
t226 = qJD(4) * t247 + t346 * t251 - t252 * t429;
t360 = g(1) * t416 - g(2) * t415 - t226 * t343 + t338 * t369;
t359 = t239 * t343 - t365;
t357 = t379 * t398 - t376;
t260 = qJDD(2) - t394 * pkin(3) + (-pkin(1) + (-t349 * pkin(3) - pkin(2)) * t345) * qJDD(1);
t224 = t237 * pkin(4) + t236 * qJ(5) - qJD(5) * t368 + t260;
t339 = -pkin(7) - t426;
t329 = -pkin(3) * t429 - pkin(4);
t326 = pkin(3) * t346 + qJ(5);
t304 = qJ(5) * t415;
t303 = qJ(5) * t416;
t301 = -qJDD(1) * t392 + qJDD(2);
t273 = t296 * t429 + t297 * t346;
t250 = qJD(4) * t274 - t346 * t291 + t367 * t429;
t249 = t291 * t429 + t296 * t389 + t297 * t400 + t346 * t367;
t245 = pkin(4) * t273 - qJ(5) * t274 + t277;
t244 = pkin(3) * t290 + t248;
t238 = t343 * qJ(5) + t240;
t234 = -t343 * pkin(4) + t399;
t227 = t250 * pkin(4) + t249 * qJ(5) - t274 * qJD(5) + t364;
t223 = t380 - t434;
t222 = t381 + t438;
t1 = [qJDD(1) * MDP(1) + t437 * MDP(2) + t376 * MDP(3) + (t379 * t435 + t357) * MDP(6) + (pkin(1) * t371 + (t404 * t435 + t357) * qJ(2)) * MDP(7) + (t272 * t297 - t290 * t291) * MDP(8) + (-t272 * t296 + t291 * t289 - t290 * t367 + t297 * t366) * MDP(9) + (-qJD(3) * t291 + qJDD(3) * t297) * MDP(10) + (-qJD(3) ^ 2 * t297 - t296 * qJDD(3)) * MDP(11) + (t392 * t366 + t301 * t296 + t382 * qJDD(3) + t437 * t335 + (t297 * t302 + t353) * qJD(3)) * MDP(13) + (-t392 * t272 + t301 * t297 - t302 * t291 - ((-qJD(3) * t311 - t402) * t347 + t408) * qJD(3) - t407 * qJDD(3) - t437 * t334) * MDP(14) + (-t236 * t274 - t249 * t368) * MDP(15) + (t236 * t273 - t237 * t274 + t249 * t264 - t250 * t368) * MDP(16) + (-t249 * t343 + t274 * t338) * MDP(17) + (-t250 * t343 - t273 * t338) * MDP(18) + (t277 * t237 + t275 * t250 + t260 * t273 + t264 * t364 + t360) * MDP(20) + (-t277 * t236 - t275 * t249 + t260 * t274 + t364 * t368 - t432) * MDP(21) + (t224 * t273 + t227 * t264 + t237 * t245 + t241 * t250 + t360) * MDP(22) + (-t222 * t273 + t223 * t274 - t225 * t264 + t226 * t368 - t234 * t249 + t236 * t369 - t237 * t247 - t238 * t250 - t376) * MDP(23) + (-t224 * t274 - t227 * t368 + t236 * t245 + t241 * t249 + t432) * MDP(24) + (t222 * t247 - t223 * t369 + t224 * t245 + t238 * t225 + t234 * t226 + t241 * t227 + (g(1) * t339 - g(2) * t370) * t350 + (g(1) * t370 + g(2) * t339) * t348) * MDP(25) + t436 * (t371 + t425); -MDP(6) * t385 + (-qJ(2) * t385 - t371) * MDP(7) + (-t366 + t401) * MDP(13) + ((-t289 - t390) * qJD(3) + t393) * MDP(14) + (-t421 - t441) * MDP(23) + (-t234 * t368 + t238 * t264 + t224 - t437) * MDP(25) + (-MDP(21) + MDP(24)) * (t236 + t420) + (MDP(20) + MDP(22)) * (t237 + t422) - t436 * qJDD(1); t290 * t289 * MDP(8) + (-t289 ^ 2 + t290 ^ 2) * MDP(9) + ((t289 - t390) * qJD(3) + t393) * MDP(10) + (t366 + t401) * MDP(11) + qJDD(3) * MDP(12) + (-g(3) * t335 - t302 * t290 + t334 * t376 + t383) * MDP(13) + (g(3) * t334 + t302 * t289 + t335 * t376 - t373) * MDP(14) + (t242 * t343 + (-t264 * t290 + t338 * t429 - t343 * t400) * pkin(3) + t433) * MDP(20) + (t243 * t343 + t445 + (-t290 * t368 - t338 * t346 - t343 * t389) * pkin(3) - t365) * MDP(21) + (-t244 * t264 - t329 * t338 - t343 * t378 - t355) * MDP(22) + (-t236 * t329 - t237 * t326 + (t238 + t378) * t368 + (t234 - t406) * t264) * MDP(23) + (t244 * t368 + t326 * t338 + t343 * t406 + t365 + t438 - t446) * MDP(24) + (t222 * t326 + t223 * t329 - t241 * t244 - g(1) * (t350 * t377 + t304) - g(2) * (t348 * t377 + t303) - g(3) * t444 + t406 * t238 + t378 * t234) * MDP(25) + t449; (t424 + t433) * MDP(20) + (t359 + t445) * MDP(21) + (-t248 * t264 + t333 - t355 + t424) * MDP(22) + (pkin(4) * t236 - qJ(5) * t237 + (t238 - t240) * t368 + (t234 - t399) * t264) * MDP(23) + (t248 * t368 + 0.2e1 * t330 + 0.2e1 * t331 - t359 - t446) * MDP(24) + (t222 * qJ(5) - t223 * pkin(4) - t241 * t248 - t234 * t240 - g(1) * (-pkin(4) * t417 + t304) - g(2) * (-pkin(4) * t418 + t303) - g(3) * t405 + t399 * t238) * MDP(25) + t449; (-t338 + t448) * MDP(22) + t228 * MDP(23) + (-t343 ^ 2 - t441) * MDP(24) + (-t238 * t343 + t355) * MDP(25);];
tau = t1;
