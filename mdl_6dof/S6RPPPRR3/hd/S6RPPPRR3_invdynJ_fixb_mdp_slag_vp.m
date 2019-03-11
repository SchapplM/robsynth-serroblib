% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:34:06
% EndTime: 2019-03-09 01:34:11
% DurationCPUTime: 3.29s
% Computational Cost: add. (2174->358), mult. (4133->456), div. (0->0), fcn. (2974->12), ass. (0->169)
t339 = sin(pkin(10));
t341 = cos(pkin(10));
t345 = sin(qJ(5));
t347 = cos(qJ(5));
t288 = t339 * t347 + t341 * t345;
t283 = t288 * qJD(1);
t344 = sin(qJ(6));
t346 = cos(qJ(6));
t265 = -t346 * qJD(5) - t283 * t344;
t408 = qJD(1) * t347;
t300 = t341 * t408;
t409 = qJD(1) * t345;
t282 = t339 * t409 - t300;
t400 = qJD(6) - t282;
t445 = t400 * t265;
t267 = qJD(5) * t344 - t283 * t346;
t444 = t400 * t267;
t422 = t339 * t345;
t388 = qJD(5) * t422;
t296 = qJD(1) * t388;
t355 = -t288 * qJDD(1) + t296;
t255 = -qJD(5) * t300 + t355;
t443 = qJD(5) * qJD(6) + t255;
t348 = -pkin(1) - pkin(2);
t298 = t348 * qJD(1) + qJD(2);
t340 = sin(pkin(9));
t342 = cos(pkin(9));
t410 = qJD(1) * t342;
t279 = qJ(2) * t410 + t340 * t298;
t269 = -qJD(1) * qJ(4) + t279;
t320 = t341 * qJD(3);
t248 = t320 + (pkin(7) * qJD(1) - t269) * t339;
t259 = t339 * qJD(3) + t341 * t269;
t411 = qJD(1) * t341;
t249 = -pkin(7) * t411 + t259;
t230 = t248 * t345 + t249 * t347;
t398 = qJD(1) * qJD(2);
t306 = t342 * t398;
t297 = t348 * qJDD(1) + qJDD(2);
t394 = qJDD(1) * t342;
t416 = qJ(2) * t394 + t340 * t297;
t263 = t306 + t416;
t257 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t263;
t318 = t341 * qJDD(3);
t238 = t318 + (pkin(7) * qJDD(1) - t257) * t339;
t242 = t339 * qJDD(3) + t341 * t257;
t395 = qJDD(1) * t341;
t239 = -pkin(7) * t395 + t242;
t372 = -t238 * t347 + t239 * t345;
t222 = -qJDD(5) * pkin(5) + qJD(5) * t230 + t372;
t337 = pkin(10) + qJ(5);
t321 = sin(t337);
t322 = cos(t337);
t434 = sin(qJ(1));
t435 = cos(qJ(1));
t289 = t435 * t340 - t434 * t342;
t287 = -t434 * t340 - t435 * t342;
t433 = g(1) * t287;
t378 = g(2) * t289 + t433;
t354 = g(3) * t322 - t378 * t321;
t442 = t354 - (-pkin(5) * t283 + t400 * pkin(8)) * t400 - t222;
t385 = t400 * t346;
t285 = t288 * qJD(5);
t299 = t347 * t395;
t376 = qJDD(1) * t422 - t299;
t256 = qJD(1) * t285 + t376;
t253 = -qJDD(6) + t256;
t420 = t344 * t253;
t441 = t400 * t385 - t420;
t439 = t341 * MDP(10) - t339 * MDP(11);
t438 = qJDD(1) * pkin(3) + qJDD(4);
t229 = t248 * t347 - t249 * t345;
t227 = -qJD(5) * pkin(5) - t229;
t419 = t347 * t341;
t286 = -t419 + t422;
t302 = qJD(2) * t342 - qJD(4);
t294 = t342 * qJ(2) + t340 * t348;
t290 = -qJ(4) + t294;
t429 = pkin(7) - t290;
t274 = t429 * t339;
t275 = t429 * t341;
t368 = t274 * t347 + t275 * t345;
t231 = t368 * qJD(5) - t286 * t302;
t237 = t274 * t345 - t275 * t347;
t293 = -t340 * qJ(2) + t342 * t348;
t291 = pkin(3) - t293;
t281 = t341 * pkin(4) + t291;
t240 = -pkin(5) * t286 + pkin(8) * t288 + t281;
t284 = -qJD(5) * t419 + t388;
t412 = qJD(1) * t340;
t278 = -qJ(2) * t412 + t298 * t342;
t268 = qJD(1) * pkin(3) + qJD(4) - t278;
t261 = pkin(4) * t411 + t268;
t235 = -pkin(5) * t282 + pkin(8) * t283 + t261;
t371 = t238 * t345 + t239 * t347;
t384 = qJDD(5) * pkin(8) + qJD(5) * t229 + qJD(6) * t235 + t371;
t436 = -t222 * t288 + t227 * t284 + t237 * t253 - (qJD(6) * t240 + t231) * t400 + t384 * t286 - t433;
t432 = g(2) * t287;
t431 = g(3) * t321;
t428 = pkin(1) * qJDD(1);
t390 = t344 * qJDD(5) + t443 * t346;
t403 = qJD(6) * t344;
t233 = t283 * t403 + t390;
t427 = t233 * t344;
t426 = t265 * t283;
t425 = t267 * t283;
t424 = t287 * t322;
t423 = t289 * t322;
t349 = qJD(1) ^ 2;
t421 = t342 * t349;
t244 = t346 * t253;
t418 = t340 * t285 - t286 * t410;
t277 = t286 * t340;
t417 = -qJD(5) * t277 - t342 * t283;
t415 = t435 * pkin(1) + t434 * qJ(2);
t414 = g(1) * t434 - g(2) * t435;
t413 = t339 ^ 2 + t341 ^ 2;
t407 = qJD(2) * t340;
t228 = qJD(5) * pkin(8) + t230;
t406 = qJD(6) * t228;
t405 = qJD(6) * t283;
t404 = qJD(6) * t288;
t399 = qJ(2) * qJDD(1);
t396 = qJDD(1) * t340;
t393 = 0.2e1 * t398;
t392 = t288 * t420;
t391 = t288 * t244;
t389 = t435 * pkin(2) + t415;
t304 = t340 * t398;
t386 = -qJ(2) * t396 + t297 * t342;
t382 = qJDD(1) * t413;
t381 = qJDD(2) - t428;
t380 = -t434 * pkin(1) + t435 * qJ(2);
t379 = g(1) * t289 - t432;
t262 = -t304 + t386;
t377 = -qJD(6) * t342 - t418;
t375 = -t406 + t432;
t374 = -t233 * t286 - t267 * t285;
t326 = t346 * qJDD(5);
t234 = qJD(6) * t267 + t255 * t344 - t326;
t373 = t234 * t286 + t265 * t285;
t241 = -t257 * t339 + t318;
t370 = -t241 * t339 + t242 * t341;
t369 = (-t269 * t339 + t320) * t339 - t259 * t341;
t367 = t278 * t340 - t279 * t342;
t366 = -qJD(6) * t277 + t412;
t364 = qJD(5) * t284 - qJDD(5) * t288;
t363 = qJD(5) * t285 + qJDD(5) * t286;
t362 = -t244 + (t282 * t344 - t403) * t400;
t361 = -g(1) * t423 - t240 * t253;
t360 = -t284 * t344 + t346 * t404;
t359 = t284 * t346 + t288 * t403;
t358 = g(1) * t435 + g(2) * t434;
t260 = -t262 + t438;
t357 = -t434 * pkin(2) + t380;
t251 = pkin(4) * t395 + t260;
t356 = -t379 - t386;
t353 = pkin(8) * t253 + (t227 + t229) * t400;
t352 = -g(2) * t423 - t384 - t431;
t276 = t288 * t340;
t246 = t289 * t344 - t346 * t424;
t245 = t289 * t346 + t344 * t424;
t243 = -pkin(5) * t285 - pkin(8) * t284 + t407;
t232 = t237 * qJD(5) + t288 * t302;
t226 = -pkin(5) * t256 - pkin(8) * t255 + t251;
t225 = t346 * t226;
t224 = t228 * t346 + t235 * t344;
t223 = -t228 * t344 + t235 * t346;
t1 = [(-qJD(5) * t232 + qJDD(5) * t368 - t251 * t286 - t256 * t281 - t261 * t285 - t282 * t407 - t379 * t322) * MDP(19) + (-g(2) * t246 - t223 * t285 - t225 * t286 + t232 * t265 - t368 * t234 + (t243 * t400 + (-t227 * t288 + t228 * t286 - t237 * t400) * qJD(6) + t361) * t346 + t436 * t344) * MDP(26) + (-g(2) * t245 + t224 * t285 + t232 * t267 - t368 * t233 + (-(-qJD(6) * t237 + t243) * t400 + (t226 - t406) * t286 + t227 * t404 - t361) * t344 + t436 * t346) * MDP(27) + (t359 * t400 + t374 + t391) * MDP(23) + (t360 * t400 + t373 - t392) * MDP(24) + (t253 * t286 - t285 * t400) * MDP(25) + t439 * (qJDD(1) * t291 + t260 + t304 - t379) + (-t358 + t393 + 0.2e1 * t399) * MDP(5) + ((-t265 * t346 - t267 * t344) * t284 + (t427 + t234 * t346 + (-t265 * t344 + t267 * t346) * qJD(6)) * t288) * MDP(22) + (-qJDD(2) + t414 + 0.2e1 * t428) * MDP(4) + qJDD(1) * MDP(1) + t363 * MDP(17) + t364 * MDP(16) + (-g(1) * t357 - g(2) * t389 - t367 * qJD(2) + t262 * t293 + t263 * t294) * MDP(9) + (-qJDD(1) * t293 + 0.2e1 * t304 + t356) * MDP(7) + t358 * MDP(3) + (-t233 * t288 * t346 + t359 * t267) * MDP(21) + t414 * MDP(2) + (-t381 * pkin(1) - g(1) * t380 - g(2) * t415 + (t393 + t399) * qJ(2)) * MDP(6) + (qJDD(1) * t294 + 0.2e1 * t306 + t378 + t416) * MDP(8) + (-t413 * t302 * qJD(1) - t290 * t382 - t370 - t378) * MDP(12) + (-qJD(5) * t231 - qJDD(5) * t237 - t251 * t288 + t255 * t281 + t261 * t284 - t283 * t407 + t379 * t321) * MDP(20) + (t260 * t291 + t268 * t407 - g(1) * (t289 * pkin(3) + t287 * qJ(4) + t357) - g(2) * (-pkin(3) * t287 + qJ(4) * t289 + t389) - t369 * t302 + t370 * t290) * MDP(13) + (-t255 * t288 - t283 * t284) * MDP(14) + (t255 * t286 - t256 * t288 + t282 * t284 - t283 * t285) * MDP(15); -qJDD(1) * MDP(4) - t349 * MDP(5) + (-qJ(2) * t349 + t381 - t414) * MDP(6) + (t396 - t421) * MDP(8) + (t367 * qJD(1) + t262 * t342 + t263 * t340 - t414) * MDP(9) + (-t340 * t382 + t413 * t421) * MDP(12) + (-t260 * t342 + t370 * t340 + (-t268 * t340 + t369 * t342) * qJD(1) - t414) * MDP(13) + (-t417 * qJD(5) - qJDD(5) * t276 + t256 * t342 + t282 * t412) * MDP(19) + (t418 * qJD(5) + qJDD(5) * t277 - t255 * t342 + t283 * t412) * MDP(20) + (-(t277 * t344 - t342 * t346) * t253 + t276 * t234 - (t377 * t344 + t366 * t346) * t400 + t417 * t265) * MDP(26) + ((-t277 * t346 - t342 * t344) * t253 + t276 * t233 - (-t366 * t344 + t377 * t346) * t400 + t417 * t267) * MDP(27) + (-MDP(7) - t439) * (t340 * t349 + t394); (qJDD(3) + g(3)) * MDP(9) + (t241 * t341 + t242 * t339 + g(3)) * MDP(13) - t363 * MDP(19) + t364 * MDP(20) + (t373 + t392) * MDP(26) + (-t374 + t391) * MDP(27) - (t360 * MDP(26) - t359 * MDP(27)) * t400; (-t369 * qJD(1) + t304 + t356 + t438) * MDP(13) + t299 * MDP(19) + t296 * MDP(20) + (t362 + t426) * MDP(26) + (t425 - t441) * MDP(27) - t413 * MDP(12) * t349 + ((-MDP(20) * t345 + MDP(10)) * t341 + (-MDP(19) * t345 - MDP(20) * t347 - MDP(11)) * t339) * qJDD(1) + ((-t339 * t408 - t341 * t409 - t283) * MDP(19) + (t282 - t300) * MDP(20)) * qJD(5); -t282 ^ 2 * MDP(15) + ((-t282 - t300) * qJD(5) + t355) * MDP(16) + t376 * MDP(17) + qJDD(5) * MDP(18) + (t354 - t372) * MDP(19) + (-t261 * t282 - t378 * t322 - t371 - t431) * MDP(20) + (t267 * t385 + t427) * MDP(21) + ((t233 - t445) * t346 + (-t234 - t444) * t344) * MDP(22) + (t425 + t441) * MDP(23) + (t362 - t426) * MDP(24) + (-pkin(5) * t234 - t230 * t265 + t353 * t344 + t442 * t346) * MDP(26) + (-pkin(5) * t233 - t230 * t267 - t442 * t344 + t353 * t346) * MDP(27) + (t282 * MDP(14) + MDP(15) * t283 + t261 * MDP(19) + MDP(25) * t400 + t223 * MDP(26) - t224 * MDP(27)) * t283; t267 * t265 * MDP(21) + (-t265 ^ 2 + t267 ^ 2) * MDP(22) + (t390 + t445) * MDP(23) + (t326 + t444) * MDP(24) - t253 * MDP(25) + (-g(1) * t245 + t224 * t400 - t227 * t267 + t225) * MDP(26) + (g(1) * t246 + t223 * t400 + t227 * t265) * MDP(27) + (MDP(24) * t405 + t375 * MDP(26) + t352 * MDP(27)) * t346 + (MDP(23) * t405 - t443 * MDP(24) + t352 * MDP(26) + (-t226 - t375) * MDP(27)) * t344;];
tau  = t1;
