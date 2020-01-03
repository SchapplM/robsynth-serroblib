% Calculate vector of inverse dynamics joint torques for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:17
% EndTime: 2019-12-31 19:13:22
% DurationCPUTime: 3.09s
% Computational Cost: add. (2176->328), mult. (4227->430), div. (0->0), fcn. (2890->10), ass. (0->153)
t335 = sin(qJ(4));
t339 = cos(qJ(4));
t340 = cos(qJ(3));
t393 = qJD(1) * t340;
t336 = sin(qJ(3));
t394 = qJD(1) * t336;
t275 = t335 * t394 - t339 * t393;
t327 = qJD(3) + qJD(4);
t334 = sin(qJ(5));
t338 = cos(qJ(5));
t261 = -t275 * t334 - t338 * t327;
t281 = t335 * t340 + t336 * t339;
t276 = t281 * qJD(1);
t424 = qJD(5) + t276;
t432 = t261 * t424;
t263 = -t275 * t338 + t327 * t334;
t431 = t263 * t424;
t341 = cos(qJ(1));
t325 = g(2) * t341;
t337 = sin(qJ(1));
t414 = g(1) * t337;
t421 = -t325 + t414;
t326 = qJDD(3) + qJDD(4);
t342 = -pkin(1) - pkin(6);
t296 = t342 * qJD(1) + qJD(2);
t268 = -pkin(7) * t393 + t340 * t296;
t265 = qJD(3) * pkin(3) + t268;
t267 = -pkin(7) * t394 + t296 * t336;
t403 = t267 * t339;
t246 = t265 * t335 + t403;
t294 = qJDD(1) * t342 + qJDD(2);
t282 = t340 * t294;
t392 = qJD(3) * t336;
t383 = qJDD(1) * t340;
t385 = qJD(1) * qJD(3);
t426 = t336 * t385 - t383;
t251 = qJDD(3) * pkin(3) + t426 * pkin(7) - t296 * t392 + t282;
t391 = qJD(3) * t340;
t377 = t340 * t385;
t384 = qJDD(1) * t336;
t427 = t377 + t384;
t253 = -t427 * pkin(7) + t294 * t336 + t296 * t391;
t418 = t246 * qJD(4) - t339 * t251 + t253 * t335;
t224 = -pkin(4) * t326 + t418;
t333 = qJ(3) + qJ(4);
t323 = cos(t333);
t381 = t323 * t414;
t430 = t224 + t381;
t322 = sin(t333);
t429 = g(3) * t322 + t323 * t325;
t428 = qJDD(2) - t421;
t328 = qJDD(1) * qJ(2);
t329 = qJD(1) * qJD(2);
t362 = g(1) * t341 + g(2) * t337;
t350 = -t362 + 0.2e1 * t329;
t425 = 0.2e1 * t328 + t350;
t254 = -pkin(4) * t275 + pkin(8) * t276;
t423 = t424 * (pkin(8) * qJD(5) + t254);
t422 = t276 * t327;
t412 = pkin(1) * qJDD(1);
t420 = t412 - t428;
t390 = qJD(4) * t335;
t419 = (qJD(4) * t265 + t253) * t339 + t251 * t335 - t267 * t390;
t413 = pkin(7) - t342;
t286 = t413 * t336;
t287 = t413 * t340;
t259 = -t286 * t335 + t287 * t339;
t278 = t413 * t392;
t279 = qJD(3) * t287;
t234 = -qJD(4) * t259 + t278 * t335 - t279 * t339;
t379 = t336 * t390;
t356 = -qJD(1) * t379 - t426 * t335;
t366 = t327 * t340;
t241 = (qJD(1) * t366 + t384) * t339 + t356;
t239 = qJDD(5) + t241;
t404 = t267 * t335;
t245 = t265 * t339 - t404;
t242 = -pkin(4) * t327 - t245;
t280 = t335 * t336 - t339 * t340;
t310 = t336 * pkin(3) + qJ(2);
t255 = pkin(4) * t281 + pkin(8) * t280 + t310;
t389 = qJD(4) * t339;
t257 = -t335 * t391 - t336 * t389 - t339 * t392 - t340 * t390;
t260 = -t286 * t339 - t287 * t335;
t223 = pkin(8) * t326 + t419;
t288 = pkin(3) * t394 + qJD(1) * qJ(2);
t248 = pkin(4) * t276 + pkin(8) * t275 + t288;
t370 = qJD(5) * t248 + t223;
t417 = -t224 * t280 - t260 * t239 + t242 * t257 - (qJD(5) * t255 + t234) * t424 - t281 * t370;
t312 = g(3) * t323;
t344 = qJD(1) ^ 2;
t411 = qJ(2) * t344;
t240 = -t335 * t384 + t339 * t383 - t422;
t387 = qJD(5) * t338;
t380 = t338 * t240 + t334 * t326 + t327 * t387;
t388 = qJD(5) * t334;
t229 = t275 * t388 + t380;
t410 = t229 * t334;
t409 = t239 * t334;
t408 = t239 * t338;
t407 = t242 * t276;
t406 = t242 * t280;
t405 = t255 * t239;
t402 = t280 * t338;
t401 = t334 * t337;
t400 = t334 * t341;
t399 = t337 * t338;
t398 = t338 * t341;
t397 = t257 * t327 - t280 * t326;
t332 = t340 ^ 2;
t396 = t336 ^ 2 - t332;
t343 = qJD(3) ^ 2;
t395 = -t343 - t344;
t297 = pkin(3) * t391 + qJD(2);
t382 = qJDD(3) * t336;
t372 = t338 * t424;
t266 = t427 * pkin(3) + t328 + t329;
t226 = pkin(4) * t241 - pkin(8) * t240 + t266;
t243 = pkin(8) * t327 + t246;
t369 = qJD(5) * t243 - t226;
t315 = pkin(3) * t335 + pkin(8);
t367 = pkin(3) * t393 + qJD(5) * t315 + t254;
t249 = t268 * t335 + t403;
t364 = pkin(3) * t390 - t249;
t250 = t268 * t339 - t404;
t363 = -pkin(3) * t389 + t250;
t361 = -t239 * t315 + t407;
t258 = -t335 * t392 + t339 * t366 - t379;
t360 = -t258 * t327 - t281 * t326;
t228 = t243 * t338 + t248 * t334;
t358 = -t228 * t275 + t242 * t387 + t430 * t334;
t227 = -t243 * t334 + t248 * t338;
t357 = t227 * t275 + t242 * t388 + t429 * t338;
t354 = t257 * t338 + t280 * t388;
t353 = 0.2e1 * qJ(2) * t385 + qJDD(3) * t342;
t352 = -t411 - t421;
t351 = -pkin(8) * t239 + t245 * t424 + t407;
t348 = -t342 * t343 + t425;
t307 = t338 * t326;
t230 = qJD(5) * t263 + t240 * t334 - t307;
t347 = ((t229 - t432) * t338 + (-t230 - t431) * t334) * MDP(22) + (t263 * t372 + t410) * MDP(21) + (-t334 * t424 ^ 2 - t261 * t275 + t408) * MDP(24) + (t263 * t275 + t372 * t424 + t409) * MDP(23) + (t240 + t422) * MDP(16) + (-t275 * t327 + (-t327 * t393 - t384) * t339 - t356) * MDP(17) + (t275 ^ 2 - t276 ^ 2) * MDP(15) + t326 * MDP(18) + (-MDP(14) * t276 + MDP(25) * t424) * t275;
t346 = t276 * t288 + t421 * t322 + t312 - t419;
t345 = t275 * t288 - t381 - t418 + t429;
t321 = qJDD(3) * t340;
t316 = -pkin(3) * t339 - pkin(4);
t272 = t322 * t398 - t401;
t271 = t322 * t400 + t399;
t270 = t322 * t399 + t400;
t269 = -t322 * t401 + t398;
t235 = qJD(4) * t260 - t278 * t339 - t279 * t335;
t233 = pkin(4) * t258 - pkin(8) * t257 + t297;
t225 = t338 * t226;
t1 = [(-t336 * t343 + t321) * MDP(9) + (t420 * pkin(1) + (t350 + t328) * qJ(2)) * MDP(6) + t421 * MDP(2) + (-t240 * t281 + t241 * t280 - t257 * t276 + t258 * t275) * MDP(15) + (t336 * t348 + t340 * t353) * MDP(12) + (-t336 * t353 + t340 * t348) * MDP(13) + (-t340 * t343 - t382) * MDP(10) + t397 * MDP(16) + ((-t261 * t338 - t263 * t334) * t257 + (t410 + t230 * t338 + (-t261 * t334 + t263 * t338) * qJD(5)) * t280) * MDP(22) + t362 * MDP(3) + (-t234 * t327 + t240 * t310 + t257 * t288 - t260 * t326 - t266 * t280 - t275 * t297 - t323 * t362) * MDP(20) + t360 * MDP(17) + (-t235 * t327 + t241 * t310 + t258 * t288 - t259 * t326 + t266 * t281 + t276 * t297 - t322 * t362) * MDP(19) + qJDD(1) * MDP(1) + (-t240 * t280 - t257 * t275) * MDP(14) + (t239 * t281 + t258 * t424) * MDP(25) + (g(1) * t271 - g(2) * t269 - t228 * t258 + t259 * t229 + t235 * t263 + (-(-qJD(5) * t260 + t233) * t424 - t405 + t369 * t281 + qJD(5) * t406) * t334 + t417 * t338) * MDP(27) + (-g(1) * t272 - g(2) * t270 + t225 * t281 + t227 * t258 + t259 * t230 + t235 * t261 + (t233 * t424 + t405 + (-t243 * t281 - t260 * t424 - t406) * qJD(5)) * t338 + t417 * t334) * MDP(26) + (t280 * t409 - t230 * t281 - t258 * t261 + (-t257 * t334 + t280 * t387) * t424) * MDP(24) + (t229 * t281 - t239 * t402 + t258 * t263 + t354 * t424) * MDP(23) + (qJDD(1) * t332 - 0.2e1 * t336 * t377) * MDP(7) + 0.2e1 * (-t336 * t383 + t385 * t396) * MDP(8) + (-t229 * t402 + t263 * t354) * MDP(21) + t425 * MDP(5) + (-0.2e1 * t412 + t428) * MDP(4); qJDD(1) * MDP(4) - t344 * MDP(5) + (-t411 - t420) * MDP(6) + (t336 * t395 + t321) * MDP(12) + (t340 * t395 - t382) * MDP(13) + (-qJD(1) * t276 + t397) * MDP(19) + (qJD(1) * t275 + t360) * MDP(20) + (t230 * t280 - t257 * t261 - t281 * t409) * MDP(26) + (t229 * t280 - t257 * t263 - t281 * t408) * MDP(27) + ((-qJD(1) * t338 - t258 * t334 - t281 * t387) * MDP(26) + (qJD(1) * t334 - t258 * t338 + t281 * t388) * MDP(27)) * t424; -MDP(10) * t384 + (t249 * t327 + (-t276 * t393 + t326 * t339 - t327 * t390) * pkin(3) + t345) * MDP(19) + (t250 * t327 + (t275 * t393 - t326 * t335 - t327 * t389) * pkin(3) + t346) * MDP(20) + (g(3) * t340 + (-t294 - t352) * t336) * MDP(13) + (g(3) * t336 + t340 * t352 + t282) * MDP(12) + qJDD(3) * MDP(11) + t347 + (t316 * t230 - t430 * t338 + t361 * t334 + t364 * t261 + (t334 * t363 - t338 * t367) * t424 + t357) * MDP(26) + (t316 * t229 + t361 * t338 - t429 * t334 + t364 * t263 + (t334 * t367 + t338 * t363) * t424 + t358) * MDP(27) + MDP(9) * t383 + (MDP(7) * t336 * t340 - MDP(8) * t396) * t344; (t246 * t327 + t345) * MDP(19) + (t245 * t327 + t346) * MDP(20) + t347 + (-pkin(4) * t230 - t246 * t261 + t351 * t334 + (-t430 - t423) * t338 + t357) * MDP(26) + (-pkin(4) * t229 - t246 * t263 + t351 * t338 + (-t429 + t423) * t334 + t358) * MDP(27); t263 * t261 * MDP(21) + (-t261 ^ 2 + t263 ^ 2) * MDP(22) + (t380 + t432) * MDP(23) + (t307 + t431) * MDP(24) + t239 * MDP(25) + (-g(1) * t269 - g(2) * t271 + t228 * t424 - t242 * t263 + t225) * MDP(26) + (g(1) * t270 - g(2) * t272 + t227 * t424 + t242 * t261) * MDP(27) + ((-t223 + t312) * MDP(27) + (MDP(24) * t275 - MDP(26) * t243 - MDP(27) * t248) * qJD(5)) * t338 + (qJD(5) * t275 * MDP(23) + (-qJD(5) * t327 - t240) * MDP(24) + (-t370 + t312) * MDP(26) + t369 * MDP(27)) * t334;];
tau = t1;
