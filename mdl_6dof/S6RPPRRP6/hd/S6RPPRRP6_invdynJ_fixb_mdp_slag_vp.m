% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:24
% EndTime: 2019-03-09 02:11:31
% DurationCPUTime: 4.15s
% Computational Cost: add. (2486->423), mult. (4248->509), div. (0->0), fcn. (2360->6), ass. (0->166)
t324 = sin(qJ(4));
t295 = qJD(1) * t324 + qJD(5);
t301 = qJ(2) * qJD(1) + qJD(3);
t291 = -pkin(7) * qJD(1) + t301;
t285 = t324 * t291;
t266 = qJD(4) * pkin(8) + t285;
t327 = cos(qJ(4));
t367 = pkin(4) * t327 + pkin(8) * t324;
t276 = qJD(4) * t367 + qJD(3);
t435 = pkin(8) * t327;
t366 = pkin(4) * t324 - t435;
t314 = qJDD(1) * qJ(3);
t320 = qJDD(1) * pkin(1);
t381 = qJDD(2) - t320;
t373 = -t314 + t381;
t243 = qJD(1) * t276 + qJDD(1) * t366 - t373;
t315 = qJD(1) * qJD(2);
t316 = qJ(2) * qJDD(1);
t372 = qJDD(3) + t315 + t316;
t282 = -pkin(7) * qJDD(1) + t372;
t395 = qJD(4) * t327;
t251 = qJDD(4) * pkin(8) + t282 * t324 + t291 * t395;
t322 = pkin(1) + qJ(3);
t286 = t366 + t322;
t260 = qJD(1) * t286 - qJD(2);
t323 = sin(qJ(5));
t326 = cos(qJ(5));
t392 = qJD(5) * t326;
t379 = -t243 * t323 - t251 * t326 - t260 * t392;
t393 = qJD(5) * t323;
t345 = -t266 * t393 - t379;
t389 = qJD(1) * qJD(4);
t375 = t327 * t389;
t386 = qJDD(1) * t324;
t275 = qJDD(5) + t375 + t386;
t429 = qJ(6) * t275;
t228 = qJD(6) * t295 + t345 + t429;
t370 = t243 * t326 - t251 * t323 - t260 * t393 - t266 * t392;
t437 = pkin(5) * t275;
t229 = qJDD(6) - t370 - t437;
t325 = sin(qJ(1));
t328 = cos(qJ(1));
t403 = g(1) * t325 - g(2) * t328;
t446 = -t228 * t323 + t229 * t326 - t403;
t238 = t260 * t326 - t266 * t323;
t390 = qJD(6) - t238;
t235 = -pkin(5) * t295 + t390;
t239 = t260 * t323 + t266 * t326;
t236 = qJ(6) * t295 + t239;
t355 = t235 * t323 + t236 * t326;
t397 = qJD(4) * t323;
t399 = qJD(1) * t327;
t280 = t326 * t399 + t397;
t420 = t280 * t323;
t391 = t326 * qJD(4);
t278 = t323 * t399 - t391;
t422 = t278 * t326;
t382 = MDP(23) - MDP(26);
t383 = MDP(22) + MDP(24);
t440 = t323 * t383 + t326 * t382;
t332 = t295 * t440 + (-t420 + t422) * MDP(25) - t355 * MDP(27);
t321 = -pkin(7) + qJ(2);
t445 = qJD(2) * t324 + t321 * t395;
t363 = g(1) * t328 + g(2) * t325;
t349 = t363 * t327;
t432 = g(3) * t324;
t444 = t432 - t349;
t385 = qJDD(1) * t327;
t394 = qJD(5) * t280;
t414 = t323 * t324;
t245 = -qJDD(4) * t326 + t323 * t385 - t389 * t414 + t394;
t443 = qJD(1) * t322;
t417 = t291 * t327;
t267 = -qJD(4) * pkin(4) - t417;
t237 = pkin(5) * t278 - qJ(6) * t280 + t267;
t436 = pkin(8) * t275;
t442 = t237 * t295 - t436;
t343 = t324 * t391 + t327 * t393;
t244 = qJD(1) * t343 - qJD(5) * t391 - qJDD(4) * t323 - t326 * t385;
t292 = -qJD(2) + t443;
t441 = qJD(4) * (qJD(2) + t292 + t443) + qJDD(4) * t321;
t438 = t280 ^ 2;
t309 = 0.2e1 * t315;
t431 = g(3) * t327;
t430 = pkin(8) * qJD(5);
t428 = t236 * t295;
t427 = t239 * t295;
t426 = t244 * t323;
t425 = t245 * t326;
t424 = t278 * t280;
t423 = t278 * t323;
t421 = t280 * t295;
t419 = t280 * t326;
t418 = t286 * t326;
t416 = t295 * t326;
t415 = t321 * t323;
t413 = t324 * t326;
t412 = t324 * t328;
t411 = t325 * t323;
t410 = t325 * t326;
t409 = t326 * t327;
t408 = t328 * t326;
t360 = pkin(5) * t323 - qJ(6) * t326;
t407 = -qJD(6) * t323 + t295 * t360 - t285;
t284 = t367 * qJD(1);
t406 = t284 * t323 + t291 * t409;
t405 = t286 * t323 + t321 * t413;
t404 = pkin(1) * t328 + qJ(2) * t325;
t319 = t327 ^ 2;
t402 = t324 ^ 2 - t319;
t329 = qJD(4) ^ 2;
t330 = qJD(1) ^ 2;
t401 = -t329 - t330;
t396 = qJD(4) * t324;
t388 = qJD(3) * qJD(1);
t387 = qJDD(1) * t322;
t380 = t295 * t430;
t377 = qJ(3) * t328 + t404;
t374 = qJDD(2) - t403;
t371 = t278 * t383;
t369 = t323 * t276 + t286 * t392 + t326 * t445;
t368 = -t320 + t374;
t268 = t324 * t411 - t408;
t270 = t323 * t412 + t410;
t365 = -g(1) * t268 + g(2) * t270;
t269 = t323 * t328 + t324 * t410;
t271 = t324 * t408 - t411;
t364 = g(1) * t269 - g(2) * t271;
t361 = pkin(5) * t326 + qJ(6) * t323;
t358 = t228 * t326 + t229 * t323;
t356 = t235 * t326 - t236 * t323;
t353 = -t314 + t368;
t352 = pkin(4) + t361;
t351 = -t321 + t360;
t250 = -qJDD(4) * pkin(4) - t282 * t327 + t291 * t396;
t348 = t275 * t323 + t295 * t392;
t347 = t275 * t326 - t295 * t393;
t346 = t309 + 0.2e1 * t316 - t363;
t344 = MDP(27) * t237 + t280 * t382;
t342 = t267 * t295 - t436;
t341 = qJD(1) * t292 - t282 + t363;
t340 = t323 * t382 - t326 * t383;
t338 = -t324 * t363 - t431;
t337 = g(1) * t270 + g(2) * t268 + t323 * t431 + t370;
t230 = pkin(5) * t245 + qJ(6) * t244 - qJD(6) * t280 + t250;
t336 = -t230 - t349 - t380;
t335 = t237 * t280 + qJDD(6) - t337;
t283 = -t373 + t388;
t334 = -t321 * t329 + t283 + t387 + t388 + t403;
t333 = -g(1) * t271 - g(2) * t269 - g(3) * t409 + t345;
t331 = (t419 + t423) * MDP(25) + t356 * MDP(27) + t340 * t295;
t308 = t328 * qJ(2);
t304 = qJDD(4) * t327;
t296 = g(3) * t413;
t258 = t351 * t327;
t255 = -t418 + (-pkin(5) + t415) * t324;
t254 = qJ(6) * t324 + t405;
t252 = pkin(5) * t280 + qJ(6) * t278;
t249 = -t284 * t326 + (-pkin(5) * qJD(1) + t291 * t323) * t327;
t246 = qJ(6) * t399 + t406;
t242 = t244 * t326;
t234 = -t351 * t396 + (qJD(5) * t361 - qJD(6) * t326 - qJD(2)) * t327;
t233 = t278 * t295 - t244;
t232 = -pkin(5) * t395 + (qJD(5) * t321 * t324 - t276) * t326 + (qJD(5) * t286 + t445) * t323;
t231 = qJ(6) * t395 + (-t321 * t393 + qJD(6)) * t324 + t369;
t1 = [(qJDD(3) + t346) * MDP(7) + t346 * MDP(5) + (-t231 * t278 + t232 * t280 - t244 * t255 - t245 * t254 - t356 * t396 + (-qJD(5) * t355 + t446) * t327) * MDP(25) + (-t324 * t441 + t334 * t327) * MDP(16) + (t334 * t324 + t327 * t441) * MDP(15) + t363 * MDP(3) + (-t324 * t329 + t304) * MDP(12) + (-qJDD(4) * t324 - t327 * t329) * MDP(13) + (-0.2e1 * t320 + t374) * MDP(4) + (qJDD(1) * t319 - 0.2e1 * t324 * t375) * MDP(10) + qJDD(1) * MDP(1) + (-t353 + t387 + 0.2e1 * t388) * MDP(8) + ((-t295 * t391 - t244) * t324 + (qJD(4) * t280 + t347) * t327) * MDP(19) + (t231 * t295 - t234 * t280 + t244 * t258 + t254 * t275 + (t237 * t391 + t228) * t324 + (qJD(4) * t236 - t230 * t326 + t237 * t393) * t327 - t365) * MDP(26) + (t275 * t324 + t295 * t395) * MDP(21) + (-t232 * t295 + t234 * t278 + t245 * t258 - t255 * t275 + (-t237 * t397 - t229) * t324 + (-qJD(4) * t235 + t230 * t323 + t237 * t392) * t327 + t364) * MDP(24) + ((t295 * t397 - t245) * t324 + (-qJD(4) * t278 - t348) * t327) * MDP(20) + 0.2e1 * (-t324 * t385 + t389 * t402) * MDP(11) + t403 * MDP(2) + (-t369 * t295 - t405 * t275 + ((t295 * t321 + t266) * t393 + (-t267 * t326 + t280 * t321) * qJD(4) + t379) * t324 + (-qJD(2) * t280 - qJD(4) * t239 + t244 * t321 + t250 * t326 - t267 * t393) * t327 + t365) * MDP(23) + (-t244 * t409 - t280 * t343) * MDP(17) + ((t276 * t326 - t286 * t393) * t295 + t275 * t418 + ((-qJD(2) * t323 - t321 * t392) * t295 - t275 * t415 + (-t267 * t323 + t278 * t321) * qJD(4) + t370) * t324 + (t267 * t392 - qJD(2) * t278 - t321 * t245 + t250 * t323 + (-t295 * t415 + t238) * qJD(4)) * t327 + t364) * MDP(22) + ((t420 + t422) * t396 + (t426 - t425 + (-t419 + t423) * qJD(5)) * t327) * MDP(18) + (t283 * t322 + t292 * qJD(3) + t372 * qJ(2) + t301 * qJD(2) - g(1) * (-t322 * t325 + t308) - g(2) * t377) * MDP(9) + (t228 * t254 + t236 * t231 + t230 * t258 + t237 * t234 + t229 * t255 + t235 * t232 - g(1) * (-pkin(5) * t269 - pkin(7) * t328 - qJ(6) * t268 + t308) - g(2) * (pkin(4) * t412 + pkin(5) * t271 + qJ(6) * t270 - t328 * t435 + t377) + (g(2) * pkin(7) + g(1) * t286) * t325) * MDP(27) + (-t381 * pkin(1) - g(1) * (-pkin(1) * t325 + t308) - g(2) * t404 + (t309 + t316) * qJ(2)) * MDP(6); t368 * MDP(6) + t353 * MDP(9) + (t245 * t323 - t242) * MDP(25) + t446 * MDP(27) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t330 + t371 * t399 + t340 * t275 + t332 * qJD(5) + (-MDP(15) * t324 - MDP(16) * t327 + MDP(4) - MDP(8)) * qJDD(1) + ((-qJD(3) - t301) * MDP(9) + (-0.2e1 * MDP(15) * qJD(4) + t344) * t327 + (0.2e1 * qJD(4) * MDP(16) + t332) * t324) * qJD(1); qJDD(1) * MDP(7) - t330 * MDP(8) + (-t363 + t372) * MDP(9) + t304 * MDP(15) - t363 * MDP(27) + t371 * t396 + (-t292 * MDP(9) + t331) * qJD(1) + (t401 * MDP(16) - t230 * MDP(27) - qJD(4) * t332 + t382 * t244 - t383 * t245) * t327 + (t401 * MDP(15) - qJDD(4) * MDP(16) + (-t425 - t426) * MDP(25) + t358 * MDP(27) + t344 * qJD(4) - t440 * t275 + t331 * qJD(5)) * t324; MDP(12) * t385 - MDP(13) * t386 + qJDD(4) * MDP(14) + (-t327 * t341 + t432) * MDP(15) + (t324 * t341 + t431) * MDP(16) + (t280 * t416 - t426) * MDP(17) + (-t242 - t278 * t416 + (-t245 - t421) * t323) * MDP(18) + ((-t280 * t327 + t295 * t413) * qJD(1) + t348) * MDP(19) + ((t278 * t327 - t295 * t414) * qJD(1) + t347) * MDP(20) - t295 * MDP(21) * t399 + (-t238 * t399 - t278 * t285 - pkin(4) * t245 + t296 + (t295 * t417 + t342) * t323 + (-t250 - t349 + (-t284 - t430) * t295) * t326) * MDP(22) + (pkin(4) * t244 + t406 * t295 + t239 * t399 - t280 * t285 + t342 * t326 + (t250 + t380 - t444) * t323) * MDP(23) + (t235 * t399 - t245 * t352 + t249 * t295 + t407 * t278 + t323 * t442 + t336 * t326 + t296) * MDP(24) + (t246 * t278 - t249 * t280 + (t228 + t295 * t235 + (-t245 + t394) * pkin(8)) * t326 + (t229 - t428 + (qJD(5) * t278 - t244) * pkin(8)) * t323 + t338) * MDP(25) + (-t236 * t399 - t244 * t352 - t246 * t295 - t407 * t280 - t442 * t326 + (t336 + t432) * t323) * MDP(26) + (-t235 * t249 - t236 * t246 + t407 * t237 + (qJD(5) * t356 + t338 + t358) * pkin(8) + (-t230 + t444) * t352) * MDP(27) + (MDP(10) * t324 * t327 - MDP(11) * t402) * t330; MDP(17) * t424 + (-t278 ^ 2 + t438) * MDP(18) + t233 * MDP(19) + (t421 - t245) * MDP(20) + t275 * MDP(21) + (-t267 * t280 + t337 + t427) * MDP(22) + (t238 * t295 + t267 * t278 - t333) * MDP(23) + (-t252 * t278 - t335 + t427 + 0.2e1 * t437) * MDP(24) + (pkin(5) * t244 - qJ(6) * t245 + (t236 - t239) * t280 + (t235 - t390) * t278) * MDP(25) + (0.2e1 * t429 - t237 * t278 + t252 * t280 + (0.2e1 * qJD(6) - t238) * t295 + t333) * MDP(26) + (t228 * qJ(6) - t229 * pkin(5) - t237 * t252 - t235 * t239 - g(1) * (-pkin(5) * t270 + qJ(6) * t271) - g(2) * (-pkin(5) * t268 + qJ(6) * t269) + t360 * t431 + t390 * t236) * MDP(27); (-t275 + t424) * MDP(24) + t233 * MDP(25) + (-t295 ^ 2 - t438) * MDP(26) + (t335 - t428 - t437) * MDP(27);];
tau  = t1;
