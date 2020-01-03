% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:35
% EndTime: 2019-12-31 20:07:41
% DurationCPUTime: 3.35s
% Computational Cost: add. (2830->358), mult. (7079->484), div. (0->0), fcn. (4769->6), ass. (0->155)
t354 = sin(pkin(8));
t357 = sin(qJ(2));
t413 = qJD(1) * t357;
t389 = t354 * t413;
t355 = cos(pkin(8));
t408 = qJD(2) * t355;
t317 = -t389 + t408;
t393 = t355 * t413;
t409 = qJD(2) * t354;
t318 = t393 + t409;
t356 = sin(qJ(4));
t358 = cos(qJ(4));
t272 = -t358 * t317 + t318 * t356;
t359 = cos(qJ(2));
t412 = qJD(1) * t359;
t344 = -qJD(4) + t412;
t442 = t272 * t344;
t425 = t358 * t355;
t320 = t354 * t356 - t425;
t404 = qJD(4) * t358;
t405 = qJD(4) * t356;
t435 = -t354 * t405 + t355 * t404;
t416 = t320 * t412 + t435;
t321 = t354 * t358 + t355 * t356;
t308 = t321 * qJD(4);
t368 = t321 * t359;
t415 = -qJD(1) * t368 + t308;
t400 = qJD(1) * qJD(2);
t441 = -0.2e1 * t400;
t440 = t357 * MDP(4);
t439 = (t357 ^ 2 - t359 ^ 2) * MDP(5);
t378 = pkin(2) * t357 - qJ(3) * t359;
t323 = t378 * qJD(1);
t288 = pkin(6) * t389 + t355 * t323;
t427 = t355 * t359;
t372 = pkin(3) * t357 - pkin(7) * t427;
t265 = qJD(1) * t372 + t288;
t309 = t354 * t323;
t428 = t355 * t357;
t429 = t354 * t359;
t370 = -pkin(6) * t428 - pkin(7) * t429;
t277 = qJD(1) * t370 + t309;
t433 = pkin(7) + qJ(3);
t331 = t433 * t354;
t332 = t433 * t355;
t373 = -t331 * t358 - t332 * t356;
t438 = -qJD(3) * t320 + qJD(4) * t373 - t356 * t265 - t358 * t277;
t286 = -t331 * t356 + t332 * t358;
t437 = qJD(3) * t321 + qJD(4) * t286 + t265 * t358 - t277 * t356;
t328 = -pkin(2) * t359 - qJ(3) * t357 - pkin(1);
t316 = t355 * t328;
t278 = -pkin(7) * t428 + t316 + (-pkin(6) * t354 - pkin(3)) * t359;
t342 = pkin(6) * t427;
t293 = t354 * t328 + t342;
t430 = t354 * t357;
t284 = -pkin(7) * t430 + t293;
t436 = t356 * t278 + t358 * t284;
t399 = MDP(20) + MDP(22);
t398 = MDP(21) - MDP(24);
t387 = t359 * t400;
t381 = t354 * t387;
t406 = qJD(2) * t359;
t382 = t406 * t425;
t417 = qJD(1) * t382 + t317 * t404;
t242 = t356 * (qJD(4) * t318 + t381) - t417;
t305 = qJD(2) * t378 - qJD(3) * t357;
t407 = qJD(2) * t357;
t395 = pkin(6) * t407;
t282 = t355 * t305 + t354 * t395;
t365 = t372 * qJD(2);
t258 = t365 + t282;
t297 = t354 * t305;
t266 = qJD(2) * t370 + t297;
t434 = -qJD(4) * t436 + t258 * t358 - t266 * t356;
t432 = qJD(2) * pkin(2);
t311 = t328 * qJD(1);
t349 = pkin(6) * t412;
t334 = qJD(2) * qJ(3) + t349;
t280 = t354 * t311 + t355 * t334;
t254 = pkin(7) * t317 + t280;
t431 = t254 * t356;
t360 = qJD(2) ^ 2;
t426 = t357 * t360;
t424 = t359 * t360;
t361 = qJD(1) ^ 2;
t423 = t359 * t361;
t422 = qJ(5) * t413 - t438;
t421 = pkin(4) * t413 + t437;
t396 = pkin(3) * t412;
t312 = t354 * t396 + t349;
t420 = -pkin(4) * t415 + qJ(5) * t416 + qJD(5) * t321 + t312;
t296 = t305 * qJD(1);
t348 = pkin(6) * t413;
t326 = (qJD(3) - t348) * qJD(2);
t264 = t354 * t296 + t355 * t326;
t343 = pkin(6) * t387;
t304 = pkin(3) * t381 + t343;
t350 = pkin(6) * t406;
t392 = t354 * t406;
t313 = pkin(3) * t392 + t350;
t324 = pkin(3) * t430 + t357 * pkin(6);
t411 = qJD(2) * t373;
t410 = qJD(2) * t286;
t403 = qJD(5) * t344;
t402 = t357 * MDP(19);
t279 = t355 * t311 - t334 * t354;
t249 = -pkin(7) * t318 + t279 - t396;
t230 = t249 * t358 - t431;
t401 = qJD(5) - t230;
t397 = pkin(6) * t429;
t263 = t355 * t296 - t326 * t354;
t246 = qJD(1) * t365 + t263;
t253 = -pkin(7) * t381 + t264;
t394 = t356 * t246 + t249 * t404 + t358 * t253;
t347 = -pkin(3) * t355 - pkin(2);
t388 = t357 * t400;
t386 = pkin(1) * t441;
t385 = qJD(3) - t432;
t384 = t399 * t356;
t383 = t358 * t246 - t249 * t405 - t356 * t253 - t254 * t404;
t327 = t348 + t385;
t380 = -t327 + t385;
t231 = t249 * t356 + t254 * t358;
t375 = t278 * t358 - t284 * t356;
t275 = t317 * t356 + t318 * t358;
t371 = -t230 * t344 - t394;
t367 = -t254 * t405 + t394;
t366 = t356 * t258 + t358 * t266 + t278 * t404 - t284 * t405;
t287 = -pkin(3) * t317 + t327;
t363 = qJD(2) * t368;
t223 = -pkin(4) * t388 - t383;
t362 = -t356 * t398 + t358 * t399;
t243 = qJD(1) * t363 + qJD(4) * t275;
t224 = pkin(4) * t243 + qJ(5) * t242 - qJD(5) * t275 + t304;
t302 = t320 * t357;
t301 = t321 * t357;
t292 = t316 - t397;
t289 = -pkin(6) * t393 + t309;
t283 = -t355 * t395 + t297;
t270 = pkin(4) * t320 - qJ(5) * t321 + t347;
t261 = t435 * t357 + t363;
t260 = t308 * t357 + t356 * t392 - t382;
t250 = pkin(4) * t301 + qJ(5) * t302 + t324;
t238 = pkin(4) * t275 + qJ(5) * t272;
t237 = pkin(4) * t359 - t375;
t236 = -qJ(5) * t359 + t436;
t233 = pkin(4) * t272 - qJ(5) * t275 + t287;
t232 = -t242 - t442;
t229 = -qJ(5) * t344 + t231;
t228 = pkin(4) * t344 + t401;
t227 = pkin(4) * t261 + qJ(5) * t260 + qJD(5) * t302 + t313;
t226 = -pkin(4) * t407 - t434;
t225 = qJ(5) * t407 - qJD(5) * t359 + t366;
t222 = qJ(5) * t388 + t367 - t403;
t1 = [0.2e1 * t387 * t440 + t439 * t441 + MDP(6) * t424 - MDP(7) * t426 + (-pkin(6) * t424 + t357 * t386) * MDP(9) + (pkin(6) * t426 + t359 * t386) * MDP(10) + ((-qJD(1) * t282 - t263) * t359 + ((-pkin(6) * t317 + t327 * t354) * t359 + (t279 + (t292 + 0.2e1 * t397) * qJD(1)) * t357) * qJD(2)) * MDP(11) + ((qJD(1) * t283 + t264) * t359 + ((pkin(6) * t318 + t327 * t355) * t359 + (-t280 + (-t293 + 0.2e1 * t342) * qJD(1)) * t357) * qJD(2)) * MDP(12) + (-t282 * t318 + t283 * t317 + (-t263 * t355 - t264 * t354) * t357 + (-t279 * t355 - t280 * t354 + (-t292 * t355 - t293 * t354) * qJD(1)) * t406) * MDP(13) + (t263 * t292 + t264 * t293 + t279 * t282 + t280 * t283 + (t327 + t348) * t350) * MDP(14) + (t242 * t302 - t260 * t275) * MDP(15) + (t242 * t301 + t243 * t302 + t260 * t272 - t261 * t275) * MDP(16) + (t242 * t359 + t260 * t344 + (-qJD(1) * t302 + t275) * t407) * MDP(17) + (t243 * t359 + t261 * t344 + (-qJD(1) * t301 - t272) * t407) * MDP(18) + (-t344 - t412) * qJD(2) * t402 + (-t434 * t344 - t383 * t359 + t313 * t272 + t324 * t243 + t304 * t301 + t287 * t261 + (qJD(1) * t375 + t230) * t407) * MDP(20) + (t366 * t344 + t367 * t359 + t313 * t275 - t324 * t242 - t304 * t302 - t287 * t260 + (-qJD(1) * t436 - t231) * t407) * MDP(21) + (t223 * t359 + t224 * t301 + t226 * t344 + t227 * t272 + t233 * t261 + t243 * t250 + (-qJD(1) * t237 - t228) * t407) * MDP(22) + (-t222 * t301 - t223 * t302 - t225 * t272 + t226 * t275 - t228 * t260 - t229 * t261 - t236 * t243 - t237 * t242) * MDP(23) + (-t222 * t359 + t224 * t302 - t225 * t344 - t227 * t275 + t233 * t260 + t242 * t250 + (qJD(1) * t236 + t229) * t407) * MDP(24) + (t222 * t236 + t223 * t237 + t224 * t250 + t225 * t229 + t226 * t228 + t227 * t233) * MDP(25); -t423 * t440 + t361 * t439 + (t288 * t318 - t289 * t317 + (qJD(3) * t317 + t279 * t412 + t264) * t355 + (qJD(3) * t318 + t280 * t412 - t263) * t354) * MDP(13) + (-t279 * t288 - t280 * t289 + (-t279 * t354 + t280 * t355) * qJD(3) + (-t263 * t354 + t264 * t355) * qJ(3) + (-t327 - t432) * t349) * MDP(14) + (-t242 * t321 + t275 * t416) * MDP(15) + (t242 * t320 - t243 * t321 - t272 * t416 - t275 * t415) * MDP(16) + (-t416 * t344 + (qJD(2) * t321 - t275) * t413) * MDP(17) + (t415 * t344 + (-qJD(2) * t320 + t272) * t413) * MDP(18) + (t347 * t243 - t312 * t272 + t304 * t320 + t437 * t344 + t415 * t287 + (-t230 + t411) * t413) * MDP(20) + (-t347 * t242 - t312 * t275 + t304 * t321 + t438 * t344 + t416 * t287 + (t231 - t410) * t413) * MDP(21) + (t224 * t320 + t243 * t270 + t421 * t344 - t420 * t272 + t415 * t233 + (t228 + t411) * t413) * MDP(22) + (-t222 * t320 + t223 * t321 + t228 * t416 - t229 * t415 + t242 * t373 - t243 * t286 + t272 * t422 + t275 * t421) * MDP(23) + (-t224 * t321 + t242 * t270 + t422 * t344 + t420 * t275 - t416 * t233 + (-t229 + t410) * t413) * MDP(24) + (t222 * t286 - t223 * t373 + t224 * t270 + t228 * t421 - t229 * t422 - t233 * t420) * MDP(25) + (t361 * t357 * MDP(9) + MDP(10) * t423) * pkin(1) + (((-qJ(3) * t409 - t279) * t357 + (t288 + (t317 - t408) * pkin(6) + t380 * t354) * t359) * MDP(11) + ((-qJ(3) * t408 + t280) * t357 + (-t289 + (-t318 + t409) * pkin(6) + t380 * t355) * t359) * MDP(12) + t344 * t402) * qJD(1); (-t317 ^ 2 - t318 ^ 2) * MDP(13) + (t279 * t318 - t280 * t317 + t343) * MDP(14) - t272 ^ 2 * MDP(23) + (t229 * t272 + t224) * MDP(25) + (-t275 * MDP(23) - t228 * MDP(25) - t344 * t399) * t275 + (t317 * t384 + t318 * t362) * qJD(4) + (-t318 * MDP(11) - t317 * MDP(12) + ((MDP(12) + t384) * t355 + (MDP(11) + t362) * t354) * qJD(2)) * t412 - t398 * (-t417 - t442); t232 * MDP(17) + t371 * MDP(21) + (pkin(4) * t242 - qJ(5) * t243) * MDP(23) + (-t371 - 0.2e1 * t403) * MDP(24) + (-pkin(4) * t223 + qJ(5) * t222 - t228 * t231 + t229 * t401 - t233 * t238) * MDP(25) + (-t275 * MDP(18) + t398 * t431) * qJD(4) + (-t344 * MDP(18) - t287 * MDP(20) - t233 * MDP(22) + (t229 - t231) * MDP(23) + t238 * MDP(24) + MDP(16) * t275) * t275 + (-MDP(18) * t368 + (0.2e1 * pkin(4) * MDP(22) + 0.2e1 * qJ(5) * MDP(24) + MDP(19)) * t357) * t400 + (t275 * MDP(15) + t287 * MDP(21) - t238 * MDP(22) + (t228 - t401) * MDP(23) - t233 * MDP(24) - MDP(16) * t272) * t272 + t399 * (-t231 * t344 + t383); (t272 * t275 - t388) * MDP(22) + t232 * MDP(23) + (-t275 ^ 2 - t344 ^ 2) * MDP(24) + (t229 * t344 + t233 * t275 + t223) * MDP(25);];
tauc = t1;
