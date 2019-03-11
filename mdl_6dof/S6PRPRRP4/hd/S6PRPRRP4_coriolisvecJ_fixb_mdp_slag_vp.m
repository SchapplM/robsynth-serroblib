% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:25
% EndTime: 2019-03-08 20:12:33
% DurationCPUTime: 4.56s
% Computational Cost: add. (3824->367), mult. (9748->483), div. (0->0), fcn. (7647->10), ass. (0->147)
t338 = sin(pkin(11));
t340 = cos(pkin(11));
t343 = sin(qJ(4));
t426 = cos(qJ(4));
t362 = -t343 * t338 + t426 * t340;
t423 = pkin(8) + qJ(3);
t326 = t423 * t338;
t327 = t423 * t340;
t430 = -t426 * t326 - t343 * t327;
t268 = t362 * qJD(3) + qJD(4) * t430;
t339 = sin(pkin(6));
t346 = cos(qJ(2));
t408 = t339 * t346;
t351 = t362 * t408;
t294 = qJD(1) * t351;
t440 = -t268 + t294;
t344 = sin(qJ(2));
t395 = qJD(1) * t339;
t379 = t344 * t395;
t324 = qJD(2) * qJ(3) + t379;
t341 = cos(pkin(6));
t394 = qJD(1) * t341;
t330 = t340 * t394;
t422 = pkin(8) * qJD(2);
t288 = t330 + (-t324 - t422) * t338;
t300 = t340 * t324 + t338 * t394;
t289 = t340 * t422 + t300;
t253 = t343 * t288 + t426 * t289;
t439 = qJD(4) * t253;
t438 = t362 * qJD(2);
t316 = t362 * qJD(4);
t350 = qJD(2) * t316;
t436 = qJD(4) * qJD(5) + t350;
t322 = t426 * t338 + t343 * t340;
t309 = qJD(5) - t438;
t396 = t338 ^ 2 + t340 ^ 2;
t435 = MDP(7) * t396;
t378 = t346 * t395;
t320 = (qJD(3) + t378) * qJD(2);
t434 = t362 * t320;
t431 = t426 * t288 - t343 * t289;
t238 = qJD(4) * t431 + t434;
t249 = qJD(4) * pkin(9) + t253;
t333 = -pkin(3) * t340 - pkin(2);
t370 = qJD(3) - t378;
t308 = qJD(2) * t333 + t370;
t315 = qJD(2) * t322;
t257 = -pkin(4) * t438 - pkin(9) * t315 + t308;
t317 = t322 * qJD(4);
t307 = qJD(2) * t317;
t393 = qJD(2) * t344;
t377 = t339 * t393;
t328 = qJD(1) * t377;
t264 = t307 * pkin(4) - pkin(9) * t350 + t328;
t342 = sin(qJ(5));
t345 = cos(qJ(5));
t391 = qJD(5) * t345;
t392 = qJD(5) * t342;
t373 = t342 * t238 + t249 * t391 + t257 * t392 - t345 * t264;
t425 = pkin(5) * t307;
t222 = t373 - t425;
t232 = t249 * t345 + t257 * t342;
t227 = qJ(6) * t309 + t232;
t418 = t227 * t309;
t433 = -t222 + t418;
t292 = -t343 * t326 + t327 * t426;
t352 = t322 * t408;
t398 = -qJD(1) * t352 + qJD(3) * t322 + qJD(4) * t292;
t280 = pkin(4) * t317 - pkin(9) * t316;
t281 = -pkin(4) * t362 - pkin(9) * t322 + t333;
t432 = -t281 * t391 + t292 * t392 + t440 * t345 + (-t280 + t379) * t342;
t397 = t342 * t281 + t345 * t292;
t248 = -qJD(4) * pkin(4) - t431;
t295 = -t345 * qJD(4) + t315 * t342;
t297 = qJD(4) * t342 + t315 * t345;
t235 = t295 * pkin(5) - t297 * qJ(6) + t248;
t424 = pkin(9) * t307;
t429 = t235 * t309 - t424;
t428 = t297 ^ 2;
t427 = t309 ^ 2;
t421 = qJD(2) * pkin(2);
t420 = qJ(6) * t307;
t239 = t322 * t320 + t439;
t265 = t315 * t392 - t436 * t345;
t266 = t315 * t391 + t436 * t342;
t223 = pkin(5) * t266 + qJ(6) * t265 - qJD(6) * t297 + t239;
t419 = t223 * t342;
t417 = t232 * t309;
t416 = t265 * t342;
t415 = t295 * t438;
t414 = t295 * t315;
t413 = t297 * t295;
t375 = t297 * t309;
t412 = t297 * t315;
t411 = t309 * t342;
t410 = t322 * t345;
t409 = t339 * t344;
t302 = t342 * t307;
t304 = t345 * t307;
t404 = qJ(6) * t317 - qJD(6) * t362 - t432;
t270 = t294 * t342 - t345 * t379;
t403 = -pkin(5) * t317 + qJD(5) * t397 + t268 * t342 - t280 * t345 - t270;
t371 = pkin(5) * t342 - qJ(6) * t345;
t372 = pkin(5) * t345 + qJ(6) * t342;
t402 = t371 * t316 + (qJD(5) * t372 - qJD(6) * t345) * t322 + t398;
t401 = qJD(6) * t342 - t309 * t371 + t253;
t279 = pkin(4) * t315 - pkin(9) * t438;
t400 = t342 * t279 + t345 * t431;
t399 = -t342 * t266 - t295 * t391;
t231 = -t249 * t342 + t257 * t345;
t390 = qJD(6) - t231;
t386 = pkin(9) * t392;
t376 = t396 * t320;
t226 = -pkin(5) * t309 + t390;
t369 = t226 * t345 - t227 * t342;
t368 = (-t324 * t338 + t330) * t338 - t300 * t340;
t357 = t345 * t238 - t249 * t392 + t257 * t391 + t342 * t264;
t221 = qJD(6) * t309 + t357 + t420;
t367 = t226 * t309 + t221;
t366 = t302 + (-t345 * t438 + t391) * t309;
t365 = -t309 * t392 + t411 * t438 + t304;
t311 = -t338 * t409 + t340 * t341;
t312 = t338 * t341 + t340 * t409;
t273 = t343 * t311 + t312 * t426;
t261 = t273 * t342 + t345 * t408;
t262 = t273 * t345 - t342 * t408;
t363 = t311 * t426 - t343 * t312;
t361 = t316 * t342 + t322 * t391;
t360 = -t316 * t345 + t322 * t392;
t359 = t248 * t309 - t424;
t358 = t235 * t297 + t373;
t349 = qJD(2) * t351;
t347 = qJD(2) ^ 2;
t325 = -pkin(4) - t372;
t323 = t370 - t421;
t260 = pkin(5) * t297 + qJ(6) * t295;
t254 = t322 * t371 - t430;
t251 = qJD(2) * t352 + qJD(4) * t273;
t250 = qJD(4) * t363 + t349;
t245 = pkin(5) * t362 - t281 * t345 + t292 * t342;
t244 = -qJ(6) * t362 + t397;
t242 = t295 * t309 - t265;
t234 = -pkin(5) * t315 - t279 * t345 + t342 * t431;
t233 = qJ(6) * t315 + t400;
t230 = qJD(5) * t262 + t250 * t342 - t345 * t377;
t229 = -qJD(5) * t261 + t250 * t345 + t342 * t377;
t1 = [((-t311 * t338 + t312 * t340) * t320 + (t323 * t344 + (-t368 - t379) * t346) * t339 * qJD(2)) * MDP(8) + (-qJD(4) * t251 + (-t307 * t346 - t393 * t438) * t339) * MDP(14) + (t315 * t377 + (-t250 - t349) * qJD(4)) * MDP(15) + (-t229 * t295 + t230 * t297 - t261 * t265 - t262 * t266) * MDP(24) + (t221 * t262 + t222 * t261 - t223 * t363 + t226 * t230 + t227 * t229 + t235 * t251) * MDP(26) + (MDP(21) + MDP(23)) * (-t230 * t309 + t251 * t295 - t261 * t307 - t266 * t363) + (-MDP(22) + MDP(25)) * (t229 * t309 - t251 * t297 + t262 * t307 - t265 * t363) + ((-MDP(4) + t435) * t346 + (-MDP(5) * t340 + MDP(6) * t338 - MDP(3)) * t344) * t339 * t347; (qJD(2) * t370 * t396 + t376) * MDP(7) + (-t368 * qJD(3) + qJ(3) * t376 + (t368 * t346 + (-t323 - t421) * t344) * t395) * MDP(8) + (t315 * t316 + t322 * t350) * MDP(9) + (-t322 * t307 - t315 * t317 + t316 * t438 + t350 * t362) * MDP(10) + (t307 * t333 + t308 * t317) * MDP(14) + (t308 * t316 - t315 * t379 + t322 * t328 + t333 * t350) * MDP(15) + (-t265 * t410 - t297 * t360) * MDP(16) + ((-t295 * t345 - t297 * t342) * t316 + (t416 - t266 * t345 + (t295 * t342 - t297 * t345) * qJD(5)) * t322) * MDP(17) + (t265 * t362 + t297 * t317 + t304 * t322 - t309 * t360) * MDP(18) + (t266 * t362 - t295 * t317 - t302 * t322 - t309 * t361) * MDP(19) + (-t307 * t362 + t309 * t317) * MDP(20) + (t373 * t362 + t231 * t317 - t430 * t266 + t270 * t309 + t398 * t295 + ((-qJD(5) * t292 + t280) * t309 + t281 * t307 + t248 * qJD(5) * t322) * t345 + ((-qJD(5) * t281 - t268) * t309 - t292 * t307 + t239 * t322 + t248 * t316) * t342) * MDP(21) + (-t232 * t317 + t239 * t410 - t360 * t248 + t265 * t430 + t398 * t297 - t397 * t307 + t309 * t432 + t357 * t362) * MDP(22) + (t222 * t362 - t226 * t317 + t235 * t361 - t245 * t307 + t254 * t266 + t295 * t402 - t309 * t403 + t322 * t419) * MDP(23) + (-t244 * t266 - t245 * t265 + t369 * t316 + t403 * t297 - t404 * t295 + (-t221 * t342 + t222 * t345 + (-t226 * t342 - t227 * t345) * qJD(5)) * t322) * MDP(24) + (-t221 * t362 - t223 * t410 + t227 * t317 + t235 * t360 + t244 * t307 + t254 * t265 - t297 * t402 + t309 * t404) * MDP(25) + (t221 * t244 + t222 * t245 + t223 * t254 + t226 * t403 + t227 * t404 + t235 * t402) * MDP(26) + (t316 * MDP(11) - t317 * MDP(12) - t398 * MDP(14) + t440 * MDP(15)) * qJD(4); -t347 * t435 + (qJD(2) * t368 + t328) * MDP(8) + (t365 - t414) * MDP(21) + (-t345 * t427 - t302 - t412) * MDP(22) + (-t309 * t411 + t304 - t414) * MDP(23) + ((t265 + t415) * t345 + t342 * t375 + t399) * MDP(24) + (t366 + t412) * MDP(25) + (-t235 * t315 + t367 * t342 + t345 * t433) * MDP(26) + (0.2e1 * t315 * MDP(14) + 0.2e1 * t438 * MDP(15)) * qJD(4); -t438 ^ 2 * MDP(10) + (-t239 + t439) * MDP(14) + (-t308 * t438 - t434) * MDP(15) + (t345 * t375 - t416) * MDP(16) + ((-t265 + t415) * t345 - t297 * t411 + t399) * MDP(17) + (t366 - t412) * MDP(18) + (t365 + t414) * MDP(19) + (-pkin(4) * t266 - t253 * t295 + (-t239 + (-pkin(9) * qJD(5) - t279) * t309) * t345 + (t309 * t431 + t359) * t342) * MDP(21) + (pkin(4) * t265 + t239 * t342 - t253 * t297 + (t386 + t400) * t309 + t359 * t345) * MDP(22) + (-t223 * t345 + t266 * t325 + (-pkin(9) * t391 + t234) * t309 - t401 * t295 + t429 * t342) * MDP(23) + (t233 * t295 - t234 * t297 + ((qJD(5) * t297 - t266) * pkin(9) + t367) * t345 + ((qJD(5) * t295 - t265) * pkin(9) - t433) * t342) * MDP(24) + (-t419 + t265 * t325 + (-t233 - t386) * t309 + t401 * t297 - t429 * t345) * MDP(25) + (t223 * t325 - t226 * t234 - t227 * t233 - t401 * t235 + (qJD(5) * t369 + t221 * t345 + t222 * t342) * pkin(9)) * MDP(26) + (MDP(10) * t315 - t308 * MDP(14) - t309 * MDP(20) - t231 * MDP(21) + t232 * MDP(22) + t226 * MDP(23) - t227 * MDP(25) - MDP(9) * t438) * t315; MDP(16) * t413 + (-t295 ^ 2 + t428) * MDP(17) + t242 * MDP(18) + (-t266 + t375) * MDP(19) + t307 * MDP(20) + (-t248 * t297 - t373 + t417) * MDP(21) + (t231 * t309 + t248 * t295 - t357) * MDP(22) + (-t260 * t295 - t358 + t417 + 0.2e1 * t425) * MDP(23) + (pkin(5) * t265 - qJ(6) * t266 + (t227 - t232) * t297 + (t226 - t390) * t295) * MDP(24) + (0.2e1 * t420 - t235 * t295 + t260 * t297 + (0.2e1 * qJD(6) - t231) * t309 + t357) * MDP(25) + (-pkin(5) * t222 + qJ(6) * t221 - t226 * t232 + t227 * t390 - t235 * t260) * MDP(26); (-qJD(4) * t315 + t413) * MDP(23) + t242 * MDP(24) + (-t427 - t428) * MDP(25) + (t358 - t418 - t425) * MDP(26);];
tauc  = t1;
