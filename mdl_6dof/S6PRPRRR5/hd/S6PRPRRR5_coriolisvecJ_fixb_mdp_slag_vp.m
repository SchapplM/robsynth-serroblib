% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:41
% EndTime: 2019-03-08 20:43:48
% DurationCPUTime: 2.88s
% Computational Cost: add. (2227->310), mult. (4983->428), div. (0->0), fcn. (3655->10), ass. (0->144)
t307 = qJD(4) + qJD(5);
t316 = cos(qJ(6));
t313 = sin(qJ(5));
t317 = cos(qJ(5));
t318 = cos(qJ(4));
t376 = qJD(2) * t318;
t357 = t317 * t376;
t314 = sin(qJ(4));
t377 = qJD(2) * t314;
t273 = t313 * t377 - t357;
t312 = sin(qJ(6));
t391 = t273 * t312;
t263 = -t316 * t307 - t391;
t281 = t313 * t318 + t314 * t317;
t274 = t281 * qJD(2);
t407 = qJD(6) + t274;
t410 = t263 * t407;
t335 = t273 * t316 - t307 * t312;
t409 = t335 * t407;
t373 = qJD(5) * t313;
t375 = qJD(4) * t314;
t408 = -t313 * t375 - t314 * t373;
t319 = cos(qJ(2));
t310 = sin(pkin(6));
t380 = qJD(1) * t310;
t363 = t319 * t380;
t339 = qJD(3) - t363;
t406 = t274 * t307;
t315 = sin(qJ(2));
t360 = qJD(2) * t310 * t315;
t347 = qJD(1) * t360;
t320 = -pkin(2) - pkin(8);
t272 = t320 * qJD(2) + t339;
t354 = pkin(9) * qJD(2) - t272;
t311 = cos(pkin(6));
t379 = qJD(1) * t311;
t362 = t314 * t379;
t241 = t314 * t347 + (-t354 * t318 - t362) * qJD(4);
t253 = -pkin(9) * t376 + t318 * t272 - t362;
t248 = qJD(4) * pkin(4) + t253;
t405 = (qJD(5) * t248 + t241) * t317;
t404 = (t314 ^ 2 - t318 ^ 2) * MDP(9);
t374 = qJD(4) * t318;
t341 = -pkin(4) * t374 - t339;
t403 = t314 * MDP(13) + t318 * MDP(14);
t283 = t318 * t347;
t361 = t318 * t379;
t242 = t283 + (t354 * t314 - t361) * qJD(4);
t254 = -pkin(9) * t377 + t272 * t314 + t361;
t352 = t242 * t313 - t254 * t373;
t212 = t352 + t405;
t392 = t254 * t317;
t228 = t248 * t313 + t392;
t353 = t241 * t313 - t317 * t242;
t213 = t228 * qJD(5) + t353;
t393 = t254 * t313;
t227 = t248 * t317 - t393;
t225 = -pkin(5) * t307 - t227;
t364 = t315 * t380;
t378 = qJD(2) * qJ(3);
t282 = t364 + t378;
t268 = pkin(4) * t377 + t282;
t239 = pkin(5) * t274 + pkin(10) * t273 + t268;
t387 = t317 * t318;
t340 = t307 * t387;
t383 = t408 * qJD(2);
t250 = qJD(2) * t340 + t383;
t280 = t313 * t314 - t387;
t299 = t314 * pkin(4) + qJ(3);
t252 = pkin(5) * t281 + pkin(10) * t280 + t299;
t372 = qJD(5) * t317;
t255 = -t313 * t374 - t314 * t372 - t317 * t375 - t318 * t373;
t401 = pkin(9) - t320;
t284 = t401 * t314;
t285 = t401 * t318;
t260 = -t284 * t317 - t285 * t313;
t259 = -t284 * t313 + t285 * t317;
t276 = t401 * t375;
t277 = qJD(4) * t285;
t386 = t259 * qJD(5) - t276 * t313 + t277 * t317 + t281 * t364;
t402 = -(qJD(6) * t239 + t212) * t281 + t225 * t255 + (-qJD(6) * t252 + t386) * t407 - t213 * t280 - t260 * t250;
t400 = qJD(2) * pkin(2);
t399 = t225 * t274;
t398 = t225 * t280;
t371 = qJD(6) * t312;
t370 = qJD(6) * t316;
t384 = t307 * t370 - t316 * t406;
t232 = t273 * t371 + t384;
t397 = t232 * t280;
t396 = t232 * t312;
t395 = t406 * t312;
t394 = t252 * t250;
t390 = t310 * t319;
t389 = t312 * t250;
t388 = t316 * t250;
t385 = t260 * qJD(5) - t276 * t317 - t277 * t313 - t280 * t364;
t366 = qJD(2) * qJD(4);
t365 = pkin(4) * t376;
t356 = t318 * t366;
t355 = -pkin(4) * t307 - t248;
t351 = t316 * t407;
t251 = -pkin(5) * t273 + pkin(10) * t274;
t301 = pkin(4) * t313 + pkin(10);
t348 = qJD(6) * t301 + t251 + t365;
t229 = t253 * t313 + t392;
t346 = pkin(4) * t373 - t229;
t230 = t253 * t317 - t393;
t345 = -pkin(4) * t372 + t230;
t226 = pkin(10) * t307 + t228;
t218 = t226 * t316 + t239 * t312;
t344 = t213 * t312 - t218 * t273 + t225 * t370;
t256 = t340 + t408;
t343 = pkin(5) * t256 - pkin(10) * t255 - t341;
t342 = -t282 + t364;
t338 = -t250 * t301 + t399;
t337 = t226 * t312 - t239 * t316;
t270 = -t311 * t314 - t318 * t390;
t333 = -t311 * t318 + t314 * t390;
t336 = t270 * t317 + t313 * t333;
t245 = t270 * t313 - t317 * t333;
t334 = -t274 * MDP(20) + t273 * MDP(21);
t332 = -t213 * t316 + t225 * t371 - t273 * t337;
t331 = t268 * t273 - t353;
t330 = t268 * t274 - t352;
t329 = -t370 * t407 - t389;
t328 = t371 * t407 - t388;
t327 = t316 * t255 + t280 * t371;
t278 = (qJD(3) + t363) * qJD(2);
t326 = t342 - t378;
t266 = pkin(4) * t356 + t278;
t321 = qJD(4) ^ 2;
t325 = t339 * qJD(2) - t320 * t321 + t278;
t233 = -t335 * qJD(6) - t395;
t324 = ((t232 - t410) * t316 + (-t233 + t409) * t312) * MDP(23) + (-t335 * t351 + t396) * MDP(22) + (-t312 * t407 ^ 2 - t263 * t273 + t388) * MDP(25) + (-t273 * t335 + t351 * t407 + t389) * MDP(24) + (-t383 + (-t273 - t357) * t307) * MDP(18) + (t273 ^ 2 - t274 ^ 2) * MDP(16) + (-MDP(15) * t274 + MDP(26) * t407) * t273;
t322 = qJD(2) ^ 2;
t302 = -pkin(4) * t317 - pkin(5);
t279 = t339 - t400;
t258 = t333 * qJD(4) + t318 * t360;
t257 = t270 * qJD(4) + t314 * t360;
t222 = pkin(5) * t250 + pkin(10) * t406 + t266;
t221 = t316 * t222;
t220 = t245 * qJD(5) + t257 * t313 - t258 * t317;
t219 = t336 * qJD(5) + t257 * t317 + t258 * t313;
t1 = [((-t219 * t312 - t245 * t370) * t407 - t245 * t389 + t220 * t263 - t336 * t233) * MDP(27) + (-(t219 * t316 - t245 * t371) * t407 - t245 * t388 - t220 * t335 - t336 * t232) * MDP(28) + (-MDP(20) * t220 - MDP(21) * t219) * t307 + (MDP(13) * t258 - MDP(14) * t257) * qJD(4) + (((-MDP(4) + MDP(6) + t403) * t322 + (MDP(7) * t282 + (t316 * MDP(27) - t312 * MDP(28)) * t407 - t334) * qJD(2)) * t319 + (t278 * MDP(7) + t250 * MDP(20) - t406 * MDP(21) - t328 * MDP(27) + t329 * MDP(28) + (-MDP(3) + MDP(5)) * t322 + ((t318 * MDP(13) - t314 * MDP(14)) * qJD(4) + (t279 - t363) * MDP(7)) * qJD(2)) * t315) * t310; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t278 + qJD(3) * t282 + (-t282 * t319 + (-t279 - t400) * t315) * t380) * MDP(7) + 0.2e1 * t366 * t404 - t321 * t318 * MDP(11) - t326 * t374 * MDP(13) + (t325 * t318 + t326 * t375) * MDP(14) + (-t255 * t273 + t280 * t406) * MDP(15) + (t250 * t280 - t255 * t274 + t256 * t273 + t281 * t406) * MDP(16) + (t250 * t299 + t256 * t268 + t266 * t281 - t341 * t274) * MDP(20) + (t255 * t268 - t266 * t280 + t341 * t273 - t299 * t406) * MDP(21) + (-t316 * t397 - t327 * t335) * MDP(22) + ((-t263 * t316 + t312 * t335) * t255 + (t396 + t233 * t316 + (-t263 * t312 - t316 * t335) * qJD(6)) * t280) * MDP(23) + (t232 * t281 - t256 * t335 - t280 * t388 + t327 * t407) * MDP(24) + (t280 * t389 - t233 * t281 - t256 * t263 + (-t312 * t255 + t280 * t370) * t407) * MDP(25) + (t250 * t281 + t256 * t407) * MDP(26) + (-t337 * t256 + t221 * t281 + t259 * t233 + t385 * t263 + (t394 + t343 * t407 + (-t226 * t281 - t260 * t407 - t398) * qJD(6)) * t316 + t402 * t312) * MDP(27) + (-t218 * t256 + t259 * t232 - t385 * t335 + (-t394 - (-qJD(6) * t226 + t222) * t281 + qJD(6) * t398 + (qJD(6) * t260 - t343) * t407) * t312 + t402 * t316) * MDP(28) + (-t321 * MDP(10) + t325 * MDP(13) - 0.2e1 * MDP(8) * t356) * t314 + (t255 * MDP(17) - t256 * MDP(18) - t385 * MDP(20) + t386 * MDP(21)) * t307; -t322 * MDP(6) + (t233 * t280 - t255 * t263 - t281 * t389) * MDP(27) + (t255 * t335 - t281 * t388 + t397) * MDP(28) + (MDP(20) * t255 - MDP(21) * t256) * t307 + (t342 * MDP(7) + t334) * qJD(2) + ((-qJD(2) * t316 - t256 * t312 - t281 * t370) * MDP(27) + (qJD(2) * t312 - t256 * t316 + t281 * t371) * MDP(28)) * t407 + t403 * (-t321 - t322); t324 + (-t274 * t365 + t229 * t307 + (t355 * t313 - t392) * qJD(5) + t331) * MDP(20) + (t273 * t365 + t230 * t307 + (t355 * qJD(5) - t241) * t317 + t330) * MDP(21) + (t302 * t233 + t338 * t312 + t346 * t263 + (t345 * t312 - t348 * t316) * t407 + t332) * MDP(27) + (t302 * t232 + t338 * t316 - t346 * t335 + (t348 * t312 + t345 * t316) * t407 + t344) * MDP(28) - t342 * t377 * MDP(14) + (-t282 * t376 + t283) * MDP(13) + (t318 * t314 * MDP(8) - t404) * t322; (t331 + (-qJD(5) + t307) * t228) * MDP(20) + (t227 * t307 + t330 - t405) * MDP(21) + (-pkin(5) * t233 - (-t227 * t312 + t251 * t316) * t407 - t228 * t263 + t312 * t399 + t329 * pkin(10) + t332) * MDP(27) + (-pkin(5) * t232 + (t227 * t316 + t251 * t312) * t407 + t228 * t335 + t316 * t399 + t328 * pkin(10) + t344) * MDP(28) + t324; -t335 * t263 * MDP(22) + (-t263 ^ 2 + t335 ^ 2) * MDP(23) + (t384 + t410) * MDP(24) + (t395 - t409) * MDP(25) + t250 * MDP(26) + (-t212 * t312 + t218 * t407 + t225 * t335 + t221) * MDP(27) + (-t212 * t316 - t222 * t312 + t225 * t263 - t337 * t407) * MDP(28) + (MDP(24) * t391 + t335 * MDP(25) - t218 * MDP(27) + t337 * MDP(28)) * qJD(6);];
tauc  = t1;
