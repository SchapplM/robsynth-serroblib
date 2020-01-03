% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:52
% EndTime: 2019-12-31 19:11:00
% DurationCPUTime: 3.93s
% Computational Cost: add. (2760->335), mult. (7244->450), div. (0->0), fcn. (5616->8), ass. (0->146)
t339 = cos(pkin(9));
t345 = cos(qJ(3));
t381 = qJD(1) * t345;
t329 = t339 * t381;
t338 = sin(pkin(9));
t342 = sin(qJ(3));
t382 = qJD(1) * t342;
t370 = t338 * t382;
t308 = t329 - t370;
t412 = qJD(4) + qJD(5);
t424 = t308 - t412;
t409 = pkin(6) + qJ(2);
t322 = t409 * t338;
t315 = qJD(1) * t322;
t323 = t409 * t339;
t316 = qJD(1) * t323;
t280 = -t342 * t315 + t345 * t316;
t423 = qJD(3) * t280;
t314 = t338 * t345 + t339 * t342;
t309 = t314 * qJD(1);
t341 = sin(qJ(4));
t344 = cos(qJ(4));
t375 = t344 * qJD(3);
t291 = t309 * t341 - t375;
t343 = cos(qJ(5));
t293 = qJD(3) * t341 + t309 * t344;
t340 = sin(qJ(5));
t398 = t293 * t340;
t238 = t343 * t291 + t398;
t304 = qJD(4) - t308;
t299 = qJD(5) + t304;
t422 = t238 * t299;
t357 = t291 * t340 - t343 * t293;
t421 = t299 * t357;
t420 = (t338 ^ 2 + t339 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t318 = t340 * t344 + t341 * t343;
t384 = t424 * t318;
t419 = t338 * t381 + t339 * t382;
t380 = qJD(4) * t341;
t395 = t308 * t341;
t418 = t380 - t395;
t331 = -pkin(2) * t339 - pkin(1);
t321 = qJD(1) * t331 + qJD(2);
t254 = -pkin(3) * t308 - pkin(7) * t309 + t321;
t275 = qJD(3) * pkin(7) + t280;
t232 = t254 * t341 + t275 * t344;
t226 = -pkin(8) * t291 + t232;
t378 = qJD(5) * t340;
t224 = t226 * t378;
t413 = -t315 * t345 - t342 * t316;
t274 = -qJD(3) * pkin(3) - t413;
t235 = pkin(4) * t291 + t274;
t417 = t235 * t238 + t224;
t328 = qJD(3) * t329;
t300 = -qJD(3) * t370 + t328;
t249 = qJD(4) * t375 + t344 * t300 - t309 * t380;
t311 = t314 * qJD(3);
t301 = qJD(1) * t311;
t313 = t338 * t342 - t345 * t339;
t348 = t313 * qJD(2);
t242 = -qJD(1) * t348 + qJD(3) * t413;
t262 = pkin(3) * t301 - pkin(7) * t300;
t259 = t344 * t262;
t347 = -qJD(4) * t232 - t242 * t341 + t259;
t213 = pkin(4) * t301 - pkin(8) * t249 + t347;
t250 = t293 * qJD(4) + t300 * t341;
t379 = qJD(4) * t344;
t350 = t344 * t242 + t254 * t379 + t341 * t262 - t275 * t380;
t216 = -pkin(8) * t250 + t350;
t365 = t343 * t213 - t340 * t216;
t416 = t235 * t357 + t365;
t415 = t301 * MDP(26) + (-t238 ^ 2 + t357 ^ 2) * MDP(23) - t238 * t357 * MDP(22);
t268 = t318 * t314;
t289 = t322 * t345 + t342 * t323;
t310 = t313 * qJD(3);
t388 = t344 * t310;
t352 = -t314 * t380 - t388;
t317 = t340 * t341 - t343 * t344;
t385 = t424 * t317;
t411 = -t299 * t385 - t301 * t318;
t364 = t249 * t340 + t343 * t250;
t220 = -qJD(5) * t357 + t364;
t410 = pkin(7) + pkin(8);
t231 = t344 * t254 - t275 * t341;
t225 = -pkin(8) * t293 + t231;
t221 = pkin(4) * t304 + t225;
t407 = t221 * t343;
t406 = t226 * t343;
t405 = t238 * t309;
t404 = t357 * t309;
t403 = t249 * t341;
t402 = t291 * t304;
t401 = t291 * t309;
t400 = t293 * t304;
t399 = t293 * t309;
t394 = t314 * t341;
t393 = t314 * t344;
t390 = t341 * t301;
t389 = t341 * t310;
t290 = -t322 * t342 + t323 * t345;
t283 = t344 * t290;
t295 = t344 * t301;
t276 = pkin(3) * t309 - pkin(7) * t308;
t387 = t341 * t276 + t344 * t413;
t278 = pkin(3) * t313 - pkin(7) * t314 + t331;
t386 = t341 * t278 + t283;
t377 = qJD(5) * t343;
t374 = qJD(1) * qJD(2);
t372 = t343 * t249 - t340 * t250 - t291 * t377;
t371 = qJD(4) * t410;
t366 = t314 * t379;
t363 = t304 * t344;
t362 = qJD(5) * t221 + t216;
t243 = t419 * qJD(2) + t423;
t361 = t418 * pkin(4) - t280;
t360 = t384 * t299 - t317 * t301;
t266 = t344 * t276;
t325 = t410 * t344;
t359 = pkin(4) * t309 + qJD(5) * t325 - t413 * t341 + t266 + (-pkin(8) * t308 + t371) * t344;
t324 = t410 * t341;
t358 = -pkin(8) * t395 + qJD(5) * t324 + t341 * t371 + t387;
t215 = t221 * t340 + t406;
t354 = -t418 * t304 + t295;
t353 = t366 - t389;
t351 = -pkin(7) * t301 + t274 * t304;
t255 = -t289 * qJD(3) - t348;
t277 = pkin(3) * t311 + pkin(7) * t310;
t349 = t344 * t255 + t341 * t277 + t278 * t379 - t290 * t380;
t219 = -t293 * t378 + t372;
t256 = qJD(2) * t314 + qJD(3) * t290;
t334 = -pkin(4) * t344 - pkin(3);
t281 = t301 * t313;
t271 = t344 * t278;
t269 = t317 * t314;
t267 = t344 * t277;
t257 = pkin(4) * t394 + t289;
t234 = pkin(4) * t353 + t256;
t233 = -pkin(8) * t394 + t386;
t229 = pkin(4) * t313 - pkin(8) * t393 - t290 * t341 + t271;
t228 = pkin(4) * t250 + t243;
t223 = -t378 * t394 + (t412 * t393 - t389) * t343 + t352 * t340;
t222 = -t412 * t268 + t317 * t310;
t218 = -pkin(8) * t353 + t349;
t217 = pkin(8) * t388 + pkin(4) * t311 - t255 * t341 + t267 + (-t283 + (pkin(8) * t314 - t278) * t341) * qJD(4);
t214 = -t226 * t340 + t407;
t1 = [(t300 * t314 - t309 * t310) * MDP(8) + (-t300 * t313 - t301 * t314 - t308 * t310 - t309 * t311) * MDP(9) + (t301 * t331 + t311 * t321) * MDP(13) + (t300 * t331 - t310 * t321) * MDP(14) + (t249 * t393 + t293 * t352) * MDP(15) + (-(-t291 * t344 - t293 * t341) * t310 + (-t403 - t250 * t344 + (t291 * t341 - t293 * t344) * qJD(4)) * t314) * MDP(16) + (t249 * t313 + t293 * t311 + t295 * t314 + t304 * t352) * MDP(17) + (-t250 * t313 - t291 * t311 - t304 * t353 - t314 * t390) * MDP(18) + (t304 * t311 + t281) * MDP(19) + ((-t290 * t379 + t267) * t304 + t271 * t301 + (-t275 * t379 + t259) * t313 + t231 * t311 + t256 * t291 + t289 * t250 + t274 * t366 + ((-qJD(4) * t278 - t255) * t304 - t290 * t301 + (-qJD(4) * t254 - t242) * t313 + t243 * t314 - t274 * t310) * t341) * MDP(20) + (-t232 * t311 + t243 * t393 + t289 * t249 + t256 * t293 + t274 * t352 - t301 * t386 - t304 * t349 - t350 * t313) * MDP(21) + (-t219 * t269 - t222 * t357) * MDP(22) + (-t219 * t268 + t220 * t269 - t222 * t238 + t223 * t357) * MDP(23) + (t219 * t313 + t222 * t299 - t269 * t301 - t311 * t357) * MDP(24) + (-t220 * t313 - t223 * t299 - t238 * t311 - t268 * t301) * MDP(25) + (t299 * t311 + t281) * MDP(26) + ((t217 * t343 - t218 * t340) * t299 + (t229 * t343 - t233 * t340) * t301 + t365 * t313 + t214 * t311 + t234 * t238 + t257 * t220 + t228 * t268 + t235 * t223 + ((-t229 * t340 - t233 * t343) * t299 - t215 * t313) * qJD(5)) * MDP(27) + (-t215 * t311 + t257 * t219 + t235 * t222 + t224 * t313 - t228 * t269 - t234 * t357 + (-(-qJD(5) * t233 + t217) * t299 - t229 * t301 - t213 * t313) * t340 + (-(qJD(5) * t229 + t218) * t299 - t233 * t301 - t362 * t313) * t343) * MDP(28) + 0.2e1 * t374 * t420 + (-t310 * MDP(10) - t311 * MDP(11) - t256 * MDP(13) - t255 * MDP(14)) * qJD(3); t328 * MDP(14) + (t354 - t401) * MDP(20) + (-t304 ^ 2 * t344 - t390 - t399) * MDP(21) + (t360 - t405) * MDP(27) + (t404 + t411) * MDP(28) + ((t309 + t419) * MDP(13) + (t308 - t370) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t420; -t308 ^ 2 * MDP(9) + (t328 + (-t308 - t370) * qJD(3)) * MDP(10) + (-t243 + t423) * MDP(13) + (-t308 * t321 + t313 * t374) * MDP(14) + (t293 * t363 + t403) * MDP(15) + ((t249 - t402) * t344 + (-t250 - t400) * t341) * MDP(16) + (t304 * t363 + t390 - t399) * MDP(17) + (t354 + t401) * MDP(18) + (-pkin(3) * t250 - t243 * t344 - t280 * t291 + (-pkin(7) * t379 - t266) * t304 + (t304 * t413 + t351) * t341) * MDP(20) + (-pkin(3) * t249 + t243 * t341 - t280 * t293 + (pkin(7) * t380 + t387) * t304 + t351 * t344) * MDP(21) + (t219 * t318 - t357 * t385) * MDP(22) + (-t219 * t317 - t220 * t318 - t238 * t385 - t357 * t384) * MDP(23) + (t404 - t411) * MDP(24) + (t360 + t405) * MDP(25) + ((-t324 * t343 - t325 * t340) * t301 + t334 * t220 + t228 * t317 + (t340 * t358 - t343 * t359) * t299 + t361 * t238 - t384 * t235) * MDP(27) + (-(-t324 * t340 + t325 * t343) * t301 + t334 * t219 + t228 * t318 + (t340 * t359 + t343 * t358) * t299 - t361 * t357 + t385 * t235) * MDP(28) + (-t321 * MDP(13) - t304 * MDP(19) - t231 * MDP(20) + t232 * MDP(21) - t299 * MDP(26) - t214 * MDP(27) + t215 * MDP(28) - MDP(8) * t308 + t309 * MDP(9)) * t309; t293 * t291 * MDP(15) + (-t291 ^ 2 + t293 ^ 2) * MDP(16) + (t249 + t402) * MDP(17) + (-t250 + t400) * MDP(18) + t301 * MDP(19) + (t232 * t304 - t274 * t293 + t347) * MDP(20) + (t231 * t304 + t274 * t291 - t350) * MDP(21) + (t219 + t422) * MDP(24) + (-t220 - t421) * MDP(25) + (-(-t225 * t340 - t406) * t299 - t215 * qJD(5) + (-t238 * t293 - t299 * t378 + t343 * t301) * pkin(4) + t416) * MDP(27) + ((-t226 * t299 - t213) * t340 + (t225 * t299 - t362) * t343 + (t293 * t357 - t299 * t377 - t340 * t301) * pkin(4) + t417) * MDP(28) + t415; (t372 + t422) * MDP(24) + (-t364 - t421) * MDP(25) + (t215 * t299 + t416) * MDP(27) + (-t340 * t213 + t214 * t299 - t343 * t216 + t417) * MDP(28) + (-MDP(24) * t398 + MDP(25) * t357 - MDP(27) * t215 - MDP(28) * t407) * qJD(5) + t415;];
tauc = t1;
