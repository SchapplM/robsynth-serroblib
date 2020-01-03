% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:43
% EndTime: 2019-12-31 20:30:49
% DurationCPUTime: 2.86s
% Computational Cost: add. (1709->287), mult. (3947->396), div. (0->0), fcn. (2566->6), ass. (0->136)
t303 = sin(qJ(4));
t306 = cos(qJ(4));
t307 = cos(qJ(2));
t356 = qJD(1) * t307;
t304 = sin(qJ(2));
t357 = qJD(1) * t304;
t253 = -t303 * t356 + t306 * t357;
t302 = sin(qJ(5));
t351 = qJD(5) * t302;
t258 = t303 * t304 + t306 * t307;
t320 = t258 * qJD(4);
t348 = qJD(1) * qJD(2);
t342 = t307 * t348;
t343 = t304 * t348;
t362 = t303 * t343 + t306 * t342;
t222 = -qJD(1) * t320 + t362;
t296 = qJD(2) - qJD(4);
t305 = cos(qJ(5));
t350 = qJD(5) * t305;
t365 = t305 * t222 - t296 * t350;
t210 = -t253 * t351 + t365;
t232 = t253 * t305 - t296 * t302;
t374 = t222 * t302;
t211 = t232 * qJD(5) + t374;
t352 = qJD(4) * t306;
t353 = qJD(4) * t303;
t354 = qJD(2) * t307;
t404 = t303 * t354 + t304 * t352 - t307 * t353;
t223 = t404 * qJD(1) - t306 * t343;
t371 = t253 * t302;
t230 = t305 * t296 + t371;
t368 = t305 * t223;
t370 = t302 * t223;
t377 = t210 * t302;
t387 = t258 * qJD(1);
t398 = qJD(5) + t387;
t389 = t398 * t232;
t390 = t398 ^ 2;
t405 = t230 * t398;
t406 = -((t211 + t389) * t302 + (-t210 + t405) * t305) * MDP(23) + (t305 * t389 + t377) * MDP(22) + (-t232 * t253 + t305 * t390 + t370) * MDP(24) - (-t230 * t253 + t302 * t390 - t368) * MDP(25) - (t253 * t296 + t223) * MDP(18) - (-t253 ^ 2 + t387 ^ 2) * MDP(16) + (MDP(15) * t387 - MDP(26) * t398) * t253;
t288 = pkin(6) * t357;
t399 = -pkin(7) * t357 + qJD(3) + t288;
t393 = -0.2e1 * t348;
t299 = t304 ^ 2;
t392 = MDP(5) * (-t307 ^ 2 + t299);
t308 = -pkin(2) - pkin(3);
t346 = t308 * qJD(2);
t243 = t346 + t399;
t289 = pkin(6) * t356;
t265 = -pkin(7) * t356 + t289;
t298 = qJD(2) * qJ(3);
t254 = t265 + t298;
t219 = t243 * t306 - t254 * t303;
t216 = pkin(4) * t296 - t219;
t391 = t216 * t398;
t355 = qJD(2) * t304;
t381 = pkin(6) - pkin(7);
t264 = t381 * t355;
t297 = qJD(2) * qJD(3);
t246 = -qJD(1) * t264 + t297;
t282 = pkin(6) * t342;
t256 = -pkin(7) * t342 + t282;
t209 = t243 * t353 + t303 * t246 + t254 * t352 - t256 * t306;
t224 = pkin(4) * t253 + pkin(8) * t387;
t286 = qJ(3) * t356;
t248 = t308 * t357 + t286;
t327 = qJ(3) * t306 + t303 * t308;
t262 = -pkin(8) + t327;
t386 = t398 * (qJD(5) * t262 - t224 + t248) - t209;
t385 = t209 + (pkin(8) * qJD(5) + t224) * t398;
t208 = t243 * t352 + t306 * t246 - t254 * t353 + t303 * t256;
t272 = t381 * t307;
t266 = qJD(2) * t272;
t271 = t381 * t304;
t323 = t271 * t306 - t272 * t303;
t212 = qJD(4) * t323 - t264 * t306 + t266 * t303;
t255 = -qJD(1) * pkin(1) - pkin(2) * t356 - qJ(3) * t357;
t241 = pkin(3) * t356 - t255;
t214 = pkin(4) * t387 - pkin(8) * t253 + t241;
t269 = -t307 * pkin(2) - t304 * qJ(3) - pkin(1);
t257 = t307 * pkin(3) - t269;
t259 = -t303 * t307 + t304 * t306;
t218 = pkin(4) * t258 - pkin(8) * t259 + t257;
t228 = qJD(2) * t258 - t320;
t234 = t271 * t303 + t272 * t306;
t382 = -(qJD(5) * t218 + t212) * t398 - t234 * t223 - (qJD(5) * t214 + t208) * t258 + t209 * t259 + t216 * t228;
t380 = qJD(2) * pkin(2);
t220 = t243 * t303 + t254 * t306;
t217 = -pkin(8) * t296 + t220;
t203 = t214 * t305 - t217 * t302;
t379 = t203 * t253;
t204 = t214 * t302 + t217 * t305;
t378 = t204 * t253;
t376 = t216 * t259;
t375 = t218 * t223;
t373 = t387 * t296;
t309 = qJD(2) ^ 2;
t369 = t304 * t309;
t367 = t307 * t309;
t310 = qJD(1) ^ 2;
t366 = t307 * t310;
t326 = -qJ(3) * t303 + t306 * t308;
t364 = -qJD(4) * t326 + t265 * t303 - t306 * t399;
t363 = qJD(4) * t327 + t265 * t306 + t303 * t399;
t292 = t304 * qJD(3);
t360 = qJ(3) * t342 + qJD(1) * t292;
t359 = qJ(3) * t354 + t292;
t347 = t304 * t366;
t340 = pkin(1) * t393;
t339 = qJD(3) - t380;
t334 = qJD(1) * t269 + t255;
t330 = t296 ^ 2;
t329 = t306 * t296;
t328 = t304 * t346;
t240 = pkin(2) * t343 - t360;
t249 = pkin(2) * t355 - t359;
t322 = -pkin(6) * t309 - qJD(1) * t249 - t240;
t321 = t228 * t305 - t259 * t351;
t237 = t328 + t359;
t229 = qJD(1) * t328 + t360;
t315 = t241 * t253 + t209;
t314 = -pkin(8) * t223 + t219 * t398 + t391;
t313 = -t241 * t387 + t208;
t312 = -t262 * t223 + t364 * t398 - t391;
t267 = -pkin(6) * t343 + t297;
t268 = t288 + t339;
t270 = t289 + t298;
t311 = t267 * t307 + (t268 * t307 + (-t270 + t289) * t304) * qJD(2);
t261 = pkin(4) - t326;
t260 = pkin(2) * t357 - t286;
t227 = -t306 * t355 + t404;
t213 = qJD(4) * t234 - t264 * t303 - t266 * t306;
t207 = pkin(4) * t227 - pkin(8) * t228 + t237;
t206 = pkin(4) * t223 - pkin(8) * t222 + t229;
t205 = t305 * t206;
t1 = [0.2e1 * t304 * MDP(4) * t342 + t392 * t393 + MDP(6) * t367 - MDP(7) * t369 + (-pkin(6) * t367 + t304 * t340) * MDP(9) + (pkin(6) * t369 + t307 * t340) * MDP(10) + (t307 * t322 + t334 * t355) * MDP(11) + t311 * MDP(12) + (t304 * t322 - t334 * t354) * MDP(13) + (pkin(6) * t311 + t240 * t269 + t249 * t255) * MDP(14) + (t222 * t259 + t228 * t253) * MDP(15) + (-t222 * t258 - t223 * t259 - t227 * t253 - t228 * t387) * MDP(16) + (t223 * t257 + t227 * t241 + t229 * t258 + t237 * t387) * MDP(20) + (t222 * t257 + t228 * t241 + t229 * t259 + t237 * t253) * MDP(21) + (t210 * t259 * t305 + t321 * t232) * MDP(22) + ((-t230 * t305 - t232 * t302) * t228 + (-t377 - t211 * t305 + (t230 * t302 - t232 * t305) * qJD(5)) * t259) * MDP(23) + (t210 * t258 + t227 * t232 + t259 * t368 + t321 * t398) * MDP(24) + (-t259 * t370 - t211 * t258 - t227 * t230 + (-t228 * t302 - t259 * t350) * t398) * MDP(25) + (t223 * t258 + t227 * t398) * MDP(26) + (t203 * t227 + t205 * t258 - t323 * t211 + t213 * t230 + (t207 * t398 + t375 + (-t217 * t258 - t234 * t398 + t376) * qJD(5)) * t305 + t382 * t302) * MDP(27) + (-t204 * t227 - t323 * t210 + t213 * t232 + (-(-qJD(5) * t234 + t207) * t398 - t375 - (-qJD(5) * t217 + t206) * t258 - qJD(5) * t376) * t302 + t382 * t305) * MDP(28) + (-MDP(17) * t228 + MDP(18) * t227 + MDP(20) * t213 + MDP(21) * t212) * t296; -MDP(4) * t347 + t310 * t392 + 0.2e1 * t297 * MDP(13) + (qJ(3) * t267 + qJD(3) * t270 - t255 * t260) * MDP(14) + (qJD(4) * t387 - t362 + t373) * MDP(17) + (-t248 * t387 + t296 * t363 + t315) * MDP(20) + (-t248 * t253 - t296 * t364 + t313) * MDP(21) + (t261 * t211 + t363 * t230 + t312 * t302 - t386 * t305 + t379) * MDP(27) + (t261 * t210 + t363 * t232 + t386 * t302 + t312 * t305 - t378) * MDP(28) + ((-t255 * t304 + t260 * t307) * MDP(11) + ((t270 - t298) * t304 + (-t268 + t339) * t307) * MDP(12) + (t255 * t307 + t260 * t304) * MDP(13) + (t270 * t304 + (-t268 - t380) * t307) * pkin(6) * MDP(14)) * qJD(1) + (MDP(9) * t304 * t310 + MDP(10) * t366) * pkin(1) - t406; -MDP(11) * t347 + (-t299 * t310 - t309) * MDP(13) + (-qJD(2) * t270 + t255 * t357 + t282) * MDP(14) + (-t303 * t330 - t357 * t387) * MDP(20) + (-t253 * t357 - t306 * t330) * MDP(21) + (-t306 * t211 + (t302 * t329 - t305 * t357) * t398 + (-t230 * t296 - t350 * t398 - t370) * t303) * MDP(27) + (-t306 * t210 + (t302 * t357 + t305 * t329) * t398 + (-t232 * t296 + t351 * t398 - t368) * t303) * MDP(28); (t222 - t373) * MDP(17) + (-t220 * t296 - t315) * MDP(20) + (-t219 * t296 - t313) * MDP(21) + (-pkin(4) * t211 - t220 * t230 + t314 * t302 - t385 * t305 - t379) * MDP(27) + (-pkin(4) * t210 - t220 * t232 + t385 * t302 + t314 * t305 + t378) * MDP(28) + t406; t232 * t230 * MDP(22) + (-t230 ^ 2 + t232 ^ 2) * MDP(23) + (t365 + t405) * MDP(24) + (-t374 + t389) * MDP(25) + t223 * MDP(26) + (t204 * t398 - t208 * t302 - t216 * t232 + t205) * MDP(27) + (t203 * t398 - t206 * t302 - t208 * t305 + t216 * t230) * MDP(28) + (-MDP(24) * t371 - MDP(25) * t232 - MDP(27) * t204 - MDP(28) * t203) * qJD(5);];
tauc = t1;
