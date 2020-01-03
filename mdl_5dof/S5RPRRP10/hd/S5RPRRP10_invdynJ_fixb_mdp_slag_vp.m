% Calculate vector of inverse dynamics joint torques for
% S5RPRRP10
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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:12
% EndTime: 2019-12-31 18:52:17
% DurationCPUTime: 3.32s
% Computational Cost: add. (2506->334), mult. (5912->428), div. (0->0), fcn. (4337->10), ass. (0->150)
t316 = sin(pkin(8));
t364 = qJD(1) * qJD(2);
t397 = pkin(6) + qJ(2);
t406 = t397 * qJDD(1) + t364;
t269 = t406 * t316;
t317 = cos(pkin(8));
t270 = t406 * t317;
t321 = sin(qJ(3));
t324 = cos(qJ(3));
t293 = t397 * t316;
t289 = qJD(1) * t293;
t294 = t397 * t317;
t290 = qJD(1) * t294;
t258 = -t321 * t289 + t324 * t290;
t420 = qJD(3) * t258;
t339 = -t269 * t324 - t321 * t270 - t420;
t227 = -qJDD(3) * pkin(3) - t339;
t315 = pkin(8) + qJ(3);
t307 = sin(t315);
t308 = cos(t315);
t322 = sin(qJ(1));
t325 = cos(qJ(1));
t350 = g(1) * t325 + g(2) * t322;
t330 = -g(3) * t308 + t350 * t307;
t369 = qJD(1) * t324;
t301 = t317 * t369;
t370 = qJD(1) * t321;
t358 = t316 * t370;
t283 = t301 - t358;
t273 = qJD(4) - t283;
t368 = qJD(4) * t273;
t421 = -pkin(7) * t368 - t227 + t330;
t288 = t316 * t324 + t317 * t321;
t284 = t288 * qJD(1);
t303 = pkin(2) * t317 + pkin(1);
t292 = -t303 * qJD(1) + qJD(2);
t236 = -pkin(3) * t283 - pkin(7) * t284 + t292;
t252 = qJD(3) * pkin(7) + t258;
t320 = sin(qJ(4));
t323 = cos(qJ(4));
t221 = t236 * t320 + t252 * t323;
t286 = t288 * qJD(3);
t362 = qJDD(1) * t324;
t300 = t317 * t362;
t363 = qJDD(1) * t321;
t347 = -t316 * t363 + t300;
t256 = qJD(1) * t286 - t347;
t291 = -t303 * qJDD(1) + qJDD(2);
t359 = qJD(3) * t301 + t316 * t362 + t317 * t363;
t333 = qJD(3) * t358 - t359;
t228 = t256 * pkin(3) + t333 * pkin(7) + t291;
t225 = t323 * t228;
t346 = -t269 * t321 + t270 * t324;
t411 = -t289 * t324 - t321 * t290;
t226 = qJDD(3) * pkin(7) + qJD(3) * t411 + t346;
t365 = t323 * qJD(3);
t367 = qJD(4) * t320;
t230 = -qJD(4) * t365 - t320 * qJDD(3) + t284 * t367 + t323 * t333;
t250 = qJDD(4) + t256;
t266 = qJD(3) * t320 + t284 * t323;
t209 = pkin(4) * t250 + qJ(5) * t230 - t221 * qJD(4) - qJD(5) * t266 - t226 * t320 + t225;
t264 = t284 * t320 - t365;
t216 = -qJ(5) * t264 + t221;
t419 = t273 * t216 + t209;
t231 = qJD(4) * t266 - t323 * qJDD(3) - t320 * t333;
t366 = qJD(4) * t323;
t334 = t323 * t226 + t320 * t228 + t236 * t366 - t252 * t367;
t210 = -qJ(5) * t231 - qJD(5) * t264 + t334;
t220 = t323 * t236 - t252 * t320;
t215 = -qJ(5) * t266 + t220;
t214 = pkin(4) * t273 + t215;
t418 = -t273 * t214 + t210;
t396 = qJ(5) + pkin(7);
t417 = qJ(5) * t283 - qJD(4) * t396;
t410 = g(1) * t322 - g(2) * t325;
t416 = qJDD(2) - t410;
t413 = t410 * t307;
t379 = t323 * t325;
t382 = t320 * t322;
t274 = t308 * t382 + t379;
t380 = t322 * t323;
t381 = t320 * t325;
t276 = -t308 * t381 + t380;
t409 = -g(1) * t276 + g(2) * t274;
t408 = qJ(2) * qJDD(1);
t399 = g(3) * t307;
t407 = t350 * t308 + t399;
t405 = t266 ^ 2;
t394 = qJDD(1) * pkin(1);
t393 = t230 * t320;
t392 = t264 * t283;
t391 = t264 * t284;
t390 = t266 * t273;
t389 = t266 * t284;
t388 = t288 * t320;
t387 = t288 * t323;
t385 = t317 * MDP(4);
t384 = t320 * t250;
t383 = t320 * t273;
t240 = t323 * t250;
t262 = -t293 * t321 + t294 * t324;
t259 = t323 * t262;
t378 = -t215 + t214;
t377 = -t320 * t231 - t264 * t366;
t376 = t283 * t383 + t240;
t253 = pkin(3) * t284 - pkin(7) * t283;
t375 = t320 * t253 + t323 * t411;
t287 = t316 * t321 - t324 * t317;
t255 = pkin(3) * t287 - pkin(7) * t288 - t303;
t374 = t320 * t255 + t259;
t373 = qJD(5) * t323 + t417 * t320 - t375;
t243 = t323 * t253;
t372 = -pkin(4) * t284 - t243 + t417 * t323 + (-qJD(5) + t411) * t320;
t371 = t316 ^ 2 + t317 ^ 2;
t344 = -t293 * t324 - t294 * t321;
t237 = -t287 * qJD(2) + t344 * qJD(3);
t285 = t287 * qJD(3);
t254 = pkin(3) * t286 + pkin(7) * t285;
t360 = t323 * t237 + t320 * t254 + t255 * t366;
t357 = t288 * t366;
t356 = pkin(4) * t320 + t397;
t353 = t273 * t323;
t352 = -qJD(4) * t236 - t226;
t351 = 0.2e1 * t371;
t348 = -t252 * t366 + t225;
t305 = pkin(4) * t323 + pkin(3);
t343 = t305 * t308 + t307 * t396;
t341 = qJ(5) * t285 - qJD(5) * t288;
t340 = t394 - t416;
t251 = -qJD(3) * pkin(3) - t411;
t338 = t303 + t343;
t337 = -t285 * t320 + t357;
t336 = -t285 * t323 - t288 * t367;
t335 = -pkin(7) * t250 + t273 * t251;
t328 = t351 * t364 - t350;
t213 = pkin(4) * t231 + qJDD(5) + t227;
t238 = t288 * qJD(2) + t262 * qJD(3);
t296 = t396 * t323;
t295 = t396 * t320;
t277 = t308 * t379 + t382;
t275 = -t308 * t380 + t381;
t263 = t264 ^ 2;
t246 = t323 * t255;
t244 = t323 * t254;
t232 = pkin(4) * t264 + qJD(5) + t251;
t222 = -qJ(5) * t388 + t374;
t218 = pkin(4) * t287 - qJ(5) * t387 - t262 * t320 + t246;
t212 = -qJ(5) * t357 + (-qJD(4) * t262 + t341) * t320 + t360;
t211 = pkin(4) * t286 - t237 * t320 + t244 + t341 * t323 + (-t259 + (qJ(5) * t288 - t255) * t320) * qJD(4);
t1 = [qJDD(1) * MDP(1) + t410 * MDP(2) + t350 * MDP(3) + (t351 * t408 + t328) * MDP(6) + (t340 * pkin(1) + (t371 * t408 + t328) * qJ(2)) * MDP(7) + (-t284 * t285 - t333 * t288) * MDP(8) + (-t288 * t256 - t285 * t283 - t284 * t286 + t333 * t287) * MDP(9) + (-qJD(3) * t285 + qJDD(3) * t288) * MDP(10) + (-qJD(3) * t286 - qJDD(3) * t287) * MDP(11) + (-qJD(3) * t238 + qJDD(3) * t344 - t256 * t303 + t286 * t292 + t287 * t291 + t308 * t410) * MDP(13) + (-t237 * qJD(3) - t262 * qJDD(3) - t292 * t285 + t291 * t288 + t303 * t333 - t413) * MDP(14) + (-t230 * t387 + t336 * t266) * MDP(15) + (-(-t264 * t323 - t266 * t320) * t285 + (t393 - t231 * t323 + (t264 * t320 - t266 * t323) * qJD(4)) * t288) * MDP(16) + (-t230 * t287 + t288 * t240 + t266 * t286 + t336 * t273) * MDP(17) + (-t231 * t287 - t264 * t286 - t337 * t273 - t288 * t384) * MDP(18) + (t250 * t287 + t273 * t286) * MDP(19) + ((-t262 * t366 + t244) * t273 + t246 * t250 + t348 * t287 + t220 * t286 + t238 * t264 - t344 * t231 + t251 * t357 - g(1) * t275 - g(2) * t277 + ((-qJD(4) * t255 - t237) * t273 - t262 * t250 + t352 * t287 + t227 * t288 - t251 * t285) * t320) * MDP(20) + (-(-t262 * t367 + t360) * t273 - t374 * t250 - t334 * t287 - t221 * t286 + t238 * t266 + t344 * t230 + t227 * t387 - g(1) * t274 - g(2) * t276 + t336 * t251) * MDP(21) + (-t211 * t266 - t212 * t264 + t218 * t230 - t222 * t231 + t413 - (-t214 * t323 - t216 * t320) * t285 + (-t209 * t323 - t210 * t320 + (t214 * t320 - t216 * t323) * qJD(4)) * t288) * MDP(22) + (t210 * t222 + t216 * t212 + t209 * t218 + t214 * t211 + t213 * (pkin(4) * t388 - t344) + t232 * (t337 * pkin(4) + t238) + (-g(1) * t356 - g(2) * t338) * t325 + (g(1) * t338 - g(2) * t356) * t322) * MDP(23) + (-MDP(5) * t316 + t385) * (t340 + t394); t416 * MDP(7) - t300 * MDP(13) + t359 * MDP(14) + (t376 - t391) * MDP(20) - MDP(21) * t389 + t377 * MDP(22) + (-t232 * t284 - t410) * MDP(23) + (-t385 - pkin(1) * MDP(7) + (MDP(13) * t321 + MDP(5)) * t316) * qJDD(1) + (-MDP(20) * t368 - t250 * MDP(21) + MDP(22) * t390 + t418 * MDP(23)) * t320 + ((t316 * t369 + t317 * t370 + t284) * MDP(13) + (t283 - t358) * MDP(14)) * qJD(3) + ((t230 + t392) * MDP(22) + t419 * MDP(23) - t273 ^ 2 * MDP(21)) * t323 + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t371; -t283 ^ 2 * MDP(9) + ((-t283 - t358) * qJD(3) + t359) * MDP(10) + t347 * MDP(11) + qJDD(3) * MDP(12) + (t330 + t339 + t420) * MDP(13) + (-t283 * t292 - t346 + t407) * MDP(14) + (t266 * t353 - t393) * MDP(15) + ((-t230 + t392) * t323 - t266 * t383 + t377) * MDP(16) + (t273 * t353 + t384 - t389) * MDP(17) + (-t273 * t367 + t376 + t391) * MDP(18) + (-pkin(3) * t231 - t243 * t273 - t258 * t264 + (t273 * t411 + t335) * t320 + t421 * t323) * MDP(20) + (pkin(3) * t230 - t258 * t266 + t375 * t273 - t421 * t320 + t335 * t323) * MDP(21) + (-t230 * t295 - t231 * t296 - t373 * t264 - t372 * t266 - t419 * t320 + t418 * t323 - t407) * MDP(22) + (t210 * t296 - t209 * t295 - t213 * t305 - g(3) * t343 + (pkin(4) * t383 - t258) * t232 + t373 * t216 + t372 * t214 + t350 * (t305 * t307 - t308 * t396)) * MDP(23) + (-t292 * MDP(13) - t273 * MDP(19) - t220 * MDP(20) + t221 * MDP(21) - t283 * MDP(8) + MDP(9) * t284) * t284; t266 * t264 * MDP(15) + (-t263 + t405) * MDP(16) + (t264 * t273 - t230) * MDP(17) + (-t231 + t390) * MDP(18) + t250 * MDP(19) + (t221 * t273 - t251 * t266 + (t352 + t399) * t320 + t348 + t409) * MDP(20) + (g(1) * t277 - g(2) * t275 + t220 * t273 + t251 * t264 + t323 * t399 - t334) * MDP(21) + (pkin(4) * t230 - t378 * t264) * MDP(22) + (t378 * t216 + (-t232 * t266 + t320 * t399 + t209 + t409) * pkin(4)) * MDP(23); (-t263 - t405) * MDP(22) + (t214 * t266 + t216 * t264 + t213 - t330) * MDP(23);];
tau = t1;
