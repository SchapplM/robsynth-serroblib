% Calculate vector of inverse dynamics joint torques for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR14_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:12
% EndTime: 2021-01-15 12:17:21
% DurationCPUTime: 3.48s
% Computational Cost: add. (1782->328), mult. (3498->425), div. (0->0), fcn. (2320->10), ass. (0->145)
t304 = cos(qJ(3));
t374 = cos(pkin(8));
t336 = t374 * t304;
t325 = qJD(1) * t336;
t299 = sin(pkin(8));
t301 = sin(qJ(3));
t354 = qJD(1) * t301;
t339 = t299 * t354;
t255 = t325 - t339;
t300 = sin(qJ(5));
t303 = cos(qJ(5));
t234 = -t303 * qJD(3) + t255 * t300;
t319 = -t299 * t304 - t301 * t374;
t382 = t319 * qJD(1);
t386 = qJD(5) - t382;
t395 = t234 * t386;
t236 = qJD(3) * t300 + t255 * t303;
t394 = t236 * t386;
t302 = sin(qJ(1));
t305 = cos(qJ(1));
t383 = g(1) * t302 - g(2) * t305;
t393 = qJDD(2) - t383;
t333 = qJDD(1) * t374;
t344 = qJDD(1) * t301;
t322 = -t299 * t344 + t304 * t333;
t387 = qJD(3) * t382;
t231 = -t322 - t387;
t392 = -qJD(3) * qJD(5) + t231;
t306 = -pkin(1) - pkin(6);
t265 = qJDD(1) * t306 + qJDD(2);
t260 = t304 * t265;
t266 = qJD(1) * t306 + qJD(2);
t347 = qJD(1) * qJD(3);
t338 = t301 * t347;
t343 = qJDD(1) * t304;
t346 = qJD(1) * qJD(4);
t352 = qJD(3) * t301;
t214 = -t304 * t346 - t266 * t352 + qJDD(3) * pkin(3) + t260 + (t338 - t343) * qJ(4);
t331 = -qJ(4) * qJD(1) + t266;
t351 = qJD(3) * t304;
t220 = t331 * t351 + (-qJ(4) * qJDD(1) + t265 - t346) * t301;
t203 = t374 * t214 - t299 * t220;
t199 = -qJDD(3) * pkin(4) - t203;
t277 = pkin(3) * t299 + pkin(7);
t293 = qJ(3) + pkin(8);
t284 = sin(t293);
t285 = cos(t293);
t312 = g(3) * t284 - t285 * t383;
t353 = qJD(1) * t304;
t391 = t312 - t386 * (pkin(3) * t353 + pkin(4) * t255 - pkin(7) * t382 + qJD(5) * t277) - t199;
t294 = qJDD(1) * qJ(2);
t295 = qJD(1) * qJD(2);
t324 = g(1) * t305 + g(2) * t302;
t316 = -t324 + 0.2e1 * t295;
t390 = 0.2e1 * t294 + t316;
t372 = qJDD(1) * pkin(1);
t389 = t372 - t393;
t332 = t303 * t386;
t341 = qJD(3) * t325 + t299 * t343 + t301 * t333;
t230 = t299 * t338 - t341;
t228 = -qJDD(5) + t230;
t364 = t300 * t228;
t388 = t332 * t386 - t364;
t384 = qJ(4) - t306;
t245 = -qJ(4) * t353 + t304 * t266;
t242 = qJD(3) * pkin(3) + t245;
t244 = t331 * t301;
t365 = t299 * t244;
t215 = t242 * t374 - t365;
t211 = -qJD(3) * pkin(4) - t215;
t334 = t384 * t304;
t243 = -qJD(3) * t334 - qJD(4) * t301;
t314 = -qJD(4) * t304 + t352 * t384;
t219 = t243 * t374 + t299 * t314;
t258 = -t299 * t301 + t336;
t278 = pkin(3) * t301 + qJ(2);
t229 = -pkin(4) * t319 - pkin(7) * t258 + t278;
t263 = t384 * t301;
t233 = -t263 * t374 - t299 * t334;
t335 = qJD(3) * t374;
t254 = -t299 * t351 - t301 * t335;
t204 = t299 * t214 + t374 * t220;
t200 = qJDD(3) * pkin(7) + t204;
t262 = pkin(3) * t354 + qJD(1) * qJ(2) + qJD(4);
t217 = -pkin(4) * t382 - pkin(7) * t255 + t262;
t330 = qJD(5) * t217 + t200;
t380 = t199 * t258 + t211 * t254 + t233 * t228 - (qJD(5) * t229 + t219) * t386 + t319 * t330;
t376 = g(3) * t285;
t375 = g(3) * t301;
t308 = qJD(1) ^ 2;
t373 = qJ(2) * t308;
t340 = t300 * qJDD(3) - t392 * t303;
t350 = qJD(5) * t300;
t207 = -t255 * t350 + t340;
t371 = t207 * t258;
t370 = t207 * t300;
t369 = t211 * t258;
t368 = t229 * t228;
t367 = t234 * t255;
t366 = t236 * t255;
t363 = t300 * t302;
t362 = t300 * t305;
t361 = t302 * t303;
t225 = t303 * t228;
t360 = t303 * t305;
t240 = t374 * t244;
t216 = t299 * t242 + t240;
t298 = t304 ^ 2;
t357 = t301 ^ 2 - t298;
t307 = qJD(3) ^ 2;
t356 = -t307 - t308;
t355 = qJD(1) * t262;
t349 = qJD(5) * t303;
t270 = pkin(3) * t351 + qJD(2);
t342 = qJDD(3) * t301;
t337 = t304 * t347;
t239 = qJDD(4) + t294 + t295 + (t337 + t344) * pkin(3);
t206 = -pkin(4) * t230 + pkin(7) * t231 + t239;
t212 = qJD(3) * pkin(7) + t216;
t329 = qJD(5) * t212 - t206;
t321 = -t225 + (t300 * t382 - t350) * t386;
t320 = t254 * t303 - t258 * t350;
t318 = 0.2e1 * qJ(2) * t347 + qJDD(3) * t306;
t317 = -t383 - t373;
t222 = t245 * t374 - t365;
t311 = t277 * t228 + (t211 + t222) * t386;
t253 = t299 * t352 - t304 * t335;
t310 = t203 * t258 - t204 * t319 + t215 * t254 - t216 * t253 - t383;
t309 = -t306 * t307 + t390;
t288 = qJDD(3) * t304;
t287 = t303 * qJDD(3);
t279 = -pkin(3) * t374 - pkin(4);
t250 = t284 * t360 - t363;
t249 = t284 * t362 + t361;
t248 = t284 * t361 + t362;
t247 = -t284 * t363 + t360;
t232 = -t263 * t299 + t334 * t374;
t223 = -pkin(4) * t253 - pkin(7) * t254 + t270;
t221 = t245 * t299 + t240;
t218 = t243 * t299 - t314 * t374;
t208 = qJD(5) * t236 - t231 * t300 - t287;
t205 = t303 * t206;
t202 = t212 * t303 + t217 * t300;
t201 = -t212 * t300 + t217 * t303;
t1 = [qJDD(1) * MDP(1) + t383 * MDP(2) + t324 * MDP(3) + (-0.2e1 * t372 + t393) * MDP(4) + t390 * MDP(5) + (t389 * pkin(1) + (t316 + t294) * qJ(2)) * MDP(6) + (qJDD(1) * t298 - 0.2e1 * t301 * t337) * MDP(7) + 0.2e1 * (-t301 * t343 + t347 * t357) * MDP(8) + (-t301 * t307 + t288) * MDP(9) + (-t304 * t307 - t342) * MDP(10) + (t301 * t309 + t304 * t318) * MDP(12) + (-t301 * t318 + t304 * t309) * MDP(13) + (-qJD(3) * t218 - qJDD(3) * t232 - t230 * t278 - t239 * t319 - t253 * t262 - t270 * t382 - t284 * t324) * MDP(14) + (-qJD(3) * t219 - qJDD(3) * t233 - t231 * t278 + t239 * t258 + t254 * t262 + t255 * t270 - t285 * t324) * MDP(15) + (t218 * t255 + t219 * t382 + t230 * t233 - t231 * t232 - t310) * MDP(16) + (t204 * t233 + t216 * t219 - t203 * t232 - t215 * t218 + t239 * t278 + t262 * t270 - g(1) * (t278 * t305 - t302 * t384) - g(2) * (t278 * t302 + t305 * t384)) * MDP(17) + (t236 * t320 + t303 * t371) * MDP(18) + ((-t234 * t303 - t236 * t300) * t254 + (-t370 - t208 * t303 + (t234 * t300 - t236 * t303) * qJD(5)) * t258) * MDP(19) + (-t207 * t319 - t225 * t258 - t236 * t253 + t320 * t386) * MDP(20) + (t258 * t364 + t208 * t319 + t234 * t253 + (-t254 * t300 - t258 * t349) * t386) * MDP(21) + (t228 * t319 - t253 * t386) * MDP(22) + (-g(1) * t250 - g(2) * t248 - t201 * t253 - t205 * t319 + t232 * t208 + t218 * t234 + (t223 * t386 - t368 + (t212 * t319 - t233 * t386 + t369) * qJD(5)) * t303 + t380 * t300) * MDP(23) + (g(1) * t249 - g(2) * t247 + t202 * t253 + t232 * t207 + t218 * t236 + (-(-qJD(5) * t233 + t223) * t386 + t368 - t329 * t319 - qJD(5) * t369) * t300 + t380 * t303) * MDP(24); qJDD(1) * MDP(4) - t308 * MDP(5) + (-t373 - t389) * MDP(6) + (t301 * t356 + t288) * MDP(12) + (t304 * t356 - t342) * MDP(13) + (qJD(1) * t382 + qJD(3) * t254 + qJDD(3) * t258) * MDP(14) + (-qJD(1) * t255 + qJD(3) * t253 + qJDD(3) * t319) * MDP(15) + (-t230 * t319 + t231 * t258 - t253 * t382 - t254 * t255) * MDP(16) + (t310 - t355) * MDP(17) + (-t208 * t258 - t234 * t254 - t319 * t364) * MDP(23) + (-t225 * t319 - t236 * t254 - t371) * MDP(24) + ((-qJD(1) * t303 + t253 * t300 + t319 * t349) * MDP(23) + (qJD(1) * t300 + t253 * t303 - t319 * t350) * MDP(24)) * t386; MDP(9) * t343 - MDP(10) * t344 + qJDD(3) * MDP(11) + (t304 * t317 + t260 + t375) * MDP(12) + (g(3) * t304 + (-t265 - t317) * t301) * MDP(13) + (t221 * qJD(3) - t262 * t255 + (qJDD(3) * t374 + t353 * t382) * pkin(3) + t312 + t203) * MDP(14) + (t376 + qJD(3) * t222 - t382 * t262 + t383 * t284 + (-qJDD(3) * t299 - t255 * t353) * pkin(3) - t204) * MDP(15) + ((t216 - t221) * t255 - (-t215 + t222) * t382 + (t230 * t299 + t231 * t374) * pkin(3)) * MDP(16) + (t215 * t221 - t216 * t222 + (t374 * t203 + t375 + t204 * t299 + (-t383 - t355) * t304) * pkin(3)) * MDP(17) + (t236 * t332 + t370) * MDP(18) + ((t207 - t395) * t303 + (-t208 - t394) * t300) * MDP(19) + (-t366 + t388) * MDP(20) + (t321 + t367) * MDP(21) - t386 * t255 * MDP(22) + (-t201 * t255 + t279 * t208 - t221 * t234 + t311 * t300 + t303 * t391) * MDP(23) + (t202 * t255 + t279 * t207 - t221 * t236 - t300 * t391 + t311 * t303) * MDP(24) + (MDP(7) * t301 * t304 - MDP(8) * t357) * t308; ((t255 - t339) * qJD(3) + t341) * MDP(14) + (t322 + 0.2e1 * t387) * MDP(15) + (-t255 ^ 2 - t382 ^ 2) * MDP(16) + (t215 * t255 - t216 * t382 + t239 - t324) * MDP(17) + (t321 - t367) * MDP(23) + (-t366 - t388) * MDP(24); t236 * t234 * MDP(18) + (-t234 ^ 2 + t236 ^ 2) * MDP(19) + (t340 + t395) * MDP(20) + (t287 + t394) * MDP(21) - t228 * MDP(22) + (-g(1) * t247 - g(2) * t249 + t202 * t386 - t211 * t236 + t205) * MDP(23) + (g(1) * t248 - g(2) * t250 + t201 * t386 + t211 * t234) * MDP(24) + ((-t200 + t376) * MDP(24) + (-MDP(21) * t255 - MDP(23) * t212 - MDP(24) * t217) * qJD(5)) * t303 + (-qJD(5) * t255 * MDP(20) + t392 * MDP(21) + (-t330 + t376) * MDP(23) + t329 * MDP(24)) * t300;];
tau = t1;
