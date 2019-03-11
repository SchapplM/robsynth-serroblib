% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:09
% EndTime: 2019-03-09 01:32:14
% DurationCPUTime: 3.16s
% Computational Cost: add. (1809->314), mult. (3395->397), div. (0->0), fcn. (2438->14), ass. (0->146)
t301 = sin(pkin(10));
t306 = sin(qJ(5));
t358 = qJD(1) * t306;
t340 = t301 * t358;
t303 = cos(pkin(10));
t309 = cos(qJ(5));
t357 = qJD(1) * t309;
t342 = t303 * t357;
t253 = -t340 + t342;
t305 = sin(qJ(6));
t308 = cos(qJ(6));
t232 = -t308 * qJD(5) + t253 * t305;
t259 = t301 * t309 + t303 * t306;
t251 = t259 * qJD(1);
t388 = qJD(6) + t251;
t394 = t232 * t388;
t234 = qJD(5) * t305 + t253 * t308;
t393 = t234 * t388;
t302 = sin(pkin(9));
t272 = pkin(1) * t302 + qJ(3);
t361 = qJDD(1) * t272;
t348 = qJDD(1) * t303;
t264 = t309 * t348;
t349 = qJDD(1) * t301;
t327 = -t306 * t349 + t264;
t222 = -qJD(5) * t251 + t327;
t392 = qJD(5) * qJD(6) + t222;
t304 = cos(pkin(9));
t276 = -pkin(1) * t304 - pkin(2);
t269 = -qJ(4) + t276;
t257 = qJD(1) * t269 + qJD(3);
t235 = -qJD(2) * t301 + t303 * t257;
t228 = -pkin(7) * qJD(1) * t303 + t235;
t236 = t303 * qJD(2) + t301 * t257;
t359 = qJD(1) * t301;
t229 = -pkin(7) * t359 + t236;
t208 = t228 * t306 + t229 * t309;
t382 = qJD(1) * qJD(4) - qJDD(1) * t269;
t244 = qJDD(3) - t382;
t230 = -qJDD(2) * t301 + t303 * t244;
t226 = -pkin(7) * t348 + t230;
t231 = t303 * qJDD(2) + t301 * t244;
t227 = -pkin(7) * t349 + t231;
t325 = -t226 * t309 + t227 * t306;
t197 = -qJDD(5) * pkin(5) + qJD(5) * t208 + t325;
t295 = pkin(10) + qJ(5);
t284 = sin(t295);
t286 = cos(t295);
t296 = qJ(1) + pkin(9);
t285 = sin(t296);
t287 = cos(t296);
t386 = g(1) * t285 - g(2) * t287;
t315 = g(3) * t284 - t286 * t386;
t391 = t315 - t388 * (pkin(5) * t253 + t388 * pkin(8)) - t197;
t336 = t308 * t388;
t355 = qJD(5) * t309;
t341 = t303 * t355;
t263 = qJD(5) * t340;
t383 = -t259 * qJDD(1) + t263;
t223 = qJD(1) * t341 - t383;
t220 = qJDD(6) + t223;
t363 = t305 * t220;
t390 = -t336 * t388 - t363;
t389 = t301 * MDP(8) + t303 * MDP(9);
t385 = qJDD(1) * t276;
t207 = t228 * t309 - t229 * t306;
t205 = -qJD(5) * pkin(5) - t207;
t376 = -pkin(7) + t269;
t249 = t376 * t301;
t250 = t376 * t303;
t214 = t249 * t306 - t250 * t309;
t209 = -qJD(4) * t259 - qJD(5) * t214;
t258 = t301 * t306 - t309 * t303;
t261 = t301 * pkin(4) + t272;
t213 = pkin(5) * t259 + pkin(8) * t258 + t261;
t215 = t249 * t309 + t250 * t306;
t356 = qJD(5) * t306;
t254 = -t301 * t355 - t303 * t356;
t324 = t226 * t306 + t227 * t309;
t196 = qJDD(5) * pkin(8) + qJD(5) * t207 + t324;
t262 = t272 * qJD(1);
t260 = qJD(4) + t262;
t248 = pkin(4) * t359 + t260;
t211 = pkin(5) * t251 - pkin(8) * t253 + t248;
t334 = qJD(6) * t211 + t196;
t381 = -t197 * t258 + t205 * t254 - t215 * t220 - (qJD(6) * t213 + t209) * t388 - t259 * t334;
t377 = g(3) * t286;
t345 = t305 * qJDD(5) + t392 * t308;
t354 = qJD(6) * t305;
t203 = -t253 * t354 + t345;
t375 = t203 * t258;
t374 = t203 * t305;
t373 = t205 * t258;
t372 = t213 * t220;
t371 = t232 * t253;
t370 = t234 * t253;
t369 = t285 * t305;
t368 = t285 * t308;
t367 = t287 * t305;
t366 = t287 * t308;
t216 = t308 * t220;
t255 = -t301 * t356 + t341;
t362 = t203 * t259 + t234 * t255;
t360 = t301 ^ 2 + t303 ^ 2;
t353 = qJD(6) * t308;
t298 = qJD(3) * qJD(1);
t347 = t258 * t363;
t346 = t258 * t216;
t310 = cos(qJ(1));
t344 = t310 * pkin(1) + t287 * pkin(2) + t285 * qJ(3);
t343 = -t298 - t361;
t307 = sin(qJ(1));
t339 = -pkin(1) * t307 + t287 * qJ(3);
t335 = t360 * MDP(10);
t256 = qJDD(4) - t343;
t243 = pkin(4) * t349 + t256;
t201 = pkin(5) * t223 - pkin(8) * t222 + t243;
t206 = qJD(5) * pkin(8) + t208;
t333 = qJD(6) * t206 - t201;
t331 = -g(1) * t287 - g(2) * t285;
t329 = g(1) * t307 - g(2) * t310;
t328 = qJDD(3) - t386;
t289 = t308 * qJDD(5);
t204 = qJD(6) * t234 + t222 * t305 - t289;
t326 = -t204 * t259 - t232 * t255;
t323 = t230 * t303 + t231 * t301;
t322 = t235 * t303 + t236 * t301;
t224 = qJD(5) * t254 - qJDD(5) * t258;
t225 = -qJD(5) * t255 - qJDD(5) * t259;
t321 = t216 + (-t251 * t305 - t354) * t388;
t320 = -t254 * t305 + t258 * t353;
t319 = t254 * t308 + t258 * t354;
t318 = t331 + t361;
t316 = -t323 + t386;
t314 = -pkin(8) * t220 + (t205 + t207) * t388;
t311 = qJD(1) ^ 2;
t242 = t284 * t366 - t369;
t241 = t284 * t367 + t368;
t240 = t284 * t368 + t367;
t239 = -t284 * t369 + t366;
t218 = pkin(5) * t255 - pkin(8) * t254 + qJD(3);
t210 = -qJD(4) * t258 + qJD(5) * t215;
t200 = t308 * t201;
t199 = t206 * t308 + t211 * t305;
t198 = -t206 * t305 + t211 * t308;
t1 = [(-g(1) * t242 - g(2) * t240 + t198 * t255 + t200 * t259 + t214 * t204 + t210 * t232 + (t372 + t218 * t388 + (-t206 * t259 - t215 * t388 - t373) * qJD(6)) * t308 + t381 * t305) * MDP(24) + (g(1) * t241 - g(2) * t239 - t199 * t255 + t214 * t203 + t210 * t234 + (-(-qJD(6) * t215 + t218) * t388 - t372 + t333 * t259 + qJD(6) * t373) * t305 + t381 * t308) * MDP(25) + qJDD(1) * MDP(1) + (t329 + (t302 ^ 2 + t304 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t328 + 0.2e1 * t385) * MDP(5) + (0.2e1 * t298 + t318 + t361) * MDP(6) + (-t343 * t272 + t262 * qJD(3) + (qJDD(3) + t385) * t276 - g(1) * (-pkin(2) * t285 + t339) - g(2) * t344) * MDP(7) + (t360 * t382 + t316) * MDP(10) + (t256 * t272 + t260 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t285 + t339) - g(2) * (qJ(4) * t287 + t344) + t323 * t269 - t322 * qJD(4)) * MDP(11) + (-t222 * t258 + t253 * t254) * MDP(12) + (-t222 * t259 + t223 * t258 - t251 * t254 - t253 * t255) * MDP(13) + t224 * MDP(14) + t225 * MDP(15) + (qJD(3) * t251 - qJD(5) * t210 - qJDD(5) * t214 + t223 * t261 + t243 * t259 + t248 * t255 + t284 * t331) * MDP(17) + (qJD(3) * t253 - qJD(5) * t209 - qJDD(5) * t215 + t222 * t261 - t243 * t258 + t248 * t254 + t286 * t331) * MDP(18) + (t234 * t319 - t308 * t375) * MDP(19) + ((-t232 * t308 - t234 * t305) * t254 + (t374 + t204 * t308 + (-t232 * t305 + t234 * t308) * qJD(6)) * t258) * MDP(20) + (t319 * t388 - t346 + t362) * MDP(21) + (t320 * t388 + t326 + t347) * MDP(22) + (t220 * t259 + t255 * t388) * MDP(23) + t329 * MDP(2) + (g(1) * t310 + g(2) * t307) * MDP(3) + t389 * (t256 + t318 + t298); (-t230 * t301 + t231 * t303 - g(3)) * MDP(11) + t225 * MDP(17) - t224 * MDP(18) + (-t326 + t347) * MDP(24) + (t346 + t362) * MDP(25) + (MDP(4) + MDP(7)) * (qJDD(2) - g(3)) + (MDP(24) * t320 - MDP(25) * t319) * t388; (-qJD(1) * t262 + t328) * MDP(7) + (-qJD(1) * t260 - t316) * MDP(11) + (-qJD(1) * t251 + t224) * MDP(17) + (-qJD(1) * t253 + t225) * MDP(18) + (t204 * t258 - t232 * t254 - t259 * t363) * MDP(24) + (-t259 * t216 - t234 * t254 + t375) * MDP(25) + (-MDP(6) - t389) * t311 + ((-qJD(1) * t308 - t255 * t305 - t259 * t353) * MDP(24) + (qJD(1) * t305 - t255 * t308 + t259 * t354) * MDP(25)) * t388 + (MDP(7) * t276 + MDP(5) - t335) * qJDD(1); (qJD(1) * t322 + t256 + t331) * MDP(11) - t263 * MDP(17) + t264 * MDP(18) + (t321 - t371) * MDP(24) + (-t370 + t390) * MDP(25) - t311 * t335 + ((MDP(17) * t306 + MDP(9)) * t303 + (MDP(17) * t309 - MDP(18) * t306 + MDP(8)) * t301) * qJDD(1) + ((t253 + t342) * MDP(17) + (-t301 * t357 - t303 * t358 - t251) * MDP(18)) * qJD(5); t253 * t251 * MDP(12) + (-t251 ^ 2 + t253 ^ 2) * MDP(13) + t327 * MDP(14) + ((t253 - t342) * qJD(5) + t383) * MDP(15) + qJDD(5) * MDP(16) + (-t248 * t253 + t315 - t325) * MDP(17) + (t248 * t251 + t284 * t386 - t324 + t377) * MDP(18) + (t234 * t336 + t374) * MDP(19) + ((t203 - t394) * t308 + (-t204 - t393) * t305) * MDP(20) + (-t370 - t390) * MDP(21) + (t321 + t371) * MDP(22) - t388 * t253 * MDP(23) + (-pkin(5) * t204 - t198 * t253 - t208 * t232 + t314 * t305 + t391 * t308) * MDP(24) + (-pkin(5) * t203 + t199 * t253 - t208 * t234 - t391 * t305 + t314 * t308) * MDP(25); t234 * t232 * MDP(19) + (-t232 ^ 2 + t234 ^ 2) * MDP(20) + (t345 + t394) * MDP(21) + (t289 + t393) * MDP(22) + t220 * MDP(23) + (-g(1) * t239 - g(2) * t241 + t199 * t388 - t205 * t234 + t200) * MDP(24) + (g(1) * t240 - g(2) * t242 + t198 * t388 + t205 * t232) * MDP(25) + ((-t196 + t377) * MDP(25) + (-t253 * MDP(22) - MDP(24) * t206 - MDP(25) * t211) * qJD(6)) * t308 + (-qJD(6) * t253 * MDP(21) - t392 * MDP(22) + (-t334 + t377) * MDP(24) + t333 * MDP(25)) * t305;];
tau  = t1;
