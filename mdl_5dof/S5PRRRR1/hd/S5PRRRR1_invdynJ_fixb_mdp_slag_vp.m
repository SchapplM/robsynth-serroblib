% Calculate vector of inverse dynamics joint torques for
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:57
% EndTime: 2019-07-18 13:29:03
% DurationCPUTime: 3.38s
% Computational Cost: add. (1321->290), mult. (3253->417), div. (0->0), fcn. (2642->10), ass. (0->128)
t264 = qJ(3) + qJ(4);
t258 = sin(t264);
t259 = cos(t264);
t268 = sin(qJ(2));
t272 = cos(qJ(2));
t347 = g(1) * t272;
t354 = g(3) * t268 + t347;
t352 = g(2) * t259 + t354 * t258;
t266 = sin(qJ(4));
t270 = cos(qJ(4));
t271 = cos(qJ(3));
t313 = qJDD(2) * t271;
t267 = sin(qJ(3));
t314 = qJDD(2) * t267;
t229 = t266 * t271 + t267 * t270;
t261 = qJD(3) + qJD(4);
t355 = t261 * t229;
t196 = qJD(2) * t355 + t266 * t314 - t270 * t313;
t315 = qJDD(1) * t268;
t317 = qJD(1) * qJD(2);
t285 = t272 * t317 + t315;
t353 = t285 * t271;
t194 = qJDD(5) + t196;
t323 = qJD(2) * t271;
t251 = t270 * t323;
t325 = qJD(2) * t267;
t308 = t266 * t325;
t224 = -t251 + t308;
t223 = qJD(5) + t224;
t322 = qJD(3) * t267;
t326 = qJD(1) * t268;
t351 = (-t194 * t271 + t223 * t322) * pkin(2) - t223 * t326;
t311 = t267 * t326;
t235 = qJD(3) * pkin(2) - t311;
t321 = qJD(4) * t235;
t278 = t321 + t353;
t307 = qJD(4) * t326;
t301 = t271 * t307;
t205 = -t267 * t315 + qJDD(3) * pkin(2) + (-qJD(3) * t268 * t271 - t272 * t325) * qJD(1);
t302 = qJD(3) * t311;
t330 = -t270 * t205 - t266 * t302;
t283 = t270 * t301 + t330;
t186 = t278 * t266 + t283;
t228 = t266 * t267 - t270 * t271;
t201 = t261 * t228;
t309 = t271 * t326;
t210 = -t270 * t235 + t266 * t309;
t286 = t228 * t272;
t216 = qJD(1) * t286;
t299 = -t270 * t302 + (t205 - t301) * t266;
t185 = t278 * t270 + t299;
t236 = -pkin(2) * t323 - qJD(1) * t272;
t304 = qJD(5) * t236 + t185;
t350 = -t304 * t228 + t186 * t229 - t210 * t201 + (qJD(5) * t271 * pkin(2) - t216) * t223;
t348 = g(1) * t268;
t346 = g(2) * t258;
t226 = -t266 * t323 - t270 * t325;
t316 = qJD(2) * qJD(3);
t306 = t271 * t316;
t195 = qJD(4) * t251 - t261 * t308 + t266 * t313 + (t306 + t314) * t270;
t260 = qJDD(3) + qJDD(4);
t265 = sin(qJ(5));
t269 = cos(qJ(5));
t318 = qJD(5) * t269;
t312 = t269 * t195 + t265 * t260 + t261 * t318;
t319 = qJD(5) * t265;
t187 = t226 * t319 + t312;
t344 = t187 * t265;
t343 = t194 * t265;
t342 = t194 * t269;
t207 = -t226 * t265 - t269 * t261;
t341 = t207 * t223;
t209 = -t226 * t269 + t261 * t265;
t340 = t209 * t223;
t339 = t210 * t224;
t338 = t210 * t229;
t305 = t223 * t269;
t337 = t229 * t269;
t336 = t265 * t268;
t335 = t265 * t272;
t334 = t268 * t269;
t333 = t269 * t272;
t332 = qJDD(1) - g(3);
t331 = t186 * t265 + t210 * t318;
t262 = t267 ^ 2;
t329 = -t271 ^ 2 + t262;
t328 = MDP(12) * t224;
t327 = MDP(23) * t223;
t324 = qJD(2) * t268;
t320 = qJD(4) * t270;
t237 = t270 * t309;
t211 = t235 * t266 + t237;
t212 = t268 * t317 - qJDD(1) * t272 + (t267 * t316 - t313) * pkin(2);
t303 = qJD(5) * t211 - t212;
t300 = -g(3) * t272 + t348;
t215 = t229 * t326;
t297 = -t215 * t223 + t339;
t218 = t228 * t268;
t296 = -t218 * t269 - t335;
t295 = t218 * t265 - t333;
t294 = qJD(5) * t266 + t325;
t274 = qJD(2) ^ 2;
t293 = qJDD(2) * t272 - t268 * t274;
t291 = -t201 * t269 - t229 * t319;
t289 = -t332 * t268 + t347;
t288 = t210 * t319 + (-t186 + t352) * t269;
t287 = t229 * t272;
t284 = -qJDD(3) * t268 - 0.2e1 * t272 * t316;
t253 = t269 * t260;
t188 = t209 * qJD(5) + t195 * t265 - t253;
t282 = ((t187 - t341) * t269 + (-t188 - t340) * t265) * MDP(20) + (t209 * t305 + t344) * MDP(19) + (-t223 ^ 2 * t265 - t207 * t226 + t342) * MDP(22) + (t209 * t226 + t223 * t305 + t343) * MDP(21) + (t224 * t261 + t195) * MDP(14) + (-t226 * t261 - t196) * MDP(15) + (-t224 ^ 2 + t226 ^ 2) * MDP(13) + t260 * MDP(16);
t273 = qJD(3) ^ 2;
t281 = -t268 * t273 + t293;
t279 = t224 * t236 + t354 * t259 - t299 - t346;
t276 = -t353 + (-pkin(2) * t261 - t235) * qJD(4);
t222 = t259 * t333 + t336;
t221 = -t259 * t335 + t334;
t220 = -t259 * t334 + t335;
t219 = t259 * t336 + t333;
t217 = t229 * t268;
t214 = qJD(1) * t287;
t213 = t266 * t311 - t237;
t206 = t269 * t212;
t199 = t211 * t269 + t236 * t265;
t198 = -t211 * t265 + t236 * t269;
t192 = qJD(2) * t287 - t201 * t268;
t191 = -qJD(2) * t286 - t268 * t355;
t1 = [t332 * MDP(1) + t293 * MDP(3) + (-qJDD(2) * t268 - t272 * t274) * MDP(4) + (-t192 * t261 - t196 * t272 - t217 * t260 + t224 * t324) * MDP(17) + (-t191 * t261 - t195 * t272 + t218 * t260 - t226 * t324) * MDP(18) + ((-t296 * qJD(5) - t191 * t265 + t269 * t324) * t223 + t295 * t194 + t192 * t207 + t217 * t188) * MDP(24) + (-(t295 * qJD(5) + t191 * t269 + t265 * t324) * t223 - t296 * t194 + t192 * t209 + t217 * t187) * MDP(25) + (t281 * MDP(10) + t284 * MDP(11)) * t271 + (t284 * MDP(10) - t281 * MDP(11)) * t267; qJDD(2) * MDP(2) + t289 * MDP(4) + (qJDD(2) * t262 + 0.2e1 * t267 * t306) * MDP(5) + 0.2e1 * (t267 * t313 - t329 * t316) * MDP(6) + (qJDD(3) * t267 + t271 * t273) * MDP(7) + (qJDD(3) * t271 - t267 * t273) * MDP(8) + (t195 * t229 + t201 * t226) * MDP(12) + (-t195 * t228 - t196 * t229 + t201 * t224 + t226 * t355) * MDP(13) + (-t201 * t261 + t229 * t260) * MDP(14) + (-t228 * t260 - t261 * t355) * MDP(15) + (-t224 * t326 + t355 * t236 + t212 * t228 + t214 * t261 + t300 * t259 + (-t196 * t271 + t224 * t322) * pkin(2)) * MDP(17) + (t226 * t326 - t201 * t236 + t212 * t229 - t216 * t261 - t300 * t258 + (-t195 * t271 - t226 * t322) * pkin(2)) * MDP(18) + (t187 * t337 + t291 * t209) * MDP(19) + (-(-t207 * t269 - t209 * t265) * t201 + (-t344 - t188 * t269 + (t207 * t265 - t209 * t269) * qJD(5)) * t229) * MDP(20) + (t187 * t228 + t194 * t337 + t209 * t355 + t291 * t223) * MDP(21) + (-t229 * t343 - t188 * t228 - t355 * t207 + (t201 * t265 - t229 * t318) * t223) * MDP(22) + (t194 * t228 + t223 * t355) * MDP(23) + (-g(1) * t220 - g(3) * t222 + t198 * t355 + t206 * t228 - t214 * t207 + ((-t211 * t228 + t338) * qJD(5) + t351) * t269 + t350 * t265) * MDP(24) + (-g(1) * t219 - g(3) * t221 - t199 * t355 - t214 * t209 + t350 * t269 + (-qJD(5) * t338 + t303 * t228 - t351) * t265) * MDP(25) + (MDP(10) * t271 - MDP(11) * t267 + MDP(3)) * (t332 * t272 + t348); t282 + (-t213 * t261 + t226 * t236 + (-t224 * t325 + t260 * t270) * pkin(2) + t276 * t266 - t283 + t352) * MDP(17) + (-t215 * t261 + (t226 * t325 - t260 * t266) * pkin(2) + t276 * t270 + t279) * MDP(18) + (-t199 * t226 + t213 * t209 + t297 * t269 - t352 * t265 + (-t270 * t187 + (qJD(4) * t209 - t342) * t266 + (t294 * t265 - t269 * t320) * t223) * pkin(2) + t331) * MDP(25) + (t198 * t226 + t213 * t207 + t297 * t265 + (-t270 * t188 + (qJD(4) * t207 - t343) * t266 + (-t265 * t320 - t294 * t269) * t223) * pkin(2) + t288) * MDP(24) + qJDD(3) * MDP(9) + (g(2) * t271 + t289 * t267) * MDP(10) + (-g(2) * t267 + t289 * t271) * MDP(11) + MDP(8) * t313 + MDP(7) * t314 - t226 * t328 + t226 * t327 + (-t267 * t271 * MDP(5) + t329 * MDP(6)) * t274; (t211 * t261 - t266 * t321 - t330 + t352) * MDP(17) + (-t210 * t261 - t235 * t320 + t279) * MDP(18) + (-t207 * t211 + t288) * MDP(24) + (-t209 * t211 - t210 * t305 + t269 * t339 + t331) * MDP(25) + (MDP(17) * t236 + MDP(24) * t198 - MDP(25) * t199 + t327 - t328) * t226 + (-t285 * MDP(17) * t266 + (-MDP(17) * t307 - t285 * MDP(18)) * t270) * t271 + (-t352 * MDP(25) + (-t223 + t224) * MDP(24) * t210) * t265 + t282; t209 * t207 * MDP(19) + (-t207 ^ 2 + t209 ^ 2) * MDP(20) + (t312 + t341) * MDP(21) + (t253 + t340) * MDP(22) + t194 * MDP(23) + (-g(1) * t221 + g(3) * t219 + t199 * t223 - t209 * t210 + t206) * MDP(24) + (g(1) * t222 - g(3) * t220 + t198 * t223 + t207 * t210) * MDP(25) + ((-t185 - t346) * MDP(25) + (MDP(22) * t226 - MDP(24) * t211 - MDP(25) * t236) * qJD(5)) * t269 + (qJD(5) * t226 * MDP(21) + (-qJD(5) * t261 - t195) * MDP(22) + (-t304 - t346) * MDP(24) + t303 * MDP(25)) * t265;];
tau  = t1;
