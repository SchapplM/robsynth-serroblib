% Calculate vector of inverse dynamics joint torques for
% S5RPRRP6
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:22
% EndTime: 2019-12-31 18:43:26
% DurationCPUTime: 2.37s
% Computational Cost: add. (1561->308), mult. (3214->414), div. (0->0), fcn. (1960->10), ass. (0->143)
t284 = sin(pkin(8));
t269 = pkin(1) * t284 + pkin(6);
t260 = t269 * qJDD(1);
t386 = -qJD(2) * qJD(3) - t260;
t262 = t269 * qJD(1);
t288 = sin(qJ(3));
t291 = cos(qJ(3));
t343 = qJD(3) * t291;
t329 = -t262 * t343 + t386 * t288;
t310 = -qJDD(3) * pkin(3) - t329;
t332 = qJDD(2) * t291;
t211 = t310 - t332;
t349 = qJD(1) * t291;
t268 = -qJD(4) + t349;
t281 = qJ(1) + pkin(8);
t274 = sin(t281);
t275 = cos(t281);
t318 = g(1) * t275 + g(2) * t274;
t305 = t318 * t288;
t385 = qJD(4) * pkin(7) * t268 - g(3) * t291 - t211 + t305;
t232 = qJD(2) * t291 - t288 * t262;
t384 = qJD(3) * t232;
t341 = qJD(4) * t288;
t383 = qJD(1) * t341 - qJDD(3);
t287 = sin(qJ(4));
t290 = cos(qJ(4));
t333 = qJDD(1) * t288;
t218 = ((qJD(4) + t349) * qJD(3) + t333) * t287 + t383 * t290;
t361 = qJDD(2) - g(3);
t382 = t361 * t291;
t364 = t287 * t291;
t227 = t274 * t364 + t275 * t290;
t229 = t274 * t290 - t275 * t364;
t381 = -g(1) * t229 + g(2) * t227;
t345 = qJD(3) * t287;
t350 = qJD(1) * t288;
t253 = t290 * t350 + t345;
t379 = t253 ^ 2;
t378 = pkin(4) * t287;
t373 = g(3) * t288;
t372 = qJ(5) + pkin(7);
t331 = t291 * qJDD(1);
t335 = qJD(1) * qJD(3);
t247 = t288 * t335 + qJDD(4) - t331;
t371 = t247 * t287;
t336 = t290 * qJD(3);
t251 = t287 * t350 - t336;
t370 = t251 * t268;
t369 = t253 * t268;
t368 = t253 * t288;
t367 = t268 * t290;
t366 = t269 * t287;
t365 = t287 * t288;
t363 = t288 * t290;
t362 = t290 * t291;
t233 = t288 * qJD(2) + t291 * t262;
t224 = qJD(3) * pkin(7) + t233;
t285 = cos(pkin(8));
t270 = -pkin(1) * t285 - pkin(2);
t243 = -pkin(3) * t291 - pkin(7) * t288 + t270;
t225 = t243 * qJD(1);
t205 = -t224 * t287 + t290 * t225;
t202 = -qJ(5) * t253 + t205;
t201 = -pkin(4) * t268 + t202;
t360 = -t202 + t201;
t327 = t291 * t336;
t359 = -t218 * t363 - t251 * t327;
t319 = pkin(3) * t288 - pkin(7) * t291;
t256 = t319 * qJD(1);
t358 = t290 * t232 + t287 * t256;
t257 = t319 * qJD(3);
t340 = qJD(4) * t290;
t357 = t243 * t340 + t287 * t257;
t323 = qJD(4) * t372;
t339 = qJD(5) * t290;
t356 = t339 - t358 + (qJ(5) * t349 - t323) * t287;
t241 = t290 * t256;
t308 = pkin(4) * t288 - qJ(5) * t362;
t355 = -t308 * qJD(1) - t290 * t323 - t241 + (-qJD(5) + t232) * t287;
t344 = qJD(3) * t288;
t354 = t290 * t257 + t344 * t366;
t255 = t269 * t362;
t353 = t287 * t243 + t255;
t282 = t288 ^ 2;
t352 = -t291 ^ 2 + t282;
t351 = MDP(19) * t287;
t263 = qJD(1) * t270;
t347 = qJD(3) * t251;
t346 = qJD(3) * t269;
t342 = qJD(4) * t287;
t223 = -qJD(3) * pkin(3) - t232;
t216 = pkin(4) * t251 + qJD(5) + t223;
t338 = t216 * qJD(3);
t337 = t290 * MDP(18);
t210 = qJDD(3) * pkin(7) + qJDD(2) * t288 + t260 * t291 + t384;
t219 = qJD(1) * t257 + t243 * qJDD(1);
t330 = t290 * t210 + t287 * t219 + t225 * t340;
t328 = pkin(6) + t378;
t326 = t291 * t335;
t324 = t269 + t378;
t322 = t268 * t269 + t224;
t321 = -qJD(4) * t225 - t210;
t317 = g(1) * t274 - g(2) * t275;
t289 = sin(qJ(1));
t292 = cos(qJ(1));
t316 = g(1) * t289 - g(2) * t292;
t206 = t224 * t290 + t225 * t287;
t214 = t290 * t219;
t217 = -qJD(4) * t336 + (-t326 - t333) * t290 + t383 * t287;
t196 = pkin(4) * t247 + qJ(5) * t217 - t206 * qJD(4) - qJD(5) * t253 - t287 * t210 + t214;
t302 = -t224 * t342 + t330;
t197 = -qJ(5) * t218 - qJD(5) * t251 + t302;
t315 = -t196 * t287 + t197 * t290;
t203 = -qJ(5) * t251 + t206;
t314 = t201 * t290 + t203 * t287;
t313 = t201 * t287 - t203 * t290;
t273 = pkin(4) * t290 + pkin(3);
t312 = t273 * t291 + t288 * t372;
t309 = t316 * pkin(1);
t306 = pkin(2) + t312;
t304 = -t268 * t340 + t371;
t301 = -qJD(1) * t263 + t318;
t300 = -t287 * t341 + t327;
t299 = -pkin(7) * t247 - t268 * t223;
t298 = 0.2e1 * t263 * qJD(3) - qJDD(3) * t269;
t297 = pkin(4) * t218 + qJDD(5) + t310;
t293 = qJD(3) ^ 2;
t296 = -0.2e1 * qJDD(1) * t270 - t269 * t293 + t317;
t265 = t372 * t290;
t264 = t372 * t287;
t259 = qJDD(3) * t291 - t288 * t293;
t258 = qJDD(3) * t288 + t291 * t293;
t246 = t251 ^ 2;
t236 = t253 * t344;
t235 = t290 * t243;
t230 = t274 * t287 + t275 * t362;
t228 = -t274 * t362 + t275 * t287;
t215 = -qJ(5) * t365 + t353;
t212 = -qJ(5) * t363 + t235 + (-pkin(4) - t366) * t291;
t200 = t297 - t332;
t199 = (-qJ(5) * qJD(4) - t346) * t363 + (-qJD(5) * t288 + (-qJ(5) * qJD(3) - qJD(4) * t269) * t291) * t287 + t357;
t198 = -t288 * t339 + t308 * qJD(3) + (-t255 + (qJ(5) * t288 - t243) * t287) * qJD(4) + t354;
t1 = [qJDD(1) * MDP(1) + t316 * MDP(2) + (g(1) * t292 + g(2) * t289) * MDP(3) + ((t284 ^ 2 + t285 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t309) * MDP(4) + (qJDD(1) * t282 + 0.2e1 * t288 * t326) * MDP(5) + 0.2e1 * (t288 * t331 - t352 * t335) * MDP(6) + t258 * MDP(7) + t259 * MDP(8) + (t298 * t288 + t296 * t291) * MDP(10) + (-t296 * t288 + t298 * t291) * MDP(11) + (-t217 * t363 + t300 * t253) * MDP(12) + (-t340 * t368 + (-t253 * t343 + (qJD(4) * t251 + t217) * t288) * t287 + t359) * MDP(13) + (t217 * t291 + t247 * t363 - t300 * t268 + t236) * MDP(14) + ((t268 * t345 + t218) * t291 + (-t304 - t347) * t288) * MDP(15) + (-t247 * t291 - t268 * t344) * MDP(16) + (-(-t243 * t342 + t354) * t268 + t235 * t247 - g(1) * t228 - g(2) * t230 + (t251 * t346 - t214 + t322 * t340 + (qJD(3) * t223 - t247 * t269 - t321) * t287) * t291 + (t205 * qJD(3) + t211 * t287 + t269 * t218 + t223 * t340) * t288) * MDP(17) + (t357 * t268 - t353 * t247 - g(1) * t227 - g(2) * t229 + (-t322 * t342 + (t223 * t290 + t253 * t269) * qJD(3) + t330) * t291 + (-t223 * t342 + t211 * t290 - t269 * t217 + (-t269 * t367 - t206) * qJD(3)) * t288) * MDP(18) + (-t198 * t253 - t199 * t251 + t212 * t217 - t215 * t218 - t314 * t343 + (t313 * qJD(4) - t196 * t290 - t197 * t287 + t317) * t288) * MDP(19) + (t196 * t212 + t197 * t215 + t201 * t198 + t203 * t199 + t324 * t291 * t338 + (t216 * pkin(4) * t340 + t200 * t324) * t288 + t309 + (-g(1) * t328 - g(2) * t306) * t275 + (g(1) * t306 - g(2) * t328) * t274) * MDP(20); t361 * MDP(4) + t259 * MDP(10) - t258 * MDP(11) + t236 * MDP(18) + t359 * MDP(19) - g(3) * MDP(20) + (-t218 * MDP(17) + t217 * MDP(18) - t200 * MDP(20) + (t253 * t351 - t313 * MDP(20) + (t287 * MDP(17) + t337) * t268) * qJD(3)) * t291 + ((t347 - t371) * MDP(17) - t247 * t337 - t217 * t351 + (t315 + t338) * MDP(20) + ((t251 * t287 + t253 * t290) * MDP(19) - t314 * MDP(20) + (t290 * MDP(17) - t287 * MDP(18)) * t268) * qJD(4)) * t288; MDP(7) * t333 + MDP(8) * t331 + qJDD(3) * MDP(9) + (qJD(3) * t233 + t301 * t288 + t329 + t382) * MDP(10) + (t384 + (qJD(3) * t262 - t361) * t288 + (t301 + t386) * t291) * MDP(11) + (-t217 * t287 - t253 * t367) * MDP(12) + ((-t217 + t370) * t290 + (-t218 + t369) * t287) * MDP(13) + ((t268 * t362 - t368) * qJD(1) + t304) * MDP(14) + (t268 * t342 + t247 * t290 + (t251 * t288 - t268 * t364) * qJD(1)) * MDP(15) + t268 * MDP(16) * t350 + (-t205 * t350 - pkin(3) * t218 - t233 * t251 + t241 * t268 + (-t232 * t268 + t299) * t287 + t385 * t290) * MDP(17) + (pkin(3) * t217 + t206 * t350 - t233 * t253 - t358 * t268 - t385 * t287 + t299 * t290) * MDP(18) + (-t373 - t217 * t264 - t218 * t265 - t355 * t253 - t356 * t251 - t314 * qJD(4) + (t314 * qJD(1) - t318) * t291 + t315) * MDP(19) + (t197 * t265 - t196 * t264 - t200 * t273 - g(3) * t312 + (-t268 * t378 - t233) * t216 + t356 * t203 + t355 * t201 + t318 * (t273 * t288 - t291 * t372)) * MDP(20) + (-t288 * t291 * MDP(5) + t352 * MDP(6)) * qJD(1) ^ 2; t253 * t251 * MDP(12) + (-t246 + t379) * MDP(13) + (-t217 - t370) * MDP(14) + (-t218 - t369) * MDP(15) + t247 * MDP(16) + (-t224 * t340 - t206 * t268 - t223 * t253 + t214 + (t321 + t373) * t287 + t381) * MDP(17) + (g(1) * t230 - g(2) * t228 + g(3) * t363 - t205 * t268 + t223 * t251 - t302) * MDP(18) + (pkin(4) * t217 - t360 * t251) * MDP(19) + (t360 * t203 + (g(3) * t365 - t216 * t253 + t196 + t381) * pkin(4)) * MDP(20); (-t246 - t379) * MDP(19) + (t201 * t253 + t203 * t251 + t297 - t305 - t382) * MDP(20);];
tau = t1;
