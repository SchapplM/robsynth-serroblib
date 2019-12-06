% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:34
% EndTime: 2019-12-05 18:16:38
% DurationCPUTime: 1.49s
% Computational Cost: add. (1347->210), mult. (2211->280), div. (0->0), fcn. (1414->14), ass. (0->128)
t386 = pkin(7) + pkin(8);
t301 = sin(qJ(5));
t302 = sin(qJ(4));
t305 = cos(qJ(5));
t306 = cos(qJ(4));
t242 = t301 * t306 + t302 * t305;
t295 = qJD(1) + qJD(3);
t236 = t242 * t295;
t300 = cos(pkin(9));
t279 = pkin(1) * t300 + pkin(2);
t299 = sin(pkin(9));
t379 = pkin(1) * t299;
t350 = qJD(1) * t379;
t385 = -qJD(3) * t350 + t279 * qJDD(1);
t286 = qJ(1) + pkin(9) + qJ(3);
t277 = sin(t286);
t278 = cos(t286);
t384 = g(2) * t278 + g(3) * t277;
t382 = -g(2) * t277 + g(3) * t278;
t263 = t279 * qJD(1);
t383 = qJD(3) * t263 + qJDD(1) * t379;
t360 = t305 * t306;
t364 = t301 * t302;
t241 = -t360 + t364;
t303 = sin(qJ(3));
t307 = cos(qJ(3));
t338 = t279 * t307 - t303 * t379;
t357 = t303 * t279 + t307 * t379;
t230 = t263 * t303 + t307 * t350;
t341 = t386 * t295 + t230;
t215 = t306 * qJD(2) - t341 * t302;
t294 = qJD(4) + qJD(5);
t216 = qJD(2) * t302 + t341 * t306;
t381 = t385 * t303 + t383 * t307;
t293 = qJDD(1) + qJDD(3);
t378 = pkin(3) * t293;
t377 = pkin(3) * t295;
t376 = pkin(4) * t306;
t238 = pkin(7) + t357;
t374 = -pkin(8) - t238;
t373 = t216 * t305;
t229 = t263 * t307 - t303 * t350;
t372 = t229 * t294;
t371 = t230 * t295;
t233 = t357 * qJD(3);
t370 = t233 * t295;
t298 = qJ(4) + qJ(5);
t289 = cos(t298);
t369 = t277 * t289;
t368 = t278 * t289;
t366 = t293 * t306;
t365 = t295 * t302;
t361 = t302 * t306;
t359 = qJDD(2) - g(1);
t335 = -t383 * t303 + t385 * t307;
t211 = -t335 - t378;
t223 = -t229 - t377;
t352 = qJD(4) * t306;
t358 = t211 * t302 + t223 * t352;
t296 = t302 ^ 2;
t356 = -t306 ^ 2 + t296;
t353 = qJD(4) * t302;
t351 = qJD(5) * t301;
t349 = pkin(4) * t353;
t347 = t295 * t364;
t346 = t295 * t360;
t345 = t223 * t353 + t384 * t306;
t283 = -pkin(3) - t376;
t344 = qJD(4) * t386;
t343 = t295 * t352;
t210 = pkin(7) * t293 + t381;
t342 = pkin(8) * t293 + t210;
t340 = qJD(4) * t374;
t199 = (t295 * t353 - t366) * pkin(4) + t211;
t217 = t283 * t295 - t229;
t222 = t294 * t242;
t336 = g(2) * t368 + g(3) * t369 + t199 * t241 + t217 * t222;
t237 = -pkin(3) - t338;
t334 = -t230 + t349;
t304 = sin(qJ(1));
t308 = cos(qJ(1));
t332 = g(2) * t308 + g(3) * t304;
t331 = t241 * t293;
t214 = qJD(4) * pkin(4) + t215;
t329 = -t214 * t301 - t373;
t221 = t294 * t241;
t292 = qJDD(4) + qJDD(5);
t328 = -t221 * t294 + t242 * t292;
t225 = t374 * t302;
t290 = t306 * pkin(8);
t226 = t238 * t306 + t290;
t327 = t225 * t305 - t226 * t301;
t326 = t225 * t301 + t226 * t305;
t268 = t386 * t302;
t269 = pkin(7) * t306 + t290;
t325 = -t268 * t305 - t269 * t301;
t324 = -t268 * t301 + t269 * t305;
t323 = t335 + t384;
t309 = qJD(4) ^ 2;
t321 = -pkin(7) * t309 + t371 + t378;
t320 = -t237 * t293 - t238 * t309 - t370;
t202 = qJD(5) * t346 + t242 * t293 - t294 * t347 + t305 * t343;
t234 = -t346 + t347;
t319 = t236 * t234 * MDP(15) + (t234 * t294 + t202) * MDP(17) - t331 * MDP(18) + (-t234 ^ 2 + t236 ^ 2) * MDP(16) + t292 * MDP(19);
t318 = -pkin(7) * qJDD(4) + (t229 - t377) * qJD(4);
t317 = -t223 * t295 - t210 + t382;
t232 = t338 * qJD(3);
t316 = -qJDD(4) * t238 + (t237 * t295 - t232) * qJD(4);
t288 = sin(t298);
t315 = t199 * t242 - t217 * t221 - t288 * t384;
t203 = t222 * t295 + t331;
t212 = -t222 * t294 - t241 * t292;
t253 = qJDD(4) * t302 + t306 * t309;
t254 = qJDD(4) * t306 - t302 * t309;
t314 = (-t202 * t241 - t203 * t242 + t221 * t234 - t222 * t236) * MDP(16) + (t202 * t242 - t221 * t236) * MDP(15) + t328 * MDP(17) + t212 * MDP(18) + 0.2e1 * (-t356 * t295 * qJD(4) + t293 * t361) * MDP(9) + (t293 * t296 + 0.2e1 * t302 * t343) * MDP(8) + t253 * MDP(10) + t254 * MDP(11) + t293 * MDP(5);
t285 = t306 * qJDD(2);
t193 = qJDD(4) * pkin(4) - t216 * qJD(4) - t342 * t302 + t285;
t313 = -g(2) * t369 + t217 * t234 + t216 * t351 + g(3) * t368 + g(1) * t288 + (-t216 * t294 - t193) * t301;
t312 = -t381 + t382;
t194 = t215 * qJD(4) + t302 * qJDD(2) + t342 * t306;
t311 = -g(1) * t289 + t329 * qJD(5) + t305 * t193 - t301 * t194 - t217 * t236 + t382 * t288;
t250 = t306 * t344;
t249 = t302 * t344;
t228 = t237 - t376;
t227 = t233 + t349;
t209 = -t232 * t302 + t306 * t340;
t208 = t232 * t306 + t302 * t340;
t1 = [(t316 * t306 + (-t320 - t384) * t302 + t358) * MDP(14) + (t332 + (t299 ^ 2 + t300 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + qJDD(1) * MDP(1) + t314 + t332 * MDP(2) + (-g(2) * t304 + g(3) * t308) * MDP(3) + (t338 * t293 + t323 - t370) * MDP(6) + (t316 * t302 + (-t211 + t320) * t306 + t345) * MDP(13) + (t227 * t234 + t228 * t203 + (-t326 * qJD(5) - t208 * t301 + t209 * t305) * t294 + t327 * t292 + t336) * MDP(20) + (-t232 * t295 - t357 * t293 + t312) * MDP(7) + (t227 * t236 + t228 * t202 - (t327 * qJD(5) + t208 * t305 + t209 * t301) * t294 - t326 * t292 + t315) * MDP(21); t254 * MDP(13) - t253 * MDP(14) + t212 * MDP(20) - t328 * MDP(21) + t359 * MDP(4); (t318 * t306 + (-t321 - t384) * t302 + t358) * MDP(14) + (t318 * t302 + (-t211 + t321) * t306 + t345) * MDP(13) + t314 + (t323 + t371) * MDP(6) + (t283 * t203 + (-t324 * qJD(5) + t249 * t301 - t250 * t305) * t294 + t325 * t292 + t334 * t234 + t242 * t372 + t336) * MDP(20) + (t229 * t295 + t312) * MDP(7) + (t283 * t202 - (t325 * qJD(5) - t249 * t305 - t250 * t301) * t294 - t324 * t292 + t334 * t236 - t241 * t372 + t315) * MDP(21); t302 * t293 * MDP(10) + MDP(11) * t366 + qJDD(4) * MDP(12) + (-g(1) * t306 + t317 * t302 + t285) * MDP(13) + (-t359 * t302 + t317 * t306) * MDP(14) + (-(-t215 * t301 - t373) * t294 + (-t234 * t365 + t292 * t305 - t294 * t351) * pkin(4) + t311) * MDP(20) + ((-qJD(5) * t214 + t215 * t294 - t194) * t305 + (-qJD(5) * t294 * t305 - t236 * t365 - t292 * t301) * pkin(4) + t313) * MDP(21) + t319 + (-MDP(8) * t361 + t356 * MDP(9)) * t295 ^ 2; (-t329 * t294 + t311) * MDP(20) + ((-t194 + (-qJD(5) + t294) * t214) * t305 + t313) * MDP(21) + t319;];
tau = t1;
