% Calculate vector of inverse dynamics joint torques for
% S5RPRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:06
% EndTime: 2021-01-15 11:24:12
% DurationCPUTime: 2.46s
% Computational Cost: add. (1536->279), mult. (2776->325), div. (0->0), fcn. (1637->8), ass. (0->119)
t269 = sin(qJ(1));
t271 = cos(qJ(1));
t336 = g(1) * t269 - g(2) * t271;
t344 = qJDD(2) - t336;
t270 = cos(qJ(3));
t308 = qJD(1) * qJD(3);
t300 = t270 * t308;
t268 = sin(qJ(3));
t306 = qJDD(1) * t268;
t343 = t300 + t306;
t261 = qJDD(1) * qJ(2);
t262 = (qJD(1) * qJD(2));
t295 = g(1) * t271 + g(2) * t269;
t284 = -t295 + (2 * t262);
t342 = 0.2e1 * t261 + t284;
t324 = qJDD(1) * pkin(1);
t341 = t324 - t344;
t266 = sin(pkin(7));
t267 = cos(pkin(7));
t223 = t266 * t270 + t267 * t268;
t335 = t223 * qJD(1);
t340 = t335 * qJD(3);
t339 = MDP(16) + MDP(19);
t272 = -pkin(1) - pkin(6);
t338 = qJ(4) - t272;
t334 = MDP(14) + MDP(18);
t333 = -MDP(15) + MDP(20);
t232 = qJDD(1) * t272 + qJDD(2);
t224 = t270 * t232;
t233 = qJD(1) * t272 + qJD(2);
t301 = t268 * t308;
t305 = qJDD(1) * t270;
t307 = qJD(1) * qJD(4);
t311 = qJD(3) * t268;
t188 = -t270 * t307 - t233 * t311 + qJDD(3) * pkin(3) + t224 + (t301 - t305) * qJ(4);
t298 = -qJ(4) * qJD(1) + t233;
t310 = qJD(3) * t270;
t194 = t298 * t310 + (-qJ(4) * qJDD(1) + t232 - t307) * t268;
t176 = t188 * t266 + t194 * t267;
t259 = qJ(3) + pkin(7);
t252 = sin(t259);
t253 = cos(t259);
t332 = g(3) * t253 + t252 * t336 - t176;
t299 = t338 * t270;
t210 = -qJD(3) * t299 - qJD(4) * t268;
t283 = -qJD(4) * t270 + t311 * t338;
t193 = t267 * t210 + t266 * t283;
t228 = t338 * t268;
t203 = -t267 * t228 - t266 * t299;
t331 = qJD(3) * t193 + qJDD(3) * t203 + t253 * t295;
t192 = t210 * t266 - t267 * t283;
t303 = t266 * t305 + t267 * t343;
t199 = t266 * t301 - t303;
t292 = t266 * t306 - t267 * t305;
t200 = t292 + t340;
t202 = -t228 * t266 + t267 * t299;
t313 = qJD(1) * t268;
t302 = t266 * t313;
t312 = qJD(1) * t270;
t219 = t267 * t312 - t302;
t330 = t192 * t219 - t193 * t335 + t199 * t203 - t200 * t202;
t214 = t219 ^ 2;
t327 = pkin(3) * t268;
t326 = g(3) * t268;
t274 = qJD(1) ^ 2;
t325 = qJ(2) * t274;
t226 = pkin(3) * t313 + qJD(1) * qJ(2) + qJD(4);
t191 = pkin(4) * t335 - qJ(5) * t219 + t226;
t323 = t191 * t219;
t211 = t298 * t268;
t320 = t211 * t266;
t319 = t338 * t269;
t207 = t267 * t211;
t175 = t188 * t267 - t194 * t266;
t212 = -qJ(4) * t312 + t233 * t270;
t209 = qJD(3) * pkin(3) + t212;
t190 = t209 * t266 + t207;
t265 = t270 ^ 2;
t316 = t268 ^ 2 - t265;
t273 = qJD(3) ^ 2;
t315 = -t273 - t274;
t314 = qJD(1) * t226;
t196 = t212 * t267 - t320;
t309 = qJD(5) - t196;
t237 = pkin(3) * t310 + qJD(2);
t304 = qJDD(3) * t268;
t247 = qJ(2) + t327;
t189 = t209 * t267 - t320;
t206 = pkin(3) * t343 + qJDD(4) + t261 + t262;
t174 = -qJDD(3) * pkin(4) + qJDD(5) - t175;
t227 = pkin(4) * t267 + qJ(5) * t266 + pkin(3);
t229 = -pkin(4) * t266 + qJ(5) * t267;
t289 = t227 * t268 - t229 * t270 + qJ(2);
t288 = 0.2e1 * qJ(2) * t308 + qJDD(3) * t272;
t286 = -t336 - t325;
t195 = t212 * t266 + t207;
t245 = g(3) * t252;
t280 = qJD(3) * t195 - t253 * t336 + t175 + t245;
t279 = -pkin(4) * t199 + qJ(5) * t200 + t206;
t278 = -qJD(3) * t192 - qJDD(3) * t202 - t252 * t295;
t260 = qJDD(3) * qJ(5);
t173 = qJD(3) * qJD(5) + t176 + t260;
t183 = -qJD(3) * pkin(4) + qJD(5) - t189;
t184 = qJD(3) * qJ(5) + t190;
t217 = t266 * t311 - t267 * t310;
t218 = -t266 * t310 - t267 * t311;
t222 = -t266 * t268 + t267 * t270;
t277 = t173 * t223 - t174 * t222 - t183 * t218 - t184 * t217 - t336;
t276 = t175 * t222 + t176 * t223 + t189 * t218 - t190 * t217 - t336;
t275 = -t272 * t273 + t342;
t254 = qJDD(3) * t270;
t248 = -pkin(3) * t267 - pkin(4);
t244 = pkin(3) * t266 + qJ(5);
t243 = t338 * t271;
t198 = pkin(4) * t223 - qJ(5) * t222 + t247;
t197 = pkin(3) * t312 + pkin(4) * t219 + qJ(5) * t335;
t179 = -pkin(4) * t217 - qJ(5) * t218 - qJD(5) * t222 + t237;
t172 = -qJD(5) * t219 + t279;
t1 = [qJDD(1) * MDP(1) + t336 * MDP(2) + t295 * MDP(3) + (-0.2e1 * t324 + t344) * MDP(4) + t342 * MDP(5) + (t341 * pkin(1) + (t284 + t261) * qJ(2)) * MDP(6) + (qJDD(1) * t265 - 0.2e1 * t268 * t300) * MDP(7) + 0.2e1 * (-t268 * t305 + t308 * t316) * MDP(8) + (-t268 * t273 + t254) * MDP(9) + (-t270 * t273 - t304) * MDP(10) + (t268 * t275 + t270 * t288) * MDP(12) + (-t268 * t288 + t270 * t275) * MDP(13) + (-t199 * t247 + t206 * t223 - t217 * t226 + t237 * t335 + t278) * MDP(14) + (-t200 * t247 + t206 * t222 + t218 * t226 + t219 * t237 - t331) * MDP(15) + (-t276 + t330) * MDP(16) + (t176 * t203 + t190 * t193 - t175 * t202 - t189 * t192 + t206 * t247 + t226 * t237 - g(1) * (t247 * t271 - t319) - g(2) * (t247 * t269 + t243)) * MDP(17) + (t172 * t223 + t179 * t335 - t191 * t217 - t198 * t199 + t278) * MDP(18) + (-t277 + t330) * MDP(19) + (-t172 * t222 - t179 * t219 - t191 * t218 + t198 * t200 + t331) * MDP(20) + (t173 * t203 + t184 * t193 + t172 * t198 + t191 * t179 + t174 * t202 + t183 * t192 - g(1) * (t271 * t289 - t319) - g(2) * (t269 * t289 + t243)) * MDP(21); qJDD(1) * MDP(4) - t274 * MDP(5) + (-t325 - t341) * MDP(6) + (t268 * t315 + t254) * MDP(12) + (t270 * t315 - t304) * MDP(13) + (t276 - t314) * MDP(17) + (-qJD(1) * t191 + t277) * MDP(21) + t339 * (t199 * t223 + t200 * t222 + t217 * t335 - t218 * t219) + t334 * (-qJD(1) * t335 + qJD(3) * t218 + qJDD(3) * t222) + t333 * (qJD(1) * t219 - qJD(3) * t217 + qJDD(3) * t223); MDP(9) * t305 - MDP(10) * t306 + qJDD(3) * MDP(11) + (t270 * t286 + t224 + t326) * MDP(12) + (g(3) * t270 + (-t232 - t286) * t268) * MDP(13) + (-t219 * t226 + (qJDD(3) * t267 - t312 * t335) * pkin(3) + t280) * MDP(14) + (qJD(3) * t196 + t335 * t226 + (-qJDD(3) * t266 - t219 * t312) * pkin(3) + t332) * MDP(15) + ((t190 - t195) * t219 + (-t189 + t196) * t335 + (t199 * t266 + t200 * t267) * pkin(3)) * MDP(16) + (t189 * t195 - t190 * t196 + (t326 + t175 * t267 + t176 * t266 + (-t336 - t314) * t270) * pkin(3)) * MDP(17) + (-t323 - t197 * t335 - qJDD(5) + (pkin(4) - t248) * qJDD(3) + t280) * MDP(18) + (t199 * t244 - t200 * t248 + (t184 - t195) * t219 + (t183 - t309) * t335) * MDP(19) + (qJDD(3) * t244 - t191 * t335 + t197 * t219 + t260 + (0.2e1 * qJD(5) - t196) * qJD(3) - t332) * MDP(20) + (t173 * t244 + t174 * t248 - t191 * t197 - t183 * t195 - g(3) * (-pkin(4) * t252 + qJ(5) * t253 - t327) + t309 * t184 - t336 * (t227 * t270 + t229 * t268)) * MDP(21) + (t270 * t268 * MDP(7) - MDP(8) * t316) * t274; (t189 * t219 + t190 * t335 + t206 - t295) * MDP(17) + (t184 * t335 + (-qJD(5) - t183) * t219 + t279 - t295) * MDP(21) + t333 * (t292 + 0.2e1 * t340) + t339 * (-t335 ^ 2 - t214) + (t303 + qJD(3) * (t219 - t302)) * t334; (t219 * t335 - qJDD(3)) * MDP(18) - t292 * MDP(19) + (-t214 - t273) * MDP(20) + (t222 * t336 + t174 - t245 + t323) * MDP(21) + ((-t266 * t312 - t267 * t313 + t335) * MDP(19) - t184 * MDP(21)) * qJD(3);];
tau = t1;
