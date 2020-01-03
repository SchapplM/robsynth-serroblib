% Calculate vector of inverse dynamics joint torques for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:57
% EndTime: 2019-12-31 17:47:59
% DurationCPUTime: 2.10s
% Computational Cost: add. (895->266), mult. (1979->363), div. (0->0), fcn. (1371->8), ass. (0->134)
t342 = pkin(3) + qJ(2);
t263 = sin(pkin(7));
t259 = t263 ^ 2;
t265 = cos(pkin(7));
t261 = t265 ^ 2;
t270 = qJD(1) ^ 2;
t341 = (t259 + t261) * t270;
t262 = sin(pkin(8));
t264 = cos(pkin(8));
t296 = (-t262 ^ 2 - t264 ^ 2) * MDP(14);
t340 = MDP(9) + t296;
t269 = cos(qJ(1));
t257 = g(2) * t269;
t267 = sin(qJ(1));
t336 = g(1) * t267;
t298 = -t257 + t336;
t307 = qJDD(1) * t265;
t312 = qJD(1) * qJD(2);
t202 = t265 * t312 + t342 * t307 + qJDD(4);
t226 = t342 * t265;
t301 = t261 * t312;
t339 = (qJDD(1) * t226 + t202) * t265 + t301;
t338 = t264 * MDP(12) + MDP(8);
t337 = -MDP(12) * t262 - MDP(10);
t335 = qJDD(1) * pkin(1);
t318 = qJD(1) * t265;
t224 = t264 * t318 + qJD(5);
t334 = t224 * t264;
t333 = t261 * t270;
t332 = t263 * t264;
t266 = sin(qJ(5));
t331 = t263 * t266;
t268 = cos(qJ(5));
t330 = t263 * t268;
t329 = t263 * t269;
t328 = t263 * t270;
t327 = t264 * t265;
t326 = t265 * t266;
t325 = t265 * t267;
t324 = t265 * t268;
t323 = t265 * t269;
t300 = -qJ(3) * t263 - pkin(1);
t214 = (-pkin(2) - qJ(4)) * t265 + t300;
t222 = -qJD(3) * t263 - qJD(4) * t265;
t188 = t222 * qJD(1) + t214 * qJDD(1) + qJDD(2);
t309 = qJDD(1) * t263;
t215 = qJ(2) * t309 + t263 * t312 + qJDD(3);
t201 = pkin(3) * t309 + t215;
t176 = t264 * t188 + t262 * t201;
t198 = t214 * qJD(1) + qJD(2);
t319 = qJD(1) * t263;
t230 = qJ(2) * t319 + qJD(3);
t216 = pkin(3) * t319 + t230;
t180 = t264 * t198 + t262 * t216;
t225 = t342 * t263;
t190 = t264 * t214 + t262 * t225;
t237 = 0.2e1 * t301;
t306 = qJDD(1) * qJ(2) ^ 2;
t322 = qJ(2) * t237 + t261 * t306;
t321 = t269 * pkin(1) + t267 * qJ(2);
t317 = qJD(2) * t263;
t316 = qJD(2) * t265;
t315 = qJD(5) * t224;
t223 = t264 * t307 + qJDD(5);
t314 = t223 * MDP(20);
t313 = qJ(2) * qJDD(1);
t311 = qJD(1) * qJD(3);
t220 = -pkin(2) * t265 + t300;
t310 = qJDD(1) * t220;
t308 = qJDD(1) * t264;
t305 = t262 * t318;
t217 = t342 * t318 + qJD(4);
t304 = t263 * t311;
t303 = qJD(5) * t319;
t302 = t262 * t307;
t254 = t269 * qJ(2);
t299 = -pkin(1) * t267 + t254;
t174 = pkin(6) * t309 + t176;
t284 = (pkin(4) * t264 + pkin(6) * t262) * t265;
t185 = qJDD(1) * t284 + t202;
t297 = -t266 * t174 + t268 * t185;
t250 = qJDD(2) - t335;
t294 = pkin(2) * t323 + qJ(3) * t329 + t321;
t293 = t268 * t305;
t292 = qJD(5) * t305;
t291 = g(1) * t269 + g(2) * t267;
t290 = t268 * t174 + t266 * t185;
t178 = pkin(6) * t319 + t180;
t191 = qJD(1) * t284 + t217;
t289 = -t178 * t268 - t191 * t266;
t288 = t178 * t266 - t191 * t268;
t187 = pkin(6) * t263 + t190;
t195 = t284 + t226;
t287 = t187 * t268 + t195 * t266;
t286 = -t187 * t266 + t195 * t268;
t175 = -t188 * t262 + t201 * t264;
t179 = -t198 * t262 + t216 * t264;
t189 = -t214 * t262 + t225 * t264;
t285 = -t250 + t335 - t257;
t213 = -t267 * t263 * t262 + t264 * t269;
t283 = t213 * t266 + t267 * t324;
t282 = -t213 * t268 + t266 * t325;
t211 = t262 * t324 - t331;
t210 = t262 * t326 + t330;
t273 = qJDD(2) + t310;
t197 = t273 - t304;
t281 = -t197 - t310 - t257;
t280 = t291 * MDP(15);
t279 = -t223 * t266 - t268 * t315;
t278 = -t223 * t268 + t266 * t315;
t199 = t222 * t262 - t264 * t317;
t277 = -qJD(1) * t199 + qJDD(1) * t189 + t175;
t200 = t222 * t264 + t262 * t317;
t276 = -t200 * qJD(1) - t190 * qJDD(1) - t176;
t275 = 0.2e1 * t261 * t313 + t237 - t291;
t274 = (t312 + t313) * t259;
t183 = (-t302 + t303) * t268 + (t292 + t309) * t266;
t204 = t210 * qJD(1);
t207 = -t266 * t319 + t293;
t272 = (-t179 * t262 + t180 * t264) * MDP(15) + (-t204 * t262 - t266 * t334) * MDP(21) + (-t207 * t262 - t268 * t334) * MDP(22);
t255 = g(3) * t265;
t239 = g(1) * t325;
t228 = t266 * t303;
t212 = t262 * t329 + t264 * t267;
t209 = t220 * qJD(1) + qJD(2);
t206 = t211 * qJD(5);
t205 = t210 * qJD(5);
t193 = t212 * t268 + t266 * t323;
t192 = -t212 * t266 + t268 * t323;
t186 = -pkin(4) * t263 - t189;
t184 = t210 * qJDD(1) + t268 * t292 - t228;
t177 = -pkin(4) * t319 - t179;
t173 = -pkin(4) * t309 - t175;
t1 = [qJDD(1) * MDP(1) + t298 * MDP(2) + t291 * MDP(3) + (t285 * t265 + t239) * MDP(4) + (-t285 - t336) * t263 * MDP(5) + (0.2e1 * t274 + t275) * MDP(6) + (-t250 * pkin(1) - g(1) * t299 - g(2) * t321 + (0.2e1 * qJ(2) * t312 + t306) * t259 + t322) * MDP(7) + (t215 * t263 + t274 + t275) * MDP(8) + (-t239 + (-t281 - t304) * t265) * MDP(9) + (t259 * t311 + (t281 + t336) * t263) * MDP(10) + (t197 * t220 - g(1) * (-pkin(2) * t325 + t299) - g(2) * t294 + (qJ(2) * t215 + qJ(3) * t336 + qJD(2) * t230 - qJD(3) * t209) * t263 + t322) * MDP(11) + (-g(1) * t213 - g(2) * t212 + t277 * t263 + t339 * t264) * MDP(12) + ((t298 * t264 + t276) * t263 + (t291 - t339) * t262) * MDP(13) + (t239 + (t277 * t262 + t276 * t264 - t257) * t265) * MDP(14) + (t176 * t190 + t180 * t200 + t175 * t189 - t179 * t199 + t202 * t226 + t217 * t316 - g(1) * (pkin(3) * t269 + t254) - g(2) * (qJ(4) * t323 + t294) + (-g(2) * pkin(3) - g(1) * t214) * t267) * MDP(15) + (-t183 * t211 - t205 * t207) * MDP(16) + (t183 * t210 - t184 * t211 + t204 * t205 - t206 * t207) * MDP(17) + (t183 * t327 + t205 * t224 - t211 * t223) * MDP(18) + (t184 * t327 + t206 * t224 + t210 * t223) * MDP(19) + t314 * t327 + ((-t287 * qJD(5) - t200 * t266 + t268 * t316) * t224 + t286 * t223 + (t289 * qJD(5) + t297) * t327 - t199 * t204 - t186 * t184 - t173 * t210 - t177 * t206 + g(1) * t282 - g(2) * t193) * MDP(21) + (-(t286 * qJD(5) + t200 * t268 + t266 * t316) * t224 - t287 * t223 - (-t288 * qJD(5) + t290) * t327 - t199 * t207 + t186 * t183 - t173 * t211 + t177 * t205 + g(1) * t283 - g(2) * t192) * MDP(22); (t250 - t298) * MDP(7) + (-qJ(2) * t333 + (-qJD(3) - t230) * t319 + t273 - t298) * MDP(11) + (t262 * t341 - t263 * t308) * MDP(13) + (-t175 * t262 + t176 * t264 + (-t217 * t265 + (-t179 * t264 - t180 * t262) * t263) * qJD(1) - t298) * MDP(15) + (-t262 * t184 + t279 * t264 + (-(-t262 * t331 + t324) * t224 - t204 * t332) * qJD(1)) * MDP(21) + (t262 * t183 + t278 * t264 + ((t262 * t330 + t326) * t224 - t207 * t332) * qJD(1)) * MDP(22) + (MDP(5) + t337) * t309 - (qJ(2) * MDP(7) + MDP(6) + t338) * t341 + (-MDP(4) + t340) * t307; (t255 + t215) * MDP(11) + t255 * MDP(15) + (MDP(15) * t175 + MDP(21) * t184 - MDP(22) * t183) * t264 + (t176 * MDP(15) + t279 * MDP(21) + t278 * MDP(22)) * t262 + (-MDP(13) * t264 + t337) * t270 * t259 + (-t291 * MDP(11) - t280 + t340 * t270 * t265 + (-t262 * MDP(13) + t338) * qJDD(1) + (t209 * MDP(11) + t272) * qJD(1)) * t263; (-g(3) * t263 + t202) * MDP(15) - t278 * MDP(21) + t279 * MDP(22) + t296 * t333 + ((-t262 * t328 + t308) * MDP(12) + (-qJDD(1) * t262 - t264 * t328) * MDP(13) - t280 + t272 * qJD(1)) * t265; t207 * t204 * MDP(16) + (-t204 ^ 2 + t207 ^ 2) * MDP(17) + (-t204 * t224 + t183) * MDP(18) + (-t207 * t224 + t266 * t302 + t268 * t309 - t228) * MDP(19) + t314 + (-g(1) * t192 - g(2) * t283 - g(3) * t210 + t177 * t207 - t289 * t224 + t297) * MDP(21) + (g(1) * t193 + g(2) * t282 - g(3) * t211 - t177 * t204 - t288 * t224 - t290) * MDP(22) + (MDP(19) * t293 + t289 * MDP(21) + t288 * MDP(22)) * qJD(5);];
tau = t1;
