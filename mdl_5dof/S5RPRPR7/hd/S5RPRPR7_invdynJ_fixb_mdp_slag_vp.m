% Calculate vector of inverse dynamics joint torques for
% S5RPRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:57
% EndTime: 2019-12-31 18:20:00
% DurationCPUTime: 2.11s
% Computational Cost: add. (1653->277), mult. (3580->376), div. (0->0), fcn. (2478->14), ass. (0->139)
t267 = sin(pkin(8));
t248 = pkin(1) * t267 + pkin(6);
t324 = qJ(4) + t248;
t266 = sin(pkin(9));
t268 = cos(pkin(9));
t272 = sin(qJ(3));
t275 = cos(qJ(3));
t235 = t266 * t275 + t268 * t272;
t317 = qJD(1) * qJD(3);
t307 = t275 * t317;
t308 = t272 * t317;
t209 = qJDD(1) * t235 - t266 * t308 + t268 * t307;
t345 = qJD(3) * qJD(5) + t209;
t252 = pkin(3) * t275 + pkin(2);
t269 = cos(pkin(8));
t341 = pkin(1) * t269;
t292 = -t252 - t341;
t344 = pkin(3) * t308 + t292 * qJDD(1) + qJDD(4);
t259 = t275 * qJDD(2);
t238 = t248 * qJDD(1);
t281 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t238;
t302 = t324 * qJD(1);
t291 = t302 * qJD(3);
t189 = qJDD(3) * pkin(3) - t272 * t281 - t275 * t291 + t259;
t192 = (qJDD(2) - t291) * t272 + t281 * t275;
t176 = t189 * t268 - t192 * t266;
t174 = -qJDD(3) * pkin(4) - t176;
t218 = t275 * qJD(2) - t272 * t302;
t338 = qJD(3) * pkin(3);
t213 = t218 + t338;
t219 = qJD(2) * t272 + t275 * t302;
t331 = t219 * t266;
t190 = t213 * t268 - t331;
t186 = -qJD(3) * pkin(4) - t190;
t304 = qJD(3) * t324;
t220 = qJD(4) * t275 - t272 * t304;
t285 = -qJD(4) * t272 - t275 * t304;
t196 = t268 * t220 + t266 * t285;
t326 = t268 * t275;
t234 = t266 * t272 - t326;
t201 = pkin(4) * t234 - pkin(7) * t235 + t292;
t233 = t324 * t275;
t305 = t324 * t272;
t205 = t268 * t233 - t266 * t305;
t229 = t235 * qJD(3);
t314 = qJDD(1) * t275;
t315 = qJDD(1) * t272;
t294 = -t266 * t315 + t268 * t314;
t207 = qJD(1) * t229 + qJDD(5) - t294;
t320 = qJD(1) * t272;
t228 = qJD(1) * t326 - t266 * t320;
t226 = qJD(5) - t228;
t232 = t234 * qJD(3);
t177 = t266 * t189 + t268 * t192;
t175 = qJDD(3) * pkin(7) + t177;
t227 = qJD(1) * t292 + qJD(4);
t230 = t235 * qJD(1);
t197 = -pkin(4) * t228 - pkin(7) * t230 + t227;
t301 = qJD(5) * t197 + t175;
t343 = t174 * t235 - t186 * t232 - t205 * t207 - (qJD(5) * t201 + t196) * t226 - t234 * t301;
t247 = pkin(3) * t266 + pkin(7);
t262 = qJ(3) + pkin(9);
t253 = sin(t262);
t255 = cos(t262);
t263 = qJ(1) + pkin(8);
t254 = sin(t263);
t256 = cos(t263);
t297 = g(1) * t256 + g(2) * t254;
t342 = t226 * (pkin(3) * t320 + pkin(4) * t230 - pkin(7) * t228 + qJD(5) * t247) - t253 * t297 + g(3) * t255 + t174;
t340 = g(3) * t253;
t339 = g(3) * t275;
t271 = sin(qJ(5));
t274 = cos(qJ(5));
t309 = t271 * qJDD(3) + t345 * t274;
t318 = qJD(5) * t271;
t183 = -t230 * t318 + t309;
t337 = t183 * t271;
t336 = t201 * t207;
t215 = -t274 * qJD(3) + t230 * t271;
t335 = t215 * t226;
t334 = t215 * t230;
t217 = qJD(3) * t271 + t230 * t274;
t333 = t217 * t226;
t332 = t217 * t230;
t330 = t254 * t271;
t329 = t254 * t274;
t328 = t256 * t271;
t327 = t256 * t274;
t211 = t268 * t219;
t325 = t271 * t207;
t202 = t274 * t207;
t323 = qJDD(2) - g(3);
t322 = t183 * t234 + t217 * t229;
t191 = t266 * t213 + t211;
t264 = t272 ^ 2;
t321 = -t275 ^ 2 + t264;
t250 = -pkin(2) - t341;
t241 = qJD(1) * t250;
t319 = qJD(5) * t235;
t312 = t272 * t338;
t311 = t235 * t325;
t310 = t235 * t202;
t303 = t226 * t274;
t208 = -qJD(3) * t230 + t294;
t181 = -pkin(4) * t208 - pkin(7) * t209 + t344;
t187 = qJD(3) * pkin(7) + t191;
t300 = qJD(5) * t187 - t181;
t296 = g(1) * t254 - g(2) * t256;
t273 = sin(qJ(1));
t276 = cos(qJ(1));
t295 = g(1) * t273 - g(2) * t276;
t258 = t274 * qJDD(3);
t184 = t217 * qJD(5) + t209 * t271 - t258;
t293 = -t184 * t234 - t215 * t229;
t290 = t202 + (t228 * t271 - t318) * t226;
t289 = t232 * t271 - t274 * t319;
t288 = t232 * t274 + t235 * t318;
t284 = -qJD(1) * t241 - t238 + t297;
t283 = 0.2e1 * qJD(3) * t241 - qJDD(3) * t248;
t194 = t218 * t268 - t331;
t282 = -t247 * t207 + (t186 + t194) * t226;
t277 = qJD(3) ^ 2;
t280 = -0.2e1 * qJDD(1) * t250 - t248 * t277 + t296;
t270 = -qJ(4) - pkin(6);
t249 = -pkin(3) * t268 - pkin(4);
t237 = qJDD(3) * t275 - t272 * t277;
t236 = qJDD(3) * t272 + t275 * t277;
t224 = t255 * t327 + t330;
t223 = -t255 * t328 + t329;
t222 = -t255 * t329 + t328;
t221 = t255 * t330 + t327;
t204 = t233 * t266 + t268 * t305;
t199 = pkin(4) * t229 + pkin(7) * t232 + t312;
t195 = t220 * t266 - t268 * t285;
t193 = t218 * t266 + t211;
t180 = t274 * t181;
t179 = t187 * t274 + t197 * t271;
t178 = -t187 * t271 + t197 * t274;
t1 = [qJDD(1) * MDP(1) + t295 * MDP(2) + (g(1) * t276 + g(2) * t273) * MDP(3) + (t295 + (t267 ^ 2 + t269 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t264 + 0.2e1 * t272 * t307) * MDP(5) + 0.2e1 * (t272 * t314 - t317 * t321) * MDP(6) + t236 * MDP(7) + t237 * MDP(8) + (t272 * t283 + t275 * t280) * MDP(10) + (-t272 * t280 + t275 * t283) * MDP(11) + (-t176 * t235 - t177 * t234 + t190 * t232 - t191 * t229 + t195 * t230 + t196 * t228 + t204 * t209 + t205 * t208 - t297) * MDP(12) + (t177 * t205 + t191 * t196 - t176 * t204 - t190 * t195 + t227 * t312 - g(1) * (-pkin(1) * t273 - t252 * t254 - t256 * t270) - g(2) * (pkin(1) * t276 + t252 * t256 - t254 * t270) + t344 * t292) * MDP(13) + (t183 * t235 * t274 - t217 * t288) * MDP(14) + (-(-t215 * t274 - t217 * t271) * t232 + (-t337 - t184 * t274 + (t215 * t271 - t217 * t274) * qJD(5)) * t235) * MDP(15) + (-t226 * t288 + t310 + t322) * MDP(16) + (t226 * t289 + t293 - t311) * MDP(17) + (t207 * t234 + t226 * t229) * MDP(18) + (-g(1) * t222 - g(2) * t224 + t178 * t229 + t180 * t234 + t204 * t184 + t195 * t215 + (t199 * t226 + t336 + (t186 * t235 - t187 * t234 - t205 * t226) * qJD(5)) * t274 + t343 * t271) * MDP(19) + (-g(1) * t221 - g(2) * t223 - t179 * t229 + t204 * t183 + t195 * t217 + (-(-qJD(5) * t205 + t199) * t226 - t336 + t300 * t234 - t186 * t319) * t271 + t343 * t274) * MDP(20); t323 * MDP(4) + t237 * MDP(10) - t236 * MDP(11) + (t208 * t235 + t209 * t234 - t228 * t232 + t229 * t230) * MDP(12) + (-t176 * t234 + t177 * t235 - t190 * t229 - t191 * t232 - g(3)) * MDP(13) + (-t293 - t311) * MDP(19) + (-t310 + t322) * MDP(20) + (MDP(19) * t289 + MDP(20) * t288) * t226; MDP(7) * t315 + MDP(8) * t314 + qJDD(3) * MDP(9) + (t272 * t284 + t259 - t339) * MDP(10) + (-t272 * t323 + t275 * t284) * MDP(11) + ((t191 - t193) * t230 + (t190 - t194) * t228 + (t208 * t266 - t209 * t268) * pkin(3)) * MDP(12) + (t190 * t193 - t191 * t194 + (-t339 + t176 * t268 + t177 * t266 + (-qJD(1) * t227 + t297) * t272) * pkin(3)) * MDP(13) + (t217 * t303 + t337) * MDP(14) + ((t183 - t335) * t274 + (-t184 - t333) * t271) * MDP(15) + (t226 * t303 + t325 - t332) * MDP(16) + (t290 + t334) * MDP(17) - t226 * t230 * MDP(18) + (-t178 * t230 + t249 * t184 - t193 * t215 + t282 * t271 - t342 * t274) * MDP(19) + (t179 * t230 + t249 * t183 - t193 * t217 + t342 * t271 + t282 * t274) * MDP(20) + (-t272 * t275 * MDP(5) + t321 * MDP(6)) * qJD(1) ^ 2; (-t228 ^ 2 - t230 ^ 2) * MDP(12) + (t290 - t334) * MDP(19) + (-t226 ^ 2 * t274 - t325 - t332) * MDP(20) + (t190 * t230 - t191 * t228 - t296 + t344) * MDP(13); t217 * t215 * MDP(14) + (-t215 ^ 2 + t217 ^ 2) * MDP(15) + (t309 + t335) * MDP(16) + (t258 + t333) * MDP(17) + t207 * MDP(18) + (-g(1) * t223 + g(2) * t221 + t179 * t226 - t186 * t217 + t180) * MDP(19) + (g(1) * t224 - g(2) * t222 + t178 * t226 + t186 * t215) * MDP(20) + ((-t175 + t340) * MDP(20) + (-MDP(17) * t230 - MDP(19) * t187 - MDP(20) * t197) * qJD(5)) * t274 + (-qJD(5) * t230 * MDP(16) - t345 * MDP(17) + (-t301 + t340) * MDP(19) + t300 * MDP(20)) * t271;];
tau = t1;
