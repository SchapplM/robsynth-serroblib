% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:45
% EndTime: 2022-01-20 09:51:49
% DurationCPUTime: 1.31s
% Computational Cost: add. (1081->208), mult. (1728->271), div. (0->0), fcn. (1120->16), ass. (0->128)
t284 = cos(qJ(2));
t333 = pkin(1) * qJD(2);
t314 = qJD(1) * t333;
t281 = sin(qJ(2));
t317 = qJDD(1) * t281;
t346 = pkin(1) * t317 + t284 * t314;
t278 = cos(pkin(9));
t283 = cos(qJ(5));
t322 = t283 * t278;
t276 = sin(pkin(9));
t280 = sin(qJ(5));
t326 = t276 * t280;
t221 = -t322 + t326;
t222 = t276 * t283 + t278 * t280;
t274 = qJD(1) + qJD(2);
t208 = t222 * t274;
t319 = t276 ^ 2 + t278 ^ 2;
t345 = t319 * t274;
t340 = pkin(1) * t284;
t256 = qJDD(1) * t340;
t270 = qJDD(1) + qJDD(2);
t210 = pkin(2) * t270 - t281 * t314 + t256;
t277 = sin(pkin(8));
t279 = cos(pkin(8));
t189 = t277 * t210 + t346 * t279;
t182 = qJ(4) * t270 + qJD(4) * t274 + t189;
t259 = t278 * qJDD(3);
t176 = -t182 * t276 + t259;
t177 = t276 * qJDD(3) + t278 * t182;
t344 = -t176 * t276 + t177 * t278;
t275 = qJ(1) + qJ(2);
t265 = sin(t275);
t266 = cos(t275);
t343 = g(1) * t265 - g(2) * t266;
t334 = pkin(1) * qJD(1);
t316 = t281 * t334;
t238 = t277 * t316;
t315 = t284 * t334;
t215 = t279 * t315 - t238;
t318 = qJD(4) - t215;
t342 = pkin(1) * t281;
t282 = sin(qJ(1));
t341 = pkin(1) * t282;
t339 = pkin(2) * t265;
t338 = pkin(2) * t279;
t337 = pkin(4) * t278;
t264 = pkin(8) + t275;
t248 = sin(t264);
t336 = g(1) * t248;
t273 = pkin(9) + qJ(5);
t262 = sin(t273);
t331 = t248 * t262;
t263 = cos(t273);
t330 = t248 * t263;
t249 = cos(t264);
t329 = t249 * t262;
t328 = t249 * t263;
t327 = t270 * t278;
t323 = t279 * t281;
t225 = pkin(2) * t274 + t315;
t239 = t279 * t316;
t203 = t277 * t225 + t239;
t255 = pkin(2) + t340;
t321 = pkin(1) * t323 + t277 * t255;
t320 = g(1) * t266 + g(2) * t265;
t244 = t277 * t342;
t312 = t274 * t326;
t311 = t274 * t322;
t310 = qJD(5) * t311 + t222 * t270;
t254 = pkin(2) * t266;
t309 = t249 * pkin(3) + t248 * qJ(4) + t254;
t308 = -pkin(3) - t337;
t188 = t210 * t279 - t346 * t277;
t290 = qJDD(4) - t188;
t183 = -pkin(3) * t270 + t290;
t307 = -g(2) * t249 - t183;
t306 = t319 * t270;
t202 = t225 * t279 - t238;
t305 = t255 * t279 - t244;
t304 = qJD(1) * (-qJD(2) + t274);
t303 = qJD(2) * (-qJD(1) - t274);
t212 = -pkin(3) - t305;
t301 = t256 + t343;
t300 = t221 * t270;
t299 = qJD(4) - t202;
t186 = -qJD(5) * t312 + t310;
t217 = t222 * qJD(5);
t187 = t217 * t274 + t300;
t216 = t221 * qJD(5);
t195 = -qJD(5) * t216 + qJDD(5) * t222;
t196 = -qJD(5) * t217 - qJDD(5) * t221;
t206 = -t311 + t312;
t298 = (-t186 * t221 - t187 * t222 + t206 * t216 - t208 * t217) * MDP(12) + (t186 * t222 - t208 * t216) * MDP(11) + t195 * MDP(13) + t196 * MDP(14) + t270 * MDP(4);
t199 = qJ(4) * t274 + t203;
t190 = t278 * qJD(3) - t199 * t276;
t191 = t276 * qJD(3) + t278 * t199;
t297 = t190 * t276 - t191 * t278;
t211 = qJ(4) + t321;
t200 = (-pkin(7) - t211) * t276;
t267 = t278 * pkin(7);
t201 = t211 * t278 + t267;
t296 = t200 * t283 - t201 * t280;
t295 = t200 * t280 + t201 * t283;
t246 = pkin(2) * t277 + qJ(4);
t218 = (-pkin(7) - t246) * t276;
t219 = t246 * t278 + t267;
t294 = t218 * t283 - t219 * t280;
t293 = t218 * t280 + t219 * t283;
t292 = t279 * t284 * t333 - qJD(2) * t244;
t291 = -pkin(3) * t248 + t249 * qJ(4) - t339;
t178 = t270 * t308 + t290;
t192 = t274 * t308 + t299;
t289 = -g(1) * t331 + g(2) * t329 + t178 * t222 - t192 * t216;
t288 = g(1) * t330 - g(2) * t328 + t178 * t221 + t192 * t217;
t286 = -g(1) * t249 - g(2) * t248 + t344;
t285 = cos(qJ(1));
t268 = t285 * pkin(1);
t250 = -pkin(3) - t338;
t228 = t278 * t336;
t227 = t308 - t338;
t214 = (t277 * t284 + t323) * t333;
t213 = t277 * t315 + t239;
t209 = qJD(4) + t292;
t204 = t212 - t337;
t198 = -pkin(3) * t274 + t299;
t172 = pkin(7) * t327 + t177;
t171 = t259 + (-pkin(7) * t270 - t182) * t276;
t1 = [qJDD(1) * MDP(1) + (g(1) * t282 - g(2) * t285) * MDP(2) + (g(1) * t285 + g(2) * t282) * MDP(3) + ((t270 * t284 + t281 * t303) * pkin(1) + t301) * MDP(5) + (((-qJDD(1) - t270) * t281 + t284 * t303) * pkin(1) + t320) * MDP(6) + (t189 * t321 + t203 * t292 + t188 * t305 - t202 * t214 - g(1) * (-t339 - t341) - g(2) * (t254 + t268)) * MDP(7) + (t228 + (-t212 * t270 - t214 * t274 + t307) * t278) * MDP(8) + (t209 * t345 + t211 * t306 + t286) * MDP(9) + (t183 * t212 + t198 * t214 - g(1) * (t291 - t341) - g(2) * (t268 + t309) + t344 * t211 - t297 * t209) * MDP(10) + (t214 * t206 + t204 * t187 + t296 * qJDD(5) + (-qJD(5) * t295 - t209 * t222) * qJD(5) + t288) * MDP(16) + (t214 * t208 + t204 * t186 - t295 * qJDD(5) + (-qJD(5) * t296 + t209 * t221) * qJD(5) + t289) * MDP(17) + t298; (t304 * t342 + t301) * MDP(5) + ((t284 * t304 - t317) * pkin(1) + t320) * MDP(6) + (t202 * t213 - t203 * t215 + (t188 * t279 + t189 * t277 + t343) * pkin(2)) * MDP(7) + (t228 + (t213 * t274 - t250 * t270 + t307) * t278) * MDP(8) + (t246 * t306 + t318 * t345 + t286) * MDP(9) + (t183 * t250 - t198 * t213 - g(1) * t291 - g(2) * t309 + (t177 * t246 + t191 * t318) * t278 + (-t176 * t246 - t190 * t318) * t276) * MDP(10) + (t227 * t187 + t294 * qJDD(5) - t213 * t206 + (-qJD(5) * t293 - t318 * t222) * qJD(5) + t288) * MDP(16) + (t227 * t186 - t293 * qJDD(5) - t213 * t208 + (-qJD(5) * t294 + t318 * t221) * qJD(5) + t289) * MDP(17) + t298; (qJDD(3) - g(3)) * MDP(7) + (t176 * t278 + t177 * t276 - g(3)) * MDP(10) + t196 * MDP(16) - t195 * MDP(17); -MDP(8) * t327 + (t274 * t297 - t307 - t336) * MDP(10) + t300 * MDP(16) + t310 * MDP(17) - t319 * MDP(9) * t274 ^ 2 + (0.2e1 * t208 * MDP(16) + (-t206 - t312) * MDP(17)) * qJD(5); t208 * t206 * MDP(11) + (-t206 ^ 2 + t208 ^ 2) * MDP(12) - t300 * MDP(14) + qJDD(5) * MDP(15) + (g(1) * t329 + g(2) * t331 - g(3) * t263 + t283 * t171 - t280 * t172 - t192 * t208) * MDP(16) + (g(1) * t328 + g(2) * t330 + g(3) * t262 - t280 * t171 - t283 * t172 + t192 * t206) * MDP(17) + (t310 + (t206 - t312) * qJD(5)) * MDP(13);];
tau = t1;
