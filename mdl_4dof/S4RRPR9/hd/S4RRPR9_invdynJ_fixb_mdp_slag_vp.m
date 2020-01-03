% Calculate vector of inverse dynamics joint torques for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:10
% EndTime: 2019-12-31 17:10:14
% DurationCPUTime: 3.25s
% Computational Cost: add. (1289->306), mult. (2986->426), div. (0->0), fcn. (2043->10), ass. (0->138)
t281 = sin(qJ(2));
t323 = qJDD(1) * t281;
t267 = pkin(5) * t323;
t284 = cos(qJ(2));
t324 = qJD(1) * qJD(2);
t313 = t284 * t324;
t231 = -qJDD(2) * pkin(2) + pkin(5) * t313 + qJDD(3) + t267;
t282 = sin(qJ(1));
t285 = cos(qJ(1));
t310 = g(1) * t285 + g(2) * t282;
t346 = g(3) * t284;
t292 = t281 * t310 - t346;
t352 = t231 - t292;
t332 = qJD(1) * t284;
t264 = -qJD(4) + t332;
t278 = sin(pkin(7));
t333 = qJD(1) * t281;
t316 = t278 * t333;
t279 = cos(pkin(7));
t326 = t279 * qJD(2);
t244 = t316 - t326;
t315 = t279 * t333;
t331 = qJD(2) * t278;
t246 = t315 + t331;
t280 = sin(qJ(4));
t283 = cos(qJ(4));
t304 = t244 * t280 - t246 * t283;
t353 = t264 * t304;
t309 = g(1) * t282 - g(2) * t285;
t351 = pkin(5) * t244;
t350 = pkin(5) * t246;
t347 = g(3) * t281;
t345 = pkin(6) + qJ(3);
t343 = t246 * t280;
t197 = t283 * t244 + t343;
t344 = t197 * t264;
t342 = t278 * t280;
t341 = t278 * t284;
t340 = t279 * t281;
t339 = t279 * t284;
t338 = t282 * t284;
t337 = t284 * t285;
t308 = pkin(2) * t281 - qJ(3) * t284;
t233 = qJD(2) * t308 - qJD(3) * t281;
t302 = pkin(2) * t284 + qJ(3) * t281 + pkin(1);
t196 = qJD(1) * t233 - qJDD(1) * t302;
t321 = t284 * qJDD(1);
t294 = t281 * t324 - t321;
t222 = -pkin(5) * t294 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t185 = t278 * t196 + t279 * t222;
t249 = -t283 * t279 + t342;
t296 = t249 * t284;
t336 = qJD(1) * t296 - t249 * qJD(4);
t250 = t278 * t283 + t279 * t280;
t297 = t250 * t284;
t335 = -qJD(1) * t297 + t250 * qJD(4);
t239 = t302 * qJD(1);
t269 = pkin(5) * t332;
t259 = qJD(2) * qJ(3) + t269;
t203 = -t278 * t239 + t279 * t259;
t252 = t308 * qJD(1);
t211 = pkin(5) * t316 + t279 * t252;
t330 = qJD(2) * t281;
t319 = pkin(5) * t330;
t204 = t279 * t233 + t278 * t319;
t217 = pkin(5) * t339 - t278 * t302;
t276 = t281 ^ 2;
t334 = -t284 ^ 2 + t276;
t329 = qJD(2) * t284;
t328 = qJD(4) * t281;
t327 = qJD(4) * t283;
t254 = -qJD(2) * pkin(2) + pkin(5) * t333 + qJD(3);
t325 = qJD(3) - t254;
t322 = qJDD(2) * t278;
t320 = pkin(3) * t332;
t270 = t279 * qJDD(2);
t295 = t313 + t323;
t214 = t278 * t295 - t270;
t215 = t279 * t295 + t322;
t318 = -t280 * t214 + t283 * t215 - t244 * t327;
t317 = pkin(3) * t278 + pkin(5);
t314 = qJ(3) * t321;
t184 = t279 * t196 - t222 * t278;
t180 = pkin(3) * t294 - pkin(6) * t215 + t184;
t183 = -pkin(6) * t214 + t185;
t312 = t283 * t180 - t183 * t280;
t311 = t283 * t214 + t280 * t215;
t202 = -t279 * t239 - t259 * t278;
t307 = t180 * t280 + t183 * t283;
t186 = -pkin(6) * t246 + t202 - t320;
t187 = -pkin(6) * t244 + t203;
t177 = t186 * t283 - t187 * t280;
t178 = t186 * t280 + t187 * t283;
t243 = t279 * t302;
t201 = -pkin(6) * t340 - t243 + (-pkin(5) * t278 - pkin(3)) * t284;
t206 = -pkin(6) * t278 * t281 + t217;
t306 = t201 * t283 - t206 * t280;
t305 = t201 * t280 + t206 * t283;
t303 = pkin(3) * t281 - pkin(6) * t339;
t301 = -0.2e1 * pkin(1) * t324 - pkin(5) * qJDD(2);
t258 = t345 * t279;
t300 = qJD(1) * t303 + qJD(3) * t278 + qJD(4) * t258 + t211;
t237 = t278 * t252;
t257 = t345 * t278;
t298 = -pkin(5) * t340 - pkin(6) * t341;
t299 = -qJD(1) * t298 + qJD(3) * t279 - qJD(4) * t257 - t237;
t181 = -qJD(4) * t343 + t318;
t287 = qJD(1) ^ 2;
t293 = pkin(1) * t287 + t310;
t290 = -t284 * t310 - t347;
t286 = qJD(2) ^ 2;
t289 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t286 + t309;
t182 = -qJD(4) * t304 + t311;
t275 = pkin(7) + qJ(4);
t273 = cos(t275);
t272 = sin(t275);
t266 = -pkin(3) * t279 - pkin(2);
t253 = t317 * t281;
t248 = qJDD(4) + t294;
t241 = t317 * t329;
t240 = t278 * t320 + t269;
t230 = t249 * t281;
t229 = t250 * t281;
t226 = t272 * t282 + t273 * t337;
t225 = -t272 * t337 + t273 * t282;
t224 = t272 * t285 - t273 * t338;
t223 = t272 * t338 + t273 * t285;
t220 = t278 * t233;
t216 = -pkin(5) * t341 - t243;
t212 = -pkin(5) * t315 + t237;
t210 = pkin(3) * t244 + t254;
t205 = -t279 * t319 + t220;
t195 = qJD(2) * t298 + t220;
t191 = qJD(2) * t297 + t327 * t340 - t328 * t342;
t190 = -qJD(2) * t296 - t250 * t328;
t189 = pkin(3) * t214 + t231;
t188 = qJD(2) * t303 + t204;
t1 = [qJDD(1) * MDP(1) + t309 * MDP(2) + t310 * MDP(3) + (qJDD(1) * t276 + 0.2e1 * t281 * t313) * MDP(4) + 0.2e1 * (t281 * t321 - t324 * t334) * MDP(5) + (qJDD(2) * t281 + t284 * t286) * MDP(6) + (qJDD(2) * t284 - t281 * t286) * MDP(7) + (t281 * t301 + t284 * t289) * MDP(9) + (-t281 * t289 + t284 * t301) * MDP(10) + (-t310 * t278 + (pkin(5) * t214 + t231 * t278 + (qJD(1) * t216 + t202) * qJD(2)) * t281 + (-t204 * qJD(1) - t216 * qJDD(1) - t184 + t309 * t279 + (t254 * t278 + t351) * qJD(2)) * t284) * MDP(11) + (-t310 * t279 + (pkin(5) * t215 + t231 * t279 + (-qJD(1) * t217 - t203) * qJD(2)) * t281 + (t205 * qJD(1) + t217 * qJDD(1) + t185 - t309 * t278 + (t254 * t279 + t350) * qJD(2)) * t284) * MDP(12) + (-t204 * t246 - t205 * t244 - t214 * t217 - t215 * t216 + (-t202 * t279 - t203 * t278) * t329 + (-t184 * t279 - t185 * t278 + t309) * t281) * MDP(13) + (t184 * t216 + t185 * t217 + t202 * t204 + t203 * t205 + (t231 * t281 + t254 * t329 - t310) * pkin(5) + t309 * t302) * MDP(14) + (-t181 * t230 - t190 * t304) * MDP(15) + (-t181 * t229 + t182 * t230 - t190 * t197 + t191 * t304) * MDP(16) + (-t181 * t284 - t190 * t264 - t230 * t248 - t304 * t330) * MDP(17) + (t182 * t284 + t191 * t264 - t197 * t330 - t229 * t248) * MDP(18) + (-t248 * t284 - t264 * t330) * MDP(19) + (-(t188 * t283 - t195 * t280) * t264 + t306 * t248 - t312 * t284 + t177 * t330 + t241 * t197 + t253 * t182 + t189 * t229 + t210 * t191 - g(1) * t224 - g(2) * t226 + (t178 * t284 + t264 * t305) * qJD(4)) * MDP(20) + ((t188 * t280 + t195 * t283) * t264 - t305 * t248 + t307 * t284 - t178 * t330 - t241 * t304 + t253 * t181 - t189 * t230 + t210 * t190 - g(1) * t223 - g(2) * t225 + (t177 * t284 + t264 * t306) * qJD(4)) * MDP(21); MDP(6) * t323 + MDP(7) * t321 + qJDD(2) * MDP(8) + (t281 * t293 - t267 - t346) * MDP(9) + (t347 + (-pkin(5) * qJDD(1) + t293) * t284) * MDP(10) + (t278 * t314 - pkin(2) * t214 - t352 * t279 + ((-qJ(3) * t331 - t202) * t281 + (t278 * t325 + t211 - t351) * t284) * qJD(1)) * MDP(11) + (t279 * t314 - pkin(2) * t215 + t352 * t278 + ((-qJ(3) * t326 + t203) * t281 + (t279 * t325 - t212 - t350) * t284) * qJD(1)) * MDP(12) + (t211 * t246 + t212 * t244 + (-qJ(3) * t214 - qJD(3) * t244 + t202 * t332 + t185) * t279 + (qJ(3) * t215 + qJD(3) * t246 + t203 * t332 - t184) * t278 + t290) * MDP(13) + (-t254 * t269 - t202 * t211 - t203 * t212 + (-t202 * t278 + t203 * t279) * qJD(3) - t352 * pkin(2) + (-t184 * t278 + t185 * t279 + t290) * qJ(3)) * MDP(14) + (t181 * t250 - t304 * t336) * MDP(15) + (-t181 * t249 - t182 * t250 - t197 * t336 + t304 * t335) * MDP(16) + (t248 * t250 - t264 * t336 + t304 * t333) * MDP(17) + (t197 * t333 - t248 * t249 + t264 * t335) * MDP(18) + t264 * MDP(19) * t333 + ((-t257 * t283 - t258 * t280) * t248 + t266 * t182 + t189 * t249 - t177 * t333 - t240 * t197 + (t280 * t299 + t283 * t300) * t264 + t335 * t210 + t292 * t273) * MDP(20) + (-(-t257 * t280 + t258 * t283) * t248 + t266 * t181 + t189 * t250 + t178 * t333 + t240 * t304 + (-t280 * t300 + t283 * t299) * t264 + t336 * t210 - t292 * t272) * MDP(21) + (-MDP(4) * t281 * t284 + MDP(5) * t334) * t287; (t278 * t323 - t270 + (-t246 + t331) * t332) * MDP(11) + (t279 * t323 + t322 + (t244 + t326) * t332) * MDP(12) + (-t244 ^ 2 - t246 ^ 2) * MDP(13) + (t202 * t246 + t203 * t244 + t352) * MDP(14) + (t182 + t353) * MDP(20) + (t181 + t344) * MDP(21); -t304 * t197 * MDP(15) + (-t197 ^ 2 + t304 ^ 2) * MDP(16) + (t318 - t344) * MDP(17) + (-t311 + t353) * MDP(18) + t248 * MDP(19) + (-g(1) * t225 + g(2) * t223 - t178 * t264 + t210 * t304 + t272 * t347 + t312) * MDP(20) + (g(1) * t226 - g(2) * t224 - t177 * t264 + t197 * t210 + t273 * t347 - t307) * MDP(21) + (-MDP(17) * t343 + MDP(18) * t304 - MDP(20) * t178 - MDP(21) * t177) * qJD(4);];
tau = t1;
