% Calculate vector of inverse dynamics joint torques for
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:22
% EndTime: 2019-12-31 17:23:24
% DurationCPUTime: 1.29s
% Computational Cost: add. (968->201), mult. (1463->280), div. (0->0), fcn. (934->12), ass. (0->119)
t256 = sin(qJ(4));
t257 = sin(qJ(3));
t260 = cos(qJ(4));
t261 = cos(qJ(3));
t203 = t256 * t261 + t257 * t260;
t251 = qJD(1) + qJD(2);
t194 = t203 * t251;
t258 = sin(qJ(2));
t329 = pkin(1) * t258;
t295 = qJD(2) * t329;
t262 = cos(qJ(2));
t328 = pkin(1) * t262;
t307 = qJD(1) * t295 - qJDD(1) * t328;
t249 = qJDD(1) + qJDD(2);
t327 = pkin(2) * t249;
t196 = t307 - t327;
t255 = qJ(1) + qJ(2);
t245 = cos(t255);
t325 = g(2) * t245;
t331 = t196 + t325;
t309 = t260 * t261;
t313 = t256 * t257;
t202 = -t309 + t313;
t243 = sin(t255);
t306 = g(1) * t245 + g(2) * t243;
t250 = qJD(3) + qJD(4);
t330 = -pkin(6) - pkin(7);
t326 = pkin(2) * t251;
t233 = g(1) * t243;
t237 = pkin(6) + t329;
t324 = -pkin(7) - t237;
t323 = pkin(1) * qJD(1);
t297 = t258 * t323;
t209 = pkin(6) * t251 + t297;
t287 = pkin(7) * t251 + t209;
t189 = t287 * t261;
t322 = t189 * t260;
t254 = qJ(3) + qJ(4);
t242 = sin(t254);
t321 = t242 * t245;
t320 = t243 * t242;
t244 = cos(t254);
t319 = t243 * t244;
t318 = t245 * t244;
t317 = t249 * t257;
t316 = t249 * t261;
t315 = t250 * t262;
t314 = t251 * t257;
t310 = t257 * t261;
t304 = qJD(1) * t262;
t296 = pkin(1) * t304;
t210 = -t296 - t326;
t301 = qJD(3) * t257;
t308 = t210 * t301 + t261 * t233;
t252 = t257 ^ 2;
t305 = -t261 ^ 2 + t252;
t303 = qJD(2) * t262;
t302 = qJD(3) * t251;
t300 = qJD(3) * t261;
t299 = qJD(4) * t256;
t298 = qJDD(1) * t258;
t294 = pkin(1) * t303;
t293 = pkin(3) * t301;
t292 = t251 * t313;
t291 = t251 * t309;
t290 = t210 * t300 + t257 * t331;
t239 = -pkin(3) * t261 - pkin(2);
t289 = qJD(3) * t330;
t288 = t251 * t300;
t285 = qJD(3) * t324;
t284 = t251 * t297;
t283 = t202 * t249;
t188 = t287 * t257;
t185 = qJD(3) * pkin(3) - t188;
t282 = -t185 * t256 - t322;
t200 = t324 * t257;
t246 = t261 * pkin(7);
t201 = t237 * t261 + t246;
t281 = t200 * t260 - t201 * t256;
t280 = t200 * t256 + t201 * t260;
t226 = t330 * t257;
t227 = pkin(6) * t261 + t246;
t279 = t226 * t260 - t227 * t256;
t278 = t226 * t256 + t227 * t260;
t238 = -pkin(2) - t328;
t264 = qJD(3) ^ 2;
t277 = t237 * t264 + t238 * t249;
t276 = t233 - t307 - t325;
t275 = t251 * t301 - t316;
t180 = pkin(3) * t275 + t196;
t182 = t250 * t202;
t195 = t239 * t251 - t296;
t274 = -g(1) * t320 + g(2) * t321 + t180 * t203 - t182 * t195;
t183 = t250 * t203;
t273 = g(1) * t319 - g(2) * t318 + t180 * t202 + t183 * t195;
t272 = -qJDD(3) * t237 + t238 * t302;
t168 = qJD(4) * t291 + t203 * t249 - t250 * t292 + t260 * t288;
t192 = -t291 + t292;
t248 = qJDD(3) + qJDD(4);
t271 = t194 * t192 * MDP(14) + (t192 * t250 + t168) * MDP(16) - t283 * MDP(17) + (-t192 ^ 2 + t194 ^ 2) * MDP(15) + t248 * MDP(18);
t197 = pkin(6) * t249 + (qJD(1) * t303 + t298) * pkin(1);
t270 = -t210 * t251 - t197 + t306;
t169 = t183 * t251 + t283;
t269 = (-t168 * t202 - t169 * t203 + t182 * t192 - t183 * t194) * MDP(15) + (t168 * t203 - t182 * t194) * MDP(14) + (-t182 * t250 + t203 * t248) * MDP(16) + (-t183 * t250 - t202 * t248) * MDP(17) + 0.2e1 * (t249 * t310 - t302 * t305) * MDP(8) + (t249 * t252 + 0.2e1 * t257 * t288) * MDP(7) + (qJDD(3) * t261 - t257 * t264) * MDP(10) + (qJDD(3) * t257 + t261 * t264) * MDP(9) + t249 * MDP(4);
t268 = pkin(6) * t264 - t284 - t327;
t267 = -pkin(6) * qJDD(3) + (t296 - t326) * qJD(3);
t172 = -t209 * t300 + qJDD(3) * pkin(3) - t197 * t257 + (-t288 - t317) * pkin(7);
t266 = t195 * t192 + t189 * t299 + g(2) * t319 + g(1) * t318 + g(3) * t242 + (-t189 * t250 - t172) * t256;
t173 = -pkin(7) * t275 + t197 * t261 - t209 * t301;
t265 = g(1) * t321 + g(2) * t320 - g(3) * t244 + qJD(4) * t282 + t172 * t260 - t256 * t173 - t195 * t194;
t263 = cos(qJ(1));
t259 = sin(qJ(1));
t219 = t239 - t328;
t208 = t293 + t295;
t207 = t261 * t289;
t206 = t257 * t289;
t187 = -t257 * t294 + t261 * t285;
t186 = t257 * t285 + t261 * t294;
t1 = [((-t277 - t331) * MDP(12) + t272 * MDP(13)) * t261 + (t272 * MDP(12) + (t277 - t233) * MDP(13)) * t257 + t269 + (g(1) * t259 - g(2) * t263) * MDP(2) + (g(1) * t263 + g(2) * t259) * MDP(3) + qJDD(1) * MDP(1) + t308 * MDP(12) + t290 * MDP(13) + t276 * MDP(5) + t306 * MDP(6) + (t208 * t192 + t219 * t169 + (-qJD(4) * t280 - t186 * t256 + t187 * t260) * t250 + t281 * t248 + t273) * MDP(19) + (t208 * t194 + t219 * t168 - (qJD(4) * t281 + t186 * t260 + t187 * t256) * t250 - t280 * t248 + t274) * MDP(20) + (t249 * t262 * MDP(5) + (-qJDD(1) - t249) * MDP(6) * t258 + ((-MDP(12) * t261 + MDP(13) * t257 - MDP(5)) * t258 * t251 + ((-qJD(1) - t251) * MDP(6) + (-MDP(12) * t257 - MDP(13) * t261) * qJD(3)) * t262) * qJD(2)) * pkin(1); t269 + (t276 + t284) * MDP(5) + (t267 * t257 + (-t268 - t331) * t261 + t308) * MDP(12) + (t267 * t261 + (t268 - t233) * t257 + t290) * MDP(13) + ((-t298 + (-qJD(2) + t251) * t304) * pkin(1) + t306) * MDP(6) + (t192 * t293 + t239 * t169 + (-qJD(4) * t278 - t206 * t256 + t207 * t260) * t250 + t279 * t248 + (-t258 * t192 + t203 * t315) * t323 + t273) * MDP(19) + (t194 * t293 + t239 * t168 - (qJD(4) * t279 + t206 * t260 + t207 * t256) * t250 - t278 * t248 + (-t258 * t194 - t202 * t315) * t323 + t274) * MDP(20); MDP(9) * t317 + MDP(10) * t316 + qJDD(3) * MDP(11) + (-g(3) * t261 + t257 * t270) * MDP(12) + (g(3) * t257 + t261 * t270) * MDP(13) + (-(t188 * t256 - t322) * t250 + (-t192 * t314 + t248 * t260 - t250 * t299) * pkin(3) + t265) * MDP(19) + ((-qJD(4) * t185 - t188 * t250 - t173) * t260 + (-qJD(4) * t250 * t260 - t194 * t314 - t248 * t256) * pkin(3) + t266) * MDP(20) + t271 + (-MDP(7) * t310 + MDP(8) * t305) * t251 ^ 2; (-t250 * t282 + t265) * MDP(19) + ((-t173 + (-qJD(4) + t250) * t185) * t260 + t266) * MDP(20) + t271;];
tau = t1;
