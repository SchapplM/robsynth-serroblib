% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:55
% EndTime: 2022-01-20 10:05:58
% DurationCPUTime: 1.51s
% Computational Cost: add. (1209->231), mult. (1931->322), div. (0->0), fcn. (1162->14), ass. (0->130)
t255 = cos(qJ(2));
t315 = pkin(1) * qJD(2);
t287 = qJD(1) * t315;
t252 = sin(qJ(2));
t290 = qJDD(1) * t252;
t326 = pkin(1) * t290 + t255 * t287;
t321 = pkin(1) * t255;
t231 = qJDD(1) * t321;
t240 = qJDD(1) + qJDD(2);
t196 = pkin(2) * t240 - t252 * t287 + t231;
t248 = sin(pkin(8));
t250 = cos(pkin(8));
t180 = t248 * t196 + t250 * t326;
t243 = qJD(1) + qJD(2);
t175 = qJ(4) * t240 + qJD(4) * t243 + t180;
t247 = sin(pkin(9));
t249 = cos(pkin(9));
t174 = qJDD(3) * t247 + t175 * t249;
t173 = -qJDD(3) * t249 + t175 * t247;
t314 = t173 * t247;
t325 = t174 * t249 + t314;
t246 = qJ(1) + qJ(2);
t236 = sin(t246);
t237 = cos(t246);
t324 = g(1) * t236 - g(2) * t237;
t251 = sin(qJ(5));
t254 = cos(qJ(5));
t266 = MDP(16) * t251 + MDP(17) * t254;
t316 = pkin(1) * qJD(1);
t289 = t252 * t316;
t214 = t248 * t289;
t288 = t255 * t316;
t201 = t250 * t288 - t214;
t292 = qJD(4) - t201;
t323 = pkin(1) * t252;
t253 = sin(qJ(1));
t322 = pkin(1) * t253;
t320 = pkin(2) * t236;
t319 = pkin(2) * t250;
t235 = pkin(8) + t246;
t223 = sin(t235);
t318 = g(1) * t223;
t205 = pkin(2) * t243 + t288;
t215 = t250 * t289;
t189 = t205 * t248 + t215;
t185 = qJ(4) * t243 + t189;
t182 = -qJD(3) * t249 + t185 * t247;
t313 = t182 * t247;
t267 = -pkin(4) * t249 - pkin(7) * t247 - pkin(3);
t220 = t248 * t323;
t230 = pkin(2) + t321;
t278 = t230 * t250 - t220;
t186 = t267 - t278;
t312 = t186 * t254;
t265 = t250 * t255 * t315 - qJD(2) * t220;
t195 = qJD(4) + t265;
t311 = t195 * t243;
t310 = t195 * t251;
t309 = t240 * t249;
t308 = t240 * t251;
t307 = t243 * t249;
t306 = t247 * t251;
t305 = t249 * t251;
t304 = t249 * t254;
t303 = t250 * t252;
t302 = t251 * t254;
t301 = pkin(1) * t303 + t230 * t248;
t300 = g(1) * t237 + g(2) * t236;
t241 = t247 ^ 2;
t242 = t249 ^ 2;
t299 = t241 + t242;
t245 = t254 ^ 2;
t298 = t251 ^ 2 - t245;
t295 = qJD(5) * t251;
t294 = qJD(5) * t254;
t209 = -qJDD(5) + t309;
t293 = t209 * MDP(15);
t211 = -qJD(5) + t307;
t291 = -qJD(5) - t211;
t224 = cos(t235);
t229 = pkin(2) * t237;
t285 = t224 * pkin(3) + t223 * qJ(4) + t229;
t284 = t243 * t294;
t179 = t196 * t250 - t248 * t326;
t263 = qJDD(4) - t179;
t176 = -pkin(3) * t240 + t263;
t283 = -g(2) * t224 - t176;
t282 = t299 * t240;
t170 = t240 * t267 + t263;
t281 = t170 * t254 - t251 * t174;
t188 = t205 * t250 - t214;
t280 = t209 - t309;
t279 = t209 + t309;
t277 = t292 * t254;
t276 = qJD(1) * (-qJD(2) + t243);
t275 = qJD(2) * (-qJD(1) - t243);
t273 = t231 + t324;
t272 = qJD(4) - t188;
t271 = t251 * t170 + t254 * t174;
t177 = t243 * t267 + t272;
t183 = qJD(3) * t247 + t185 * t249;
t270 = t177 * t254 - t183 * t251;
t269 = -t177 * t251 - t183 * t254;
t268 = t183 * t249 + t313;
t264 = -pkin(3) * t223 + t224 * qJ(4) - t320;
t202 = t267 - t319;
t222 = pkin(2) * t248 + qJ(4);
t262 = t202 * t254 - t222 * t305;
t190 = t223 * t305 + t224 * t254;
t192 = t223 * t254 - t224 * t305;
t261 = -g(1) * t190 - g(2) * t192 + (qJD(5) * t270 + t271) * t249 + t254 * t314;
t191 = -t223 * t304 + t224 * t251;
t193 = t223 * t251 + t224 * t304;
t260 = -g(1) * t191 - g(2) * t193 + t173 * t306 + t294 * t313;
t259 = t222 * t294 + t251 * t292;
t258 = -g(1) * t224 - g(2) * t223 + t325;
t204 = t247 * t295 * t307;
t257 = t249 * t293 + (t204 + (t211 * t295 - t254 * t279) * t247) * MDP(13) + (t279 * t251 + (t211 + t307) * t294) * t247 * MDP(14) + t240 * MDP(4) + (0.2e1 * (qJD(5) * t243 * t298 - t240 * t302) * MDP(12) + (t240 * t245 - 0.2e1 * t251 * t284) * MDP(11)) * t241;
t256 = cos(qJ(1));
t239 = t243 ^ 2;
t238 = t256 * pkin(1);
t225 = -pkin(3) - t319;
t206 = t249 * t318;
t200 = (t248 * t255 + t303) * t315;
t199 = t248 * t288 + t215;
t198 = -pkin(3) - t278;
t197 = qJ(4) + t301;
t184 = -pkin(3) * t243 + t272;
t165 = qJD(5) * t269 + t281;
t1 = [qJDD(1) * MDP(1) + (g(1) * t253 - g(2) * t256) * MDP(2) + (g(1) * t256 + g(2) * t253) * MDP(3) + ((t240 * t255 + t252 * t275) * pkin(1) + t273) * MDP(5) + (((-qJDD(1) - t240) * t252 + t255 * t275) * pkin(1) + t300) * MDP(6) + (t180 * t301 + t189 * t265 + t179 * t278 - t188 * t200 - g(1) * (-t320 - t322) - g(2) * (t229 + t238)) * MDP(7) + (t206 + (-t198 * t240 - t200 * t243 + t283) * t249) * MDP(8) + (t197 * t282 + t299 * t311 + t258) * MDP(9) + (t176 * t198 + t184 * t200 - g(1) * (t264 - t322) - g(2) * (t238 + t285) + t325 * t197 + t268 * t195) * MDP(10) + (-(-t186 * t295 + t200 * t254) * t211 - t209 * t312 + (-(-t197 * t294 - t310) * t211 + t197 * t251 * t209 - t165) * t249 + (t243 * t310 + (t284 + t308) * t197) * t241 + t260) * MDP(16) + ((t195 * t304 + t200 * t251) * t211 + (t186 * t251 + t197 * t304) * t209 + (t197 * t240 + t311) * t254 * t241 + (t211 * t312 + (-t313 + (-t211 * t249 - t241 * t243) * t197) * t251) * qJD(5) + t261) * MDP(17) + t257; (t276 * t323 + t273) * MDP(5) + ((t255 * t276 - t290) * pkin(1) + t300) * MDP(6) + (t188 * t199 - t189 * t201 + (t179 * t250 + t180 * t248 + t324) * pkin(2)) * MDP(7) + (t206 + (t199 * t243 - t225 * t240 + t283) * t249) * MDP(8) + (t243 * t292 * t299 + t222 * t282 + t258) * MDP(9) + (t176 * t225 - t184 * t199 - g(1) * t264 - g(2) * t285 + (t174 * t222 + t183 * t292) * t249 + (t173 * t222 + t182 * t292) * t247) * MDP(10) + (-t262 * t209 - t165 * t249 + (t254 * t199 + t202 * t295 + t249 * t259) * t211 + (t222 * t308 + t243 * t259) * t241 + t260) * MDP(16) + ((t202 * t251 + t222 * t304) * t209 + (-t199 * t251 + t249 * t277) * t211 + (-t182 * t306 + t211 * t262) * qJD(5) + (t222 * t254 * t240 + (-t222 * t295 + t277) * t243) * t241 + t261) * MDP(17) + t257; (qJDD(3) - g(3)) * MDP(7) + (-t173 * t249 - g(3)) * MDP(10) + t204 * MDP(17) + (t174 * MDP(10) + (-qJD(5) * t211 * MDP(17) + MDP(16) * t280) * t251 + (t280 * MDP(17) + (t211 - t307) * MDP(16) * qJD(5)) * t254) * t247; -MDP(8) * t309 + (-t243 * t268 - t283 - t318) * MDP(10) + (-MDP(16) * t254 + MDP(17) * t251) * t209 + (-t242 * MDP(9) + (-MDP(9) - t266) * t241) * t239 - t266 * t211 ^ 2; -t293 + (-g(1) * t192 + g(2) * t190 + t211 * t269 + t281) * MDP(16) + (g(1) * t193 - g(2) * t191 - t211 * t270 - t271) * MDP(17) + (MDP(16) * t269 - MDP(17) * t270) * qJD(5) + (MDP(11) * t302 - MDP(12) * t298) * t241 * t239 + ((MDP(13) * t254 - MDP(14) * t251) * t240 + t266 * g(3) + ((MDP(14) * t291 - MDP(16) * t182) * t254 + (MDP(13) * t291 + MDP(17) * t182) * t251) * t243) * t247;];
tau = t1;
