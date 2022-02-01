% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:55
% EndTime: 2022-01-20 11:30:57
% DurationCPUTime: 1.18s
% Computational Cost: add. (1290->203), mult. (2175->267), div. (0->0), fcn. (1241->16), ass. (0->124)
t242 = qJD(1) + qJD(2);
t254 = cos(qJ(2));
t312 = pkin(1) * qJD(1);
t208 = pkin(2) * t242 + t254 * t312;
t253 = cos(qJ(3));
t295 = qJD(2) * t254;
t284 = qJD(1) * t295;
t250 = sin(qJ(2));
t288 = qJDD(1) * t250;
t266 = t284 + t288;
t263 = t266 * pkin(1);
t323 = (qJD(3) * t208 + t263) * t253;
t241 = qJDD(1) + qJDD(2);
t236 = qJDD(3) + t241;
t249 = sin(qJ(3));
t317 = pkin(1) * t254;
t233 = qJDD(1) * t317;
t287 = t250 * t312;
t194 = pkin(2) * t241 - qJD(2) * t287 + t233;
t191 = t253 * t194;
t293 = qJD(3) * t249;
t275 = -t208 * t293 + t191;
t291 = qJD(3) * t253;
t285 = t250 * t291;
t167 = pkin(3) * t236 + (-t249 * t288 + (-t249 * t295 - t285) * qJD(1)) * pkin(1) + t275;
t292 = qJD(3) * t250;
t283 = qJD(1) * t292;
t278 = pkin(1) * t283;
t211 = t249 * t278;
t281 = t194 * t249 - t211;
t173 = t281 + t323;
t246 = sin(pkin(9));
t247 = cos(pkin(9));
t162 = t167 * t247 - t173 * t246;
t245 = qJ(1) + qJ(2);
t240 = qJ(3) + t245;
t228 = pkin(9) + t240;
t215 = cos(t228);
t322 = -pkin(4) * t236 + g(2) * t215 - t162;
t229 = sin(t240);
t230 = cos(t240);
t320 = g(1) * t230 + g(2) * t229;
t319 = g(1) * t229 - g(2) * t230;
t231 = pkin(2) * t253 + pkin(3);
t309 = t246 * t249;
t268 = -pkin(2) * t309 + t231 * t247;
t195 = -pkin(4) - t268;
t307 = t247 * t249;
t298 = pkin(2) * t307 + t246 * t231;
t196 = pkin(8) + t298;
t237 = qJD(3) + t242;
t256 = qJD(5) ^ 2;
t305 = t250 * t253;
t271 = -t249 * t254 - t305;
t200 = t271 * t312;
t306 = t249 * t250;
t270 = t253 * t254 - t306;
t201 = t270 * t312;
t311 = pkin(2) * qJD(3);
t301 = -t200 * t247 + t201 * t246 - (t246 * t253 + t307) * t311;
t318 = -t195 * t236 - t196 * t256 + t301 * t237;
t316 = pkin(2) * t236;
t214 = sin(t228);
t315 = g(1) * t214;
t190 = t208 * t249 + t253 * t287;
t310 = t190 * t246;
t308 = t247 * t190;
t304 = qJDD(4) - g(3);
t189 = t253 * t208 - t249 * t287;
t186 = pkin(3) * t237 + t189;
t174 = t186 * t247 - t310;
t171 = -pkin(4) * t237 - t174;
t248 = sin(qJ(5));
t252 = cos(qJ(5));
t303 = t171 * qJD(5) * t248 + t252 * t315;
t163 = t246 * t167 + t247 * t173;
t232 = pkin(2) + t317;
t213 = t253 * t232;
t199 = -pkin(1) * t306 + pkin(3) + t213;
t202 = pkin(1) * t305 + t232 * t249;
t302 = t246 * t199 + t247 * t202;
t300 = -t200 * t246 - t201 * t247 + (t247 * t253 - t309) * t311;
t239 = cos(t245);
t299 = pkin(2) * t239 + pkin(3) * t230;
t238 = sin(t245);
t297 = g(1) * t239 + g(2) * t238;
t243 = t248 ^ 2;
t296 = -t252 ^ 2 + t243;
t290 = qJD(5) * t237;
t289 = qJD(5) * t252;
t286 = t171 * t289 + t322 * t248;
t280 = qJD(1) * (-qJD(2) + t242);
t279 = qJD(2) * (-qJD(1) - t242);
t277 = g(1) * t238 - g(2) * t239 + t233;
t276 = -pkin(2) * t238 - pkin(3) * t229;
t209 = qJDD(5) * t248 + t252 * t256;
t210 = qJDD(5) * t252 - t248 * t256;
t274 = 0.2e1 * (t236 * t248 * t252 - t296 * t290) * MDP(12) + (0.2e1 * t237 * t248 * t289 + t236 * t243) * MDP(11) + t209 * MDP(13) + t210 * MDP(14) + t236 * MDP(7);
t272 = t199 * t247 - t202 * t246;
t269 = -t281 + t320;
t267 = t241 * MDP(4) + t274;
t182 = t232 * t291 + (qJD(2) * t270 - t249 * t292) * pkin(1);
t183 = -t232 * t293 + (qJD(2) * t271 - t285) * pkin(1);
t164 = t182 * t246 - t183 * t247;
t178 = -pkin(4) - t272;
t179 = pkin(8) + t302;
t265 = t164 * t237 + t178 * t236 + t179 * t256;
t176 = t189 * t246 + t308;
t222 = pkin(3) * t246 + pkin(8);
t223 = -pkin(3) * t247 - pkin(4);
t264 = -t176 * t237 + t222 * t256 + t223 * t236;
t262 = t275 + t319;
t261 = -pkin(8) * t236 + g(1) * t215 + g(2) * t214 - t171 * t237 - t163;
t165 = t182 * t247 + t183 * t246;
t260 = -qJDD(5) * t179 + (t178 * t237 - t165) * qJD(5);
t177 = t189 * t247 - t310;
t259 = qJD(5) * t177 - qJDD(5) * t222 + t223 * t290;
t258 = -qJDD(5) * t196 + (t195 * t237 - t300) * qJD(5);
t257 = (-pkin(2) * t237 - t208) * qJD(3) - t263;
t255 = cos(qJ(1));
t251 = sin(qJ(1));
t235 = t237 ^ 2;
t175 = t246 * t186 + t308;
t1 = [qJDD(1) * MDP(1) + (g(1) * t251 - g(2) * t255) * MDP(2) + (g(1) * t255 + g(2) * t251) * MDP(3) + ((t241 * t254 + t250 * t279) * pkin(1) + t277) * MDP(5) + (((-qJDD(1) - t241) * t250 + t254 * t279) * pkin(1) + t297) * MDP(6) + (t183 * t237 + t213 * t236 + (-t253 * t283 + (-t284 + (-qJDD(1) - t236) * t250) * t249) * pkin(1) + t262) * MDP(8) + (-t182 * t237 - t202 * t236 + t269 - t323) * MDP(9) + (t163 * t302 + t175 * t165 + t162 * t272 - t174 * t164 - g(1) * (-pkin(1) * t251 + t276) - g(2) * (pkin(1) * t255 + t299)) * MDP(10) + (t260 * t248 + (-t265 - t322) * t252 + t303) * MDP(16) + (t260 * t252 + (t265 - t315) * t248 + t286) * MDP(17) + t267; (pkin(1) * t250 * t280 + t277) * MDP(5) + ((t254 * t280 - t288) * pkin(1) + t297) * MDP(6) + (-t200 * t237 + t191 + (-t278 + t316) * t253 + t257 * t249 + t319) * MDP(8) + (t201 * t237 + t211 + (-t194 - t316) * t249 + t257 * t253 + t320) * MDP(9) + (-g(1) * t276 - g(2) * t299 + t162 * t268 + t163 * t298 + t301 * t174 + t300 * t175) * MDP(10) + (t258 * t248 + (-t322 + t318) * t252 + t303) * MDP(16) + (t258 * t252 + (-t315 - t318) * t248 + t286) * MDP(17) + t267; (t190 * t237 + t262) * MDP(8) + (t189 * t237 - t208 * t291 + t269) * MDP(9) + t303 * MDP(16) + t286 * MDP(17) + (-t266 * MDP(8) * t249 + (-MDP(8) * t283 - MDP(9) * t266) * t253) * pkin(1) + (t259 * MDP(16) + (t264 - t315) * MDP(17)) * t248 + ((-t264 - t322) * MDP(16) + t259 * MDP(17)) * t252 + t274 + (t174 * t176 - t175 * t177 + (t162 * t247 + t163 * t246 + t319) * pkin(3)) * MDP(10); t304 * MDP(10) + t210 * MDP(16) - t209 * MDP(17); qJDD(5) * MDP(15) + t296 * MDP(12) * t235 + (t236 * MDP(14) + t304 * MDP(16) + t261 * MDP(17)) * t252 + (-t235 * t252 * MDP(11) + t236 * MDP(13) + t261 * MDP(16) - t304 * MDP(17)) * t248;];
tau = t1;
