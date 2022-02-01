% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:08
% EndTime: 2022-01-23 09:19:10
% DurationCPUTime: 1.12s
% Computational Cost: add. (1043->184), mult. (1711->240), div. (0->0), fcn. (1095->16), ass. (0->109)
t255 = cos(pkin(8));
t234 = pkin(1) * t255 + pkin(2);
t253 = sin(pkin(8));
t315 = pkin(1) * t253;
t292 = qJD(3) * t315;
t320 = -qJD(1) * t292 + t234 * qJDD(1);
t219 = t234 * qJD(1);
t319 = qJD(3) * t219 + qJDD(1) * t315;
t254 = cos(pkin(9));
t259 = cos(qJ(5));
t299 = t259 * t254;
t252 = sin(pkin(9));
t256 = sin(qJ(5));
t302 = t252 * t256;
t202 = -t299 + t302;
t203 = t252 * t259 + t254 * t256;
t250 = qJD(1) + qJD(3);
t195 = t203 * t250;
t295 = t252 ^ 2 + t254 ^ 2;
t318 = t295 * t250;
t246 = qJDD(1) + qJDD(3);
t257 = sin(qJ(3));
t260 = cos(qJ(3));
t316 = -t320 * t257 - t319 * t260;
t169 = qJ(4) * t246 + qJD(4) * t250 - t316;
t238 = t254 * qJDD(2);
t165 = -t169 * t252 + t238;
t166 = t252 * qJDD(2) + t254 * t169;
t317 = -t165 * t252 + t166 * t254;
t296 = t257 * t234 + t260 * t315;
t293 = qJD(1) * t315;
t190 = t260 * t219 - t257 * t293;
t278 = qJD(4) - t190;
t314 = pkin(3) * t246;
t313 = pkin(4) * t254;
t251 = qJ(1) + pkin(8);
t243 = qJ(3) + t251;
t232 = sin(t243);
t227 = g(1) * t232;
t233 = cos(t243);
t312 = g(2) * t233;
t191 = t219 * t257 + t260 * t293;
t310 = t191 * t250;
t196 = t296 * qJD(3);
t309 = t196 * t250;
t249 = pkin(9) + qJ(5);
t241 = sin(t249);
t308 = t232 * t241;
t242 = cos(t249);
t307 = t232 * t242;
t306 = t233 * t241;
t305 = t233 * t242;
t304 = t234 * t260;
t303 = t246 * t254;
t298 = t233 * pkin(3) + t232 * qJ(4);
t297 = -g(1) * t233 - g(2) * t232;
t290 = t250 * t302;
t289 = t250 * t299;
t288 = qJD(5) * t289 + t203 * t246;
t235 = -pkin(3) - t313;
t287 = -t232 * pkin(3) + t233 * qJ(4);
t281 = -t319 * t257 + t320 * t260;
t271 = qJDD(4) - t281;
t170 = t271 - t314;
t286 = -t170 - t312;
t285 = t295 * t246;
t283 = -t257 * t315 + t304;
t198 = -pkin(3) - t283;
t258 = sin(qJ(1));
t261 = cos(qJ(1));
t280 = g(1) * t258 - g(2) * t261;
t279 = t202 * t246;
t175 = -qJD(5) * t290 + t288;
t200 = t203 * qJD(5);
t176 = t250 * t200 + t279;
t199 = t202 * qJD(5);
t182 = -qJD(5) * t199 + qJDD(5) * t203;
t183 = -qJD(5) * t200 - qJDD(5) * t202;
t193 = -t289 + t290;
t277 = (-t175 * t202 - t176 * t203 + t193 * t199 - t195 * t200) * MDP(12) + (t175 * t203 - t195 * t199) * MDP(11) + t182 * MDP(13) + t183 * MDP(14) + t246 * MDP(5);
t186 = qJ(4) * t250 + t191;
t177 = t254 * qJD(2) - t186 * t252;
t178 = t252 * qJD(2) + t254 * t186;
t276 = t177 * t252 - t178 * t254;
t197 = qJ(4) + t296;
t187 = (-pkin(7) - t197) * t252;
t244 = t254 * pkin(7);
t188 = t197 * t254 + t244;
t275 = t187 * t259 - t188 * t256;
t274 = t187 * t256 + t188 * t259;
t210 = (-pkin(7) - qJ(4)) * t252;
t211 = qJ(4) * t254 + t244;
t273 = t210 * t259 - t211 * t256;
t272 = t210 * t256 + t211 * t259;
t270 = qJD(3) * t304 - t257 * t292;
t269 = t297 + t317;
t167 = t235 * t246 + t271;
t179 = t235 * t250 + t278;
t267 = -g(1) * t308 + g(2) * t306 + t167 * t203 - t179 * t199;
t266 = g(1) * t307 - g(2) * t305 + t167 * t202 + t179 * t200;
t264 = -t227 - t281 + t312;
t263 = -t297 + t316;
t209 = t254 * t227;
t192 = qJD(4) + t270;
t189 = t198 - t313;
t185 = -pkin(3) * t250 + t278;
t161 = pkin(7) * t303 + t166;
t160 = t238 + (-pkin(7) * t246 - t169) * t252;
t1 = [qJDD(1) * MDP(1) + t280 * MDP(2) + (g(1) * t261 + g(2) * t258) * MDP(3) + (t280 + (t253 ^ 2 + t255 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t283 * t246 - t264 - t309) * MDP(6) + (-t296 * t246 - t270 * t250 + t263) * MDP(7) + (t209 + (-t198 * t246 + t286 - t309) * t254) * MDP(8) + (t192 * t318 + t197 * t285 + t269) * MDP(9) + (t170 * t198 + t185 * t196 - g(1) * (-pkin(2) * sin(t251) - t258 * pkin(1) + t287) - g(2) * (pkin(2) * cos(t251) + t261 * pkin(1) + t298) + t317 * t197 - t276 * t192) * MDP(10) + (t196 * t193 + t189 * t176 + t275 * qJDD(5) + (-t274 * qJD(5) - t203 * t192) * qJD(5) + t266) * MDP(16) + (t196 * t195 + t189 * t175 - t274 * qJDD(5) + (-t275 * qJD(5) + t202 * t192) * qJD(5) + t267) * MDP(17) + t277; (qJDD(2) - g(3)) * MDP(4) + (t165 * t254 + t166 * t252 - g(3)) * MDP(10) + t183 * MDP(16) - t182 * MDP(17); (-t264 + t310) * MDP(6) + (t190 * t250 + t263) * MDP(7) + (t209 + (t286 + t310 + t314) * t254) * MDP(8) + (qJ(4) * t285 + t278 * t318 + t269) * MDP(9) + (-t170 * pkin(3) - t185 * t191 - g(1) * t287 - g(2) * t298 + (t166 * qJ(4) + t278 * t178) * t254 + (-t165 * qJ(4) - t278 * t177) * t252) * MDP(10) + (t235 * t176 + t273 * qJDD(5) - t191 * t193 + (-t272 * qJD(5) - t278 * t203) * qJD(5) + t266) * MDP(16) + (t235 * t175 - t272 * qJDD(5) - t191 * t195 + (-t273 * qJD(5) + t278 * t202) * qJD(5) + t267) * MDP(17) + t277; -MDP(8) * t303 + (t276 * t250 - t227 - t286) * MDP(10) + t279 * MDP(16) + t288 * MDP(17) - t295 * MDP(9) * t250 ^ 2 + (0.2e1 * t195 * MDP(16) + (-t193 - t290) * MDP(17)) * qJD(5); t195 * t193 * MDP(11) + (-t193 ^ 2 + t195 ^ 2) * MDP(12) - t279 * MDP(14) + qJDD(5) * MDP(15) + (g(1) * t306 + g(2) * t308 - g(3) * t242 + t259 * t160 - t256 * t161 - t179 * t195) * MDP(16) + (g(1) * t305 + g(2) * t307 + g(3) * t241 - t256 * t160 - t259 * t161 + t179 * t193) * MDP(17) + (t288 + (t193 - t290) * qJD(5)) * MDP(13);];
tau = t1;
