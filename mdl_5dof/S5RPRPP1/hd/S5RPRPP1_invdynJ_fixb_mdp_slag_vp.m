% Calculate vector of inverse dynamics joint torques for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:08
% EndTime: 2021-01-15 11:14:13
% DurationCPUTime: 2.17s
% Computational Cost: add. (1588->274), mult. (3205->333), div. (0->0), fcn. (2020->12), ass. (0->124)
t331 = 2 * qJD(3);
t255 = sin(pkin(7));
t237 = pkin(1) * t255 + pkin(6);
t310 = qJ(4) + t237;
t257 = cos(pkin(7));
t239 = -pkin(1) * t257 - pkin(2);
t261 = cos(qJ(3));
t247 = t261 * pkin(3);
t327 = t239 - t247;
t330 = qJDD(1) * t327;
t329 = MDP(14) + MDP(17);
t250 = qJ(1) + pkin(7);
t244 = cos(t250);
t242 = sin(t250);
t321 = g(1) * t242;
t285 = -g(2) * t244 + t321;
t319 = g(2) * t242;
t286 = g(1) * t244 + t319;
t254 = sin(pkin(8));
t256 = cos(pkin(8));
t259 = sin(qJ(3));
t217 = t254 * t261 + t256 * t259;
t212 = t217 * qJD(1);
t297 = MDP(12) + MDP(16);
t211 = t217 * qJD(3);
t299 = qJDD(1) * t261;
t230 = t256 * t299;
t300 = qJDD(1) * t259;
t192 = qJD(1) * t211 + t254 * t300 - t230;
t301 = qJD(1) * qJD(3);
t292 = t259 * t301;
t271 = t217 * qJDD(1) - t254 * t292;
t291 = t261 * t301;
t193 = t256 * t291 + t271;
t326 = pkin(4) * t192 - qJ(5) * t193 - qJD(5) * t212;
t215 = t310 * t261;
t290 = t310 * t259;
t191 = t256 * t215 - t254 * t290;
t249 = qJ(3) + pkin(8);
t241 = sin(t249);
t325 = -qJDD(3) * t191 - t285 * t241;
t245 = t261 * qJDD(2);
t225 = t237 * qJDD(1);
t269 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + (qJD(2) * qJD(3)) + t225;
t288 = t310 * qJD(1);
t280 = t288 * qJD(3);
t176 = qJDD(3) * pkin(3) - t269 * t259 - t261 * t280 + t245;
t179 = (qJDD(2) - t280) * t259 + t269 * t261;
t166 = t254 * t176 + t256 * t179;
t243 = cos(t249);
t324 = g(3) * t241 + t286 * t243 - t166;
t207 = t212 ^ 2;
t323 = pkin(3) * t259;
t318 = g(3) * t261;
t316 = qJDD(3) * pkin(4);
t203 = qJD(2) * t259 + t288 * t261;
t315 = t203 * t254;
t314 = t241 * t244;
t313 = t243 * t244;
t312 = t254 * t259;
t197 = t256 * t203;
t311 = t256 * t261;
t309 = qJDD(2) - g(3);
t165 = t256 * t176 - t254 * t179;
t202 = t261 * qJD(2) - t288 * t259;
t199 = qJD(3) * pkin(3) + t202;
t178 = t254 * t199 + t197;
t240 = t247 + pkin(2);
t262 = cos(qJ(1));
t308 = t262 * pkin(1) + t244 * t240;
t252 = t259 ^ 2;
t307 = -t261 ^ 2 + t252;
t228 = qJD(1) * t239;
t306 = qJD(1) * t259;
t181 = t202 * t254 + t197;
t305 = qJD(3) * t181;
t304 = qJD(3) * t259;
t182 = t202 * t256 - t315;
t302 = qJD(5) - t182;
t296 = MDP(13) - MDP(18);
t295 = pkin(3) * t292 + qJDD(4);
t294 = pkin(3) * t304;
t293 = qJD(1) * t311;
t289 = qJD(3) * t310;
t260 = sin(qJ(1));
t284 = g(1) * t260 - g(2) * t262;
t258 = -qJ(4) - pkin(6);
t283 = -pkin(1) * t260 - t244 * t258;
t282 = -pkin(4) * t243 - qJ(5) * t241;
t177 = t199 * t256 - t315;
t278 = t295 - t285;
t277 = g(1) * t314 - g(3) * t243 + t241 * t319 + t165;
t190 = t215 * t254 + t256 * t290;
t275 = -g(2) * t313 - qJDD(3) * t190 + t243 * t321;
t273 = -qJD(4) * t259 - t261 * t289;
t272 = -qJD(1) * t228 - t225 + t286;
t208 = qJD(1) * t327 + qJD(4);
t270 = -qJDD(3) * t237 + t228 * t331;
t201 = t295 + t330;
t209 = t254 * t306 - t293;
t186 = pkin(4) * t209 - qJ(5) * t212 + t208;
t268 = -t186 * t212 - qJDD(5) + t277;
t263 = qJD(3) ^ 2;
t267 = -0.2e1 * qJDD(1) * t239 - t237 * t263 + t285;
t204 = qJD(4) * t261 - t259 * t289;
t184 = t204 * t254 - t256 * t273;
t185 = t256 * t204 + t254 * t273;
t266 = t184 * t212 - t185 * t209 + t190 * t193 - t191 * t192 - t286;
t251 = qJDD(3) * qJ(5);
t238 = -pkin(3) * t256 - pkin(4);
t234 = pkin(3) * t254 + qJ(5);
t224 = qJDD(3) * t261 - t259 * t263;
t223 = qJDD(3) * t259 + t261 * t263;
t216 = -t311 + t312;
t214 = qJD(3) * t311 - t254 * t304;
t189 = pkin(4) * t216 - qJ(5) * t217 + t327;
t188 = pkin(3) * t306 + pkin(4) * t212 + qJ(5) * t209;
t180 = pkin(4) * t211 - qJ(5) * t214 - qJD(5) * t217 + t294;
t173 = qJD(3) * qJ(5) + t178;
t172 = -qJD(3) * pkin(4) + qJD(5) - t177;
t167 = t201 + t326;
t164 = qJDD(5) - t165 - t316;
t163 = qJD(3) * qJD(5) + t166 + t251;
t1 = [qJDD(1) * MDP(1) + t284 * MDP(2) + (g(1) * t262 + g(2) * t260) * MDP(3) + (t284 + (t255 ^ 2 + t257 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t252 + 0.2e1 * t259 * t291) * MDP(5) + 0.2e1 * (t259 * t299 - t307 * t301) * MDP(6) + t223 * MDP(7) + t224 * MDP(8) + (t270 * t259 + t267 * t261) * MDP(10) + (-t267 * t259 + t270 * t261) * MDP(11) + (t192 * t327 + t201 * t216 + t208 * t211 + (t209 * t323 - t184) * qJD(3) + t275) * MDP(12) + (t193 * t327 + t201 * t217 + t208 * t214 + (t212 * t323 - t185) * qJD(3) + t325) * MDP(13) + (-t165 * t217 - t166 * t216 - t177 * t214 - t178 * t211 + t266) * MDP(14) + (t166 * t191 + t178 * t185 - t165 * t190 - t177 * t184 + t201 * t327 + t208 * t294 - g(1) * (-t240 * t242 + t283) - g(2) * (-t242 * t258 + t308)) * MDP(15) + (-qJD(3) * t184 + t167 * t216 + t180 * t209 + t186 * t211 + t189 * t192 + t275) * MDP(16) + (-t163 * t216 + t164 * t217 + t172 * t214 - t173 * t211 + t266) * MDP(17) + (qJD(3) * t185 - t167 * t217 - t180 * t212 - t186 * t214 - t189 * t193 - t325) * MDP(18) + (t163 * t191 + t173 * t185 + t167 * t189 + t186 * t180 + t164 * t190 + t172 * t184 - g(1) * t283 - g(2) * (pkin(4) * t313 + qJ(5) * t314 + t308) + (-g(1) * (-t240 + t282) + g(2) * t258) * t242) * MDP(19); t309 * MDP(4) + t224 * MDP(10) - t223 * MDP(11) + (-t165 * t216 + t166 * t217 - t177 * t211 + t178 * t214 - g(3)) * MDP(15) + (t163 * t217 + t164 * t216 + t172 * t211 + t173 * t214 - g(3)) * MDP(19) - t296 * (qJD(3) * t214 + qJDD(3) * t217) + t297 * (-qJD(3) * t211 - qJDD(3) * t216) + t329 * (-t217 * t192 + t193 * t216 - t214 * t209 + t211 * t212); MDP(7) * t300 + MDP(8) * t299 + qJDD(3) * MDP(9) + (t272 * t259 + t245 - t318) * MDP(10) + (-t309 * t259 + t272 * t261) * MDP(11) + (t305 - t208 * t212 + (qJDD(3) * t256 - t209 * t306) * pkin(3) + t277) * MDP(12) + (qJD(3) * t182 + t208 * t209 + (-qJDD(3) * t254 - t212 * t306) * pkin(3) + t324) * MDP(13) + ((t178 - t181) * t212 + (-t177 + t182) * t209 + (-t192 * t254 - t193 * t256) * pkin(3)) * MDP(14) + (t177 * t181 - t178 * t182 + (-t318 + t165 * t256 + t166 * t254 + (-qJD(1) * t208 + t286) * t259) * pkin(3)) * MDP(15) + (t305 - t188 * t209 + (pkin(4) - t238) * qJDD(3) + t268) * MDP(16) + (-t192 * t234 + t193 * t238 + (t173 - t181) * t212 + (t172 - t302) * t209) * MDP(17) + (qJDD(3) * t234 - t186 * t209 + t188 * t212 + t251 + (0.2e1 * qJD(5) - t182) * qJD(3) - t324) * MDP(18) + (t163 * t234 + t164 * t238 - t186 * t188 - t172 * t181 - g(3) * (t247 - t282) + t302 * t173 + t286 * (pkin(4) * t241 - qJ(5) * t243 + t323)) * MDP(19) + (-t259 * t261 * MDP(5) + t307 * MDP(6)) * qJD(1) ^ 2; (t177 * t212 + t178 * t209 + t278) * MDP(15) + (-t172 * t212 + t173 * t209 + t278 + t326) * MDP(19) + t296 * ((-t209 + t293) * qJD(3) + t271) + (MDP(15) + MDP(19)) * t330 + t329 * (-t209 ^ 2 - t207) + (t312 * qJDD(1) + t212 * t331 - t230) * t297; (t209 * t212 - qJDD(3)) * MDP(16) + ((t209 + t293) * qJD(3) + t271) * MDP(17) + (-t207 - t263) * MDP(18) + (-qJD(3) * t173 - t268 - t316) * MDP(19);];
tau = t1;
