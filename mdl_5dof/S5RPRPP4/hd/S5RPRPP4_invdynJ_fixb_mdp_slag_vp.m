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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
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
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:54
% EndTime: 2019-12-31 18:14:57
% DurationCPUTime: 1.80s
% Computational Cost: add. (1360->258), mult. (2458->305), div. (0->0), fcn. (1433->8), ass. (0->113)
t269 = qJD(1) ^ 2;
t264 = sin(qJ(1));
t266 = cos(qJ(1));
t314 = g(1) * t264 - g(2) * t266;
t278 = -qJ(2) * t269 - t314;
t260 = sin(pkin(7));
t261 = cos(pkin(7));
t263 = sin(qJ(3));
t265 = cos(qJ(3));
t284 = t260 * t265 + t261 * t263;
t209 = t284 * qJD(1);
t310 = qJD(1) * t263;
t296 = t260 * t310;
t309 = qJD(1) * t265;
t213 = t261 * t309 - t296;
t220 = pkin(3) * t310 + qJD(1) * qJ(2) + qJD(4);
t186 = pkin(4) * t209 - qJ(5) * t213 + t220;
t253 = qJ(3) + pkin(7);
t244 = sin(t253);
t245 = cos(t253);
t267 = -pkin(1) - pkin(6);
t223 = t267 * qJDD(1) + qJDD(2);
t218 = t265 * t223;
t224 = t267 * qJD(1) + qJD(2);
t304 = qJD(1) * qJD(3);
t295 = t263 * t304;
t301 = qJDD(1) * t265;
t303 = qJD(1) * qJD(4);
t308 = qJD(3) * t263;
t183 = -t265 * t303 - t224 * t308 + qJDD(3) * pkin(3) + t218 + (t295 - t301) * qJ(4);
t290 = -qJ(4) * qJD(1) + t224;
t307 = qJD(3) * t265;
t189 = t290 * t307 + (-qJ(4) * qJDD(1) + t223 - t303) * t263;
t171 = t261 * t183 - t260 * t189;
t293 = -qJDD(5) + t171;
t331 = g(3) * t244 - t186 * t213 - t314 * t245 + t293;
t294 = t265 * t304;
t302 = qJDD(1) * t263;
t330 = t294 + t302;
t329 = MDP(14) + MDP(17);
t255 = qJDD(1) * qJ(2);
t288 = g(1) * t266 + g(2) * t264;
t256 = qJD(1) * qJD(2);
t299 = 0.2e1 * t256;
t327 = 0.2e1 * t255 + t299 - t288;
t316 = qJ(4) - t267;
t291 = t316 * t265;
t204 = -qJD(3) * t291 - qJD(4) * t263;
t275 = -qJD(4) * t265 + t316 * t308;
t187 = t204 * t260 - t261 * t275;
t188 = t261 * t204 + t260 * t275;
t297 = t260 * t301 + t261 * t330;
t194 = t260 * t295 - t297;
t285 = t260 * t302 - t261 * t301;
t195 = t284 * t304 + t285;
t221 = t316 * t263;
t197 = -t221 * t260 + t261 * t291;
t198 = -t261 * t221 - t260 * t291;
t326 = t187 * t213 - t188 * t209 + t198 * t194 - t195 * t197;
t208 = t213 ^ 2;
t324 = g(3) * t263;
t249 = t263 * pkin(3);
t323 = pkin(1) * qJDD(1);
t321 = qJDD(3) * pkin(4);
t205 = t290 * t263;
t318 = t205 * t260;
t201 = t261 * t205;
t317 = qJ(2) + t249;
t172 = t260 * t183 + t261 * t189;
t206 = -qJ(4) * t309 + t265 * t224;
t203 = qJD(3) * pkin(3) + t206;
t185 = t260 * t203 + t201;
t315 = t266 * pkin(1) + t264 * qJ(2);
t259 = t265 ^ 2;
t313 = t263 ^ 2 - t259;
t268 = qJD(3) ^ 2;
t312 = -t268 - t269;
t311 = qJD(1) * t220;
t306 = pkin(3) * t307 + qJD(2);
t191 = t206 * t261 - t318;
t305 = qJD(5) - t191;
t300 = qJDD(3) * t263;
t298 = qJDD(3) * qJ(5) + t172;
t292 = -pkin(1) * t264 + t266 * qJ(2);
t289 = qJDD(2) - t323;
t286 = pkin(4) * t244 - qJ(5) * t245;
t184 = t203 * t261 - t318;
t283 = pkin(3) * t330 + qJDD(4) + t255 + t256;
t262 = -qJ(4) - pkin(6);
t282 = t266 * t249 + t264 * t262 + t292;
t281 = t264 * t249 - t262 * t266 + t315;
t280 = 0.2e1 * qJ(2) * t304 + qJDD(3) * t267;
t276 = -t260 * t309 - t261 * t310;
t273 = -pkin(4) * t194 + qJ(5) * t195 + t283;
t169 = qJD(3) * qJD(5) + t298;
t170 = -t293 - t321;
t178 = -qJD(3) * pkin(4) + qJD(5) - t184;
t179 = qJD(3) * qJ(5) + t185;
t211 = t260 * t308 - t261 * t307;
t212 = -t260 * t307 - t261 * t308;
t216 = -t260 * t263 + t261 * t265;
t272 = t169 * t284 - t170 * t216 - t178 * t212 - t179 * t211 - t314;
t271 = t171 * t216 + t172 * t284 + t184 * t212 - t185 * t211 - t314;
t270 = -t267 * t268 + t327;
t246 = qJDD(3) * t265;
t239 = -pkin(3) * t261 - pkin(4);
t235 = pkin(3) * t260 + qJ(5);
t193 = pkin(4) * t284 - qJ(5) * t216 + t317;
t192 = pkin(3) * t309 + pkin(4) * t213 + qJ(5) * t209;
t190 = t206 * t260 + t201;
t175 = -pkin(4) * t211 - qJ(5) * t212 - qJD(5) * t216 + t306;
t168 = -qJD(5) * t213 + t273;
t1 = [qJDD(1) * MDP(1) + t314 * MDP(2) + t288 * MDP(3) + (qJDD(2) - t314 - 0.2e1 * t323) * MDP(4) + t327 * MDP(5) + (-t289 * pkin(1) - g(1) * t292 - g(2) * t315 + (t299 + t255) * qJ(2)) * MDP(6) + (qJDD(1) * t259 - 0.2e1 * t263 * t294) * MDP(7) + 0.2e1 * (-t263 * t301 + t313 * t304) * MDP(8) + (-t263 * t268 + t246) * MDP(9) + (-t265 * t268 - t300) * MDP(10) + (t270 * t263 + t280 * t265) * MDP(12) + (-t280 * t263 + t270 * t265) * MDP(13) + (-t271 + t326) * MDP(14) + (-g(1) * t282 - g(2) * t281 - t171 * t197 + t172 * t198 - t184 * t187 + t185 * t188 + t220 * t306 + t283 * t317) * MDP(15) + (-qJD(3) * t187 - qJDD(3) * t197 + t168 * t284 + t175 * t209 - t186 * t211 - t193 * t194 - t288 * t244) * MDP(16) + (-t272 + t326) * MDP(17) + (qJD(3) * t188 + qJDD(3) * t198 - t168 * t216 - t175 * t213 - t186 * t212 + t193 * t195 + t288 * t245) * MDP(18) + (t169 * t198 + t179 * t188 + t168 * t193 + t186 * t175 + t170 * t197 + t178 * t187 - g(1) * (t286 * t266 + t282) - g(2) * (t286 * t264 + t281)) * MDP(19); qJDD(1) * MDP(4) - t269 * MDP(5) + (t289 + t278) * MDP(6) + (t312 * t263 + t246) * MDP(12) + (t312 * t265 - t300) * MDP(13) + (t271 - t311) * MDP(15) + (-qJD(1) * t209 + qJD(3) * t212 + qJDD(3) * t216) * MDP(16) + (qJD(1) * t213 - qJD(3) * t211 + qJDD(3) * t284) * MDP(18) + (-qJD(1) * t186 + t272) * MDP(19) + t329 * (t194 * t284 + t195 * t216 + t211 * t209 - t212 * t213); MDP(9) * t301 - MDP(10) * t302 + qJDD(3) * MDP(11) + (t278 * t265 + t218 + t324) * MDP(12) + (g(3) * t265 + (-t223 - t278) * t263) * MDP(13) + ((t185 - t190) * t213 + (-t184 + t191) * t209 + (t194 * t260 + t195 * t261) * pkin(3)) * MDP(14) + (t184 * t190 - t185 * t191 + (t324 + t171 * t261 + t172 * t260 + (-t314 - t311) * t265) * pkin(3)) * MDP(15) + (qJD(3) * t190 - t192 * t209 + (pkin(4) - t239) * qJDD(3) + t331) * MDP(16) + (t194 * t235 - t195 * t239 + (t179 - t190) * t213 + (t178 - t305) * t209) * MDP(17) + (-g(3) * t245 + qJDD(3) * t235 - t186 * t209 + t192 * t213 - t314 * t244 + (0.2e1 * qJD(5) - t191) * qJD(3) + t298) * MDP(18) + (t169 * t235 + t170 * t239 - t186 * t192 - t178 * t190 - g(3) * (-t249 - t286) + t305 * t179 - t314 * (pkin(3) * t265 + pkin(4) * t245 + qJ(5) * t244)) * MDP(19) + (t265 * t263 * MDP(7) - t313 * MDP(8)) * t269; (t184 * t213 + t185 * t209 + t283 - t288) * MDP(15) + t297 * MDP(16) + t285 * MDP(18) + (t179 * t209 + (-qJD(5) - t178) * t213 + t273 - t288) * MDP(19) + ((t213 - t296) * MDP(16) + (t209 - t276) * MDP(18)) * qJD(3) + t329 * (-t209 ^ 2 - t208); (t209 * t213 - qJDD(3)) * MDP(16) - t285 * MDP(17) + (-t208 - t268) * MDP(18) + (-t321 - t331) * MDP(19) + ((t209 + t276) * MDP(17) - t179 * MDP(19)) * qJD(3);];
tau = t1;
