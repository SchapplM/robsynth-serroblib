% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:31
% EndTime: 2020-01-03 11:36:36
% DurationCPUTime: 1.44s
% Computational Cost: add. (1208->212), mult. (2001->295), div. (0->0), fcn. (1186->14), ass. (0->114)
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t233 = cos(pkin(8));
t217 = pkin(1) * t233 + pkin(2);
t202 = t217 * qJD(1);
t231 = sin(pkin(8));
t310 = pkin(1) * t231;
t313 = qJD(3) * t202 + qJDD(1) * t310;
t276 = qJD(3) * t310;
t314 = -qJD(1) * t276 + t217 * qJDD(1);
t263 = -t313 * t235 + t314 * t238;
t253 = qJDD(4) - t263;
t223 = qJDD(1) + qJDD(3);
t309 = pkin(3) * t223;
t171 = t253 - t309;
t227 = qJ(1) + pkin(8);
t221 = qJ(3) + t227;
t215 = sin(t221);
t216 = cos(t221);
t262 = g(2) * t216 + g(3) * t215;
t250 = -t171 - t262;
t226 = qJD(1) + qJD(3);
t305 = qJ(4) * t223;
t311 = -t314 * t235 - t313 * t238;
t169 = qJD(4) * t226 + t305 - t311;
t230 = sin(pkin(9));
t232 = cos(pkin(9));
t166 = qJDD(2) * t230 + t169 * t232;
t165 = -t232 * qJDD(2) + t169 * t230;
t304 = t165 * t230;
t312 = t166 * t232 + t304;
t287 = t235 * t217 + t238 * t310;
t234 = sin(qJ(5));
t237 = cos(qJ(5));
t252 = MDP(17) * t234 + MDP(18) * t237;
t277 = qJD(1) * t310;
t185 = t238 * t202 - t235 * t277;
t260 = qJD(4) - t185;
t306 = MDP(8) * t232;
t186 = t202 * t235 + t238 * t277;
t178 = qJ(4) * t226 + t186;
t175 = -t232 * qJD(2) + t178 * t230;
t303 = t175 * t230;
t195 = -pkin(4) * t232 - pkin(7) * t230 - pkin(3);
t297 = t217 * t238;
t266 = -t235 * t310 + t297;
t179 = t195 - t266;
t302 = t179 * t237;
t301 = t186 * t226;
t251 = qJD(3) * t297 - t235 * t276;
t187 = qJD(4) + t251;
t300 = t187 * t226;
t299 = t187 * t234;
t188 = t287 * qJD(3);
t298 = t188 * t226;
t296 = t223 * t232;
t295 = t223 * t234;
t294 = t226 * t232;
t293 = t230 * t234;
t292 = t232 * t234;
t291 = t232 * t237;
t290 = t234 * t237;
t289 = t216 * pkin(3) + t215 * qJ(4);
t288 = -g(2) * t215 + g(3) * t216;
t224 = t230 ^ 2;
t225 = t232 ^ 2;
t286 = t224 + t225;
t229 = t237 ^ 2;
t285 = t234 ^ 2 - t229;
t281 = qJD(5) * t234;
t280 = qJD(5) * t237;
t201 = -qJDD(5) + t296;
t279 = t201 * MDP(16);
t205 = -qJD(5) + t294;
t278 = -qJD(5) - t205;
t274 = t250 * t230;
t273 = t226 * t280;
t272 = t286 * t223;
t271 = t215 * pkin(3) - t216 * qJ(4);
t164 = t195 * t223 + t253;
t270 = t237 * t164 - t234 * t166;
t268 = t201 - t296;
t267 = t201 + t296;
t265 = t260 * t237;
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t261 = -g(2) * t239 - g(3) * t236;
t259 = -t301 - t309;
t258 = t234 * t164 + t237 * t166;
t172 = t195 * t226 + t260;
t176 = qJD(2) * t230 + t178 * t232;
t257 = t172 * t237 - t176 * t234;
t256 = -t172 * t234 - t176 * t237;
t255 = t176 * t232 + t303;
t190 = -pkin(3) - t266;
t254 = t190 * t223 + t298;
t249 = t288 + t312;
t248 = -qJ(4) * t292 + t195 * t237;
t181 = -t215 * t292 - t216 * t237;
t183 = -t215 * t237 + t216 * t292;
t246 = g(2) * t183 - g(3) * t181 + (t257 * qJD(5) + t258) * t232 + t237 * t304;
t182 = t215 * t291 - t216 * t234;
t184 = t215 * t234 + t216 * t291;
t245 = -g(2) * t184 - g(3) * t182 + t165 * t293 + t280 * t303;
t244 = qJ(4) * t280 + t260 * t234;
t243 = t262 - t263;
t192 = t230 * t281 * t294;
t242 = t232 * t279 + (t192 + (t205 * t281 - t267 * t237) * t230) * MDP(14) + (t267 * t234 + (t205 + t294) * t280) * t230 * MDP(15) + t223 * MDP(5) + (0.2e1 * (t285 * t226 * qJD(5) - t223 * t290) * MDP(13) + (t223 * t229 - 0.2e1 * t234 * t273) * MDP(12)) * t224;
t241 = -t288 + t311;
t222 = t226 ^ 2;
t189 = qJ(4) + t287;
t177 = -pkin(3) * t226 + t260;
t159 = t256 * qJD(5) + t270;
t1 = [(t250 - t254) * t306 + qJDD(1) * MDP(1) + (t189 * t272 + t286 * t300 + t249) * MDP(10) + t242 + (t261 + (t231 ^ 2 + t233 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t254 * t230 - t274) * MDP(9) + (t171 * t190 + t177 * t188 - g(2) * (pkin(2) * cos(t227) + t239 * pkin(1) + t289) - g(3) * (pkin(2) * sin(t227) + t236 * pkin(1) + t271) + t312 * t189 + t255 * t187) * MDP(11) + (-t287 * t223 - t251 * t226 + t241) * MDP(7) + (-(-t179 * t281 + t188 * t237) * t205 - t201 * t302 + (-(-t189 * t280 - t299) * t205 + t189 * t234 * t201 - t159) * t232 + (t226 * t299 + (t273 + t295) * t189) * t224 + t245) * MDP(17) + (t266 * t223 - t243 - t298) * MDP(6) + t261 * MDP(2) + (g(2) * t236 - g(3) * t239) * MDP(3) + ((t187 * t291 + t188 * t234) * t205 + (t179 * t234 + t189 * t291) * t201 + (t189 * t223 + t300) * t237 * t224 + (t205 * t302 + (-t303 + (-t205 * t232 - t224 * t226) * t189) * t234) * qJD(5) + t246) * MDP(18); (qJDD(2) - g(1)) * MDP(4) + (-t165 * t232 - g(1)) * MDP(11) + t192 * MDP(18) + (t166 * MDP(11) + (-qJD(5) * t205 * MDP(18) + t268 * MDP(17)) * t234 + (t268 * MDP(18) + (t205 - t294) * MDP(17) * qJD(5)) * t237) * t230; (-t243 + t301) * MDP(6) + (t185 * t226 + t241) * MDP(7) + (t250 - t259) * t306 + (t259 * t230 - t274) * MDP(9) + (t260 * t226 * t286 + qJ(4) * t272 + t249) * MDP(10) + (-t171 * pkin(3) - t177 * t186 - g(2) * t289 - g(3) * t271 + (t166 * qJ(4) + t260 * t176) * t232 + (t165 * qJ(4) + t260 * t175) * t230) * MDP(11) + (-t248 * t201 - t159 * t232 + (t186 * t237 + t195 * t281 + t244 * t232) * t205 + (qJ(4) * t295 + t244 * t226) * t224 + t245) * MDP(17) + ((qJ(4) * t291 + t195 * t234) * t201 + (-t186 * t234 + t232 * t265) * t205 + (-t175 * t293 + t248 * t205) * qJD(5) + (t237 * t305 + (-qJ(4) * t281 + t265) * t226) * t224 + t246) * MDP(18) + t242; (-t255 * t226 + qJDD(4) + t243) * MDP(11) + (-pkin(3) * MDP(11) + MDP(9) * t230 - t306) * t223 + (-MDP(17) * t237 + MDP(18) * t234) * t201 + (-t225 * MDP(10) + (-MDP(10) - t252) * t224) * t222 - t252 * t205 ^ 2; -t279 + (-g(2) * t181 - g(3) * t183 + t256 * t205 + t270) * MDP(17) + (g(2) * t182 - g(3) * t184 - t257 * t205 - t258) * MDP(18) + (t256 * MDP(17) - t257 * MDP(18)) * qJD(5) + (MDP(12) * t290 - t285 * MDP(13)) * t224 * t222 + ((MDP(14) * t237 - MDP(15) * t234) * t223 + t252 * g(1) + ((t278 * MDP(15) - t175 * MDP(17)) * t237 + (t278 * MDP(14) + t175 * MDP(18)) * t234) * t226) * t230;];
tau = t1;
