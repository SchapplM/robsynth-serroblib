% Calculate vector of inverse dynamics joint torques for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:36
% EndTime: 2019-12-05 16:04:43
% DurationCPUTime: 2.25s
% Computational Cost: add. (894->283), mult. (1920->402), div. (0->0), fcn. (1452->10), ass. (0->131)
t224 = sin(qJ(4));
t278 = qJD(2) * qJD(4);
t264 = t224 * t278;
t227 = cos(qJ(4));
t273 = qJDD(2) * t227;
t317 = -t264 + t273;
t229 = -pkin(2) - pkin(7);
t220 = sin(pkin(5));
t225 = sin(qJ(2));
t306 = t220 * t225;
t269 = qJD(2) * t306;
t207 = qJD(1) * t269;
t228 = cos(qJ(2));
t277 = qJDD(1) * t220;
t261 = t228 * t277;
t247 = qJDD(3) + t207 - t261;
t182 = t229 * qJDD(2) + t247;
t292 = qJD(1) * t220;
t265 = t228 * t292;
t252 = qJD(3) - t265;
t195 = t229 * qJD(2) + t252;
t222 = cos(pkin(5));
t291 = qJD(1) * t222;
t246 = -t195 * t227 + t224 * t291;
t276 = qJDD(1) * t222;
t259 = t227 * t276;
t162 = qJDD(4) * pkin(8) - t246 * qJD(4) + t182 * t224 + t259;
t181 = t195 * t224 + t227 * t291;
t260 = t224 * t276;
t163 = -qJDD(4) * pkin(4) + qJD(4) * t181 - t182 * t227 + t260;
t170 = -qJD(4) * pkin(4) + t246;
t204 = pkin(4) * t224 - pkin(8) * t227 + qJ(3);
t270 = t225 * t292;
t184 = t204 * qJD(2) + t270;
t263 = t227 * t278;
t274 = qJDD(2) * t224;
t242 = t263 + t274;
t196 = qJDD(5) + t242;
t211 = qJD(2) * t224 + qJD(5);
t219 = sin(pkin(9));
t221 = cos(pkin(9));
t302 = t222 * t228;
t189 = t219 * t225 - t221 * t302;
t191 = t219 * t302 + t221 * t225;
t254 = g(1) * t191 + g(2) * t189;
t283 = qJD(4) * t227;
t316 = -(qJD(5) * t204 + t229 * t283) * t211 + t163 * t227 + (-qJD(4) * t170 - qJD(5) * t184 - t196 * t229 - t162) * t224 + t254;
t304 = t220 * t228;
t193 = t222 * t224 + t227 * t304;
t307 = t220 * t224;
t241 = g(1) * (t191 * t227 - t219 * t307) + g(2) * (t189 * t227 + t221 * t307) - g(3) * t193;
t255 = pkin(4) * t227 + pkin(8) * t224;
t315 = (pkin(8) * qJD(5) + t255 * qJD(2)) * t211 + t163 + t241;
t303 = t222 * t225;
t190 = t219 * t228 + t221 * t303;
t192 = -t219 * t303 + t221 * t228;
t253 = g(1) * t192 + g(2) * t190;
t295 = qJDD(1) - g(3);
t314 = -t295 * t306 + t253;
t290 = qJD(2) * qJ(3);
t203 = t270 + t290;
t313 = (-t203 + t270 - t290) * qJD(4) - qJDD(4) * t229;
t171 = qJD(4) * pkin(8) + t181;
t312 = -(t211 * t229 + t171) * qJD(5) - t253;
t311 = qJDD(2) * pkin(2);
t223 = sin(qJ(5));
t226 = cos(qJ(5));
t279 = t226 * qJD(4);
t280 = qJD(5) * t227;
t240 = -t223 * t280 - t224 * t279;
t271 = qJD(5) * t279 + t223 * qJDD(4) + t226 * t273;
t168 = t240 * qJD(2) + t271;
t310 = t168 * t223;
t288 = qJD(2) * t227;
t268 = t223 * t288;
t199 = t268 - t279;
t309 = t199 * t211;
t284 = qJD(4) * t223;
t201 = t226 * t288 + t284;
t308 = t201 * t211;
t305 = t220 * t227;
t301 = t223 * t225;
t300 = t224 * t229;
t299 = t225 * t226;
t298 = t225 * t227;
t297 = t226 * t211;
t296 = t227 * t229;
t218 = t227 ^ 2;
t294 = t224 ^ 2 - t218;
t230 = qJD(4) ^ 2;
t231 = qJD(2) ^ 2;
t293 = -t230 - t231;
t289 = qJD(2) * t203;
t287 = qJD(2) * t228;
t286 = qJD(4) * t199;
t285 = qJD(4) * t201;
t282 = qJD(5) * t223;
t281 = qJD(5) * t226;
t275 = qJDD(2) * qJ(3);
t267 = t211 * t284;
t266 = t211 * t279;
t262 = t225 * t277;
t256 = -t182 + t289;
t165 = t171 * t226 + t184 * t223;
t251 = t171 * t223 - t184 * t226;
t250 = -t226 * qJDD(4) + t317 * t223;
t194 = t222 * t227 - t224 * t304;
t249 = t223 * t228 + t224 * t299;
t248 = t224 * t301 - t226 * t228;
t245 = t223 * t196 + t211 * t281;
t244 = t226 * t196 - t211 * t282;
t197 = t255 * qJD(4) + qJD(3);
t237 = -t170 * t280 - t204 * t196 - t197 * t211;
t236 = -g(3) * t304 + t254 + t261;
t235 = qJDD(3) - t236;
t234 = -pkin(8) * t196 + (t170 - t246) * t211;
t183 = t262 + t275 + (qJD(3) + t265) * qJD(2);
t232 = -g(3) * t306 + t252 * qJD(2) - t229 * t230 + t183 - t253 + t275;
t215 = qJDD(4) * t227;
t198 = -qJD(2) * pkin(2) + t252;
t185 = t247 - t311;
t179 = t194 * t226 + t220 * t301;
t178 = -t194 * t223 + t220 * t299;
t177 = t194 * qJD(4) - t227 * t269;
t176 = -t193 * qJD(4) + t224 * t269;
t175 = -t189 * t224 + t221 * t305;
t173 = t191 * t224 + t219 * t305;
t169 = t201 * qJD(5) + t250;
t167 = t262 + t204 * qJDD(2) + (t197 + t265) * qJD(2);
t166 = t226 * t167;
t1 = [t295 * MDP(1) + (qJDD(1) * t222 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t177 - qJDD(4) * t193) * MDP(13) + (-qJD(4) * t176 - qJDD(4) * t194) * MDP(14) + ((-t176 * t223 - t194 * t281) * t211 + t178 * t196 + t177 * t199 + t193 * t169) * MDP(20) + (-(t176 * t226 - t194 * t282) * t211 - t179 * t196 + t177 * t201 + t193 * t168) * MDP(21) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t225 + t228 * t231) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t228 + t225 * t231) + ((-t185 + t289) * MDP(7) + (t224 * MDP(13) + MDP(14) * t227) * t231) * t228 + ((qJD(2) * t198 + t183) * MDP(7) + t242 * MDP(13) + t317 * MDP(14)) * t225 + ((-t225 * t282 + t226 * t287) * MDP(20) - (t223 * t287 + t225 * t281) * MDP(21)) * t211) * t220; qJDD(2) * MDP(2) + t236 * MDP(3) + t314 * MDP(4) + (t235 - 0.2e1 * t311) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t275 - t314) * MDP(6) + (t183 * qJ(3) + t203 * qJD(3) - t185 * pkin(2) - g(1) * (-pkin(2) * t191 + qJ(3) * t192) - g(2) * (-pkin(2) * t189 + qJ(3) * t190) + (-g(3) * (pkin(2) * t228 + qJ(3) * t225) + (-t198 * t225 - t203 * t228) * qJD(1)) * t220) * MDP(7) + (qJDD(2) * t218 - 0.2e1 * t224 * t263) * MDP(8) + 0.2e1 * (-t224 * t273 + t294 * t278) * MDP(9) + (-t224 * t230 + t215) * MDP(10) + (-qJDD(4) * t224 - t227 * t230) * MDP(11) + (t232 * t224 - t227 * t313) * MDP(13) + (t224 * t313 + t232 * t227) * MDP(14) + (t168 * t226 * t227 + t240 * t201) * MDP(15) + ((t199 * t226 + t201 * t223) * t224 * qJD(4) + (-t310 - t169 * t226 + (t199 * t223 - t201 * t226) * qJD(5)) * t227) * MDP(16) + ((t168 - t266) * t224 + (t244 + t285) * t227) * MDP(17) + ((-t169 + t267) * t224 + (-t245 - t286) * t227) * MDP(18) + (t196 * t224 + t211 * t283) * MDP(19) + (-t169 * t296 + t166 * t224 + (t199 * t300 - t227 * t251) * qJD(4) + (t224 * t312 - t237) * t226 + t316 * t223 + (-g(3) * t249 + (t199 * t298 + t248 * t211) * qJD(1)) * t220) * MDP(20) + (-t168 * t296 + (-t165 * t227 + t201 * t300) * qJD(4) + ((-t167 - t312) * t224 + t237) * t223 + t316 * t226 + (g(3) * t248 + (t201 * t298 + t249 * t211) * qJD(1)) * t220) * MDP(21); qJDD(2) * MDP(5) - t231 * MDP(6) + (t207 + t235 - t311) * MDP(7) + t215 * MDP(13) + (-t203 * MDP(7) + (-MDP(20) * t226 + MDP(21) * t223) * t211) * qJD(2) + (t293 * MDP(14) + (-t169 - t267) * MDP(20) + (-t168 - t266) * MDP(21)) * t227 + (t293 * MDP(13) - qJDD(4) * MDP(14) + (-t245 + t286) * MDP(20) + (-t244 + t285) * MDP(21)) * t224; MDP(10) * t273 - MDP(11) * t274 + qJDD(4) * MDP(12) + (-t256 * t227 - t241 - t260) * MDP(13) + (g(1) * t173 - g(2) * t175 + g(3) * t194 + t256 * t224 - t259) * MDP(14) + (t201 * t297 + t310) * MDP(15) + ((t168 - t309) * t226 + (-t169 - t308) * t223) * MDP(16) + ((-t201 * t227 + t224 * t297) * qJD(2) + t245) * MDP(17) + ((-t211 * t223 * t224 + t199 * t227) * qJD(2) + t244) * MDP(18) - t211 * MDP(19) * t288 + (-pkin(4) * t169 - t181 * t199 + t234 * t223 - t226 * t315 + t251 * t288) * MDP(20) + (-pkin(4) * t168 + t165 * t288 - t181 * t201 + t223 * t315 + t234 * t226) * MDP(21) + (t227 * t224 * MDP(8) - t294 * MDP(9)) * t231; t201 * t199 * MDP(15) + (-t199 ^ 2 + t201 ^ 2) * MDP(16) + (-t226 * t264 + t271 + t309) * MDP(17) + (-t250 + t308) * MDP(18) + t196 * MDP(19) + (-t223 * t162 + t166 + t165 * t211 - t170 * t201 - g(1) * (-t173 * t223 + t192 * t226) - g(2) * (t175 * t223 + t190 * t226) - g(3) * t178) * MDP(20) + (-t226 * t162 - t223 * t167 - t251 * t211 + t170 * t199 - g(1) * (-t173 * t226 - t192 * t223) - g(2) * (t175 * t226 - t190 * t223) + g(3) * t179) * MDP(21) + (-MDP(17) * t268 - t201 * MDP(18) - t165 * MDP(20) + t251 * MDP(21)) * qJD(5);];
tau = t1;
