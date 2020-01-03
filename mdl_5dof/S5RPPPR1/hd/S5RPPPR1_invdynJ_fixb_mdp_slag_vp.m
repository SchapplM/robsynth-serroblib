% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:41
% EndTime: 2020-01-03 11:20:46
% DurationCPUTime: 1.90s
% Computational Cost: add. (904->229), mult. (1830->328), div. (0->0), fcn. (1301->14), ass. (0->125)
t242 = cos(pkin(7));
t224 = -pkin(1) * t242 - pkin(2);
t285 = qJDD(1) * t224;
t211 = qJDD(3) + t285;
t236 = qJ(1) + pkin(7);
t226 = sin(t236);
t228 = cos(t236);
t310 = g(2) * t228 + g(3) * t226;
t311 = t211 + t310;
t239 = sin(pkin(7));
t219 = pkin(1) * t239 + qJ(3);
t207 = qJD(1) * qJD(3) + qJDD(1) * t219;
t238 = sin(pkin(8));
t309 = pkin(6) * t238;
t240 = cos(pkin(9));
t245 = cos(qJ(5));
t298 = t240 * t245;
t280 = t238 * t298;
t212 = qJDD(1) * t280;
t243 = sin(qJ(5));
t299 = t240 * t243;
t237 = sin(pkin(9));
t300 = t237 * t245;
t257 = t299 + t300;
t251 = t257 * qJD(5);
t283 = qJDD(1) * t243;
t275 = t237 * t283;
t165 = t212 + (-qJD(1) * t251 - t275) * t238;
t241 = cos(pkin(8));
t308 = t165 * t241;
t294 = qJD(1) * t238;
t278 = t237 * t294;
t271 = t243 * t278;
t210 = qJD(5) * t271;
t287 = qJD(1) * qJD(5);
t277 = t245 * t287;
t270 = t240 * t277;
t166 = -t210 + (t257 * qJDD(1) + t270) * t238;
t307 = t166 * t241;
t190 = t257 * t294;
t217 = t241 * qJD(1) - qJD(5);
t306 = t190 * t217;
t192 = qJD(1) * t280 - t271;
t305 = t192 * t217;
t304 = t219 * t241;
t303 = t226 * t241;
t302 = t228 * t241;
t301 = t237 * t243;
t194 = (-t238 * t301 + t280) * qJD(5);
t198 = t257 * t238;
t282 = t241 * qJDD(1);
t216 = -qJDD(5) + t282;
t297 = t194 * t217 + t198 * t216;
t206 = -pkin(3) * t241 - qJ(4) * t238 + t224;
t291 = qJD(4) * t238;
t174 = -qJD(1) * t291 + t206 * qJDD(1) + qJDD(3);
t187 = qJDD(2) * t238 + t207 * t241;
t161 = t237 * t174 + t240 * t187;
t188 = t206 * qJD(1) + qJD(3);
t215 = t219 * qJD(1);
t197 = qJD(2) * t238 + t215 * t241;
t164 = t237 * t188 + t240 * t197;
t173 = t237 * t206 + t240 * t304;
t232 = t238 ^ 2;
t295 = t241 ^ 2 + t232;
t293 = qJD(3) * t238;
t292 = qJD(3) * t241;
t290 = t216 * MDP(17);
t289 = t245 * MDP(18);
t284 = qJDD(1) * t238;
t281 = t240 * t309;
t246 = cos(qJ(1));
t279 = t246 * pkin(1) + t228 * pkin(2) + t226 * qJ(3);
t276 = t237 * t284;
t247 = qJD(1) ^ 2;
t274 = t295 * t247;
t160 = t240 * t174 - t187 * t237;
t255 = -pkin(4) * t241 - t281;
t157 = t255 * qJDD(1) + t160;
t158 = -pkin(6) * t276 + t161;
t273 = t245 * t157 - t243 * t158;
t163 = t240 * t188 - t197 * t237;
t272 = (-t237 ^ 2 - t240 ^ 2) * MDP(11);
t196 = qJD(2) * t241 - t238 * t215;
t186 = qJDD(2) * t241 - t238 * t207;
t268 = -g(2) * t226 + g(3) * t228;
t244 = sin(qJ(1));
t267 = -g(2) * t246 - g(3) * t244;
t266 = t244 * pkin(1) + t226 * pkin(2) - qJ(3) * t228;
t195 = qJD(4) - t196;
t265 = t243 * t157 + t245 * t158;
t159 = t255 * qJD(1) + t163;
t162 = -pkin(6) * t278 + t164;
t264 = t159 * t245 - t162 * t243;
t263 = -t159 * t243 - t162 * t245;
t201 = t240 * t206;
t167 = -t281 + t201 + (-t219 * t237 - pkin(4)) * t241;
t168 = -t237 * t309 + t173;
t262 = t167 * t245 - t168 * t243;
t261 = t167 * t243 + t168 * t245;
t260 = -t186 * t238 + t187 * t241;
t193 = t238 * t251;
t256 = -t298 + t301;
t199 = t256 * t238;
t259 = t193 * t217 + t199 * t216;
t258 = t196 * t238 - t197 * t241;
t181 = qJDD(4) - t186;
t172 = -t237 * t304 + t201;
t204 = -t237 * t292 - t240 * t291;
t254 = -t204 * qJD(1) - t172 * qJDD(1) - t160;
t205 = -t237 * t291 + t240 * t292;
t253 = t205 * qJD(1) + t173 * qJDD(1) + t161;
t252 = t217 * t256;
t250 = t285 + t311;
t235 = pkin(9) + qJ(5);
t227 = cos(t235);
t225 = sin(t235);
t202 = (pkin(4) * t237 + t219) * t238;
t185 = t225 * t226 + t227 * t302;
t184 = t225 * t302 - t226 * t227;
t183 = -t225 * t228 + t227 * t303;
t182 = -t225 * t303 - t227 * t228;
t176 = pkin(4) * t278 + t195;
t171 = pkin(4) * t276 + t181;
t1 = [qJDD(1) * MDP(1) + t267 * MDP(2) + (g(2) * t244 - g(3) * t246) * MDP(3) + (t267 + (t239 ^ 2 + t242 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t207 * t295 + t260 + t268) * MDP(7) + (-g(2) * t279 - g(3) * t266 - t258 * qJD(3) + t211 * t224 + t260 * t219) * MDP(8) + (t161 * t173 + t164 * t205 + t160 * t172 + t163 * t204 - g(2) * (pkin(3) * t302 + t279) - g(3) * (pkin(3) * t303 + t266)) * MDP(12) + (-t165 * t199 - t192 * t193) * MDP(13) + (-t165 * t198 + t166 * t199 + t190 * t193 - t192 * t194) * MDP(14) + (t259 - t308) * MDP(15) + (t297 + t307) * MDP(16) + (-g(2) * t185 - g(3) * t183 + t202 * t166 + t171 * t198 + t176 * t194 + t190 * t293 - t262 * t216) * MDP(18) + (g(2) * t184 - g(3) * t182 + t202 * t165 - t171 * t199 - t176 * t193 + t192 * t293 + t261 * t216) * MDP(19) + (t240 * MDP(10) + t237 * MDP(9)) * (t181 * t238 + t207 * t232 + t268) + ((t261 * qJD(5) - t204 * t245 + t205 * t243) * MDP(18) + (t262 * qJD(5) + t204 * t243 + t205 * t245) * MDP(19)) * t217 + (-t250 * MDP(5) + (-t240 * t310 + t254) * MDP(9) + (t237 * t310 + t253) * MDP(10) + t290 + (-t263 * qJD(5) - t273) * MDP(18) + (t264 * qJD(5) + t265) * MDP(19)) * t241 + (t250 * MDP(6) + (-t253 * t237 + t254 * t240 - t310) * MDP(11) + (-qJ(4) * t310 + t195 * qJD(3) + t181 * t219) * MDP(12)) * t238; (qJDD(2) - g(1)) * MDP(4) + (t186 * t241 + t187 * t238 - g(1)) * MDP(8) + (-t181 * t241 - g(1) + (-t160 * t237 + t161 * t240) * t238) * MDP(12) + (t297 - t307) * MDP(18) + (-t259 - t308) * MDP(19); -MDP(5) * t282 - MDP(7) * t274 + (t258 * qJD(1) + t311) * MDP(8) + (-t237 * t274 - t240 * t282) * MDP(9) + (t237 * t282 - t240 * t274) * MDP(10) + (t160 * t240 + t161 * t237 + (-t195 * t238 + (t163 * t237 - t164 * t240) * t241) * qJD(1) + t310) * MDP(12) + (t256 * t216 + t217 * t251 + (-t257 * t217 * t241 - t238 * t190) * qJD(1)) * MDP(18) + (t257 * t216 - qJD(5) * t252 + (-t238 * t192 + t241 * t252) * qJD(1)) * MDP(19) + (MDP(6) + t272) * t284; (g(1) * t241 + t181) * MDP(12) + (-t210 - t305) * MDP(18) + (t212 + t306) * MDP(19) + t247 * t232 * t272 + (t268 * MDP(12) + (MDP(10) * t237 - MDP(9) * t240) * t247 * t241 + ((t243 * MDP(18) + MDP(10)) * t240 + (-t243 * MDP(19) + MDP(9) + t289) * t237) * qJDD(1) + ((t163 * t240 + t164 * t237) * MDP(12) + (-t257 * MDP(19) + t240 * t289) * qJD(5)) * qJD(1)) * t238; t192 * t190 * MDP(13) + (-t190 ^ 2 + t192 ^ 2) * MDP(14) + (t212 - t306) * MDP(15) + (t210 - t305) * MDP(16) - t290 + (-g(2) * t182 - g(3) * t184 - t176 * t192 + t263 * t217 + t273) * MDP(18) + (g(2) * t183 - g(3) * t185 + t176 * t190 - t264 * t217 - t265) * MDP(19) + (t263 * MDP(18) - t264 * MDP(19)) * qJD(5) + ((-t237 * t277 - t287 * t299 - t275) * MDP(15) + (-qJDD(1) * t300 - t240 * t283 - t270) * MDP(16) + (MDP(18) * t225 + MDP(19) * t227) * g(1)) * t238;];
tau = t1;
