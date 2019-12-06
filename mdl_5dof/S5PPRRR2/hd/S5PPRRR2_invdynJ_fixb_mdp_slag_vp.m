% Calculate vector of inverse dynamics joint torques for
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:56
% EndTime: 2019-12-05 15:14:59
% DurationCPUTime: 1.54s
% Computational Cost: add. (737->184), mult. (1567->257), div. (0->0), fcn. (1225->14), ass. (0->104)
t224 = sin(pkin(9));
t226 = cos(pkin(9));
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t303 = -t224 * t230 + t226 * t233;
t178 = t303 * qJD(1);
t232 = cos(qJ(4));
t211 = -pkin(4) * t232 - pkin(3);
t171 = t211 * qJD(3) - t178;
t305 = t171 + t178;
t186 = t224 * t233 + t226 * t230;
t181 = t186 * qJD(3);
t247 = -qJD(1) * t181 + t303 * qJDD(1);
t293 = qJDD(3) * pkin(3);
t165 = -t247 - t293;
t234 = qJD(4) ^ 2;
t179 = t186 * qJD(1);
t219 = pkin(9) + qJ(3);
t212 = sin(t219);
t213 = cos(t219);
t225 = sin(pkin(8));
t227 = cos(pkin(8));
t259 = g(1) * t227 + g(2) * t225;
t242 = -g(3) * t213 + t259 * t212;
t239 = qJD(3) * t179 + t242;
t304 = -pkin(6) * t234 - t165 + t239 + t293;
t220 = qJD(4) + qJD(5);
t298 = pkin(6) + pkin(7);
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t267 = qJDD(3) * t232;
t229 = sin(qJ(4));
t268 = qJDD(3) * t229;
t188 = t228 * t232 + t229 * t231;
t302 = t220 * t188;
t160 = qJD(3) * t302 + t228 * t268 - t231 * t267;
t187 = t228 * t229 - t231 * t232;
t248 = t220 * t187;
t262 = t298 * qJD(3) + t179;
t167 = t232 * qJD(2) - t262 * t229;
t301 = t220 * t232;
t300 = t186 * qJDD(1);
t168 = qJD(2) * t229 + t262 * t232;
t296 = g(3) * t212;
t294 = qJD(3) * pkin(3);
t292 = t168 * t231;
t218 = qJDD(4) + qJDD(5);
t291 = t187 * t218;
t290 = t188 * t218;
t223 = qJ(4) + qJ(5);
t216 = sin(t223);
t289 = t216 * t225;
t288 = t216 * t227;
t217 = cos(t223);
t287 = t217 * t225;
t286 = t217 * t227;
t221 = t229 ^ 2;
t281 = -t232 ^ 2 + t221;
t278 = qJD(3) * t303;
t277 = qJD(3) * t229;
t276 = qJD(3) * t232;
t275 = qJD(4) * t229;
t273 = qJD(5) * t228;
t272 = qJD(5) * t231;
t269 = qJD(3) * qJD(4);
t266 = qJD(4) * t298;
t265 = t228 * t277;
t264 = t231 * t276;
t263 = t232 * t269;
t164 = qJDD(3) * pkin(6) + qJD(1) * t278 + t300;
t261 = pkin(7) * qJDD(3) + t164;
t260 = pkin(4) * t275 - t179;
t258 = -g(1) * t225 + g(2) * t227;
t166 = qJD(4) * pkin(4) + t167;
t255 = -t166 * t228 - t292;
t254 = -t220 * t248 + t290;
t195 = t298 * t229;
t196 = t298 * t232;
t253 = -t195 * t231 - t196 * t228;
t252 = -t195 * t228 + t196 * t231;
t251 = -qJD(3) * t181 + qJDD(3) * t303;
t250 = -qJDD(2) - t258;
t159 = qJD(5) * t264 - t220 * t265 + t228 * t267 + (t263 + t268) * t231;
t182 = -t264 + t265;
t184 = -t228 * t276 - t231 * t277;
t246 = -t184 * t182 * MDP(13) + (t182 * t220 + t159) * MDP(15) + (-t184 * t220 - t160) * MDP(16) + (-t182 ^ 2 + t184 ^ 2) * MDP(14) + t218 * MDP(17);
t245 = t186 * t234 - t251;
t244 = -0.2e1 * t278 * qJD(4) - qJDD(4) * t186;
t243 = t259 * t213 + t296;
t176 = -t178 - t294;
t240 = -pkin(6) * qJDD(4) + (t176 + t178 - t294) * qJD(4);
t238 = -t176 * qJD(3) - t164 + t243;
t214 = t232 * qJDD(2);
t154 = qJDD(4) * pkin(4) - qJD(4) * t168 - t261 * t229 + t214;
t237 = -g(1) * (-t213 * t286 - t289) - g(2) * (-t213 * t287 + t288) + t171 * t182 + t168 * t273 + t217 * t296 + (-t168 * t220 - t154) * t228;
t155 = qJD(4) * t167 + t229 * qJDD(2) + t261 * t232;
t236 = -g(1) * (-t213 * t288 + t287) - g(2) * (-t213 * t289 - t286) + t255 * qJD(5) + t231 * t154 - t228 * t155 + t171 * t184 + t216 * t296;
t194 = qJDD(4) * t232 - t229 * t234;
t193 = qJDD(4) * t229 + t232 * t234;
t190 = t232 * t266;
t189 = t229 * t266;
t161 = (t229 * t269 - t267) * pkin(4) + t165;
t158 = -t220 * t302 - t291;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t224 ^ 2 + t226 ^ 2) * qJDD(1)) * MDP(2) + t251 * MDP(4) + (-qJD(3) * t278 - qJDD(3) * t186) * MDP(5) + (t244 * t229 - t245 * t232) * MDP(11) + (t245 * t229 + t244 * t232) * MDP(12) + (-t303 * t160 + t181 * t182 - t278 * t302 + ((t228 * t275 + t229 * t273 - t231 * t301) * t220 - t290) * t186) * MDP(18) + (-t303 * t159 - t181 * t184 + t278 * t248 + (-(-t228 * t301 - t229 * t272 - t231 * t275) * t220 + t291) * t186) * MDP(19); t194 * MDP(11) - t193 * MDP(12) + t158 * MDP(18) - t254 * MDP(19) - t250 * MDP(2); qJDD(3) * MDP(3) + (t239 + t247) * MDP(4) + (-t300 + t243) * MDP(5) + (qJDD(3) * t221 + 0.2e1 * t229 * t263) * MDP(6) + 0.2e1 * (t229 * t267 - t281 * t269) * MDP(7) + t193 * MDP(8) + t194 * MDP(9) + (t240 * t229 + t304 * t232) * MDP(11) + (-t304 * t229 + t240 * t232) * MDP(12) + (t159 * t188 + t184 * t248) * MDP(13) + (-t159 * t187 - t160 * t188 + t182 * t248 + t184 * t302) * MDP(14) + t254 * MDP(15) + t158 * MDP(16) + ((-t252 * qJD(5) + t189 * t228 - t190 * t231) * t220 + t253 * t218 + t211 * t160 + t161 * t187 + t260 * t182 + t242 * t217 + t305 * t302) * MDP(18) + (-(t253 * qJD(5) - t189 * t231 - t190 * t228) * t220 - t252 * t218 + t211 * t159 + t161 * t188 - t260 * t184 - t242 * t216 - t305 * t248) * MDP(19); MDP(8) * t268 + MDP(9) * t267 + qJDD(4) * MDP(10) + (t238 * t229 + t258 * t232 + t214) * MDP(11) + (t250 * t229 + t238 * t232) * MDP(12) + (-(-t167 * t228 - t292) * t220 + (-t182 * t277 + t231 * t218 - t220 * t273) * pkin(4) + t236) * MDP(18) + ((-qJD(5) * t166 + t167 * t220 - t155) * t231 + (t184 * t277 - t228 * t218 - t220 * t272) * pkin(4) + t237) * MDP(19) + t246 + (-t229 * t232 * MDP(6) + t281 * MDP(7)) * qJD(3) ^ 2; (-t255 * t220 + t236) * MDP(18) + ((-t155 + (-qJD(5) + t220) * t166) * t231 + t237) * MDP(19) + t246;];
tau = t1;
