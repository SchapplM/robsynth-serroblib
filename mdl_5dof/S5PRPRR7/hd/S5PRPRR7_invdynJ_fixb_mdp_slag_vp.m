% Calculate vector of inverse dynamics joint torques for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:52
% EndTime: 2019-12-05 16:00:57
% DurationCPUTime: 1.76s
% Computational Cost: add. (685->223), mult. (1293->288), div. (0->0), fcn. (878->10), ass. (0->114)
t222 = cos(qJ(2));
t209 = g(3) * t222;
t219 = sin(qJ(2));
t215 = sin(pkin(8));
t216 = cos(pkin(8));
t252 = g(1) * t216 + g(2) * t215;
t300 = t252 * t219 - t209;
t217 = sin(qJ(5));
t218 = sin(qJ(4));
t220 = cos(qJ(5));
t221 = cos(qJ(4));
t186 = t217 * t221 + t218 * t220;
t211 = qJD(4) + qJD(5);
t236 = t186 * t211;
t305 = qJD(2) * t236;
t210 = qJDD(4) + qJDD(5);
t243 = t217 * t218 - t220 * t221;
t304 = t243 * t210;
t261 = qJDD(2) * t221;
t267 = qJD(2) * qJD(4);
t303 = t218 * t267 - t261;
t239 = t252 * t222;
t286 = qJDD(1) - g(3);
t302 = -t286 * t219 + t239;
t280 = qJD(2) * qJ(3);
t282 = qJD(1) * t219;
t193 = t280 + t282;
t270 = t193 * qJD(2);
t301 = -t270 - t300;
t274 = qJD(4) * t221;
t297 = qJD(5) * t221 + t274;
t223 = -pkin(2) - pkin(6);
t296 = (t193 + t280 - t282) * qJD(4) + qJDD(4) * t223;
t293 = g(3) * t219;
t292 = pkin(7) - t223;
t291 = qJDD(2) * pkin(2);
t281 = qJD(1) * t222;
t254 = qJD(3) - t281;
t185 = qJD(2) * t223 + t254;
t277 = qJD(2) * t218;
t168 = -pkin(7) * t277 + t185 * t218;
t290 = t168 * t220;
t289 = t186 * t210;
t288 = t215 * t219;
t287 = t216 * t219;
t285 = -t211 * t236 - t304;
t213 = t221 ^ 2;
t284 = t218 ^ 2 - t213;
t224 = qJD(4) ^ 2;
t225 = qJD(2) ^ 2;
t283 = t224 + t225;
t276 = qJD(2) * t221;
t179 = t217 * t277 - t220 * t276;
t279 = qJD(2) * t179;
t180 = t186 * qJD(2);
t278 = qJD(2) * t180;
t275 = qJD(4) * t218;
t273 = qJD(5) * t217;
t272 = qJD(5) * t220;
t269 = qJD(1) * qJD(2);
t268 = qJD(2) * qJD(3);
t266 = qJDD(1) * t219;
t265 = qJDD(1) * t222;
t264 = qJDD(2) * qJ(3);
t263 = qJDD(2) * t218;
t262 = qJDD(2) * t219;
t260 = qJDD(4) * t218;
t258 = t218 * t273;
t189 = t292 * t221;
t256 = t221 * t267;
t202 = pkin(4) * t218 + qJ(3);
t255 = t211 * t221;
t253 = qJDD(3) - t265;
t169 = -pkin(7) * t276 + t185 * t221;
t197 = pkin(4) * t274 + qJD(3);
t251 = g(1) * t215 - g(2) * t216;
t232 = -t217 * t275 - t258;
t164 = t220 * t255 + t232;
t248 = -t164 * t211 - t289;
t166 = qJD(4) * pkin(4) + t169;
t247 = -t166 * t217 - t290;
t188 = t292 * t218;
t246 = -t188 * t220 - t189 * t217;
t245 = t188 * t217 - t189 * t220;
t244 = (-qJD(2) * pkin(2) + t254) * t219 + t193 * t222;
t203 = t219 * t269;
t242 = t203 + t253;
t241 = -qJD(2) * t258 - t217 * t303;
t235 = t243 * t211;
t234 = -qJDD(4) * t222 + 0.2e1 * t219 * t267;
t233 = t222 * t283 + t262;
t158 = -t217 * t263 + t220 * t261 - t305;
t231 = -t179 * t180 * MDP(15) + (t180 * t211 + t158) * MDP(17) + (-t179 * t211 + (-t211 * t276 - t263) * t220 - t241) * MDP(18) + (t179 ^ 2 - t180 ^ 2) * MDP(16) + t210 * MDP(19);
t230 = -t239 - t293;
t229 = t253 - t300;
t174 = qJDD(2) * t223 + t242;
t167 = t221 * t174;
t156 = qJDD(4) * pkin(4) + pkin(7) * t303 - t185 * t275 + t167;
t182 = qJD(2) * t202 + t282;
t214 = qJ(4) + qJ(5);
t207 = sin(t214);
t208 = cos(t214);
t228 = t168 * t273 + (-t168 * t211 - t156) * t217 - g(1) * (-t207 * t287 - t208 * t215) - g(2) * (-t207 * t288 + t208 * t216) + t182 * t180 - t207 * t209;
t157 = t185 * t274 + t218 * t174 + (-t256 - t263) * pkin(7);
t227 = -g(1) * (-t207 * t215 + t208 * t287) - g(2) * (t207 * t216 + t208 * t288) + t247 * qJD(5) + t220 * t156 - t217 * t157 + t182 * t179 + t208 * t209;
t175 = t264 + t266 + (qJD(3) + t281) * qJD(2);
t226 = -t293 + t264 + t268 - t223 * t224 + t175 + (-t252 - t269) * t222;
t206 = qJDD(4) * t221;
t184 = qJD(4) * t189;
t183 = t292 * t275;
t177 = t242 - t291;
t162 = t266 + t202 * qJDD(2) + (t197 + t281) * qJD(2);
t159 = (qJD(2) * t255 + t263) * t220 + t241;
t1 = [t286 * MDP(1) + (qJD(2) * t244 + t175 * t219 - t177 * t222 - g(3)) * MDP(7) + (t218 * t233 + t221 * t234) * MDP(13) + (-t218 * t234 + t221 * t233) * MDP(14) + ((-qJD(2) * t235 + t159) * t219 + ((t217 * t297 + t218 * t272 + t220 * t275) * t211 + t304 + t278) * t222) * MDP(20) + ((t158 - t305) * t219 + (-(-t220 * t297 - t232) * t211 + t289 - t279) * t222) * MDP(21) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t222 + t219 * t225) + (-MDP(4) + MDP(6)) * (t222 * t225 + t262); qJDD(2) * MDP(2) + (t265 + t300) * MDP(3) + t302 * MDP(4) + (t229 - 0.2e1 * t291) * MDP(5) + (0.2e1 * t264 + 0.2e1 * t268 - t302) * MDP(6) + (t175 * qJ(3) + t193 * qJD(3) - t177 * pkin(2) - g(3) * (pkin(2) * t222 + qJ(3) * t219) - t244 * qJD(1) + t252 * (pkin(2) * t219 - qJ(3) * t222)) * MDP(7) + (qJDD(2) * t213 - 0.2e1 * t218 * t256) * MDP(8) + 0.2e1 * (-t218 * t261 + t267 * t284) * MDP(9) + (-t218 * t224 + t206) * MDP(10) + (-t221 * t224 - t260) * MDP(11) + (t226 * t218 + t221 * t296) * MDP(13) + (-t218 * t296 + t226 * t221) * MDP(14) + (-t158 * t243 + t179 * t236) * MDP(15) + (-t158 * t186 + t159 * t243 + t164 * t179 + t180 * t236) * MDP(16) + t285 * MDP(17) + t248 * MDP(18) + ((-qJD(5) * t246 + t183 * t220 + t184 * t217) * t211 + t245 * t210 + t197 * t180 + t202 * t159 + t162 * t186 + t182 * t164 + t230 * t207 + (-t222 * t180 + t219 * t235) * qJD(1)) * MDP(20) + (-(qJD(5) * t245 + t183 * t217 - t184 * t220) * t211 - t246 * t210 - t197 * t179 + t202 * t158 - t162 * t243 - t182 * t236 + t230 * t208 + (t222 * t179 + t219 * t236) * qJD(1)) * MDP(21); qJDD(2) * MDP(5) - t225 * MDP(6) + (t203 + t229 - t270 - t291) * MDP(7) + (-t218 * t283 + t206) * MDP(13) + (-t221 * t283 - t260) * MDP(14) + (-t278 + t285) * MDP(20) + (t248 + t279) * MDP(21); MDP(10) * t261 - MDP(11) * t263 + qJDD(4) * MDP(12) + (t251 * t218 + t221 * t301 + t167) * MDP(13) + (t251 * t221 + (-t174 - t301) * t218) * MDP(14) + (-(-t169 * t217 - t290) * t211 + (-t180 * t276 + t220 * t210 - t211 * t273) * pkin(4) + t227) * MDP(20) + ((-qJD(5) * t166 + t169 * t211 - t157) * t220 + (t179 * t276 - t217 * t210 - t211 * t272) * pkin(4) + t228) * MDP(21) + t231 + (MDP(8) * t218 * t221 - MDP(9) * t284) * t225; (-t211 * t247 + t227) * MDP(20) + ((-t157 + (-qJD(5) + t211) * t166) * t220 + t228) * MDP(21) + t231;];
tau = t1;
