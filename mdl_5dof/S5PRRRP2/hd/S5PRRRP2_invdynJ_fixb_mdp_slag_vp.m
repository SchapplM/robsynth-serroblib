% Calculate vector of inverse dynamics joint torques for
% S5PRRRP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:42:04
% EndTime: 2019-12-05 16:42:08
% DurationCPUTime: 1.20s
% Computational Cost: add. (970->214), mult. (1343->240), div. (0->0), fcn. (662->8), ass. (0->110)
t211 = qJDD(2) + qJDD(3);
t217 = sin(qJ(3));
t255 = qJDD(2) * t217;
t219 = cos(qJ(3));
t262 = qJD(3) * t219;
t292 = pkin(7) * t211 + (qJD(2) * t262 + t255) * pkin(2) + qJD(1) * qJD(4);
t218 = cos(qJ(4));
t213 = qJD(2) + qJD(3);
t281 = pkin(2) * qJD(2);
t250 = t217 * t281;
t176 = pkin(7) * t213 + t250;
t216 = sin(qJ(4));
t278 = t176 * t216;
t163 = qJD(1) * t218 - t278;
t291 = qJD(5) - t163;
t212 = pkin(8) + qJ(2);
t209 = qJ(3) + t212;
t197 = cos(t209);
t194 = g(1) * t197;
t214 = t216 ^ 2;
t215 = t218 ^ 2;
t265 = t214 + t215;
t290 = t213 * t265;
t248 = t216 * qJDD(1) + t292 * t218;
t254 = qJDD(4) * qJ(5);
t150 = t254 + (qJD(5) - t278) * qJD(4) + t248;
t258 = qJD(4) * t218;
t240 = -t218 * qJDD(1) + t176 * t258 + t292 * t216;
t279 = (qJDD(4) * pkin(4));
t288 = qJDD(5) - t279;
t151 = t240 + t288;
t289 = t150 * t218 + t151 * t216;
t156 = -qJD(4) * pkin(4) + t291;
t164 = qJD(1) * t216 + t176 * t218;
t157 = qJD(4) * qJ(5) + t164;
t196 = sin(t209);
t268 = g(2) * t196 + t194;
t287 = pkin(2) * t219;
t286 = pkin(3) * t211;
t285 = pkin(3) * t213;
t284 = pkin(4) * t218;
t220 = qJD(4) ^ 2;
t283 = pkin(7) * t220;
t193 = g(1) * t196;
t282 = g(2) * t197;
t280 = pkin(7) * qJDD(4);
t235 = qJ(5) * t216 + t284;
t179 = -pkin(3) - t235;
t277 = t179 * t211;
t276 = t179 * t213;
t275 = t196 * t216;
t274 = t197 * t216;
t200 = pkin(2) * t217 + pkin(7);
t273 = t200 * t220;
t272 = t211 * t216;
t271 = t213 * t216;
t270 = t216 * t218;
t269 = g(1) * t275 - g(2) * t274;
t267 = -qJD(3) * t250 + qJDD(2) * t287;
t266 = t214 - t215;
t264 = qJD(2) * t219;
t263 = qJD(3) * t217;
t260 = qJD(4) * t213;
t259 = qJD(4) * t216;
t257 = qJD(5) * t216;
t253 = qJDD(4) * t200;
t252 = pkin(2) * t264;
t251 = pkin(2) * t262;
t210 = t213 ^ 2;
t249 = t210 * t270;
t186 = t218 * t193;
t241 = t213 * t250;
t247 = t218 * t241 + t252 * t259 + t186;
t246 = t213 * t263;
t165 = -t267 - t286;
t245 = -t165 - t282;
t244 = t265 * t211;
t242 = t163 + t278;
t177 = -t252 - t285;
t239 = t165 * t216 + t177 * t258 - t269;
t238 = t283 - t286;
t205 = sin(t212);
t206 = cos(t212);
t236 = g(1) * t205 - g(2) * t206;
t234 = pkin(4) * t216 - qJ(5) * t218;
t180 = qJDD(4) * t216 + t218 * t220;
t181 = qJDD(4) * t218 - t216 * t220;
t233 = 0.2e1 * (t211 * t270 - t266 * t260) * MDP(9) + (t211 * t214 + 0.2e1 * t258 * t271) * MDP(8) + t180 * MDP(10) + t181 * MDP(11) + t211 * MDP(5);
t232 = t156 * t216 + t157 * t218;
t231 = t193 + t267 - t282;
t147 = (t234 * qJD(4) - t257) * t213 + t277 - t267;
t230 = -t147 - t277 - t283;
t172 = t179 - t287;
t229 = t172 * t213 - t251;
t167 = pkin(4) * t259 - qJ(5) * t258 - t257;
t228 = g(1) * t274 + g(2) * t275 - g(3) * t218 - t240;
t158 = pkin(2) * t263 + t167;
t227 = -t158 * t213 - t172 * t211 - t147 - t273;
t226 = t156 * t258 - t157 * t259 - t268 + t289;
t201 = -pkin(3) - t287;
t225 = pkin(2) * t246 + t201 * t211 + t273;
t224 = qJD(4) * t164 + t228;
t223 = -t253 + (t201 * t213 - t251) * qJD(4);
t222 = (t156 * t218 - t157 * t216) * qJD(4) + t289;
t221 = -pkin(7) * t194 - g(2) * (t196 * pkin(7) + qJ(5) * t274 + (pkin(3) + t284) * t197) - t179 * t193;
t169 = t177 * t259;
t168 = t234 * t213;
t155 = -t252 + t276;
t152 = t155 * t259;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t232 * qJD(4) + t150 * t216 - t151 * t218 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t181 + (-MDP(14) + MDP(17)) * t180; qJDD(2) * MDP(2) + t236 * MDP(3) + (g(1) * t206 + g(2) * t205) * MDP(4) + ((t211 * t219 - t246) * pkin(2) + t231) * MDP(6) + (((-qJDD(2) - t211) * t217 + (-qJD(2) - t213) * t262) * pkin(2) + t268) * MDP(7) + (t169 + t186 + t223 * t216 + (-t225 + t245) * t218) * MDP(13) + (t225 * t216 + t223 * t218 + t239) * MDP(14) + (t152 + t186 + (t229 * qJD(4) - t253) * t216 + (t227 - t282) * t218) * MDP(15) + (t200 * t244 + t251 * t290 + t226) * MDP(16) + ((t253 + (-t155 - t229) * qJD(4)) * t218 + t227 * t216 + t269) * MDP(17) + (t147 * t172 + t155 * t158 + (t232 * t262 + t236) * pkin(2) + t222 * t200 + t221) * MDP(18) + t233; (t231 + t241) * MDP(6) + ((-t255 + (-qJD(3) + t213) * t264) * pkin(2) + t268) * MDP(7) + (t169 + (-pkin(3) * t260 - t280) * t216 + (-t238 + t245) * t218 + t247) * MDP(13) + ((-t280 + (t252 - t285) * qJD(4)) * t218 + (t238 - t241) * t216 + t239) * MDP(14) + (t152 + (t179 * t260 - t280) * t216 + (-t167 * t213 + t230 - t282) * t218 + t247) * MDP(15) + (pkin(7) * t244 - t252 * t290 + t226) * MDP(16) + ((t280 + (-t155 - t252 - t276) * qJD(4)) * t218 + ((-t167 + t250) * t213 + t230) * t216 + t269) * MDP(17) + (t147 * t179 + t155 * t167 + (-t155 * t217 - t232 * t219) * t281 + t222 * pkin(7) + t221) * MDP(18) + t233; -MDP(8) * t249 + t266 * MDP(9) * t210 + MDP(10) * t272 + qJDD(4) * MDP(12) + (-t177 * t271 + t224) * MDP(13) + (g(3) * t216 + t242 * qJD(4) + (-t177 * t213 + t268) * t218 - t248) * MDP(14) + ((2 * t279) - qJDD(5) + (-t155 * t216 + t168 * t218) * t213 + t224) * MDP(15) + (0.2e1 * t254 + (t168 * t213 - g(3)) * t216 + (t155 * t213 - t268) * t218 + (0.2e1 * qJD(5) - t242) * qJD(4) + t248) * MDP(17) + (-t151 * pkin(4) - g(3) * t235 + t150 * qJ(5) - t155 * t168 - t156 * t164 + t291 * t157 + t268 * t234) * MDP(18) + (t218 * MDP(11) - t234 * MDP(16)) * t211; (-qJDD(4) - t249) * MDP(15) + MDP(16) * t272 + (-t210 * t214 - t220) * MDP(17) + (-qJD(4) * t157 + t155 * t271 - t228 + t288) * MDP(18);];
tau = t1;
