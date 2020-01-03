% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:16
% EndTime: 2020-01-03 11:59:18
% DurationCPUTime: 1.12s
% Computational Cost: add. (1101->209), mult. (1761->273), div. (0->0), fcn. (978->12), ass. (0->120)
t301 = qJ(5) + pkin(7);
t241 = cos(qJ(2));
t299 = pkin(1) * t241;
t217 = qJDD(1) * t299;
t229 = qJDD(1) + qJDD(2);
t238 = sin(qJ(2));
t293 = pkin(1) * qJD(2);
t272 = qJD(1) * t293;
t177 = pkin(2) * t229 - t238 * t272 + t217;
t234 = sin(pkin(8));
t235 = cos(pkin(8));
t276 = qJDD(1) * t238;
t300 = pkin(1) * t276 + t241 * t272;
t161 = t235 * t177 - t300 * t234;
t233 = qJ(1) + qJ(2);
t220 = pkin(8) + t233;
t208 = sin(t220);
t209 = cos(t220);
t258 = g(2) * t209 + g(3) * t208;
t250 = t229 * pkin(3) + t161 - t258;
t257 = g(2) * t208 - g(3) * t209;
t298 = pkin(2) * t235;
t240 = cos(qJ(4));
t297 = g(1) * t240;
t292 = qJD(4) * pkin(4);
t237 = sin(qJ(4));
t291 = t229 * t237;
t290 = t234 * t238;
t289 = t235 * t238;
t288 = t237 * t240;
t216 = pkin(2) + t299;
t283 = pkin(1) * t289 + t234 * t216;
t179 = pkin(7) + t283;
t287 = -qJ(5) - t179;
t210 = pkin(2) * t234 + pkin(7);
t286 = -qJ(5) - t210;
t285 = qJDD(3) - g(1);
t230 = qJD(1) + qJD(2);
t278 = qJD(1) * t241;
t273 = pkin(1) * t278;
t191 = pkin(2) * t230 + t273;
t279 = qJD(1) * t238;
t274 = pkin(1) * t279;
t204 = t235 * t274;
t172 = t234 * t191 + t204;
t267 = t301 * t230 + t172;
t159 = t240 * qJD(3) - t267 * t237;
t156 = t159 + t292;
t284 = t156 - t159;
t231 = t237 ^ 2;
t232 = t240 ^ 2;
t282 = -t231 - t232;
t281 = t231 - t232;
t203 = t234 * t274;
t171 = t191 * t235 - t203;
t215 = pkin(4) * t240 + pkin(3);
t163 = -t215 * t230 + qJD(5) - t171;
t280 = MDP(16) * t163;
t277 = qJD(4) * t237;
t275 = qJDD(4) * t179;
t223 = sin(t233);
t213 = pkin(2) * t223;
t270 = t208 * t215 - t209 * t301 + t213;
t162 = t234 * t177 + t300 * t235;
t269 = t230 * t277;
t224 = cos(t233);
t268 = g(2) * t223 - g(3) * t224;
t265 = -pkin(1) * t290 + t216 * t235;
t178 = -pkin(3) - t265;
t183 = (t235 * t241 - t290) * t293;
t266 = t178 * t230 - t183;
t264 = qJD(4) * t287;
t263 = qJD(4) * t286;
t262 = qJD(1) * (-qJD(2) + t230);
t167 = -pkin(3) * t230 - t171;
t260 = t167 * qJD(4) * t240 - t250 * t237;
t158 = pkin(7) * t229 + t162;
t244 = qJ(5) * t229 + qJD(3) * qJD(4) + qJD(5) * t230 + t158;
t251 = t267 * qJD(4);
t151 = (qJDD(3) - t251) * t237 + t244 * t240;
t259 = t151 * t240 - t257;
t256 = -g(2) * t224 - g(3) * t223;
t214 = pkin(2) * t224;
t255 = t208 * t301 + t209 * t215 + t214;
t243 = qJD(4) ^ 2;
t192 = qJDD(4) * t237 + t240 * t243;
t193 = qJDD(4) * t240 - t237 * t243;
t254 = 0.2e1 * (-t281 * t230 * qJD(4) + t229 * t288) * MDP(9) + (t229 * t231 + 0.2e1 * t240 * t269) * MDP(8) + t192 * MDP(10) + t193 * MDP(11) + t229 * MDP(4);
t160 = t237 * qJD(3) + t267 * t240;
t253 = t156 * t237 - t160 * t240;
t252 = t234 * t241 + t289;
t249 = t217 + t256;
t181 = t252 * t293;
t248 = t178 * t229 + t179 * t243 + t181 * t230;
t180 = t234 * t273 + t204;
t211 = -pkin(3) - t298;
t247 = -t180 * t230 + t210 * t243 + t211 * t229;
t246 = -t167 * t230 - t158 + t257;
t182 = t235 * t273 - t203;
t245 = -qJDD(4) * t210 + (t211 * t230 + t182) * qJD(4);
t152 = pkin(4) * t269 - t215 * t229 + qJDD(5) - t161;
t242 = cos(qJ(1));
t239 = sin(qJ(1));
t228 = t230 ^ 2;
t227 = t242 * pkin(1);
t226 = t239 * pkin(1);
t225 = t240 * qJ(5);
t221 = t240 * qJD(5);
t219 = t240 * qJDD(3);
t186 = t210 * t240 + t225;
t185 = t286 * t237;
t175 = -qJD(5) * t237 + t240 * t263;
t174 = t237 * t263 + t221;
t170 = t179 * t240 + t225;
t169 = t287 * t237;
t164 = t167 * t277;
t154 = (-qJD(5) - t183) * t237 + t240 * t264;
t153 = t183 * t240 + t237 * t264 + t221;
t150 = qJDD(4) * pkin(4) - t244 * t237 - t240 * t251 + t219;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t242 - g(3) * t239) * MDP(2) + (g(2) * t239 - g(3) * t242) * MDP(3) + t249 * MDP(5) + t268 * MDP(6) + (t162 * t283 + t172 * t183 + t161 * t265 - t171 * t181 - g(2) * (t214 + t227) - g(3) * (t213 + t226)) * MDP(7) + t164 * MDP(13) + t260 * MDP(14) + t259 * MDP(15) + (t151 * t170 + t160 * t153 + t150 * t169 + t156 * t154 + t152 * t178 - g(2) * (t227 + t255) - g(3) * (t226 + t270)) * MDP(16) + (t229 * t241 * MDP(5) + (-qJDD(1) - t229) * MDP(6) * t238 + ((-t230 * t238 - t279) * MDP(5) + (-t230 * t241 - t278) * MDP(6) + t252 * t280) * qJD(2)) * pkin(1) + (-MDP(13) * t275 + t248 * MDP(14) + (-t154 * t230 - t169 * t229 - t150) * MDP(15) + (t266 * MDP(13) + (-t170 * t230 - t160) * MDP(15) + pkin(4) * t280) * qJD(4)) * t237 + ((-t248 + t250) * MDP(13) - MDP(14) * t275 + (t153 * t230 + t170 * t229) * MDP(15) - t152 * pkin(4) * MDP(16) + (t266 * MDP(14) + (-t169 * t230 - t156) * MDP(15)) * qJD(4)) * t240 + t254; (t238 * pkin(1) * t262 + t249) * MDP(5) + ((t241 * t262 - t276) * pkin(1) + t268) * MDP(6) + (t171 * t180 - t172 * t182 + (t161 * t235 + t162 * t234 + t256) * pkin(2)) * MDP(7) + (t164 + t245 * t237 + (-t247 + t250) * t240) * MDP(13) + (t247 * t237 + t245 * t240 + t260) * MDP(14) + ((-qJD(4) * t156 + t186 * t229) * t240 + (-qJD(4) * t160 - t185 * t229 - t150) * t237 + (t174 * t240 - t175 * t237 + t282 * t182 + (-t185 * t240 - t186 * t237) * qJD(4)) * t230 + t259) * MDP(15) + (t151 * t186 + t150 * t185 + t152 * (-t215 - t298) - g(2) * t255 - g(3) * t270 + (pkin(4) * t277 - t180) * t163 + (-t182 * t240 + t174) * t160 + (t182 * t237 + t175) * t156) * MDP(16) + t254; t285 * MDP(7) + t193 * MDP(13) - t192 * MDP(14) + (-t253 * qJD(4) + t150 * t240 + t151 * t237 - g(1)) * MDP(16); MDP(10) * t291 + t240 * t229 * MDP(11) + qJDD(4) * MDP(12) + (t246 * t237 + t219 - t297) * MDP(13) + (-t285 * t237 + t246 * t240) * MDP(14) + (-pkin(4) * t291 + (t284 - t292) * t240 * t230) * MDP(15) + (t284 * t160 + (-t297 + t150 + (-t163 * t230 + t257) * t237) * pkin(4)) * MDP(16) + (-MDP(8) * t288 + t281 * MDP(9)) * t228; (t253 * t230 + t152 + t258) * MDP(16) + t282 * MDP(15) * t228;];
tau = t1;
