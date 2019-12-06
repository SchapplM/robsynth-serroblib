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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:22:29
% EndTime: 2019-12-05 18:22:32
% DurationCPUTime: 1.26s
% Computational Cost: add. (1101->207), mult. (1761->270), div. (0->0), fcn. (978->12), ass. (0->116)
t232 = -qJ(5) - pkin(7);
t237 = cos(qJ(2));
t288 = pkin(1) * qJD(2);
t266 = qJD(1) * t288;
t234 = sin(qJ(2));
t270 = qJDD(1) * t234;
t297 = pkin(1) * t270 + t237 * t266;
t229 = qJ(1) + qJ(2);
t218 = pkin(8) + t229;
t206 = sin(t218);
t207 = cos(t218);
t296 = g(2) * t207 + g(3) * t206;
t221 = sin(t229);
t222 = cos(t229);
t295 = g(2) * t222 + g(3) * t221;
t251 = -g(2) * t206 + g(3) * t207;
t294 = pkin(1) * t237;
t231 = cos(pkin(8));
t293 = pkin(2) * t231;
t236 = cos(qJ(4));
t292 = g(1) * t236;
t287 = qJD(4) * pkin(4);
t225 = qJDD(1) + qJDD(2);
t233 = sin(qJ(4));
t286 = t225 * t233;
t230 = sin(pkin(8));
t285 = t230 * t234;
t284 = t231 * t234;
t283 = t233 * t236;
t214 = pkin(2) + t294;
t277 = pkin(1) * t284 + t230 * t214;
t179 = pkin(7) + t277;
t282 = -qJ(5) - t179;
t208 = pkin(2) * t230 + pkin(7);
t281 = -qJ(5) - t208;
t280 = qJDD(3) - g(1);
t215 = qJDD(1) * t294;
t177 = pkin(2) * t225 - t234 * t266 + t215;
t161 = t231 * t177 - t297 * t230;
t157 = -t225 * pkin(3) - t161;
t226 = qJD(1) + qJD(2);
t272 = qJD(1) * t237;
t267 = pkin(1) * t272;
t189 = pkin(2) * t226 + t267;
t273 = qJD(1) * t234;
t268 = pkin(1) * t273;
t201 = t230 * t268;
t171 = t189 * t231 - t201;
t167 = -pkin(3) * t226 - t171;
t279 = t167 * qJD(4) * t236 + t157 * t233;
t202 = t231 * t268;
t172 = t230 * t189 + t202;
t260 = -t232 * t226 + t172;
t159 = t236 * qJD(3) - t260 * t233;
t156 = t159 + t287;
t278 = t156 - t159;
t227 = t233 ^ 2;
t228 = t236 ^ 2;
t276 = -t227 - t228;
t275 = t227 - t228;
t213 = pkin(4) * t236 + pkin(3);
t163 = -t213 * t226 + qJD(5) - t171;
t274 = MDP(16) * t163;
t271 = qJD(4) * t233;
t269 = qJDD(4) * t179;
t264 = t167 * t271 + t296 * t236;
t162 = t230 * t177 + t297 * t231;
t263 = t215 + t295;
t262 = t226 * t271;
t261 = -g(2) * t221 + g(3) * t222;
t258 = -pkin(1) * t285 + t214 * t231;
t178 = -pkin(3) - t258;
t183 = (t231 * t237 - t285) * t288;
t259 = t178 * t226 - t183;
t257 = qJD(4) * t282;
t256 = qJD(4) * t281;
t255 = qJD(1) * (-qJD(2) + t226);
t158 = pkin(7) * t225 + t162;
t241 = qJ(5) * t225 + qJD(3) * qJD(4) + qJD(5) * t226 + t158;
t246 = t260 * qJD(4);
t151 = (qJDD(3) - t246) * t233 + t241 * t236;
t253 = t151 * t236 - t251;
t235 = sin(qJ(1));
t238 = cos(qJ(1));
t250 = g(2) * t238 + g(3) * t235;
t239 = qJD(4) ^ 2;
t190 = qJDD(4) * t233 + t236 * t239;
t191 = qJDD(4) * t236 - t233 * t239;
t249 = 0.2e1 * (-t275 * t226 * qJD(4) + t225 * t283) * MDP(9) + (t225 * t227 + 0.2e1 * t236 * t262) * MDP(8) + t190 * MDP(10) + t191 * MDP(11) + t225 * MDP(4);
t160 = qJD(3) * t233 + t260 * t236;
t248 = t156 * t233 - t160 * t236;
t247 = t230 * t237 + t284;
t181 = t247 * t288;
t245 = -t178 * t225 - t179 * t239 - t181 * t226;
t180 = t230 * t267 + t202;
t209 = -pkin(3) - t293;
t244 = t180 * t226 - t208 * t239 - t209 * t225;
t243 = -t167 * t226 - t158 + t251;
t182 = t231 * t267 - t201;
t242 = -qJDD(4) * t208 + (t209 * t226 + t182) * qJD(4);
t152 = pkin(4) * t262 - t213 * t225 + qJDD(5) - t161;
t240 = -g(2) * (-pkin(2) * t222 + t206 * t232 - t207 * t213) - g(3) * (-pkin(2) * t221 - t206 * t213 - t207 * t232);
t224 = t226 ^ 2;
t223 = t236 * qJ(5);
t219 = t236 * qJD(5);
t217 = t236 * qJDD(3);
t186 = t208 * t236 + t223;
t185 = t281 * t233;
t175 = -qJD(5) * t233 + t236 * t256;
t174 = t233 * t256 + t219;
t170 = t179 * t236 + t223;
t169 = t282 * t233;
t154 = (-qJD(5) - t183) * t233 + t236 * t257;
t153 = t183 * t236 + t233 * t257 + t219;
t150 = qJDD(4) * pkin(4) - t241 * t233 - t236 * t246 + t217;
t1 = [qJDD(1) * MDP(1) + t250 * MDP(2) + (-g(2) * t235 + g(3) * t238) * MDP(3) + t263 * MDP(5) + t261 * MDP(6) + (pkin(2) * t295 + t161 * t258 + t162 * t277 - t171 * t181 + t172 * t183) * MDP(7) + t264 * MDP(13) + t279 * MDP(14) + t253 * MDP(15) + (t150 * t169 + t151 * t170 + t152 * t178 + t160 * t153 + t156 * t154 + t240) * MDP(16) + (t225 * t237 * MDP(5) + (-qJDD(1) - t225) * MDP(6) * t234 + ((-t226 * t234 - t273) * MDP(5) + (-t226 * t237 - t272) * MDP(6) + t247 * t274) * qJD(2) + (MDP(7) + MDP(16)) * t250) * pkin(1) + ((-t157 + t245) * MDP(13) - MDP(14) * t269 + (t153 * t226 + t170 * t225) * MDP(15) - pkin(4) * t152 * MDP(16) + (t259 * MDP(14) + (-t169 * t226 - t156) * MDP(15)) * qJD(4)) * t236 + (-MDP(13) * t269 + (-t245 - t296) * MDP(14) + (-t154 * t226 - t169 * t225 - t150) * MDP(15) + (t259 * MDP(13) + (-t170 * t226 - t160) * MDP(15) + pkin(4) * t274) * qJD(4)) * t233 + t249; (t234 * pkin(1) * t255 + t263) * MDP(5) + ((t237 * t255 - t270) * pkin(1) + t261) * MDP(6) + (t171 * t180 - t172 * t182 + (t161 * t231 + t162 * t230 + t295) * pkin(2)) * MDP(7) + (t242 * t233 + (-t157 + t244) * t236 + t264) * MDP(13) + (t242 * t236 + (-t244 - t296) * t233 + t279) * MDP(14) + ((-qJD(4) * t156 + t186 * t225) * t236 + (-qJD(4) * t160 - t185 * t225 - t150) * t233 + (t174 * t236 - t175 * t233 + t276 * t182 + (-t185 * t236 - t186 * t233) * qJD(4)) * t226 + t253) * MDP(15) + (t151 * t186 + t150 * t185 + t152 * (-t213 - t293) + (pkin(4) * t271 - t180) * t163 + (-t182 * t236 + t174) * t160 + (t182 * t233 + t175) * t156 + t240) * MDP(16) + t249; t280 * MDP(7) + t191 * MDP(13) - t190 * MDP(14) + (-t248 * qJD(4) + t150 * t236 + t151 * t233 - g(1)) * MDP(16); MDP(10) * t286 + t236 * t225 * MDP(11) + qJDD(4) * MDP(12) + (t243 * t233 + t217 - t292) * MDP(13) + (-t280 * t233 + t243 * t236) * MDP(14) + (-pkin(4) * t286 + (t278 - t287) * t236 * t226) * MDP(15) + (t278 * t160 + (-t292 + t150 + (-t163 * t226 + t251) * t233) * pkin(4)) * MDP(16) + (-MDP(8) * t283 + t275 * MDP(9)) * t224; (t248 * t226 + t152 - t296) * MDP(16) + t276 * MDP(15) * t224;];
tau = t1;
