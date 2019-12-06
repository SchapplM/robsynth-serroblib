% Calculate vector of inverse dynamics joint torques for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:21
% EndTime: 2019-12-05 15:11:24
% DurationCPUTime: 1.47s
% Computational Cost: add. (729->206), mult. (1639->275), div. (0->0), fcn. (1228->8), ass. (0->112)
t207 = sin(qJ(3));
t209 = cos(qJ(3));
t202 = sin(pkin(8));
t266 = qJD(1) * t202;
t241 = t209 * t266;
t183 = qJD(2) * t207 + t241;
t179 = qJD(3) * pkin(6) + t183;
t206 = sin(qJ(4));
t208 = cos(qJ(4));
t204 = cos(pkin(8));
t265 = qJD(1) * t204;
t222 = t179 * t206 + t208 * t265;
t296 = qJD(5) + t222;
t157 = -qJD(4) * pkin(4) + t296;
t205 = cos(pkin(7));
t272 = t205 * t209;
t273 = t205 * t207;
t203 = sin(pkin(7));
t274 = t203 * t209;
t275 = t203 * t207;
t231 = -g(2) * (-t204 * t275 - t272) - g(1) * (-t204 * t273 + t274);
t278 = t202 * t207;
t245 = g(3) * t278;
t295 = t245 + t231;
t259 = qJD(3) * t208;
t291 = qJD(2) * t209 - t207 * t266;
t294 = qJD(3) * t291;
t224 = pkin(4) * t208 + qJ(5) * t206 + pkin(3);
t293 = t224 * qJD(3);
t292 = t224 * qJDD(3);
t289 = MDP(11) + MDP(13);
t240 = t206 * t265;
t168 = t179 * t208 - t240;
t160 = qJD(4) * qJ(5) + t168;
t254 = qJDD(1) * t202;
t260 = qJD(3) * t207;
t223 = qJD(2) * t260 + qJD(3) * t241 - qJDD(2) * t209 + t207 * t254;
t290 = 0.2e1 * qJDD(3) * pkin(3) - t223;
t280 = qJDD(4) * pkin(4);
t288 = qJDD(5) - t280;
t210 = qJD(4) ^ 2;
t287 = pkin(6) * t210;
t284 = qJD(3) * pkin(3);
t283 = pkin(6) * qJDD(4);
t281 = qJDD(3) * pkin(6);
t279 = t202 * t206;
t277 = t202 * t208;
t276 = t202 * t209;
t271 = qJDD(1) - g(3);
t229 = pkin(4) * t206 - qJ(5) * t208;
t173 = t229 * qJD(4) - qJD(5) * t206;
t270 = t173 - t183;
t200 = t206 ^ 2;
t201 = t208 ^ 2;
t269 = t200 - t201;
t268 = t200 + t201;
t211 = qJD(3) ^ 2;
t267 = t210 + t211;
t169 = -t291 - t293;
t263 = qJD(3) * t169;
t262 = qJD(3) * t173;
t261 = qJD(3) * t206;
t257 = qJD(4) * t206;
t256 = qJD(4) * t208;
t255 = qJD(3) * qJD(4);
t253 = qJDD(1) * t204;
t252 = qJDD(2) * t207;
t251 = qJDD(3) * t208;
t250 = qJDD(3) * t209;
t249 = qJDD(4) * qJ(5);
t248 = qJDD(4) * t206;
t247 = MDP(12) - MDP(15);
t246 = MDP(14) * qJDD(3);
t244 = t208 * t276;
t243 = t206 * t208 * t211;
t242 = t183 * t259 + t208 * t245 + t257 * t291;
t239 = t202 * t260;
t238 = t206 * t255;
t237 = t208 * t255;
t236 = -g(1) * t203 + g(2) * t205;
t178 = -t291 - t284;
t235 = t178 - t284;
t233 = t169 - t293;
t159 = t209 * t254 + t252 + t281 + t294;
t232 = -qJD(4) * t240 + t206 * t159 + t179 * t256 + t208 * t253;
t175 = t204 * t274 - t273;
t177 = t204 * t272 + t275;
t230 = g(1) * t177 + g(2) * t175;
t228 = t208 * t159 - t206 * t253;
t150 = t249 + (qJD(5) - t222) * qJD(4) + t228;
t151 = t232 + t288;
t227 = t150 * t208 + t151 * t206;
t226 = t268 * MDP(14) - MDP(5);
t225 = -qJDD(3) * t207 - t209 * t211;
t180 = t204 * t208 + t206 * t276;
t219 = t231 - t287;
t217 = -qJD(3) * t183 - t295;
t152 = t223 + t262 - t292;
t216 = -t152 + t219 + t292;
t161 = t175 * t206 - t203 * t277;
t163 = t177 * t206 - t205 * t277;
t215 = g(1) * t163 + g(2) * t161 + g(3) * t180 - t232;
t162 = t175 * t208 + t203 * t279;
t164 = t177 * t208 + t205 * t279;
t181 = -t204 * t206 + t244;
t214 = g(1) * t164 + g(2) * t162 + g(3) * t181 - t228;
t213 = qJD(4) * t168 + t215;
t212 = -g(3) * t276 + (t157 * t208 - t160 * t206) * qJD(4) + t227 - t230;
t184 = t229 * qJD(3);
t166 = -qJD(4) * t180 - t208 * t239;
t165 = -qJD(4) * t244 + (qJD(4) * t204 + t239) * t206;
t1 = [t271 * MDP(1) + (qJDD(1) * t204 ^ 2 - g(3)) * MDP(2) + (t150 * t181 + t151 * t180 - t157 * t165 + t160 * t166 - g(3)) * MDP(16) - t247 * (qJD(4) * t166 + qJDD(4) * t181 + (t225 * t206 - t207 * t237) * t202) + (t180 * t206 + t181 * t208) * t246 + (-t165 * t206 + t166 * t208 + (t180 * t208 - t181 * t206) * qJD(4)) * MDP(14) * qJD(3) + ((t207 * t211 - t250) * MDP(5) + (t152 * t207 + t209 * t263) * MDP(16) + MDP(2) * t254 + (t208 * t289 + MDP(4)) * t225) * t202 + t289 * (qJD(4) * t165 - qJDD(4) * t180 + t238 * t278); (qJDD(2) + t236) * MDP(2) + t236 * MDP(16) + t289 * ((-0.2e1 * t238 + t251) * t209 + (-t267 * t208 - t248) * t207) + t247 * ((-qJDD(4) * t207 - 0.2e1 * t209 * t255) * t208 + (t267 * t207 - t250) * t206) + (qJDD(3) * MDP(4) + (t157 * t261 + t160 * t259 - t152) * MDP(16) + t226 * t211) * t209 + (-t211 * MDP(4) + (t157 * t256 - t160 * t257 + t227 + t263) * MDP(16) + t226 * qJDD(3)) * t207; qJDD(3) * MDP(3) + (-t217 - t223) * MDP(4) + (-t271 * t276 + t230 - t252) * MDP(5) + (qJDD(3) * t200 + 0.2e1 * t206 * t237) * MDP(6) + 0.2e1 * (t206 * t251 - t269 * t255) * MDP(7) + (t208 * t210 + t248) * MDP(8) + (qJDD(4) * t208 - t206 * t210) * MDP(9) + ((t235 * qJD(4) - t283) * t206 + (t219 + t290) * t208 + t242) * MDP(11) + ((-t283 + (t291 + t235) * qJD(4)) * t208 + (t217 + t287 - t290) * t206) * MDP(12) + ((t233 * qJD(4) - t283) * t206 + (t216 - t262) * t208 + t242) * MDP(13) + (t212 + (t281 - t294) * t268) * MDP(14) + ((t283 + (-t291 - t233) * qJD(4)) * t208 + (-t270 * qJD(3) + t216 + t245) * t206) * MDP(15) + ((-t157 * t206 - t160 * t208) * t291 + t270 * t169 + t212 * pkin(6) + (-t152 + t295) * t224) * MDP(16); -MDP(6) * t243 + t269 * t211 * MDP(7) + MDP(9) * t251 + qJDD(4) * MDP(10) + (-t178 * t261 + t213) * MDP(11) + (-t178 * t259 + t214) * MDP(12) + (0.2e1 * t280 - qJDD(5) + (-t169 * t206 + t184 * t208) * qJD(3) + t213) * MDP(13) + (0.2e1 * t249 + (t169 * t208 + t184 * t206) * qJD(3) + 0.2e1 * qJD(4) * qJD(5) - t214) * MDP(15) + (t150 * qJ(5) - t151 * pkin(4) - t169 * t184 - t157 * t168 - g(1) * (-pkin(4) * t163 + qJ(5) * t164) - g(2) * (-pkin(4) * t161 + qJ(5) * t162) - g(3) * (-pkin(4) * t180 + qJ(5) * t181) + t296 * t160) * MDP(16) + (-t229 * MDP(14) + t206 * MDP(8)) * qJDD(3); (-qJDD(4) - t243) * MDP(13) + t206 * t246 + (-t200 * t211 - t210) * MDP(15) + (-qJD(4) * t160 + t169 * t261 - t215 + t288) * MDP(16);];
tau = t1;
