% Calculate vector of inverse dynamics joint torques for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:19
% EndTime: 2021-01-15 17:13:22
% DurationCPUTime: 1.22s
% Computational Cost: add. (853->226), mult. (1434->269), div. (0->0), fcn. (777->6), ass. (0->112)
t216 = sin(pkin(7));
t217 = cos(pkin(7));
t221 = sin(qJ(1));
t223 = cos(qJ(1));
t177 = -t221 * t216 - t223 * t217;
t178 = t216 * t223 - t221 * t217;
t291 = -g(1) * t177 - g(2) * t178;
t226 = qJD(1) ^ 2;
t220 = sin(qJ(4));
t213 = t220 ^ 2;
t222 = cos(qJ(4));
t214 = t222 ^ 2;
t244 = MDP(17) * (t213 + t214);
t293 = t226 * t244;
t224 = pkin(1) + pkin(2);
t187 = -qJD(1) * t224 + qJD(2);
t269 = qJD(1) * qJ(2);
t169 = t216 * t187 + t217 * t269;
t164 = -qJD(1) * pkin(6) + t169;
t263 = qJD(4) * t220;
t160 = t164 * t263;
t186 = -qJDD(1) * t224 + qJDD(2);
t260 = qJD(1) * qJD(2);
t261 = qJ(2) * qJDD(1);
t162 = t216 * t186 + (t260 + t261) * t217;
t158 = -qJDD(1) * pkin(6) + t162;
t238 = qJD(3) * qJD(4) + t158;
t230 = qJ(5) * qJDD(1) - t238;
t259 = qJD(1) * qJD(4);
t234 = qJ(5) * t259 + qJDD(3);
t258 = qJD(1) * qJD(5);
t148 = -t160 + t234 * t220 + (-t230 - t258) * t222;
t207 = t222 * qJD(3);
t267 = qJD(1) * t220;
t155 = qJ(5) * t267 - t164 * t220 + t207;
t285 = qJD(4) * pkin(4);
t152 = t155 + t285;
t292 = -qJD(4) * t152 + t148;
t247 = t222 * t259;
t257 = qJDD(1) * t220;
t290 = t220 * t258 + (t247 + t257) * qJ(5);
t219 = qJ(5) + pkin(6);
t288 = pkin(4) * t222;
t249 = pkin(3) + t288;
t289 = t216 * t249 - t217 * t219 + qJ(2);
t210 = g(3) * t222;
t284 = qJDD(1) * pkin(1);
t283 = qJDD(4) * pkin(4);
t282 = t164 * t222;
t281 = t177 * t220;
t280 = t178 * t220;
t279 = t222 * t226;
t274 = t217 * qJ(2) - t216 * t224;
t180 = -pkin(6) + t274;
t278 = qJ(5) - t180;
t277 = qJDD(3) + g(3);
t276 = -t155 + t152;
t275 = g(1) * t280 - g(2) * t281;
t273 = g(1) * t221 - g(2) * t223;
t272 = t213 - t214;
t225 = qJD(4) ^ 2;
t270 = t225 + t226;
t168 = t187 * t217 - t216 * t269;
t163 = qJD(1) * pkin(3) - t168;
t268 = qJD(1) * t163;
t266 = qJD(1) * t222;
t265 = qJD(2) * t217;
t159 = pkin(4) * t266 + qJD(5) + t163;
t262 = t159 * qJD(1);
t256 = qJDD(1) * t222;
t255 = qJDD(4) * t220;
t254 = qJDD(4) * t222;
t253 = MDP(13) + MDP(15);
t252 = MDP(14) + MDP(16);
t251 = 0.2e1 * t260;
t250 = t291 * t222 + t160;
t248 = t220 * t259;
t190 = t216 * t260;
t245 = -t216 * qJ(2) - t217 * t224;
t243 = qJD(4) * t278;
t179 = pkin(3) - t245;
t176 = t179 + t288;
t242 = -qJD(1) * t176 - t159;
t241 = 0.2e1 * t247;
t240 = qJDD(2) - t284;
t239 = -qJD(5) + t265;
t206 = t222 * qJDD(3);
t237 = -g(1) * t281 - g(2) * t280 + t206 + t210;
t236 = -g(1) * t178 + g(2) * t177;
t235 = g(1) * t223 + g(2) * t221;
t161 = t186 * t217 - t216 * t261 - t190;
t232 = -qJD(3) * t220 - t282;
t156 = -qJ(5) * t266 - t232;
t233 = -t152 * t220 + t156 * t222;
t157 = qJDD(1) * pkin(3) - t161;
t229 = pkin(4) * t256 + qJDD(5) + t157;
t149 = -pkin(4) * t248 + t229;
t183 = -pkin(4) * t263 + qJD(2) * t216;
t231 = -qJD(1) * t183 - qJDD(1) * t176 - t149;
t228 = -qJDD(1) * t179 + t180 * t225 - t157 - t190;
t227 = -qJDD(4) * t180 + (-qJD(1) * t179 - t163 - t265) * qJD(4);
t209 = t223 * qJ(2);
t208 = t221 * qJ(2);
t185 = -t220 * t225 + t254;
t184 = -t222 * t225 - t255;
t167 = t216 * t219 + t217 * t249 + t224;
t166 = t278 * t222;
t165 = t278 * t220;
t154 = -t220 * t239 + t222 * t243;
t153 = t220 * t243 + t222 * t239;
t147 = qJD(4) * t232 - t220 * t158 + t206 + t283 + t290;
t1 = [qJDD(1) * MDP(1) + t273 * MDP(2) + t235 * MDP(3) + (-qJDD(2) + t273 + 0.2e1 * t284) * MDP(4) + (-t235 + t251 + 0.2e1 * t261) * MDP(5) + (-t240 * pkin(1) - g(1) * (-pkin(1) * t221 + t209) - g(2) * (pkin(1) * t223 + t208) + (t251 + t261) * qJ(2)) * MDP(6) + (t162 * t274 + t161 * t245 - g(1) * (-t221 * t224 + t209) - g(2) * (t223 * t224 + t208) + (-t168 * t216 + t169 * t217) * qJD(2)) * MDP(7) + (qJDD(1) * t213 + t220 * t241) * MDP(8) + 0.2e1 * (t220 * t256 - t259 * t272) * MDP(9) + t184 * MDP(10) - t185 * MDP(11) + (t227 * t220 + (-t228 + t236) * t222) * MDP(13) + (t220 * t228 + t222 * t227 + t275) * MDP(14) + (qJDD(4) * t165 + (t220 * t242 + t154) * qJD(4) + (-t231 + t236) * t222) * MDP(15) + (qJDD(4) * t166 + t231 * t220 + (t222 * t242 - t153) * qJD(4) + t275) * MDP(16) + ((qJDD(1) * t166 + (qJD(4) * t165 - t153) * qJD(1) - t292) * t222 + (qJD(4) * t156 + qJDD(1) * t165 + t147 + (-qJD(4) * t166 + t154) * qJD(1)) * t220 + t291) * MDP(17) + (-t148 * t166 + t156 * t153 + t147 * t165 + t152 * t154 + t149 * t176 + t159 * t183 - g(1) * (-t167 * t221 + t289 * t223) - g(2) * (t167 * t223 + t289 * t221)) * MDP(18); -qJDD(1) * MDP(4) - t226 * MDP(5) + (-qJ(2) * t226 + t240) * MDP(6) + t252 * ((t241 + t257) * t217 + (t220 * t270 - t254) * t216) + t253 * ((0.2e1 * t248 - t256) * t217 + (-t222 * t270 - t255) * t216) + ((-qJD(1) * t169 + t161) * MDP(7) + (t152 * t267 - t156 * t266 - t149) * MDP(18) + t293) * t217 + ((qJD(1) * t168 + t162) * MDP(7) + (-t147 * t220 - t156 * t263 + t292 * t222 - t262) * MDP(18) - qJDD(1) * t244) * t216 + (-MDP(6) - MDP(7) - MDP(18)) * t273; t277 * MDP(7) + (qJD(4) * t233 + t147 * t222 + t148 * t220 + g(3)) * MDP(18) + t253 * t185 + t252 * t184; -t220 * MDP(8) * t279 + t272 * t226 * MDP(9) - MDP(10) * t257 - MDP(11) * t256 + qJDD(4) * MDP(12) + ((-t158 + t268) * t220 + t237) * MDP(13) + (t207 * qJD(4) + (-qJD(4) * t164 - t277) * t220 + (-t238 + t268) * t222 + t250) * MDP(14) + (0.2e1 * t283 + (t156 - t282) * qJD(4) + (pkin(4) * t279 - t238 + t262) * t220 + t237 + t290) * MDP(15) + (-pkin(4) * t213 * t226 + qJD(4) * t155 + (-g(3) - t234) * t220 + ((qJD(5) + t159) * qJD(1) + t230) * t222 + t250) * MDP(16) + (pkin(4) * t257 + (-t276 + t285) * t266) * MDP(17) + (t276 * t156 + (t210 + t147 + (t262 + t291) * t220) * pkin(4)) * MDP(18); (t229 + t236) * MDP(18) - t293 + (MDP(15) * t222 - MDP(16) * t220) * qJDD(1) + (t233 * MDP(18) + (-0.2e1 * t222 * MDP(16) + (-MDP(18) * pkin(4) - 0.2e1 * MDP(15)) * t220) * qJD(4)) * qJD(1);];
tau = t1;
