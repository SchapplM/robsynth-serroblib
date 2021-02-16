% Calculate vector of inverse dynamics joint torques for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:56
% EndTime: 2021-01-15 15:05:01
% DurationCPUTime: 2.12s
% Computational Cost: add. (940->250), mult. (1910->326), div. (0->0), fcn. (1234->6), ass. (0->122)
t221 = sin(pkin(8));
t225 = cos(qJ(4));
t263 = qJD(2) * qJD(4);
t249 = t225 * t263;
t224 = sin(qJ(4));
t260 = qJDD(2) * t224;
t308 = t249 + t260;
t309 = t221 * t308;
t218 = pkin(7) + qJ(2);
t214 = sin(t218);
t215 = cos(t218);
t307 = g(1) * t214 - g(2) * t215;
t264 = qJD(2) * qJD(3);
t266 = qJ(3) * qJDD(2);
t232 = t264 + t266;
t222 = cos(pkin(8));
t306 = -t221 * (-qJ(5) - pkin(6)) + (pkin(4) * t225 + pkin(3)) * t222;
t286 = t222 * t224;
t176 = t214 * t286 + t215 * t225;
t178 = t214 * t225 - t215 * t286;
t305 = -g(1) * t178 + g(2) * t176;
t276 = qJD(2) * t222;
t203 = -qJD(4) + t276;
t304 = qJD(4) + t203;
t289 = t221 * t224;
t303 = g(3) * t289 + t305;
t302 = pkin(4) * t309 + qJDD(5);
t261 = qJDD(2) * t222;
t201 = -qJDD(4) + t261;
t301 = pkin(4) * t201;
t299 = pkin(4) * t224;
t295 = qJDD(2) * pkin(2);
t212 = t222 * qJDD(1);
t180 = t232 * t221 - t212;
t294 = t180 * t221;
t292 = t215 * t224;
t216 = t221 ^ 2;
t226 = qJD(2) ^ 2;
t291 = t216 * t226;
t288 = t221 * t225;
t287 = t222 * MDP(5);
t285 = t222 * t225;
t194 = -pkin(3) * t222 - pkin(6) * t221 - pkin(2);
t183 = t194 * qJD(2) + qJD(3);
t192 = qJ(3) * t276 + qJD(1) * t221;
t247 = t225 * t183 - t192 * t224;
t274 = qJD(2) * t225;
t252 = t221 * t274;
t165 = -qJ(5) * t252 + t247;
t158 = -pkin(4) * t203 + t165;
t284 = -t165 + t158;
t270 = qJD(4) * t225;
t273 = qJD(3) * t222;
t283 = t194 * t270 + t225 * t273;
t202 = qJ(3) * t285;
t282 = t224 * t194 + t202;
t281 = t215 * pkin(2) + t214 * qJ(3);
t280 = t222 ^ 2 + t216;
t219 = t224 ^ 2;
t220 = t225 ^ 2;
t279 = t219 - t220;
t278 = MDP(17) * t221;
t190 = (qJ(3) + t299) * t221;
t277 = qJD(2) * t190;
t275 = qJD(2) * t224;
t272 = qJD(4) * t192;
t271 = qJD(4) * t224;
t269 = qJD(5) * t221;
t268 = t201 * MDP(12);
t213 = t222 * qJD(1);
t171 = qJD(5) - t213 + t277;
t267 = qJD(5) + t171;
t265 = qJ(5) * qJDD(2);
t262 = qJD(2) * qJD(5);
t259 = qJDD(2) * t225;
t258 = MDP(13) + MDP(15);
t257 = MDP(14) + MDP(16);
t256 = qJ(3) * t286;
t255 = qJ(5) * t288;
t254 = t224 * t291;
t253 = t221 * t275;
t251 = t203 * t271;
t250 = t224 * t263;
t246 = MDP(17) * (-t219 - t220);
t245 = t171 + t277;
t181 = qJDD(1) * t221 + t232 * t222;
t244 = -qJD(4) * t183 - t181;
t243 = t201 - t261;
t242 = t201 + t261;
t182 = t194 * qJDD(2) + qJDD(3);
t241 = t225 * t181 + t224 * t182 + t183 * t270 - t192 * t271;
t240 = t221 * t250;
t239 = -g(1) * t176 - g(2) * t178;
t177 = -t214 * t285 + t292;
t179 = t214 * t224 + t215 * t285;
t238 = -g(1) * t177 - g(2) * t179;
t237 = g(1) * t215 + g(2) * t214;
t235 = t181 * t222 + t294;
t234 = -t183 * t224 - t192 * t225;
t191 = qJ(3) * qJD(2) * t221 - t213;
t233 = t191 * t221 + t192 * t222;
t167 = t180 + t302;
t187 = (pkin(4) * t270 + qJD(3)) * t221;
t231 = qJD(2) * t187 + qJDD(2) * t190 + t167;
t174 = t225 * t182;
t229 = qJ(5) * t240 + t244 * t224 + t174;
t228 = -t203 ^ 2 - t291;
t227 = g(1) * t179 - g(2) * t177 + g(3) * t288 - t241;
t211 = qJDD(3) - t295;
t207 = t215 * qJ(3);
t199 = t221 * t259;
t193 = t222 * t240;
t189 = t225 * t194;
t186 = t203 * t252;
t169 = -qJ(5) * t289 + t282;
t168 = -t255 + t189 + (-qJ(3) * t224 - pkin(4)) * t222;
t166 = -qJ(5) * t253 - t234;
t162 = -t224 * t273 - t225 * t269 + (-t202 + (qJ(5) * t221 - t194) * t224) * qJD(4);
t161 = -t224 * t269 + (-t255 - t256) * qJD(4) + t283;
t157 = (-qJ(5) * t308 - t224 * t262) * t221 + t241;
t156 = -t301 + (-t272 + (-t262 - t265) * t221) * t225 + t229;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-t180 * t222 - g(3)) * MDP(7) + (-t167 * t222 - g(3)) * MDP(18) + t257 * t193 + (t258 * (t243 * t224 + (t203 - t276) * t270) + t257 * (t243 * t225 - t251) + t181 * MDP(7) + (-t156 * t224 + t157 * t225 - t158 * t270 - t166 * t271) * MDP(18)) * t221; qJDD(2) * MDP(2) + t307 * MDP(3) + t237 * MDP(4) + (-t211 + t307 + t295) * t287 + (t232 * t280 + t235 - t237) * MDP(6) + (-t211 * pkin(2) - g(1) * (-pkin(2) * t214 + t207) - g(2) * t281 + t233 * qJD(3) + t235 * qJ(3)) * MDP(7) + (qJDD(2) * t220 - 0.2e1 * t224 * t249) * t216 * MDP(8) + 0.2e1 * (-t224 * t259 + t279 * t263) * t216 * MDP(9) + (t193 + (-t242 * t225 + t251) * t221) * MDP(10) + (t242 * t224 + (t203 + t276) * t270) * t221 * MDP(11) + t222 * t268 + (-t174 * t222 - t189 * t201 + ((qJD(2) * t216 + t203 * t222) * qJ(3) + t233) * t270 + (-(-qJD(4) * t194 - t273) * t203 - t244 * t222 + t216 * t264 + t294 + (qJDD(2) * t216 + t201 * t222) * qJ(3)) * t224 + t238) * MDP(13) + ((-qJD(4) * t256 + t283) * t203 + t282 * t201 + t241 * t222 + (t180 * t225 - t191 * t271) * t221 + (t225 * t264 + (-t250 + t259) * qJ(3)) * t216 + t239) * MDP(14) + (-t156 * t222 - t162 * t203 - t168 * t201 + (t231 * t224 + t245 * t270) * t221 + t238) * MDP(15) + (t157 * t222 + t161 * t203 + t169 * t201 + (t231 * t225 - t245 * t271) * t221 + t239) * MDP(16) + ((-qJD(4) * t166 - qJDD(2) * t168 - t156 + (-qJD(4) * t169 - t162) * qJD(2)) * t225 + (qJD(4) * t158 - qJDD(2) * t169 - t157 + (qJD(4) * t168 - t161) * qJD(2)) * t224 + t307) * t278 + (t157 * t169 + t166 * t161 + t156 * t168 + t158 * t162 + t167 * t190 + t171 * t187 - g(1) * (pkin(4) * t292 + t207) - g(2) * (t306 * t215 + t281) + (-g(1) * (-pkin(2) - t306) - g(2) * t299) * t214) * MDP(18); (qJDD(3) - t307) * MDP(7) + (t156 * t225 + t157 * t224 - t158 * t271 + t166 * t270 - t307) * MDP(18) - t280 * MDP(6) * t226 + t257 * (t224 * t201 + t228 * t225) + t258 * (-t201 * t225 + t228 * t224) + (-pkin(2) * MDP(7) + t221 * t246 - t287) * qJDD(2) + (-t233 * MDP(7) + (t158 * t286 - t166 * t285 - t171 * t221) * MDP(18)) * qJD(2); t225 * MDP(8) * t254 - t279 * MDP(9) * t291 + (-t304 * t253 + t199) * MDP(10) + (-t186 - t309) * MDP(11) - t268 + (-t224 * t181 - t191 * t252 + t304 * t234 + t174 + t303) * MDP(13) + (t191 * t253 - t247 * t203 + t227) * MDP(14) + (-0.2e1 * t301 - t166 * t203 + (-pkin(4) * t254 - t272 + (-t267 * qJD(2) - t265) * t221) * t225 + t229 + t303) * MDP(15) + (-pkin(4) * t220 * t291 - t165 * t203 + (qJ(5) * t260 + (qJ(5) * t270 + t267 * t224) * qJD(2)) * t221 + t227) * MDP(16) + (-pkin(4) * t259 + (pkin(4) * qJD(4) - t284) * t275) * t278 + (t284 * t166 + (t156 + (g(3) * t224 - t171 * t274) * t221 + t305) * pkin(4)) * MDP(18); -t186 * MDP(15) + t199 * MDP(16) + (g(3) * t222 - t212 + t302) * MDP(18) + t246 * t291 + (MDP(15) * t260 + (-t237 + t266) * MDP(18) + (MDP(15) * t270 + (t158 * t225 + qJD(3)) * MDP(18) + ((-qJD(4) + t203) * MDP(16) + t166 * MDP(18)) * t224) * qJD(2)) * t221;];
tau = t1;
