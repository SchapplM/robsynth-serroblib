% Calculate vector of inverse dynamics joint torques for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:58
% EndTime: 2019-12-31 17:05:00
% DurationCPUTime: 1.35s
% Computational Cost: add. (1037->208), mult. (2476->289), div. (0->0), fcn. (1785->10), ass. (0->107)
t248 = sin(pkin(7));
t249 = cos(pkin(7));
t252 = sin(qJ(2));
t255 = cos(qJ(2));
t221 = -t248 * t252 + t249 * t255;
t211 = t221 * qJD(1);
t254 = cos(qJ(4));
t203 = t254 * t211;
t222 = t248 * t255 + t249 * t252;
t213 = t222 * qJD(1);
t251 = sin(qJ(4));
t287 = t213 * t251;
t177 = t203 - t287;
t245 = qJD(2) + qJD(4);
t288 = t177 * t245;
t253 = sin(qJ(1));
t256 = cos(qJ(1));
t273 = g(1) * t256 + g(2) * t253;
t291 = qJ(3) + pkin(5);
t275 = qJD(2) * t291;
t209 = -qJD(3) * t252 - t255 * t275;
t233 = t291 * t252;
t183 = qJDD(2) * pkin(2) + qJD(1) * t209 - qJDD(1) * t233;
t208 = qJD(3) * t255 - t252 * t275;
t234 = t291 * t255;
t190 = qJD(1) * t208 + qJDD(1) * t234;
t156 = t249 * t183 - t190 * t248;
t281 = qJD(1) * qJD(2);
t276 = t255 * t281;
t277 = t252 * t281;
t187 = qJDD(1) * t222 - t248 * t277 + t249 * t276;
t154 = qJDD(2) * pkin(3) - pkin(6) * t187 + t156;
t157 = t248 * t183 + t249 * t190;
t212 = t222 * qJD(2);
t186 = -qJD(1) * t212 + qJDD(1) * t221;
t155 = pkin(6) * t186 + t157;
t226 = qJD(1) * t233;
t290 = qJD(2) * pkin(2);
t220 = -t226 + t290;
t227 = qJD(1) * t234;
t286 = t249 * t227;
t185 = t248 * t220 + t286;
t296 = pkin(6) * t211;
t165 = t185 + t296;
t242 = pkin(2) * t255 + pkin(1);
t228 = -qJD(1) * t242 + qJD(3);
t193 = -pkin(3) * t211 + t228;
t243 = qJ(2) + pkin(7) + qJ(4);
t238 = sin(t243);
t239 = cos(t243);
t282 = qJD(4) * t251;
t301 = g(3) * t238 - t251 * t154 - t254 * t155 + t165 * t282 - t193 * t177 + t273 * t239;
t244 = qJDD(2) + qJDD(4);
t268 = t211 * t251 + t254 * t213;
t300 = t244 * MDP(17) + (-t177 ^ 2 + t268 ^ 2) * MDP(14) - t177 * MDP(13) * t268;
t289 = t268 * t245;
t298 = -g(3) * t239 + t254 * t154 - t251 * t155 - t193 * t268 + t273 * t238;
t274 = -t254 * t186 + t187 * t251;
t153 = qJD(4) * t268 + t274;
t297 = pkin(2) * t248;
t295 = pkin(6) * t213;
t292 = g(3) * t255;
t216 = t248 * t227;
t184 = t249 * t220 - t216;
t163 = qJD(2) * pkin(3) + t184 - t295;
t285 = t254 * t163;
t169 = t249 * t208 + t248 * t209;
t192 = -t249 * t226 - t216;
t195 = -t248 * t233 + t249 * t234;
t246 = t252 ^ 2;
t284 = -t255 ^ 2 + t246;
t280 = qJDD(1) * t255;
t279 = t252 * t290;
t278 = qJD(4) * t203 + t251 * t186 + t254 * t187;
t168 = -t208 * t248 + t249 * t209;
t191 = t226 * t248 - t286;
t194 = -t249 * t233 - t234 * t248;
t272 = g(1) * t253 - g(2) * t256;
t271 = -t251 * t163 - t254 * t165;
t170 = -pkin(6) * t222 + t194;
t171 = pkin(6) * t221 + t195;
t270 = t170 * t254 - t171 * t251;
t269 = t170 * t251 + t171 * t254;
t267 = t254 * t221 - t222 * t251;
t189 = t221 * t251 + t222 * t254;
t240 = pkin(2) * t249 + pkin(3);
t266 = t240 * t251 + t254 * t297;
t265 = t240 * t254 - t251 * t297;
t264 = -0.2e1 * pkin(1) * t281 - pkin(5) * qJDD(2);
t152 = -t213 * t282 + t278;
t263 = pkin(2) * t277 - qJDD(1) * t242 + qJDD(3);
t257 = qJD(2) ^ 2;
t260 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t257 + t272;
t258 = qJD(1) ^ 2;
t259 = pkin(1) * t258 - pkin(5) * qJDD(1) + t273;
t215 = t221 * qJD(2);
t199 = -pkin(3) * t221 - t242;
t197 = pkin(3) * t212 + t279;
t196 = pkin(2) * qJD(1) * t252 + pkin(3) * t213;
t167 = t192 - t295;
t166 = t191 - t296;
t164 = -pkin(3) * t186 + t263;
t161 = -pkin(6) * t212 + t169;
t160 = -pkin(6) * t215 + t168;
t159 = qJD(4) * t189 + t254 * t212 + t215 * t251;
t158 = qJD(4) * t267 - t212 * t251 + t215 * t254;
t1 = [qJDD(1) * MDP(1) + t272 * MDP(2) + t273 * MDP(3) + (qJDD(1) * t246 + 0.2e1 * t252 * t276) * MDP(4) + 0.2e1 * (t252 * t280 - t281 * t284) * MDP(5) + (qJDD(2) * t252 + t255 * t257) * MDP(6) + (qJDD(2) * t255 - t252 * t257) * MDP(7) + (t252 * t264 + t255 * t260) * MDP(9) + (-t252 * t260 + t255 * t264) * MDP(10) + (-t156 * t222 + t157 * t221 - t168 * t213 + t169 * t211 - t184 * t215 - t185 * t212 + t186 * t195 - t187 * t194 - t273) * MDP(11) + (t157 * t195 + t185 * t169 + t156 * t194 + t184 * t168 - t263 * t242 + t228 * t279 - g(1) * (-t242 * t253 + t256 * t291) - g(2) * (t242 * t256 + t253 * t291)) * MDP(12) + (t152 * t189 + t158 * t268) * MDP(13) + (t152 * t267 - t153 * t189 + t158 * t177 - t159 * t268) * MDP(14) + (t158 * t245 + t189 * t244) * MDP(15) + (-t159 * t245 + t244 * t267) * MDP(16) + (-t197 * t177 + t199 * t153 - t164 * t267 + t193 * t159 + (-qJD(4) * t269 + t160 * t254 - t161 * t251) * t245 + t270 * t244 + t272 * t239) * MDP(18) + (t197 * t268 + t199 * t152 + t164 * t189 + t193 * t158 - (qJD(4) * t270 + t160 * t251 + t161 * t254) * t245 - t269 * t244 - t272 * t238) * MDP(19); t252 * qJDD(1) * MDP(6) + MDP(7) * t280 + qJDD(2) * MDP(8) + (t252 * t259 - t292) * MDP(9) + (g(3) * t252 + t255 * t259) * MDP(10) + ((t185 + t191) * t213 + (t184 - t192) * t211 + (t186 * t248 - t187 * t249) * pkin(2)) * MDP(11) + (-t184 * t191 - t185 * t192 + (-t292 + t156 * t249 + t157 * t248 + (-qJD(1) * t228 + t273) * t252) * pkin(2)) * MDP(12) + (t152 - t288) * MDP(15) + (-t153 + t289) * MDP(16) + (t265 * t244 + t196 * t177 - (t166 * t254 - t167 * t251) * t245 + (-t245 * t266 + t271) * qJD(4) + t298) * MDP(18) + (-t266 * t244 - t196 * t268 + (t166 * t251 + t167 * t254) * t245 + (-t245 * t265 - t285) * qJD(4) + t301) * MDP(19) + (-MDP(4) * t252 * t255 + MDP(5) * t284) * t258 + t300; (-t211 ^ 2 - t213 ^ 2) * MDP(11) + (t184 * t213 - t185 * t211 + t263 - t272) * MDP(12) + (t153 + t289) * MDP(18) + (t152 + t288) * MDP(19); (t278 - t288) * MDP(15) + (-t274 + t289) * MDP(16) + (-t245 * t271 + t298) * MDP(18) + ((-t165 * t251 + t285) * t245 + t301) * MDP(19) + (-MDP(15) * t287 - t268 * MDP(16) + t271 * MDP(18) - MDP(19) * t285) * qJD(4) + t300;];
tau = t1;
