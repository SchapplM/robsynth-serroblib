% Calculate vector of inverse dynamics joint torques for
% S4RRPR8
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:27
% EndTime: 2019-12-31 17:08:29
% DurationCPUTime: 1.70s
% Computational Cost: add. (609->231), mult. (1323->302), div. (0->0), fcn. (802->6), ass. (0->124)
t213 = cos(qJ(2));
t259 = qJD(1) * qJD(2);
t251 = t213 * t259;
t210 = sin(qJ(2));
t257 = qJDD(1) * t210;
t288 = t251 + t257;
t214 = cos(qJ(1));
t278 = g(2) * t214;
t211 = sin(qJ(1));
t281 = g(1) * t211;
t287 = -t278 + t281;
t279 = g(2) * t211;
t244 = g(1) * t214 + t279;
t203 = qJD(2) - qJD(4);
t286 = qJD(4) + t203;
t209 = sin(qJ(4));
t212 = cos(qJ(4));
t171 = t209 * t210 + t212 * t213;
t165 = t171 * qJD(1);
t263 = qJD(1) * t213;
t264 = qJD(1) * t210;
t167 = -t209 * t263 + t212 * t264;
t202 = qJDD(2) - qJDD(4);
t285 = t167 * t165 * MDP(15) - (t165 ^ 2 - t167 ^ 2) * MDP(16) - t202 * MDP(19);
t201 = t213 * pkin(2);
t200 = t210 * qJ(3);
t248 = pkin(1) + t200;
t235 = t248 + t201;
t169 = t235 * qJD(1);
t267 = t201 + t200;
t179 = -pkin(1) - t267;
t277 = pkin(5) * qJDD(2);
t284 = qJD(2) * (qJD(1) * t179 - t169) - t277;
t283 = pkin(2) + pkin(3);
t282 = pkin(5) - pkin(6);
t276 = qJ(3) * t213;
t275 = qJDD(2) * pkin(2);
t274 = t165 * t203;
t273 = t167 * t203;
t272 = t209 * t213;
t271 = t210 * t212;
t270 = t210 * t214;
t217 = qJD(1) ^ 2;
t269 = t210 * t217;
t258 = qJD(1) * qJD(4);
t250 = t213 * t258;
t252 = t210 * t259;
t268 = -t209 * t250 - t212 * t252;
t207 = t210 ^ 2;
t208 = t213 ^ 2;
t266 = t207 - t208;
t243 = pkin(2) * t210 - t276;
t261 = qJD(3) * t210;
t164 = qJD(2) * t243 - t261;
t265 = qJD(1) * t164;
t206 = qJD(2) * qJ(3);
t262 = qJD(2) * t210;
t198 = pkin(5) * t264;
t260 = -pkin(6) * t264 + qJD(3) + t198;
t256 = qJDD(1) * t213;
t182 = t282 * t213;
t255 = t213 * t269;
t254 = t209 * t252 + t288 * t212;
t192 = pkin(5) * t251;
t196 = pkin(5) * t257;
t249 = qJDD(3) + t192 + t196;
t247 = -qJD(2) * pkin(2) + qJD(3);
t199 = pkin(5) * t263;
t176 = -pkin(6) * t263 + t199;
t216 = qJD(2) ^ 2;
t245 = pkin(5) * t216 + t278;
t242 = qJ(3) * t212 - t209 * t283;
t241 = -qJ(3) * t209 - t212 * t283;
t178 = t198 + t247;
t180 = t199 + t206;
t240 = t178 * t213 - t180 * t210;
t181 = t282 * t210;
t239 = t181 * t212 - t182 * t209;
t238 = t181 * t209 + t182 * t212;
t237 = -t271 + t272;
t236 = g(1) * t270 - g(3) * t213 + t210 * t279 - t196;
t234 = -t210 * t283 + t276;
t233 = -0.2e1 * pkin(1) * t259 - t277;
t150 = -pkin(6) * t288 - qJDD(2) * t283 + t249;
t168 = t176 + t206;
t232 = t286 * t168 - t150;
t197 = pkin(5) * t256;
t204 = qJDD(2) * qJ(3);
t205 = qJD(2) * qJD(3);
t156 = -pkin(5) * t252 + t197 + t204 + t205;
t151 = (t252 - t256) * pkin(6) + t156;
t157 = -qJD(2) * t283 + t260;
t231 = -t286 * t157 - t151;
t229 = -qJDD(3) + t236;
t227 = -t210 * t258 - t256;
t226 = t171 * qJD(4);
t225 = t213 * t283 + t248;
t224 = -t209 * t256 + t254;
t223 = 0.2e1 * qJDD(1) * pkin(1) - t245;
t222 = qJD(2) * t272 + qJD(4) * t271;
t155 = t225 * qJD(1);
t160 = t237 * t211;
t162 = -t212 * t270 + t214 * t272;
t221 = g(1) * t162 + g(2) * t160 + g(3) * t171 - t155 * t167;
t161 = t171 * t211;
t163 = t171 * t214;
t220 = g(1) * t163 + g(2) * t161 - g(3) * t237 + t155 * t165;
t154 = qJD(2) * t234 + t261;
t149 = -qJDD(1) * t235 + t265;
t219 = -qJDD(1) * t179 - t149 - t245 - t265;
t158 = t249 - t275;
t218 = qJD(2) * t240 + t156 * t213 + t158 * t210 - t244;
t147 = qJD(1) * t222 + qJDD(1) * t171 + t268;
t194 = t213 * t281;
t177 = qJD(2) * t182;
t175 = t282 * t262;
t173 = t243 * qJD(1);
t170 = pkin(3) * t213 - t179;
t159 = t234 * qJD(1);
t153 = qJD(2) * t171 - t226;
t152 = -qJD(4) * t272 - t212 * t262 + t222;
t148 = qJD(1) * t154 + qJDD(1) * t225;
t146 = -qJD(1) * t226 + t224;
t1 = [qJDD(1) * MDP(1) + t287 * MDP(2) + t244 * MDP(3) + (qJDD(1) * t207 + 0.2e1 * t210 * t251) * MDP(4) + 0.2e1 * (t210 * t256 - t259 * t266) * MDP(5) + (qJDD(2) * t210 + t213 * t216) * MDP(6) + (qJDD(2) * t213 - t210 * t216) * MDP(7) + (t210 * t233 + t213 * t223 + t194) * MDP(9) + (t233 * t213 + (-t223 - t281) * t210) * MDP(10) + (t284 * t210 + t219 * t213 + t194) * MDP(11) + ((t207 + t208) * qJDD(1) * pkin(5) + t218) * MDP(12) + (-t284 * t213 + (t219 + t281) * t210) * MDP(13) + (pkin(5) * t218 + t149 * t179 - t169 * t164 + t287 * t235) * MDP(14) + (-t146 * t237 + t153 * t167) * MDP(15) + (-t146 * t171 + t147 * t237 - t152 * t167 - t153 * t165) * MDP(16) + (-t153 * t203 + t202 * t237) * MDP(17) + (t152 * t203 + t171 * t202) * MDP(18) + (t154 * t165 + t170 * t147 + t148 * t171 + t155 * t152 - (-qJD(4) * t238 + t175 * t209 + t177 * t212) * t203 - t239 * t202 + g(1) * t161 - g(2) * t163) * MDP(20) + (t154 * t167 + t170 * t146 - t148 * t237 + t155 * t153 + (qJD(4) * t239 - t175 * t212 + t177 * t209) * t203 + t238 * t202 - g(1) * t160 + g(2) * t162) * MDP(21); -MDP(4) * t255 + t266 * MDP(5) * t217 + MDP(6) * t257 + MDP(7) * t256 + qJDD(2) * MDP(8) + (pkin(1) * t269 + t236) * MDP(9) + (g(3) * t210 - t197 + (pkin(1) * t217 + t244) * t213) * MDP(10) + (0.2e1 * t275 + (t169 * t210 + t173 * t213) * qJD(1) + t229) * MDP(11) + (-t243 * qJDD(1) + ((t180 - t206) * t210 + (-t178 + t247) * t213) * qJD(1)) * MDP(12) + (t197 + 0.2e1 * t204 + 0.2e1 * t205 + (qJD(1) * t173 - g(3)) * t210 + (-qJD(1) * t169 - t244) * t213) * MDP(13) + (-pkin(5) * qJD(1) * t240 - t158 * pkin(2) - g(3) * t267 + t156 * qJ(3) + t180 * qJD(3) + t169 * t173 + t244 * t243) * MDP(14) + (t171 * t258 - t224 + t274) * MDP(17) + (t147 + t273) * MDP(18) + (-t241 * t202 + t209 * t151 - t212 * t150 - t159 * t165 + (t212 * t176 + t209 * t260) * t203 + (t209 * t157 + t212 * t168 + t203 * t242) * qJD(4) - t221) * MDP(20) + (t242 * t202 + t212 * t151 + t209 * t150 - t159 * t167 + (-t209 * t176 + t212 * t260) * t203 + (t212 * t157 - t209 * t168 + t203 * t241) * qJD(4) - t220) * MDP(21) - t285; (-qJDD(2) - t255) * MDP(11) + MDP(12) * t257 + (-t207 * t217 - t216) * MDP(13) + (-qJD(2) * t180 - t169 * t264 + t192 - t229 - t275) * MDP(14) + (-t165 * t264 - t212 * t202) * MDP(20) + (-t167 * t264 + t209 * t202) * MDP(21) + (-MDP(20) * t209 - MDP(21) * t212) * t203 ^ 2; (t254 - t274) * MDP(17) + (-t268 - t273) * MDP(18) + t221 * MDP(20) + t220 * MDP(21) + (-MDP(17) * t250 + MDP(18) * t227 - MDP(20) * t232 + MDP(21) * t231) * t212 + (MDP(17) * t227 - MDP(18) * t288 + MDP(20) * t231 + MDP(21) * t232) * t209 + t285;];
tau = t1;
