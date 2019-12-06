% Calculate vector of inverse dynamics joint torques for
% S5PRPRP4
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:11
% EndTime: 2019-12-05 15:36:15
% DurationCPUTime: 1.24s
% Computational Cost: add. (841->201), mult. (1636->259), div. (0->0), fcn. (1135->10), ass. (0->111)
t214 = cos(qJ(2));
t202 = t214 * qJDD(1);
t212 = sin(qJ(2));
t253 = qJD(1) * qJD(2);
t175 = qJDD(2) * pkin(2) - t212 * t253 + t202;
t207 = sin(pkin(8));
t209 = cos(pkin(8));
t283 = qJDD(1) * t212 + t214 * t253;
t152 = t207 * t175 + t283 * t209;
t284 = qJDD(2) * pkin(6) + qJD(3) * qJD(4) + t152;
t213 = cos(qJ(4));
t261 = qJD(1) * t214;
t187 = qJD(2) * pkin(2) + t261;
t262 = qJD(1) * t212;
t191 = t209 * t262;
t163 = t207 * t187 + t191;
t159 = qJD(2) * pkin(6) + t163;
t211 = sin(qJ(4));
t272 = t159 * t211;
t155 = qJD(3) * t213 - t272;
t282 = qJD(5) - t155;
t204 = qJ(2) + pkin(8);
t198 = sin(t204);
t208 = sin(pkin(7));
t210 = cos(pkin(7));
t232 = g(1) * t210 + g(2) * t208;
t281 = t198 * t232;
t150 = -qJD(4) * pkin(4) + t282;
t156 = qJD(3) * t211 + t159 * t213;
t153 = qJD(4) * qJ(5) + t156;
t205 = t211 ^ 2;
t206 = t213 ^ 2;
t263 = t205 + t206;
t280 = qJD(2) * t263;
t273 = qJDD(4) * pkin(4);
t279 = qJDD(5) - t273;
t278 = pkin(2) * t209;
t275 = g(3) * t198;
t199 = cos(t204);
t274 = g(3) * t199;
t271 = t208 * t211;
t270 = t208 * t213;
t269 = t210 * t211;
t268 = t210 * t213;
t267 = t211 * t213;
t266 = qJDD(1) - g(3);
t231 = pkin(4) * t211 - qJ(5) * t213;
t169 = qJD(4) * t231 - qJD(5) * t211;
t171 = t207 * t261 + t191;
t265 = t169 - t171;
t264 = t205 - t206;
t260 = qJD(2) * t169;
t259 = qJD(2) * t171;
t176 = t207 * t212 - t209 * t214;
t258 = qJD(2) * t176;
t257 = qJD(2) * t211;
t255 = qJD(4) * t211;
t254 = qJD(4) * t213;
t252 = qJD(2) * qJD(4);
t249 = qJDD(2) * t211;
t248 = qJDD(2) * t213;
t247 = qJDD(4) * qJ(5);
t194 = pkin(2) * t207 + pkin(6);
t246 = qJDD(4) * t194;
t245 = MDP(11) + MDP(13);
t244 = MDP(12) - MDP(15);
t216 = qJD(2) ^ 2;
t243 = t216 * t267;
t242 = t211 * qJDD(3) + t284 * t213;
t240 = -g(1) * t208 + g(2) * t210;
t238 = t155 + t272;
t190 = t207 * t262;
t162 = t187 * t209 - t190;
t227 = pkin(4) * t213 + qJ(5) * t211 + pkin(3);
t154 = -qJD(2) * t227 - t162;
t174 = -t227 - t278;
t237 = qJD(2) * t174 + t154;
t158 = -qJD(2) * pkin(3) - t162;
t195 = -pkin(3) - t278;
t236 = qJD(2) * t195 + t158;
t235 = qJDD(2) * t263;
t234 = -t213 * qJDD(3) + t159 * t254 + t284 * t211;
t173 = t209 * t261 - t190;
t233 = t173 * t255 + t213 * t259 + (g(1) * t268 + g(2) * t270) * t198;
t151 = t175 * t209 - t283 * t207;
t215 = qJD(4) ^ 2;
t230 = t194 * t215 + t274;
t143 = t247 + (qJD(5) - t272) * qJD(4) + t242;
t144 = t234 + t279;
t229 = t143 * t213 + t144 * t211;
t228 = t150 * t211 + t153 * t213;
t177 = t207 * t214 + t209 * t212;
t166 = t199 * t270 - t269;
t168 = t199 * t268 + t271;
t226 = g(1) * t168 + g(2) * t166 - t242;
t170 = t177 * qJD(2);
t225 = qJD(2) * t170 + qJDD(2) * t176 + t177 * t215;
t145 = -qJDD(2) * t227 - t151 + t260;
t224 = -qJDD(2) * t174 - t145 - t230;
t223 = -t151 + t230 + (-pkin(3) + t195) * qJDD(2);
t222 = 0.2e1 * t258 * qJD(4) - qJDD(4) * t177;
t221 = -g(3) * t214 + t212 * t232;
t165 = t199 * t271 + t268;
t167 = t199 * t269 - t270;
t220 = g(1) * t167 + g(2) * t165 + t211 * t275 - t234;
t218 = qJD(4) * t156 + t220;
t217 = (t150 * t213 - t153 * t211) * qJD(4) + t229;
t183 = qJDD(4) * t213 - t211 * t215;
t182 = qJDD(4) * t211 + t213 * t215;
t178 = t231 * qJD(2);
t1 = [t266 * MDP(1) + (qJDD(2) * t214 - t212 * t216) * MDP(3) + (-qJDD(2) * t212 - t214 * t216) * MDP(4) + (-t151 * t176 - t162 * t170 - g(3)) * MDP(5) + (t145 * t176 + t154 * t170 - g(3)) * MDP(16) + t245 * (t211 * t222 - t213 * t225) + t244 * (t211 * t225 + t213 * t222) - (MDP(14) * t280 + MDP(16) * t228 + t163 * MDP(5)) * t258 + (t152 * MDP(5) + (t150 * t254 - t153 * t255 + t229) * MDP(16) + MDP(14) * t235) * t177; qJDD(2) * MDP(2) + (t202 + t221) * MDP(3) + (-t212 * t266 + t214 * t232) * MDP(4) + (t162 * t171 - t163 * t173 + (t151 * t209 + t152 * t207 + t221) * pkin(2)) * MDP(5) + (qJDD(2) * t205 + 0.2e1 * t252 * t267) * MDP(6) + 0.2e1 * (t211 * t248 - t252 * t264) * MDP(7) + t182 * MDP(8) + t183 * MDP(9) + ((qJD(4) * t236 - t246) * t211 - t223 * t213 + t233) * MDP(11) + ((-t246 + (t173 + t236) * qJD(4)) * t213 + (t223 - t259 - t281) * t211) * MDP(12) + ((qJD(4) * t237 - t246) * t211 + (t224 - t260) * t213 + t233) * MDP(13) + (-t173 * t280 + t194 * t235 - t199 * t232 + t217 - t275) * MDP(14) + ((t246 + (-t173 - t237) * qJD(4)) * t213 + (-qJD(2) * t265 + t224 + t281) * t211) * MDP(15) + (t145 * t174 - g(3) * (pkin(2) * t214 + pkin(6) * t198) - t227 * t274 - t228 * t173 + t265 * t154 + t217 * t194 + t232 * (pkin(2) * t212 - pkin(6) * t199 + t198 * t227)) * MDP(16); (qJDD(3) + t240) * MDP(5) + (qJD(4) * t228 + t143 * t211 - t144 * t213 + t240) * MDP(16) + t245 * t183 - t244 * t182; -MDP(6) * t243 + t264 * t216 * MDP(7) + MDP(8) * t249 + MDP(9) * t248 + qJDD(4) * MDP(10) + (-t158 * t257 + t218) * MDP(11) + ((-qJD(2) * t158 + t275) * t213 + t238 * qJD(4) + t226) * MDP(12) + (0.2e1 * t273 - qJDD(5) + (-t154 * t211 + t178 * t213) * qJD(2) + t218) * MDP(13) - t231 * qJDD(2) * MDP(14) + (-t213 * t275 + 0.2e1 * t247 + (t154 * t213 + t178 * t211) * qJD(2) + (0.2e1 * qJD(5) - t238) * qJD(4) - t226) * MDP(15) + (t143 * qJ(5) - t144 * pkin(4) - t154 * t178 - t150 * t156 - g(1) * (-pkin(4) * t167 + qJ(5) * t168) - g(2) * (-pkin(4) * t165 + qJ(5) * t166) + t231 * t275 + t282 * t153) * MDP(16); (-qJDD(4) - t243) * MDP(13) + MDP(14) * t249 + (-t205 * t216 - t215) * MDP(15) + (-qJD(4) * t153 + t154 * t257 - t220 + t279) * MDP(16);];
tau = t1;
