% Calculate vector of inverse dynamics joint torques for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:21
% EndTime: 2019-12-31 16:59:23
% DurationCPUTime: 1.50s
% Computational Cost: add. (529->224), mult. (1100->271), div. (0->0), fcn. (526->4), ass. (0->112)
t204 = cos(qJ(2));
t262 = qJ(3) * t204;
t202 = sin(qJ(2));
t266 = pkin(2) + pkin(3);
t269 = t266 * t202;
t217 = t262 - t269;
t241 = qJD(1) * qJD(2);
t233 = t204 * t241;
t240 = qJDD(1) * t202;
t271 = t233 + t240;
t270 = t266 * qJDD(2);
t251 = qJD(1) * t202;
t182 = pkin(5) * t251;
t155 = qJ(4) * t251 - t182;
t243 = qJD(3) - t155;
t146 = -t266 * qJD(2) + t243;
t250 = qJD(1) * t204;
t183 = pkin(5) * t250;
t157 = -qJ(4) * t250 + t183;
t198 = qJD(2) * qJ(3);
t152 = t157 + t198;
t203 = sin(qJ(1));
t205 = cos(qJ(1));
t268 = g(1) * t205 + g(2) * t203;
t189 = t204 * pkin(2);
t184 = t202 * qJ(3);
t231 = pkin(1) + t184;
t218 = -t231 - t189;
t153 = t218 * qJD(1);
t255 = t189 + t184;
t159 = -pkin(1) - t255;
t263 = pkin(5) * qJDD(2);
t267 = (qJD(1) * t159 + t153) * qJD(2) - t263;
t194 = g(1) * t203;
t265 = g(2) * t205;
t188 = t204 * pkin(3);
t264 = pkin(5) - qJ(4);
t201 = qJDD(1) * pkin(1);
t261 = qJDD(2) * pkin(2);
t260 = t202 * t203;
t259 = t202 * t205;
t208 = qJD(1) ^ 2;
t258 = t202 * t208;
t257 = t203 * t204;
t256 = t204 * t205;
t199 = t202 ^ 2;
t200 = t204 ^ 2;
t253 = t199 - t200;
t252 = t199 + t200;
t249 = qJD(2) * t152;
t163 = t264 * t204;
t248 = qJD(2) * t163;
t247 = qJD(2) * t202;
t246 = qJD(3) * t202;
t245 = qJD(4) * t202;
t244 = qJD(4) * t204;
t144 = qJD(4) + (t266 * t204 + t231) * qJD(1);
t242 = qJD(4) + t144;
t239 = qJDD(1) * t204;
t238 = t204 * t258;
t179 = pkin(5) * t239;
t196 = qJDD(2) * qJ(3);
t197 = qJD(2) * qJD(3);
t237 = t179 + 0.2e1 * t196 + 0.2e1 * t197;
t236 = t179 + t196 + t197;
t235 = t188 + t255;
t234 = t202 * t241;
t170 = pkin(5) * t233;
t178 = pkin(5) * t240;
t232 = qJDD(3) + t170 + t178;
t230 = t194 - t265;
t229 = -qJD(2) * pkin(2) + qJD(3);
t154 = pkin(1) + t235;
t228 = qJD(1) * t154 + t144;
t226 = t205 * pkin(1) + pkin(2) * t256 + t203 * pkin(5) + qJ(3) * t259;
t225 = g(1) * t259 + g(2) * t260 - g(3) * t204 - t178;
t207 = qJD(2) ^ 2;
t224 = pkin(5) * t207 + t265;
t222 = pkin(2) * t202 - t262;
t221 = pkin(2) * t239 + t271 * qJ(3) + qJD(1) * t246 + t201;
t158 = t182 + t229;
t160 = t183 + t198;
t220 = t158 * t204 - t160 * t202;
t219 = -qJDD(3) + t225;
t216 = -0.2e1 * pkin(1) * t241 - t263;
t215 = -t170 + t219;
t214 = -t224 + 0.2e1 * t201;
t213 = pkin(3) * t239 + qJDD(4) + t221;
t139 = -t266 * t234 + t213;
t143 = t217 * qJD(2) + t246;
t212 = qJD(1) * t143 + qJDD(1) * t154 + t139 - t265;
t211 = -t215 - t261;
t142 = pkin(2) * t234 - t221;
t151 = t222 * qJD(2) - t246;
t210 = -qJD(1) * t151 - qJDD(1) * t159 - t142 - t224;
t145 = -pkin(5) * t234 + t236;
t147 = t232 - t261;
t209 = t220 * qJD(2) + t145 * t204 + t147 * t202;
t190 = t205 * pkin(5);
t174 = g(1) * t257;
t173 = g(1) * t260;
t169 = qJ(3) * t256;
t167 = qJ(3) * t257;
t165 = qJ(4) * t234;
t162 = t264 * t202;
t156 = t222 * qJD(1);
t150 = -t245 + t248;
t149 = -t264 * t247 - t244;
t148 = t217 * qJD(1);
t141 = -qJ(4) * t239 + t165 + (-pkin(5) * t247 - t244) * qJD(1) + t236;
t140 = -t271 * qJ(4) - qJD(1) * t245 + t232 - t270;
t1 = [qJDD(1) * MDP(1) + t230 * MDP(2) + t268 * MDP(3) + (qJDD(1) * t199 + 0.2e1 * t202 * t233) * MDP(4) + 0.2e1 * (t202 * t239 - t253 * t241) * MDP(5) + (qJDD(2) * t202 + t204 * t207) * MDP(6) + (qJDD(2) * t204 - t202 * t207) * MDP(7) + (t216 * t202 + t214 * t204 + t174) * MDP(9) + (-t214 * t202 + t216 * t204 - t173) * MDP(10) + (t267 * t202 + t210 * t204 + t174) * MDP(11) + (t252 * qJDD(1) * pkin(5) + t209 - t268) * MDP(12) + (t210 * t202 - t267 * t204 + t173) * MDP(13) + (t209 * pkin(5) - g(1) * t190 - g(2) * t226 + t142 * t159 + t153 * t151 - t218 * t194) * MDP(14) + (-qJDD(2) * t162 + t174 + (-t228 * t202 - t150) * qJD(2) + t212 * t204) * MDP(15) + (qJDD(2) * t163 + t173 + (t228 * t204 + t149) * qJD(2) + t212 * t202) * MDP(16) + ((-qJD(2) * t146 - qJDD(1) * t163 - t141 + (-qJD(2) * t162 - t149) * qJD(1)) * t204 + (t249 - qJDD(1) * t162 - t140 + (-t150 + t248) * qJD(1)) * t202 + t268) * MDP(17) + (t141 * t163 + t152 * t149 + t140 * t162 + t146 * t150 + t139 * t154 + t144 * t143 - g(1) * (-qJ(4) * t205 + t190) - g(2) * (pkin(3) * t256 + t226) + (-g(1) * (t218 - t188) + g(2) * qJ(4)) * t203) * MDP(18); -MDP(4) * t238 + t253 * t208 * MDP(5) + MDP(6) * t240 + MDP(7) * t239 + qJDD(2) * MDP(8) + (pkin(1) * t258 + t225) * MDP(9) + (g(3) * t202 - t179 + (pkin(1) * t208 + t268) * t204) * MDP(10) + (0.2e1 * t261 + (-t153 * t202 + t156 * t204) * qJD(1) + t219) * MDP(11) + (-t222 * qJDD(1) + ((t160 - t198) * t202 + (-t158 + t229) * t204) * qJD(1)) * MDP(12) + ((qJD(1) * t156 - g(3)) * t202 + (qJD(1) * t153 - t268) * t204 + t237) * MDP(13) + (t145 * qJ(3) + t160 * qJD(3) - t147 * pkin(2) - t153 * t156 - g(1) * (-pkin(2) * t259 + t169) - g(2) * (-pkin(2) * t260 + t167) - g(3) * t255 - t220 * qJD(1) * pkin(5)) * MDP(14) + (qJ(4) * t240 + qJD(2) * t157 + 0.2e1 * t270 + ((qJ(4) * qJD(2) - t148) * t204 + t242 * t202) * qJD(1) + t215) * MDP(15) + (-qJD(2) * t155 + t165 + (-g(3) + (-pkin(5) * qJD(2) - t148) * qJD(1)) * t202 + (-qJ(4) * qJDD(1) - t242 * qJD(1) - t268) * t204 + t237) * MDP(16) - t217 * qJDD(1) * MDP(17) + (-g(1) * t169 - g(2) * t167 - g(3) * t235 + t141 * qJ(3) - t140 * t266 - t144 * t148 - t146 * t157 + t243 * t152 + t268 * t269) * MDP(18); (-qJD(2) * t160 + t211) * MDP(14) + (-qJDD(2) * pkin(3) - qJ(4) * t233 + t211 - t249) * MDP(18) + (MDP(13) + MDP(16)) * (-t199 * t208 - t207) + (MDP(11) + MDP(15)) * (-qJDD(2) - t238) + ((-MDP(18) * qJ(4) + MDP(12) - MDP(17)) * qJDD(1) + (t153 * MDP(14) - t242 * MDP(18)) * qJD(1)) * t202; (t213 + t230) * MDP(18) - t252 * MDP(17) * t208 + (MDP(15) * t204 + t202 * MDP(16)) * qJDD(1) + ((t146 * t202 + t152 * t204) * MDP(18) + (0.2e1 * t204 * MDP(16) + (-t266 * MDP(18) - 0.2e1 * MDP(15)) * t202) * qJD(2)) * qJD(1);];
tau = t1;
