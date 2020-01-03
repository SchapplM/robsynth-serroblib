% Calculate vector of inverse dynamics joint torques for
% S4RRPP5
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
%   see S4RRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:40
% EndTime: 2019-12-31 17:00:42
% DurationCPUTime: 1.61s
% Computational Cost: add. (530->219), mult. (1092->250), div. (0->0), fcn. (518->4), ass. (0->108)
t197 = sin(qJ(2));
t184 = t197 * qJ(3);
t199 = cos(qJ(2));
t242 = t199 * pkin(2) + t184;
t268 = -pkin(1) - t242;
t148 = t268 * qJD(1);
t269 = qJDD(1) * t268;
t262 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t239 = qJD(1) * t199;
t253 = pkin(2) + qJ(4);
t267 = qJD(2) * t253;
t180 = pkin(5) * t239;
t153 = pkin(3) * t239 + t180;
t266 = -qJD(4) - t153;
t265 = t253 * qJDD(2);
t235 = qJD(1) * qJD(2);
t229 = t199 * t235;
t234 = qJDD(1) * t197;
t264 = (-t229 - t234) * pkin(3);
t263 = 0.2e1 * t262;
t238 = qJD(2) * qJ(3);
t147 = t238 - t266;
t198 = sin(qJ(1));
t200 = cos(qJ(1));
t221 = g(1) * t200 + g(2) * t198;
t247 = t197 * t200;
t248 = t197 * t198;
t261 = -g(1) * t247 - g(2) * t248 + g(3) * t199;
t210 = -t199 * t253 - pkin(1) - t184;
t140 = t210 * qJD(1);
t260 = t140 * t239 - t221 * t199;
t259 = pkin(3) + pkin(5);
t258 = g(1) * t198;
t255 = g(2) * t200;
t254 = g(3) * t197;
t252 = qJD(2) * pkin(2);
t251 = pkin(5) * qJDD(2);
t250 = qJ(4) * t199;
t249 = qJDD(2) * pkin(2);
t202 = qJD(1) ^ 2;
t246 = t197 * t202;
t245 = t198 * t199;
t244 = t199 * t200;
t233 = qJDD(1) * t199;
t243 = qJ(3) * t233 + qJD(3) * t239;
t194 = t197 ^ 2;
t195 = t199 ^ 2;
t241 = t194 - t195;
t240 = qJD(1) * t197;
t160 = t259 * t199;
t154 = qJD(2) * t160;
t237 = qJD(3) * t197;
t179 = pkin(5) * t240;
t151 = -pkin(3) * t240 - t179;
t236 = qJD(3) - t151;
t232 = t199 * t246;
t231 = t259 * qJD(2);
t230 = t197 * t235;
t165 = pkin(5) * t229;
t175 = pkin(5) * t234;
t228 = qJDD(3) + t165 + t175;
t177 = pkin(5) * t233;
t227 = pkin(3) * t233 + qJDD(4) + t177;
t220 = t242 + t250;
t149 = -pkin(1) - t220;
t226 = qJD(1) * t149 + t140;
t224 = pkin(1) * t200 + pkin(2) * t244 + pkin(5) * t198 + qJ(3) * t247;
t223 = t175 + t261;
t152 = t197 * t231;
t201 = qJD(2) ^ 2;
t222 = pkin(5) * t201 + t255;
t218 = -qJ(3) * t199 + qJ(4) * t197;
t157 = -t180 - t238;
t217 = (qJD(3) + t179 - t252) * t199 + t157 * t197;
t216 = t227 + t262;
t215 = -qJDD(3) - t223;
t213 = -0.2e1 * pkin(1) * t235 - t251;
t212 = -t165 + t215;
t211 = -t199 * t238 - t237;
t209 = 0.2e1 * qJDD(1) * pkin(1) - t222;
t208 = -0.2e1 * qJD(2) * t148 + t251;
t166 = pkin(2) * t230;
t203 = qJD(2) * t218 - qJD(4) * t199 - t237;
t134 = qJD(1) * t203 + qJDD(1) * t210 + t166;
t182 = t197 * t252;
t138 = t182 + t203;
t207 = -qJD(1) * t138 - qJDD(1) * t149 - t134 - t255;
t206 = t228 - t264 - t265;
t136 = qJD(1) * t211 + t166 + t269;
t146 = t182 + t211;
t205 = qJD(1) * t146 + t136 + t222 + t269;
t141 = pkin(5) * t230 - t177 - t262;
t144 = t228 - t249;
t204 = qJD(2) * t217 - t141 * t199 + t144 * t197;
t189 = t200 * pkin(5);
t183 = pkin(2) * t240;
t171 = g(1) * t245;
t170 = g(1) * t248;
t164 = qJ(3) * t244;
t162 = qJ(3) * t245;
t159 = t259 * t197;
t150 = -qJ(3) * t239 + t183;
t145 = qJD(1) * t218 + t183;
t143 = t236 - t267;
t142 = t148 * t240;
t137 = -qJD(1) * t152 + t216;
t135 = -qJD(2) * qJD(4) + t206;
t1 = [qJDD(1) * MDP(1) + (-t255 + t258) * MDP(2) + t221 * MDP(3) + (qJDD(1) * t194 + 0.2e1 * t197 * t229) * MDP(4) + 0.2e1 * (t197 * t233 - t235 * t241) * MDP(5) + (qJDD(2) * t197 + t199 * t201) * MDP(6) + (qJDD(2) * t199 - t197 * t201) * MDP(7) + (t197 * t213 + t199 * t209 + t171) * MDP(9) + (-t197 * t209 + t199 * t213 - t170) * MDP(10) + ((t194 + t195) * qJDD(1) * pkin(5) + t204 - t221) * MDP(11) + (t197 * t208 + t199 * t205 - t171) * MDP(12) + (-t197 * t205 + t199 * t208 + t170) * MDP(13) + (pkin(5) * t204 - g(1) * t189 - g(2) * t224 + t148 * t146 + (t136 - t258) * t268) * MDP(14) + ((qJD(2) * t143 + qJDD(1) * t160 + t137 + (qJD(2) * t159 - t152) * qJD(1)) * t199 + (-qJD(2) * t147 + qJDD(1) * t159 + t135) * t197 - t221) * MDP(15) + (qJDD(2) * t160 + t170 + (-t199 * t226 - t152) * qJD(2) + t207 * t197) * MDP(16) + (-qJDD(2) * t159 + t171 + (t197 * t226 - t154) * qJD(2) + t207 * t199) * MDP(17) + (t134 * t149 + t140 * t138 + t135 * t159 + t143 * t154 + t137 * t160 - t147 * t152 - g(1) * (pkin(3) * t200 + t189) - g(2) * (qJ(4) * t244 + t224) + (-g(1) * (t268 - t250) - g(2) * pkin(3)) * t198) * MDP(18); -MDP(4) * t232 + t241 * t202 * MDP(5) + MDP(6) * t234 + MDP(7) * t233 + qJDD(2) * MDP(8) + (pkin(1) * t246 - t223) * MDP(9) + (t254 - t177 + (pkin(1) * t202 + t221) * t199) * MDP(10) + (-pkin(2) * t234 + (-qJD(2) * t242 - t217) * qJD(1) + t243) * MDP(11) + (-t150 * t239 + t142 - t215 - 0.2e1 * t249) * MDP(12) + (t177 + (qJD(1) * t150 - g(3)) * t197 + (qJD(1) * t148 - t221) * t199 + t263) * MDP(13) + (-t141 * qJ(3) - t157 * qJD(3) - t144 * pkin(2) - t148 * t150 - g(1) * (-pkin(2) * t247 + t164) - g(2) * (-pkin(2) * t248 + t162) - g(3) * t242 - t217 * qJD(1) * pkin(5)) * MDP(14) + (-t253 * t234 + (-t143 - t151 - t267) * t239 + t243) * MDP(15) + (-qJD(2) * t151 + (-g(3) + (t145 - t231) * qJD(1)) * t197 + t227 + t260 + t263) * MDP(16) + ((0.2e1 * qJD(4) + t153) * qJD(2) + (-t140 * t197 + t145 * t199) * qJD(1) + 0.2e1 * t265 + t212 + t264) * MDP(17) + (-g(1) * t164 - g(2) * t162 - g(3) * t220 + t137 * qJ(3) - t140 * t145 + t266 * t143 + t236 * t147 + (t197 * t221 - t135) * t253) * MDP(18); (qJD(2) * t157 + t142 - t212 - t249) * MDP(14) + (t140 * t240 + (-qJD(4) - t147) * qJD(2) + t206 + t261) * MDP(18) + (MDP(13) + MDP(16)) * (-t194 * t202 - t201) + (MDP(12) - MDP(17)) * (qJDD(2) + t232) + (MDP(11) + MDP(15)) * t234; qJDD(2) * MDP(16) + (-t195 * t202 - t201) * MDP(17) + (MDP(15) * qJDD(1) - MDP(16) * t246) * t199 + (t216 - t254 + (-t240 * t259 + t143) * qJD(2) + t260) * MDP(18);];
tau = t1;
