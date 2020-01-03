% Calculate vector of inverse dynamics joint torques for
% S4RRPR7
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
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:48
% EndTime: 2019-12-31 17:06:50
% DurationCPUTime: 1.66s
% Computational Cost: add. (1116->239), mult. (2683->335), div. (0->0), fcn. (1873->10), ass. (0->115)
t218 = sin(pkin(7));
t219 = cos(pkin(7));
t222 = sin(qJ(2));
t225 = cos(qJ(2));
t195 = t218 * t225 + t219 * t222;
t256 = qJD(1) * qJD(2);
t248 = t225 * t256;
t249 = t222 * t256;
t169 = qJDD(1) * t195 - t218 * t249 + t219 * t248;
t280 = qJD(2) * qJD(4) + t169;
t275 = qJ(3) + pkin(5);
t246 = qJD(2) * t275;
t234 = -qJD(3) * t222 - t225 * t246;
t250 = t275 * t222;
t163 = qJDD(2) * pkin(2) + qJD(1) * t234 - qJDD(1) * t250;
t185 = qJD(3) * t225 - t222 * t246;
t200 = t275 * t225;
t170 = qJD(1) * t185 + qJDD(1) * t200;
t147 = t163 * t219 - t170 * t218;
t145 = -qJDD(2) * pkin(3) - t147;
t153 = t219 * t185 + t218 * t234;
t197 = qJD(1) * t250;
t274 = qJD(2) * pkin(2);
t193 = -t197 + t274;
t198 = qJD(1) * t200;
t266 = t198 * t218;
t166 = t193 * t219 - t266;
t158 = -qJD(2) * pkin(3) - t166;
t187 = t195 * qJD(2);
t253 = qJDD(1) * t225;
t254 = qJDD(1) * t222;
t238 = -t218 * t254 + t219 * t253;
t164 = qJD(1) * t187 + qJDD(4) - t238;
t265 = t219 * t225;
t194 = t218 * t222 - t265;
t209 = pkin(2) * t225 + pkin(1);
t165 = pkin(3) * t194 - pkin(6) * t195 - t209;
t174 = t219 * t200 - t218 * t250;
t259 = qJD(1) * t222;
t186 = qJD(1) * t265 - t218 * t259;
t180 = qJD(4) - t186;
t190 = t194 * qJD(2);
t148 = t218 * t163 + t219 * t170;
t146 = qJDD(2) * pkin(6) + t148;
t188 = t195 * qJD(1);
t199 = -qJD(1) * t209 + qJD(3);
t151 = -pkin(3) * t186 - pkin(6) * t188 + t199;
t243 = qJD(4) * t151 + t146;
t279 = t145 * t195 - t158 * t190 - t174 * t164 - (qJD(4) * t165 + t153) * t180 - t194 * t243;
t206 = pkin(2) * t218 + pkin(6);
t215 = qJ(2) + pkin(7);
t210 = sin(t215);
t211 = cos(t215);
t223 = sin(qJ(1));
t226 = cos(qJ(1));
t240 = g(1) * t226 + g(2) * t223;
t278 = t180 * (pkin(2) * t259 + pkin(3) * t188 - pkin(6) * t186 + qJD(4) * t206) - t210 * t240 + g(3) * t211 + t145;
t277 = g(3) * t210;
t276 = g(3) * t225;
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t251 = t221 * qJDD(2) + t224 * t280;
t257 = qJD(4) * t221;
t149 = -t188 * t257 + t251;
t273 = t149 * t221;
t272 = t164 * t221;
t271 = t165 * t164;
t175 = -t224 * qJD(2) + t188 * t221;
t270 = t175 * t180;
t269 = t175 * t188;
t177 = qJD(2) * t221 + t188 * t224;
t268 = t177 * t180;
t267 = t177 * t188;
t191 = t219 * t198;
t264 = t221 * t223;
t263 = t221 * t226;
t262 = t223 * t224;
t157 = t224 * t164;
t261 = t224 * t226;
t167 = t218 * t193 + t191;
t216 = t222 ^ 2;
t260 = -t225 ^ 2 + t216;
t258 = qJD(4) * t195;
t252 = t222 * t274;
t245 = t180 * t224;
t168 = -qJD(2) * t188 + t238;
t233 = pkin(2) * t249 - qJDD(1) * t209 + qJDD(3);
t144 = -pkin(3) * t168 - pkin(6) * t169 + t233;
t159 = qJD(2) * pkin(6) + t167;
t244 = qJD(4) * t159 - t144;
t239 = g(1) * t223 - g(2) * t226;
t237 = t157 + (t186 * t221 - t257) * t180;
t236 = -0.2e1 * pkin(1) * t256 - pkin(5) * qJDD(2);
t235 = -t190 * t224 - t195 * t257;
t172 = -t197 * t219 - t266;
t231 = -t206 * t164 + (t158 + t172) * t180;
t227 = qJD(2) ^ 2;
t230 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t227 + t239;
t228 = qJD(1) ^ 2;
t229 = pkin(1) * t228 - pkin(5) * qJDD(1) + t240;
t213 = t224 * qJDD(2);
t207 = -pkin(2) * t219 - pkin(3);
t184 = t211 * t261 + t264;
t183 = -t211 * t263 + t262;
t182 = -t211 * t262 + t263;
t181 = t211 * t264 + t261;
t173 = t200 * t218 + t219 * t250;
t171 = -t197 * t218 + t191;
t155 = pkin(3) * t187 + pkin(6) * t190 + t252;
t152 = t185 * t218 - t219 * t234;
t150 = t177 * qJD(4) + t169 * t221 - t213;
t143 = t151 * t221 + t159 * t224;
t142 = t151 * t224 - t159 * t221;
t141 = t224 * t144;
t1 = [qJDD(1) * MDP(1) + t239 * MDP(2) + t240 * MDP(3) + (qJDD(1) * t216 + 0.2e1 * t222 * t248) * MDP(4) + 0.2e1 * (t222 * t253 - t256 * t260) * MDP(5) + (qJDD(2) * t222 + t225 * t227) * MDP(6) + (qJDD(2) * t225 - t222 * t227) * MDP(7) + (t222 * t236 + t225 * t230) * MDP(9) + (-t222 * t230 + t225 * t236) * MDP(10) + (-t147 * t195 - t148 * t194 + t152 * t188 + t153 * t186 + t166 * t190 - t167 * t187 + t168 * t174 + t169 * t173 - t240) * MDP(11) + (t148 * t174 + t167 * t153 - t147 * t173 - t166 * t152 - t233 * t209 + t199 * t252 - g(1) * (-t209 * t223 + t226 * t275) - g(2) * (t209 * t226 + t223 * t275)) * MDP(12) + (t149 * t195 * t224 + t177 * t235) * MDP(13) + (-(-t175 * t224 - t177 * t221) * t190 + (-t273 - t150 * t224 + (t175 * t221 - t177 * t224) * qJD(4)) * t195) * MDP(14) + (t149 * t194 + t157 * t195 + t177 * t187 + t180 * t235) * MDP(15) + (-t195 * t272 - t150 * t194 - t175 * t187 + (t190 * t221 - t224 * t258) * t180) * MDP(16) + (t164 * t194 + t180 * t187) * MDP(17) + (-g(1) * t182 - g(2) * t184 + t141 * t194 + t142 * t187 + t173 * t150 + t152 * t175 + (t155 * t180 + t271 + (t158 * t195 - t159 * t194 - t174 * t180) * qJD(4)) * t224 + t279 * t221) * MDP(18) + (-g(1) * t181 - g(2) * t183 - t143 * t187 + t173 * t149 + t152 * t177 + (-(-qJD(4) * t174 + t155) * t180 - t271 + t244 * t194 - t158 * t258) * t221 + t279 * t224) * MDP(19); MDP(6) * t254 + MDP(7) * t253 + qJDD(2) * MDP(8) + (t222 * t229 - t276) * MDP(9) + (g(3) * t222 + t229 * t225) * MDP(10) + ((t167 - t171) * t188 + (t166 - t172) * t186 + (t168 * t218 - t169 * t219) * pkin(2)) * MDP(11) + (t166 * t171 - t167 * t172 + (-t276 + t147 * t219 + t148 * t218 + (-qJD(1) * t199 + t240) * t222) * pkin(2)) * MDP(12) + (t177 * t245 + t273) * MDP(13) + ((t149 - t270) * t224 + (-t150 - t268) * t221) * MDP(14) + (t180 * t245 - t267 + t272) * MDP(15) + (t237 + t269) * MDP(16) - t180 * t188 * MDP(17) + (-t142 * t188 + t207 * t150 - t171 * t175 + t231 * t221 - t224 * t278) * MDP(18) + (t143 * t188 + t207 * t149 - t171 * t177 + t221 * t278 + t231 * t224) * MDP(19) + (-t222 * t225 * MDP(4) + t260 * MDP(5)) * t228; (-t186 ^ 2 - t188 ^ 2) * MDP(11) + (t166 * t188 - t167 * t186 + t233 - t239) * MDP(12) + (t237 - t269) * MDP(18) + (-t180 ^ 2 * t224 - t267 - t272) * MDP(19); t177 * t175 * MDP(13) + (-t175 ^ 2 + t177 ^ 2) * MDP(14) + (t251 + t270) * MDP(15) + (t213 + t268) * MDP(16) + t164 * MDP(17) + (-g(1) * t183 + g(2) * t181 + t143 * t180 - t158 * t177 + t141) * MDP(18) + (g(1) * t184 - g(2) * t182 + t142 * t180 + t158 * t175) * MDP(19) + ((-t146 + t277) * MDP(19) + (-MDP(16) * t188 - MDP(18) * t159 - MDP(19) * t151) * qJD(4)) * t224 + (-qJD(4) * t188 * MDP(15) - t280 * MDP(16) + (-t243 + t277) * MDP(18) + t244 * MDP(19)) * t221;];
tau = t1;
