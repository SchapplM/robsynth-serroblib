% Calculate vector of inverse dynamics joint torques for
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:58
% EndTime: 2019-12-31 18:06:02
% DurationCPUTime: 1.68s
% Computational Cost: add. (700->246), mult. (1253->315), div. (0->0), fcn. (713->6), ass. (0->110)
t191 = cos(qJ(4));
t226 = qJDD(1) * t191;
t267 = qJD(4) * qJD(5) + t226;
t189 = sin(qJ(1));
t192 = cos(qJ(1));
t244 = g(1) * t189 - g(2) * t192;
t186 = pkin(1) + qJ(3);
t266 = qJD(1) * t186;
t179 = qJDD(1) * qJ(2);
t180 = qJD(1) * qJD(2);
t215 = qJDD(3) + t179 + t180;
t156 = -pkin(6) * qJDD(1) + t215;
t166 = qJ(2) * qJD(1) + qJD(3);
t161 = -pkin(6) * qJD(1) + t166;
t188 = sin(qJ(4));
t238 = qJD(4) * t188;
t142 = -qJDD(4) * pkin(4) - t156 * t191 + t161 * t238;
t253 = t161 * t191;
t146 = -qJD(4) * pkin(4) - t253;
t231 = qJD(1) * qJD(4);
t218 = t191 * t231;
t227 = qJDD(1) * t188;
t151 = qJDD(5) + t218 + t227;
t208 = pkin(4) * t188 - pkin(7) * t191;
t159 = t208 + t186;
t164 = qJD(1) * t188 + qJD(5);
t185 = -pkin(6) + qJ(2);
t144 = t159 * qJD(1) - qJD(2);
t237 = qJD(4) * t191;
t211 = -qJDD(4) * pkin(7) - qJD(5) * t144 - t156 * t188 - t161 * t237;
t265 = -(qJD(5) * t159 + t185 * t237) * t164 + t142 * t191 + (-qJD(2) * t164 - t146 * qJD(4) - t185 * t151 + t211) * t188;
t162 = -qJD(2) + t266;
t264 = (qJD(2) + t162 + t266) * qJD(4) + qJDD(4) * t185;
t207 = g(1) * t192 + g(2) * t189;
t209 = pkin(4) * t191 + pkin(7) * t188;
t262 = g(3) * t188;
t263 = (pkin(7) * qJD(5) + t209 * qJD(1)) * t164 + t207 * t191 + t142 - t262;
t174 = 0.2e1 * t180;
t261 = g(3) * t191;
t187 = sin(qJ(5));
t190 = cos(qJ(5));
t232 = t190 * qJD(4);
t234 = qJD(5) * t191;
t199 = -t187 * t234 - t188 * t232;
t221 = t187 * qJDD(4) + t190 * t267;
t140 = t199 * qJD(1) + t221;
t260 = t140 * t187;
t241 = qJD(1) * t191;
t153 = t187 * t241 - t232;
t259 = t153 * t164;
t258 = t153 * t191;
t233 = t187 * qJD(4);
t155 = t190 * t241 + t233;
t257 = t155 * t164;
t256 = t155 * t191;
t255 = t159 * t151;
t254 = t161 * t188;
t252 = t187 * t188;
t251 = t187 * t189;
t250 = t187 * t192;
t249 = t189 * t190;
t248 = t190 * t164;
t247 = t190 * t192;
t246 = -t190 * qJDD(4) - t231 * t252;
t245 = t192 * pkin(1) + t189 * qJ(2);
t183 = t191 ^ 2;
t243 = t188 ^ 2 - t183;
t193 = qJD(4) ^ 2;
t194 = qJD(1) ^ 2;
t242 = -t193 - t194;
t240 = qJD(4) * t153;
t239 = qJD(4) * t155;
t145 = qJD(4) * pkin(7) + t254;
t236 = qJD(5) * t145;
t235 = qJD(5) * t164;
t230 = qJD(3) * qJD(1);
t228 = qJDD(1) * t186;
t184 = qJDD(1) * pkin(1);
t224 = qJDD(2) - t184;
t223 = t164 * t252;
t222 = t188 * t248;
t220 = t164 * t233;
t219 = t164 * t232;
t217 = qJDD(2) - t244;
t178 = qJDD(1) * qJ(3);
t216 = -t178 + t224;
t152 = t209 * qJD(4) + qJD(3);
t139 = t152 * qJD(1) + t208 * qJDD(1) - t216;
t212 = -t139 + t236;
t210 = -t184 + t217;
t205 = -t178 + t210;
t204 = t211 + t261;
t203 = t187 * t151 + t190 * t235;
t202 = -t190 * t151 + t187 * t235;
t200 = 0.2e1 * t179 + t174 - t207;
t198 = qJD(1) * t162 - t156 + t207;
t197 = -pkin(7) * t151 + (t146 + t253) * t164;
t157 = -t216 + t230;
t196 = -t185 * t193 + t157 + t228 + t230 + t244;
t172 = t192 * qJ(2);
t169 = qJDD(4) * t191;
t150 = t188 * t247 - t251;
t149 = -t188 * t250 - t249;
t148 = -t188 * t249 - t250;
t147 = t188 * t251 - t247;
t141 = t155 * qJD(5) + t187 * t226 + t246;
t138 = t190 * t139;
t137 = t187 * t144 + t190 * t145;
t136 = t190 * t144 - t187 * t145;
t1 = [qJDD(1) * MDP(1) + t244 * MDP(2) + t207 * MDP(3) + (-0.2e1 * t184 + t217) * MDP(4) + t200 * MDP(5) + (-t224 * pkin(1) - g(1) * (-pkin(1) * t189 + t172) - g(2) * t245 + (t179 + t174) * qJ(2)) * MDP(6) + (qJDD(3) + t200) * MDP(7) + (-t205 + t228 + 0.2e1 * t230) * MDP(8) + (t157 * t186 + t162 * qJD(3) + t215 * qJ(2) + t166 * qJD(2) - g(1) * (-t186 * t189 + t172) - g(2) * (qJ(3) * t192 + t245)) * MDP(9) + (qJDD(1) * t183 - 0.2e1 * t188 * t218) * MDP(10) + 0.2e1 * (-t188 * t226 + t243 * t231) * MDP(11) + (-t188 * t193 + t169) * MDP(12) + (-qJDD(4) * t188 - t191 * t193) * MDP(13) + (t196 * t188 + t264 * t191) * MDP(15) + (-t264 * t188 + t196 * t191) * MDP(16) + (t140 * t190 * t191 + t199 * t155) * MDP(17) + ((t153 * t190 + t155 * t187) * t238 + (-t260 - t141 * t190 + (t153 * t187 - t155 * t190) * qJD(5)) * t191) * MDP(18) + ((t140 - t219) * t188 + (-t202 + t239) * t191) * MDP(19) + ((-t141 + t220) * t188 + (-t203 - t240) * t191) * MDP(20) + (t151 * t188 + t164 * t237) * MDP(21) + (-g(1) * t148 - g(2) * t150 + (t185 * t240 + t138) * t188 + (-qJD(2) * t153 + qJD(4) * t136 - t141 * t185) * t191 + (t255 + t152 * t164 + (t146 * t191 + (-t164 * t185 - t145) * t188) * qJD(5)) * t190 + t265 * t187) * MDP(22) + (t185 * t155 * t238 - g(1) * t147 - g(2) * t149 + (-qJD(2) * t155 - qJD(4) * t137 - t140 * t185) * t191 + (-(-qJD(5) * t185 * t188 + t152) * t164 - t255 + t212 * t188 - t146 * t234) * t187 + t265 * t190) * MDP(23); t210 * MDP(6) + t205 * MDP(9) + t202 * MDP(22) + t203 * MDP(23) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t194 + (-MDP(15) * t188 - MDP(16) * t191 + MDP(4) - MDP(8)) * qJDD(1) + ((-qJD(3) - t166) * MDP(9) + (t223 + t258) * MDP(22) + (t222 + t256) * MDP(23) + 0.2e1 * (-MDP(15) * t191 + MDP(16) * t188) * qJD(4)) * qJD(1); qJDD(1) * MDP(7) - t194 * MDP(8) + (-t207 + t215) * MDP(9) + t169 * MDP(15) + (-t162 * MDP(9) + (-MDP(22) * t190 + MDP(23) * t187) * t164) * qJD(1) + (t242 * MDP(16) + (-t141 - t220) * MDP(22) + (-t140 - t219) * MDP(23)) * t191 + (t242 * MDP(15) - qJDD(4) * MDP(16) + (-t203 + t240) * MDP(22) + (t202 + t239) * MDP(23)) * t188; MDP(12) * t226 - MDP(13) * t227 + qJDD(4) * MDP(14) + (-t198 * t191 + t262) * MDP(15) + (t198 * t188 + t261) * MDP(16) + (t155 * t248 + t260) * MDP(17) + ((t140 - t259) * t190 + (-t141 - t257) * t187) * MDP(18) + ((t222 - t256) * qJD(1) + t203) * MDP(19) + ((-t223 + t258) * qJD(1) - t202) * MDP(20) - t164 * MDP(21) * t241 + (-pkin(4) * t141 - t136 * t241 - t153 * t254 + t197 * t187 - t190 * t263) * MDP(22) + (-pkin(4) * t140 + t137 * t241 - t155 * t254 + t187 * t263 + t197 * t190) * MDP(23) + (t191 * t188 * MDP(10) - t243 * MDP(11)) * t194; t155 * t153 * MDP(17) + (-t153 ^ 2 + t155 ^ 2) * MDP(18) + (t221 + t259) * MDP(19) + (-t246 + t257) * MDP(20) + t151 * MDP(21) + (-g(1) * t149 + g(2) * t147 + t137 * t164 - t146 * t155 + t138) * MDP(22) + (g(1) * t150 - g(2) * t148 + t136 * t164 + t146 * t153) * MDP(23) + (-MDP(22) * t236 + t204 * MDP(23) + (-MDP(19) * t238 - MDP(20) * t234) * qJD(1)) * t190 + (-qJD(1) * MDP(19) * t234 - MDP(20) * t267 + t204 * MDP(22) + t212 * MDP(23)) * t187;];
tau = t1;
