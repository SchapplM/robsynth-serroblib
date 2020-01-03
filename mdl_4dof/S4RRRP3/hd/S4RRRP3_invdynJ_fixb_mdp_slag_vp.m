% Calculate vector of inverse dynamics joint torques for
% S4RRRP3
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:18
% EndTime: 2019-12-31 17:14:19
% DurationCPUTime: 0.96s
% Computational Cost: add. (775->199), mult. (1151->230), div. (0->0), fcn. (552->8), ass. (0->108)
t192 = qJ(1) + qJ(2);
t186 = cos(t192);
t179 = g(1) * t186;
t189 = qJD(1) + qJD(2);
t193 = sin(qJ(3));
t190 = t193 ^ 2;
t196 = cos(qJ(3));
t191 = t196 ^ 2;
t241 = t190 + t191;
t268 = t189 * t241;
t188 = qJDD(1) + qJDD(2);
t194 = sin(qJ(2));
t232 = qJDD(1) * t194;
t197 = cos(qJ(2));
t238 = qJD(2) * t197;
t151 = pkin(6) * t188 + (qJD(1) * t238 + t232) * pkin(1);
t146 = t196 * t151;
t231 = qJDD(3) * qJ(4);
t259 = pkin(1) * qJD(1);
t229 = t194 * t259;
t162 = pkin(6) * t189 + t229;
t255 = t162 * t193;
t138 = t231 + t146 + (qJD(4) - t255) * qJD(3);
t145 = t193 * t151;
t234 = qJD(3) * t196;
t256 = qJDD(3) * pkin(3);
t266 = t162 * t234 - t256;
t139 = qJDD(4) + t145 + t266;
t267 = t138 * t196 + t139 * t193;
t185 = sin(t192);
t244 = g(2) * t185 + t179;
t265 = pkin(1) * t197;
t264 = pkin(2) * t188;
t263 = pkin(2) * t189;
t262 = pkin(3) * t196;
t199 = qJD(3) ^ 2;
t261 = pkin(6) * t199;
t178 = g(1) * t185;
t260 = g(2) * t186;
t258 = qJD(3) * pkin(3);
t257 = pkin(6) * qJDD(3);
t254 = t162 * t196;
t214 = qJ(4) * t193 + t262;
t164 = -pkin(2) - t214;
t253 = t164 * t188;
t252 = t164 * t189;
t180 = pkin(1) * t194 + pkin(6);
t251 = t180 * t199;
t250 = t185 * t193;
t249 = t186 * t193;
t248 = t188 * t193;
t247 = t189 * t193;
t246 = t193 * t196;
t245 = g(1) * t250 - g(2) * t249;
t239 = qJD(2) * t194;
t227 = pkin(1) * t239;
t243 = -qJD(1) * t227 + qJDD(1) * t265;
t242 = t190 - t191;
t240 = qJD(1) * t197;
t237 = qJD(3) * qJ(4);
t236 = qJD(3) * t189;
t235 = qJD(3) * t193;
t233 = qJD(4) * t193;
t230 = qJDD(3) * t180;
t228 = pkin(1) * t240;
t226 = pkin(1) * t238;
t187 = t189 ^ 2;
t225 = t187 * t246;
t171 = t196 * t178;
t219 = t189 * t229;
t224 = t196 * t219 + t228 * t235 + t171;
t223 = t189 * t239;
t150 = -t243 - t264;
t222 = -t150 - t260;
t221 = t241 * t188;
t220 = qJD(4) + t255;
t163 = -t228 - t263;
t218 = t150 * t193 + t163 * t234 - t245;
t217 = t261 - t264;
t195 = sin(qJ(1));
t198 = cos(qJ(1));
t215 = g(1) * t195 - g(2) * t198;
t213 = pkin(3) * t193 - qJ(4) * t196;
t212 = 0.2e1 * (t188 * t246 - t242 * t236) * MDP(8) + (t188 * t190 + 0.2e1 * t234 * t247) * MDP(7) + (qJDD(3) * t196 - t193 * t199) * MDP(10) + (qJDD(3) * t193 + t196 * t199) * MDP(9) + t188 * MDP(4);
t149 = t220 - t258;
t152 = t237 + t254;
t211 = t149 * t193 + t152 * t196;
t210 = t178 + t243 - t260;
t209 = g(1) * t249 + g(2) * t250 - g(3) * t196 - t145;
t135 = (t213 * qJD(3) - t233) * t189 + t253 - t243;
t208 = -t135 - t253 - t261;
t158 = t164 - t265;
t207 = t158 * t189 - t226;
t206 = -qJDD(4) + t209;
t153 = pkin(3) * t235 - qJ(4) * t234 - t233;
t144 = t153 + t227;
t205 = -t144 * t189 - t158 * t188 - t135 - t251;
t204 = t149 * t234 - t152 * t235 - t244 + t267;
t181 = -pkin(2) - t265;
t203 = pkin(1) * t223 + t181 * t188 + t251;
t202 = -t230 + (t181 * t189 - t226) * qJD(3);
t201 = (t149 * t196 - t152 * t193) * qJD(3) + t267;
t200 = -pkin(6) * t179 - g(2) * (t185 * pkin(6) + qJ(4) * t249 + (pkin(2) + t262) * t186) - t164 * t178;
t156 = t163 * t235;
t154 = t213 * t189;
t142 = -t228 + t252;
t140 = t142 * t235;
t1 = [qJDD(1) * MDP(1) + t215 * MDP(2) + (g(1) * t198 + g(2) * t195) * MDP(3) + ((t188 * t197 - t223) * pkin(1) + t210) * MDP(5) + (((-qJDD(1) - t188) * t194 + (-qJD(1) - t189) * t238) * pkin(1) + t244) * MDP(6) + (t156 + t171 + t202 * t193 + (-t203 + t222) * t196) * MDP(12) + (t203 * t193 + t202 * t196 + t218) * MDP(13) + (t140 + t171 + (t207 * qJD(3) - t230) * t193 + (t205 - t260) * t196) * MDP(14) + (t180 * t221 + t226 * t268 + t204) * MDP(15) + ((t230 + (-t142 - t207) * qJD(3)) * t196 + t205 * t193 + t245) * MDP(16) + (t135 * t158 + t142 * t144 + (t211 * t238 + t215) * pkin(1) + t201 * t180 + t200) * MDP(17) + t212; (t210 + t219) * MDP(5) + ((-t232 + (-qJD(2) + t189) * t240) * pkin(1) + t244) * MDP(6) + (t156 + (-pkin(2) * t236 - t257) * t193 + (-t217 + t222) * t196 + t224) * MDP(12) + ((-t257 + (t228 - t263) * qJD(3)) * t196 + (t217 - t219) * t193 + t218) * MDP(13) + (t140 + (t164 * t236 - t257) * t193 + (-t153 * t189 + t208 - t260) * t196 + t224) * MDP(14) + (pkin(6) * t221 - t228 * t268 + t204) * MDP(15) + ((t257 + (-t142 - t228 - t252) * qJD(3)) * t196 + ((-t153 + t229) * t189 + t208) * t193 + t245) * MDP(16) + (t135 * t164 + t142 * t153 + (-t142 * t194 - t211 * t197) * t259 + t201 * pkin(6) + t200) * MDP(17) + t212; -MDP(7) * t225 + t242 * MDP(8) * t187 + MDP(9) * t248 + t196 * t188 * MDP(10) + qJDD(3) * MDP(11) + (-t163 * t247 + t209) * MDP(12) + (g(3) * t193 - t146 + (-t163 * t189 + t244) * t196) * MDP(13) + (0.2e1 * t256 + (-t142 * t193 + t154 * t196) * t189 + t206) * MDP(14) + (-t213 * t188 + ((t152 - t237) * t193 + (qJD(4) - t149 - t258) * t196) * t189) * MDP(15) + (0.2e1 * t231 + 0.2e1 * qJD(3) * qJD(4) + t146 + (t154 * t189 - g(3)) * t193 + (t142 * t189 - t244) * t196) * MDP(16) + (-t139 * pkin(3) - g(3) * t214 + t138 * qJ(4) - t142 * t154 - t149 * t254 + t220 * t152 + t244 * t213) * MDP(17); (-qJDD(3) - t225) * MDP(14) + MDP(15) * t248 + (-t187 * t190 - t199) * MDP(16) + (-qJD(3) * t152 + t142 * t247 - t206 + t266) * MDP(17);];
tau = t1;
