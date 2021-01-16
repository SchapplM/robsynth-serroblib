% Calculate vector of inverse dynamics joint torques for
% S4RRRP2
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
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:31
% EndTime: 2021-01-15 11:04:33
% DurationCPUTime: 0.97s
% Computational Cost: add. (824->215), mult. (1213->269), div. (0->0), fcn. (597->8), ass. (0->115)
t203 = cos(qJ(3));
t184 = pkin(3) * t203 + pkin(2);
t195 = qJD(1) + qJD(2);
t271 = t184 * t195;
t200 = sin(qJ(3));
t201 = sin(qJ(2));
t265 = pkin(1) * qJD(1);
t240 = t201 * t265;
t164 = pkin(6) * t195 + t240;
t225 = qJ(4) * t195 + t164;
t147 = t225 * t200;
t144 = qJD(3) * pkin(3) - t147;
t198 = qJ(1) + qJ(2);
t189 = sin(t198);
t180 = g(2) * t189;
t190 = cos(t198);
t182 = g(1) * t190;
t252 = t180 + t182;
t244 = qJD(3) * t200;
t230 = t195 * t244;
t194 = qJDD(1) + qJDD(2);
t260 = t184 * t194;
t270 = -pkin(3) * t230 + t260;
t215 = qJD(3) * t225;
t204 = cos(qJ(2));
t269 = pkin(1) * t204;
t196 = t200 ^ 2;
t268 = pkin(3) * t196;
t181 = g(1) * t189;
t267 = g(2) * t190;
t266 = t194 * pkin(2);
t199 = -qJ(4) - pkin(6);
t263 = qJDD(3) * pkin(3);
t249 = qJD(1) * t204;
t239 = pkin(1) * t249;
t151 = qJD(4) - t239 - t271;
t262 = t151 * t195;
t167 = -t184 - t269;
t261 = t167 * t195;
t259 = t189 * t203;
t258 = t190 * t200;
t193 = t195 ^ 2;
t257 = t193 * t200;
t256 = t194 * t203;
t183 = pkin(1) * t201 + pkin(6);
t255 = -qJ(4) - t183;
t254 = g(1) * t258 + t200 * t180;
t248 = qJD(2) * t201;
t238 = pkin(1) * t248;
t253 = qJD(1) * t238 - qJDD(1) * t269;
t197 = t203 ^ 2;
t251 = -t196 - t197;
t250 = t196 - t197;
t247 = qJD(2) * t204;
t148 = t225 * t203;
t246 = qJD(3) * t148;
t245 = qJD(3) * t195;
t243 = qJD(3) * t203;
t242 = qJDD(1) * t201;
t241 = qJDD(3) * t183;
t237 = pkin(1) * t247;
t210 = qJD(1) * t247 + t242;
t153 = t210 * pkin(1) + pkin(6) * t194;
t226 = -qJ(4) * t194 - t153;
t211 = qJD(4) * t195 - t226;
t138 = -t200 * t215 + t211 * t203;
t236 = t138 * t203 - t252;
t140 = qJDD(4) + t253 - t270;
t175 = g(2) * t258;
t235 = t140 * t200 + t151 * t243 + t175;
t152 = t253 - t266;
t165 = -pkin(2) * t195 - t239;
t234 = t152 * t200 + t165 * t243 + t175;
t176 = g(1) * t259;
t221 = qJD(3) * t239;
t222 = t195 * t240;
t233 = t200 * t221 + t203 * t222 + t176;
t232 = g(2) * t259 + g(3) * t200 + t203 * t182;
t231 = t195 * t248;
t145 = t151 * t244;
t229 = -t140 - t267;
t228 = -t152 - t267;
t227 = qJD(3) * t199;
t224 = qJD(3) * t255;
t223 = 0.2e1 * t195 * t243;
t206 = qJD(3) ^ 2;
t220 = -pkin(6) * t206 + t266;
t202 = sin(qJ(1));
t205 = cos(qJ(1));
t219 = g(1) * t202 - g(2) * t205;
t218 = 0.2e1 * (t200 * t256 - t250 * t245) * MDP(8) + (t194 * t196 + t200 * t223) * MDP(7) + (qJDD(3) * t203 - t200 * t206) * MDP(10) + (qJDD(3) * t200 + t203 * t206) * MDP(9) + t194 * MDP(4);
t217 = t144 * t200 - t148 * t203;
t162 = pkin(3) * t244 + t238;
t216 = t162 * t195 + t167 * t194;
t214 = -t181 + t253 + t267;
t185 = -pkin(2) - t269;
t213 = t185 * t195 - t237;
t212 = -pkin(2) * t245 - pkin(6) * qJDD(3);
t209 = -t222 - t181;
t208 = pkin(1) * t231 + t183 * t206 + t185 * t194;
t207 = -g(1) * (-t184 * t189 - t190 * t199) - g(2) * (t190 * t184 - t189 * t199);
t191 = t203 * qJ(4);
t188 = t203 * qJD(4);
t171 = pkin(6) * t203 + t191;
t170 = t199 * t200;
t169 = t203 * t221;
t159 = t183 * t203 + t191;
t158 = t255 * t200;
t156 = t165 * t244;
t155 = -qJD(4) * t200 + t203 * t227;
t154 = t200 * t227 + t188;
t143 = (-qJD(4) - t237) * t200 + t203 * t224;
t142 = t200 * t224 + t203 * t237 + t188;
t137 = -t211 * t200 - t203 * t215 + t263;
t1 = [(t156 + t176) * MDP(12) + (-qJD(3) * t142 - qJDD(3) * t159 + t235) * MDP(15) + t252 * MDP(6) + t234 * MDP(13) + (qJD(3) * t143 + qJDD(3) * t158 + t145 + t176) * MDP(14) + t218 + ((-t208 + t228) * MDP(12) - MDP(13) * t241 + (-t216 + t229) * MDP(14) + (t142 * t195 + t159 * t194) * MDP(16) + (t213 * MDP(13) + MDP(15) * t261 + (-t158 * t195 - t144) * MDP(16)) * qJD(3)) * t203 + (-MDP(12) * t241 + (t208 - t181) * MDP(13) + (t216 - t181) * MDP(15) + (-t143 * t195 - t158 * t194 - t137) * MDP(16) + (t213 * MDP(12) + MDP(14) * t261 + (-t159 * t195 - t148) * MDP(16)) * qJD(3)) * t200 - t214 * MDP(5) + ((t194 * t204 - t231) * MDP(5) + (-t194 * t201 - t195 * t247 - t210) * MDP(6) + t219 * MDP(17)) * pkin(1) + (t137 * t158 + t138 * t159 + t140 * t167 + t148 * t142 + t144 * t143 + t151 * t162 + t207) * MDP(17) + t236 * MDP(16) + qJDD(1) * MDP(1) + t219 * MDP(2) + (g(1) * t205 + g(2) * t202) * MDP(3); (-t214 + t222) * MDP(5) + ((-t242 + (-qJD(2) + t195) * t249) * pkin(1) + t252) * MDP(6) + (t156 + t212 * t200 + (t220 + t228) * t203 + t233) * MDP(12) + (t169 + t212 * t203 + (t209 - t220) * t200 + t234) * MDP(13) + (qJDD(3) * t170 + t145 + (-t200 * t271 + t155) * qJD(3) + (t229 + t270) * t203 + t233) * MDP(14) + (-qJDD(3) * t171 + t169 + (t209 - t260) * t200 + (-t154 + (-t184 * t203 + t268) * t195) * qJD(3) + t235) * MDP(15) + ((-qJD(3) * t144 + t171 * t194) * t203 + (-t170 * t194 - t137 - t246) * t200 + (t154 * t203 - t155 * t200 + (-t170 * t203 - t171 * t200) * qJD(3) + t251 * t239) * t195 + t236) * MDP(16) + (t138 * t171 + t148 * t154 + t137 * t170 + t144 * t155 - t140 * t184 + pkin(3) * t145 + (-t151 * t201 + t217 * t204) * t265 + t207) * MDP(17) + t218; qJDD(3) * MDP(11) + t254 * MDP(12) + t232 * MDP(13) + (t246 + t254 + 0.2e1 * t263) * MDP(14) + (-qJD(3) * t147 + t232) * MDP(15) + (pkin(3) * t137 + (t144 + t147) * t148) * MDP(17) + (-MDP(15) * t268 + t250 * MDP(8)) * t193 + (t194 * MDP(9) + (-t165 * t195 - t153) * MDP(12) + (-t211 - t262) * MDP(14) + MDP(15) * t215 + (-t194 * MDP(16) + (-t262 + t252) * MDP(17)) * pkin(3)) * t200 + (-MDP(7) * t257 + t194 * MDP(10) - t153 * MDP(13) + (pkin(3) * t257 - qJD(3) * t164) * MDP(14) + t226 * MDP(15) + (-MDP(17) * pkin(3) - MDP(12) - MDP(14)) * g(3) + (-t165 * MDP(13) - qJ(4) * qJD(3) * MDP(14) + (-qJD(4) - t151) * MDP(15)) * t195) * t203; (0.2e1 * t230 - t256) * MDP(14) + (t194 * t200 + t223) * MDP(15) + (t217 * t195 - t181 - t229) * MDP(17) + t251 * MDP(16) * t193;];
tau = t1;
