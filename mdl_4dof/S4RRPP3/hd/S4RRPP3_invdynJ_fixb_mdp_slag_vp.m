% Calculate vector of inverse dynamics joint torques for
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:08
% EndTime: 2021-01-15 10:36:12
% DurationCPUTime: 1.79s
% Computational Cost: add. (963->234), mult. (2194->282), div. (0->0), fcn. (1373->8), ass. (0->97)
t210 = sin(qJ(1));
t212 = cos(qJ(1));
t256 = g(1) * t210 - g(2) * t212;
t226 = g(1) * t212 + g(2) * t210;
t211 = cos(qJ(2));
t253 = pkin(2) * t211;
t198 = pkin(1) + t253;
t224 = t198 * qJDD(1);
t237 = (qJD(1) * qJD(2));
t259 = -2 * t237;
t209 = sin(qJ(2));
t204 = t209 ^ 2;
t258 = (-t211 ^ 2 + t204) * MDP(5);
t206 = sin(pkin(6));
t207 = cos(pkin(6));
t178 = t206 * t211 + t207 * t209;
t171 = t178 * qJD(1);
t255 = MDP(12) - MDP(17);
t208 = -qJ(3) - pkin(5);
t228 = qJD(2) * t208;
t222 = -qJD(3) * t209 + t211 * t228;
t231 = t208 * t209;
t149 = qJDD(2) * pkin(2) + t222 * qJD(1) + qJDD(1) * t231;
t165 = qJD(3) * t211 + t209 * t228;
t185 = t208 * t211;
t155 = t165 * qJD(1) - qJDD(1) * t185;
t134 = t206 * t149 + t207 * t155;
t202 = qJ(2) + pkin(6);
t199 = sin(t202);
t200 = cos(t202);
t254 = g(3) * t199 + t226 * t200 - t134;
t167 = t171 ^ 2;
t249 = g(3) * t200;
t248 = g(3) * t211;
t247 = MDP(4) * t211;
t243 = t207 * t211;
t232 = qJD(1) * t243;
t240 = qJD(1) * t209;
t168 = t206 * t240 - t232;
t183 = -t198 * qJD(1) + qJD(3);
t139 = pkin(3) * t168 - qJ(4) * t171 + t183;
t246 = t139 * t171;
t181 = qJD(1) * t185;
t245 = t181 * t206;
t244 = t206 * t209;
t174 = t207 * t181;
t242 = t208 * t210;
t133 = t207 * t149 - t206 * t155;
t180 = qJD(1) * t231;
t176 = qJD(2) * pkin(2) + t180;
t152 = t206 * t176 - t174;
t239 = qJD(2) * t209;
t157 = t180 * t207 + t245;
t238 = qJD(4) - t157;
t236 = qJDD(1) * t209;
t235 = qJDD(1) * t211;
t230 = t209 * t237;
t233 = pkin(2) * t230 + qJDD(3);
t151 = t176 * t207 + t245;
t142 = t207 * t165 + t206 * t222;
t161 = -t207 * t185 + t206 * t231;
t223 = -qJD(2) * t142 - qJDD(2) * t161;
t132 = -qJDD(2) * pkin(3) + qJDD(4) - t133;
t170 = t178 * qJD(2);
t221 = t178 * qJDD(1) - t206 * t230;
t156 = t180 * t206 - t174;
t220 = qJD(2) * t156 + t226 * t199 + t133 - t249;
t213 = qJD(2) ^ 2;
t219 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t213 + t256;
t214 = qJD(1) ^ 2;
t218 = pkin(1) * t214 - pkin(5) * qJDD(1) + t226;
t141 = t165 * t206 - t207 * t222;
t160 = -t185 * t206 - t207 * t231;
t217 = -qJD(2) * t141 - qJDD(2) * t160 + t256 * t200;
t190 = t207 * t235;
t153 = qJD(1) * t170 + t206 * t236 - t190;
t154 = qJD(2) * t232 + t221;
t216 = pkin(3) * t153 - qJ(4) * t154 - qJD(4) * t171 + t233;
t215 = t141 * t171 - t142 * t168 - t161 * t153 + t154 * t160 - t226;
t203 = qJDD(2) * qJ(4);
t197 = t208 * t212;
t196 = -pkin(2) * t207 - pkin(3);
t194 = pkin(2) * t206 + qJ(4);
t184 = -pkin(3) * t206 + qJ(4) * t207;
t182 = pkin(3) * t207 + qJ(4) * t206 + pkin(2);
t177 = -t243 + t244;
t173 = qJD(2) * t243 - t206 * t239;
t164 = -t224 + t233;
t158 = t182 * t211 + t184 * t209 + pkin(1);
t150 = pkin(3) * t177 - qJ(4) * t178 - t198;
t146 = qJD(2) * qJ(4) + t152;
t143 = -qJD(2) * pkin(3) + qJD(4) - t151;
t140 = pkin(2) * t240 + pkin(3) * t171 + qJ(4) * t168;
t137 = pkin(2) * t239 + pkin(3) * t170 - qJ(4) * t173 - qJD(4) * t178;
t131 = qJD(2) * qJD(4) + t134 + t203;
t130 = -t224 + t216;
t1 = [t226 * MDP(3) + (-t153 * t198 + t164 * t177 + t170 * t183 + t217) * MDP(11) + (-t154 * t198 + t164 * t178 + t173 * t183 + t223) * MDP(12) + (-t133 * t178 - t134 * t177 - t151 * t173 - t152 * t170 + t215) * MDP(13) + (t134 * t161 + t152 * t142 - t133 * t160 - t151 * t141 - t164 * t198 - g(1) * (-t198 * t210 - t197) - g(2) * (t198 * t212 - t242)) * MDP(14) + (t130 * t177 + t137 * t168 + t139 * t170 + t150 * t153 + t217) * MDP(15) + (-t131 * t177 + t132 * t178 + t143 * t173 - t146 * t170 + t215) * MDP(16) + (-t130 * t178 - t137 * t171 - t139 * t173 - t150 * t154 - t223) * MDP(17) + (t131 * t161 + t146 * t142 + t130 * t150 + t139 * t137 + t132 * t160 + t143 * t141 - g(1) * (-t158 * t210 - t197) - g(2) * (t158 * t212 - t242)) * MDP(18) + (t204 * MDP(4) + MDP(1)) * qJDD(1) + t258 * t259 + (t213 * MDP(6) + qJDD(2) * MDP(7) + t219 * MDP(9) + (pkin(1) * t259 - pkin(5) * qJDD(2)) * MDP(10)) * t211 + (0.2e1 * MDP(5) * t235 - t213 * MDP(7) - t219 * MDP(10) + (-pkin(5) * MDP(9) + MDP(6)) * qJDD(2) + (0.2e1 * (-MDP(9) * pkin(1) + t247) * qJD(1) + (MDP(11) * t168 + MDP(12) * t171 + MDP(14) * t183) * pkin(2)) * qJD(2)) * t209 + (-t255 * t199 + MDP(2)) * t256; MDP(6) * t236 + MDP(7) * t235 + qJDD(2) * MDP(8) + (t218 * t209 - t248) * MDP(9) + (g(3) * t209 + t218 * t211) * MDP(10) + (-t171 * t183 + (qJDD(2) * t207 - t168 * t240) * pkin(2) + t220) * MDP(11) + (qJD(2) * t157 + t168 * t183 + (-qJDD(2) * t206 - t171 * t240) * pkin(2) + t254) * MDP(12) + ((t152 - t156) * t171 + (-t151 + t157) * t168 + (-t153 * t206 - t154 * t207) * pkin(2)) * MDP(13) + (t151 * t156 - t152 * t157 + (-t248 + t133 * t207 + t134 * t206 + (-qJD(1) * t183 + t226) * t209) * pkin(2)) * MDP(14) + (-t246 - t140 * t168 - qJDD(4) + (pkin(3) - t196) * qJDD(2) + t220) * MDP(15) + (-t153 * t194 + t154 * t196 + (t146 - t156) * t171 + (t143 - t238) * t168) * MDP(16) + (qJDD(2) * t194 - t139 * t168 + t140 * t171 + t203 + (0.2e1 * qJD(4) - t157) * qJD(2) - t254) * MDP(17) + (t131 * t194 + t132 * t196 - t139 * t140 - t143 * t156 - g(3) * (pkin(3) * t200 + qJ(4) * t199 + t253) - t226 * (-t182 * t209 + t184 * t211) + t238 * t146) * MDP(18) + (-t209 * t247 + t258) * t214; (t151 * t171 + t152 * t168 + t233 - t256) * MDP(14) + (-t143 * t171 + t146 * t168 + t216 - t256) * MDP(18) + t255 * ((-t168 + t232) * qJD(2) + t221) + (-MDP(14) - MDP(18)) * t224 + (MDP(13) + MDP(16)) * (-t168 ^ 2 - t167) + (0.2e1 * qJD(2) * t171 + qJDD(1) * t244 - t190) * (MDP(11) + MDP(15)); (t168 * t171 - qJDD(2)) * MDP(15) + ((t168 + t232) * qJD(2) + t221) * MDP(16) + (-t167 - t213) * MDP(17) + (-qJD(2) * t146 - t226 * t178 + t132 + t246 + t249) * MDP(18);];
tau = t1;
