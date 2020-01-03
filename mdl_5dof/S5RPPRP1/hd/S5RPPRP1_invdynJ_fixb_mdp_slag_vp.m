% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:57
% EndTime: 2020-01-03 11:26:01
% DurationCPUTime: 1.51s
% Computational Cost: add. (891->222), mult. (1740->297), div. (0->0), fcn. (1071->10), ass. (0->113)
t205 = cos(qJ(4));
t238 = qJD(1) * qJD(4);
t226 = t205 * t238;
t203 = sin(qJ(4));
t233 = qJDD(1) * t203;
t279 = t226 + t233;
t201 = cos(pkin(7));
t185 = -pkin(1) * t201 - pkin(2);
t198 = sin(pkin(8));
t200 = cos(pkin(8));
t167 = -pkin(3) * t200 - pkin(6) * t198 + t185;
t157 = t167 * qJD(1) + qJD(3);
t199 = sin(pkin(7));
t180 = pkin(1) * t199 + qJ(3);
t172 = t180 * qJD(1);
t164 = qJD(2) * t198 + t172 * t200;
t250 = qJD(1) * t198;
t229 = qJ(5) * t250;
t144 = t164 * t205 + (t157 - t229) * t203;
t278 = qJD(4) * t144;
t195 = qJ(1) + pkin(7);
t189 = sin(t195);
t190 = cos(t195);
t276 = g(2) * t190 + g(3) * t189;
t239 = qJD(1) * qJD(3);
t168 = qJDD(1) * t180 + t239;
t261 = t200 * t203;
t158 = -t189 * t261 - t190 * t205;
t160 = -t189 * t205 + t190 * t261;
t275 = -g(2) * t158 - g(3) * t160;
t271 = qJ(5) * t198;
t187 = t200 * qJDD(2);
t268 = t168 * t198;
t155 = -t187 + t268;
t270 = t155 * t198;
t267 = (pkin(4) * t205 + pkin(3)) * t200;
t266 = t189 * t203;
t265 = t190 * t203;
t193 = t198 ^ 2;
t207 = qJD(1) ^ 2;
t264 = t193 * t207;
t263 = t198 * t205;
t262 = t200 * MDP(5);
t260 = t200 * t205;
t225 = t205 * t157 - t164 * t203;
t143 = -t205 * t229 + t225;
t249 = qJD(1) * t200;
t178 = -qJD(4) + t249;
t140 = -pkin(4) * t178 + t143;
t259 = -t143 + t140;
t243 = qJD(4) * t205;
t246 = qJD(3) * t205;
t258 = t167 * t243 + t200 * t246;
t169 = t180 * t260;
t257 = t203 * t167 + t169;
t255 = t200 ^ 2 + t193;
t196 = t203 ^ 2;
t197 = t205 ^ 2;
t254 = t196 - t197;
t253 = MDP(10) * t193;
t252 = MDP(12) * t198;
t251 = MDP(16) * t198;
t248 = qJD(1) * t203;
t247 = qJD(3) * t200;
t245 = qJD(4) * t178;
t244 = qJD(4) * t203;
t242 = qJD(5) * t198;
t234 = qJDD(1) * t200;
t177 = -qJDD(4) + t234;
t241 = t177 * MDP(13);
t240 = -qJD(4) - t178;
t237 = qJD(1) * qJD(5);
t235 = qJDD(1) * t185;
t232 = qJDD(1) * t205;
t231 = qJ(5) * t263;
t206 = cos(qJ(1));
t230 = t206 * pkin(1) + t190 * pkin(2) + t189 * qJ(3);
t228 = t200 * t244;
t188 = t200 * qJD(2);
t163 = t172 * t198 - t188;
t227 = t163 * t250;
t224 = (-t196 - t197) * MDP(16);
t154 = t167 * qJDD(1) + qJDD(3);
t151 = t205 * t154;
t156 = qJDD(2) * t198 + t168 * t200;
t138 = -pkin(4) * t177 - t203 * t156 + t151 + (-qJ(5) * qJDD(1) - t237) * t263 - t278;
t223 = -t138 - t278;
t211 = -t203 * t154 - t205 * t156 - t157 * t243 + t164 * t244;
t139 = (-t279 * qJ(5) - t203 * t237) * t198 - t211;
t222 = qJD(4) * t140 - t139;
t221 = qJD(1) * t240;
t220 = t177 - t234;
t219 = t177 + t234;
t218 = t279 * pkin(4) * t198 + qJDD(5) - t187;
t216 = -g(2) * t189 + g(3) * t190;
t204 = sin(qJ(1));
t215 = -g(2) * t206 - g(3) * t204;
t214 = t204 * pkin(1) + t189 * pkin(2) - qJ(3) * t190;
t213 = t156 * t200 + t270;
t212 = t163 * t198 + t164 * t200;
t209 = -t178 * t249 + t245 - t264;
t171 = qJDD(3) + t235;
t170 = t228 * t250;
t166 = t205 * t167;
t161 = t190 * t260 + t266;
t159 = t189 * t260 - t265;
t148 = qJD(5) - t188 + (pkin(4) * t248 + t172) * t198;
t147 = -t203 * t271 + t257;
t146 = t218 + t268;
t145 = -t231 + t166 + (-t180 * t203 - pkin(4)) * t200;
t142 = -t203 * t247 - t205 * t242 + (-t169 + (-t167 + t271) * t203) * qJD(4);
t141 = -t203 * t242 + (-t180 * t261 - t231) * qJD(4) + t258;
t1 = [qJDD(1) * MDP(1) + t215 * MDP(2) + (g(2) * t204 - g(3) * t206) * MDP(3) + (t215 + (t199 ^ 2 + t201 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t168 * t255 + t213 + t216) * MDP(7) + (-g(2) * t230 - g(3) * t214 + t212 * qJD(3) + t171 * t185 + t213 * t180) * MDP(8) + (qJDD(1) * t197 - 0.2e1 * t203 * t226) * t193 * MDP(9) + 0.2e1 * (-t203 * t232 + t254 * t238) * t253 + (t170 + (t178 * t244 - t219 * t205) * t198) * MDP(11) + (t219 * t203 + (t178 + t249) * t243) * t252 + t200 * t241 + (-g(2) * t161 - g(3) * t159 - t151 * t200 - t166 * t177 + ((qJD(1) * t193 + t178 * t200) * t180 + t212) * t243 + (-(-qJD(4) * t167 - t247) * t178 - (-qJD(4) * t157 - t156) * t200 + t193 * t239 + t270 + (qJDD(1) * t193 + t177 * t200) * t180) * t203) * MDP(14) + ((-t180 * t228 + t258) * t178 + t257 * t177 - t211 * t200 + g(2) * t160 - g(3) * t158 + (t155 * t205 - t163 * t244) * t198 + (t180 * t232 + (-t180 * t244 + t246) * qJD(1)) * t193) * MDP(15) + ((-qJDD(1) * t145 + (-qJD(4) * t147 - t142) * qJD(1) + t223) * t205 + (-qJDD(1) * t147 + (qJD(4) * t145 - t141) * qJD(1) + t222) * t203 - t276) * t251 + (t139 * t147 + t144 * t141 + t138 * t145 + t140 * t142 - g(2) * (pkin(4) * t266 + t190 * t267 + t230) - g(3) * (-pkin(4) * t265 + t189 * t267 + t214) + (t146 * (pkin(4) * t203 + t180) + t148 * (pkin(4) * t243 + qJD(3)) + t276 * (-qJ(5) - pkin(6))) * t198) * MDP(17) + (t198 * MDP(6) - t262) * (t171 + t276 + t235); (qJDD(2) - g(1)) * MDP(4) + (-t155 * t200 - g(1)) * MDP(8) + t170 * MDP(15) + (-t146 * t200 - g(1)) * MDP(17) + (t156 * MDP(8) + (t220 * MDP(14) - MDP(15) * t245 + t223 * MDP(17)) * t203 + (t220 * MDP(15) + t139 * MDP(17) + ((t178 - t249) * MDP(14) - t140 * MDP(17)) * qJD(4)) * t205) * t198; (-t164 * t249 + qJDD(3) - t227 + t276) * MDP(8) + (-t148 * t250 + t276) * MDP(17) - t255 * MDP(7) * t207 + (-t177 * MDP(14) + t209 * MDP(15) + (-t144 * t249 - t223) * MDP(17)) * t205 + (t209 * MDP(14) + t177 * MDP(15) + (t140 * t249 - t222) * MDP(17)) * t203 + (-t262 + t185 * MDP(8) + (MDP(6) + t224) * t198) * qJDD(1); t205 * t203 * MDP(9) * t264 - t254 * t207 * t253 + (t203 * t221 + t232) * t198 * MDP(11) + (t205 * t221 - t233) * t252 - t241 + (t151 + (t240 * t164 - t227) * t205 + (g(1) * t198 + t240 * t157 - t156) * t203 + t275) * MDP(14) + (-t225 * t178 + g(2) * t159 - g(3) * t161 + (g(1) * t205 + t163 * t248) * t198 + t211) * MDP(15) + (-pkin(4) * t232 + (pkin(4) * qJD(4) - t259) * t248) * t251 + (t259 * t144 + (t138 + (-qJD(1) * t148 * t205 + g(1) * t203) * t198 + t275) * pkin(4)) * MDP(17); t224 * t264 + (g(1) * t200 + t218 + (t168 + (t140 * t205 + t144 * t203) * qJD(1) + t216) * t198) * MDP(17);];
tau = t1;
