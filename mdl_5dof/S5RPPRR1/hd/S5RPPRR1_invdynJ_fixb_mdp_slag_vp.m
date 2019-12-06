% Calculate vector of inverse dynamics joint torques for
% S5RPPRR1
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
%   see S5RPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:21
% EndTime: 2019-12-05 17:38:25
% DurationCPUTime: 1.28s
% Computational Cost: add. (612->207), mult. (1069->258), div. (0->0), fcn. (625->8), ass. (0->99)
t204 = sin(qJ(1));
t207 = cos(qJ(1));
t225 = g(1) * t207 + g(2) * t204;
t262 = g(1) * t204 - g(2) * t207;
t201 = pkin(1) + qJ(3);
t261 = qJD(1) * t201;
t191 = qJD(4) + qJD(5);
t168 = -qJD(2) + t261;
t200 = -pkin(6) + qJ(2);
t260 = (qJD(2) + t168 + t261) * qJD(4) + qJDD(4) * t200;
t194 = (qJD(1) * qJD(2));
t186 = 2 * t194;
t257 = pkin(7) - t200;
t179 = qJ(2) * qJD(1) + qJD(3);
t167 = -pkin(6) * qJD(1) + t179;
t203 = sin(qJ(4));
t246 = qJD(1) * t203;
t148 = -pkin(7) * t246 + t167 * t203;
t205 = cos(qJ(5));
t256 = t148 * t205;
t202 = sin(qJ(5));
t206 = cos(qJ(4));
t245 = qJD(1) * t206;
t152 = t202 * t246 - t205 * t245;
t255 = t152 * t191;
t156 = t202 * t206 + t203 * t205;
t153 = t156 * qJD(1);
t254 = t153 * t191;
t141 = t191 * t156;
t157 = -t202 * t203 + t205 * t206;
t190 = qJDD(4) + qJDD(5);
t253 = -t141 * t191 + t157 * t190;
t241 = qJD(1) * qJD(4);
t232 = t203 * t241;
t242 = qJD(5) * t202;
t233 = t203 * t242;
t252 = qJD(1) * t233 + t202 * t232;
t251 = t207 * pkin(1) + t204 * qJ(2);
t197 = t206 ^ 2;
t249 = t203 ^ 2 - t197;
t208 = qJD(4) ^ 2;
t209 = qJD(1) ^ 2;
t248 = -t208 - t209;
t247 = MDP(22) * t205;
t244 = qJD(4) * t203;
t243 = qJD(4) * t206;
t240 = qJD(3) * qJD(1);
t193 = qJDD(1) * qJ(2);
t239 = qJDD(1) * t201;
t238 = qJDD(1) * t203;
t237 = qJDD(1) * t206;
t235 = qJDD(4) * t203;
t199 = qJDD(1) * pkin(1);
t234 = qJDD(2) - t199;
t162 = t257 * t206;
t231 = t206 * t241;
t230 = qJDD(2) - t262;
t229 = qJDD(3) + t193 + t194;
t228 = t191 * t206;
t227 = MDP(23) * t191;
t174 = pkin(4) * t203 + t201;
t226 = -t199 + t230;
t149 = -pkin(7) * t245 + t206 * t167;
t142 = -t202 * t244 + t205 * t228 - t233;
t222 = -t142 * t191 - t156 * t190;
t145 = qJD(4) * pkin(4) + t149;
t221 = -t145 * t202 - t256;
t161 = t257 * t203;
t220 = -t161 * t205 - t162 * t202;
t219 = t161 * t202 - t162 * t205;
t192 = qJDD(1) * qJ(3);
t218 = -t192 + t226;
t217 = t202 * t237 - t252;
t159 = t192 - t234 + t240;
t216 = 0.2e1 * t193 + t186 - t225;
t215 = t231 + t238;
t214 = qJD(1) * t168 + t225;
t170 = t205 * t237;
t134 = -t141 * qJD(1) - t202 * t238 + t170;
t213 = -t152 * t153 * MDP(17) + (t134 + t254) * MDP(19) + (-t255 + (-t191 * t245 - t238) * t205 - t217) * MDP(20) + (t152 ^ 2 - t153 ^ 2) * MDP(18) + t190 * MDP(21);
t158 = -pkin(6) * qJDD(1) + t229;
t151 = t206 * t158;
t138 = -t167 * t244 + qJDD(4) * pkin(4) + t151 + (t232 - t237) * pkin(7);
t155 = t174 * qJD(1) - qJD(2);
t198 = qJ(4) + qJ(5);
t181 = sin(t198);
t182 = cos(t198);
t212 = t148 * t242 + g(3) * t182 + (-t148 * t191 - t138) * t202 + t155 * t153 + t225 * t181;
t211 = -t200 * t208 + t159 + t239 + t240 + t262;
t139 = -t215 * pkin(7) + t158 * t203 + t167 * t243;
t210 = g(3) * t181 + t221 * qJD(5) + t205 * t138 - t202 * t139 + t155 * t152 - t225 * t182;
t184 = t207 * qJ(2);
t180 = qJDD(4) * t206;
t169 = pkin(4) * t243 + qJD(3);
t147 = qJD(2) * t203 - qJD(4) * t162;
t146 = qJD(2) * t206 + t257 * t244;
t143 = t215 * pkin(4) + t159;
t135 = (qJD(1) * t228 + t238) * t205 + t217;
t1 = [qJDD(1) * MDP(1) + t262 * MDP(2) + t225 * MDP(3) + (-0.2e1 * t199 + t230) * MDP(4) + t216 * MDP(5) + (-t234 * pkin(1) - g(1) * (-pkin(1) * t204 + t184) - g(2) * t251 + (t193 + t186) * qJ(2)) * MDP(6) + (qJDD(3) + t216) * MDP(7) + (-t218 + t239 + 0.2e1 * t240) * MDP(8) + (t159 * t201 + t168 * qJD(3) + t229 * qJ(2) + t179 * qJD(2) - g(1) * (-t201 * t204 + t184) - g(2) * (qJ(3) * t207 + t251)) * MDP(9) + (qJDD(1) * t197 - 0.2e1 * t203 * t231) * MDP(10) + 0.2e1 * (-t203 * t237 + t249 * t241) * MDP(11) + (-t203 * t208 + t180) * MDP(12) + (-t206 * t208 - t235) * MDP(13) + (t211 * t203 + t260 * t206) * MDP(15) + (-t260 * t203 + t211 * t206) * MDP(16) + (t134 * t157 + t141 * t152) * MDP(17) + (-t134 * t156 - t135 * t157 + t141 * t153 + t142 * t152) * MDP(18) + t253 * MDP(19) + t222 * MDP(20) + (t169 * t153 + t174 * t135 + t143 * t156 + t155 * t142 + (-t220 * qJD(5) + t146 * t205 - t147 * t202) * t191 + t219 * t190 + t262 * t181) * MDP(22) + (-t169 * t152 + t174 * t134 + t143 * t157 - t155 * t141 - (t219 * qJD(5) + t146 * t202 + t147 * t205) * t191 - t220 * t190 + t262 * t182) * MDP(23); t226 * MDP(6) + t218 * MDP(9) + (t252 + t255) * MDP(22) + (-t170 + t254) * MDP(23) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t209 + (MDP(4) - MDP(8) + (-MDP(22) * t202 - MDP(16)) * t206 + (MDP(23) * t202 - MDP(15) - t247) * t203) * qJDD(1) + ((-qJD(3) - t179) * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + t205 * t227) * t203 + (-0.2e1 * qJD(4) * MDP(15) - t191 * t247 + t202 * t227) * t206) * qJD(1); qJDD(1) * MDP(7) - t209 * MDP(8) + (-t214 + t229) * MDP(9) + (t248 * t203 + t180) * MDP(15) + (t248 * t206 - t235) * MDP(16) + (-qJD(1) * t153 + t253) * MDP(22) + (qJD(1) * t152 + t222) * MDP(23); MDP(12) * t237 - MDP(13) * t238 + qJDD(4) * MDP(14) + (g(3) * t203 - t214 * t206 + t151) * MDP(15) + (g(3) * t206 + (-t158 + t214) * t203) * MDP(16) + (-(-t149 * t202 - t256) * t191 + (-t153 * t245 + t190 * t205 - t191 * t242) * pkin(4) + t210) * MDP(22) + ((-qJD(5) * t145 + t149 * t191 - t139) * t205 + (-qJD(5) * t191 * t205 + t152 * t245 - t190 * t202) * pkin(4) + t212) * MDP(23) + t213 + (t206 * t203 * MDP(10) - t249 * MDP(11)) * t209; (-t221 * t191 + t210) * MDP(22) + ((-t139 + (-qJD(5) + t191) * t145) * t205 + t212) * MDP(23) + t213;];
tau = t1;
