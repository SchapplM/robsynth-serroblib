% Calculate vector of inverse dynamics joint torques for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:58
% EndTime: 2019-12-31 18:00:00
% DurationCPUTime: 1.41s
% Computational Cost: add. (800->230), mult. (1434->315), div. (0->0), fcn. (862->10), ass. (0->105)
t184 = sin(pkin(8));
t169 = pkin(1) * t184 + qJ(3);
t240 = qJDD(1) * t169;
t190 = cos(qJ(4));
t222 = qJDD(1) * t190;
t261 = qJD(4) * qJD(5) + t222;
t179 = qJ(1) + pkin(8);
t174 = sin(t179);
t175 = cos(t179);
t209 = g(1) * t174 - g(2) * t175;
t185 = cos(pkin(8));
t171 = -pkin(1) * t185 - pkin(2);
t260 = qJDD(1) * t171;
t167 = -pkin(6) + t171;
t152 = t167 * qJDD(1) + qJDD(3);
t187 = sin(qJ(4));
t154 = t167 * qJD(1) + qJD(3);
t144 = qJD(2) * t190 + t154 * t187;
t235 = qJD(4) * t144;
t135 = -qJDD(4) * pkin(4) + qJDD(2) * t187 - t152 * t190 + t235;
t143 = -qJD(2) * t187 + t154 * t190;
t141 = -qJD(4) * pkin(4) - t143;
t210 = pkin(4) * t187 - pkin(7) * t190;
t151 = t169 + t210;
t226 = qJD(1) * qJD(4);
t217 = t190 * t226;
t223 = qJDD(1) * t187;
t156 = qJDD(5) + t217 + t223;
t168 = qJD(1) * t187 + qJD(5);
t145 = t151 * qJD(1);
t236 = qJD(4) * t143;
t213 = -qJDD(4) * pkin(7) - qJD(5) * t145 - qJDD(2) * t190 - t152 * t187 - t236;
t231 = qJD(4) * t190;
t259 = -(qJD(5) * t151 + t167 * t231) * t168 + t135 * t190 + (-qJD(4) * t141 - t156 * t167 + t213) * t187;
t164 = qJD(1) * t169;
t258 = 0.2e1 * qJD(4) * t164 + qJDD(4) * t167;
t211 = pkin(4) * t190 + pkin(7) * t187;
t257 = (pkin(7) * qJD(5) + t211 * qJD(1)) * t168 + t209 * t190 - g(3) * t187 + t135;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t227 = t189 * qJD(4);
t218 = t187 * t227;
t228 = qJD(5) * t190;
t201 = -t186 * t228 - t218;
t220 = t186 * qJDD(4) + t261 * t189;
t139 = t201 * qJD(1) + t220;
t254 = t139 * t186;
t233 = qJD(4) * t186;
t237 = qJD(1) * t190;
t160 = t189 * t237 + t233;
t232 = qJD(4) * t187;
t219 = t186 * t232;
t241 = -qJD(1) * t219 - t189 * qJDD(4);
t140 = qJD(5) * t160 + t186 * t222 + t241;
t253 = t140 * t187;
t252 = t151 * t156;
t251 = t156 * t189;
t158 = t186 * t237 - t227;
t250 = t158 * t168;
t249 = t160 * t168;
t248 = t160 * t189;
t247 = t167 * t187;
t246 = t186 * t187;
t245 = t187 * t189;
t244 = t190 * t139;
t243 = qJDD(2) - g(3);
t242 = t139 * t187 + t160 * t231;
t182 = t190 ^ 2;
t239 = t187 ^ 2 - t182;
t192 = qJD(4) ^ 2;
t193 = qJD(1) ^ 2;
t238 = -t192 - t193;
t234 = qJD(4) * t158;
t142 = qJD(4) * pkin(7) + t144;
t230 = qJD(5) * t142;
t229 = qJD(5) * t168;
t225 = qJD(3) * qJD(1);
t157 = t211 * qJD(4) + qJD(3);
t138 = t157 * qJD(1) + t210 * qJDD(1) + t240;
t212 = -t138 + t230;
t188 = sin(qJ(1));
t191 = cos(qJ(1));
t208 = g(1) * t188 - g(2) * t191;
t207 = qJDD(3) - t209;
t206 = -qJD(4) * t154 - t243;
t205 = g(3) * t190 + t213;
t204 = t156 * t186 + t189 * t229;
t203 = t186 * t229 - t251;
t200 = -g(1) * t175 - g(2) * t174 + t240;
t199 = (-t204 + t234) * MDP(20);
t198 = -pkin(7) * t156 + (t141 + t143) * t168;
t197 = qJD(1) * t164 + qJD(2) * qJD(4) - t152 + t209;
t155 = t225 + t240;
t195 = -t167 * t192 + t155 + t200 + t225;
t163 = qJDD(4) * t190 - t187 * t192;
t162 = -qJDD(4) * t187 - t190 * t192;
t153 = t168 * t219;
t149 = -t174 * t186 + t175 * t245;
t148 = t174 * t189 + t175 * t246;
t147 = t174 * t245 + t175 * t186;
t146 = -t174 * t246 + t175 * t189;
t136 = t189 * t138;
t133 = t142 * t189 + t145 * t186;
t132 = -t142 * t186 + t145 * t189;
t1 = [qJDD(1) * MDP(1) + t208 * MDP(2) + (g(1) * t191 + g(2) * t188) * MDP(3) + (t208 + (t184 ^ 2 + t185 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t207 + 0.2e1 * t260) * MDP(5) + (t200 + 0.2e1 * t225 + t240) * MDP(6) + (t155 * t169 + t164 * qJD(3) + (qJDD(3) + t260) * t171 - g(1) * (-pkin(1) * t188 - pkin(2) * t174 + qJ(3) * t175) - g(2) * (pkin(1) * t191 + pkin(2) * t175 + qJ(3) * t174)) * MDP(7) + (qJDD(1) * t182 - 0.2e1 * t187 * t217) * MDP(8) + 0.2e1 * (-t187 * t222 + t239 * t226) * MDP(9) + t163 * MDP(10) + t162 * MDP(11) + (t195 * t187 + t190 * t258) * MDP(13) + (-t187 * t258 + t195 * t190) * MDP(14) + (t201 * t160 + t189 * t244) * MDP(15) + ((t158 * t189 + t160 * t186) * t232 + (-t254 - t140 * t189 + (t158 * t186 - t248) * qJD(5)) * t190) * MDP(16) + (t201 * t168 + t190 * t251 + t242) * MDP(17) + (-t253 + t153 + (-t204 - t234) * t190) * MDP(18) + (t156 * t187 + t168 * t231) * MDP(19) + (-t190 * t167 * t140 - g(1) * t149 - g(2) * t147 + t136 * t187 + (t132 * t190 + t158 * t247) * qJD(4) + (t252 + t157 * t168 + (t141 * t190 + (-t167 * t168 - t142) * t187) * qJD(5)) * t189 + t259 * t186) * MDP(20) + (-t167 * t244 + g(1) * t148 - g(2) * t146 + (-t133 * t190 + t160 * t247) * qJD(4) + (-(-qJD(5) * t247 + t157) * t168 - t252 + t212 * t187 - t141 * t228) * t186 + t259 * t189) * MDP(21); t162 * MDP(13) - t163 * MDP(14) + (t153 + t253) * MDP(20) + (t168 * t218 + t242) * MDP(21) + (MDP(4) + MDP(7)) * t243 + (t203 * MDP(21) + t199) * t190; -t193 * MDP(6) + t207 * MDP(7) + (t171 * MDP(7) + MDP(5)) * qJDD(1) + (-t164 * MDP(7) + (-t189 * MDP(20) + t186 * MDP(21)) * t168) * qJD(1) + (qJDD(4) * MDP(13) + t238 * MDP(14) + (-t168 * t233 - t140) * MDP(20) + (-t168 * t227 - t139) * MDP(21)) * t190 + (t238 * MDP(13) - qJDD(4) * MDP(14) + t199 + (qJD(4) * t160 + t203) * MDP(21)) * t187; MDP(10) * t222 - MDP(11) * t223 + qJDD(4) * MDP(12) + (t206 * t187 - t197 * t190 + t235) * MDP(13) + (t197 * t187 + t206 * t190 + t236) * MDP(14) + (t168 * t248 + t254) * MDP(15) + ((t139 - t250) * t189 + (-t140 - t249) * t186) * MDP(16) + ((-t160 * t190 + t168 * t245) * qJD(1) + t204) * MDP(17) + ((t158 * t190 - t168 * t246) * qJD(1) - t203) * MDP(18) - t168 * MDP(19) * t237 + (-pkin(4) * t140 - t132 * t237 - t144 * t158 + t198 * t186 - t189 * t257) * MDP(20) + (-pkin(4) * t139 + t133 * t237 - t144 * t160 + t186 * t257 + t198 * t189) * MDP(21) + (t190 * t187 * MDP(8) - t239 * MDP(9)) * t193; t160 * t158 * MDP(15) + (-t158 ^ 2 + t160 ^ 2) * MDP(16) + (t220 + t250) * MDP(17) + (-t241 + t249) * MDP(18) + t156 * MDP(19) + (-g(1) * t146 - g(2) * t148 + t133 * t168 - t141 * t160 + t136) * MDP(20) + (g(1) * t147 - g(2) * t149 + t132 * t168 + t141 * t158) * MDP(21) + (-MDP(20) * t230 + t205 * MDP(21) + (-MDP(17) * t232 - MDP(18) * t228) * qJD(1)) * t189 + (-qJD(1) * MDP(17) * t228 - t261 * MDP(18) + t205 * MDP(20) + t212 * MDP(21)) * t186;];
tau = t1;
