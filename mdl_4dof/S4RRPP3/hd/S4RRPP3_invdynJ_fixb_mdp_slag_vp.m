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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
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
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:54
% EndTime: 2019-12-31 16:57:56
% DurationCPUTime: 1.17s
% Computational Cost: add. (835->209), mult. (1910->265), div. (0->0), fcn. (1185->8), ass. (0->92)
t190 = sin(qJ(1));
t192 = cos(qJ(1));
t237 = g(1) * t190 - g(2) * t192;
t207 = g(1) * t192 + g(2) * t190;
t186 = sin(pkin(6));
t187 = cos(pkin(6));
t191 = cos(qJ(2));
t223 = qJD(1) * t191;
t214 = t187 * t223;
t189 = sin(qJ(2));
t224 = qJD(1) * t189;
t153 = t186 * t224 - t214;
t163 = t186 * t191 + t187 * t189;
t156 = t163 * qJD(1);
t235 = pkin(2) * t191;
t178 = pkin(1) + t235;
t167 = -qJD(1) * t178 + qJD(3);
t127 = pkin(3) * t153 - qJ(4) * t156 + t167;
t182 = qJ(2) + pkin(6);
t179 = sin(t182);
t180 = cos(t182);
t230 = qJ(3) + pkin(5);
t208 = qJD(2) * t230;
t201 = -qJD(3) * t189 - t191 * t208;
t213 = t230 * t189;
t137 = qJDD(2) * pkin(2) + qJD(1) * t201 - qJDD(1) * t213;
t150 = qJD(3) * t191 - t189 * t208;
t168 = t230 * t191;
t143 = qJD(1) * t150 + qJDD(1) * t168;
t122 = t187 * t137 - t186 * t143;
t210 = -qJDD(4) + t122;
t236 = -g(3) * t180 - t127 * t156 + t207 * t179 + t210;
t152 = t156 ^ 2;
t231 = g(3) * t191;
t229 = qJDD(2) * pkin(3);
t166 = qJD(1) * t168;
t228 = t166 * t186;
t227 = t186 * t189;
t159 = t187 * t166;
t226 = t187 * t191;
t123 = t186 * t137 + t187 * t143;
t165 = qJD(1) * t213;
t161 = qJD(2) * pkin(2) - t165;
t140 = t186 * t161 + t159;
t184 = t189 ^ 2;
t225 = -t191 ^ 2 + t184;
t222 = qJD(2) * t189;
t145 = -t165 * t187 - t228;
t221 = qJD(4) - t145;
t220 = qJD(1) * qJD(2);
t219 = qJDD(1) * t189;
t218 = qJDD(1) * t191;
t212 = t189 * t220;
t217 = pkin(2) * t212 + qJDD(3);
t216 = pkin(2) * t222;
t215 = qJDD(2) * qJ(4) + t123;
t211 = t191 * t220;
t205 = t178 * qJDD(1);
t204 = pkin(3) * t180 + qJ(4) * t179;
t139 = t161 * t187 - t228;
t203 = -0.2e1 * pkin(1) * t220 - pkin(5) * qJDD(2);
t155 = t163 * qJD(2);
t169 = t186 * t212;
t199 = t163 * qJDD(1) - t169;
t193 = qJD(2) ^ 2;
t198 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t193 + t237;
t194 = qJD(1) ^ 2;
t197 = pkin(1) * t194 - pkin(5) * qJDD(1) + t207;
t171 = t187 * t218;
t141 = qJD(1) * t155 + t186 * t219 - t171;
t142 = t187 * t211 + t199;
t196 = pkin(3) * t141 - qJ(4) * t142 - qJD(4) * t156 + t217;
t129 = t150 * t186 - t187 * t201;
t130 = t187 * t150 + t186 * t201;
t146 = t168 * t186 + t187 * t213;
t147 = t187 * t168 - t186 * t213;
t195 = t129 * t156 - t130 * t153 - t147 * t141 + t142 * t146 - t207;
t177 = -pkin(2) * t187 - pkin(3);
t175 = pkin(2) * t186 + qJ(4);
t170 = t192 * t178;
t162 = -t226 + t227;
t158 = qJD(2) * t226 - t186 * t222;
t144 = -t165 * t186 + t159;
t138 = pkin(3) * t162 - qJ(4) * t163 - t178;
t134 = qJD(2) * qJ(4) + t140;
t131 = -qJD(2) * pkin(3) + qJD(4) - t139;
t128 = pkin(2) * t224 + pkin(3) * t156 + qJ(4) * t153;
t126 = pkin(3) * t155 - qJ(4) * t158 - qJD(4) * t163 + t216;
t121 = -t210 - t229;
t120 = qJD(2) * qJD(4) + t215;
t119 = -t205 + t196;
t1 = [qJDD(1) * MDP(1) + t237 * MDP(2) + t207 * MDP(3) + (qJDD(1) * t184 + 0.2e1 * t189 * t211) * MDP(4) + 0.2e1 * (t189 * t218 - t220 * t225) * MDP(5) + (qJDD(2) * t189 + t191 * t193) * MDP(6) + (qJDD(2) * t191 - t189 * t193) * MDP(7) + (t189 * t203 + t191 * t198) * MDP(9) + (-t189 * t198 + t191 * t203) * MDP(10) + (-t122 * t163 - t123 * t162 - t139 * t158 - t140 * t155 + t195) * MDP(11) + (t123 * t147 + t140 * t130 - t122 * t146 - t139 * t129 - (-t205 + t217) * t178 + t167 * t216 - g(1) * (-t178 * t190 + t192 * t230) - g(2) * (t190 * t230 + t170)) * MDP(12) + (-qJD(2) * t129 - qJDD(2) * t146 + t119 * t162 + t126 * t153 + t127 * t155 + t138 * t141 + t180 * t237) * MDP(13) + (-t120 * t162 + t121 * t163 + t131 * t158 - t134 * t155 + t195) * MDP(14) + (qJD(2) * t130 + qJDD(2) * t147 - t119 * t163 - t126 * t156 - t127 * t158 - t138 * t142 + t179 * t237) * MDP(15) + (-g(2) * t170 + t119 * t138 + t120 * t147 + t121 * t146 + t127 * t126 + t131 * t129 + t134 * t130 + (-g(1) * t230 - g(2) * t204) * t192 + (-g(1) * (-t178 - t204) - g(2) * t230) * t190) * MDP(16); MDP(6) * t219 + MDP(7) * t218 + qJDD(2) * MDP(8) + (t189 * t197 - t231) * MDP(9) + (g(3) * t189 + t191 * t197) * MDP(10) + ((t140 - t144) * t156 + (-t139 + t145) * t153 + (-t141 * t186 - t142 * t187) * pkin(2)) * MDP(11) + (t139 * t144 - t140 * t145 + (-t231 + t122 * t187 + t123 * t186 + (-qJD(1) * t167 + t207) * t189) * pkin(2)) * MDP(12) + (qJD(2) * t144 - t128 * t153 + (pkin(3) - t177) * qJDD(2) + t236) * MDP(13) + (-t141 * t175 + t142 * t177 + (t134 - t144) * t156 + (t131 - t221) * t153) * MDP(14) + (-g(3) * t179 + qJDD(2) * t175 - t127 * t153 + t128 * t156 - t207 * t180 + (0.2e1 * qJD(4) - t145) * qJD(2) + t215) * MDP(15) + (t120 * t175 + t121 * t177 - t127 * t128 - t131 * t144 - g(3) * (t204 + t235) + t221 * t134 + t207 * (pkin(2) * t189 + pkin(3) * t179 - qJ(4) * t180)) * MDP(16) + (-MDP(4) * t189 * t191 + MDP(5) * t225) * t194; (t139 * t156 + t140 * t153 + t217 - t237) * MDP(12) - t171 * MDP(13) + t169 * MDP(15) + (-t131 * t156 + t134 * t153 + t196 - t237) * MDP(16) + (MDP(13) * t227 - MDP(15) * t163 + (-MDP(12) - MDP(16)) * t178) * qJDD(1) + ((t186 * t223 + t187 * t224 + t156) * MDP(13) + (t153 - t214) * MDP(15)) * qJD(2) + (MDP(11) + MDP(14)) * (-t153 ^ 2 - t152); (t153 * t156 - qJDD(2)) * MDP(13) + ((t153 + t214) * qJD(2) + t199) * MDP(14) + (-t152 - t193) * MDP(15) + (-qJD(2) * t134 - t229 - t236) * MDP(16);];
tau = t1;
