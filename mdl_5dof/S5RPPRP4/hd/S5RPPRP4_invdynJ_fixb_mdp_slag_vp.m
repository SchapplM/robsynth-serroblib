% Calculate vector of inverse dynamics joint torques for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:24
% EndTime: 2019-12-31 17:52:26
% DurationCPUTime: 0.86s
% Computational Cost: add. (698->191), mult. (1174->227), div. (0->0), fcn. (632->6), ass. (0->94)
t185 = sin(qJ(4));
t178 = t185 ^ 2;
t186 = cos(qJ(4));
t179 = t186 ^ 2;
t206 = (t178 + t179) * MDP(17);
t238 = MDP(8) - t206;
t212 = qJDD(1) * t186;
t237 = pkin(4) * t212 + qJDD(5);
t187 = -pkin(1) - pkin(2);
t236 = cos(qJ(1));
t235 = sin(qJ(1));
t234 = g(3) * t186;
t233 = qJD(4) * pkin(4);
t232 = pkin(1) * qJDD(1);
t152 = t187 * qJD(1) + qJD(2);
t181 = sin(pkin(7));
t182 = cos(pkin(7));
t221 = qJ(2) * qJD(1);
t140 = t181 * t152 + t182 * t221;
t136 = -qJD(1) * pkin(6) + t140;
t219 = qJD(1) * t186;
t129 = -qJ(5) * t219 + qJD(3) * t185 + t136 * t186;
t231 = t129 * t186;
t147 = t182 * qJ(2) + t181 * t187;
t144 = -pkin(6) + t147;
t230 = qJ(5) - t144;
t229 = qJDD(3) + g(3);
t220 = qJD(1) * t185;
t128 = qJ(5) * t220 + t186 * qJD(3) - t136 * t185;
t125 = t128 + t233;
t228 = t125 - t128;
t151 = t187 * qJDD(1) + qJDD(2);
t216 = qJ(2) * qJDD(1);
t227 = t181 * t151 + t182 * t216;
t226 = t236 * pkin(1) + t235 * qJ(2);
t225 = g(1) * t235 - g(2) * t236;
t224 = t178 - t179;
t188 = qJD(4) ^ 2;
t189 = qJD(1) ^ 2;
t222 = t188 + t189;
t218 = qJD(2) * t182;
t217 = MDP(18) * qJD(1);
t215 = qJD(1) * qJD(2);
t214 = qJD(1) * qJD(4);
t213 = qJDD(1) * t185;
t211 = 0.2e1 * t215;
t156 = t182 * t215;
t134 = t156 + t227;
t210 = t236 * pkin(2) + t226;
t209 = t185 * t214;
t154 = t181 * t215;
t208 = t151 * t182 - t181 * t216;
t139 = t152 * t182 - t181 * t221;
t146 = -t181 * qJ(2) + t182 * t187;
t207 = qJD(4) * t230;
t169 = t186 * qJDD(3);
t131 = -qJDD(1) * pkin(6) + t134;
t192 = qJ(5) * qJDD(1) + qJD(1) * qJD(5) - qJD(3) * qJD(4) - t131;
t197 = (qJ(5) * qJD(1) - t136) * qJD(4);
t122 = qJDD(4) * pkin(4) + t192 * t185 + t186 * t197 + t169;
t205 = -qJD(4) * t129 - t122;
t123 = (qJDD(3) + t197) * t185 - t192 * t186;
t204 = -qJD(4) * t125 + t123;
t203 = 0.2e1 * t186 * t214;
t202 = qJDD(2) - t232;
t201 = -qJD(5) + t218;
t143 = pkin(3) - t146;
t200 = -t235 * pkin(1) + t236 * qJ(2);
t141 = -t235 * t181 - t236 * t182;
t142 = t236 * t181 - t235 * t182;
t199 = g(1) * t142 - g(2) * t141;
t198 = -g(1) * t141 - g(2) * t142;
t133 = -t154 + t208;
t135 = qJD(1) * pkin(3) - t139;
t180 = qJDD(1) * pkin(3);
t130 = -t133 + t180;
t196 = g(1) * t236 + g(2) * t235;
t195 = -t235 * pkin(2) + t200;
t194 = -t199 - t208;
t193 = qJD(1) * t135 - t131 + t198;
t191 = -qJDD(4) * t144 + (-qJD(1) * t143 - t135 - t218) * qJD(4);
t190 = qJDD(1) * t143 - t144 * t188 + t130 + t154 - t199;
t184 = -qJ(5) - pkin(6);
t173 = t186 * pkin(4);
t165 = t173 + pkin(3);
t150 = qJDD(4) * t186 - t185 * t188;
t149 = -qJDD(4) * t185 - t186 * t188;
t138 = t230 * t186;
t137 = t230 * t185;
t132 = pkin(4) * t219 + qJD(5) + t135;
t127 = -t201 * t185 + t186 * t207;
t126 = t185 * t207 + t201 * t186;
t124 = -pkin(4) * t209 + t130 + t237;
t1 = [qJDD(1) * MDP(1) + t225 * MDP(2) + t196 * MDP(3) + (-qJDD(2) + t225 + 0.2e1 * t232) * MDP(4) + (-t196 + t211 + 0.2e1 * t216) * MDP(5) + (-t202 * pkin(1) - g(1) * t200 - g(2) * t226 + (t211 + t216) * qJ(2)) * MDP(6) + (-qJDD(1) * t146 + 0.2e1 * t154 + t194) * MDP(7) + (qJDD(1) * t147 + 0.2e1 * t156 - t198 + t227) * MDP(8) + (t134 * t147 + t133 * t146 - g(1) * t195 - g(2) * t210 + (-t139 * t181 + t140 * t182) * qJD(2)) * MDP(9) + (qJDD(1) * t178 + t185 * t203) * MDP(10) + 0.2e1 * (t185 * t212 - t224 * t214) * MDP(11) + t149 * MDP(12) - t150 * MDP(13) + (t191 * t185 + t190 * t186) * MDP(15) + (-t190 * t185 + t191 * t186) * MDP(16) + ((qJDD(1) * t138 + (qJD(4) * t137 - t126) * qJD(1) - t204) * t186 + (qJDD(1) * t137 + (-qJD(4) * t138 + t127) * qJD(1) - t205) * t185 + t198) * MDP(17) + (-t123 * t138 + t129 * t126 + t122 * t137 + t125 * t127 + t124 * (t143 + t173) + t132 * (qJD(2) * t181 - t185 * t233) - g(1) * (-t141 * t184 + t142 * t165 + t195) - g(2) * (-t141 * t165 - t142 * t184 + t210)) * MDP(18); -qJDD(1) * MDP(4) - t189 * MDP(5) + (-qJ(2) * t189 + t202) * MDP(6) + (-qJDD(1) * MDP(7) + (-qJD(1) * t140 + t133) * MDP(9) + (0.2e1 * t209 - t212) * MDP(15) + (t203 + t213) * MDP(16) + (t125 * t220 - t129 * t219 - t124) * MDP(18) - t238 * t189) * t182 + (-t189 * MDP(7) + (qJD(1) * t139 + t134) * MDP(9) - t132 * t217 + (-t222 * MDP(15) - qJDD(4) * MDP(16) + t204 * MDP(18)) * t186 + (-qJDD(4) * MDP(15) + t222 * MDP(16) + t205 * MDP(18)) * t185 + t238 * qJDD(1)) * t181 + (-MDP(6) - MDP(9) - MDP(18)) * t225; t229 * MDP(9) + t150 * MDP(15) + t149 * MDP(16) + (t122 * t186 + t123 * t185 + g(3) + (-t125 * t185 + t231) * qJD(4)) * MDP(18); -MDP(12) * t213 - MDP(13) * t212 + qJDD(4) * MDP(14) + (t193 * t185 + t169 + t234) * MDP(15) + (-t229 * t185 + t193 * t186) * MDP(16) + (pkin(4) * t213 + (-t228 + t233) * t219) * MDP(17) + (t228 * t129 + (t234 + t122 + (qJD(1) * t132 + t198) * t185) * pkin(4)) * MDP(18) + (-t185 * t186 * MDP(10) + t224 * MDP(11)) * t189; (t154 + t180 + t194 + t237) * MDP(18) - t189 * t206 + (t231 + (-t125 - t233) * t185) * t217;];
tau = t1;
