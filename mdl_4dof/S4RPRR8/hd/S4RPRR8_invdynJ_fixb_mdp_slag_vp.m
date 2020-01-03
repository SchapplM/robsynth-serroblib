% Calculate vector of inverse dynamics joint torques for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:18
% EndTime: 2019-12-31 16:55:19
% DurationCPUTime: 0.90s
% Computational Cost: add. (483->153), mult. (919->209), div. (0->0), fcn. (565->8), ass. (0->81)
t164 = sin(qJ(1));
t167 = cos(qJ(1));
t214 = -g(1) * t164 + g(2) * t167;
t215 = qJDD(2) + t214;
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t193 = qJDD(1) * t166;
t196 = qJD(1) * qJD(3);
t213 = t163 * t196 - t193;
t208 = pkin(1) * qJDD(1);
t211 = t208 - t215;
t158 = qJD(3) + qJD(4);
t168 = -pkin(1) - pkin(5);
t209 = pkin(6) - t168;
t170 = qJD(1) ^ 2;
t207 = qJ(2) * t170;
t142 = t168 * qJD(1) + qJD(2);
t202 = qJD(1) * t163;
t122 = -pkin(6) * t202 + t142 * t163;
t165 = cos(qJ(4));
t206 = t122 * t165;
t162 = sin(qJ(4));
t130 = t162 * t166 + t163 * t165;
t117 = t158 * t130;
t131 = -t162 * t163 + t165 * t166;
t157 = qJDD(3) + qJDD(4);
t205 = -t117 * t158 + t131 * t157;
t160 = t166 ^ 2;
t204 = t163 ^ 2 - t160;
t169 = qJD(3) ^ 2;
t203 = -t169 - t170;
t201 = qJD(1) * t166;
t200 = qJD(3) * t163;
t199 = qJD(3) * t166;
t198 = qJD(4) * t162;
t197 = qJD(1) * qJD(2);
t195 = qJDD(1) * qJ(2);
t194 = qJDD(1) * t163;
t192 = qJDD(3) * t163;
t191 = t163 * t198;
t135 = t209 * t166;
t189 = t166 * t196;
t148 = pkin(3) * t163 + qJ(2);
t187 = t158 * t166;
t123 = -pkin(6) * t201 + t166 * t142;
t186 = g(1) * t167 + g(2) * t164;
t118 = -t162 * t200 + t165 * t187 - t191;
t185 = -t118 * t158 - t130 * t157;
t120 = qJD(3) * pkin(3) + t123;
t184 = -t120 * t162 - t206;
t134 = t209 * t163;
t183 = -t134 * t165 - t135 * t162;
t182 = t134 * t162 - t135 * t165;
t181 = -qJD(1) * t191 - t213 * t162;
t179 = 0.2e1 * qJ(2) * t196 + qJDD(3) * t168;
t178 = t189 + t194;
t177 = -t207 + t214;
t176 = -t186 + 0.2e1 * t197;
t110 = -t117 * qJD(1) - t162 * t194 + t165 * t193;
t125 = t162 * t202 - t165 * t201;
t126 = t130 * qJD(1);
t175 = -t125 * t126 * MDP(14) + (t126 * t158 + t110) * MDP(16) + (-t125 * t158 + (-t158 * t201 - t194) * t165 - t181) * MDP(17) + (t125 ^ 2 - t126 ^ 2) * MDP(15) + t157 * MDP(18);
t174 = t176 + 0.2e1 * t195;
t173 = -t168 * t169 + t174;
t141 = t168 * qJDD(1) + qJDD(2);
t132 = t166 * t141;
t114 = qJDD(3) * pkin(3) + t213 * pkin(6) - t142 * t200 + t132;
t136 = t148 * qJD(1);
t161 = qJ(3) + qJ(4);
t154 = sin(t161);
t155 = cos(t161);
t172 = t122 * t198 + g(3) * t155 + (-t122 * t158 - t114) * t162 + t136 * t126 - t214 * t154;
t115 = -t178 * pkin(6) + t141 * t163 + t142 * t199;
t171 = g(3) * t154 + t184 * qJD(4) + t165 * t114 - t162 * t115 + t136 * t125 + t214 * t155;
t153 = qJDD(3) * t166;
t143 = pkin(3) * t199 + qJD(2);
t129 = qJD(3) * t135;
t128 = t209 * t200;
t121 = t178 * pkin(3) + t195 + t197;
t111 = (qJD(1) * t187 + t194) * t165 + t181;
t1 = [qJDD(1) * MDP(1) - t214 * MDP(2) + t186 * MDP(3) + (-0.2e1 * t208 + t215) * MDP(4) + t174 * MDP(5) + (t211 * pkin(1) + (t176 + t195) * qJ(2)) * MDP(6) + (qJDD(1) * t160 - 0.2e1 * t163 * t189) * MDP(7) + 0.2e1 * (-t163 * t193 + t204 * t196) * MDP(8) + (-t163 * t169 + t153) * MDP(9) + (-t166 * t169 - t192) * MDP(10) + (t173 * t163 + t179 * t166) * MDP(12) + (-t179 * t163 + t173 * t166) * MDP(13) + (t110 * t131 + t117 * t125) * MDP(14) + (-t110 * t130 - t111 * t131 + t117 * t126 + t118 * t125) * MDP(15) + t205 * MDP(16) + t185 * MDP(17) + (t143 * t126 + t148 * t111 + t121 * t130 + t136 * t118 + (-t183 * qJD(4) + t128 * t165 + t129 * t162) * t158 + t182 * t157 - t186 * t154) * MDP(19) + (-t143 * t125 + t148 * t110 + t121 * t131 - t136 * t117 - (t182 * qJD(4) + t128 * t162 - t129 * t165) * t158 - t183 * t157 - t186 * t155) * MDP(20); qJDD(1) * MDP(4) - t170 * MDP(5) + (-t207 - t211) * MDP(6) + (t203 * t163 + t153) * MDP(12) + (t203 * t166 - t192) * MDP(13) + (-qJD(1) * t126 + t205) * MDP(19) + (qJD(1) * t125 + t185) * MDP(20); MDP(9) * t193 - MDP(10) * t194 + qJDD(3) * MDP(11) + (g(3) * t163 + t177 * t166 + t132) * MDP(12) + (g(3) * t166 + (-t141 - t177) * t163) * MDP(13) + (-(-t123 * t162 - t206) * t158 + (-t126 * t201 + t165 * t157 - t158 * t198) * pkin(3) + t171) * MDP(19) + ((-qJD(4) * t120 + t123 * t158 - t115) * t165 + (-qJD(4) * t165 * t158 + t125 * t201 - t162 * t157) * pkin(3) + t172) * MDP(20) + t175 + (t166 * t163 * MDP(7) - t204 * MDP(8)) * t170; (-t184 * t158 + t171) * MDP(19) + ((-t115 + (-qJD(4) + t158) * t120) * t165 + t172) * MDP(20) + t175;];
tau = t1;
