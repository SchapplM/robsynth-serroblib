% Calculate vector of inverse dynamics joint torques for
% S4RPRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:21
% EndTime: 2019-12-31 16:49:22
% DurationCPUTime: 0.74s
% Computational Cost: add. (510->138), mult. (1034->204), div. (0->0), fcn. (667->12), ass. (0->79)
t171 = sin(pkin(7));
t154 = pkin(1) * t171 + pkin(5);
t214 = pkin(6) + t154;
t167 = qJ(1) + pkin(7);
t159 = sin(t167);
t160 = cos(t167);
t217 = g(1) * t160 + g(2) * t159;
t174 = sin(qJ(3));
t177 = cos(qJ(3));
t197 = t214 * qJD(1);
t122 = t177 * qJD(2) - t197 * t174;
t166 = qJD(3) + qJD(4);
t123 = qJD(2) * t174 + t177 * t197;
t213 = qJD(3) * pkin(3);
t176 = cos(qJ(4));
t212 = t123 * t176;
t211 = qJDD(2) - g(3);
t168 = t174 ^ 2;
t210 = -t177 ^ 2 + t168;
t172 = cos(pkin(7));
t155 = -pkin(1) * t172 - pkin(2);
t145 = qJD(1) * t155;
t209 = qJD(1) * t174;
t208 = qJD(1) * t177;
t173 = sin(qJ(4));
t206 = qJD(4) * t173;
t205 = qJD(1) * qJD(3);
t204 = qJDD(1) * t174;
t203 = qJDD(1) * t177;
t202 = t174 * t213;
t201 = t173 * t209;
t200 = t176 * t208;
t199 = t177 * t205;
t198 = qJD(3) * t214;
t142 = t154 * qJDD(1);
t196 = pkin(6) * qJDD(1) + t142;
t195 = g(1) * t159 - g(2) * t160;
t175 = sin(qJ(1));
t178 = cos(qJ(1));
t194 = g(1) * t175 - g(2) * t178;
t193 = t173 * t204 - t176 * t203;
t133 = t173 * t174 - t176 * t177;
t118 = t166 * t133;
t134 = t173 * t177 + t174 * t176;
t165 = qJDD(3) + qJDD(4);
t191 = -t118 * t166 + t134 * t165;
t121 = t122 + t213;
t190 = -t121 * t173 - t212;
t131 = t214 * t174;
t132 = t214 * t177;
t189 = -t131 * t176 - t132 * t173;
t188 = -t131 * t173 + t132 * t176;
t139 = -pkin(3) * t177 + t155;
t115 = qJD(4) * t200 - t166 * t201 + t173 * t203 + (t199 + t204) * t176;
t127 = -t200 + t201;
t129 = -t173 * t208 - t176 * t209;
t187 = -t129 * t127 * MDP(12) + (t127 * t166 + t115) * MDP(14) + (-t193 + (-qJD(1) * t134 - t129) * t166) * MDP(15) + (-t127 ^ 2 + t129 ^ 2) * MDP(13) + t165 * MDP(16);
t186 = -qJD(1) * t145 - t142 + t217;
t185 = 0.2e1 * qJD(3) * t145 - qJDD(3) * t154;
t179 = qJD(3) ^ 2;
t184 = -0.2e1 * qJDD(1) * t155 - t154 * t179 + t195;
t161 = t177 * qJDD(2);
t112 = qJDD(3) * pkin(3) - qJD(3) * t123 - t196 * t174 + t161;
t130 = t139 * qJD(1);
t170 = qJ(3) + qJ(4);
t163 = sin(t170);
t164 = cos(t170);
t183 = t130 * t127 + t123 * t206 + g(3) * t163 + (-t123 * t166 - t112) * t173 + t217 * t164;
t119 = t166 * t134;
t113 = qJD(3) * t122 + t174 * qJDD(2) + t196 * t177;
t182 = -g(3) * t164 + qJD(4) * t190 + t176 * t112 - t173 * t113 + t130 * t129 + t163 * t217;
t141 = qJDD(3) * t177 - t174 * t179;
t140 = qJDD(3) * t174 + t177 * t179;
t126 = t177 * t198;
t125 = t174 * t198;
t124 = qJD(1) * t202 + qJDD(1) * t139;
t116 = qJD(1) * t119 + t193;
t114 = -t119 * t166 - t133 * t165;
t1 = [qJDD(1) * MDP(1) + t194 * MDP(2) + (g(1) * t178 + g(2) * t175) * MDP(3) + (t194 + (t171 ^ 2 + t172 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t168 + 0.2e1 * t174 * t199) * MDP(5) + 0.2e1 * (t174 * t203 - t205 * t210) * MDP(6) + t140 * MDP(7) + t141 * MDP(8) + (t174 * t185 + t177 * t184) * MDP(10) + (-t174 * t184 + t177 * t185) * MDP(11) + (t115 * t134 + t118 * t129) * MDP(12) + (-t115 * t133 - t116 * t134 + t118 * t127 + t119 * t129) * MDP(13) + t191 * MDP(14) + t114 * MDP(15) + (t127 * t202 + t139 * t116 + t124 * t133 + t130 * t119 + (-qJD(4) * t188 + t125 * t173 - t126 * t176) * t166 + t189 * t165 + t195 * t164) * MDP(17) + (-t129 * t202 + t139 * t115 + t124 * t134 - t130 * t118 - (qJD(4) * t189 - t125 * t176 - t126 * t173) * t166 - t188 * t165 - t195 * t163) * MDP(18); MDP(10) * t141 - MDP(11) * t140 + MDP(17) * t114 - MDP(18) * t191 + MDP(4) * t211; MDP(7) * t204 + MDP(8) * t203 + qJDD(3) * MDP(9) + (-g(3) * t177 + t174 * t186 + t161) * MDP(10) + (-t174 * t211 + t177 * t186) * MDP(11) + (-(-t122 * t173 - t212) * t166 + (-t127 * t209 + t165 * t176 - t166 * t206) * pkin(3) + t182) * MDP(17) + ((-qJD(4) * t121 + t122 * t166 - t113) * t176 + (-qJD(4) * t166 * t176 + t129 * t209 - t165 * t173) * pkin(3) + t183) * MDP(18) + t187 + (-t174 * t177 * MDP(5) + MDP(6) * t210) * qJD(1) ^ 2; (-t166 * t190 + t182) * MDP(17) + ((-t113 + (-qJD(4) + t166) * t121) * t176 + t183) * MDP(18) + t187;];
tau = t1;
