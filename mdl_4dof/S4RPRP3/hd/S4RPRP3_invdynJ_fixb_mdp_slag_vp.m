% Calculate vector of inverse dynamics joint torques for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:39
% EndTime: 2021-01-15 10:20:40
% DurationCPUTime: 0.79s
% Computational Cost: add. (416->154), mult. (801->197), div. (0->0), fcn. (413->8), ass. (0->80)
t179 = qJDD(2) - g(3);
t129 = sin(qJ(3));
t124 = t129 ^ 2;
t131 = cos(qJ(3));
t125 = t131 ^ 2;
t178 = (t124 - t125) * MDP(6);
t126 = sin(pkin(6));
t115 = pkin(1) * t126 + pkin(5);
t106 = t115 * qJD(1);
t157 = qJ(4) * qJD(1);
t147 = t106 + t157;
t142 = t147 * t131;
t127 = cos(pkin(6));
t176 = pkin(1) * t127;
t175 = pkin(3) * t124;
t174 = pkin(3) * t131;
t123 = qJ(1) + pkin(6);
t118 = sin(t123);
t173 = g(1) * t118;
t119 = cos(t123);
t172 = g(1) * t119;
t170 = g(2) * t118;
t169 = g(2) * t119;
t167 = g(3) * t131;
t165 = qJD(3) * pkin(3);
t121 = t131 * qJD(2);
t93 = -t129 * t147 + t121;
t91 = t93 + t165;
t166 = t91 - t93;
t164 = qJDD(3) * pkin(3);
t163 = t106 * t129;
t162 = t118 * t131;
t161 = t119 * t129;
t134 = qJD(1) ^ 2;
t160 = t131 * t134;
t159 = qJ(4) + t115;
t117 = pkin(2) + t174;
t101 = -t117 - t176;
t156 = qJD(1) * t101;
t116 = -pkin(2) - t176;
t107 = qJD(1) * t116;
t155 = qJD(1) * t107;
t154 = qJD(1) * qJD(3);
t153 = qJDD(1) * t101;
t152 = qJDD(1) * t129;
t151 = qJDD(1) * t131;
t150 = qJDD(3) * t115;
t149 = t129 * pkin(3) * t154 + qJDD(4);
t148 = qJD(3) * t159;
t92 = t149 + t153;
t146 = t92 + t153;
t104 = t115 * qJDD(1);
t145 = -qJD(3) * qJD(2) - t104;
t144 = -t170 - t172;
t94 = qJD(2) * t129 + t142;
t143 = t129 * t91 - t131 * t94;
t120 = t131 * qJDD(2);
t141 = g(1) * t161 + t129 * t170 + t120 - t167;
t133 = qJD(3) ^ 2;
t140 = 0.2e1 * qJDD(1) * t116 + t115 * t133;
t139 = -qJ(4) * qJDD(1) + t145;
t100 = qJD(3) * t163;
t138 = g(2) * t162 - t179 * t129 + t131 * t172 + t100;
t137 = qJD(1) * qJD(4) - t139;
t97 = qJD(4) + t156;
t136 = (-qJD(4) - t97) * qJD(1) + t139;
t132 = cos(qJ(1));
t130 = sin(qJ(1));
t128 = -qJ(4) - pkin(5);
t111 = g(1) * t162;
t110 = g(2) * t161;
t103 = qJDD(3) * t131 - t129 * t133;
t102 = qJDD(3) * t129 + t131 * t133;
t99 = t159 * t131;
t98 = t159 * t129;
t96 = -qJD(4) * t129 - t131 * t148;
t95 = qJD(4) * t131 - t129 * t148;
t90 = -t100 + (-qJ(4) * t154 + qJDD(2)) * t129 + t137 * t131;
t89 = -qJD(3) * t142 - t129 * t137 + t120 + t164;
t1 = [(g(1) * t132 + g(2) * t130) * MDP(3) + t102 * MDP(7) + t103 * MDP(8) + t111 * MDP(10) + t110 * MDP(11) + (-qJDD(3) * t98 + t111) * MDP(12) + (-qJDD(3) * t99 + t110) * MDP(13) + t144 * MDP(14) + (t90 * t99 + t94 * t95 - t89 * t98 + t91 * t96 + t92 * t101 - g(1) * (-pkin(1) * t130 - t117 * t118 - t119 * t128) - g(2) * (pkin(1) * t132 + t117 * t119 - t118 * t128)) * MDP(15) + (t124 * MDP(5) + MDP(1) + (t126 ^ 2 + t127 ^ 2) * MDP(4) * pkin(1) ^ 2) * qJDD(1) + (t96 * MDP(12) - t95 * MDP(13) + (MDP(13) * t175 - 0.2e1 * t178) * qJD(1)) * qJD(3) + ((-t140 - t169) * MDP(10) - MDP(11) * t150 + (-t146 - t169) * MDP(12) + (qJD(1) * t95 + qJDD(1) * t99 + t90) * MDP(14) + (0.2e1 * t107 * MDP(11) + (t97 + t156) * MDP(13) + (qJD(1) * t98 - t91) * MDP(14)) * qJD(3)) * t131 + (0.2e1 * MDP(6) * t151 - MDP(10) * t150 + (t140 - t173) * MDP(11) + (t146 - t173) * MDP(13) + (-qJD(1) * t96 + qJDD(1) * t98 - t89) * MDP(14) + (t107 * MDP(10) - t94 * MDP(14) + (pkin(3) * MDP(15) + MDP(12)) * t97 + (0.2e1 * t131 * MDP(5) + t116 * MDP(10) + (t101 - t174) * MDP(12) - t99 * MDP(14)) * qJD(1)) * qJD(3)) * t129 + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t130 - g(2) * t132); t179 * MDP(4) + (-qJD(3) * t143 + t129 * t90 + t131 * t89 - g(3)) * MDP(15) + (MDP(10) + MDP(12)) * t103 + (-MDP(11) - MDP(13)) * t102; -t129 * MDP(5) * t160 + t134 * t178 + MDP(7) * t152 + MDP(8) * t151 + qJDD(3) * MDP(9) + ((-t104 - t155) * t129 + t141) * MDP(10) + ((t121 - t163) * qJD(3) + (t145 - t155) * t131 + t138) * MDP(11) + (0.2e1 * t164 + (t94 - t142) * qJD(3) + (pkin(3) * t160 + t136) * t129 + t141) * MDP(12) + (-t134 * t175 + (t129 * t157 + t93) * qJD(3) + t136 * t131 + t138) * MDP(13) + (-pkin(3) * t152 + (-t165 + t166) * t131 * qJD(1)) * MDP(14) + (t166 * t94 + (-t167 + t89 + (-qJD(1) * t97 - t144) * t129) * pkin(3)) * MDP(15); (t149 + t169 - t173) * MDP(15) + (-t124 - t125) * MDP(14) * t134 + (t143 * MDP(15) + 0.2e1 * (MDP(12) * t129 + MDP(13) * t131) * qJD(3)) * qJD(1) + (-t131 * MDP(12) + t129 * MDP(13) + MDP(15) * t101) * qJDD(1);];
tau = t1;
