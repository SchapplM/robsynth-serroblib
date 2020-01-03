% Calculate vector of inverse dynamics joint torques for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:19
% EndTime: 2019-12-31 16:47:20
% DurationCPUTime: 0.65s
% Computational Cost: add. (358->139), mult. (610->169), div. (0->0), fcn. (262->4), ass. (0->70)
t107 = sin(qJ(3));
t109 = cos(qJ(3));
t133 = qJDD(3) * qJ(4);
t111 = -pkin(1) - pkin(5);
t95 = t111 * qJD(1) + qJD(2);
t152 = t109 * t95;
t157 = t111 * qJDD(1);
t94 = qJDD(2) + t157;
t91 = t107 * t94;
t81 = t133 + t91 + (qJD(4) + t152) * qJD(3);
t148 = qJDD(3) * pkin(3);
t153 = t107 * t95;
t92 = t109 * t94;
t82 = qJD(3) * t153 + qJDD(4) - t148 - t92;
t126 = qJD(4) - t152;
t154 = qJD(3) * pkin(3);
t84 = t126 - t154;
t138 = qJD(3) * qJ(4);
t86 = t138 + t153;
t116 = t81 * t107 - t82 * t109 + (t107 * t84 + t109 * t86) * qJD(3);
t110 = cos(qJ(1));
t103 = g(2) * t110;
t108 = sin(qJ(1));
t104 = g(1) * t108;
t142 = t104 - t103;
t159 = t116 - t142;
t147 = t109 * t110;
t158 = g(2) * t147 + g(3) * t107;
t124 = pkin(3) * t109 + qJ(4) * t107;
t83 = t124 * qJD(3) - t109 * qJD(4) + qJD(2);
t155 = t107 * pkin(3);
t123 = -t109 * qJ(4) + t155;
t93 = qJ(2) + t123;
t80 = qJD(1) * t83 + qJDD(1) * t93;
t156 = 0.2e1 * qJD(3);
t151 = t110 * pkin(1) + t108 * qJ(2);
t150 = pkin(1) * qJDD(1);
t85 = qJD(1) * t93;
t112 = qJD(3) ^ 2;
t146 = t111 * t112;
t113 = qJD(1) ^ 2;
t145 = t113 * qJ(2);
t144 = t113 * t107;
t143 = t85 * qJD(1);
t105 = t107 ^ 2;
t106 = t109 ^ 2;
t141 = -t105 - t106;
t140 = t105 - t106;
t139 = t112 + t113;
t136 = qJD(1) * qJD(3);
t135 = qJDD(1) * qJ(2);
t134 = qJDD(1) * t109;
t132 = qJDD(3) * t107;
t131 = qJDD(3) * t111;
t130 = t92 + t158;
t129 = 0.2e1 * qJD(1) * qJD(2);
t128 = t109 * t136;
t127 = qJDD(2) - t142;
t125 = g(1) * t110 + g(2) * t108;
t122 = -t145 - t104;
t121 = t143 + t104;
t118 = t85 * t156;
t117 = -t125 + t129 + 0.2e1 * t135;
t115 = t117 - t146;
t114 = t125 + t146 - 0.2e1 * t80;
t100 = t110 * qJ(2);
t98 = qJDD(3) * t109;
t96 = t109 * t131;
t90 = t124 * qJD(1);
t1 = [qJDD(1) * MDP(1) + t142 * MDP(2) + t125 * MDP(3) + (t127 - 0.2e1 * t150) * MDP(4) + t117 * MDP(5) + (-(qJDD(2) - t150) * pkin(1) - g(1) * (-t108 * pkin(1) + t100) - g(2) * t151 + (t129 + t135) * qJ(2)) * MDP(6) + (t106 * qJDD(1) - 0.2e1 * t107 * t128) * MDP(7) + 0.2e1 * (-t107 * t134 + t140 * t136) * MDP(8) + (-t112 * t107 + t98) * MDP(9) + (-t112 * t109 - t132) * MDP(10) + (0.2e1 * qJ(2) * t128 + t115 * t107 + t96) * MDP(12) + ((-0.2e1 * qJ(2) * t136 - t131) * t107 + t115 * t109) * MDP(13) + (-t114 * t107 + t109 * t118 + t96) * MDP(14) + (t141 * t157 - t159) * MDP(15) + ((t118 + t131) * t107 + t114 * t109) * MDP(16) + (t80 * t93 + t85 * t83 - g(1) * (-qJ(4) * t147 + t110 * t155 + t100) - g(2) * (t110 * pkin(5) + t151) + (-g(1) * t111 - g(2) * t123) * t108 + t116 * t111) * MDP(17); -t113 * MDP(5) + (t127 - t145) * MDP(6) + (-t143 + t159) * MDP(17) + (MDP(12) + MDP(14)) * (-t139 * t107 + t98) + (-MDP(13) + MDP(16)) * (t139 * t109 + t132) + (t141 * MDP(15) - pkin(1) * MDP(6) + MDP(4)) * qJDD(1); t109 * MDP(7) * t144 - t140 * MDP(8) * t113 + MDP(9) * t134 - t107 * qJDD(1) * MDP(10) + qJDD(3) * MDP(11) + (t122 * t109 + t130) * MDP(12) + (g(3) * t109 - t91 + (-t122 - t103) * t107) * MDP(13) + (-t109 * t104 + 0.2e1 * t148 - qJDD(4) + (-t107 * t90 - t109 * t85) * qJD(1) + t130) * MDP(14) + (-t124 * qJDD(1) + ((t86 - t138) * t109 + (-qJD(4) + t84 + t154) * t107) * qJD(1)) * MDP(15) + (0.2e1 * t133 + qJD(4) * t156 + t91 + (qJD(1) * t90 - g(3)) * t109 + (-t121 + t103) * t107) * MDP(16) + (-t82 * pkin(3) + g(3) * t123 + t81 * qJ(4) - t142 * t124 + t126 * t86 - t84 * t153 - t85 * t90) * MDP(17); -qJDD(3) * MDP(14) + (-t106 * t113 - t112) * MDP(16) + (-t86 * qJD(3) - t158 + t82) * MDP(17) + (MDP(14) * t144 + qJDD(1) * MDP(15) + t121 * MDP(17)) * t109;];
tau = t1;
