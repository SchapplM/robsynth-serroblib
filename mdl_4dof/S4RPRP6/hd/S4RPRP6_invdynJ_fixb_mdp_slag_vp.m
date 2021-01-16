% Calculate vector of inverse dynamics joint torques for
% S4RPRP6
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
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:43
% EndTime: 2021-01-15 10:27:45
% DurationCPUTime: 0.71s
% Computational Cost: add. (385->148), mult. (673->177), div. (0->0), fcn. (294->4), ass. (0->70)
t132 = -pkin(1) - pkin(5);
t109 = t132 * qJD(1) + qJD(2);
t175 = -qJ(4) * qJD(1) + t109;
t124 = (qJDD(1) * qJ(2));
t125 = (qJD(1) * qJD(2));
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t141 = g(1) * t131 + g(2) * t129;
t137 = -t141 + (2 * t125);
t178 = (2 * t124) + t137;
t121 = g(2) * t131;
t122 = g(1) * t129;
t164 = t122 - t121;
t147 = qJDD(2) - t164;
t130 = cos(qJ(3));
t154 = qJD(1) * qJD(3);
t148 = t130 * t154;
t128 = sin(qJ(3));
t153 = qJDD(1) * t128;
t95 = qJDD(4) + t124 + t125 + (t148 + t153) * pkin(3);
t177 = t95 - t141;
t176 = qJ(4) - t132;
t100 = t175 * t130;
t96 = qJD(3) * pkin(3) + t100;
t171 = pkin(1) * qJDD(1);
t134 = qJD(1) ^ 2;
t170 = qJ(2) * t134;
t99 = t175 * t128;
t169 = qJD(3) * t99;
t168 = qJDD(3) * pkin(3);
t167 = t130 * t134;
t165 = g(3) * t130 + t128 * t122;
t126 = t128 ^ 2;
t127 = t130 ^ 2;
t163 = t126 - t127;
t133 = qJD(3) ^ 2;
t162 = -t133 - t134;
t115 = pkin(3) * t128 + qJ(2);
t159 = qJD(1) * t115;
t105 = qJD(4) + t159;
t160 = qJD(1) * t105;
t107 = t176 * t130;
t158 = qJD(3) * t107;
t157 = qJD(3) * t128;
t156 = qJD(3) * t130;
t155 = qJ(4) * qJDD(1);
t152 = qJDD(1) * t130;
t151 = qJDD(3) * t128;
t108 = t132 * qJDD(1) + qJDD(2);
t103 = t130 * t108;
t150 = g(3) * t128 + t130 * t121 + t103;
t149 = t128 * t154;
t146 = -t108 - t121;
t145 = (-t126 - t127) * MDP(16);
t143 = t105 + t159;
t142 = -0.2e1 * t149;
t140 = -t160 - t122;
t139 = -qJD(1) * qJD(4) - t155;
t138 = 0.2e1 * qJ(2) * t154 + qJDD(3) * t132;
t110 = pkin(3) * t156 + qJD(2);
t136 = qJD(1) * t110 + qJDD(1) * t115 + t177;
t135 = -t132 * t133 + t178;
t118 = qJDD(3) * t130;
t111 = qJ(4) * t149;
t106 = t176 * t128;
t98 = -qJD(4) * t128 - t158;
t97 = -qJD(4) * t130 + t157 * t176;
t94 = t175 * t156 + (t108 + t139) * t128;
t93 = -t109 * t157 + t139 * t130 + t103 + t111 + t168;
t1 = [qJDD(1) * MDP(1) + t164 * MDP(2) + t141 * MDP(3) + (t147 - 0.2e1 * t171) * MDP(4) + t178 * MDP(5) + ((t171 - t147) * pkin(1) + (t137 + t124) * qJ(2)) * MDP(6) + (qJDD(1) * t127 + t130 * t142) * MDP(7) + 0.2e1 * (-t128 * t152 + t163 * t154) * MDP(8) + (-t128 * t133 + t118) * MDP(9) + (-t130 * t133 - t151) * MDP(10) + (t135 * t128 + t138 * t130) * MDP(12) + (-t138 * t128 + t135 * t130) * MDP(13) + (-qJDD(3) * t107 + (t143 * t130 + t97) * qJD(3) + t136 * t128) * MDP(14) + (qJDD(3) * t106 + (-t143 * t128 - t98) * qJD(3) + t136 * t130) * MDP(15) + ((-t169 + qJDD(1) * t107 - t93 + (qJD(3) * t106 - t97) * qJD(1)) * t130 + (qJD(3) * t96 + qJDD(1) * t106 - t94 + (-t98 - t158) * qJD(1)) * t128 + t164) * MDP(16) + (-t94 * t106 + t99 * t98 - t93 * t107 + t96 * t97 + t95 * t115 + t105 * t110 - g(1) * (t115 * t131 - t129 * t176) - g(2) * (t115 * t129 + t131 * t176)) * MDP(17); -t134 * MDP(5) + (t147 - t170) * MDP(6) + (-t160 + t128 * t94 + t130 * t93 + (-t128 * t96 + t130 * t99) * qJD(3) - t164) * MDP(17) + (MDP(12) + MDP(14)) * (t162 * t128 + t118) + (MDP(13) + MDP(15)) * (t162 * t130 - t151) + (-pkin(1) * MDP(6) + MDP(4) + t145) * qJDD(1); qJDD(3) * MDP(11) + t150 * MDP(12) + t165 * MDP(13) + (t111 + t150 + 0.2e1 * t168 + t169) * MDP(14) + (qJD(3) * t100 + t165) * MDP(15) + (pkin(3) * t93 + (-t100 + t96) * t99) * MDP(17) + (-pkin(3) * t127 * MDP(15) - t163 * MDP(8)) * t134 + (qJDD(1) * MDP(9) + (-t170 - t122) * MDP(12) + (t139 + t140) * MDP(14) - t175 * MDP(15) * qJD(3) + (-qJDD(1) * MDP(16) + (t140 + t121) * MDP(17)) * pkin(3)) * t130 + (MDP(7) * t167 - qJDD(1) * MDP(10) + (t146 + t170) * MDP(13) + (-pkin(3) * t167 - qJD(3) * t109) * MDP(14) + pkin(3) * g(3) * MDP(17) + (t146 + t155 + (qJD(4) + t105) * qJD(1)) * MDP(15)) * t128; (0.2e1 * t148 + t153) * MDP(14) + (t142 + t152) * MDP(15) + ((t128 * t99 + t130 * t96) * qJD(1) + t177) * MDP(17) + t134 * t145;];
tau = t1;
