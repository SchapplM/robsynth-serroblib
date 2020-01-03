% Calculate vector of inverse dynamics joint torques for
% S4RPRP4
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
%   see S4RPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:56
% EndTime: 2019-12-31 16:43:57
% DurationCPUTime: 0.63s
% Computational Cost: add. (414->126), mult. (767->165), div. (0->0), fcn. (391->8), ass. (0->68)
t172 = 2 * qJD(3);
t125 = cos(qJ(3));
t121 = sin(pkin(6));
t110 = pkin(1) * t121 + pkin(5);
t105 = t110 * qJD(1);
t123 = sin(qJ(3));
t161 = t105 * t123;
t93 = qJD(2) * t125 - t161;
t171 = qJD(4) - t93;
t103 = t110 * qJDD(1);
t170 = qJD(2) * qJD(3) + t103;
t118 = qJ(1) + pkin(6);
t114 = sin(t118);
t115 = cos(t118);
t142 = g(1) * t115 + g(2) * t114;
t90 = -(qJD(3) * pkin(3)) + t171;
t160 = t105 * t125;
t94 = qJD(2) * t123 + t160;
t91 = (qJD(3) * qJ(4)) + t94;
t139 = pkin(3) * t123 - qJ(4) * t125;
t95 = t139 * qJD(3) - qJD(4) * t123;
t140 = pkin(3) * t125 + qJ(4) * t123;
t137 = pkin(2) + t140;
t122 = cos(pkin(6));
t167 = pkin(1) * t122;
t96 = -t137 - t167;
t89 = qJD(1) * t95 + qJDD(1) * t96;
t162 = qJDD(3) * pkin(3);
t169 = qJDD(4) - t162;
t149 = qJDD(3) * t110;
t92 = qJD(1) * t96;
t168 = t92 * t172 - t149;
t166 = g(1) * t114;
t159 = t123 * t125;
t119 = t123 ^ 2;
t120 = t125 ^ 2;
t158 = t119 - t120;
t111 = -pkin(2) - t167;
t106 = qJD(1) * t111;
t157 = qJD(1) * t123;
t154 = qJD(1) * qJD(3);
t152 = qJDD(1) * t123;
t151 = qJDD(1) * t125;
t150 = qJDD(3) * qJ(4);
t148 = t123 * qJDD(2) + t125 * t170;
t128 = qJD(1) ^ 2;
t147 = t128 * t159;
t145 = t93 + t161;
t143 = qJD(3) * t160 - t125 * qJDD(2) + t170 * t123;
t124 = sin(qJ(1));
t126 = cos(qJ(1));
t141 = g(1) * t124 - g(2) * t126;
t127 = qJD(3) ^ 2;
t138 = g(2) * t115 + t110 * t127;
t136 = t141 * pkin(1);
t135 = -g(3) * t125 + t142 * t123 - t143;
t134 = -0.2e1 * qJDD(1) * t111 - t138;
t133 = t106 * t172 - t149;
t132 = qJD(3) * t94 + t135;
t131 = -t138 - 0.2e1 * t89;
t87 = t150 + (qJD(4) - t161) * qJD(3) + t148;
t88 = t143 + t169;
t130 = t88 * t123 + t87 * t125 + (-t123 * t91 + t125 * t90) * qJD(3);
t108 = t125 * t166;
t102 = qJDD(3) * t125 - t123 * t127;
t101 = qJDD(3) * t123 + t125 * t127;
t100 = t139 * qJD(1);
t1 = [qJDD(1) * MDP(1) + t141 * MDP(2) + (g(1) * t126 + g(2) * t124) * MDP(3) + ((t121 ^ 2 + t122 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t136) * MDP(4) + (qJDD(1) * t119 + 0.2e1 * t154 * t159) * MDP(5) + 0.2e1 * (t123 * t151 - t158 * t154) * MDP(6) + t101 * MDP(7) + t102 * MDP(8) + (t133 * t123 + t134 * t125 + t108) * MDP(10) + (t133 * t125 + (-t134 - t166) * t123) * MDP(11) + (t168 * t123 + t131 * t125 + t108) * MDP(12) + ((t119 + t120) * t103 + t130 - t142) * MDP(13) + (-t168 * t125 + (t131 + t166) * t123) * MDP(14) + (t89 * t96 + t92 * t95 + t136 + (-g(1) * pkin(5) - g(2) * t137) * t115 + (-g(2) * pkin(5) + g(1) * t137) * t114 + t130 * t110) * MDP(15); (qJDD(2) - g(3)) * MDP(4) + (t123 * t87 - t125 * t88 - g(3) + (t123 * t90 + t125 * t91) * qJD(3)) * MDP(15) + (MDP(10) + MDP(12)) * t102 + (-MDP(11) + MDP(14)) * t101; -MDP(5) * t147 + t158 * MDP(6) * t128 + MDP(7) * t152 + MDP(8) * t151 + (qJDD(3) * MDP(9)) + (-t106 * t157 + t132) * MDP(10) + (g(3) * t123 + t145 * qJD(3) + (-qJD(1) * t106 + t142) * t125 - t148) * MDP(11) + ((2 * t162) - qJDD(4) + (t100 * t125 - t123 * t92) * qJD(1) + t132) * MDP(12) - t139 * qJDD(1) * MDP(13) + ((2 * t150) + (qJD(1) * t100 - g(3)) * t123 + (qJD(1) * t92 - t142) * t125 + (0.2e1 * qJD(4) - t145) * qJD(3) + t148) * MDP(14) + (-t88 * pkin(3) - g(3) * t140 + t87 * qJ(4) - t92 * t100 + t142 * t139 + t171 * t91 - t90 * t94) * MDP(15); (-qJDD(3) - t147) * MDP(12) + MDP(13) * t152 + (-t119 * t128 - t127) * MDP(14) + (-qJD(3) * t91 + t92 * t157 - t135 + t169) * MDP(15);];
tau = t1;
