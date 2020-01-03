% Calculate vector of inverse dynamics joint torques for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:05
% EndTime: 2019-12-31 17:51:07
% DurationCPUTime: 0.86s
% Computational Cost: add. (512->150), mult. (799->185), div. (0->0), fcn. (380->8), ass. (0->71)
t127 = qJD(3) * qJD(1);
t125 = qJ(1) + pkin(7);
t122 = sin(t125);
t123 = cos(t125);
t131 = sin(pkin(7));
t115 = pkin(1) * t131 + qJ(3);
t167 = qJDD(1) * t115;
t144 = -g(1) * t123 - g(2) * t122 + t167;
t180 = t127 + t144;
t179 = -g(1) * t122 + g(2) * t123;
t134 = sin(qJ(4));
t128 = t134 ^ 2;
t136 = cos(qJ(4));
t129 = t136 ^ 2;
t178 = MDP(9) * (t128 - t129);
t132 = cos(pkin(7));
t118 = -pkin(1) * t132 - pkin(2);
t111 = -pkin(6) + t118;
t169 = qJ(5) - t111;
t103 = t169 * t134;
t102 = qJD(5) + (pkin(4) * t134 + t115) * qJD(1);
t177 = -qJD(1) * t102 + t179;
t110 = qJD(1) * t115;
t143 = -qJD(1) * t110 + t179;
t160 = qJDD(1) * t134;
t161 = qJD(1) * qJD(4);
t176 = qJDD(5) + (t136 * t161 + t160) * pkin(4);
t135 = sin(qJ(1));
t173 = g(1) * t135;
t137 = cos(qJ(1));
t124 = t137 * pkin(1);
t170 = qJD(4) * pkin(4);
t106 = qJD(1) * t111 + qJD(3);
t101 = t136 * t106;
t164 = qJ(5) * qJD(1);
t96 = -qJD(2) * t134 - t136 * t164 + t101;
t95 = t96 + t170;
t172 = t95 - t96;
t171 = MDP(16) * pkin(4);
t168 = qJDD(2) - g(3);
t138 = qJD(4) ^ 2;
t139 = qJD(1) ^ 2;
t165 = -t138 - t139;
t163 = qJD(1) * t134;
t162 = qJD(2) * t136;
t159 = qJDD(1) * t136;
t158 = qJDD(4) * MDP(13);
t157 = qJDD(4) * MDP(14);
t156 = t123 * pkin(2) + t122 * qJ(3) + t124;
t107 = t127 + t167;
t104 = t169 * t136;
t154 = (-t128 - t129) * MDP(15);
t153 = -t106 + t164;
t152 = 0.2e1 * t110;
t151 = t118 * MDP(7);
t97 = -qJ(5) * t163 + t106 * t134 + t162;
t149 = t134 * t97 + t136 * t95;
t146 = -pkin(1) * t135 - pkin(2) * t122 + t123 * qJ(3);
t145 = -qJ(5) * qJDD(1) - qJD(1) * qJD(5);
t142 = t176 + t180;
t141 = -t111 * t138 + t107 + t180;
t133 = -qJ(5) - pkin(6);
t109 = qJDD(4) * t136 - t134 * t138;
t108 = -qJDD(4) * t134 - t136 * t138;
t105 = qJDD(1) * t111 + qJDD(3);
t100 = t136 * t105;
t99 = -qJD(4) * t104 - qJD(5) * t134;
t98 = qJD(4) * t103 - qJD(5) * t136;
t94 = (-qJD(4) * t153 + qJDD(2)) * t136 + (-qJD(2) * qJD(4) + t105 + t145) * t134;
t93 = qJDD(4) * pkin(4) - t134 * qJDD(2) + t100 + t145 * t136 + (t134 * t153 - t162) * qJD(4);
t1 = [(-g(2) * t137 + t173) * MDP(2) + (g(1) * t137 + g(2) * t135) * MDP(3) + (pkin(1) * t173 - g(2) * t124) * MDP(4) + (qJDD(3) + t179) * MDP(5) + (0.2e1 * t127 + t144) * MDP(6) + (-g(1) * t146 - g(2) * t156 + t110 * qJD(3) + qJDD(3) * t118 + t107 * t115) * MDP(7) + t109 * MDP(10) + t108 * MDP(11) - t179 * MDP(15) + (-t94 * t103 + t97 * t99 - t93 * t104 + t95 * t98 + (t107 + t176) * t115 + t102 * qJD(3) - g(1) * (t122 * t133 + t146) - g(2) * (-t123 * t133 + t156)) * MDP(16) + 0.2e1 * t161 * t178 + (MDP(1) + t115 * MDP(6) + t129 * MDP(8) + (t131 ^ 2 + t132 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * MDP(5) + t151) * t118) * qJDD(1) + (t141 * MDP(13) - t111 * t157 + (-qJD(1) * t99 + qJDD(1) * t103 - t94) * MDP(15) + (-t152 * MDP(14) + (-qJD(1) * t104 + t95) * MDP(15)) * qJD(4) + t142 * t171) * t134 + (-0.2e1 * MDP(9) * t160 + t111 * t158 + t141 * MDP(14) + (-qJD(1) * t98 + qJDD(1) * t104 - t93) * MDP(15) + (-0.2e1 * MDP(8) * t163 + t152 * MDP(13) + (qJD(1) * t103 - t97) * MDP(15) + t102 * t171) * qJD(4)) * t136; t108 * MDP(13) - t109 * MDP(14) + (-qJD(4) * t149 - t134 * t93 + t136 * t94 - g(3)) * MDP(16) + (MDP(4) + MDP(7)) * t168; -t139 * MDP(6) + (qJDD(3) + t143) * MDP(7) + t177 * MDP(16) + (t158 + t165 * MDP(14) + (qJD(4) * t97 + t93) * MDP(16)) * t136 + (t165 * MDP(13) - t157 + (-qJD(4) * t95 + t94) * MDP(16)) * t134 + (MDP(5) + t151 + t154) * qJDD(1); MDP(10) * t159 - MDP(11) * t160 + qJDD(4) * MDP(12) + (-t134 * t168 + t136 * t143 + t100) * MDP(13) + (t101 * qJD(4) + (-qJD(4) * t106 - t168) * t136 + (-t105 - t143) * t134) * MDP(14) + (-pkin(4) * t159 + (t170 - t172) * t163) * MDP(15) + (t172 * t97 + (g(3) * t134 + t177 * t136 + t93) * pkin(4)) * MDP(16) + (MDP(8) * t134 * t136 - t178) * t139; (qJD(1) * t149 + t142) * MDP(16) + t139 * t154;];
tau = t1;
