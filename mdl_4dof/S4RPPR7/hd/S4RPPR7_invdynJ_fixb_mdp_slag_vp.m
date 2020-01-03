% Calculate vector of inverse dynamics joint torques for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:43
% EndTime: 2019-12-31 16:41:44
% DurationCPUTime: 0.70s
% Computational Cost: add. (329->131), mult. (607->166), div. (0->0), fcn. (372->8), ass. (0->67)
t127 = -pkin(1) - qJ(3);
t125 = sin(pkin(6));
t126 = cos(pkin(6));
t154 = t125 ^ 2 + t126 ^ 2;
t170 = (qJD(1) * t127 + qJD(2)) * t154;
t122 = qJDD(1) * qJ(2);
t123 = qJD(1) * qJD(2);
t168 = t122 + t123;
t104 = qJDD(3) + t168;
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t141 = g(1) * t131 + g(2) * t129;
t169 = t104 - t141;
t155 = g(1) * t129 - g(2) * t131;
t167 = t125 * MDP(7) + t126 * MDP(8);
t166 = -qJD(1) * qJD(3) + qJDD(1) * t127;
t165 = 0.2e1 * t123;
t164 = pkin(3) * t125;
t163 = -pkin(5) + t127;
t130 = cos(qJ(4));
t157 = t126 * t130;
t128 = sin(qJ(4));
t159 = t125 * t128;
t137 = -t157 + t159;
t98 = t125 * t130 + t126 * t128;
t94 = t98 * qJD(4);
t162 = -t94 * qJD(4) - qJDD(4) * t137;
t161 = pkin(1) * qJDD(1);
t156 = t131 * pkin(1) + t129 * qJ(2);
t153 = qJD(1) * t128;
t152 = qJD(1) * t130;
t111 = qJD(1) * qJ(2) + qJD(3);
t150 = qJDD(1) * t128;
t149 = qJDD(1) * t130;
t148 = t126 * t152;
t147 = qJD(4) * t157;
t146 = t125 * t153;
t145 = qJDD(2) - t155;
t144 = t154 * MDP(9);
t100 = qJDD(2) + t166;
t143 = t154 * t100;
t142 = -pkin(5) * qJDD(1) + t100;
t107 = t126 * t149;
t140 = -t125 * t150 + t107;
t101 = t163 * t125;
t102 = t163 * t126;
t139 = t101 * t130 + t102 * t128;
t138 = t101 * t128 - t102 * t130;
t95 = -qJD(4) * t159 + t147;
t136 = -qJD(4) * t95 - qJDD(4) * t98;
t134 = -t125 * t152 - t126 * t153;
t132 = qJD(1) ^ 2;
t121 = pkin(6) + qJ(4);
t115 = t131 * qJ(2);
t113 = cos(t121);
t112 = sin(t121);
t110 = qJ(2) + t164;
t105 = qJD(4) * t146;
t103 = qJD(1) * t164 + t111;
t97 = qJDD(1) * t164 + t104;
t93 = -t146 + t148;
t91 = t98 * qJD(1);
t87 = t142 * t126;
t86 = t142 * t125;
t85 = qJD(1) * t147 + qJDD(1) * t98 - t105;
t84 = -qJD(1) * t94 + t140;
t1 = [qJDD(1) * MDP(1) + t155 * MDP(2) + t141 * MDP(3) + (t145 - 0.2e1 * t161) * MDP(4) + (0.2e1 * t122 + t165 - t141) * MDP(5) + (-(qJDD(2) - t161) * pkin(1) - g(1) * (-pkin(1) * t129 + t115) - g(2) * t156 + (t122 + t165) * qJ(2)) * MDP(6) + (t155 + t154 * (-t100 - t166)) * MDP(9) + (t104 * qJ(2) + t111 * qJD(2) - g(1) * (t127 * t129 + t115) - g(2) * (qJ(3) * t131 + t156) + t127 * t143 - qJD(3) * t170) * MDP(10) + (-t137 * t84 - t93 * t94) * MDP(11) + (t137 * t85 - t84 * t98 + t91 * t94 - t93 * t95) * MDP(12) + t162 * MDP(13) + t136 * MDP(14) + (qJD(2) * t91 + t110 * t85 + t97 * t98 + t103 * t95 - t138 * qJDD(4) - t141 * t112 + (t137 * qJD(3) - t139 * qJD(4)) * qJD(4)) * MDP(16) + (qJD(2) * t93 + t110 * t84 - t97 * t137 - t103 * t94 - t139 * qJDD(4) - t141 * t113 + (qJD(3) * t98 + qJD(4) * t138) * qJD(4)) * MDP(17) + t167 * (t168 + t169); t145 * MDP(6) + (-qJD(1) * t111 + t143 - t155) * MDP(10) + (-qJD(1) * t91 + t162) * MDP(16) + (-qJD(1) * t93 + t136) * MDP(17) + (-MDP(6) * qJ(2) - MDP(5) - t167) * t132 + (-pkin(1) * MDP(6) + MDP(4) - t144) * qJDD(1); (qJD(1) * t170 + t169) * MDP(10) - t105 * MDP(16) + t107 * MDP(17) - t132 * t144 + ((t128 * MDP(16) + MDP(8)) * t126 + (t130 * MDP(16) - t128 * MDP(17) + MDP(7)) * t125) * qJDD(1) + ((t93 + t148) * MDP(16) + (t134 - t91) * MDP(17)) * qJD(4); t93 * t91 * MDP(11) + (-t91 ^ 2 + t93 ^ 2) * MDP(12) + t140 * MDP(13) + (-t125 * t149 - t126 * t150 + t105) * MDP(14) + qJDD(4) * MDP(15) + (g(3) * t112 - t103 * t93 - t155 * t113 - t128 * t86 + t130 * t87) * MDP(16) + (g(3) * t113 + t103 * t91 + t155 * t112 - t128 * t87 - t130 * t86) * MDP(17) + ((t134 + t91) * MDP(13) + (t93 - t148) * MDP(14)) * qJD(4);];
tau = t1;
