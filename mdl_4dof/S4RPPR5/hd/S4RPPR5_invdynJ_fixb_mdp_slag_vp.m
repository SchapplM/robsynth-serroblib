% Calculate vector of inverse dynamics joint torques for
% S4RPPR5
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:50
% EndTime: 2019-12-31 16:39:51
% DurationCPUTime: 0.41s
% Computational Cost: add. (255->109), mult. (438->139), div. (0->0), fcn. (243->6), ass. (0->54)
t114 = -pkin(1) - pkin(2);
t108 = sin(pkin(6));
t109 = cos(pkin(6));
t132 = qJ(2) * qJDD(1);
t92 = t114 * qJDD(1) + qJDD(2);
t142 = t108 * t92 + t109 * t132;
t93 = t114 * qJD(1) + qJD(2);
t141 = t109 * t93;
t88 = t109 * qJ(2) + t108 * t114;
t140 = (pkin(1) * qJDD(1));
t139 = t108 * qJ(2);
t112 = cos(qJ(4));
t116 = qJD(1) ^ 2;
t138 = t112 * t116;
t137 = qJDD(3) + g(3);
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t136 = t113 * pkin(1) + t111 * qJ(2);
t135 = g(1) * t111 - g(2) * t113;
t110 = sin(qJ(4));
t106 = t110 ^ 2;
t134 = -t112 ^ 2 + t106;
t133 = qJ(2) * qJD(1);
t131 = qJD(1) * qJD(2);
t130 = qJD(1) * qJD(4);
t129 = qJDD(1) * t112;
t96 = t109 * t131;
t77 = t96 + t142;
t128 = 0.2e1 * t131;
t127 = 0.2e1 * t130;
t126 = t108 * t131;
t125 = -t108 * t132 + t109 * t92;
t124 = t112 * t127;
t123 = qJDD(2) - t140;
t82 = t113 * t108 - t111 * t109;
t83 = t111 * t108 + t113 * t109;
t122 = g(1) * t83 - g(2) * t82;
t121 = -g(1) * t82 - g(2) * t83;
t120 = g(1) * t113 + g(2) * t111;
t87 = t109 * t114 - t139;
t115 = qJD(4) ^ 2;
t91 = qJDD(4) * t112 - t115 * t110;
t90 = -qJDD(4) * t110 - t115 * t112;
t76 = t125 - t126;
t78 = -t141 + (pkin(3) + t139) * qJD(1);
t119 = qJDD(1) * pkin(5) + t78 * qJD(1) + t122 - t77;
t84 = pkin(3) - t87;
t85 = -pkin(5) + t88;
t118 = -qJDD(4) * t85 + (-qJD(1) * t84 - qJD(2) * t109 - t78) * qJD(4);
t117 = -t115 * t85 + t121 + t126 - t76 + (pkin(3) + t84) * qJDD(1);
t102 = t113 * qJ(2);
t81 = t108 * t93 + t109 * t133;
t80 = -t108 * t133 + t141;
t1 = [qJDD(1) * MDP(1) + t135 * MDP(2) + t120 * MDP(3) + (-qJDD(2) + t135 + (2 * t140)) * MDP(4) + (-t120 + t128 + 0.2e1 * t132) * MDP(5) + (-t123 * pkin(1) - g(1) * (-t111 * pkin(1) + t102) - g(2) * t136 + (t128 + t132) * qJ(2)) * MDP(6) + (-t87 * qJDD(1) + t121 - t125 + 0.2e1 * t126) * MDP(7) + (t88 * qJDD(1) - t122 + t142 + 0.2e1 * t96) * MDP(8) + (t77 * t88 + t76 * t87 - g(1) * (t114 * t111 + t102) - g(2) * (t113 * pkin(2) + t136) + (-t80 * t108 + t81 * t109) * qJD(2)) * MDP(9) + (t106 * qJDD(1) + t110 * t124) * MDP(10) + 0.2e1 * (t110 * t129 - t134 * t130) * MDP(11) + t90 * MDP(12) - t91 * MDP(13) + (t110 * t118 + t112 * t117) * MDP(15) + (-t110 * t117 + t112 * t118) * MDP(16); -qJDD(1) * MDP(4) - t116 * MDP(5) + (-t116 * qJ(2) + t123 - t135) * MDP(6) - t135 * MDP(9) + (-qJDD(1) * MDP(7) - t116 * MDP(8) + (-qJD(1) * t81 + t76) * MDP(9) + (t110 * t127 - t129) * MDP(15) + (qJDD(1) * t110 + t124) * MDP(16)) * t109 + (-t116 * MDP(7) + qJDD(1) * MDP(8) + (qJD(1) * t80 + t77) * MDP(9) + (t90 - t138) * MDP(15) + (t110 * t116 - t91) * MDP(16)) * t108; t91 * MDP(15) + t90 * MDP(16) + t137 * MDP(9); qJDD(4) * MDP(14) + t134 * MDP(11) * t116 + (-qJDD(1) * MDP(13) + t137 * MDP(15) + t119 * MDP(16)) * t112 + (-MDP(10) * t138 - qJDD(1) * MDP(12) + t119 * MDP(15) - t137 * MDP(16)) * t110;];
tau = t1;
