% Calculate vector of inverse dynamics joint torques for
% S4RRPR1
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x6] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:47
% EndTime: 2019-01-31 13:16:49
% DurationCPUTime: 0.41s
% Computational Cost: add. (402->108), mult. (683->145), div. (0->0), fcn. (386->12), ass. (0->59)
t115 = sin(qJ(2));
t132 = qJDD(1) * t115;
t118 = cos(qJ(2));
t134 = qJD(1) * t118;
t144 = pkin(1) * (qJD(2) * t134 + t132);
t111 = qJ(1) + qJ(2);
t107 = sin(t111);
t108 = cos(t111);
t143 = g(1) * t107 - g(2) * t108;
t110 = qJD(1) + qJD(2);
t106 = qJD(4) + t110;
t142 = qJD(4) - t106;
t141 = pkin(1) * t118;
t112 = sin(pkin(7));
t140 = pkin(2) * t112;
t138 = g(1) * t108 + g(2) * t107;
t109 = qJDD(1) + qJDD(2);
t105 = qJDD(4) + t109;
t96 = t105 * MDP(8);
t137 = t109 * MDP(4) + t96;
t136 = t112 * t115;
t113 = cos(pkin(7));
t135 = t113 * t115;
t114 = sin(qJ(4));
t130 = pkin(1) * qJD(1) * t115;
t89 = pkin(1) * t134 + pkin(2) * t110;
t79 = t112 * t89 + t113 * t130;
t101 = pkin(7) + qJ(4) + t111;
t94 = sin(t101);
t95 = cos(t101);
t131 = qJD(4) * t114 * t79 + g(1) * t95 + g(2) * t94;
t129 = qJD(2) * (-qJD(1) - t110);
t102 = pkin(2) + t141;
t128 = -pkin(1) * t136 + t113 * t102;
t103 = qJDD(1) * t141;
t127 = t103 + t143;
t117 = cos(qJ(4));
t78 = -t112 * t130 + t113 * t89;
t76 = pkin(3) * t110 + t78;
t126 = -t114 * t76 - t117 * t79;
t97 = pkin(2) * t113 + pkin(3);
t125 = t114 * t97 + t117 * t140;
t124 = -t114 * t140 + t117 * t97;
t81 = pkin(2) * t109 - qJD(2) * t130 + t103;
t74 = -t112 * t144 + t113 * t81;
t73 = pkin(3) * t109 + t74;
t75 = t112 * t81 + t113 * t144;
t123 = g(1) * t94 - g(2) * t95 - t114 * t75 + t117 * t73;
t122 = pkin(1) * (-t112 * t118 - t135);
t121 = pkin(1) * (t113 * t118 - t136);
t119 = cos(qJ(1));
t116 = sin(qJ(1));
t87 = pkin(1) * t135 + t102 * t112;
t86 = qJD(2) * t121;
t85 = qJD(1) * t121;
t84 = qJD(2) * t122;
t83 = qJD(1) * t122;
t82 = pkin(3) + t128;
t1 = [qJDD(1) * MDP(1) + (g(1) * t116 - g(2) * t119) * MDP(2) + (g(1) * t119 + g(2) * t116) * MDP(3) + ((t109 * t118 + t115 * t129) * pkin(1) + t127) * MDP(5) + (((-qJDD(1) - t109) * t115 + t118 * t129) * pkin(1) + t138) * MDP(6) + (t75 * t87 + t79 * t86 + t74 * t128 + t78 * t84 - g(1) * (-pkin(1) * t116 - pkin(2) * t107) - g(2) * (pkin(1) * t119 + pkin(2) * t108)) * MDP(7) + ((-t114 * t86 + t117 * t84) * t106 + (-t114 * t87 + t117 * t82) * t105 + ((-t114 * t82 - t117 * t87) * t106 + t126) * qJD(4) + t123) * MDP(9) + ((-(-qJD(4) * t87 + t84) * t106 - t82 * t105 - t73) * t114 + (-(qJD(4) * t82 + t86) * t106 - t87 * t105 - t75 - qJD(4) * t76) * t117 + t131) * MDP(10) + t137; t127 * MDP(5) + t138 * MDP(6) + (-t78 * t83 - t79 * t85 + (t112 * t75 + t113 * t74 + t143) * pkin(2)) * MDP(7) + (t124 * t105 - (-t114 * t85 + t117 * t83) * t106 + t123) * MDP(9) + (-t125 * t105 - t117 * t75 - t114 * t73 + (t114 * t83 + t117 * t85) * t106 + t131) * MDP(10) + (-MDP(6) * t132 + (t115 * MDP(5) + t118 * MDP(6)) * qJD(1) * (-qJD(2) + t110)) * pkin(1) + ((-t125 * t106 + t126) * MDP(9) + (-t124 * t106 - t117 * t76) * MDP(10)) * qJD(4) + t137; (qJDD(3) - g(3)) * MDP(7); t96 + (t126 * t142 + t123) * MDP(9) + ((-t79 * t106 - t73) * t114 + (-t142 * t76 - t75) * t117 + t131) * MDP(10);];
tau  = t1;
