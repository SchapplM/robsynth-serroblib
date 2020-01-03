% Calculate vector of inverse dynamics joint torques for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:00
% EndTime: 2019-12-31 16:28:01
% DurationCPUTime: 0.58s
% Computational Cost: add. (323->116), mult. (599->149), div. (0->0), fcn. (296->4), ass. (0->56)
t101 = cos(qJ(3));
t100 = sin(qJ(3));
t129 = qJD(2) * t100;
t121 = pkin(5) * t129;
t80 = qJD(1) * t101 - t121;
t142 = qJD(4) - t80;
t97 = pkin(6) + qJ(2);
t93 = sin(t97);
t94 = cos(t97);
t115 = g(1) * t94 + g(2) * t93;
t113 = pkin(3) * t101 + qJ(4) * t100;
t110 = pkin(2) + t113;
t77 = t110 * qJD(2);
t112 = pkin(3) * t100 - qJ(4) * t101;
t75 = t112 * qJD(3) - qJD(4) * t100;
t72 = qJD(2) * t75 - t110 * qJDD(2);
t76 = -qJD(3) * pkin(3) + t142;
t81 = pkin(5) * qJD(2) * t101 + qJD(1) * t100;
t78 = qJD(3) * qJ(4) + t81;
t135 = g(2) * t94;
t138 = g(1) * t93;
t141 = -t135 + t138;
t131 = qJDD(3) * pkin(3);
t140 = qJDD(4) - t131;
t133 = pkin(5) * qJDD(3);
t139 = -0.2e1 * t77 * qJD(3) - t133;
t98 = t100 ^ 2;
t99 = t101 ^ 2;
t134 = t98 - t99;
t103 = qJD(2) ^ 2;
t130 = t100 * t103;
t127 = qJD(1) * qJD(3);
t126 = qJD(2) * qJD(3);
t125 = qJDD(2) * t100;
t124 = qJDD(2) * t101;
t123 = qJDD(3) * qJ(4);
t122 = pkin(5) * t124 + t100 * qJDD(1) + t101 * t127;
t120 = t101 * t130;
t119 = t101 * t126;
t118 = -t101 * qJDD(1) + t100 * t127 + (t119 + t125) * pkin(5);
t102 = qJD(3) ^ 2;
t114 = pkin(5) * t102 + t135;
t111 = g(3) * t100 - t122;
t109 = -0.2e1 * pkin(2) * t126 - t133;
t108 = 0.2e1 * qJDD(2) * pkin(2) - t114;
t107 = -g(3) * t101 + t115 * t100 - t118;
t106 = qJD(3) * t81 + t107;
t105 = -t114 - 0.2e1 * t72;
t73 = t123 + (qJD(4) - t121) * qJD(3) + t122;
t74 = t118 + t140;
t104 = t74 * t100 + t73 * t101 + (-t100 * t78 + t101 * t76) * qJD(3) - t115;
t86 = t101 * t138;
t84 = qJDD(3) * t101 - t100 * t102;
t83 = qJDD(3) * t100 + t101 * t102;
t79 = t112 * qJD(2);
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t100 * t73 - t101 * t74 - g(3) + (t100 * t76 + t101 * t78) * qJD(3)) * MDP(15) + (MDP(10) + MDP(12)) * t84 + (-MDP(11) + MDP(14)) * t83; qJDD(2) * MDP(2) + t141 * MDP(3) + t115 * MDP(4) + (qJDD(2) * t98 + 0.2e1 * t100 * t119) * MDP(5) + 0.2e1 * (t100 * t124 - t134 * t126) * MDP(6) + t83 * MDP(7) + t84 * MDP(8) + (t109 * t100 + t108 * t101 + t86) * MDP(10) + (t109 * t101 + (-t108 - t138) * t100) * MDP(11) + (t100 * t139 + t105 * t101 + t86) * MDP(12) + ((t98 + t99) * qJDD(2) * pkin(5) + t104) * MDP(13) + (-t139 * t101 + (t105 + t138) * t100) * MDP(14) + (t104 * pkin(5) - t77 * t75 + (t141 - t72) * t110) * MDP(15); -MDP(5) * t120 + t134 * MDP(6) * t103 + MDP(7) * t125 + MDP(8) * t124 + qJDD(3) * MDP(9) + (pkin(2) * t130 + t106) * MDP(10) + ((t80 + t121) * qJD(3) + (pkin(2) * t103 + t115) * t101 + t111) * MDP(11) + (0.2e1 * t131 - qJDD(4) + (t100 * t77 + t101 * t79) * qJD(2) + t106) * MDP(12) - t112 * qJDD(2) * MDP(13) + (0.2e1 * t123 - t115 * t101 + (0.2e1 * qJD(4) - t80) * qJD(3) + (-t101 * t77 + (-pkin(5) * qJD(3) + t79) * t100) * qJD(2) - t111) * MDP(14) + (-t74 * pkin(3) - g(3) * t113 + t73 * qJ(4) + t115 * t112 + t142 * t78 - t76 * t81 + t77 * t79) * MDP(15); (-qJDD(3) - t120) * MDP(12) + MDP(13) * t125 + (-t103 * t98 - t102) * MDP(14) + (-qJD(3) * t78 - t77 * t129 - t107 + t140) * MDP(15);];
tau = t1;
