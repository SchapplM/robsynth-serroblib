% Calculate vector of inverse dynamics joint torques for
% S4PRRP5
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
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:24
% EndTime: 2019-12-31 16:29:25
% DurationCPUTime: 0.67s
% Computational Cost: add. (299->124), mult. (660->167), div. (0->0), fcn. (378->6), ass. (0->60)
t135 = qJ(4) + pkin(5);
t120 = (qJD(2) * qJD(3));
t142 = -2 * t120;
t97 = sin(qJ(3));
t92 = t97 ^ 2;
t99 = cos(qJ(3));
t93 = t99 ^ 2;
t132 = t92 + t93;
t141 = (t92 - t93) * MDP(6);
t94 = sin(pkin(6));
t95 = cos(pkin(6));
t114 = g(1) * t95 + g(2) * t94;
t140 = t132 * MDP(12);
t100 = cos(qJ(2));
t101 = qJD(3) ^ 2;
t119 = t100 * qJDD(1);
t121 = qJD(1) * qJD(2);
t98 = sin(qJ(2));
t90 = t98 * t121;
t139 = (2 * qJDD(2) * pkin(2)) - pkin(5) * t101 - g(3) * t100 + (t114 + t121) * t98 + t119 - t90;
t130 = qJD(3) * pkin(3);
t115 = t98 * qJD(1) + t135 * qJD(2);
t78 = t115 * t97;
t77 = -t78 + t130;
t136 = t77 * t97;
t134 = t78 + t77;
t131 = qJD(2) * pkin(2);
t129 = qJD(2) * t99;
t122 = t100 * qJD(1);
t91 = t99 * pkin(3) + pkin(2);
t82 = -t91 * qJD(2) + qJD(4) - t122;
t127 = t82 * qJD(2);
t126 = qJDD(1) - g(3);
t102 = qJD(2) ^ 2;
t125 = t101 + t102;
t124 = qJDD(2) * t97;
t123 = qJDD(2) * t99;
t118 = t97 * t120;
t117 = MDP(13) * t122;
t116 = qJD(3) * t135;
t113 = g(1) * t94 - g(2) * t95;
t112 = -MDP(4) + t140;
t111 = qJD(3) * t115;
t110 = t113 * t99;
t88 = -t122 - t131;
t109 = t122 + t88 - t131;
t84 = qJDD(2) * pkin(5) + t98 * qJDD(1) + t100 * t121;
t108 = (qJ(4) * qJDD(2)) + qJD(2) * qJD(4) + t84;
t105 = pkin(3) * t118 - t91 * qJDD(2) + qJDD(4) + t90;
t104 = g(3) * t98 + t114 * t100;
t103 = -t88 * qJD(2) + t104 - t84;
t86 = t135 * t99;
t85 = t135 * t97;
t81 = -t97 * qJD(4) - t99 * t116;
t80 = t99 * qJD(4) - t97 * t116;
t79 = t115 * t99;
t76 = t105 - t119;
t75 = t108 * t99 - t97 * t111;
t74 = qJDD(3) * pkin(3) - t108 * t97 - t99 * t111;
t1 = [t126 * MDP(1) - g(3) * MDP(13) + ((qJDD(2) * MDP(3)) + (-0.2e1 * t118 + t123) * MDP(10) + (t99 * t142 - t124) * MDP(11) + (-qJD(2) * t136 + t79 * t129 - t76) * MDP(13) + t112 * t102) * t100 + (MDP(13) * t127 - t102 * MDP(3) + (-t125 * MDP(10) - qJDD(3) * MDP(11) + (-qJD(3) * t77 + t75) * MDP(13)) * t99 + (-qJDD(3) * MDP(10) + t125 * MDP(11) + (-qJD(3) * t79 - t74) * MDP(13)) * t97 + t112 * qJDD(2)) * t98; (-t74 * t85 + t75 * t86 - t76 * t91 + t77 * t81 + t79 * t80) * MDP(13) + (t92 * MDP(5) + MDP(2)) * qJDD(2) + t141 * t142 + (t114 * MDP(3) - t126 * MDP(4) - g(3) * MDP(12) + (-g(3) * t135 - t82 * qJD(1) + t114 * t91) * MDP(13)) * t98 + (t126 * MDP(3) + t114 * MDP(4) + (-t132 * t121 - t114) * MDP(12) + (-g(3) * t91 - t114 * t135) * MDP(13)) * t100 + (t101 * MDP(7) + (qJD(2) * t80 + qJDD(2) * t86 + t75) * MDP(12) - t79 * t117 + (-pkin(5) * MDP(11) + MDP(8)) * qJDD(3) + (t109 * MDP(11) + (qJD(2) * t85 - t77) * MDP(12)) * qJD(3) + t139 * MDP(10)) * t99 + (0.2e1 * MDP(6) * t123 - t101 * MDP(8) + (-qJD(2) * t81 + qJDD(2) * t85 - t74) * MDP(12) + t77 * t117 + (-pkin(5) * MDP(10) + MDP(7)) * qJDD(3) + (0.2e1 * MDP(5) * t129 + t109 * MDP(10) + (-qJD(2) * t86 - t79) * MDP(12) + pkin(3) * t82 * MDP(13)) * qJD(3) - t139 * MDP(11)) * t97; MDP(7) * t124 + MDP(8) * t123 + qJDD(3) * MDP(9) + (t103 * t97 - t110) * MDP(10) + (t103 * t99 + t113 * t97) * MDP(11) + (-pkin(3) * t124 + (-t130 + t134) * t129) * MDP(12) + (t134 * t79 + (t74 - t110 + (t104 - t127) * t97) * pkin(3)) * MDP(13) + (-t97 * t99 * MDP(5) + t141) * t102; (-t114 * t98 - t126 * t100 + (-t79 * t99 + t136) * qJD(2) + t105) * MDP(13) - t102 * t140;];
tau = t1;
