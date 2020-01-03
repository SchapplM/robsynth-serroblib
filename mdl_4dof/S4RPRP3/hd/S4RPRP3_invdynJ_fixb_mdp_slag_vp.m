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
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
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
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4RPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:50
% EndTime: 2019-12-31 16:42:51
% DurationCPUTime: 0.47s
% Computational Cost: add. (307->103), mult. (603->145), div. (0->0), fcn. (309->8), ass. (0->57)
t90 = sin(pkin(6));
t80 = pkin(1) * t90 + pkin(5);
t118 = qJ(4) + t80;
t91 = cos(pkin(6));
t123 = pkin(1) * t91;
t95 = cos(qJ(3));
t82 = pkin(3) * t95 + pkin(2);
t105 = -t82 - t123;
t113 = qJD(1) * qJD(3);
t93 = sin(qJ(3));
t112 = t93 * t113;
t124 = pkin(3) * t112 + t105 * qJDD(1) + qJDD(4);
t122 = g(3) * t95;
t119 = qJD(3) * pkin(3);
t110 = t118 * qJD(1);
t66 = t95 * qJD(2) - t110 * t93;
t65 = t66 + t119;
t121 = t65 - t66;
t88 = t93 ^ 2;
t89 = t95 ^ 2;
t120 = t88 - t89;
t81 = -pkin(2) - t123;
t78 = qJD(1) * t81;
t117 = qJDD(2) - g(3);
t116 = qJDD(1) * t93;
t115 = qJDD(1) * t95;
t111 = qJD(3) * t118;
t87 = qJ(1) + pkin(6);
t83 = sin(t87);
t84 = cos(t87);
t109 = g(1) * t84 + g(2) * t83;
t108 = -g(1) * t83 + g(2) * t84;
t94 = sin(qJ(1));
t96 = cos(qJ(1));
t107 = g(1) * t94 - g(2) * t96;
t67 = qJD(2) * t93 + t110 * t95;
t106 = t65 * t93 - t67 * t95;
t104 = t110 * qJD(3);
t75 = t80 * qJDD(1);
t103 = -qJD(1) * t78 + t109 - t75;
t102 = 0.2e1 * qJD(3) * t78 - qJDD(3) * t80;
t101 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t75;
t97 = qJD(3) ^ 2;
t100 = -0.2e1 * qJDD(1) * t81 - t80 * t97 - t108;
t98 = qJD(1) ^ 2;
t92 = -qJ(4) - pkin(5);
t85 = t95 * qJDD(2);
t74 = qJDD(3) * t95 - t93 * t97;
t73 = qJDD(3) * t93 + t95 * t97;
t72 = t118 * t95;
t71 = t118 * t93;
t70 = t105 * qJD(1) + qJD(4);
t69 = -qJD(4) * t93 - t95 * t111;
t68 = qJD(4) * t95 - t93 * t111;
t64 = (qJDD(2) - t104) * t93 + t101 * t95;
t63 = qJDD(3) * pkin(3) - t101 * t93 - t95 * t104 + t85;
t1 = [qJDD(1) * MDP(1) + t107 * MDP(2) + (g(1) * t96 + g(2) * t94) * MDP(3) + (t107 + (t90 ^ 2 + t91 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t88 + 0.2e1 * t95 * t112) * MDP(5) + 0.2e1 * (-t120 * t113 + t93 * t115) * MDP(6) + t73 * MDP(7) + t74 * MDP(8) + (t100 * t95 + t102 * t93) * MDP(10) + (-t100 * t93 + t102 * t95) * MDP(11) + ((-qJD(3) * t65 + qJDD(1) * t72 + t64 + (qJD(3) * t71 + t68) * qJD(1)) * t95 + (-qJD(3) * t67 + qJDD(1) * t71 - t63 + (-qJD(3) * t72 - t69) * qJD(1)) * t93 - t109) * MDP(12) + (t64 * t72 + t67 * t68 - t63 * t71 + t65 * t69 + t70 * t93 * t119 - g(1) * (-pkin(1) * t94 - t82 * t83 - t84 * t92) - g(2) * (pkin(1) * t96 + t82 * t84 - t83 * t92) + t124 * t105) * MDP(13); t117 * MDP(4) + t74 * MDP(10) - t73 * MDP(11) + (-t106 * qJD(3) + t63 * t95 + t64 * t93 - g(3)) * MDP(13); MDP(7) * t116 + MDP(8) * t115 + qJDD(3) * MDP(9) + (t103 * t93 - t122 + t85) * MDP(10) + (t103 * t95 - t117 * t93) * MDP(11) + (-pkin(3) * t116 + (-t119 + t121) * t95 * qJD(1)) * MDP(12) + (t121 * t67 + (-t122 + t63 + (-qJD(1) * t70 + t109) * t93) * pkin(3)) * MDP(13) + (-t93 * t95 * MDP(5) + t120 * MDP(6)) * t98; (-t88 - t89) * MDP(12) * t98 + (t106 * qJD(1) + t108 + t124) * MDP(13);];
tau = t1;
