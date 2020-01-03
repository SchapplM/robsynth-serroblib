% Calculate vector of inverse dynamics joint torques for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:20
% EndTime: 2019-12-31 16:23:21
% DurationCPUTime: 0.37s
% Computational Cost: add. (213->84), mult. (452->127), div. (0->0), fcn. (328->10), ass. (0->50)
t101 = sin(qJ(2));
t103 = cos(qJ(2));
t117 = qJD(1) * qJD(2);
t128 = qJDD(1) * t101 + t103 * t117;
t102 = cos(qJ(4));
t100 = sin(qJ(4));
t94 = t100 ^ 2;
t127 = (-t102 ^ 2 + t94) * MDP(7);
t97 = sin(pkin(6));
t99 = cos(pkin(6));
t113 = g(1) * t99 + g(2) * t97;
t92 = t103 * qJDD(1);
t78 = qJDD(2) * pkin(2) - t101 * t117 + t92;
t96 = sin(pkin(7));
t98 = cos(pkin(7));
t67 = -t128 * t96 + t78 * t98;
t122 = t102 * MDP(6);
t121 = qJDD(1) - g(3);
t120 = qJD(1) * t101;
t119 = qJD(1) * t103;
t88 = pkin(2) * t96 + pkin(5);
t118 = qJDD(4) * t88;
t116 = qJD(2) * qJD(4);
t68 = t128 * t98 + t96 * t78;
t85 = qJD(2) * pkin(2) + t119;
t71 = -t96 * t120 + t85 * t98;
t69 = -qJD(2) * pkin(3) - t71;
t79 = t101 * t96 - t103 * t98;
t77 = t79 * qJD(1);
t89 = -pkin(2) * t98 - pkin(3);
t112 = qJD(2) * t89 + t69 - t77;
t80 = t101 * t98 + t103 * t96;
t111 = g(1) * t97 - g(2) * t99 - qJDD(3);
t104 = qJD(4) ^ 2;
t74 = t80 * qJD(2);
t109 = qJD(2) * t74 + qJDD(2) * t79 + t104 * t80;
t76 = qJD(2) * t79;
t108 = qJD(4) * t76 - qJDD(4) * t80 + t79 * t116;
t93 = qJ(2) + pkin(7);
t90 = sin(t93);
t91 = cos(t93);
t107 = -qJDD(2) * pkin(5) + g(3) * t90 - t69 * qJD(2) + t113 * t91 - t68;
t87 = t98 * t120;
t75 = t96 * t119 + t87;
t106 = -g(3) * t91 + qJD(2) * t75 - t104 * t88 + t113 * t90 + t67 + (pkin(3) - t89) * qJDD(2);
t105 = qJD(2) ^ 2;
t83 = qJDD(4) * t102 - t100 * t104;
t82 = qJDD(4) * t100 + t102 * t104;
t72 = t96 * t85 + t87;
t1 = [t121 * MDP(1) + (qJDD(2) * t103 - t101 * t105) * MDP(3) + (-qJDD(2) * t101 - t103 * t105) * MDP(4) + (-t67 * t79 + t68 * t80 - t71 * t74 - t72 * t76 - g(3)) * MDP(5) + (-t109 * MDP(11) + t108 * MDP(12)) * t102 + (t108 * MDP(11) + t109 * MDP(12)) * t100; t92 * MDP(3) + t82 * MDP(8) + t83 * MDP(9) + (-g(3) * MDP(3) + t113 * MDP(4)) * t103 + (t113 * MDP(3) - t121 * MDP(4)) * t101 + (t94 * MDP(6) + MDP(2)) * qJDD(2) - 0.2e1 * t116 * t127 + (t106 * MDP(11) + (t112 * qJD(4) - t118) * MDP(12)) * t102 + (0.2e1 * qJDD(2) * t102 * MDP(7) - MDP(11) * t118 - t106 * MDP(12) + (t112 * MDP(11) + 0.2e1 * qJD(2) * t122) * qJD(4)) * t100 + (t71 * t75 + t72 * t77 + (-g(3) * t103 + t113 * t101 + t67 * t98 + t68 * t96) * pkin(2)) * MDP(5); t83 * MDP(11) - t82 * MDP(12) - t111 * MDP(5); qJDD(4) * MDP(10) + t105 * t127 + (-t111 * MDP(11) + t107 * MDP(12) + qJDD(2) * MDP(9)) * t102 + (t107 * MDP(11) + t111 * MDP(12) + qJDD(2) * MDP(8) - t105 * t122) * t100;];
tau = t1;
