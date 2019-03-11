% Calculate vector of inverse dynamics joint torques for
% S4RPPR2
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
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:33
% EndTime: 2019-03-08 18:28:34
% DurationCPUTime: 0.46s
% Computational Cost: add. (292->106), mult. (417->133), div. (0->0), fcn. (238->8), ass. (0->58)
t137 = qJD(1) - qJD(4);
t105 = -pkin(1) - pkin(2);
t100 = cos(pkin(6));
t127 = t100 * qJDD(1);
t85 = qJDD(1) * t105 + qJDD(2);
t99 = sin(pkin(6));
t136 = qJ(2) * t127 + t85 * t99;
t102 = sin(qJ(1));
t104 = cos(qJ(1));
t135 = pkin(1) * t104 + qJ(2) * t102;
t134 = g(1) * t102 - g(2) * t104;
t133 = t99 * qJ(2);
t132 = pkin(1) * qJDD(1);
t97 = qJDD(1) - qJDD(4);
t131 = t97 * MDP(10);
t130 = qJ(2) * qJD(1);
t129 = qJ(2) * qJDD(1);
t128 = qJD(1) * qJD(2);
t126 = pkin(6) + qJ(4);
t125 = 0.2e1 * t128;
t124 = t99 * t128;
t123 = -pkin(3) - t133;
t122 = t100 * t128;
t82 = t100 * t105 - t133;
t121 = cos(t126);
t120 = sin(t126);
t119 = qJDD(2) - t132;
t118 = g(1) * t104 + g(2) * t102;
t89 = qJD(1) * t105 + qJD(2);
t74 = t100 * t130 + t89 * t99;
t84 = t100 * t89;
t117 = t74 * t100 - (-t130 * t99 + t84) * t99;
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t72 = qJD(1) * t123 + t84;
t116 = t101 * t74 - t103 * t72;
t115 = -t101 * t72 - t103 * t74;
t79 = -pkin(3) + t82;
t83 = qJ(2) * t100 + t105 * t99;
t114 = t101 * t83 - t103 * t79;
t113 = t101 * t79 + t103 * t83;
t112 = t100 * t103 - t101 * t99;
t111 = t100 * t101 + t103 * t99;
t110 = t112 * t137;
t109 = t111 * t137;
t81 = t100 * t85;
t69 = qJDD(1) * t123 - t124 + t81;
t71 = t122 + t136;
t75 = -t102 * t120 - t104 * t121;
t76 = -t102 * t121 + t104 * t120;
t108 = g(1) * t76 - g(2) * t75 - t101 * t71 + t103 * t69;
t107 = -g(1) * t75 - g(2) * t76 - t101 * t69 - t103 * t71;
t106 = qJD(1) ^ 2;
t93 = t104 * qJ(2);
t78 = t100 * t104 + t102 * t99;
t77 = -t100 * t102 + t104 * t99;
t70 = t81 + (-t128 - t129) * t99;
t1 = [qJDD(1) * MDP(1) + t134 * MDP(2) + t118 * MDP(3) + (-qJDD(2) + 0.2e1 * t132 + t134) * MDP(4) + (-t118 + t125 + 0.2e1 * t129) * MDP(5) + (-t119 * pkin(1) - g(1) * (-pkin(1) * t102 + t93) - g(2) * t135 + (t125 + t129) * qJ(2)) * MDP(6) + (0.2e1 * t124 - g(1) * t77 - g(2) * t78 - t81 + (-t82 + t133) * qJDD(1)) * MDP(7) + (-g(1) * t78 + g(2) * t77 + qJDD(1) * t83 + 0.2e1 * t122 + t136) * MDP(8) + (t71 * t83 + t70 * t82 - g(1) * (t102 * t105 + t93) - g(2) * (pkin(2) * t104 + t135) + t117 * qJD(2)) * MDP(9) + t131 + (t114 * t97 + qJD(2) * t109 + (t113 * t137 - t115) * qJD(4) - t108) * MDP(11) + (t113 * t97 + qJD(2) * t110 + (-t114 * t137 - t116) * qJD(4) - t107) * MDP(12); -qJDD(1) * MDP(4) - t106 * MDP(5) + (-qJ(2) * t106 + t119 - t134) * MDP(6) + (-t106 * t99 - t127) * MDP(7) + (qJDD(1) * t99 - t100 * t106) * MDP(8) + (-qJD(1) * t117 + t70 * t100 + t71 * t99 - t134) * MDP(9) + (-t109 * t137 - t112 * t97) * MDP(11) + (-t110 * t137 + t111 * t97) * MDP(12); (qJDD(3) + g(3)) * MDP(9); -t131 + (t115 * t137 + t108) * MDP(11) + (t116 * t137 + t107) * MDP(12) + (t115 * MDP(11) + MDP(12) * t116) * qJD(4);];
tau  = t1;
