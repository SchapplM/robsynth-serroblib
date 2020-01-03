% Calculate vector of inverse dynamics joint torques for
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:16
% EndTime: 2019-12-31 16:46:17
% DurationCPUTime: 0.51s
% Computational Cost: add. (288->117), mult. (507->148), div. (0->0), fcn. (217->4), ass. (0->57)
t101 = sin(qJ(1));
t103 = cos(qJ(1));
t134 = g(1) * t101 - g(2) * t103;
t100 = sin(qJ(3));
t136 = pkin(3) * t100;
t116 = qJ(2) + t136;
t81 = t116 * qJD(1) + qJD(4);
t140 = (-qJD(1) * t81 - t134) * MDP(15);
t104 = -pkin(1) - pkin(5);
t85 = t104 * qJD(1) + qJD(2);
t139 = -qJ(4) * qJD(1) + t85;
t102 = cos(qJ(3));
t78 = t139 * t102;
t74 = qJD(3) * pkin(3) + t78;
t112 = g(1) * t103 + g(2) * t101;
t96 = qJD(1) * qJD(2);
t119 = 0.2e1 * t96;
t95 = qJDD(1) * qJ(2);
t138 = -t112 + 0.2e1 * t95 + t119;
t135 = t103 * pkin(1) + t101 * qJ(2);
t97 = t100 ^ 2;
t98 = t102 ^ 2;
t133 = t97 - t98;
t131 = pkin(1) * qJDD(1);
t106 = qJD(1) ^ 2;
t130 = qJ(2) * t106;
t127 = qJ(4) - t104;
t83 = t127 * t102;
t128 = qJD(3) * t83;
t105 = qJD(3) ^ 2;
t126 = -t105 - t106;
t124 = qJD(3) * t100;
t123 = qJD(3) * t102;
t122 = qJD(1) * qJD(3);
t121 = qJD(1) * qJD(4);
t120 = qJDD(1) * t102;
t118 = qJDD(2) - t134;
t117 = t102 * t122;
t115 = (-t97 - t98) * MDP(14);
t84 = t104 * qJDD(1) + qJDD(2);
t79 = t102 * t84;
t72 = -t102 * t121 - t85 * t124 + qJDD(3) * pkin(3) + t79 + (t100 * t122 - t120) * qJ(4);
t77 = t139 * t100;
t114 = qJD(3) * t77 + t72;
t73 = t139 * t123 + (-qJ(4) * qJDD(1) - t121 + t84) * t100;
t113 = -qJD(3) * t74 + t73;
t110 = pkin(3) * t117 + qJDD(1) * t136 + qJDD(4) + t95 + t96;
t109 = 0.2e1 * qJ(2) * t122 + qJDD(3) * t104;
t108 = -t134 - t130;
t107 = -t104 * t105 + t138;
t99 = -qJ(4) - pkin(5);
t91 = t103 * qJ(2);
t89 = qJDD(3) * t102;
t82 = t127 * t100;
t76 = -qJD(4) * t100 - t128;
t75 = -qJD(4) * t102 + t127 * t124;
t1 = [qJDD(1) * MDP(1) + t134 * MDP(2) + t112 * MDP(3) + (t118 - 0.2e1 * t131) * MDP(4) + t138 * MDP(5) + (-(qJDD(2) - t131) * pkin(1) - g(1) * (-pkin(1) * t101 + t91) - g(2) * t135 + (t119 + t95) * qJ(2)) * MDP(6) + (qJDD(1) * t98 - 0.2e1 * t100 * t117) * MDP(7) + 0.2e1 * (-t100 * t120 + t133 * t122) * MDP(8) + (-t100 * t105 + t89) * MDP(9) + (-qJDD(3) * t100 - t102 * t105) * MDP(10) + (t107 * t100 + t109 * t102) * MDP(12) + (-t109 * t100 + t107 * t102) * MDP(13) + ((qJDD(1) * t83 + (qJD(3) * t82 - t75) * qJD(1) - t114) * t102 + (qJDD(1) * t82 + (-t76 - t128) * qJD(1) - t113) * t100 + t134) * MDP(14) + (-t73 * t82 + t77 * t76 - t72 * t83 + t74 * t75 + t110 * t116 + t81 * (pkin(3) * t123 + qJD(2)) - g(1) * (t103 * t136 + t91 + (-pkin(1) + t99) * t101) - g(2) * (t101 * t136 - t103 * t99 + t135)) * MDP(15); -t106 * MDP(5) + (t118 - t130) * MDP(6) + t89 * MDP(12) + t140 + (t126 * MDP(13) + t114 * MDP(15)) * t102 + (t126 * MDP(12) - qJDD(3) * MDP(13) + t113 * MDP(15)) * t100 + (-pkin(1) * MDP(6) + MDP(4) + t115) * qJDD(1); qJDD(3) * MDP(11) + t79 * MDP(12) + (pkin(3) * t72 + (t74 - t78) * t77) * MDP(15) - t133 * MDP(8) * t106 + (qJDD(1) * MDP(9) + t108 * MDP(12) + g(3) * MDP(13) + (-qJDD(1) * MDP(14) + t140) * pkin(3)) * t102 + (t102 * t106 * MDP(7) - qJDD(1) * MDP(10) + (-t108 - t84) * MDP(13) + (MDP(15) * pkin(3) + MDP(12)) * g(3)) * t100; ((t100 * t77 + t102 * t74) * qJD(1) + t110 - t112) * MDP(15) + t106 * t115;];
tau = t1;
