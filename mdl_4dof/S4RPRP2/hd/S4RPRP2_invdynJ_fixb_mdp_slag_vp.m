% Calculate vector of inverse dynamics joint torques for
% S4RPRP2
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:54
% EndTime: 2019-03-08 18:30:55
% DurationCPUTime: 0.31s
% Computational Cost: add. (246->91), mult. (358->106), div. (0->0), fcn. (169->4), ass. (0->46)
t111 = (qJD(1) * qJD(2));
t98 = -pkin(1) - pkin(2);
t79 = t98 * qJD(1) + qJD(2);
t103 = qJD(3) * t79 + t111;
t112 = (qJ(2) * qJDD(1));
t101 = -t103 - t112;
t96 = cos(qJ(3));
t86 = t96 * t98;
t94 = sin(qJ(3));
t128 = -t94 * qJ(2) + t86;
t95 = sin(qJ(1));
t121 = t95 * t94;
t97 = cos(qJ(1));
t72 = -t97 * t96 - t121;
t120 = t97 * t94;
t73 = -t95 * t96 + t120;
t127 = -g(1) * t72 - g(2) * t73;
t126 = -g(1) * t73 + g(2) * t72;
t110 = qJD(1) - qJD(3);
t119 = t97 * pkin(1) + t95 * qJ(2);
t118 = g(1) * t95 - g(2) * t97;
t92 = qJDD(1) - qJDD(3);
t117 = t92 * MDP(7);
t115 = pkin(1) * qJDD(1);
t113 = qJ(2) * qJD(1);
t109 = 2 * t111;
t108 = qJD(3) * t113;
t78 = t98 * qJDD(1) + qJDD(2);
t107 = -t101 * t96 + t94 * t78;
t106 = t110 ^ 2;
t105 = qJDD(2) - t115;
t104 = g(1) * t97 + g(2) * t95;
t69 = -t94 * t113 + t96 * t79;
t77 = t96 * qJ(2) + t94 * t98;
t102 = (-t108 + t78) * t96;
t65 = -t94 * t108 + t107;
t100 = -t102 + t126;
t99 = qJD(1) ^ 2;
t88 = t97 * qJ(2);
t85 = t96 * pkin(3) + pkin(2);
t70 = t96 * t113 + t94 * t79;
t68 = -t94 * qJD(2) - t77 * qJD(3);
t67 = t96 * qJD(2) + qJD(3) * t128;
t66 = -pkin(3) * t110 + t69;
t64 = -t92 * pkin(3) + t101 * t94 + t102;
t1 = [qJDD(1) * MDP(1) + t118 * MDP(2) + t104 * MDP(3) + (-qJDD(2) + 0.2e1 * t115 + t118) * MDP(4) + (-t104 + t109 + (2 * t112)) * MDP(5) + (-t105 * pkin(1) - g(1) * (-t95 * pkin(1) + t88) - g(2) * t119 + (t109 + t112) * qJ(2)) * MDP(6) + t117 + (-t68 * t110 - t86 * t92 + ((qJDD(1) + t92) * qJ(2) + t103) * t94 + t100) * MDP(8) + (t110 * t67 + t77 * t92 - t127 + t65) * MDP(9) + (t65 * t77 + t70 * t67 + t64 * (-pkin(3) + t128) + t66 * t68 - g(1) * (pkin(3) * t120 + t88 + (-pkin(1) - t85) * t95) - g(2) * (pkin(3) * t121 + t97 * t85 + t119)) * MDP(10); -qJDD(1) * MDP(4) - t99 * MDP(5) + (-t99 * qJ(2) + t105 - t118) * MDP(6) - t118 * MDP(10) + (-t92 * MDP(8) + (-t110 * t70 + t64) * MDP(10) - MDP(9) * t106) * t96 + (t92 * MDP(9) + (t110 * t66 + t65) * MDP(10) - MDP(8) * t106) * t94; -t117 - t100 * MDP(8) + (-t110 * t69 - t107 + t127) * MDP(9) + (t101 * MDP(8) + MDP(9) * t108) * t94 + (-t110 * MDP(8) + (-t69 + t66) * MDP(10)) * t70 + (t64 - t126) * MDP(10) * pkin(3); (qJDD(4) + g(3)) * MDP(10);];
tau  = t1;
