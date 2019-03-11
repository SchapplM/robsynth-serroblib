% Calculate vector of inverse dynamics joint torques for
% S4RPRR1
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:51
% EndTime: 2019-03-08 18:31:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (368->94), mult. (640->114), div. (0->0), fcn. (356->12), ass. (0->53)
t106 = sin(qJ(3));
t103 = sin(pkin(7));
t130 = pkin(1) * t103;
t109 = cos(qJ(3));
t104 = cos(pkin(7));
t95 = pkin(1) * t104 + pkin(2);
t83 = t109 * t95;
t135 = -t106 * t130 + t83;
t100 = qJ(1) + pkin(7) + qJ(3);
t93 = sin(t100);
t94 = cos(t100);
t134 = g(1) * t93 - g(2) * t94;
t102 = qJD(1) + qJD(3);
t99 = qJD(4) + t102;
t133 = qJD(4) - t99;
t123 = qJD(1) * t130;
t119 = t106 * t123;
t82 = t95 * qJD(1);
t126 = qJD(3) * t82;
t80 = t95 * qJDD(1);
t69 = (qJDD(1) * t130 + t126) * t109 - qJD(3) * t119 + t106 * t80;
t101 = qJDD(1) + qJDD(3);
t98 = qJDD(4) + t101;
t132 = pkin(3) * t98;
t92 = t98 * MDP(8);
t129 = t101 * MDP(5) + t92;
t105 = sin(qJ(4));
t73 = t106 * t82 + t109 * t123;
t128 = t105 * t73;
t108 = cos(qJ(4));
t127 = t108 * t73;
t96 = qJ(4) + t100;
t90 = sin(t96);
t91 = cos(t96);
t125 = g(1) * t91 + g(2) * t90 + qJD(4) * t128;
t72 = t109 * t82 - t119;
t70 = pkin(3) * t102 + t72;
t122 = -pkin(3) * t99 - t70;
t121 = qJD(1) * qJD(3) * t109;
t118 = -t106 * t126 + t109 * t80;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t117 = g(1) * t107 - g(2) * t110;
t116 = -t105 * t70 - t127;
t77 = t106 * t95 + t109 * t130;
t112 = (-qJDD(1) * t106 - t121) * t130 + t118;
t68 = pkin(3) * t101 + t112;
t114 = g(1) * t90 - g(2) * t91 - t105 * t69 + t108 * t68;
t113 = g(1) * t94 + g(2) * t93 - t69;
t76 = pkin(3) + t135;
t75 = t77 * qJD(3);
t74 = t135 * qJD(3);
t1 = [qJDD(1) * MDP(1) + t117 * MDP(2) + (g(1) * t110 + g(2) * t107) * MDP(3) + (t117 + (t103 ^ 2 + t104 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t83 * t101 - t75 * t102 + (-t121 + (-qJDD(1) - t101) * t106) * t130 + t118 + t134) * MDP(6) + (-t101 * t77 - t102 * t74 + t113) * MDP(7) + ((-t105 * t74 - t108 * t75) * t99 + (-t105 * t77 + t108 * t76) * t98 + ((-t105 * t76 - t108 * t77) * t99 + t116) * qJD(4) + t114) * MDP(9) + ((-(-qJD(4) * t77 - t75) * t99 - t76 * t98 - t68) * t105 + (-(qJD(4) * t76 + t74) * t99 - t77 * t98 - t69 - qJD(4) * t70) * t108 + t125) * MDP(10) + t129; (qJDD(2) - g(3)) * MDP(4); (t102 * t73 + t112 + t134) * MDP(6) + (t102 * t72 + t113) * MDP(7) + (t108 * t132 - (-t105 * t72 - t127) * t99 + t114) * MDP(9) + (-t108 * t69 + (t108 * t72 - t128) * t99 + t125 + (-t132 - t68) * t105) * MDP(10) + (t122 * MDP(9) * t105 + (t122 * MDP(10) - t73 * MDP(9)) * t108) * qJD(4) + t129; t92 + (t133 * t116 + t114) * MDP(9) + ((-t73 * t99 - t68) * t105 + (-t133 * t70 - t69) * t108 + t125) * MDP(10);];
tau  = t1;
