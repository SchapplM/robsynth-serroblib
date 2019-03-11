% Calculate vector of inverse dynamics joint torques for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:21
% EndTime: 2019-03-08 18:25:22
% DurationCPUTime: 0.37s
% Computational Cost: add. (228->78), mult. (293->100), div. (0->0), fcn. (148->10), ass. (0->45)
t83 = cos(qJ(4));
t101 = t83 * MDP(10);
t81 = sin(qJ(4));
t111 = t81 * MDP(9) + t101;
t84 = cos(qJ(3));
t110 = t84 * pkin(2);
t78 = qJDD(2) + qJDD(3);
t73 = qJDD(4) + t78;
t109 = t81 * t73;
t93 = -qJD(3) - qJD(4);
t76 = qJD(2) - t93;
t108 = t81 * t76;
t107 = t83 * t73;
t106 = pkin(2) * qJD(2);
t104 = qJD(3) * t84;
t103 = qJD(4) * t81;
t102 = qJD(4) * t83;
t80 = qJD(2) + qJD(3);
t99 = -qJD(2) - t80;
t97 = -qJD(3) + t80;
t96 = -qJD(4) + t76;
t82 = sin(qJ(3));
t95 = qJDD(2) * t82;
t94 = -t73 - qJDD(2);
t79 = pkin(7) + qJ(2);
t92 = t82 * t106;
t77 = qJ(3) + t79;
t90 = t96 * MDP(9);
t89 = t93 * t76;
t71 = qJDD(2) * t110;
t55 = t78 * pkin(3) - qJD(3) * t92 + t71;
t69 = qJ(4) + t77;
t64 = sin(t69);
t65 = cos(t69);
t88 = g(1) * t64 - g(2) * t65 + t83 * t55;
t87 = g(1) * t65 + g(2) * t64 + t92 * t103 - t81 * t55;
t66 = t73 * MDP(8);
t67 = sin(t77);
t68 = cos(t77);
t85 = t66 + t78 * MDP(5) + (g(1) * t68 + g(2) * t67) * MDP(7) + (g(1) * t67 - g(2) * t68 + t71) * MDP(6);
t75 = cos(t79);
t74 = sin(t79);
t70 = pkin(3) + t110;
t56 = t80 * pkin(3) + t84 * t106;
t1 = [(qJDD(1) - g(3)) * MDP(1); qJDD(2) * MDP(2) + (g(1) * t74 - g(2) * t75) * MDP(3) + (g(1) * t75 + g(2) * t74) * MDP(4) + (t70 * t107 + t88) * MDP(9) + (-t70 * t109 + t87) * MDP(10) + ((t78 * MDP(6) + (t99 * MDP(7) + t111 * (-qJD(2) - t76)) * qJD(3)) * t84 + ((-qJDD(2) - t78) * MDP(7) + t99 * MDP(6) * qJD(3) + (-MDP(10) * t89 + MDP(9) * t94) * t81 + ((-qJD(2) * qJD(4) + t89) * MDP(9) + t94 * MDP(10)) * t83) * t82) * pkin(2) + t85 + t111 * qJD(4) * (-t70 * t76 - t56); (-t56 * t103 + t88) * MDP(9) + (-t56 * t102 + t87) * MDP(10) + ((-t76 * t103 + t107) * MDP(9) + (-t76 * t102 - t109) * MDP(10)) * pkin(3) + ((-MDP(7) - t111) * t95 + ((-MDP(10) * t108 + MDP(6) * t97 + t83 * t90) * t82 + (MDP(7) * t97 + t111 * (-qJD(3) + t76)) * t84) * qJD(2)) * pkin(2) + t85; t66 + t88 * MDP(9) + t87 * MDP(10) + (t96 * t101 + t81 * t90) * t56 + (-t111 * t95 + ((-t81 * t104 + (t83 * t76 - t102) * t82) * MDP(9) + (-t83 * t104 - t82 * t108) * MDP(10)) * qJD(2)) * pkin(2);];
tau  = t1;
