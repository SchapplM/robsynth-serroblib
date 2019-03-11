% Calculate vector of inverse dynamics joint torques for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:09
% DurationCPUTime: 0.33s
% Computational Cost: add. (203->76), mult. (292->99), div. (0->0), fcn. (148->10), ass. (0->44)
t80 = cos(qJ(3));
t102 = t80 * MDP(9);
t77 = sin(qJ(3));
t109 = t77 * MDP(8) + t102;
t81 = cos(qJ(2));
t108 = t81 * pkin(1);
t74 = qJDD(1) + qJDD(2);
t69 = qJDD(3) + t74;
t107 = t77 * t69;
t91 = -qJD(2) - qJD(3);
t70 = qJD(1) - t91;
t106 = t77 * t70;
t105 = t80 * t69;
t104 = pkin(1) * qJD(1);
t101 = qJD(2) * t81;
t100 = qJD(3) * t77;
t99 = qJD(3) * t80;
t76 = qJ(1) + qJ(2);
t75 = qJD(1) + qJD(2);
t97 = -qJD(1) - t75;
t95 = -qJD(2) + t75;
t94 = -qJD(3) + t70;
t78 = sin(qJ(2));
t93 = qJDD(1) * t78;
t92 = -t69 - qJDD(1);
t90 = t78 * t104;
t88 = t94 * MDP(8);
t87 = t91 * t70;
t67 = qJDD(1) * t108;
t54 = t74 * pkin(2) - qJD(2) * t90 + t67;
t73 = qJ(3) + t76;
t64 = sin(t73);
t65 = cos(t73);
t86 = g(1) * t64 - g(2) * t65 + t80 * t54;
t84 = g(1) * t65 + g(2) * t64 + t90 * t100 - t77 * t54;
t57 = t69 * MDP(7);
t71 = sin(t76);
t72 = cos(t76);
t83 = t57 + t74 * MDP(4) + (g(1) * t72 + g(2) * t71) * MDP(6) + (g(1) * t71 - g(2) * t72 + t67) * MDP(5);
t82 = cos(qJ(1));
t79 = sin(qJ(1));
t66 = pkin(2) + t108;
t55 = t75 * pkin(2) + t104 * t81;
t1 = [qJDD(1) * MDP(1) + (g(1) * t79 - g(2) * t82) * MDP(2) + (g(1) * t82 + g(2) * t79) * MDP(3) + (t105 * t66 + t86) * MDP(8) + (-t107 * t66 + t84) * MDP(9) + ((t74 * MDP(5) + (MDP(6) * t97 + t109 * (-qJD(1) - t70)) * qJD(2)) * t81 + ((-qJDD(1) - t74) * MDP(6) + t97 * MDP(5) * qJD(2) + (MDP(8) * t92 - MDP(9) * t87) * t77 + ((-qJD(1) * qJD(3) + t87) * MDP(8) + t92 * MDP(9)) * t80) * t78) * pkin(1) + t83 + t109 * qJD(3) * (-t66 * t70 - t55); (-t55 * t100 + t86) * MDP(8) + (-t55 * t99 + t84) * MDP(9) + ((-t100 * t70 + t105) * MDP(8) + (-t70 * t99 - t107) * MDP(9)) * pkin(2) + ((-MDP(6) - t109) * t93 + ((MDP(5) * t95 - MDP(9) * t106 + t80 * t88) * t78 + (MDP(6) * t95 + t109 * (-qJD(2) + t70)) * t81) * qJD(1)) * pkin(1) + t83; t57 + t86 * MDP(8) + t84 * MDP(9) + (t102 * t94 + t77 * t88) * t55 + (-t109 * t93 + ((-t101 * t77 + (t80 * t70 - t99) * t78) * MDP(8) + (-t101 * t80 - t106 * t78) * MDP(9)) * qJD(1)) * pkin(1);];
tau  = t1;
