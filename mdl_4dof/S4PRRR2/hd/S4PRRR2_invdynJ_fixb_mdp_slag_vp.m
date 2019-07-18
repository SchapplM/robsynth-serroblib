% Calculate vector of inverse dynamics joint torques for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:26
% EndTime: 2019-07-18 13:27:26
% DurationCPUTime: 0.36s
% Computational Cost: add. (204->77), mult. (293->101), div. (0->0), fcn. (148->10), ass. (0->44)
t80 = cos(qJ(4));
t102 = MDP(10) * t80;
t77 = sin(qJ(4));
t109 = MDP(9) * t77 + t102;
t81 = cos(qJ(3));
t108 = pkin(1) * t81;
t91 = -qJD(3) - qJD(4);
t70 = qJD(2) - t91;
t78 = sin(qJ(3));
t107 = t70 * t78;
t74 = qJDD(2) + qJDD(3);
t69 = qJDD(4) + t74;
t106 = t77 * t69;
t105 = t80 * t69;
t103 = pkin(1) * qJD(2);
t101 = qJD(3) * t81;
t100 = qJD(4) * t77;
t99 = qJD(4) * t80;
t76 = qJ(2) + qJ(3);
t75 = qJD(2) + qJD(3);
t97 = -qJD(2) - t75;
t95 = -qJD(3) + t75;
t94 = -qJD(4) + t70;
t93 = qJDD(2) * t78;
t92 = -t69 - qJDD(2);
t67 = qJDD(2) * t108;
t89 = t78 * t103;
t54 = pkin(2) * t74 - qJD(3) * t89 + t67;
t73 = qJ(4) + t76;
t64 = sin(t73);
t65 = cos(t73);
t90 = g(1) * t65 + g(3) * t64 + t80 * t54;
t87 = MDP(9) * t94;
t86 = t91 * t70;
t84 = -g(1) * t64 + g(3) * t65 + t89 * t100 - t77 * t54;
t57 = t69 * MDP(8);
t71 = sin(t76);
t72 = cos(t76);
t83 = t57 + t74 * MDP(5) + (g(1) * t72 + g(3) * t71 + t67) * MDP(6) + (-g(1) * t71 + g(3) * t72) * MDP(7);
t82 = cos(qJ(2));
t79 = sin(qJ(2));
t66 = pkin(2) + t108;
t55 = pkin(2) * t75 + t81 * t103;
t1 = [(qJDD(1) + g(2)) * MDP(1); qJDD(2) * MDP(2) + (g(1) * t82 + g(3) * t79) * MDP(3) + (-g(1) * t79 + g(3) * t82) * MDP(4) + (t66 * t105 + t90) * MDP(9) + (-t66 * t106 + t84) * MDP(10) + ((t74 * MDP(6) + (t97 * MDP(7) + t109 * (-qJD(2) - t70)) * qJD(3)) * t81 + ((-qJDD(2) - t74) * MDP(7) + t97 * MDP(6) * qJD(3) + (-MDP(10) * t86 + t92 * MDP(9)) * t77 + ((-qJD(2) * qJD(4) + t86) * MDP(9) + t92 * MDP(10)) * t80) * t78) * pkin(1) + t83 + t109 * qJD(4) * (-t66 * t70 - t55); (-t55 * t100 + t90) * MDP(9) + (-t55 * t99 + t84) * MDP(10) + ((-t70 * t100 + t105) * MDP(9) + (-t70 * t99 - t106) * MDP(10)) * pkin(2) + ((-MDP(7) - t109) * t93 + ((-t77 * t70 * MDP(10) + t95 * MDP(6) + t80 * t87) * t78 + (t95 * MDP(7) + t109 * (-qJD(3) + t70)) * t81) * qJD(2)) * pkin(1) + t83; t57 + t90 * MDP(9) + t84 * MDP(10) + (t94 * t102 + t77 * t87) * t55 + (-t109 * t93 + ((-t77 * t101 + t80 * t107 - t78 * t99) * MDP(9) + (-t80 * t101 - t77 * t107) * MDP(10)) * qJD(2)) * pkin(1);];
tau  = t1;
