% Calculate vector of inverse dynamics joint torques for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:57
% EndTime: 2019-12-31 16:25:57
% DurationCPUTime: 0.39s
% Computational Cost: add. (164->95), mult. (325->116), div. (0->0), fcn. (190->6), ass. (0->47)
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t73 = sin(pkin(6));
t74 = cos(pkin(6));
t92 = g(1) * t74 + g(2) * t73;
t121 = -g(3) * t78 + t92 * t76;
t104 = qJD(2) * qJ(3);
t107 = t76 * qJD(1);
t67 = t104 + t107;
t108 = t67 * qJD(2);
t79 = -pkin(2) - pkin(5);
t98 = qJD(1) * qJD(2);
t68 = t76 * t98;
t99 = t78 * qJDD(1);
t93 = qJDD(3) - t99;
t89 = t68 + t93;
t120 = -t79 * qJDD(2) + t108 + t121 - t89;
t77 = cos(qJ(4));
t69 = qJDD(4) * t77;
t75 = sin(qJ(4));
t80 = qJD(4) ^ 2;
t118 = -t80 * t75 + t69;
t117 = (t104 + t67 - t107) * qJD(4) + qJDD(4) * t79;
t105 = qJDD(1) - g(3);
t116 = -t105 * t76 + t92 * t78;
t81 = qJD(2) ^ 2;
t112 = t81 * t77;
t72 = t77 ^ 2;
t111 = t75 ^ 2 - t72;
t110 = -t80 - t81;
t109 = qJDD(2) * pkin(2);
t106 = t78 * qJD(1);
t103 = qJDD(2) * t75;
t102 = qJDD(2) * t77;
t101 = qJDD(4) * t75;
t97 = qJD(2) * qJD(3);
t96 = qJD(2) * qJD(4);
t95 = qJDD(2) * qJ(3);
t94 = -0.2e1 * t75 * t96;
t91 = g(1) * t73 - g(2) * t74;
t88 = -t80 * t77 - t101;
t83 = t93 - t121;
t62 = t95 + t76 * qJDD(1) + (qJD(3) + t106) * qJD(2);
t82 = -g(3) * t76 + t95 + t97 - t79 * t80 + t62 + (-t92 - t98) * t78;
t66 = -qJD(2) * pkin(2) + qJD(3) - t106;
t63 = t89 - t109;
t1 = [t105 * MDP(1) - g(3) * MDP(7) + (-MDP(3) + MDP(5)) * (-t78 * qJDD(2) + t81 * t76) + (-MDP(4) + MDP(6)) * (qJDD(2) * t76 + t81 * t78) + ((qJD(2) * t66 + t62) * MDP(7) + (0.2e1 * t77 * t96 + t103) * MDP(13) + (t94 + t102) * MDP(14)) * t76 + ((-t63 + t108) * MDP(7) + (t81 * t75 - t118) * MDP(13) + (-t88 + t112) * MDP(14)) * t78; qJDD(2) * MDP(2) + (t99 + t121) * MDP(3) + t116 * MDP(4) + (t83 - 0.2e1 * t109) * MDP(5) + (-t116 + 0.2e1 * t95 + 0.2e1 * t97) * MDP(6) + (t62 * qJ(3) + t67 * qJD(3) - t63 * pkin(2) - g(3) * (t78 * pkin(2) + t76 * qJ(3)) + (-t66 * t76 - t67 * t78) * qJD(1) + t92 * (pkin(2) * t76 - qJ(3) * t78)) * MDP(7) + (t72 * qJDD(2) + t77 * t94) * MDP(8) + 0.2e1 * (-t75 * t102 + t111 * t96) * MDP(9) + t118 * MDP(10) + t88 * MDP(11) + (t117 * t77 + t82 * t75) * MDP(13) + (-t117 * t75 + t82 * t77) * MDP(14); qJDD(2) * MDP(5) - t81 * MDP(6) + (-t108 + t68 + t83 - t109) * MDP(7) + (t110 * t75 + t69) * MDP(13) + (t110 * t77 - t101) * MDP(14); t75 * MDP(8) * t112 - t111 * MDP(9) * t81 + MDP(10) * t102 - MDP(11) * t103 + qJDD(4) * MDP(12) + (-t120 * t77 + t91 * t75) * MDP(13) + (t120 * t75 + t91 * t77) * MDP(14);];
tau = t1;
