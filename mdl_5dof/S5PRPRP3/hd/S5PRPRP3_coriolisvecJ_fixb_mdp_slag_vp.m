% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:43
% EndTime: 2019-12-05 15:33:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (410->96), mult. (980->151), div. (0->0), fcn. (635->6), ass. (0->61)
t102 = sin(qJ(4));
t104 = cos(qJ(4));
t100 = sin(pkin(8));
t105 = cos(qJ(2));
t123 = t105 * qJD(1);
t91 = qJD(2) * pkin(2) + t123;
t101 = cos(pkin(8));
t103 = sin(qJ(2));
t127 = qJD(1) * t103;
t93 = t101 * t127;
t79 = t100 * t91 + t93;
t115 = t79 + (pkin(6) + qJ(5)) * qJD(2);
t70 = t104 * qJD(3) - t115 * t102;
t71 = t102 * qJD(3) + t115 * t104;
t139 = qJD(4) * t71;
t98 = t102 ^ 2;
t99 = t104 ^ 2;
t138 = (t98 - t99) * MDP(7);
t137 = t102 * MDP(11) + t104 * MDP(12);
t131 = t98 + t99;
t136 = t131 * MDP(13);
t106 = qJD(4) ^ 2;
t88 = t100 * t103 - t101 * t105;
t84 = t88 * qJD(2);
t81 = qJD(1) * t84;
t114 = qJD(2) * qJD(5) - t81;
t68 = -t114 * t102 - t139;
t135 = (t68 + t139) * MDP(14) - t106 * MDP(12);
t134 = t101 * pkin(2);
t130 = qJD(4) * pkin(4);
t69 = t70 + t130;
t133 = t69 - t70;
t95 = t100 * pkin(2) + pkin(6);
t129 = qJ(5) + t95;
t128 = qJD(4) * t88;
t121 = qJD(2) * qJD(4);
t120 = -t104 * pkin(4) - pkin(3);
t119 = t102 * t121;
t92 = t100 * t127;
t78 = t101 * t91 - t92;
t116 = qJD(4) * t129;
t113 = -t102 * t69 + t104 * t71;
t111 = t100 * t105 + t101 * t103;
t82 = t111 * qJD(2);
t80 = qJD(1) * t82;
t83 = t100 * t123 + t93;
t110 = qJD(2) * t83 - t106 * t95 - t80;
t74 = -qJD(2) * pkin(3) - t78;
t85 = t101 * t123 - t92;
t109 = qJD(4) * (qJD(2) * (-pkin(3) - t134) + t74 + t85);
t67 = t70 * qJD(4) + t114 * t104;
t108 = -t106 * MDP(11) + (-qJD(4) * t69 + t67) * MDP(14);
t107 = qJD(2) ^ 2;
t94 = pkin(4) * t119;
t87 = t129 * t104;
t86 = t129 * t102;
t77 = -t102 * qJD(5) - t104 * t116;
t76 = t104 * qJD(5) - t102 * t116;
t73 = t94 + t80;
t72 = t120 * qJD(2) + qJD(5) - t78;
t1 = [(-t78 * t82 + t80 * t88) * MDP(5) + (t72 * t82 + t73 * t88) * MDP(14) + (-t103 * MDP(3) - t105 * MDP(4)) * t107 - (t113 * MDP(14) + t79 * MDP(5) - t137 * qJD(4)) * t84 + (-t81 * MDP(5) - t135 * t102 + t108 * t104) * t111 + ((t102 * t128 - t104 * t82) * MDP(11) + (t102 * t82 + t104 * t128) * MDP(12) - t84 * t136) * qJD(2); (t78 * t83 - t79 * t85 + (-t100 * t81 - t101 * t80) * pkin(2)) * MDP(5) + 0.2e1 * t104 * MDP(6) * t119 - 0.2e1 * t121 * t138 + (t102 * t109 + t110 * t104) * MDP(11) + (-t110 * t102 + t104 * t109) * MDP(12) + (-t68 * t102 + t67 * t104 + (-t102 * t71 - t104 * t69) * qJD(4) + (-t102 * t77 + t104 * t76 - t131 * t85 + (-t102 * t87 + t104 * t86) * qJD(4)) * qJD(2)) * MDP(13) + (t67 * t87 - t68 * t86 + t69 * t77 + t73 * (t120 - t134) - t72 * t83 + (-t104 * t85 + t76) * t71 + (t72 * t130 + t69 * t85) * t102) * MDP(14) + (t104 * MDP(8) - t102 * MDP(9)) * t106; t108 * t102 + t135 * t104; t107 * t138 + (t133 * t71 + (-qJD(2) * t102 * t72 + t68) * pkin(4)) * MDP(14) + t137 * (-t74 * qJD(2) + t81) + (-t102 * t107 * MDP(6) + (-t130 + t133) * qJD(2) * MDP(13)) * t104; -t107 * t136 + (t94 + (t111 * qJD(1) - t113) * qJD(2)) * MDP(14);];
tauc = t1;
