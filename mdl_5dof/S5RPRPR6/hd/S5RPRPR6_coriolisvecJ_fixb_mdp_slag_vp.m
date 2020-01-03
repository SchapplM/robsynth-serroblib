% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:50
% EndTime: 2019-12-31 18:17:51
% DurationCPUTime: 0.29s
% Computational Cost: add. (312->72), mult. (646->96), div. (0->0), fcn. (304->6), ass. (0->50)
t92 = qJD(1) + qJD(3);
t128 = t92 * qJ(4);
t100 = cos(qJ(3));
t134 = pkin(1) * sin(pkin(8));
t118 = qJD(1) * t134;
t89 = cos(pkin(8)) * pkin(1) + pkin(2);
t87 = t89 * qJD(1);
t98 = sin(qJ(3));
t130 = t98 * t87;
t75 = t100 * t118 + t130;
t72 = t75 + t128;
t117 = qJD(3) * t134;
t113 = qJD(1) * t117;
t73 = qJD(3) * t130 + t100 * t113;
t138 = -t72 * t92 + t73;
t99 = cos(qJ(5));
t137 = MDP(17) * t99;
t136 = MDP(6) - MDP(8);
t97 = sin(qJ(5));
t135 = (t97 ^ 2 - t99 ^ 2) * MDP(12);
t91 = t92 ^ 2;
t74 = -t100 * t87 + t98 * t118;
t132 = t74 * t92;
t123 = qJD(3) * t100;
t108 = -t98 * t117 + t89 * t123;
t76 = -qJD(4) - t108;
t131 = t76 * t92;
t90 = t92 * qJD(4);
t126 = t99 * MDP(16);
t125 = t99 * qJD(5);
t124 = qJD(4) + t74;
t102 = qJD(5) ^ 2;
t122 = t102 * MDP(13);
t121 = t102 * MDP(14);
t120 = MDP(16) * qJD(5);
t119 = MDP(17) * qJD(5);
t107 = t100 * t134 + t98 * t89;
t78 = t107 * qJD(3);
t80 = qJ(4) + t107;
t115 = t80 * t92 + t78;
t114 = -t75 + t128;
t112 = -t100 * t89 + t98 * t134 - pkin(3);
t109 = -t102 * (-pkin(7) + t112) - t131;
t106 = t98 * t113 - t87 * t123;
t70 = t106 - t90;
t105 = -(-pkin(3) - pkin(7)) * t102 + t132 + t90;
t104 = (t72 * t125 - t70 * t97) * MDP(16) - t70 * t137 + 0.2e1 * (-t97 * MDP(11) * t125 + qJD(5) * t135) * t92;
t103 = t106 - t132;
t71 = -t92 * pkin(3) + t124;
t1 = [(-t108 * t92 + t106) * MDP(7) + (-t70 - t131) * MDP(9) + (t73 * t112 - t70 * t80 + t71 * t78 - t72 * t76) * MDP(10) + (t109 * MDP(17) + t115 * t120 - t121) * t99 + (-t122 + t109 * MDP(16) + (-t115 - t72) * t119) * t97 + t104 - t136 * (t78 * t92 + t73); (t97 * MDP(17) - t126) * t102; t103 * MDP(7) + (-t103 + 0.2e1 * t90) * MDP(9) + (-t73 * pkin(3) - t70 * qJ(4) + t124 * t72 - t71 * t75) * MDP(10) + (t105 * MDP(17) + t114 * t120 - t121) * t99 + (-t122 + t105 * MDP(16) + (-t114 - t72) * t119) * t97 + t104 + t136 * (t75 * t92 - t73); -t91 * MDP(9) + t138 * MDP(10) + (MDP(16) * t97 + t137) * (-t102 - t91); t138 * t126 - t91 * t135 + (t99 * t91 * MDP(11) - MDP(17) * t138) * t97;];
tauc = t1;
