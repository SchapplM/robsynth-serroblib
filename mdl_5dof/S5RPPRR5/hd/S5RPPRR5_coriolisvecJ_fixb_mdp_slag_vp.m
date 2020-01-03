% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:39
% EndTime: 2019-12-31 17:56:40
% DurationCPUTime: 0.28s
% Computational Cost: add. (321->63), mult. (551->96), div. (0->0), fcn. (245->6), ass. (0->44)
t89 = cos(qJ(4));
t111 = qJD(5) * t89;
t88 = cos(qJ(5));
t108 = t88 * MDP(17);
t86 = sin(qJ(5));
t94 = t86 * MDP(16) + t108;
t124 = t94 * t111;
t102 = qJD(1) - qJD(4);
t103 = qJD(3) * qJD(1);
t112 = qJD(4) * t89;
t87 = sin(qJ(4));
t113 = qJD(4) * t87;
t76 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t74 = t76 * qJD(1) + qJD(3);
t77 = sin(pkin(8)) * pkin(1) + qJ(3);
t75 = qJD(1) * t77;
t61 = t87 * t103 + t75 * t112 + t74 * t113;
t123 = -t61 - (t87 * t74 + t89 * t75) * t102;
t95 = t87 * t76 + t89 * t77;
t122 = t61 + (t87 * qJD(3) + t95 * qJD(4)) * t102;
t121 = MDP(6) * qJD(1) + t75 * MDP(7);
t120 = (t86 ^ 2 - t88 ^ 2) * MDP(12);
t119 = t102 * pkin(4);
t90 = qJD(5) ^ 2;
t116 = t87 * t90;
t110 = t102 * qJD(5) * t120;
t109 = t88 * MDP(11);
t107 = t90 * MDP(13);
t106 = t90 * MDP(14);
t104 = MDP(17) * qJD(5);
t101 = t102 * t109;
t100 = qJD(4) * t111;
t64 = t89 * t74 - t87 * t75;
t62 = -t64 + t119;
t93 = -t89 * t103 - t74 * t112 + t75 * t113;
t99 = t102 * t62 + t93;
t98 = t62 + t64 + t119;
t96 = t89 * t76 - t87 * t77;
t66 = t89 * qJD(3) + t96 * qJD(4);
t97 = -(pkin(4) - t96) * t102 - t62 - t66;
t92 = pkin(7) * t90 - t123;
t91 = (-pkin(7) + t95) * t90 - t122;
t80 = t102 ^ 2;
t1 = [t122 * MDP(9) + (t102 * t66 - t93) * MDP(10) - 0.2e1 * t110 + 0.2e1 * t121 * qJD(3) + (-t91 * MDP(16) + t97 * t104 - t107) * t88 + (t106 + t91 * MDP(17) + (t97 * MDP(16) + 0.2e1 * t101) * qJD(5)) * t86; t94 * t90; (-t86 * t100 - t88 * t116) * MDP(16) + (-t88 * t100 + t86 * t116) * MDP(17) + (-t121 + t124) * qJD(1) - (-t124 + (t89 * MDP(10) + (t88 * MDP(16) - t86 * MDP(17) + MDP(9)) * t87) * t102) * t102; t123 * MDP(9) + (-t102 * t64 + t93) * MDP(10) + 0.2e1 * t110 + (-t92 * MDP(16) + t98 * t104 + t107) * t88 + (-t106 + t92 * MDP(17) + (t98 * MDP(16) - 0.2e1 * t101) * qJD(5)) * t86; t99 * t108 + t80 * t120 + (t99 * MDP(16) - t80 * t109) * t86;];
tauc = t1;
