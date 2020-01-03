% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:38
% EndTime: 2019-12-31 17:38:40
% DurationCPUTime: 0.33s
% Computational Cost: add. (248->86), mult. (532->128), div. (0->0), fcn. (307->6), ass. (0->46)
t104 = sin(pkin(8));
t105 = cos(pkin(8));
t107 = sin(qJ(2));
t129 = t107 * qJD(1);
t123 = qJD(2) * t129;
t109 = cos(qJ(2));
t125 = t109 * qJD(1);
t91 = (qJD(3) + t125) * qJD(2);
t77 = t104 * t91 - t105 * t123;
t80 = t104 * t125 - t105 * t129;
t143 = (t104 * qJD(3) - t80) * qJD(2) + t77;
t111 = qJD(5) ^ 2;
t110 = -pkin(2) - pkin(3);
t133 = t105 * qJ(3) + t104 * t110;
t142 = -t111 * (-pkin(6) + t133) + t143;
t86 = -t109 * t104 + t107 * t105;
t106 = sin(qJ(5));
t108 = cos(qJ(5));
t141 = (t106 ^ 2 - t108 ^ 2) * MDP(12);
t78 = t104 * t123 + t105 * t91;
t139 = qJD(2) * pkin(2);
t95 = qJD(2) * qJ(3) + t129;
t138 = t109 * t95;
t137 = t106 * t111;
t135 = t108 * t111;
t130 = t106 * qJD(5);
t128 = t108 * MDP(11);
t127 = t108 * MDP(17);
t126 = t108 * qJD(5);
t118 = qJD(3) - t125;
t87 = t110 * qJD(2) + t118;
t75 = -t104 * t95 + t105 * t87;
t73 = qJD(2) * pkin(4) - t75;
t121 = t73 * qJD(2) - t78;
t116 = t107 * t104 + t109 * t105;
t81 = t116 * qJD(1);
t119 = t105 * qJD(3) - t81;
t117 = -t104 * qJ(3) + t105 * t110;
t114 = t106 * MDP(16) + t127;
t113 = qJD(5) * (-qJD(2) * (pkin(4) - t117) - t119 - t73);
t112 = qJD(2) ^ 2;
t92 = t118 - t139;
t83 = t116 * qJD(2);
t82 = t86 * qJD(2);
t76 = t104 * t87 + t105 * t95;
t1 = [t91 * t107 * MDP(7) + (t116 * t77 + t75 * t82 + t76 * t83 + t78 * t86) * MDP(10) + (-t83 * t130 - t86 * t135) * MDP(16) + (-t83 * t126 + t86 * t137) * MDP(17) + ((-MDP(4) + MDP(6)) * t109 + (-MDP(3) - MDP(5)) * t107) * t112 + ((t138 + (-t125 + t92) * t107) * MDP(7) - t82 * MDP(8) + t83 * MDP(9) + (-t108 * t82 - t116 * t130) * MDP(16) + (t106 * t82 - t116 * t126) * MDP(17)) * qJD(2); (t91 * qJ(3) + t95 * qJD(3) + (-t138 + (-t92 - t139) * t107) * qJD(1)) * MDP(7) + t143 * MDP(8) + (t119 * qJD(2) + t78) * MDP(9) + (t78 * t133 - t77 * t117 - t76 * t81 + t75 * t80 + (-t75 * t104 + t76 * t105) * qJD(3)) * MDP(10) - MDP(13) * t135 + MDP(14) * t137 + (t106 * t113 + t142 * t108) * MDP(16) + (-t142 * t106 + t108 * t113) * MDP(17) + 0.2e1 * (qJD(3) * MDP(6) + (t128 * t106 - t141) * qJD(5)) * qJD(2); -t112 * MDP(6) + (-t95 + t129) * MDP(7) * qJD(2) + (-t77 * MDP(10) - t112 * MDP(9) + (-t76 * MDP(10) + 0.2e1 * qJD(5) * t114) * qJD(2)) * t105 + (-t112 * MDP(8) + (qJD(2) * t75 + t78) * MDP(10) + (-t108 * MDP(16) + t106 * MDP(17)) * (t111 + t112)) * t104; -t114 * t111; t112 * t141 + t121 * t127 + (t121 * MDP(16) - t112 * t128) * t106;];
tauc = t1;
