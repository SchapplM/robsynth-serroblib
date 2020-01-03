% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:16
% EndTime: 2019-12-31 17:45:17
% DurationCPUTime: 0.40s
% Computational Cost: add. (214->66), mult. (486->104), div. (0->0), fcn. (292->6), ass. (0->38)
t89 = sin(pkin(7)) * pkin(1) + qJ(3);
t109 = qJD(1) * t89;
t123 = t109 * MDP(7);
t100 = cos(qJ(5));
t95 = sin(pkin(8));
t97 = cos(pkin(8));
t99 = sin(qJ(5));
t80 = t100 * t95 + t99 * t97;
t110 = qJD(1) * t80;
t122 = t110 * MDP(17);
t114 = t95 ^ 2 + t97 ^ 2;
t88 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t121 = qJD(1) * t88;
t113 = t100 * t97;
t81 = -t99 * t95 + t113;
t120 = MDP(10) * t114;
t119 = t95 * MDP(8) + t97 * MDP(9) + MDP(6);
t118 = 0.2e1 * qJD(1);
t116 = -pkin(6) + t88;
t108 = qJD(1) * t95;
t107 = t99 * t108;
t82 = qJD(4) + t109;
t105 = qJD(1) * t113;
t104 = t114 * (qJD(3) + t121);
t103 = t81 * qJD(4);
t77 = t80 * qJD(5);
t102 = t80 * qJD(4);
t101 = qJD(1) ^ 2;
t85 = qJD(5) * t107;
t83 = t95 * pkin(4) + t89;
t78 = t81 * qJD(5);
t76 = t105 - t107;
t73 = t116 * t97;
t72 = t116 * t95;
t71 = pkin(4) * t108 + t82;
t68 = qJD(5) * t105 - t85;
t67 = qJD(1) * t77;
t1 = [(-t67 * t81 - t76 * t77) * MDP(12) + (t110 * t77 + t67 * t80 - t81 * t68 - t76 * t78) * MDP(13) - t77 * qJD(5) * MDP(14) - t78 * qJD(5) * MDP(15) + (t83 * t68 + t71 * t78 + ((-t100 * t72 - t73 * t99) * qJD(5) - t103) * qJD(5)) * MDP(17) + (-t83 * t67 - t71 * t77 + ((-t100 * t73 + t72 * t99) * qJD(5) + t102) * qJD(5)) * MDP(18) + (t118 * t120 + (-t114 * t121 - t104) * MDP(11)) * qJD(4) + (0.2e1 * t123 + (t82 + t109) * MDP(11) + 0.2e1 * t122 + (qJD(1) * t81 + t76) * MDP(18) + t119 * t118) * qJD(3); (-t78 * MDP(17) + t77 * MDP(18)) * qJD(5); (-t77 * MDP(17) - t78 * MDP(18)) * qJD(5) - t119 * t101 + (-t123 + (-t114 * qJD(4) - t82) * MDP(11) - t122 - t76 * MDP(18)) * qJD(1); (t76 * qJD(5) - t85) * MDP(17) - t110 * qJD(5) * MDP(18) - t101 * t120 + ((qJD(3) + t104) * MDP(11) + (MDP(17) * t113 - MDP(18) * t80) * qJD(5)) * qJD(1); t76 * t110 * MDP(12) + (-t110 ^ 2 + t76 ^ 2) * MDP(13) + (t85 + (t76 - t105) * qJD(5)) * MDP(15) + (-qJD(1) * t103 - t71 * t76) * MDP(17) + (qJD(1) * t102 + t110 * t71) * MDP(18);];
tauc = t1;
