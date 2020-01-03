% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:02
% EndTime: 2019-12-31 17:44:03
% DurationCPUTime: 0.38s
% Computational Cost: add. (217->76), mult. (546->123), div. (0->0), fcn. (336->6), ass. (0->44)
t102 = sin(pkin(8));
t100 = t102 ^ 2;
t104 = cos(pkin(8));
t134 = t104 ^ 2 + t100;
t113 = cos(pkin(7)) * pkin(1) + t102 * qJ(4) + pkin(2);
t133 = (pkin(3) * t104 + t113) * qJD(1);
t132 = t134 * (MDP(7) + MDP(10));
t98 = sin(pkin(7)) * pkin(1) + qJ(3);
t130 = -pkin(6) + t98;
t94 = t98 * qJD(1);
t129 = t104 * (t102 * qJD(2) + t104 * t94);
t106 = sin(qJ(5));
t107 = cos(qJ(5));
t89 = t102 * t106 + t104 * t107;
t128 = qJD(1) * t89;
t127 = t102 * t107;
t126 = t104 * MDP(9);
t125 = t104 * t106;
t123 = qJD(3) * t102;
t122 = qJD(4) * t102;
t121 = t100 * MDP(11);
t120 = qJD(1) * qJD(3);
t119 = t134 * t98 * t120 + qJD(3) * t129;
t116 = qJD(1) * t127;
t115 = qJD(5) * t127;
t114 = qJD(1) * t125;
t75 = qJD(2) * t104 - t102 * t94;
t90 = -t125 + t127;
t112 = t90 * qJD(3);
t111 = t89 * qJD(3);
t79 = t89 * qJD(5);
t74 = (pkin(3) + pkin(4)) * t104 + t113;
t109 = qJD(1) ^ 2;
t95 = qJD(5) * t114;
t85 = t130 * t104;
t84 = t130 * t102;
t83 = -t114 + t116;
t80 = -qJD(5) * t125 + t115;
t78 = qJD(1) * t115 - t95;
t77 = qJD(1) * t79;
t73 = qJD(4) - t75;
t72 = qJD(3) - t133;
t68 = t74 * qJD(1) - qJD(3);
t1 = [(-t75 * t123 + t119) * MDP(8) + (t73 * t123 + (-t72 + t133) * t122 + t119) * MDP(12) + (-t77 * t90 - t79 * t83) * MDP(13) + (t128 * t79 + t77 * t89 - t78 * t90 - t80 * t83) * MDP(14) - t79 * qJD(5) * MDP(15) - t80 * qJD(5) * MDP(16) + (t68 * t80 + t74 * t78 + ((-t106 * t84 - t107 * t85) * qJD(5) + t112) * qJD(5) + 0.2e1 * t128 * t122) * MDP(18) + (-t68 * t79 - t74 * t77 + ((t106 * t85 - t107 * t84) * qJD(5) - t111) * qJD(5) + (qJD(1) * t90 + t83) * t122) * MDP(19) + 0.2e1 * t120 * t132 + 0.2e1 * (t126 * t102 + t121) * qJD(1) * qJD(4); (-MDP(18) * t80 + MDP(19) * t79) * qJD(5); (-qJD(5) * t83 + t95) * MDP(18) + t128 * qJD(5) * MDP(19) - t109 * t132 + ((t102 * t75 - t129) * MDP(8) + (-t102 * t73 - t122 - t129) * MDP(12) + (-MDP(18) * t127 + t89 * MDP(19)) * qJD(5)) * qJD(1); -t109 * t121 + (-MDP(18) * t106 - MDP(19) * t107) * qJD(5) ^ 2 + (-t109 * t126 + ((qJD(3) + t72) * MDP(12) - t128 * MDP(18) - t83 * MDP(19)) * qJD(1)) * t102; t83 * t128 * MDP(13) + (-t128 ^ 2 + t83 ^ 2) * MDP(14) + (t95 + (t83 - t116) * qJD(5)) * MDP(16) + (qJD(1) * t112 - t68 * t83) * MDP(18) + (-qJD(1) * t111 + t128 * t68) * MDP(19);];
tauc = t1;
