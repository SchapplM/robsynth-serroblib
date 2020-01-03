% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:46
% EndTime: 2019-12-31 17:35:47
% DurationCPUTime: 0.27s
% Computational Cost: add. (229->69), mult. (501->116), div. (0->0), fcn. (310->6), ass. (0->49)
t95 = qJD(3) + qJD(4);
t99 = sin(qJ(4));
t128 = qJD(4) * t99;
t100 = sin(qJ(3));
t102 = cos(qJ(4));
t103 = cos(qJ(3));
t87 = t102 * t100 + t99 * t103;
t138 = -pkin(3) * t128 + t87 * qJD(2);
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t125 = t98 * MDP(14);
t137 = t101 * MDP(15) + t125;
t136 = (-t101 ^ 2 + t98 ^ 2) * MDP(10);
t135 = t95 * pkin(4);
t123 = qJD(2) * t100;
t119 = t103 * qJD(2);
t90 = qJD(3) * pkin(3) + t119;
t134 = (t102 * t123 + t99 * t90) * t95;
t120 = t101 * qJD(5);
t107 = t87 * qJD(3);
t122 = qJD(4) * t102;
t113 = t100 * t122;
t115 = t90 * t128;
t73 = t115 + (t107 + t113) * qJD(2);
t114 = t99 * t123;
t81 = t102 * t90 - t114;
t78 = -t81 - t135;
t133 = t78 * t120 + t73 * t98;
t132 = t95 * t114;
t130 = MDP(9) * t98;
t104 = qJD(5) ^ 2;
t129 = t104 * t98;
t127 = t101 * t104;
t126 = t102 * t103;
t124 = t98 * qJD(5);
t118 = t104 * MDP(12);
t112 = -t78 * t95 - (qJD(3) * t119 + qJD(4) * t90) * t102 + t132;
t110 = pkin(7) * t104 - t134;
t109 = qJD(5) * (t81 - t135);
t86 = t99 * t100 - t126;
t108 = -t90 * t122 + t132;
t106 = MDP(11) * t127 + (-MDP(7) * t113 + (-t87 * MDP(7) - MDP(8) * t126) * qJD(3)) * qJD(2) + 0.2e1 * (t130 * t101 - t136) * qJD(5) * t95;
t94 = t95 ^ 2;
t93 = -t102 * pkin(3) - pkin(4);
t92 = t99 * pkin(3) + pkin(7);
t76 = t78 * t124;
t75 = t87 * qJD(4) + t107;
t74 = t86 * t95;
t1 = [t137 * t104; (t74 * t124 - t87 * t127) * MDP(14) + (t74 * t120 + t87 * t129) * MDP(15) + (-t100 * MDP(4) - t103 * MDP(5)) * qJD(3) ^ 2 + (-t75 * MDP(7) + t74 * MDP(8) + (-t101 * t75 + t86 * t124) * MDP(14) + (t86 * t120 + t75 * t98) * MDP(15)) * t95; -MDP(7) * t115 + t108 * MDP(8) - t98 * t118 + (-t73 * t101 - t92 * t127 + t76) * MDP(14) + (t92 * t129 + t133) * MDP(15) + (t138 * MDP(7) + (t138 * t101 + t93 * t124) * MDP(14) + (t93 * t120 - t138 * t98) * MDP(15)) * t95 + t106 + (t95 * MDP(8) + t137 * qJD(5)) * (-pkin(3) * t122 - t86 * qJD(2)); (-t115 + t134) * MDP(7) + (t81 * t95 + t108) * MDP(8) + t76 * MDP(14) + t133 * MDP(15) + (MDP(14) * t109 + t110 * MDP(15) - t118) * t98 + ((-t110 - t73) * MDP(14) + MDP(15) * t109) * t101 + t106; t112 * t125 + t94 * t136 + (t112 * MDP(15) - t94 * t130) * t101;];
tauc = t1;
