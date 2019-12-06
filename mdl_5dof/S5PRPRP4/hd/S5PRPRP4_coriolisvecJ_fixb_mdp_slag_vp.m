% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP4
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:08
% EndTime: 2019-12-05 15:36:11
% DurationCPUTime: 0.57s
% Computational Cost: add. (507->120), mult. (1194->174), div. (0->0), fcn. (774->6), ass. (0->65)
t128 = cos(qJ(4));
t126 = sin(qJ(4));
t129 = cos(qJ(2));
t150 = qJD(1) * t129;
t115 = qJD(2) * pkin(2) + t150;
t125 = cos(pkin(8));
t127 = sin(qJ(2));
t151 = qJD(1) * t127;
t117 = t125 * t151;
t124 = sin(pkin(8));
t102 = t124 * t115 + t117;
t96 = qJD(2) * pkin(6) + t102;
t159 = t126 * t96;
t92 = qJD(3) * t128 - t159;
t164 = qJD(5) - t92;
t111 = t124 * t127 - t125 * t129;
t108 = t111 * qJD(2);
t104 = qJD(1) * t108;
t163 = qJD(3) * qJD(4) - t104;
t89 = -qJD(4) * pkin(4) + t164;
t93 = qJD(3) * t126 + t128 * t96;
t90 = qJD(4) * qJ(5) + t93;
t122 = t126 ^ 2;
t123 = t128 ^ 2;
t162 = (t122 - t123) * MDP(7);
t161 = MDP(14) * (t122 + t123);
t143 = MDP(11) + MDP(13);
t142 = MDP(12) - MDP(15);
t160 = pkin(2) * t125;
t157 = t163 * t128;
t112 = t124 * t129 + t125 * t127;
t130 = qJD(4) ^ 2;
t156 = t112 * t130;
t118 = pkin(2) * t124 + pkin(6);
t155 = t118 * t130;
t154 = t126 * t128;
t107 = t124 * t150 + t117;
t149 = qJD(2) * t107;
t148 = qJD(2) * t128;
t146 = qJD(4) * t126;
t145 = qJD(4) * t128;
t87 = t163 * t126 + t96 * t145;
t141 = t92 + t159;
t133 = -pkin(4) * t128 - qJ(5) * t126 - pkin(3);
t110 = t133 - t160;
t116 = t124 * t151;
t101 = t115 * t125 - t116;
t91 = t133 * qJD(2) - t101;
t139 = qJD(2) * t110 + t91;
t95 = -qJD(2) * pkin(3) - t101;
t138 = qJD(2) * (-pkin(3) - t160) + t95;
t106 = t112 * qJD(2);
t103 = qJD(1) * t106;
t137 = -t103 - t155;
t136 = MDP(16) * t118 + MDP(14);
t135 = pkin(4) * t126 - qJ(5) * t128;
t134 = t126 * t89 + t128 * t90;
t105 = t135 * qJD(4) - t126 * qJD(5);
t88 = (t112 * qJD(1) + t105) * qJD(2);
t132 = -qJD(2) * t105 - t155 - t88;
t131 = qJD(2) ^ 2;
t113 = t135 * qJD(2);
t109 = t125 * t150 - t116;
t86 = (qJD(5) - t159) * qJD(4) + t157;
t1 = [(-t101 * t106 + t103 * t111) * MDP(5) + (t106 * t91 + t111 * t88) * MDP(16) + t143 * (t108 * t146 - t128 * t156 + (-t106 * t128 + t111 * t146) * qJD(2)) + t142 * (t108 * t145 + t126 * t156 + (t106 * t126 + t111 * t145) * qJD(2)) + (-MDP(3) * t127 - MDP(4) * t129) * t131 + (-t104 * MDP(5) + (t126 * t87 + t128 * t86 + t89 * t145 - t90 * t146) * MDP(16)) * t112 - (t134 * MDP(16) + t102 * MDP(5) + qJD(2) * t161) * t108; (t101 * t107 - t102 * t109 + (-t103 * t125 - t104 * t124) * pkin(2)) * MDP(5) + (t110 * t88 + (t105 - t107) * t91) * MDP(16) + (-0.2e1 * qJD(4) * t162 - t109 * t161) * qJD(2) + (t130 * MDP(8) + t137 * MDP(11) + t132 * MDP(13) + t86 * MDP(14) + (-t109 * t90 + t118 * t86) * MDP(16) + ((t109 + t138) * MDP(12) + (-t109 - t139) * MDP(15) + t136 * t89) * qJD(4)) * t128 + (-t130 * MDP(9) + (-t137 - t149) * MDP(12) + t87 * MDP(14) + (t132 + t149) * MDP(15) + (-t109 * t89 + t118 * t87) * MDP(16) + (t138 * MDP(11) + t139 * MDP(13) + 0.2e1 * MDP(6) * t148 - t136 * t90) * qJD(4)) * t126 + t143 * (t107 * t148 + t109 * t146); (t134 * qJD(4) + t86 * t126 - t87 * t128) * MDP(16) + (-t143 * t126 - t142 * t128) * t130; (qJ(5) * t86 - t113 * t91 + t164 * t90 - t89 * t93) * MDP(16) + (-MDP(6) * t154 + t162) * t131 + (t141 * MDP(12) + (0.2e1 * qJD(5) - t141) * MDP(15) + t143 * t93) * qJD(4) + ((-t95 * MDP(11) - t91 * MDP(13) + t113 * MDP(15)) * t126 + (-t95 * MDP(12) + t113 * MDP(13) + t91 * MDP(15)) * t128) * qJD(2) + (-pkin(4) * MDP(16) - t143) * t87 - t142 * t157; -t131 * MDP(13) * t154 + (-t122 * t131 - t130) * MDP(15) + (qJD(2) * t126 * t91 - qJD(4) * t90 + t87) * MDP(16);];
tauc = t1;
