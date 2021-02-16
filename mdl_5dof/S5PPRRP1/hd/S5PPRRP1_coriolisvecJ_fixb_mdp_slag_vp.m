% Calculate Coriolis joint torque vector for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:19
% EndTime: 2021-01-15 14:48:22
% DurationCPUTime: 0.76s
% Computational Cost: add. (432->122), mult. (1120->177), div. (0->0), fcn. (742->6), ass. (0->65)
t119 = sin(pkin(8));
t120 = cos(pkin(8));
t122 = sin(qJ(3));
t124 = cos(qJ(3));
t164 = -t119 * t122 + t120 * t124;
t108 = t119 * t124 + t120 * t122;
t106 = t108 * qJD(3);
t121 = sin(qJ(4));
t123 = cos(qJ(4));
t104 = t108 * qJD(1);
t98 = qJD(3) * pkin(6) + t104;
t132 = qJ(5) * qJD(3) + t98;
t128 = t132 * t123;
t90 = t121 * qJD(2) + t128;
t117 = t121 ^ 2;
t118 = t123 ^ 2;
t163 = (t117 - t118) * MDP(7);
t103 = t164 * qJD(1);
t147 = t117 + t118;
t162 = t147 * MDP(15);
t161 = pkin(4) * t117;
t160 = pkin(4) * t123;
t159 = -qJ(5) - pkin(6);
t156 = qJD(4) * pkin(4);
t116 = t123 * qJD(2);
t89 = -t132 * t121 + t116;
t88 = t89 + t156;
t158 = t88 - t89;
t143 = qJD(4) * t121;
t144 = qJD(3) * t123;
t157 = t103 * t143 + t104 * t144;
t97 = -qJD(3) * pkin(3) - t103;
t155 = qJD(3) * t97;
t125 = qJD(4) ^ 2;
t154 = t108 * t125;
t126 = qJD(3) ^ 2;
t149 = t123 * t126;
t100 = qJD(1) * t106;
t146 = MDP(16) * t103;
t145 = qJD(3) * t121;
t142 = qJD(4) * t123;
t139 = qJD(4) * qJD(2);
t138 = MDP(11) + MDP(13);
t137 = MDP(12) + MDP(14);
t135 = qJD(3) * t143;
t92 = pkin(4) * t135 + t100;
t115 = -pkin(3) - t160;
t134 = pkin(6) * t125 + t100;
t133 = qJD(4) * t159;
t105 = t164 * qJD(3);
t99 = qJD(1) * t105;
t131 = t99 + t139;
t130 = -qJD(3) * qJD(5) - t99;
t129 = t121 * t88 - t123 * t90;
t91 = t115 * qJD(3) + qJD(5) - t103;
t127 = (-qJD(5) - t91) * qJD(3) - t131;
t110 = t159 * t123;
t109 = t159 * t121;
t102 = -t121 * qJD(5) + t123 * t133;
t101 = t123 * qJD(5) + t121 * t133;
t95 = t103 * t142;
t93 = t98 * t143;
t87 = -t90 * qJD(4) + t130 * t121;
t86 = -qJ(5) * t135 - t93 + (-t130 + t139) * t123;
t1 = [t138 * (-t105 * t143 - t123 * t154 + (-t106 * t123 - t143 * t164) * qJD(3)) + t137 * (-t105 * t142 + t121 * t154 + (t106 * t121 - t142 * t164) * qJD(3)) + (-t106 * MDP(4) + (-MDP(5) + t162) * t105) * qJD(3) + (-t129 * t105 + t91 * t106 - t92 * t164 + (-t87 * t121 + t86 * t123 + (-t121 * t90 - t123 * t88) * qJD(4)) * t108) * MDP(16); (-t129 * qJD(4) + t86 * t121 + t87 * t123) * MDP(16) + (-t138 * t121 - t137 * t123) * t125; -t100 * MDP(4) + t157 * MDP(11) + t95 * MDP(12) + (t102 * qJD(4) + t157) * MDP(13) + (-t101 * qJD(4) + t95) * MDP(14) + (t101 * t90 + t102 * t88 - t104 * t91 + t109 * t87 - t110 * t86 + t115 * t92) * MDP(16) + (t125 * MDP(8) - t134 * MDP(11) - t92 * MDP(13) + t86 * MDP(15) - t90 * t146 + (t97 * MDP(12) + MDP(14) * t91 - MDP(15) * t88) * qJD(4)) * t123 + (-t125 * MDP(9) + t134 * MDP(12) + t92 * MDP(14) - t87 * MDP(15) + t88 * t146 + (t97 * MDP(11) - t90 * MDP(15) + (MDP(16) * pkin(4) + MDP(13)) * t91) * qJD(4)) * t121 + ((t101 * t123 - t102 * t121 - t147 * t103) * MDP(15) + (-t137 * t121 + MDP(4)) * t104 + (-0.2e1 * t163 + MDP(14) * t161 + (-MDP(12) * pkin(3) + MDP(14) * t115 - MDP(15) * t109) * t123 + (0.2e1 * t123 * MDP(6) - pkin(3) * MDP(11) + (t115 - t160) * MDP(13) + t110 * MDP(15)) * t121) * qJD(4)) * qJD(3); -t121 * MDP(6) * t149 + t126 * t163 + (-t99 - t155) * t121 * MDP(11) + (t93 + (-t121 * t98 + t116) * qJD(4) + (-t131 - t155) * t123) * MDP(12) + ((t90 - t128) * qJD(4) + (pkin(4) * t149 + t127) * t121) * MDP(13) + (-t126 * t161 + t93 + (qJ(5) * t145 + t89) * qJD(4) + t127 * t123) * MDP(14) + (-t156 + t158) * MDP(15) * t144 + (t158 * t90 + (-t91 * t145 + t87) * pkin(4)) * MDP(16); t92 * MDP(16) - t126 * t162 + (t129 * MDP(16) + 0.2e1 * (t121 * MDP(13) + MDP(14) * t123) * qJD(4)) * qJD(3);];
tauc = t1;
