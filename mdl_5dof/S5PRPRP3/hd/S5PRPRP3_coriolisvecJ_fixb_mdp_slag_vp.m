% Calculate Coriolis joint torque vector for
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:02
% EndTime: 2021-01-15 15:14:06
% DurationCPUTime: 0.68s
% Computational Cost: add. (516->136), mult. (1234->193), div. (0->0), fcn. (795->6), ass. (0->76)
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t130 = cos(qJ(2));
t152 = qJD(1) * t130;
t116 = qJD(2) * pkin(2) + t152;
t126 = cos(pkin(8));
t128 = sin(qJ(2));
t153 = qJD(1) * t128;
t118 = t126 * t153;
t125 = sin(pkin(8));
t103 = t125 * t116 + t118;
t96 = qJD(2) * pkin(6) + t103;
t140 = qJ(5) * qJD(2) + t96;
t134 = t140 * t129;
t91 = qJD(3) * t127 + t134;
t123 = t127 ^ 2;
t124 = t129 ^ 2;
t168 = (t123 - t124) * MDP(7);
t167 = t103 * MDP(5);
t144 = MDP(11) + MDP(13);
t166 = pkin(2) * t126;
t165 = pkin(4) * t123;
t164 = pkin(4) * t129;
t162 = qJD(4) * pkin(4);
t122 = t129 * qJD(3);
t90 = -t140 * t127 + t122;
t89 = t90 + t162;
t163 = t89 - t90;
t117 = t125 * t153;
t102 = t116 * t126 - t117;
t95 = -qJD(2) * pkin(3) - t102;
t160 = qJD(2) * t95;
t113 = t125 * t130 + t126 * t128;
t131 = qJD(4) ^ 2;
t159 = t113 * t131;
t132 = qJD(2) ^ 2;
t158 = t129 * t132;
t120 = pkin(2) * t125 + pkin(6);
t157 = qJ(5) + t120;
t155 = t123 + t124;
t109 = t126 * t152 - t117;
t154 = MDP(16) * t109;
t107 = t125 * t152 + t118;
t151 = qJD(2) * t107;
t150 = qJD(2) * t127;
t148 = qJD(4) * t127;
t147 = qJD(4) * t129;
t146 = MDP(15) * qJD(2);
t145 = qJD(2) * qJD(5);
t143 = MDP(12) + MDP(14);
t142 = -pkin(3) - t164;
t141 = qJD(2) * t148;
t106 = t113 * qJD(2);
t104 = qJD(1) * t106;
t139 = -t120 * t131 - t104;
t138 = MDP(15) * t155;
t137 = qJD(4) * t157;
t112 = t125 * t128 - t126 * t130;
t108 = t112 * qJD(2);
t105 = qJD(1) * t108;
t136 = qJD(4) * qJD(3) - t105;
t135 = -t127 * t89 + t129 * t91;
t92 = t142 * qJD(2) + qJD(5) - t102;
t133 = (-qJD(5) - t92) * qJD(2) - t136;
t121 = -pkin(3) - t166;
t119 = pkin(4) * t141;
t115 = t142 - t166;
t111 = t157 * t129;
t110 = t157 * t127;
t98 = -qJD(5) * t127 - t129 * t137;
t97 = qJD(5) * t129 - t127 * t137;
t94 = t119 + t104;
t93 = t96 * t148;
t88 = (t105 - t145) * t127 - qJD(4) * t91;
t87 = -qJ(5) * t141 - t93 + (t136 + t145) * t129;
t1 = [(-t102 * t106 + t104 * t112) * MDP(5) + (t106 * t92 + t112 * t94) * MDP(16) + t144 * (t108 * t148 - t129 * t159 + (-t106 * t129 + t112 * t148) * qJD(2)) + t143 * (t108 * t147 + t127 * t159 + (t106 * t127 + t112 * t147) * qJD(2)) + (-t128 * MDP(3) - t130 * MDP(4)) * t132 + (-t105 * MDP(5) + (-t127 * t88 + t129 * t87 - t89 * t147 - t91 * t148) * MDP(16)) * t113 - (t135 * MDP(16) + t155 * t146 + t167) * t108; (t102 * t107 + (-t104 * t126 - t105 * t125) * pkin(2)) * MDP(5) + (-t107 * t92 - t110 * t88 + t111 * t87 + t115 * t94 + t89 * t98 + t91 * t97) * MDP(16) + (t131 * MDP(8) + t139 * MDP(11) - t94 * MDP(13) + (qJD(2) * t97 + t87) * MDP(15) - t91 * t154) * t129 + (-t131 * MDP(9) + (-t139 - t151) * MDP(12) + (t94 - t151) * MDP(14) + (-qJD(2) * t98 - t88) * MDP(15) + t89 * t154) * t127 + (t98 * MDP(13) - t97 * MDP(14) + (MDP(14) * t165 - 0.2e1 * t168) * qJD(2) + ((qJD(2) * t121 + t95) * MDP(12) + (qJD(2) * t115 + t92) * MDP(14) + (qJD(2) * t110 - t89) * MDP(15)) * t129 + (t95 * MDP(11) - t91 * MDP(15) + (pkin(4) * MDP(16) + MDP(13)) * t92 + (0.2e1 * t129 * MDP(6) + t121 * MDP(11) + (t115 - t164) * MDP(13) - t111 * MDP(15)) * qJD(2)) * t127) * qJD(4) + t144 * (t109 * t148 + t129 * t151) + (-qJD(2) * t138 + t143 * t147 - t167) * t109; (qJD(4) * t135 + t87 * t127 + t88 * t129) * MDP(16) + (-t127 * t144 - t129 * t143) * t131; -t127 * MDP(6) * t158 + t132 * t168 + (t105 - t160) * t127 * MDP(11) + (t93 + (-t127 * t96 + t122) * qJD(4) + (-t136 - t160) * t129) * MDP(12) + ((t91 - t134) * qJD(4) + (pkin(4) * t158 + t133) * t127) * MDP(13) + (-t132 * t165 + t93 + (qJ(5) * t150 + t90) * qJD(4) + t133 * t129) * MDP(14) + (-t162 + t163) * t129 * t146 + (t163 * t91 + (-t150 * t92 + t88) * pkin(4)) * MDP(16); t119 * MDP(16) - t132 * t138 + ((qJD(1) * t113 - t135) * MDP(16) + 0.2e1 * (t127 * MDP(13) + t129 * MDP(14)) * qJD(4)) * qJD(2);];
tauc = t1;
