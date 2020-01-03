% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:37
% EndTime: 2019-12-31 16:50:39
% DurationCPUTime: 0.87s
% Computational Cost: add. (392->139), mult. (971->217), div. (0->0), fcn. (568->6), ass. (0->76)
t119 = sin(qJ(3));
t175 = MDP(5) * t119;
t114 = t119 ^ 2;
t121 = cos(qJ(3));
t174 = (-t121 ^ 2 + t114) * MDP(6);
t173 = -t119 * MDP(10) - t121 * MDP(11);
t118 = sin(qJ(4));
t147 = qJD(4) * t118;
t138 = t119 * t147;
t120 = cos(qJ(4));
t141 = qJD(1) * qJD(3);
t136 = t121 * t141;
t140 = qJD(3) * qJD(4);
t155 = (t136 + t140) * t120;
t88 = -qJD(1) * t138 + t155;
t172 = t118 * t88;
t110 = sin(pkin(7)) * pkin(1) + pkin(5);
t106 = t110 * qJD(1);
t95 = qJD(2) * t121 - t106 * t119;
t90 = -qJD(3) * pkin(3) - t95;
t171 = t118 * t90;
t111 = -cos(pkin(7)) * pkin(1) - pkin(2);
t100 = -pkin(3) * t121 - pkin(6) * t119 + t111;
t94 = t100 * qJD(1);
t170 = t118 * t94;
t169 = t120 * t90;
t168 = t120 * t94;
t146 = qJD(4) * t120;
t137 = t119 * t146;
t148 = qJD(3) * t121;
t89 = t118 * t140 + (t118 * t148 + t137) * qJD(1);
t167 = t121 * t89;
t96 = qJD(2) * t119 + t106 * t121;
t93 = qJD(3) * t96;
t166 = t93 * t118;
t165 = t93 * t120;
t152 = qJD(1) * t119;
t139 = t118 * t152;
t143 = t120 * qJD(3);
t101 = t139 - t143;
t151 = qJD(1) * t121;
t109 = -qJD(4) + t151;
t164 = t101 * t109;
t163 = t101 * t119;
t150 = qJD(3) * t118;
t103 = t120 * t152 + t150;
t162 = t103 * t109;
t161 = t109 * t118;
t160 = t109 * t120;
t159 = t118 * t121;
t122 = qJD(3) ^ 2;
t158 = t119 * t122;
t157 = t120 * t121;
t156 = t121 * t122;
t107 = qJD(1) * t111;
t149 = qJD(3) * t119;
t144 = t119 * MDP(16);
t135 = qJD(3) * t144;
t92 = qJD(3) * t95;
t129 = pkin(3) * t119 - pkin(6) * t121;
t105 = t129 * qJD(3);
t99 = qJD(1) * t105;
t134 = t118 * t92 - t120 * t99;
t133 = t103 * t149 - t121 * t88;
t91 = qJD(3) * pkin(6) + t96;
t132 = t109 * t110 + t91;
t131 = t109 * t138;
t130 = t109 * t137;
t87 = t120 * t91 + t170;
t128 = t118 * t91 - t168;
t127 = t118 * t99 + t120 * t92;
t126 = qJD(1) * t114 - t109 * t121;
t125 = 0.2e1 * qJD(3) * t107;
t124 = t126 * t118;
t104 = t129 * qJD(1);
t1 = [0.2e1 * t136 * t175 - 0.2e1 * t141 * t174 + MDP(7) * t156 - MDP(8) * t158 + (-t110 * t156 + t119 * t125) * MDP(10) + (t110 * t158 + t121 * t125) * MDP(11) + (t119 * t120 * t88 + (t121 * t143 - t138) * t103) * MDP(12) + ((-t101 * t120 - t103 * t118) * t148 + (-t172 - t120 * t89 + (t101 * t118 - t103 * t120) * qJD(4)) * t119) * MDP(13) + (t126 * t143 + t131 + t133) * MDP(14) + (t130 + t167 + (-t124 - t163) * qJD(3)) * MDP(15) + (-t109 - t151) * t135 + (-(-t100 * t147 + t105 * t120) * t109 + ((t101 * t110 + t171) * qJD(3) + (t132 * t120 + t170) * qJD(4) + t134) * t121 + (t90 * t146 + t110 * t89 + t166 + (-t110 * t161 + (t100 * t120 - t110 * t159) * qJD(1) - t128) * qJD(3)) * t119) * MDP(17) + ((t100 * t146 + t105 * t118) * t109 + ((t103 * t110 + t169) * qJD(3) + (-t132 * t118 + t168) * qJD(4) + t127) * t121 + (-t90 * t147 + t110 * t88 + t165 + (-t110 * t160 - (t100 * t118 + t110 * t157) * qJD(1) - t87) * qJD(3)) * t119) * MDP(18); (t130 - t167) * MDP(17) + (-t131 + t133) * MDP(18) + t173 * t122 + ((-t124 + t163) * MDP(17) - t126 * MDP(18) * t120) * qJD(3); (-t103 * t160 + t172) * MDP(12) + ((t88 + t164) * t120 + (-t89 + t162) * t118) * MDP(13) + (-pkin(3) * t89 - t165 - t96 * t101 + (pkin(6) * t160 + t171) * qJD(4)) * MDP(17) + (-pkin(3) * t88 + t166 - t96 * t103 + (-pkin(6) * t161 + t169) * qJD(4)) * MDP(18) + (-t146 * MDP(14) + t147 * MDP(15) + (t104 * t120 - t118 * t95) * MDP(17) + (-t104 * t118 - t120 * t95) * MDP(18)) * t109 + ((t109 * t157 + (-t103 + t150) * t119) * MDP(14) + (-t109 * t159 + (t101 + t143) * t119) * MDP(15) + t109 * t144 + (t128 * t119 + (-pkin(6) * t149 - t121 * t90) * t118) * MDP(17) + (-t90 * t157 + (-pkin(6) * t143 + t87) * t119) * MDP(18) + t173 * t107 + (-t121 * t175 + t174) * qJD(1)) * qJD(1); t103 * t101 * MDP(12) + (-t101 ^ 2 + t103 ^ 2) * MDP(13) + (t155 - t164) * MDP(14) + (-t118 * t136 - t162) * MDP(15) + qJD(1) * t135 + (-t103 * t90 - t109 * t87 - t134) * MDP(17) + (t101 * t90 + t109 * t128 - t127) * MDP(18) + (-MDP(14) * t139 - t103 * MDP(15) - t87 * MDP(17) + t128 * MDP(18)) * qJD(4);];
tauc = t1;
