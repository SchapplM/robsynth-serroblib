% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:23
% EndTime: 2020-01-03 11:52:24
% DurationCPUTime: 0.41s
% Computational Cost: add. (659->94), mult. (1447->133), div. (0->0), fcn. (798->8), ass. (0->62)
t115 = cos(pkin(9)) * pkin(1) + pkin(2);
t114 = t115 * qJD(1);
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t169 = pkin(1) * sin(pkin(9));
t151 = qJD(1) * t169;
t105 = t131 * t114 - t128 * t151;
t127 = sin(qJ(4));
t159 = qJD(4) * t127;
t106 = t128 * t114 + t131 * t151;
t130 = cos(qJ(4));
t162 = t130 * t106;
t172 = -pkin(3) * t159 + t127 * t105 + t162;
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t155 = t129 * MDP(17);
t171 = t126 * MDP(16) + t155;
t170 = (t126 ^ 2 - t129 ^ 2) * MDP(12);
t102 = t105 * qJD(3);
t103 = t106 * qJD(3);
t147 = t127 * t102 + t130 * t103;
t121 = qJD(1) + qJD(3);
t100 = t121 * pkin(3) + t105;
t94 = t127 * t100 + t162;
t86 = t94 * qJD(4) + t147;
t146 = -t127 * t103 - t106 * t159;
t135 = -(qJD(4) * t100 + t102) * t130 - t146;
t120 = qJD(4) + t121;
t168 = t120 * pkin(4);
t154 = t129 * qJD(5);
t164 = t127 * t106;
t93 = t130 * t100 - t164;
t91 = -t93 - t168;
t167 = t86 * t126 + t91 * t154;
t136 = t131 * t115 - t128 * t169;
t107 = t136 * qJD(3);
t110 = t128 * t115 + t131 * t169;
t108 = t110 * qJD(3);
t109 = pkin(3) + t136;
t137 = t127 * t109 + t130 * t110;
t166 = (t137 * qJD(4) + t127 * t107 + t130 * t108) * t120;
t165 = t94 * t120;
t132 = qJD(5) ^ 2;
t163 = t129 * t132;
t160 = MDP(11) * t129;
t158 = qJD(4) * t130;
t156 = t126 * qJD(5);
t153 = t132 * MDP(14);
t149 = MDP(13) * t163 + 0.2e1 * (t160 * t126 - t170) * qJD(5) * t120;
t148 = -t91 * t120 + t135;
t142 = pkin(8) * t132 - t165;
t141 = qJD(5) * (t93 - t168);
t140 = t132 * (pkin(8) + t137) + t166;
t138 = t130 * t109 - t127 * t110;
t87 = t138 * qJD(4) + t130 * t107 - t127 * t108;
t139 = qJD(5) * (t120 * (-pkin(4) - t138) - t87);
t89 = t91 * t156;
t134 = t89 * MDP(16) + t167 * MDP(17) + t149;
t119 = t120 ^ 2;
t118 = -t130 * pkin(3) - pkin(4);
t117 = t127 * pkin(3) + pkin(8);
t1 = [(-t108 * t121 - t103) * MDP(6) + (-t107 * t121 - t102) * MDP(7) + (-t86 - t166) * MDP(9) + (-t87 * t120 + t135) * MDP(10) + ((-t140 - t86) * MDP(16) + MDP(17) * t139) * t129 + (MDP(16) * t139 + t140 * MDP(17) - t153) * t126 + t134; -t171 * t132; (t106 * t121 - t103) * MDP(6) + (t105 * t121 - t102) * MDP(7) + (-t100 * t159 - t106 * t158 - t147) * MDP(9) + (-t100 * t158 - t130 * t102 - t146) * MDP(10) - t126 * t153 + (-t117 * t163 - t86 * t129 + t89) * MDP(16) + (t132 * t126 * t117 + t167) * MDP(17) + (t172 * MDP(9) + (t118 * t156 + t172 * t129) * MDP(16) + (t118 * t154 - t172 * t126) * MDP(17)) * t120 + t149 + (t120 * MDP(10) + t171 * qJD(5)) * (-pkin(3) * t158 + t130 * t105 - t164); (-t86 + t165) * MDP(9) + (t93 * t120 + t135) * MDP(10) + ((-t142 - t86) * MDP(16) + MDP(17) * t141) * t129 + (MDP(16) * t141 + t142 * MDP(17) - t153) * t126 + t134; t148 * t155 + t119 * t170 + (t148 * MDP(16) - t119 * t160) * t126;];
tauc = t1;
