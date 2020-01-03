% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:25
% EndTime: 2019-12-31 17:49:27
% DurationCPUTime: 0.57s
% Computational Cost: add. (709->139), mult. (1688->186), div. (0->0), fcn. (1134->6), ass. (0->66)
t148 = sin(pkin(8));
t150 = cos(pkin(8));
t182 = qJD(1) * (t148 ^ 2 + t150 ^ 2);
t142 = sin(pkin(7)) * pkin(1) + qJ(3);
t135 = qJD(1) * t142;
t118 = t148 * qJD(2) + t150 * t135;
t145 = t150 * qJD(2);
t181 = ((-t135 * t148 + t145) * t148 - t118 * t150) * MDP(8);
t152 = sin(qJ(4));
t153 = cos(qJ(4));
t133 = t148 * t153 + t150 * t152;
t125 = t133 * qJD(1);
t168 = MDP(14) + MDP(16);
t180 = t125 ^ 2;
t179 = pkin(6) + t142;
t178 = pkin(6) * qJD(1);
t116 = t150 * t178 + t118;
t177 = t116 * t152;
t115 = t145 + (-t135 - t178) * t148;
t98 = t115 * t153 - t177;
t176 = qJD(5) - t98;
t128 = t133 * qJD(4);
t120 = qJD(1) * t128;
t172 = qJD(1) * t153;
t164 = t150 * t172;
t173 = qJD(1) * t152;
t166 = t148 * t173;
t123 = -t164 + t166;
t171 = qJD(4) * t152;
t165 = t148 * t171;
t170 = qJD(4) * t153;
t127 = -t150 * t170 + t165;
t175 = -t133 * t120 + t127 * t123;
t169 = qJD(1) * qJD(3);
t167 = MDP(15) - MDP(18);
t163 = t152 * t169;
t162 = t153 * t169;
t161 = t98 + t177;
t94 = t115 * t171 + t116 * t170 + t148 * t162 + t150 * t163;
t139 = qJD(4) * t164;
t119 = qJD(1) * t165 - t139;
t160 = pkin(4) * t120 + qJ(5) * t119;
t99 = t115 * t152 + t116 * t153;
t132 = t148 * t152 - t153 * t150;
t158 = -t119 * t132 + t125 * t128;
t129 = t179 * t148;
t130 = t179 * t150;
t157 = -t129 * t153 - t130 * t152;
t106 = -t129 * t152 + t130 * t153;
t134 = -cos(pkin(7)) * pkin(1) - pkin(3) * t150 - pkin(2);
t121 = t134 * qJD(1) + qJD(3);
t102 = pkin(4) * t123 - qJ(5) * t125 + t121;
t156 = -t102 * t125 - t94;
t155 = -t115 * t170 + t148 * t163 - t150 * t162;
t122 = t123 ^ 2;
t109 = pkin(4) * t125 + qJ(5) * t123;
t108 = t139 + (t123 - t166) * qJD(4);
t104 = pkin(4) * t132 - qJ(5) * t133 + t134;
t103 = pkin(4) * t128 + qJ(5) * t127 - qJD(5) * t133;
t101 = t133 * qJD(3) + t106 * qJD(4);
t100 = -t132 * qJD(3) + t157 * qJD(4);
t97 = -qJD(5) * t125 + t160;
t96 = qJD(4) * qJ(5) + t99;
t95 = -qJD(4) * pkin(4) + t176;
t93 = (qJD(5) - t177) * qJD(4) - t155;
t1 = [(-t119 * t133 - t125 * t127) * MDP(9) + (-t158 + t175) * MDP(10) + (t120 * t134 + t121 * t128) * MDP(14) + (-t119 * t134 - t121 * t127) * MDP(15) + (t102 * t128 + t103 * t123 + t104 * t120 + t97 * t132) * MDP(16) + (-t100 * t123 + t101 * t125 - t106 * t120 + t119 * t157 - t127 * t95 - t128 * t96 - t132 * t93 + t133 * t94) * MDP(17) + (t102 * t127 - t103 * t125 + t104 * t119 - t133 * t97) * MDP(18) + (t100 * t96 + t101 * t95 + t102 * t103 + t104 * t97 + t106 * t93 - t157 * t94) * MDP(19) + (-MDP(11) * t127 - MDP(12) * t128 - t167 * t100 - t168 * t101) * qJD(4) + (-t181 + (t142 * MDP(8) + (2 * MDP(7))) * t182) * qJD(3); (t158 + t175) * MDP(17) + (-t127 * t96 + t128 * t95 + t132 * t94 + t133 * t93) * MDP(19) + (t167 * t127 - t168 * t128) * qJD(4); (-t122 - t180) * MDP(17) + (t123 * t96 + (-qJD(5) - t95) * t125 + t160) * MDP(19) - t167 * (-t139 + (t123 + t166) * qJD(4)) + 0.2e1 * t168 * t125 * qJD(4) + (-MDP(7) * t182 + t181) * qJD(1); (-t122 + t180) * MDP(10) + t108 * MDP(11) + (-t121 * t125 - t94) * MDP(14) + t155 * MDP(15) + t156 * MDP(16) + (pkin(4) * t119 - qJ(5) * t120 - (-t96 + t99) * t125) * MDP(17) + (t109 * t125 - t155) * MDP(18) + (-pkin(4) * t94 + qJ(5) * t93 - t102 * t109 + t176 * t96 - t95 * t99) * MDP(19) + (t125 * MDP(9) + t121 * MDP(15) - t109 * MDP(16) + (t95 - t176) * MDP(17) - t102 * MDP(18)) * t123 + ((-t148 * t172 - t150 * t173 + t125) * MDP(12) + t161 * MDP(15) + (0.2e1 * qJD(5) - t161) * MDP(18) + t168 * t99) * qJD(4); t125 * t123 * MDP(16) + t108 * MDP(17) + (-qJD(4) ^ 2 - t180) * MDP(18) + (-qJD(4) * t96 - t156) * MDP(19);];
tauc = t1;
