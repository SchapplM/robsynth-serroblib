% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:12
% EndTime: 2019-12-05 15:22:15
% DurationCPUTime: 0.73s
% Computational Cost: add. (356->125), mult. (1036->210), div. (0->0), fcn. (734->6), ass. (0->74)
t131 = cos(pkin(8));
t156 = qJD(2) * t131;
t124 = -qJD(5) + t156;
t169 = qJD(5) + t124;
t129 = sin(pkin(8));
t168 = pkin(6) * t129;
t128 = sin(pkin(9));
t130 = cos(pkin(9));
t132 = sin(qJ(5));
t133 = cos(qJ(5));
t138 = t128 * t133 + t130 * t132;
t157 = qJD(2) * t129;
t97 = t138 * t157;
t167 = t124 * t97;
t148 = t128 * t157;
t144 = t132 * t148;
t160 = t130 * t133;
t149 = t129 * t160;
t99 = qJD(2) * t149 - t144;
t166 = t124 * t99;
t135 = t138 * qJD(5);
t100 = t129 * t135;
t93 = qJD(2) * t100;
t165 = t131 * t93;
t118 = qJD(5) * t144;
t145 = qJD(5) * t149;
t94 = qJD(2) * t145 - t118;
t164 = t131 * t94;
t163 = qJ(3) * t131;
t162 = t100 * t124;
t161 = t128 * t132;
t119 = -pkin(3) * t131 - qJ(4) * t129 - pkin(2);
t108 = t119 * qJD(2) + qJD(3);
t117 = qJ(3) * t156 + qJD(1) * t129;
t89 = t128 * t108 + t130 * t117;
t159 = t128 * t119 + t130 * t163;
t126 = t129 ^ 2;
t127 = t131 ^ 2;
t158 = t126 + t127;
t155 = qJD(3) * t129;
t154 = qJD(3) * t131;
t153 = qJD(4) * t129;
t152 = qJD(2) * qJD(3);
t151 = t130 * t168;
t150 = 0.2e1 * qJD(3) * t126;
t109 = -t128 * t154 - t130 * t153;
t104 = qJD(2) * t109;
t110 = -t128 * t153 + t130 * t154;
t105 = qJD(2) * t110;
t147 = t133 * t104 - t132 * t105;
t88 = t130 * t108 - t117 * t128;
t116 = -qJ(3) * t157 + qJD(1) * t131;
t86 = (-pkin(4) * t131 - t151) * qJD(2) + t88;
t87 = -pkin(6) * t148 + t89;
t143 = t132 * t87 - t133 * t86;
t142 = -t132 * t86 - t133 * t87;
t114 = qJD(4) - t116;
t141 = -t104 * t130 - t105 * t128;
t140 = t132 * t104 + t133 * t105;
t139 = t116 * t129 - t117 * t131;
t137 = t160 - t161;
t136 = t137 * t124;
t134 = qJD(2) ^ 2;
t123 = t126 * qJ(3) * t152;
t115 = (pkin(4) * t128 + qJ(3)) * t129;
t113 = t130 * t119;
t107 = t137 * t129;
t106 = t138 * t129;
t101 = -qJD(5) * t129 * t161 + t145;
t96 = pkin(4) * t148 + t114;
t92 = t101 * t124;
t91 = -t128 * t168 + t159;
t90 = -t151 + t113 + (-qJ(3) * t128 - pkin(4)) * t131;
t1 = [(t92 - t164) * MDP(18) + (-t162 + t165) * MDP(19) + (-t104 * t128 + t105 * t130 - t131 * t152) * MDP(12) * t129; 0.2e1 * t158 * MDP(7) * t152 + (t123 + (qJ(3) * qJD(2) * t127 - t139) * qJD(3)) * MDP(8) + (-t104 * t131 + (-t109 * t131 + t128 * t150) * qJD(2)) * MDP(9) + (t105 * t131 + (t110 * t131 + t130 * t150) * qJD(2)) * MDP(10) + ((-t109 * t130 - t110 * t128) * qJD(2) + t141) * MDP(11) * t129 + (t105 * t159 + t89 * t110 + t104 * (-t128 * t163 + t113) + t88 * t109 + t123 + t114 * t155) * MDP(12) + (-t100 * t99 - t107 * t93) * MDP(13) + (t100 * t97 - t101 * t99 + t106 * t93 - t107 * t94) * MDP(14) + (t162 + t165) * MDP(15) + (t92 + t164) * MDP(16) + (-(t109 * t133 - t110 * t132) * t124 - t147 * t131 + t115 * t94 + t96 * t101 + (-(-t132 * t90 - t133 * t91) * t124 - t142 * t131) * qJD(5) + (qJD(2) * t106 + t97) * t155) * MDP(18) + ((t109 * t132 + t110 * t133) * t124 + t140 * t131 - t115 * t93 - t96 * t100 + ((-t132 * t91 + t133 * t90) * t124 - t143 * t131) * qJD(5) + (qJD(2) * t107 + t99) * t155) * MDP(19); t139 * qJD(2) * MDP(8) + ((-t114 * t129 + (t128 * t88 - t130 * t89) * t131) * qJD(2) - t141) * MDP(12) + (t124 * t135 + (-t124 * t131 * t138 - t129 * t97) * qJD(2)) * MDP(18) + (qJD(5) * t136 + (-t129 * t99 - t131 * t136) * qJD(2)) * MDP(19) - (t130 * MDP(10) + t128 * MDP(9) + MDP(7)) * t158 * t134; (-t118 - t166) * MDP(18) + MDP(19) * t167 + (-t128 ^ 2 - t130 ^ 2) * t134 * MDP(11) * t126 + ((MDP(10) * t128 - MDP(9) * t130) * t134 * t131 + ((t128 * t89 + t130 * t88 + qJD(3)) * MDP(12) + (MDP(18) * t160 - MDP(19) * t138) * qJD(5)) * qJD(2)) * t129; t99 * t97 * MDP(13) + (-t97 ^ 2 + t99 ^ 2) * MDP(14) + (-t93 - t167) * MDP(15) + (-t94 - t166) * MDP(16) + (t169 * t142 - t96 * t99 + t147) * MDP(18) + (t169 * t143 + t96 * t97 - t140) * MDP(19);];
tauc = t1;
