% Calculate Coriolis joint torque vector for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:46
% EndTime: 2022-01-20 09:12:49
% DurationCPUTime: 0.75s
% Computational Cost: add. (409->122), mult. (1103->202), div. (0->0), fcn. (769->8), ass. (0->74)
t136 = cos(pkin(8));
t128 = qJD(1) * t136 - qJD(5);
t173 = qJD(5) + t128;
t133 = sin(pkin(8));
t172 = pkin(6) * t133;
t132 = sin(pkin(9));
t135 = cos(pkin(9));
t138 = sin(qJ(5));
t139 = cos(qJ(5));
t144 = t132 * t139 + t135 * t138;
t141 = t144 * qJD(5);
t106 = t133 * t141;
t98 = qJD(1) * t106;
t171 = t136 * t98;
t161 = qJD(1) * t133;
t153 = t132 * t161;
t149 = t138 * t153;
t124 = qJD(5) * t149;
t164 = t135 * t139;
t154 = t133 * t164;
t150 = qJD(5) * t154;
t99 = qJD(1) * t150 - t124;
t170 = t136 * t99;
t120 = -cos(pkin(7)) * pkin(1) - pkin(3) * t136 - qJ(4) * t133 - pkin(2);
t100 = qJD(1) * t120 + qJD(3);
t129 = sin(pkin(7)) * pkin(1) + qJ(3);
t127 = qJD(1) * t129;
t112 = qJD(2) * t133 + t127 * t136;
t91 = t132 * t100 + t135 * t112;
t103 = t144 * t161;
t169 = t103 * t128;
t105 = qJD(1) * t154 - t149;
t168 = t105 * t128;
t167 = t106 * t128;
t166 = t129 * t136;
t165 = t132 * t138;
t163 = t132 * t120 + t135 * t166;
t130 = t133 ^ 2;
t131 = t136 ^ 2;
t162 = t130 + t131;
t160 = qJD(3) * t136;
t159 = qJD(4) * t133;
t158 = t133 * qJD(3);
t157 = qJD(1) * qJD(3);
t156 = t135 * t172;
t155 = 0.2e1 * qJD(3) * t130;
t90 = t135 * t100 - t112 * t132;
t118 = -t132 * t160 - t135 * t159;
t109 = qJD(1) * t118;
t119 = -t132 * t159 + t135 * t160;
t110 = qJD(1) * t119;
t151 = t139 * t109 - t138 * t110;
t111 = qJD(2) * t136 - t133 * t127;
t88 = (-pkin(4) * t136 - t156) * qJD(1) + t90;
t89 = -pkin(6) * t153 + t91;
t148 = t138 * t89 - t139 * t88;
t147 = -t138 * t88 - t139 * t89;
t108 = qJD(4) - t111;
t146 = t138 * t109 + t139 * t110;
t145 = t111 * t133 - t112 * t136;
t143 = t164 - t165;
t142 = t143 * t128;
t140 = qJD(1) ^ 2;
t123 = t130 * t129 * t157;
t117 = (pkin(4) * t132 + t129) * t133;
t116 = t135 * t120;
t114 = t143 * t133;
t113 = t144 * t133;
t107 = -qJD(5) * t133 * t165 + t150;
t95 = pkin(4) * t153 + t108;
t94 = t107 * t128;
t93 = -t132 * t172 + t163;
t92 = -t156 + t116 + (-t129 * t132 - pkin(4)) * t136;
t1 = [0.2e1 * t162 * MDP(6) * t157 + (t123 + (t127 * t131 - t145) * qJD(3)) * MDP(7) + (-t109 * t136 + (-t118 * t136 + t132 * t155) * qJD(1)) * MDP(8) + (t110 * t136 + (t119 * t136 + t135 * t155) * qJD(1)) * MDP(9) + (t110 * t163 + t91 * t119 + t109 * (-t132 * t166 + t116) + t90 * t118 + t123 + t108 * t158) * MDP(10) + (-t105 * t106 - t114 * t98) * MDP(11) + (t103 * t106 - t105 * t107 + t113 * t98 - t114 * t99) * MDP(12) + (t167 + t171) * MDP(13) + (t94 + t170) * MDP(14) + (-(t118 * t139 - t119 * t138) * t128 - t151 * t136 + t117 * t99 + t95 * t107 + (-(-t138 * t92 - t139 * t93) * t128 - t147 * t136) * qJD(5) + (qJD(1) * t113 + t103) * t158) * MDP(16) + ((t118 * t138 + t119 * t139) * t128 + t146 * t136 - t117 * t98 - t95 * t106 + ((-t138 * t93 + t139 * t92) * t128 - t148 * t136) * qJD(5) + (qJD(1) * t114 + t105) * t158) * MDP(17); (t94 - t170) * MDP(16) + (-t167 + t171) * MDP(17) + (-t109 * t132 + t110 * t135 - t136 * t157) * MDP(10) * t133; t145 * qJD(1) * MDP(7) + (t109 * t135 + t110 * t132 + (-t108 * t133 + (t132 * t90 - t135 * t91) * t136) * qJD(1)) * MDP(10) + (t128 * t141 + (-t128 * t136 * t144 - t133 * t103) * qJD(1)) * MDP(16) + (qJD(5) * t142 + (-t133 * t105 - t136 * t142) * qJD(1)) * MDP(17) - (MDP(8) * t132 + MDP(9) * t135 + MDP(6)) * t162 * t140; (-t124 - t168) * MDP(16) + MDP(17) * t169 + ((-MDP(8) * t135 + MDP(9) * t132) * t140 * t136 + ((t132 * t91 + t135 * t90 + qJD(3)) * MDP(10) + (MDP(16) * t164 - MDP(17) * t144) * qJD(5)) * qJD(1)) * t133; t105 * t103 * MDP(11) + (-t103 ^ 2 + t105 ^ 2) * MDP(12) + (-t98 - t169) * MDP(13) + (-t99 - t168) * MDP(14) + (-t95 * t105 + t173 * t147 + t151) * MDP(16) + (t95 * t103 + t173 * t148 - t146) * MDP(17);];
tauc = t1;
