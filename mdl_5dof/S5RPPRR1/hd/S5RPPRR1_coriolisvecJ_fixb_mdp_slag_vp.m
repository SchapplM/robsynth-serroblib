% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:17
% EndTime: 2019-12-05 17:38:20
% DurationCPUTime: 0.66s
% Computational Cost: add. (380->125), mult. (810->179), div. (0->0), fcn. (450->4), ass. (0->66)
t132 = sin(qJ(5));
t133 = sin(qJ(4));
t154 = qJD(5) * t132;
t156 = qJD(4) * t133;
t175 = t132 * t156 + t133 * t154;
t135 = cos(qJ(4));
t153 = t133 * MDP(16);
t174 = t135 * MDP(15) - t153;
t161 = MDP(16) * t135;
t162 = MDP(15) * t133;
t173 = t161 + t162;
t172 = (t133 ^ 2 - t135 ^ 2) * MDP(11);
t126 = qJD(4) + qJD(5);
t171 = qJD(5) - t126;
t131 = pkin(1) + qJ(3);
t130 = -pkin(6) + qJ(2);
t170 = pkin(7) - t130;
t134 = cos(qJ(5));
t159 = qJD(1) * t135;
t160 = qJD(1) * t133;
t104 = t132 * t160 - t134 * t159;
t169 = t104 * t126;
t108 = t132 * t135 + t133 * t134;
t105 = t108 * qJD(1);
t168 = t105 * t126;
t125 = qJD(1) * qJ(2) + qJD(3);
t116 = -pkin(6) * qJD(1) + t125;
t102 = -pkin(7) * t160 + t116 * t133;
t167 = t134 * t102;
t166 = t134 * t135;
t165 = t175 * qJD(1);
t158 = qJD(2) * t133;
t157 = qJD(2) * t135;
t155 = qJD(4) * t135;
t117 = qJD(1) * t131 - qJD(2);
t151 = qJD(2) - t117;
t149 = pkin(4) * t159;
t146 = t135 * t133 * MDP(10);
t103 = -pkin(7) * t159 + t135 * t116;
t99 = qJD(4) * pkin(4) + t103;
t145 = -pkin(4) * t126 - t99;
t113 = t170 * t135;
t144 = t134 * t126;
t143 = MDP(23) * t126;
t122 = pkin(4) * t133 + t131;
t118 = pkin(4) * t155 + qJD(3);
t142 = t126 * t166;
t107 = qJD(1) * t122 - qJD(2);
t123 = qJD(1) * t157;
t96 = t123 + (pkin(7) * qJD(1) - t116) * t156;
t97 = t116 * t155 + (-pkin(7) * t155 + t158) * qJD(1);
t140 = t107 * t104 - t132 * t97 + t134 * t96;
t94 = t126 * t108;
t91 = t94 * qJD(1);
t139 = -t104 * t105 * MDP(17) + (-t91 + t168) * MDP(19) + (-t144 * t159 + t165 - t169) * MDP(20) + (t104 ^ 2 - t105 ^ 2) * MDP(18);
t138 = t102 * t154 + (-t102 * t126 - t96) * t132 + t107 * t105;
t137 = qJD(1) ^ 2;
t136 = qJD(4) ^ 2;
t112 = t170 * t133;
t111 = t118 * qJD(1);
t109 = -t132 * t133 + t166;
t101 = -qJD(4) * t113 + t158;
t100 = t156 * t170 + t157;
t95 = t142 - t175;
t92 = qJD(1) * t142 - t165;
t1 = [(qJD(2) * t125 + qJD(3) * t117) * MDP(9) + (t104 * t94 - t109 * t91) * MDP(17) + (t104 * t95 + t105 * t94 + t108 * t91 - t109 * t92) * MDP(18) + (t118 * t105 + t107 * t95 + t111 * t108 + t122 * t92) * MDP(22) + (-t118 * t104 - t107 * t94 + t111 * t109 - t122 * t91) * MDP(23) + (-t94 * MDP(19) - t95 * MDP(20) + (t100 * t134 - t101 * t132) * MDP(22) + (-t100 * t132 - t101 * t134) * MDP(23) + ((t112 * t134 + t113 * t132) * MDP(22) + (-t112 * t132 + t113 * t134) * MDP(23)) * qJD(5)) * t126 + (-t133 * MDP(12) - t135 * MDP(13) - t173 * t130) * t136 + ((t131 * MDP(9) + (2 * MDP(8)) + 0.2e1 * t161 + 0.2e1 * t162) * qJD(3) + (t174 * t131 - 0.2e1 * t146 + 0.2e1 * t172) * qJD(4) + ((2 * MDP(5)) + (2 * MDP(7)) + ((2 * MDP(6)) + MDP(9)) * qJ(2)) * qJD(2)) * qJD(1) + t174 * qJD(4) * (qJD(2) + t117); (t165 + t169) * MDP(22) + MDP(23) * t168 + (-qJ(2) * MDP(6) - MDP(5) - MDP(7)) * t137 + ((-qJD(3) - t125) * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + t134 * t143) * t133 + (-0.2e1 * qJD(4) * MDP(15) - MDP(22) * t144 + t132 * t143) * t135) * qJD(1); -t137 * MDP(8) + (-MDP(22) * t94 - MDP(23) * t95) * t126 + (-t105 * MDP(22) + t104 * MDP(23) + MDP(9) * t151) * qJD(1) + t173 * (-t136 - t137); (-t117 * t159 + t123) * MDP(15) - t151 * qJD(1) * t153 + (-t105 * t149 - (-t103 * t132 - t167) * t126 + (t132 * t145 - t167) * qJD(5) + t140) * MDP(22) + (t104 * t149 + (qJD(5) * t145 + t103 * t126 - t97) * t134 + t138) * MDP(23) + t139 + (t146 - t172) * t137; (t140 + t171 * (-t132 * t99 - t167)) * MDP(22) + ((-t171 * t99 - t97) * t134 + t138) * MDP(23) + t139;];
tauc = t1;
