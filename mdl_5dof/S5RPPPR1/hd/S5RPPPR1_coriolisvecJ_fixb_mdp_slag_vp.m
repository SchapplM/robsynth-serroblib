% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:08
% EndTime: 2019-12-05 17:29:11
% DurationCPUTime: 0.75s
% Computational Cost: add. (419->127), mult. (1145->211), div. (0->0), fcn. (797->8), ass. (0->75)
t137 = cos(pkin(8));
t129 = qJD(1) * t137 - qJD(5);
t175 = qJD(5) + t129;
t134 = sin(pkin(8));
t174 = pkin(6) * t134;
t133 = sin(pkin(9));
t136 = cos(pkin(9));
t139 = sin(qJ(5));
t140 = cos(qJ(5));
t145 = t133 * t140 + t136 * t139;
t142 = t145 * qJD(5);
t107 = t134 * t142;
t99 = qJD(1) * t107;
t173 = t137 * t99;
t121 = -cos(pkin(7)) * pkin(1) - pkin(3) * t137 - qJ(4) * t134 - pkin(2);
t101 = t121 * qJD(1) + qJD(3);
t130 = sin(pkin(7)) * pkin(1) + qJ(3);
t128 = qJD(1) * t130;
t113 = qJD(2) * t134 + t128 * t137;
t92 = t133 * t101 + t136 * t113;
t163 = qJD(1) * t134;
t155 = t133 * t163;
t152 = t139 * t155;
t125 = qJD(5) * t152;
t166 = t136 * t140;
t156 = t134 * t166;
t151 = qJD(1) * t156;
t100 = qJD(5) * t151 - t125;
t172 = t100 * t137;
t104 = t145 * t163;
t171 = t104 * t129;
t106 = t151 - t152;
t170 = t106 * t129;
t169 = t107 * t129;
t168 = t130 * t137;
t167 = t133 * t139;
t165 = t133 * t121 + t136 * t168;
t131 = t134 ^ 2;
t132 = t137 ^ 2;
t164 = t131 + t132;
t162 = qJD(3) * t134;
t161 = qJD(3) * t137;
t160 = qJD(4) * t134;
t159 = qJD(1) * qJD(3);
t158 = t136 * t174;
t157 = 0.2e1 * qJD(3) * t131;
t91 = t136 * t101 - t113 * t133;
t119 = -t133 * t161 - t136 * t160;
t110 = qJD(1) * t119;
t120 = -t133 * t160 + t136 * t161;
t111 = qJD(1) * t120;
t153 = t140 * t110 - t139 * t111;
t112 = qJD(2) * t137 - t134 * t128;
t89 = (-pkin(4) * t137 - t158) * qJD(1) + t91;
t90 = -pkin(6) * t155 + t92;
t150 = t139 * t90 - t140 * t89;
t149 = -t139 * t89 - t140 * t90;
t109 = qJD(4) - t112;
t148 = -t110 * t136 - t111 * t133;
t147 = t139 * t110 + t140 * t111;
t146 = t112 * t134 - t113 * t137;
t144 = t166 - t167;
t143 = t144 * t129;
t141 = qJD(1) ^ 2;
t124 = t131 * t130 * t159;
t118 = (pkin(4) * t133 + t130) * t134;
t117 = t136 * t121;
t115 = t144 * t134;
t114 = t145 * t134;
t108 = (-t134 * t167 + t156) * qJD(5);
t96 = pkin(4) * t155 + t109;
t95 = t108 * t129;
t94 = -t133 * t174 + t165;
t93 = -t158 + t117 + (-t130 * t133 - pkin(4)) * t137;
t1 = [0.2e1 * t164 * MDP(7) * t159 + (t124 + (t132 * t128 - t146) * qJD(3)) * MDP(8) + (-t110 * t137 + (-t119 * t137 + t133 * t157) * qJD(1)) * MDP(9) + (t111 * t137 + (t120 * t137 + t136 * t157) * qJD(1)) * MDP(10) + ((-t119 * t136 - t120 * t133) * qJD(1) + t148) * MDP(11) * t134 + (t111 * t165 + t92 * t120 + t110 * (-t133 * t168 + t117) + t91 * t119 + t124 + t109 * t162) * MDP(12) + (-t106 * t107 - t115 * t99) * MDP(13) + (-t100 * t115 + t104 * t107 - t106 * t108 + t114 * t99) * MDP(14) + (t169 + t173) * MDP(15) + (t95 + t172) * MDP(16) + (-(t119 * t140 - t120 * t139) * t129 - t153 * t137 + t118 * t100 + t96 * t108 + (-(-t139 * t93 - t140 * t94) * t129 - t149 * t137) * qJD(5) + (qJD(1) * t114 + t104) * t162) * MDP(18) + ((t119 * t139 + t120 * t140) * t129 + t147 * t137 - t118 * t99 - t96 * t107 + ((-t139 * t94 + t140 * t93) * t129 - t150 * t137) * qJD(5) + (qJD(1) * t115 + t106) * t162) * MDP(19); (t95 - t172) * MDP(18) + (-t169 + t173) * MDP(19) + (-t110 * t133 + t111 * t136 - t137 * t159) * MDP(12) * t134; t146 * qJD(1) * MDP(8) + ((-t109 * t134 + (t133 * t91 - t136 * t92) * t137) * qJD(1) - t148) * MDP(12) + (t129 * t142 + (-t145 * t129 * t137 - t134 * t104) * qJD(1)) * MDP(18) + (qJD(5) * t143 + (-t134 * t106 - t137 * t143) * qJD(1)) * MDP(19) - (t136 * MDP(10) + t133 * MDP(9) + MDP(7)) * t164 * t141; (-t125 - t170) * MDP(18) + MDP(19) * t171 + (-t133 ^ 2 - t136 ^ 2) * t141 * MDP(11) * t131 + ((MDP(10) * t133 - MDP(9) * t136) * t141 * t137 + ((t133 * t92 + t136 * t91 + qJD(3)) * MDP(12) + (MDP(18) * t166 - t145 * MDP(19)) * qJD(5)) * qJD(1)) * t134; t106 * t104 * MDP(13) + (-t104 ^ 2 + t106 ^ 2) * MDP(14) + (-t99 - t171) * MDP(15) + (-t100 - t170) * MDP(16) + (-t96 * t106 + t175 * t149 + t153) * MDP(18) + (t96 * t104 + t175 * t150 - t147) * MDP(19);];
tauc = t1;
