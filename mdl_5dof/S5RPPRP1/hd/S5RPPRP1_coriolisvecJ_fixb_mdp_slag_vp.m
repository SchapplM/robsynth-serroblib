% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP1
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:27
% EndTime: 2019-12-05 17:36:30
% DurationCPUTime: 0.67s
% Computational Cost: add. (519->113), mult. (1207->187), div. (0->0), fcn. (705->6), ass. (0->66)
t119 = cos(pkin(8));
t149 = qJD(1) * t119;
t164 = qJD(4) - t149;
t110 = sin(pkin(7)) * pkin(1) + qJ(3);
t121 = sin(qJ(4));
t122 = cos(qJ(4));
t117 = sin(pkin(8));
t102 = -cos(pkin(7)) * pkin(1) - pkin(3) * t119 - pkin(6) * t117 - pkin(2);
t132 = qJ(5) * t117 - t102;
t163 = -t110 * t119 * t122 + t121 * t132;
t113 = t117 ^ 2;
t115 = t121 ^ 2;
t116 = t122 ^ 2;
t162 = (-t122 * t121 * MDP(9) + MDP(10) * (t115 - t116)) * t113;
t114 = t119 ^ 2;
t104 = t110 * qJD(1);
t100 = qJD(2) * t117 + t104 * t119;
t150 = qJD(1) * t117;
t137 = t122 * t150;
t97 = qJD(1) * t102 + qJD(3);
t96 = t122 * t97;
t90 = -qJ(5) * t137 - t100 * t121 + t96;
t87 = pkin(4) * t164 + t90;
t161 = t87 - t90;
t146 = qJD(3) * t122;
t108 = t119 * t146;
t144 = qJD(4) * t122;
t160 = -qJD(1) * t108 - t97 * t144;
t159 = t102 * t144 + t108;
t158 = t164 * t119;
t123 = qJD(1) ^ 2;
t157 = t113 * t123;
t156 = t117 * t122;
t155 = t119 * t121;
t140 = qJD(1) * qJD(4);
t133 = t117 * t140;
t101 = t122 * pkin(4) * t133 + qJD(3) * t150;
t154 = t113 + t114;
t151 = qJD(1) * t113;
t148 = qJD(1) * t121;
t147 = qJD(3) * t121;
t145 = qJD(4) * t121;
t143 = qJD(5) * t121;
t142 = qJD(5) * t122;
t141 = qJD(4) - t164;
t138 = t117 * t148;
t136 = t119 * t147;
t135 = t110 * t145;
t134 = t117 * t145;
t131 = qJD(1) * t154;
t112 = t119 * qJD(2);
t99 = t104 * t117 - t112;
t128 = t100 * t119 + t117 * t99;
t127 = -t100 * t122 - t121 * t97;
t126 = -t100 * t145 - t160;
t125 = t127 * qJD(4);
t124 = -t164 ^ 2 - t157;
t103 = t133 * t155;
t94 = qJD(5) - t112 + (pkin(4) * t148 + t104) * t117;
t92 = -t132 * t122 + (-t110 * t121 - pkin(4)) * t119;
t91 = -qJ(5) * t138 - t127;
t89 = t163 * qJD(4) - t117 * t142 - t136;
t88 = -t117 * t143 + (-qJ(5) * t156 - t110 * t155) * qJD(4) + t159;
t86 = t125 + (-t136 + (qJ(5) * t145 - t142) * t117) * qJD(1);
t85 = (-qJ(5) * t144 - t143) * t150 + t126;
t1 = [(-t134 * t164 + t103) * MDP(11) + ((-t158 + (0.2e1 * t113 + t114) * qJD(1)) * t147 + ((-t102 * t164 + t119 * t97) * t121 + ((t151 - t158) * t110 + t128) * t122) * qJD(4)) * MDP(14) + (-(-t119 * t135 + t159) * t164 + t126 * t119 - t99 * t134 + (-t135 + 0.2e1 * t146) * t151) * MDP(15) + (-t163 * t85 + t86 * t92 + t87 * t89 + t91 * t88) * MDP(17) + 0.2e1 * t162 * t140 + (0.2e1 * MDP(7) * t131 + (t110 * t131 + t128) * MDP(8)) * qJD(3) + ((-t164 + t149) * MDP(12) * t144 + (-t121 * t85 - t122 * t86 + (t121 * t87 - t122 * t91) * qJD(4) + (-t121 * t88 - t122 * t89 + (t121 * t92 + t122 * t163) * qJD(4)) * qJD(1)) * MDP(16) + (t101 * (pkin(4) * t121 + t110) + t94 * (pkin(4) * t144 + qJD(3))) * MDP(17)) * t117; -t101 * t119 * MDP(17) + t103 * MDP(15) + ((-t121 * t86 + t122 * t85) * MDP(17) + ((MDP(15) * t164 - t91 * MDP(17)) * t121 + ((-t164 - t149) * MDP(14) - t87 * MDP(17)) * t122) * qJD(4)) * t117; -t154 * MDP(7) * t123 + (-t94 * t117 * MDP(17) - t128 * MDP(8)) * qJD(1) + (t124 * MDP(15) + (t164 * t91 + t86) * MDP(17)) * t122 + (t124 * MDP(14) + (-t164 * t87 + t85) * MDP(17)) * t121; (-t127 * t164 + t125 + (-t156 * t99 - t136) * qJD(1)) * MDP(14) + (t96 * t164 + (t100 * t141 + t150 * t99) * t121 + t160) * MDP(15) + (pkin(4) * qJD(4) - t161) * MDP(16) * t138 + (t161 * t91 + (-t137 * t94 + t86) * pkin(4)) * MDP(17) - (MDP(11) * t121 + MDP(12) * t122) * t141 * t150 - t162 * t123; ((t121 * t91 + t122 * t87) * t150 + t101) * MDP(17) + (-t115 - t116) * MDP(16) * t157;];
tauc = t1;
