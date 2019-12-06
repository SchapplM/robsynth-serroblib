% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:37
% EndTime: 2019-12-05 16:17:40
% DurationCPUTime: 0.69s
% Computational Cost: add. (419->105), mult. (799->174), div. (0->0), fcn. (401->6), ass. (0->61)
t120 = cos(pkin(9));
t122 = sin(qJ(3));
t124 = cos(qJ(3));
t178 = t124 * MDP(7) + (MDP(8) * t120 + MDP(6)) * t122;
t116 = qJD(2) + qJD(3);
t165 = pkin(2) * qJD(3);
t146 = qJD(2) * t165;
t103 = t116 * qJD(4) + t124 * t146;
t119 = sin(pkin(9));
t114 = t119 ^ 2;
t155 = t120 ^ 2 + t114;
t177 = t155 * t103;
t121 = sin(qJ(5));
t123 = cos(qJ(5));
t161 = t114 * t123;
t176 = t121 * MDP(12) * t161 - MDP(13) * t114 * (t121 ^ 2 - t123 ^ 2);
t171 = t155 * t116;
t166 = pkin(2) * qJD(2);
t133 = -t124 * t166 + qJD(4);
t170 = pkin(2) * t122;
t169 = pkin(2) * t124;
t138 = t122 * t146;
t157 = t120 * t123;
t106 = -t120 * pkin(4) - t119 * pkin(7) - pkin(3);
t91 = t106 * t116 + t133;
t105 = qJ(4) * t116 + t122 * t166;
t95 = qJD(1) * t119 + t105 * t120;
t168 = (t121 * t138 + t103 * t157 + (-t121 * t95 + t123 * t91) * qJD(5)) * t120 + t103 * t161;
t94 = -t120 * qJD(1) + t105 * t119;
t164 = t119 * t94;
t159 = t116 * t120;
t109 = -qJD(5) + t159;
t163 = t109 * t123;
t162 = t114 * t121;
t160 = t116 * t119;
t158 = t120 * t121;
t156 = t120 * t124;
t151 = qJD(5) * t119;
t150 = qJD(5) + t109;
t148 = t122 * t165;
t147 = t94 * t160;
t145 = t116 * t162;
t144 = t116 * t161;
t143 = t123 * t151;
t142 = t121 * t151;
t102 = t142 * t159;
t140 = (t109 * t142 + t102) * MDP(14) + (t109 + t159) * MDP(15) * t143 - 0.2e1 * t176 * qJD(5) * t116;
t139 = t160 * t170;
t130 = -t121 * t91 - t123 * t95;
t88 = qJD(5) * t130 - t103 * t158 + t123 * t138;
t136 = t103 * t162 - t88 * t120 + t94 * t143;
t131 = t120 * t95 + t164;
t129 = t109 * t120 + t114 * t116;
t126 = qJD(4) * t129;
t113 = t116 ^ 2;
t111 = qJ(4) + t170;
t110 = t124 * t165 + qJD(4);
t107 = t119 * t138;
t104 = -t116 * pkin(3) + t133;
t101 = t106 - t169;
t1 = [t102 * MDP(18) + (-t109 * t121 * MDP(18) + (t109 - t159) * MDP(17) * t123) * t151; (qJD(3) * t139 + t107) * MDP(9) + (t110 * t171 + t177) * MDP(10) + (t131 * t110 + t111 * t177 + (t104 + (-pkin(3) - t169) * qJD(2)) * t148) * MDP(11) + (-(-t110 * t158 + t123 * t148) * t109 + t110 * t145 + (-(-t101 * t121 - t111 * t157) * t109 + t111 * t144) * qJD(5) + t136) * MDP(17) + ((t110 * t157 + t121 * t148) * t109 + t110 * t144 + (t101 * t163 + (-t111 * t129 - t164) * t121) * qJD(5) + t168) * MDP(18) + t140 + t178 * (-qJD(2) - t116) * t165; (-qJD(2) * t139 + t107) * MDP(9) + (t133 * t171 + t177) * MDP(10) + (t131 * qJD(4) + qJ(4) * t177 + (-t131 * t124 + (-pkin(3) * qJD(3) - t104) * t122) * t166) * MDP(11) + (t121 * t126 + (-(-qJ(4) * t157 - t106 * t121) * t109 + qJ(4) * t144) * qJD(5) + ((-t121 * t156 + t122 * t123) * t109 - t124 * t145) * t166 + t136) * MDP(17) + (t123 * t126 + (t106 * t163 + (-qJ(4) * t129 - t164) * t121) * qJD(5) + (-(t121 * t122 + t123 * t156) * t109 - t124 * t144) * t166 + t168) * MDP(18) + t140 + t178 * (-qJD(3) + t116) * t166; (-t116 * t131 + t138) * MDP(11) - t155 * MDP(10) * t113 + (t121 * MDP(17) + t123 * MDP(18)) * (-t109 ^ 2 - t114 * t113); (t109 * t130 - t123 * t147 + t88) * MDP(17) + ((-t120 * t103 - t150 * t91) * t123 + (t150 * t95 - t138 + t147) * t121) * MDP(18) - (MDP(14) * t121 + MDP(15) * t123) * t150 * t160 + t176 * t113;];
tauc = t1;
