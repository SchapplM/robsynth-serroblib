% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:14
% EndTime: 2019-12-05 15:45:16
% DurationCPUTime: 0.35s
% Computational Cost: add. (443->92), mult. (1013->149), div. (0->0), fcn. (746->8), ass. (0->62)
t126 = sin(pkin(9));
t127 = cos(pkin(9));
t130 = sin(qJ(2));
t133 = cos(qJ(2));
t116 = t126 * t133 + t127 * t130;
t112 = t116 * qJD(1);
t115 = -t126 * t130 + t127 * t133;
t114 = t115 * qJD(1);
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t120 = t127 * pkin(2) + pkin(3);
t163 = pkin(2) * t126;
t137 = t129 * t120 + t132 * t163;
t167 = t137 * qJD(4) - t132 * t112 - t129 * t114;
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t149 = t131 * MDP(15);
t166 = t128 * MDP(14) + t149;
t165 = (t128 ^ 2 - t131 ^ 2) * MDP(10);
t111 = t116 * qJD(2);
t107 = qJD(1) * t111;
t113 = t115 * qJD(2);
t108 = qJD(1) * t113;
t142 = t132 * t107 + t129 * t108;
t119 = qJD(2) * pkin(2) + t133 * qJD(1);
t154 = qJD(1) * t130;
t103 = t127 * t119 - t126 * t154;
t100 = qJD(2) * pkin(3) + t103;
t104 = t126 * t119 + t127 * t154;
t94 = t129 * t100 + t132 * t104;
t86 = t94 * qJD(4) + t142;
t153 = qJD(4) * t129;
t143 = -t104 * t153 - t129 * t107;
t164 = -(qJD(4) * t100 + t108) * t132 - t143;
t123 = qJD(2) + qJD(4);
t162 = t123 * pkin(4);
t148 = t131 * qJD(5);
t93 = t132 * t100 - t129 * t104;
t91 = -t93 - t162;
t161 = t86 * t128 + t91 * t148;
t160 = t94 * t123;
t158 = MDP(9) * t131;
t134 = qJD(5) ^ 2;
t157 = t128 * t134;
t156 = t131 * t134;
t152 = qJD(4) * t132;
t150 = t128 * qJD(5);
t147 = t134 * MDP(12);
t145 = MDP(11) * t156 + 0.2e1 * (t158 * t128 - t165) * qJD(5) * t123;
t144 = -t91 * t123 + t164;
t140 = pkin(7) * t134 - t160;
t139 = qJD(5) * (t93 - t162);
t138 = t132 * t115 - t129 * t116;
t98 = t129 * t115 + t132 * t116;
t136 = t132 * t120 - t129 * t163;
t122 = t123 ^ 2;
t110 = pkin(7) + t137;
t109 = -pkin(4) - t136;
t89 = t91 * t150;
t88 = t98 * qJD(4) + t132 * t111 + t129 * t113;
t87 = t138 * qJD(4) - t129 * t111 + t132 * t113;
t1 = [(-t103 * t111 + t104 * t113 - t107 * t115 + t108 * t116) * MDP(5) + (-t87 * t150 - t98 * t156) * MDP(14) + (-t87 * t148 + t98 * t157) * MDP(15) + (-t130 * MDP(3) - t133 * MDP(4)) * qJD(2) ^ 2 + (-t88 * MDP(7) - t87 * MDP(8) + (-t131 * t88 - t138 * t150) * MDP(14) + (t128 * t88 - t138 * t148) * MDP(15)) * t123; (t103 * t112 - t104 * t114 + (-t107 * t127 + t108 * t126) * pkin(2)) * MDP(5) + (-t100 * t153 - t104 * t152 - t142) * MDP(7) + (-t100 * t152 - t132 * t108 - t143) * MDP(8) - t128 * t147 + (-t110 * t156 - t86 * t131 + t89) * MDP(14) + (t110 * t157 + t161) * MDP(15) + (-t167 * MDP(7) + (t109 * t150 - t167 * t131) * MDP(14) + (t109 * t148 + t167 * t128) * MDP(15)) * t123 + t145 + (t123 * MDP(8) + t166 * qJD(5)) * (-t136 * qJD(4) - t129 * t112 + t132 * t114); -t166 * t134; (t160 - t86) * MDP(7) + (t93 * t123 + t164) * MDP(8) + t89 * MDP(14) + t161 * MDP(15) + ((-t140 - t86) * MDP(14) + MDP(15) * t139) * t131 + (MDP(14) * t139 + t140 * MDP(15) - t147) * t128 + t145; t144 * t149 + t122 * t165 + (t144 * MDP(14) - t122 * t158) * t128;];
tauc = t1;
