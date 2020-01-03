% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:15
% EndTime: 2019-12-31 18:01:16
% DurationCPUTime: 0.37s
% Computational Cost: add. (473->81), mult. (825->131), div. (0->0), fcn. (432->6), ass. (0->55)
t114 = -qJD(1) + qJD(4);
t117 = sin(pkin(8));
t118 = cos(pkin(8));
t120 = sin(qJ(4));
t122 = cos(qJ(4));
t130 = t122 * t117 + t120 * t118;
t166 = t114 * t130;
t123 = -pkin(1) - pkin(2);
t111 = t123 * qJD(1) + qJD(2);
t107 = t118 * t111;
t154 = t117 * qJ(2);
t95 = t107 + (-pkin(3) - t154) * qJD(1);
t159 = t120 * t95;
t142 = qJD(1) * qJD(2);
t137 = t122 * t142;
t138 = t118 * t142;
t149 = qJD(4) * t122;
t150 = qJD(1) * qJ(2);
t97 = t117 * t111 + t118 * t150;
t83 = qJD(4) * t159 + t117 * t137 + t120 * t138 + t97 * t149;
t165 = -t83 + (t122 * t97 + t159) * t114;
t135 = t118 * t123 - t154;
t105 = -pkin(3) + t135;
t106 = t118 * qJ(2) + t117 * t123;
t131 = t120 * t105 + t122 * t106;
t164 = t83 - (t130 * qJD(2) + t131 * qJD(4)) * t114;
t163 = qJ(2) * MDP(6) + MDP(5);
t119 = sin(qJ(5));
t121 = cos(qJ(5));
t145 = t121 * MDP(19);
t162 = t119 * MDP(18) + t145;
t161 = MDP(14) * (t119 ^ 2 - t121 ^ 2);
t139 = t117 * t142;
t82 = -(qJD(4) * t97 + t139) * t120 + t118 * t137 + t95 * t149;
t160 = t114 * pkin(4);
t153 = t114 * qJD(5) * t161;
t147 = t119 * qJD(5);
t146 = t121 * MDP(13);
t124 = qJD(5) ^ 2;
t144 = t124 * MDP(15);
t143 = t124 * MDP(16);
t140 = t114 * t146;
t88 = -t120 * t97 + t122 * t95;
t84 = -t88 - t160;
t136 = -t84 * t114 - t82;
t134 = t84 + t88 - t160;
t133 = (-t117 * t150 + t107) * t117 - t97 * t118;
t132 = t122 * t105 - t120 * t106;
t103 = t120 * t117 - t122 * t118;
t129 = pkin(7) * t124 - t165;
t128 = t124 * (-pkin(7) + t131) - t164;
t86 = -t103 * qJD(2) + t132 * qJD(4);
t127 = qJD(5) * (t114 * (pkin(4) - t132) - t84 - t86);
t113 = t114 ^ 2;
t1 = [0.2e1 * MDP(7) * t139 + 0.2e1 * MDP(8) * t138 + ((t118 * t106 - t117 * t135) * qJD(1) - t133) * qJD(2) * MDP(9) + t164 * MDP(11) + (-t86 * t114 + t82) * MDP(12) - 0.2e1 * t140 * t147 + 0.2e1 * t153 - t121 * t144 + t119 * t143 + (t119 * t127 - t128 * t121) * MDP(18) + (t128 * t119 + t121 * t127) * MDP(19) + 0.2e1 * t163 * t142; (-t121 * MDP(18) + t119 * MDP(19)) * t124 * t130 + (-t166 * MDP(11) + (t103 * t147 - t166 * t121) * MDP(18) + (qJD(5) * t103 * t121 + t166 * t119) * MDP(19)) * t114 + (t114 * MDP(12) + t162 * qJD(5)) * t114 * t103 + (t133 * MDP(9) + (-t117 * MDP(7) - t118 * MDP(8) - t163) * qJD(1)) * qJD(1); -t162 * t124; t165 * MDP(11) + (t88 * t114 - t82) * MDP(12) - 0.2e1 * t153 + (t134 * MDP(19) * qJD(5) - t129 * MDP(18) + t144) * t121 + (-t143 + t129 * MDP(19) + (t134 * MDP(18) + 0.2e1 * t140) * qJD(5)) * t119; t136 * t145 + t113 * t161 + (t136 * MDP(18) - t113 * t146) * t119;];
tauc = t1;
