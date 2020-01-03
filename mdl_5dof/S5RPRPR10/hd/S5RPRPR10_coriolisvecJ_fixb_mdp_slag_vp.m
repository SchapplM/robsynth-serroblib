% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:06
% EndTime: 2019-12-31 18:26:07
% DurationCPUTime: 0.40s
% Computational Cost: add. (585->96), mult. (998->146), div. (0->0), fcn. (494->6), ass. (0->63)
t149 = qJD(1) - qJD(3);
t128 = sin(pkin(8));
t129 = cos(pkin(8));
t131 = sin(qJ(3));
t133 = cos(qJ(3));
t140 = t128 * t133 + t129 * t131;
t172 = t149 * t140;
t130 = sin(qJ(5));
t132 = cos(qJ(5));
t154 = t132 * MDP(17);
t171 = t130 * MDP(16) + t154;
t170 = (t130 ^ 2 - t132 ^ 2) * MDP(12);
t134 = -pkin(1) - pkin(2);
t161 = qJ(2) * qJD(1);
t120 = t134 * qJD(1) + qJD(2);
t166 = t131 * t120;
t106 = t133 * t161 + t166;
t169 = t128 * t106;
t168 = t129 * t106;
t167 = t131 * MDP(8);
t114 = t128 * t131 - t129 * t133;
t165 = t149 * t114;
t142 = -t131 * qJ(2) + t133 * t134;
t117 = -pkin(3) + t142;
t119 = t133 * qJ(2) + t131 * t134;
t164 = t128 * t117 + t129 * t119;
t153 = t133 * qJD(2);
t160 = qJD(3) * t133;
t163 = qJD(1) * t153 + t120 * t160;
t159 = qJD(5) * t114;
t158 = t149 * qJD(5) * t170;
t156 = t131 * qJD(2);
t155 = t132 * MDP(11);
t135 = qJD(5) ^ 2;
t152 = t135 * MDP(13);
t151 = t135 * MDP(14);
t150 = MDP(17) * qJD(5);
t148 = qJD(3) * t166;
t147 = t149 * t155;
t146 = t131 * t161;
t137 = -t148 + (-qJ(2) * t160 - t156) * qJD(1);
t99 = -qJD(3) * t146 + t163;
t87 = t128 * t137 + t129 * t99;
t105 = t133 * t120 - t146;
t101 = -pkin(3) * t149 + t105;
t90 = t129 * t101 - t169;
t88 = pkin(4) * t149 - t90;
t145 = t149 * t88 - t87;
t141 = t129 * t117 - t128 * t119;
t103 = t142 * qJD(3) + t153;
t104 = -t119 * qJD(3) - t156;
t93 = t129 * t103 + t128 * t104;
t144 = -t149 * (pkin(4) - t141) - t88 - t93;
t95 = t129 * t105 - t169;
t143 = -(-t129 * pkin(3) - pkin(4)) * t149 + t88 + t95;
t86 = t128 * t99 - t129 * t137;
t92 = t128 * t103 - t129 * t104;
t139 = -t149 * t92 + t135 * (-pkin(7) + t164) - t86;
t94 = t128 * t105 + t168;
t138 = (t128 * pkin(3) + pkin(7)) * t135 + t149 * t94 + t86;
t124 = t149 ^ 2;
t91 = t128 * t101 + t168;
t1 = [(-t104 * t149 + t148) * MDP(8) + (t103 * t149 + t163) * MDP(9) + (-t86 * t141 + t87 * t164 - t90 * t92 + t91 * t93) * MDP(10) - 0.2e1 * t158 + (((2 * MDP(5)) + t167) * qJD(2) + (0.2e1 * qJD(2) * MDP(6) + (MDP(8) * t133 - MDP(9) * t131) * qJD(3)) * qJ(2)) * qJD(1) + (-t139 * MDP(16) + t144 * t150 - t152) * t132 + (t151 + t139 * MDP(17) + (t144 * MDP(16) + 0.2e1 * t147) * qJD(5)) * t130; (t86 * t114 + t165 * t91 + t172 * t90) * MDP(10) + (-qJ(2) * MDP(6) - MDP(5)) * qJD(1) ^ 2 + (t87 * MDP(10) + (-MDP(16) * t132 + MDP(17) * t130) * t135) * t140 - t171 * t165 * qJD(5) - ((t130 * t159 + t172 * t132) * MDP(16) + (-t172 * t130 + t132 * t159) * MDP(17) + (t133 * MDP(9) + t167) * t149) * t149; (-t106 * t149 + t137) * MDP(8) + (-t105 * t149 - t99) * MDP(9) + (t90 * t94 - t91 * t95 + (t128 * t87 - t129 * t86) * pkin(3)) * MDP(10) + 0.2e1 * t158 + (-t138 * MDP(16) + t143 * t150 + t152) * t132 + (-t151 + t138 * MDP(17) + (t143 * MDP(16) - 0.2e1 * t147) * qJD(5)) * t130; -t171 * t135; t145 * t154 + t124 * t170 + (t145 * MDP(16) - t124 * t155) * t130;];
tauc = t1;
