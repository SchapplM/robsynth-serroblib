% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:34
% EndTime: 2019-12-31 17:02:36
% DurationCPUTime: 0.50s
% Computational Cost: add. (402->94), mult. (765->144), div. (0->0), fcn. (448->6), ass. (0->53)
t115 = sin(pkin(7));
t116 = cos(pkin(7));
t138 = t115 ^ 2 + t116 ^ 2;
t114 = qJD(1) + qJD(2);
t120 = cos(qJ(2));
t143 = pkin(1) * qJD(2);
t135 = qJD(1) * t143;
t96 = qJD(3) * t114 + t120 * t135;
t153 = t138 * t96;
t118 = sin(qJ(2));
t152 = t120 * MDP(6) + t118 * (MDP(7) * t116 + MDP(5));
t117 = sin(qJ(4));
t119 = cos(qJ(4));
t98 = t115 * t119 + t116 * t117;
t85 = t98 * t114;
t93 = t98 * qJD(4);
t149 = pkin(1) * t118;
t148 = pkin(1) * t120;
t126 = t118 * t135;
t109 = -pkin(3) * t116 - pkin(2);
t144 = pkin(1) * qJD(1);
t123 = -t120 * t144 + qJD(3);
t82 = t109 * t114 + t123;
t139 = t116 * t119;
t142 = t115 * t117;
t97 = -t139 + t142;
t92 = t97 * qJD(4);
t147 = t98 * t126 - t82 * t92;
t146 = t97 * t126 + t82 * t93;
t136 = t118 * t143;
t134 = t114 * t142;
t133 = t114 * t139;
t101 = qJD(4) * t133;
t80 = -qJD(4) * t134 + t101;
t81 = t114 * t93;
t83 = -t133 + t134;
t131 = (-t80 * t97 - t81 * t98 + t83 * t92 - t85 * t93) * MDP(12) + (t80 * t98 - t85 * t92) * MDP(11) + (-MDP(13) * t92 - MDP(14) * t93) * qJD(4);
t130 = t138 * t120;
t107 = t120 * t143 + qJD(3);
t129 = t138 * t107;
t128 = t138 * qJD(3);
t127 = t114 * t115 * t149;
t110 = t116 * pkin(6);
t108 = qJ(3) + t149;
t105 = t115 * t126;
t104 = qJ(3) * t116 + t110;
t103 = (-pkin(6) - qJ(3)) * t115;
t102 = t109 - t148;
t100 = qJ(3) * t114 + t118 * t144;
t99 = -pkin(2) * t114 + t123;
t95 = t108 * t116 + t110;
t94 = (-pkin(6) - t108) * t115;
t1 = [(qJD(2) * t127 + t105) * MDP(8) + (t114 * t129 + t153) * MDP(9) + (t108 * t153 + t100 * t129 + (t99 + (-pkin(2) - t148) * qJD(1)) * t136) * MDP(10) + (t83 * t136 + t102 * t81 + ((-t117 * t94 - t119 * t95) * qJD(4) - t98 * t107) * qJD(4) + t146) * MDP(16) + (t85 * t136 + t102 * t80 + ((t117 * t95 - t119 * t94) * qJD(4) + t97 * t107) * qJD(4) + t147) * MDP(17) + t131 + t152 * (-qJD(1) - t114) * t143; (-qJD(1) * t127 + t105) * MDP(8) + ((-t130 * t144 + t128) * t114 + t153) * MDP(9) + (t100 * t128 + qJ(3) * t153 + ((-pkin(2) * qJD(2) - t99) * t118 - t100 * t130) * t144) * MDP(10) + (t109 * t81 + ((-t103 * t117 - t104 * t119) * qJD(4) - t98 * qJD(3)) * qJD(4) + (-t118 * t83 + t120 * t93) * t144 + t146) * MDP(16) + (t109 * t80 + ((-t103 * t119 + t104 * t117) * qJD(4) + t97 * qJD(3)) * qJD(4) + (-t118 * t85 - t120 * t92) * t144 + t147) * MDP(17) + t131 + t152 * (-qJD(2) + t114) * t144; (-t100 * t114 * t138 + t126) * MDP(10) + t101 * MDP(17) - t138 * MDP(9) * t114 ^ 2 + (0.2e1 * t85 * MDP(16) + (-t83 - t134) * MDP(17)) * qJD(4); t85 * t83 * MDP(11) + (-t83 ^ 2 + t85 ^ 2) * MDP(12) + (t101 + (t83 - t134) * qJD(4)) * MDP(13) + (-t82 * t85 - t96 * t98) * MDP(16) + (t82 * t83 + t96 * t97) * MDP(17);];
tauc = t1;
