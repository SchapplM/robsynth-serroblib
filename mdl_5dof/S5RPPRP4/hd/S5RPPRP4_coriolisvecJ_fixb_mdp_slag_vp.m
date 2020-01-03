% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:21
% EndTime: 2019-12-31 17:52:23
% DurationCPUTime: 0.38s
% Computational Cost: add. (415->96), mult. (797->145), div. (0->0), fcn. (374->4), ass. (0->60)
t157 = -qJ(2) * MDP(6) - MDP(5);
t113 = sin(qJ(4));
t114 = cos(qJ(4));
t110 = sin(pkin(7));
t111 = cos(pkin(7));
t140 = qJD(1) * qJ(2);
t115 = -pkin(1) - pkin(2);
t98 = t115 * qJD(1) + qJD(2);
t93 = t110 * t98 + t111 * t140;
t89 = -qJD(1) * pkin(6) + t93;
t156 = t114 * qJD(3) - t113 * t89;
t155 = t113 * MDP(15) + t114 * MDP(16);
t154 = MDP(10) * t114;
t138 = qJD(2) * t111;
t125 = -qJD(5) + t138;
t121 = t125 * t113;
t123 = -t113 * qJD(3) - t114 * t89;
t81 = t123 * qJD(4) + (qJD(4) * t114 * qJ(5) - t121) * qJD(1);
t139 = qJD(1) * t114;
t86 = -qJ(5) * t139 - t123;
t153 = (qJD(4) * t86 + t81) * MDP(18);
t108 = t113 ^ 2;
t109 = t114 ^ 2;
t152 = (t108 + t109) * MDP(17);
t151 = (t108 - t109) * MDP(11);
t149 = qJD(4) * pkin(4);
t136 = t113 * qJD(1);
t85 = qJ(5) * t136 + t156;
t82 = t85 + t149;
t150 = t82 - t85;
t147 = t114 * t86;
t144 = t111 * qJ(2) + t110 * t115;
t96 = -pkin(6) + t144;
t146 = qJ(5) - t96;
t116 = qJD(4) ^ 2;
t117 = qJD(1) ^ 2;
t141 = t116 + t117;
t135 = t113 * qJD(4);
t133 = qJD(1) * qJD(2);
t132 = qJD(1) * qJD(4);
t130 = t113 * t132;
t99 = t110 * t133;
t92 = -t110 * t140 + t111 * t98;
t128 = qJD(4) * t146;
t127 = -t110 * qJ(2) + t111 * t115;
t126 = 0.2e1 * t99;
t95 = pkin(3) - t127;
t88 = qJD(1) * pkin(3) - t92;
t120 = t125 * t114;
t80 = t156 * qJD(4) + (qJ(5) * t135 + t120) * qJD(1);
t124 = (-qJD(4) * t82 + t80) * MDP(18);
t119 = -t116 * t96 + t126;
t118 = qJD(4) * (-qJD(1) * t95 - t138 - t88);
t94 = -pkin(4) * t130 + t99;
t91 = t146 * t114;
t90 = t146 * t113;
t87 = pkin(4) * t139 + qJD(5) + t88;
t84 = t114 * t128 - t121;
t83 = t113 * t128 + t120;
t1 = [MDP(7) * t126 + (-t92 * t110 + t93 * t111 + (-t110 * t127 + t111 * t144) * qJD(1)) * qJD(2) * MDP(9) + 0.2e1 * t130 * t154 - 0.2e1 * t132 * t151 + (t113 * t118 + t119 * t114) * MDP(15) + (-t119 * t113 + t114 * t118) * MDP(16) + (t81 * t113 - t80 * t114 + (t113 * t86 + t114 * t82) * qJD(4) + (t113 * t84 - t114 * t83 + (-t113 * t91 + t114 * t90) * qJD(4)) * qJD(1)) * MDP(17) + (-t80 * t91 + t86 * t83 + t81 * t90 + t82 * t84 + t94 * (t114 * pkin(4) + t95) + t87 * (-pkin(4) * t135 + t110 * qJD(2))) * MDP(18) + (-t114 * MDP(12) + t113 * MDP(13)) * t116 + 0.2e1 * (t111 * MDP(8) - t157) * t133; t157 * t117 + (-t94 * MDP(18) + (-MDP(8) + t152) * t117 + (-t93 * MDP(9) + (t113 * t82 - t147) * MDP(18) + 0.2e1 * t155 * qJD(4)) * qJD(1)) * t111 + (-t117 * MDP(7) + (-t87 * MDP(18) + t92 * MDP(9)) * qJD(1) + (-t141 * MDP(15) + t124) * t114 + (t141 * MDP(16) - t153) * t113) * t110; (-t116 * MDP(16) + t153) * t114 + (-t116 * MDP(15) + t124) * t113; (t149 - t150) * MDP(17) * t139 + (t150 * t86 + (t87 * t136 + t81) * pkin(4)) * MDP(18) + t155 * (t88 - t138) * qJD(1) + (-t113 * t154 + t151) * t117; -t117 * t152 + (t99 + (t147 + (-t82 - t149) * t113) * qJD(1)) * MDP(18);];
tauc = t1;
