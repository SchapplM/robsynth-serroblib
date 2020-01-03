% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:11
% EndTime: 2019-12-31 17:13:12
% DurationCPUTime: 0.48s
% Computational Cost: add. (368->104), mult. (676->163), div. (0->0), fcn. (301->4), ass. (0->66)
t155 = qJ(4) + pkin(6);
t114 = sin(qJ(3));
t112 = t114 ^ 2;
t116 = cos(qJ(3));
t113 = t116 ^ 2;
t154 = (t112 - t113) * MDP(8);
t117 = cos(qJ(2));
t153 = t117 * pkin(1);
t148 = qJD(3) * pkin(3);
t111 = qJD(1) + qJD(2);
t115 = sin(qJ(2));
t150 = pkin(1) * qJD(1);
t134 = t115 * t150;
t126 = t155 * t111 + t134;
t86 = t126 * t114;
t85 = -t86 + t148;
t151 = t85 + t86;
t149 = pkin(1) * qJD(2);
t132 = qJD(1) * t149;
t104 = t115 * t132;
t139 = qJD(3) * t116;
t133 = t117 * t150;
t99 = -t111 * pkin(2) - t133;
t147 = t114 * t104 + t99 * t139;
t146 = MDP(7) * t114;
t118 = qJD(3) ^ 2;
t145 = t114 * t118;
t144 = t116 * t118;
t106 = t115 * pkin(1) + pkin(6);
t143 = -qJ(4) - t106;
t140 = qJD(3) * t114;
t130 = t111 * t140;
t92 = pkin(3) * t130 + t104;
t142 = -t112 - t113;
t138 = qJD(3) * t117;
t137 = -qJD(2) + t111;
t136 = t117 * t149;
t135 = pkin(3) * t140;
t131 = -t116 * pkin(3) - pkin(2);
t129 = t111 * t139;
t123 = t117 * t132;
t119 = qJD(4) * t111 + t123;
t121 = qJD(3) * t126;
t81 = -t114 * t121 + t119 * t116;
t82 = -t119 * t114 - t116 * t121;
t128 = -t82 * t114 + t81 * t116;
t127 = qJD(3) * t155;
t125 = qJD(3) * t143;
t124 = qJD(3) * t136;
t87 = t126 * t116;
t122 = t114 * t85 - t116 * t87;
t120 = -0.2e1 * t111 * qJD(3) * t154 - MDP(10) * t145 + MDP(9) * t144 + 0.2e1 * t129 * t146;
t110 = t111 ^ 2;
t109 = t116 * qJ(4);
t108 = t116 * qJD(4);
t102 = t116 * pkin(6) + t109;
t101 = t155 * t114;
t96 = t116 * t106 + t109;
t95 = t143 * t114;
t93 = t99 * t140;
t91 = -t114 * qJD(4) - t116 * t127;
t90 = -t114 * t127 + t108;
t89 = t131 * t111 + qJD(4) - t133;
t84 = (-qJD(4) - t136) * t114 + t116 * t125;
t83 = t114 * t125 + t116 * t136 + t108;
t1 = [-t104 * MDP(5) - MDP(6) * t123 + (-t116 * t104 - t106 * t144 - t114 * t124 + t93) * MDP(12) + (t106 * t145 - t116 * t124 + t147) * MDP(13) + (-t85 * t139 - t87 * t140 + t128) * MDP(14) + (t81 * t96 + t87 * t83 + t82 * t95 + t85 * t84 + t92 * (t131 - t153) + t89 * (t115 * t149 + t135)) * MDP(15) + ((-t114 * t84 + t116 * t83) * MDP(14) + ((-t114 * t96 - t116 * t95) * MDP(14) + (t114 * MDP(12) + t116 * MDP(13)) * (-pkin(2) - t153)) * qJD(3) + (-t117 * MDP(6) + (-t116 * MDP(12) + t114 * MDP(13) - MDP(5)) * t115) * t149) * t111 + t120; (t111 * t134 - t104) * MDP(5) + t137 * MDP(6) * t133 + (-pkin(2) * t130 - pkin(6) * t144 + t93 + (t137 * t116 * t115 + t114 * t138) * t150) * MDP(12) + (-pkin(2) * t129 + pkin(6) * t145 + (-t111 * t114 * t115 + t116 * t138) * t150 + t147) * MDP(13) + ((-t114 * t87 - t116 * t85) * qJD(3) + (-t114 * t91 + t116 * t90 + (t101 * t116 - t102 * t114) * qJD(3) + t142 * t133) * t111 + t128) * MDP(14) + (t81 * t102 + t87 * t90 - t82 * t101 + t85 * t91 + t92 * t131 + t89 * t135 + (-t115 * t89 + t122 * t117) * t150) * MDP(15) + t120; (t82 * pkin(3) + t151 * t87) * MDP(15) + ((-t111 * t99 - t123) * MDP(12) - pkin(3) * t111 * t89 * MDP(15)) * t114 + t110 * t154 + (-MDP(13) * t123 - t110 * t146 + (-t99 * MDP(13) + (-t148 + t151) * MDP(14)) * t111) * t116; (t122 * t111 + t92) * MDP(15) + t142 * MDP(14) * t110;];
tauc = t1;
