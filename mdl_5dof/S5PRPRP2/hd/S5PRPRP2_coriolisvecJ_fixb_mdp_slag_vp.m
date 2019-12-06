% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:07
% EndTime: 2019-12-05 15:31:10
% DurationCPUTime: 0.69s
% Computational Cost: add. (413->110), mult. (1020->182), div. (0->0), fcn. (599->4), ass. (0->64)
t113 = cos(pkin(8));
t142 = qJD(2) * t113;
t157 = qJD(4) - t142;
t114 = sin(qJ(4));
t115 = cos(qJ(4));
t112 = sin(pkin(8));
t100 = -pkin(3) * t113 - pkin(6) * t112 - pkin(2);
t125 = qJ(5) * t112 - t100;
t156 = -qJ(3) * t113 * t115 + t125 * t114;
t108 = t112 ^ 2;
t110 = t114 ^ 2;
t111 = t115 ^ 2;
t155 = (-t115 * t114 * MDP(9) + MDP(10) * (t110 - t111)) * t108;
t109 = t113 ^ 2;
t143 = qJD(2) * t112;
t130 = t115 * t143;
t95 = t100 * qJD(2) + qJD(3);
t93 = t115 * t95;
t98 = qJ(3) * t142 + qJD(1) * t112;
t87 = -qJ(5) * t130 - t114 * t98 + t93;
t84 = pkin(4) * t157 + t87;
t154 = t84 - t87;
t140 = qJD(3) * t115;
t104 = t113 * t140;
t138 = qJD(4) * t115;
t153 = -qJD(2) * t104 - t95 * t138;
t152 = t100 * t138 + t104;
t151 = t157 * t113;
t116 = qJD(2) ^ 2;
t150 = t108 * t116;
t149 = t112 * t115;
t148 = t113 * t114;
t134 = qJD(2) * qJD(4);
t127 = t112 * t134;
t94 = t115 * pkin(4) * t127 + qJD(3) * t143;
t147 = t108 + t109;
t144 = qJD(2) * t108;
t141 = qJD(3) * t114;
t139 = qJD(4) * t114;
t137 = qJD(5) * t114;
t136 = qJD(5) * t115;
t135 = qJD(4) - t157;
t132 = qJ(3) * t139;
t131 = t114 * t143;
t129 = t113 * t141;
t128 = t112 * t139;
t126 = pkin(4) * t114 + qJ(3);
t124 = t147 * qJD(2);
t107 = t113 * qJD(1);
t97 = qJ(3) * t143 - t107;
t121 = t112 * t97 + t113 * t98;
t120 = -t114 * t95 - t115 * t98;
t119 = -t98 * t139 - t153;
t118 = t120 * qJD(4);
t117 = -t157 ^ 2 - t150;
t99 = t127 * t148;
t91 = t126 * t143 + qJD(5) - t107;
t89 = -t125 * t115 + (-qJ(3) * t114 - pkin(4)) * t113;
t88 = -qJ(5) * t131 - t120;
t86 = t156 * qJD(4) - t112 * t136 - t129;
t85 = -t112 * t137 + (-qJ(3) * t148 - qJ(5) * t149) * qJD(4) + t152;
t83 = t118 + (-t129 + (qJ(5) * t139 - t136) * t112) * qJD(2);
t82 = (-qJ(5) * t138 - t137) * t143 + t119;
t1 = [-t94 * t113 * MDP(17) + t99 * MDP(15) + ((-t114 * t83 + t115 * t82) * MDP(17) + ((MDP(15) * t157 - MDP(17) * t88) * t114 + ((-t157 - t142) * MDP(14) - t84 * MDP(17)) * t115) * qJD(4)) * t112; (-t128 * t157 + t99) * MDP(11) + ((-t151 + (0.2e1 * t108 + t109) * qJD(2)) * t141 + ((-t100 * t157 + t113 * t95) * t114 + ((t144 - t151) * qJ(3) + t121) * t115) * qJD(4)) * MDP(14) + (-(-t113 * t132 + t152) * t157 + t119 * t113 - t97 * t128 + (-t132 + 0.2e1 * t140) * t144) * MDP(15) + (-t156 * t82 + t83 * t89 + t84 * t86 + t88 * t85) * MDP(17) + 0.2e1 * t155 * t134 + (0.2e1 * MDP(7) * t124 + (qJ(3) * t124 + t121) * MDP(8)) * qJD(3) + ((-t157 + t142) * MDP(12) * t138 + (-t114 * t82 - t115 * t83 + (t114 * t84 - t115 * t88) * qJD(4) + (-t114 * t85 - t115 * t86 + (t114 * t89 + t115 * t156) * qJD(4)) * qJD(2)) * MDP(16) + (t94 * t126 + t91 * (pkin(4) * t138 + qJD(3))) * MDP(17)) * t112; -t147 * MDP(7) * t116 + (-t91 * t112 * MDP(17) - t121 * MDP(8)) * qJD(2) + (t117 * MDP(15) + (t157 * t88 + t83) * MDP(17)) * t115 + (t117 * MDP(14) + (-t157 * t84 + t82) * MDP(17)) * t114; (-t120 * t157 + t118 + (-t97 * t149 - t129) * qJD(2)) * MDP(14) + (t93 * t157 + (t135 * t98 + t97 * t143) * t114 + t153) * MDP(15) + (pkin(4) * qJD(4) - t154) * MDP(16) * t131 + (t154 * t88 + (-t91 * t130 + t83) * pkin(4)) * MDP(17) - (t114 * MDP(11) + t115 * MDP(12)) * t135 * t143 - t155 * t116; ((t114 * t88 + t115 * t84) * t143 + t94) * MDP(17) + (-t110 - t111) * MDP(16) * t150;];
tauc = t1;
