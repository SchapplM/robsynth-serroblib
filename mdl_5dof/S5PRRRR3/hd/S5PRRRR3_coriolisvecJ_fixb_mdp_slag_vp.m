% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:21
% DurationCPUTime: 0.43s
% Computational Cost: add. (371->93), mult. (749->139), div. (0->0), fcn. (366->6), ass. (0->60)
t135 = -qJD(3) - qJD(4);
t108 = sin(qJ(4));
t112 = cos(qJ(3));
t109 = sin(qJ(3));
t111 = cos(qJ(4));
t148 = t109 * t111;
t119 = t108 * t112 + t148;
t143 = qJD(4) * t108;
t152 = pkin(2) * qJD(2);
t158 = -pkin(3) * t143 + t119 * t152;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t138 = t110 * MDP(17);
t157 = t107 * MDP(16) + t138;
t156 = (t107 ^ 2 - t110 ^ 2) * MDP(12);
t103 = qJD(2) - t135;
t155 = pkin(4) * t103;
t140 = qJD(5) * t110;
t142 = qJD(4) * t111;
t128 = t109 * t142;
t114 = (t119 * qJD(3) + t128) * pkin(2);
t104 = qJD(2) + qJD(3);
t131 = t112 * t152;
t95 = pkin(3) * t104 + t131;
t129 = t95 * t143;
t80 = qJD(2) * t114 + t129;
t132 = t109 * t152;
t126 = t108 * t132;
t88 = t111 * t95 - t126;
t85 = -t88 - t155;
t154 = t80 * t107 + t85 * t140;
t153 = t135 * t126;
t101 = pkin(2) * t112 + pkin(3);
t151 = t103 * (t101 * t143 + t114);
t150 = t103 * (t108 * t95 + t111 * t132);
t149 = t108 * t109;
t113 = qJD(5) ^ 2;
t147 = t110 * t113;
t146 = t111 * MDP(9);
t144 = MDP(11) * t110;
t141 = qJD(5) * t107;
t137 = t111 * MDP(10);
t136 = t113 * MDP(14);
t134 = MDP(13) * t147 + 0.2e1 * (t144 * t107 - t156) * qJD(5) * t103;
t127 = -t103 * t85 - (qJD(3) * t131 + qJD(4) * t95) * t111 - t153;
t125 = MDP(9) * t128;
t123 = pkin(8) * t113 - t150;
t122 = qJD(5) * (t88 - t155);
t121 = t113 * (pkin(2) * t148 + t101 * t108 + pkin(8)) + t151;
t118 = t111 * t112 - t149;
t81 = t101 * t142 + (t118 * qJD(3) - t109 * t143) * pkin(2);
t120 = qJD(5) * (t103 * (pkin(2) * t149 - t101 * t111 - pkin(4)) - t81);
t117 = -t95 * t142 - t153;
t116 = -t108 * MDP(9) - MDP(7) - t137;
t83 = t85 * t141;
t115 = t83 * MDP(16) + t154 * MDP(17) + t134;
t102 = t103 ^ 2;
t100 = -pkin(3) * t111 - pkin(4);
t99 = pkin(3) * t108 + pkin(8);
t1 = [-t157 * t113; (-t129 - t151) * MDP(9) + (-t103 * t81 + t117) * MDP(10) + ((-t121 - t80) * MDP(16) + MDP(17) * t120) * t110 + (MDP(16) * t120 + t121 * MDP(17) - t136) * t107 + (-qJD(2) * t125 + ((-t109 * MDP(6) - t112 * MDP(7)) * t104 + ((-MDP(6) - t146) * t109 + t116 * t112) * qJD(2)) * qJD(3)) * pkin(2) + t115; -MDP(9) * t129 + t117 * MDP(10) - t107 * t136 + (-t110 * t80 - t99 * t147 + t83) * MDP(16) + (t107 * t113 * t99 + t154) * MDP(17) + (t158 * MDP(9) + (t100 * t141 + t158 * t110) * MDP(16) + (t100 * t140 - t158 * t107) * MDP(17)) * t103 + ((t104 * MDP(7) + t116 * qJD(3)) * t112 + ((-qJD(3) + t104) * MDP(6) + t135 * t146) * t109) * t152 + t134 + (t103 * MDP(10) + t157 * qJD(5)) * (-pkin(3) * t142 + t118 * t152); (-t129 + t150) * MDP(9) + (t103 * t88 + t117) * MDP(10) + (-t125 + (-t119 * MDP(9) - t112 * t137) * qJD(3)) * t152 + ((-t123 - t80) * MDP(16) + MDP(17) * t122) * t110 + (MDP(16) * t122 + t123 * MDP(17) - t136) * t107 + t115; t127 * t138 + t102 * t156 + (t127 * MDP(16) - t102 * t144) * t107;];
tauc = t1;
