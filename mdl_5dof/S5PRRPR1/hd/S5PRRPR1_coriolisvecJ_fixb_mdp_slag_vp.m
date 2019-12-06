% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR1
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
%   see S5PRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:52
% EndTime: 2019-12-05 16:15:55
% DurationCPUTime: 0.57s
% Computational Cost: add. (421->94), mult. (793->144), div. (0->0), fcn. (480->6), ass. (0->51)
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t147 = t125 ^ 2 + t126 ^ 2;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t164 = t130 * MDP(7) + (t126 * MDP(8) + MDP(6)) * t128;
t124 = qJD(2) + qJD(3);
t152 = pkin(2) * qJD(3);
t143 = qJD(2) * t152;
t104 = t124 * qJD(4) + t130 * t143;
t163 = t147 * t104;
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t106 = t125 * t129 + t126 * t127;
t92 = t106 * t124;
t159 = t147 * t124;
t100 = t106 * qJD(5);
t153 = pkin(2) * qJD(2);
t134 = -t130 * t153 + qJD(4);
t158 = pkin(2) * t128;
t157 = pkin(2) * t130;
t137 = t128 * t143;
t117 = -t126 * pkin(4) - pkin(3);
t89 = t117 * t124 + t134;
t148 = t126 * t129;
t151 = t125 * t127;
t105 = -t148 + t151;
t99 = t105 * qJD(5);
t156 = t106 * t137 - t89 * t99;
t155 = t89 * t100 + t105 * t137;
t144 = t128 * t152;
t142 = t124 * t151;
t141 = t124 * t148;
t109 = qJD(5) * t141;
t85 = -qJD(5) * t142 + t109;
t86 = t124 * t100;
t90 = -t141 + t142;
t140 = (-t100 * t92 - t105 * t85 - t106 * t86 + t90 * t99) * MDP(13) + (t106 * t85 - t92 * t99) * MDP(12) + (-t99 * MDP(14) - t100 * MDP(15)) * qJD(5);
t138 = t124 * t125 * t158;
t133 = t147 * (qJ(4) * t124 + t128 * t153);
t120 = t126 * pkin(7);
t116 = qJ(4) + t158;
t115 = t130 * t152 + qJD(4);
t113 = t125 * t137;
t112 = qJ(4) * t126 + t120;
t111 = (-pkin(7) - qJ(4)) * t125;
t110 = t117 - t157;
t107 = -t124 * pkin(3) + t134;
t103 = t116 * t126 + t120;
t102 = (-pkin(7) - t116) * t125;
t1 = [(-MDP(17) * t100 + MDP(18) * t99) * qJD(5); (qJD(3) * t138 + t113) * MDP(9) + (t115 * t159 + t163) * MDP(10) + (t133 * t115 + t116 * t163 + (t107 + (-pkin(3) - t157) * qJD(2)) * t144) * MDP(11) + (t90 * t144 + t110 * t86 + ((-t102 * t127 - t103 * t129) * qJD(5) - t106 * t115) * qJD(5) + t155) * MDP(17) + (t92 * t144 + t110 * t85 + ((-t102 * t129 + t103 * t127) * qJD(5) + t105 * t115) * qJD(5) + t156) * MDP(18) + t140 + t164 * (-qJD(2) - t124) * t152; (-qJD(2) * t138 + t113) * MDP(9) + (t134 * t159 + t163) * MDP(10) + (t133 * qJD(4) + qJ(4) * t163 + (-t133 * t130 + (-pkin(3) * qJD(3) - t107) * t128) * t153) * MDP(11) + (t117 * t86 + ((-t111 * t127 - t112 * t129) * qJD(5) - t106 * qJD(4)) * qJD(5) + (t130 * t100 - t128 * t90) * t153 + t155) * MDP(17) + (t117 * t85 + ((-t111 * t129 + t112 * t127) * qJD(5) + t105 * qJD(4)) * qJD(5) + (-t128 * t92 - t130 * t99) * t153 + t156) * MDP(18) + t140 + t164 * (-qJD(3) + t124) * t153; (-t133 * t124 + t137) * MDP(11) + t109 * MDP(18) - t147 * MDP(10) * t124 ^ 2 + (0.2e1 * t92 * MDP(17) + (-t90 - t142) * MDP(18)) * qJD(5); t92 * t90 * MDP(12) + (-t90 ^ 2 + t92 ^ 2) * MDP(13) + (t109 + (t90 - t142) * qJD(5)) * MDP(14) + (-t106 * t104 - t89 * t92) * MDP(17) + (t105 * t104 + t89 * t90) * MDP(18);];
tauc = t1;
