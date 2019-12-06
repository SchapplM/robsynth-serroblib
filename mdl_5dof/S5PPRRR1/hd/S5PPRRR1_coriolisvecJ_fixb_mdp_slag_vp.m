% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:52
% EndTime: 2019-12-05 15:12:54
% DurationCPUTime: 0.44s
% Computational Cost: add. (355->82), mult. (879->130), div. (0->0), fcn. (680->8), ass. (0->56)
t114 = sin(pkin(9));
t115 = cos(pkin(9));
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t157 = -t118 * t114 + t121 * t115;
t117 = sin(qJ(4));
t141 = qJD(4) * t117;
t104 = t121 * t114 + t118 * t115;
t100 = t104 * qJD(1);
t120 = cos(qJ(4));
t144 = t120 * t100;
t99 = t157 * qJD(1);
t156 = -pkin(3) * t141 + t117 * t99 + t144;
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t137 = t119 * MDP(15);
t155 = t116 * MDP(14) + t137;
t154 = (t116 ^ 2 - t119 ^ 2) * MDP(10);
t101 = t157 * qJD(3);
t97 = qJD(1) * t101;
t102 = t104 * qJD(3);
t98 = qJD(1) * t102;
t130 = t117 * t97 + t120 * t98;
t96 = qJD(3) * pkin(3) + t99;
t88 = t117 * t96 + t144;
t80 = t88 * qJD(4) + t130;
t129 = -t100 * t141 - t117 * t98;
t153 = -(qJD(4) * t96 + t97) * t120 - t129;
t111 = qJD(3) + qJD(4);
t152 = t111 * pkin(4);
t136 = t119 * qJD(5);
t147 = t117 * t100;
t87 = t120 * t96 - t147;
t85 = -t87 - t152;
t151 = t80 * t116 + t85 * t136;
t150 = t88 * t111;
t149 = MDP(9) * t119;
t122 = qJD(5) ^ 2;
t148 = t116 * t122;
t145 = t119 * t122;
t140 = qJD(4) * t120;
t138 = t116 * qJD(5);
t135 = t122 * MDP(12);
t132 = MDP(11) * t145 + 0.2e1 * (t149 * t116 - t154) * qJD(5) * t111;
t131 = -t85 * t111 + t153;
t126 = pkin(7) * t122 - t150;
t125 = qJD(5) * (t87 - t152);
t124 = -t117 * t104 + t120 * t157;
t92 = t120 * t104 + t117 * t157;
t110 = t111 ^ 2;
t109 = -t120 * pkin(3) - pkin(4);
t108 = t117 * pkin(3) + pkin(7);
t83 = t85 * t138;
t82 = t92 * qJD(4) + t117 * t101 + t120 * t102;
t81 = t124 * qJD(4) + t120 * t101 - t117 * t102;
t1 = [(-t81 * t138 - t92 * t145) * MDP(14) + (-t81 * t136 + t92 * t148) * MDP(15) + (-t102 * MDP(4) - t101 * MDP(5)) * qJD(3) + (-t82 * MDP(7) - t81 * MDP(8) + (-t119 * t82 - t124 * t138) * MDP(14) + (t116 * t82 - t124 * t136) * MDP(15)) * t111; -t155 * t122; (-t100 * t140 - t96 * t141 - t130) * MDP(7) + (-t120 * t97 - t96 * t140 - t129) * MDP(8) - t116 * t135 + (-t108 * t145 - t80 * t119 + t83) * MDP(14) + (t108 * t148 + t151) * MDP(15) + (t156 * MDP(7) + (t109 * t138 + t156 * t119) * MDP(14) + (t109 * t136 - t156 * t116) * MDP(15)) * t111 + (t100 * MDP(4) + t99 * MDP(5) + (-t104 * MDP(4) - MDP(5) * t157) * qJD(1)) * qJD(3) + t132 + (t111 * MDP(8) + t155 * qJD(5)) * (-pkin(3) * t140 + t120 * t99 - t147); (t150 - t80) * MDP(7) + (t87 * t111 + t153) * MDP(8) + t83 * MDP(14) + t151 * MDP(15) + ((-t126 - t80) * MDP(14) + MDP(15) * t125) * t119 + (MDP(14) * t125 + t126 * MDP(15) - t135) * t116 + t132; t131 * t137 + t110 * t154 + (t131 * MDP(14) - t110 * t149) * t116;];
tauc = t1;
