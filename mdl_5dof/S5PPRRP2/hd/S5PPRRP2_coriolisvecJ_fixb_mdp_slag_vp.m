% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:12
% EndTime: 2019-12-05 15:09:14
% DurationCPUTime: 0.64s
% Computational Cost: add. (424->106), mult. (1081->150), div. (0->0), fcn. (725->6), ass. (0->57)
t119 = sin(pkin(8));
t120 = cos(pkin(8));
t122 = sin(qJ(3));
t124 = cos(qJ(3));
t162 = -t119 * t122 + t120 * t124;
t104 = t162 * qJD(1);
t109 = t119 * t124 + t120 * t122;
t107 = t109 * qJD(3);
t123 = cos(qJ(4));
t105 = t109 * qJD(1);
t100 = qJD(3) * pkin(6) + t105;
t121 = sin(qJ(4));
t155 = t100 * t121;
t92 = qJD(2) * t123 - t155;
t164 = qJD(5) - t92;
t106 = t162 * qJD(3);
t163 = qJD(1) * t106 + qJD(2) * qJD(4);
t89 = -qJD(4) * pkin(4) + t164;
t93 = qJD(2) * t121 + t100 * t123;
t90 = qJD(4) * qJ(5) + t93;
t117 = t121 ^ 2;
t118 = t123 ^ 2;
t161 = (t117 - t118) * MDP(7);
t137 = MDP(11) + MDP(13);
t136 = MDP(12) - MDP(15);
t160 = -MDP(5) + (t117 + t118) * MDP(14);
t125 = qJD(4) ^ 2;
t159 = pkin(6) * t125;
t157 = qJD(3) * pkin(3);
t156 = t163 * t123;
t154 = t109 * t125;
t149 = t121 * t123;
t102 = qJD(1) * t107;
t129 = pkin(4) * t121 - qJ(5) * t123;
t103 = t129 * qJD(4) - qJD(5) * t121;
t146 = qJD(3) * t103;
t145 = qJD(3) * t105;
t111 = -pkin(4) * t123 - qJ(5) * t121 - pkin(3);
t144 = qJD(3) * t111;
t143 = qJD(3) * t123;
t141 = qJD(4) * t121;
t140 = qJD(4) * t123;
t87 = t100 * t140 + t163 * t121;
t99 = -t104 - t157;
t135 = t99 - t157;
t134 = -t102 - t159;
t133 = t92 + t155;
t132 = pkin(6) * MDP(16) + MDP(14);
t91 = -t104 + t144;
t130 = t91 + t144;
t128 = t121 * t89 + t123 * t90;
t88 = t146 + t102;
t127 = -t146 - t88 - t159;
t126 = qJD(3) ^ 2;
t110 = t129 * qJD(3);
t86 = (qJD(5) - t155) * qJD(4) + t156;
t1 = [t137 * (-t106 * t141 - t123 * t154 + (-t107 * t123 - t141 * t162) * qJD(3)) + t136 * (-t106 * t140 + t121 * t154 + (t107 * t121 - t140 * t162) * qJD(3)) + (-t107 * MDP(4) + t160 * t106) * qJD(3) + (t128 * t106 + t91 * t107 - t88 * t162 + (t87 * t121 + t86 * t123 + (-t121 * t90 + t123 * t89) * qJD(4)) * t109) * MDP(16); (t128 * qJD(4) + t86 * t121 - t87 * t123) * MDP(16) + (-t137 * t121 - t136 * t123) * t125; -t102 * MDP(4) + (t111 * t88 + (t103 - t105) * t91) * MDP(16) + (t105 * MDP(4) - 0.2e1 * qJD(4) * t161 + (-MDP(5) - t160) * t104) * qJD(3) + (t125 * MDP(8) + t134 * MDP(11) + t127 * MDP(13) + t86 * MDP(14) + (pkin(6) * t86 - t104 * t90) * MDP(16) + ((t104 + t135) * MDP(12) + (-t104 - t130) * MDP(15) + t132 * t89) * qJD(4)) * t123 + (-t125 * MDP(9) + (-t134 - t145) * MDP(12) + t87 * MDP(14) + (t127 + t145) * MDP(15) + (pkin(6) * t87 - t104 * t89) * MDP(16) + (t135 * MDP(11) + t130 * MDP(13) + 0.2e1 * MDP(6) * t143 - t132 * t90) * qJD(4)) * t121 + t137 * (t104 * t141 + t105 * t143); (qJ(5) * t86 - t110 * t91 + t164 * t90 - t89 * t93) * MDP(16) + (-MDP(6) * t149 + t161) * t126 + (t133 * MDP(12) + (0.2e1 * qJD(5) - t133) * MDP(15) + t137 * t93) * qJD(4) + ((-t99 * MDP(11) - t91 * MDP(13) + t110 * MDP(15)) * t121 + (-t99 * MDP(12) + t110 * MDP(13) + t91 * MDP(15)) * t123) * qJD(3) + (-pkin(4) * MDP(16) - t137) * t87 - t136 * t156; -t126 * MDP(13) * t149 + (-t117 * t126 - t125) * MDP(15) + (qJD(3) * t121 * t91 - qJD(4) * t90 + t87) * MDP(16);];
tauc = t1;
