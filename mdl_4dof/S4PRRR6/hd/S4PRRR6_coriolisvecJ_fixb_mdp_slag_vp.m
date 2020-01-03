% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:08
% EndTime: 2019-12-31 16:35:10
% DurationCPUTime: 0.66s
% Computational Cost: add. (312->104), mult. (821->177), div. (0->0), fcn. (538->6), ass. (0->63)
t121 = sin(qJ(4));
t122 = sin(qJ(3));
t124 = cos(qJ(4));
t125 = cos(qJ(3));
t105 = t121 * t125 + t122 * t124;
t118 = qJD(3) + qJD(4);
t168 = t105 * t118;
t89 = qJD(2) * t168;
t163 = qJD(4) - t118;
t147 = qJD(2) * qJD(3);
t169 = -0.2e1 * t147;
t104 = t121 * t122 - t124 * t125;
t131 = t104 * t118;
t123 = sin(qJ(2));
t127 = qJD(3) ^ 2;
t128 = qJD(2) ^ 2;
t167 = (t127 + t128) * t123;
t166 = (t122 ^ 2 - t125 ^ 2) * MDP(6);
t151 = qJD(3) * t125;
t165 = -qJD(4) * t125 - t151;
t152 = qJD(3) * t122;
t164 = qJD(4) * t122 + t152;
t162 = pkin(5) + pkin(6);
t161 = qJD(2) * pkin(2);
t148 = t123 * qJD(1);
t112 = qJD(2) * pkin(5) + t148;
t139 = pkin(6) * qJD(2) + t112;
t98 = t139 * t125;
t160 = t124 * t98;
t159 = t122 * t127;
t158 = t125 * t127;
t126 = cos(qJ(2));
t155 = qJD(1) * t126;
t154 = qJD(2) * t122;
t153 = qJD(2) * t125;
t146 = pkin(3) * t152;
t145 = pkin(3) * t154;
t117 = -pkin(3) * t125 - pkin(2);
t144 = qJD(3) * t162;
t143 = t121 * t154;
t142 = t124 * t153;
t97 = t139 * t122;
t96 = qJD(3) * pkin(3) - t97;
t141 = -pkin(3) * t118 - t96;
t140 = t125 * t147;
t137 = t126 * t169;
t135 = qJD(2) * t161;
t101 = -t121 * t153 - t124 * t154;
t102 = t117 * qJD(2) - t155;
t93 = -t112 * t152 + (-pkin(6) * t152 + t125 * t155) * qJD(2);
t94 = -t112 * t151 + (-pkin(6) * t151 - t122 * t155) * qJD(2);
t134 = t102 * t101 - t121 * t93 + t124 * t94;
t88 = qJD(4) * t142 - t118 * t143 + t124 * t140;
t99 = -t142 + t143;
t133 = -t101 * t99 * MDP(12) + (t118 * t99 + t88) * MDP(14) + (-t101 * t118 - t89) * MDP(15) + (t101 ^ 2 - t99 ^ 2) * MDP(13);
t130 = -0.2e1 * qJD(3) * t161;
t129 = t102 * t99 + (t163 * t98 - t94) * t121;
t109 = t162 * t125;
t108 = t162 * t122;
t107 = t125 * t144;
t106 = t122 * t144;
t103 = (t146 + t148) * qJD(2);
t1 = [(t122 * t137 - t125 * t167) * MDP(10) + (t122 * t167 + t125 * t137) * MDP(11) + (-0.2e1 * t126 * t89 + ((t164 * t121 + t165 * t124) * t118 + qJD(2) * t99) * t123) * MDP(17) + ((qJD(2) * t131 - t88) * t126 + (-(t165 * t121 - t164 * t124) * t118 - qJD(2) * t101) * t123) * MDP(18) + (-t123 * MDP(3) - t126 * MDP(4)) * t128; 0.2e1 * t122 * MDP(5) * t140 + t166 * t169 + MDP(7) * t158 - MDP(8) * t159 + (-pkin(5) * t158 + t122 * t130) * MDP(10) + (pkin(5) * t159 + t125 * t130) * MDP(11) + (t101 * t131 + t105 * t88) * MDP(12) + (t101 * t168 - t104 * t88 - t105 * t89 + t131 * t99) * MDP(13) + (t99 * t146 + t117 * t89 + t103 * t104 + t102 * t168 + (-t123 * t99 + t126 * t168) * qJD(1)) * MDP(17) + (-t101 * t146 + t117 * t88 + t103 * t105 - t102 * t131 + (t123 * t101 - t126 * t131) * qJD(1)) * MDP(18) + (-t131 * MDP(14) - t168 * MDP(15) + (t106 * t121 - t107 * t124 + (t108 * t121 - t109 * t124) * qJD(4)) * MDP(17) + (t106 * t124 + t107 * t121 - (-t108 * t124 - t109 * t121) * qJD(4)) * MDP(18)) * t118; t128 * t166 + t125 * MDP(11) * t135 + (-(t121 * t97 - t160) * t118 - t99 * t145 + (t141 * t121 - t160) * qJD(4) + t134) * MDP(17) + (t101 * t145 + (t141 * qJD(4) - t97 * t118 - t93) * t124 + t129) * MDP(18) + t133 + (-t128 * t125 * MDP(5) + MDP(10) * t135) * t122; (t134 + t163 * (-t121 * t96 - t160)) * MDP(17) + ((-t163 * t96 - t93) * t124 + t129) * MDP(18) + t133;];
tauc = t1;
