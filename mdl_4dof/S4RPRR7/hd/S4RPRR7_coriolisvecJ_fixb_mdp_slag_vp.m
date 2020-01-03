% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:07
% EndTime: 2019-12-31 16:54:09
% DurationCPUTime: 1.05s
% Computational Cost: add. (728->160), mult. (1992->230), div. (0->0), fcn. (1456->6), ass. (0->73)
t140 = sin(pkin(7));
t141 = cos(pkin(7));
t185 = (t140 ^ 2 + t141 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t145 = cos(qJ(3));
t167 = qJD(1) * t145;
t133 = t141 * t167;
t143 = sin(qJ(3));
t168 = qJD(1) * t143;
t161 = t140 * t168;
t121 = t133 - t161;
t120 = qJD(4) - t121;
t126 = t140 * t145 + t141 * t143;
t122 = t126 * qJD(1);
t181 = pkin(5) + qJ(2);
t130 = t181 * t140;
t127 = qJD(1) * t130;
t131 = t181 * t141;
t128 = qJD(1) * t131;
t110 = -t127 * t143 + t128 * t145;
t150 = t126 * qJD(2);
t95 = qJD(1) * t150 + qJD(3) * t110;
t183 = (pkin(3) * t122 + t120 * pkin(6)) * t120 + t95;
t109 = -t127 * t145 - t128 * t143;
t104 = -qJD(3) * pkin(3) - t109;
t125 = t140 * t143 - t145 * t141;
t135 = -pkin(2) * t141 - pkin(1);
t108 = pkin(3) * t125 - pkin(6) * t126 + t135;
t112 = -t130 * t143 + t131 * t145;
t124 = t126 * qJD(3);
t119 = qJD(1) * t124;
t123 = t125 * qJD(3);
t149 = t125 * qJD(2);
t94 = -qJD(1) * t149 + qJD(3) * t109;
t129 = t135 * qJD(1) + qJD(2);
t98 = -pkin(3) * t121 - pkin(6) * t122 + t129;
t155 = -t130 * t145 - t131 * t143;
t99 = t155 * qJD(3) - t149;
t182 = -t104 * t123 - t112 * t119 - (qJD(4) * t108 + t99) * t120 - (qJD(4) * t98 + t94) * t125 + t95 * t126;
t142 = sin(qJ(4));
t165 = qJD(4) * t142;
t132 = qJD(3) * t133;
t118 = -qJD(3) * t161 + t132;
t144 = cos(qJ(4));
t164 = t144 * qJD(3);
t170 = qJD(4) * t164 + t144 * t118;
t96 = -t122 * t165 + t170;
t180 = t142 * t96;
t178 = t108 * t119;
t171 = t122 * t142;
t113 = -t164 + t171;
t177 = t113 * t120;
t176 = t113 * t122;
t115 = qJD(3) * t142 + t122 * t144;
t175 = t115 * t120;
t174 = t115 * t122;
t173 = t118 * t142;
t172 = t119 * t142;
t117 = t144 * t119;
t166 = qJD(4) * t126;
t163 = qJD(1) * qJD(2);
t157 = t144 * t120;
t105 = qJD(3) * pkin(6) + t110;
t93 = t105 * t144 + t142 * t98;
t156 = t105 * t142 - t144 * t98;
t153 = t117 + (t121 * t142 - t165) * t120;
t152 = -t144 * t123 - t126 * t165;
t147 = -pkin(6) * t119 + (t104 + t109) * t120;
t107 = pkin(3) * t124 + pkin(6) * t123;
t102 = pkin(3) * t119 - pkin(6) * t118;
t101 = t144 * t102;
t100 = t112 * qJD(3) + t150;
t97 = qJD(4) * t115 + t173;
t1 = [(t118 * t126 - t122 * t123) * MDP(8) + (-t118 * t125 - t119 * t126 - t121 * t123 - t122 * t124) * MDP(9) + (t119 * t135 + t124 * t129) * MDP(13) + (t118 * t135 - t123 * t129) * MDP(14) + (t126 * t144 * t96 + t152 * t115) * MDP(15) + (-(-t113 * t144 - t115 * t142) * t123 + (-t180 - t144 * t97 + (t113 * t142 - t115 * t144) * qJD(4)) * t126) * MDP(16) + (t115 * t124 + t126 * t117 + t152 * t120 + t125 * t96) * MDP(17) + (-t126 * t172 - t113 * t124 - t125 * t97 + (t142 * t123 - t144 * t166) * t120) * MDP(18) + (t119 * t125 + t120 * t124) * MDP(19) + (t100 * t113 + t101 * t125 - t155 * t97 - t156 * t124 + (t107 * t120 + t178 + (t104 * t126 - t105 * t125 - t112 * t120) * qJD(4)) * t144 + t182 * t142) * MDP(20) + (t100 * t115 - t155 * t96 - t93 * t124 + (-(-qJD(4) * t112 + t107) * t120 - t178 - (-qJD(4) * t105 + t102) * t125 - t104 * t166) * t142 + t182 * t144) * MDP(21) + 0.2e1 * t163 * t185 + (-t123 * MDP(10) - t124 * MDP(11) - t100 * MDP(13) - t99 * MDP(14)) * qJD(3); t132 * MDP(14) + (t153 - t176) * MDP(20) + (-t120 ^ 2 * t144 - t172 - t174) * MDP(21) + ((t140 * t167 + t141 * t168 + t122) * MDP(13) + (t121 - t161) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t185; -t121 ^ 2 * MDP(9) + (t132 + (-t121 - t161) * qJD(3)) * MDP(10) + (-t121 * t129 + t125 * t163) * MDP(14) + (t115 * t157 + t180) * MDP(15) + ((t96 - t177) * t144 + (-t97 - t175) * t142) * MDP(16) + (t120 * t157 + t172 - t174) * MDP(17) + (t153 + t176) * MDP(18) + (-pkin(3) * t97 - t110 * t113 + t147 * t142 - t183 * t144) * MDP(20) + (-pkin(3) * t96 - t110 * t115 + t183 * t142 + t147 * t144) * MDP(21) + (-t120 * MDP(19) + t156 * MDP(20) + t93 * MDP(21) - t121 * MDP(8) + MDP(9) * t122 + (-qJD(2) - t129) * MDP(13)) * t122; t115 * t113 * MDP(15) + (-t113 ^ 2 + t115 ^ 2) * MDP(16) + (t170 + t177) * MDP(17) + (-t173 + t175) * MDP(18) + t119 * MDP(19) + (-t104 * t115 + t120 * t93 - t142 * t94 + t101) * MDP(20) + (-t102 * t142 + t104 * t113 - t120 * t156 - t144 * t94) * MDP(21) + (-MDP(17) * t171 - t115 * MDP(18) - t93 * MDP(20) + t156 * MDP(21)) * qJD(4);];
tauc = t1;
