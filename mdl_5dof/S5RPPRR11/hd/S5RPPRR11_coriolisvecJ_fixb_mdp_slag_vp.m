% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:56
% EndTime: 2019-12-31 18:05:59
% DurationCPUTime: 1.02s
% Computational Cost: add. (464->175), mult. (980->256), div. (0->0), fcn. (518->4), ass. (0->84)
t118 = pkin(1) + qJ(3);
t180 = qJD(1) * t118;
t122 = cos(qJ(4));
t116 = t122 ^ 2;
t120 = sin(qJ(4));
t179 = (t120 ^ 2 - t116) * MDP(11);
t178 = qJ(2) * MDP(6) + MDP(5) + MDP(7);
t119 = sin(qJ(5));
t121 = cos(qJ(5));
t152 = t121 * qJD(4);
t111 = qJD(5) * t152;
t157 = qJD(5) * t119;
t142 = t122 * t157;
t126 = -t120 * t152 - t142;
t91 = t126 * qJD(1) + t111;
t177 = t91 * t119;
t112 = qJD(1) * qJ(2) + qJD(3);
t108 = -pkin(6) * qJD(1) + t112;
t149 = qJD(1) * qJD(2);
t159 = qJD(4) * t120;
t94 = t108 * t159 - t122 * t149;
t176 = t94 * t119;
t175 = t94 * t121;
t161 = qJD(1) * t122;
t140 = t119 * t161;
t101 = t140 - t152;
t154 = t120 * qJD(1);
t110 = qJD(5) + t154;
t173 = t101 * t110;
t172 = t101 * t122;
t139 = t121 * t161;
t155 = t119 * qJD(4);
t103 = t139 + t155;
t171 = t103 * t110;
t170 = t103 * t122;
t169 = t110 * t119;
t168 = t110 * t121;
t167 = t119 * t120;
t166 = t120 * t108;
t165 = t120 * t121;
t164 = t122 * t108;
t123 = qJD(4) ^ 2;
t124 = qJD(1) ^ 2;
t162 = -t123 - t124;
t117 = -pkin(6) + qJ(2);
t160 = qJD(4) * t117;
t158 = qJD(4) * t122;
t156 = qJD(5) * t121;
t153 = t121 * MDP(22);
t151 = t122 * MDP(21);
t109 = -qJD(2) + t180;
t150 = qJD(2) - t109;
t148 = qJD(1) * qJD(4);
t146 = 0.2e1 * qJD(3) * qJD(1);
t145 = t110 * t167;
t144 = t110 * t165;
t143 = t110 * t157;
t141 = t110 * t156;
t138 = t120 * t148;
t137 = t122 * t148;
t136 = qJD(4) * t151;
t98 = qJD(4) * pkin(7) + t166;
t135 = t110 * t117 + t98;
t134 = pkin(4) * t122 + pkin(7) * t120;
t133 = t119 * t137 + t141;
t105 = t120 * pkin(4) - t122 * pkin(7) + t118;
t96 = t105 * qJD(1) - qJD(2);
t132 = t119 * t98 - t121 * t96;
t90 = t119 * t96 + t121 * t98;
t131 = qJD(2) + t109 + t180;
t130 = qJD(1) * t116 - t110 * t120;
t129 = -t117 * t123 + t146;
t99 = -qJD(4) * pkin(4) - t164;
t128 = -pkin(7) * t158 + t120 * t99;
t127 = (t119 * MDP(23) - t153) * t110;
t100 = t134 * qJD(4) + qJD(3);
t95 = t108 * t158 + t120 * t149;
t125 = -qJD(2) * t110 - t99 * qJD(4) - qJD(5) * t96 - t95;
t107 = t119 * t138;
t104 = t134 * qJD(1);
t97 = t100 * qJD(1);
t93 = t121 * t97;
t92 = t103 * qJD(5) - t107;
t1 = [MDP(8) * t146 + (t112 * qJD(2) + t109 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t118) * qJD(1)) * MDP(9) + 0.2e1 * t148 * t179 - t123 * t122 * MDP(13) + t131 * t158 * MDP(15) + (t129 * t122 - t131 * t159) * MDP(16) + (t91 * t121 * t122 + t126 * t103) * MDP(17) + ((t101 * t121 + t103 * t119) * t159 + (-t177 - t121 * t92 + (t101 * t119 - t103 * t121) * qJD(5)) * t122) * MDP(18) + (-t110 * t142 + (t130 * t121 + t170) * qJD(4)) * MDP(19) + (-t122 * t141 + (-t130 * t119 - t172) * qJD(4)) * MDP(20) + (t110 + t154) * t136 + ((t121 * t100 - t105 * t157) * t110 + (t99 * t156 - qJD(2) * t101 - t117 * t92 + t176 + (-t117 * t169 + (t121 * t105 - t117 * t167) * qJD(1) - t132) * qJD(4)) * t122) * MDP(22) + (-(t119 * t100 + t105 * t156) * t110 + (-t99 * t157 - qJD(2) * t103 - t117 * t91 + t175 + (-t117 * t168 - (t119 * t105 + t117 * t165) * qJD(1) - t90) * qJD(4)) * t122) * MDP(23) + (-0.2e1 * MDP(10) * t137 - t123 * MDP(12) + t129 * MDP(15) + t91 * MDP(19) - t92 * MDP(20) + (t101 * t160 + t125 * t119 - t135 * t156 + t93) * MDP(22) + (t103 * t160 + (t135 * qJD(5) - t97) * t119 + t125 * t121) * MDP(23)) * t120 + 0.2e1 * t178 * t149; MDP(22) * t143 + t133 * MDP(23) - t178 * t124 + ((-qJD(3) - t112) * MDP(9) + (t145 + t172) * MDP(22) + (t144 + t170) * MDP(23) + (0.2e1 * t120 * MDP(16) + (-0.2e1 * MDP(15) - t153) * t122) * qJD(4)) * qJD(1); -t124 * MDP(8) + (t150 * MDP(9) + t127) * qJD(1) + (t162 * MDP(16) + (-t110 * t155 - t92) * MDP(22) + (-t110 * t152 - t91) * MDP(23)) * t122 + (t162 * MDP(15) + qJD(5) * t127 + ((t101 - t140) * MDP(22) + (t103 - t139) * MDP(23)) * qJD(4)) * t120; (t103 * t168 + t177) * MDP(17) + ((t91 - t173) * t121 + (-t92 - t171) * t119) * MDP(18) + ((t144 - t170) * qJD(1) + t133) * MDP(19) + (-t143 + (-t145 + (t101 + t152) * t122) * qJD(1)) * MDP(20) - t110 * qJD(1) * t151 + (-pkin(4) * t92 - t175 - (t121 * t104 - t119 * t164) * t110 - t101 * t166 + (-pkin(7) * t168 + t99 * t119) * qJD(5) + (t128 * t119 + t122 * t132) * qJD(1)) * MDP(22) + (-pkin(4) * t91 + t176 + (t119 * t104 + t121 * t164) * t110 - t103 * t166 + (pkin(7) * t169 + t99 * t121) * qJD(5) + (t128 * t121 + t90 * t122) * qJD(1)) * MDP(23) + (MDP(15) * t161 - MDP(16) * t154) * t150 + (t122 * t120 * MDP(10) - t179) * t124; t103 * t101 * MDP(17) + (-t101 ^ 2 + t103 ^ 2) * MDP(18) + (-t121 * t138 + t111 + t173) * MDP(19) + (t107 + t171) * MDP(20) + qJD(1) * t136 + (-t99 * t103 + t90 * t110 - t119 * t95 + t93) * MDP(22) + (t99 * t101 - t110 * t132 - t119 * t97 - t121 * t95) * MDP(23) + (-MDP(19) * t140 - t103 * MDP(20) - t90 * MDP(22) + t132 * MDP(23)) * qJD(5);];
tauc = t1;
