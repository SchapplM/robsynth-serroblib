% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:25
% EndTime: 2019-12-31 17:08:28
% DurationCPUTime: 0.98s
% Computational Cost: add. (384->139), mult. (956->210), div. (0->0), fcn. (520->4), ass. (0->79)
t151 = sin(qJ(2));
t152 = cos(qJ(4));
t150 = sin(qJ(4));
t153 = cos(qJ(2));
t185 = t150 * t153;
t123 = t151 * t152 - t185;
t201 = qJD(2) * t185 + qJD(4) * t123;
t174 = qJD(1) * qJD(2);
t170 = t151 * t174;
t104 = qJD(1) * t201 - t152 * t170;
t122 = t150 * t151 + t152 * t153;
t115 = t122 * qJD(1);
t178 = qJD(1) * t153;
t179 = qJD(1) * t151;
t117 = -t150 * t178 + t152 * t179;
t145 = qJD(2) - qJD(4);
t200 = -(t117 * t145 + t104) * MDP(18) - (t115 ^ 2 - t117 ^ 2) * MDP(16) + t117 * t115 * MDP(15);
t199 = -0.2e1 * qJD(1);
t198 = pkin(5) * MDP(14);
t148 = t151 ^ 2;
t197 = (-t153 ^ 2 + t148) * MDP(5);
t160 = t122 * qJD(4);
t195 = qJD(4) + t145;
t193 = pkin(2) + pkin(3);
t192 = pkin(5) - pkin(6);
t191 = (pkin(1) * MDP(9));
t190 = pkin(1) * MDP(10);
t189 = qJD(2) * pkin(2);
t188 = qJ(3) * t153;
t187 = t115 * t145;
t156 = qJD(1) ^ 2;
t183 = t153 * t156;
t169 = t153 * t174;
t182 = t150 * t170 + t152 * t169;
t164 = pkin(2) * t151 - t188;
t176 = qJD(3) * t151;
t114 = qJD(2) * t164 - t176;
t109 = qJD(1) * t114;
t168 = qJ(3) * t151 + pkin(1);
t131 = -pkin(2) * t153 - t168;
t119 = qJD(1) * t131;
t147 = qJD(2) * qJ(3);
t177 = qJD(2) * t151;
t143 = pkin(5) * t179;
t175 = -pkin(6) * t179 + qJD(3) + t143;
t134 = t192 * t153;
t167 = MDP(12) + t198;
t166 = qJD(3) - t189;
t165 = -0.2e1 * t109;
t126 = t192 * t177;
t144 = pkin(5) * t178;
t127 = -pkin(6) * t178 + t144;
t111 = -qJD(2) * t193 + t175;
t118 = t127 + t147;
t163 = t152 * t111 - t150 * t118;
t162 = -t150 * t111 - t152 * t118;
t161 = -t151 * t193 + t188;
t121 = t153 * t193 + t168;
t110 = t121 * qJD(1);
t146 = qJD(2) * qJD(3);
t112 = -qJD(1) * t126 + t146;
t142 = pkin(5) * t169;
t120 = -pkin(6) * t169 + t142;
t158 = t110 * t115 - t152 * t112 - t150 * t120;
t157 = -t110 * t117 - t150 * t112 + t152 * t120;
t103 = -qJD(1) * t160 + t182;
t108 = qJD(2) * t161 + t176;
t155 = qJD(2) ^ 2;
t133 = t192 * t151;
t132 = t144 + t147;
t130 = t143 + t166;
t129 = -pkin(5) * t170 + t146;
t128 = qJD(2) * t134;
t124 = t164 * qJD(1);
t113 = t161 * qJD(1);
t107 = t108 * qJD(1);
t106 = qJD(2) * t122 - t160;
t105 = -t152 * t177 + t201;
t1 = [(t109 * t131 + t114 * t119) * MDP(14) + (t103 * t123 + t106 * t117) * MDP(15) + (-t103 * t122 - t104 * t123 - t105 * t117 - t106 * t115) * MDP(16) + (t121 * t104 + t110 * t105 + t107 * t122 + t108 * t115) * MDP(20) + (t121 * t103 + t110 * t106 + t107 * t123 + t108 * t117) * MDP(21) + (-t106 * MDP(17) + t105 * MDP(18) + (-t126 * t150 - t128 * t152) * MDP(20) + (-t126 * t152 + t128 * t150) * MDP(21) + ((t133 * t150 + t134 * t152) * MDP(20) + (t133 * t152 - t134 * t150) * MDP(21)) * qJD(4)) * t145 + (t165 * MDP(13) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(5)) * t155) * t151 + (t155 * MDP(6) + t165 * MDP(11) + t129 * MDP(12) + (t129 * MDP(14) + (-MDP(11) - MDP(9)) * t155) * pkin(5)) * t153 + (t197 * t199 + (-0.2e1 * t119 * MDP(13) + t167 * t130 + t190 * t199) * t153 + (t119 * MDP(11) - t167 * t132 + (t131 * MDP(11) - (2 * t191) + (pkin(5) * t167 + (2 * MDP(4))) * t153) * qJD(1)) * t151) * qJD(2); t156 * t197 + t183 * t190 + 0.2e1 * t146 * MDP(13) + (qJ(3) * t129 + qJD(3) * t132 - t119 * t124) * MDP(14) + (-t182 + t187) * MDP(17) + (-t113 * t115 + (t152 * t127 + t175 * t150) * t145 + (-(-qJ(3) * t152 + t150 * t193) * t145 - t162) * qJD(4) - t157) * MDP(20) + (-t113 * t117 + (-t150 * t127 + t175 * t152) * t145 + ((-qJ(3) * t150 - t152 * t193) * t145 + t163) * qJD(4) - t158) * MDP(21) + (-MDP(4) * t183 + t156 * t191) * t151 + ((-t119 * t151 + t124 * t153) * MDP(11) + ((t132 - t147) * t151 + (-t130 + t166) * t153) * MDP(12) + (t119 * t153 + t124 * t151) * MDP(13) + (t132 * t151 + (-t130 - t189) * t153) * t198 + MDP(17) * t160) * qJD(1) - t200; (-t148 * t156 - t155) * MDP(13) + (-qJD(2) * t132 + t142) * MDP(14) + (-MDP(11) * t183 + (t119 * MDP(14) - MDP(20) * t115 - MDP(21) * t117) * qJD(1)) * t151 + (-MDP(20) * t150 - MDP(21) * t152) * t145 ^ 2; (t103 - t187) * MDP(17) + (t162 * t195 + t157) * MDP(20) + (-t163 * t195 + t158) * MDP(21) + t200;];
tauc = t1;
