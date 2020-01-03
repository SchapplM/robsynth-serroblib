% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR8
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:43
% EndTime: 2019-12-31 17:42:45
% DurationCPUTime: 0.40s
% Computational Cost: add. (535->104), mult. (1214->166), div. (0->0), fcn. (872->8), ass. (0->73)
t144 = qJD(2) + qJD(3);
t150 = sin(qJ(3));
t151 = sin(qJ(2));
t153 = cos(qJ(3));
t154 = cos(qJ(2));
t132 = t150 * t154 + t153 * t151;
t158 = t132 * qJD(2);
t176 = qJD(3) * t153;
t167 = t151 * t176;
t157 = (-t158 - t167) * qJD(1);
t172 = t154 * qJD(1);
t137 = qJD(2) * pkin(2) + t172;
t177 = qJD(3) * t137;
t168 = t150 * t177;
t196 = t157 - t168;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t195 = (t149 ^ 2 - t152 ^ 2) * MDP(10);
t155 = qJD(5) ^ 2;
t142 = t153 * pkin(2) + pkin(3);
t147 = sin(pkin(9));
t148 = cos(pkin(9));
t186 = t148 * t150;
t180 = pkin(2) * t186 + t147 * t142;
t128 = t132 * qJD(1);
t184 = t153 * t154;
t131 = -t150 * t151 + t184;
t129 = t131 * qJD(1);
t193 = pkin(2) * qJD(3);
t183 = t148 * t128 + t147 * t129 - (t147 * t153 + t186) * t193;
t194 = -t183 * t144 + (pkin(7) + t180) * t155;
t178 = qJD(1) * t151;
t166 = t150 * t178;
t122 = t153 * t137 - t166;
t119 = t144 * pkin(3) + t122;
t123 = t150 * t137 + t153 * t178;
t189 = t147 * t123;
t106 = t148 * t119 - t189;
t104 = -t144 * pkin(4) - t106;
t173 = t152 * qJD(5);
t165 = qJD(2) * t172;
t181 = t144 * t166;
t111 = (t165 + t177) * t153 - t181;
t98 = t147 * t111 - t148 * t196;
t192 = t104 * t173 + t98 * t149;
t191 = MDP(9) * t152;
t188 = t147 * t150;
t187 = t148 * t123;
t185 = t152 * t155;
t182 = t147 * t128 - t148 * t129 + (t148 * t153 - t188) * t193;
t175 = t149 * qJD(5);
t174 = t152 * MDP(15);
t171 = t155 * MDP(12);
t169 = MDP(11) * t185 + 0.2e1 * (t191 * t149 - t195) * qJD(5) * t144;
t99 = t148 * t111 + t196 * t147;
t164 = -t104 * t144 - t99;
t163 = (-pkin(2) * t144 - t137) * qJD(3);
t108 = t147 * t122 + t187;
t162 = -t108 * t144 + (t147 * pkin(3) + pkin(7)) * t155;
t109 = t148 * t122 - t189;
t161 = qJD(5) * ((-t148 * pkin(3) - pkin(4)) * t144 + t109);
t160 = -pkin(2) * t188 + t148 * t142;
t159 = qJD(5) * ((-pkin(4) - t160) * t144 - t182);
t143 = t144 ^ 2;
t117 = -t132 * qJD(3) - t158;
t116 = t144 * t131;
t115 = t147 * t131 + t148 * t132;
t114 = -t148 * t131 + t147 * t132;
t107 = t147 * t119 + t187;
t102 = t104 * t175;
t101 = t148 * t116 + t147 * t117;
t100 = t147 * t116 - t148 * t117;
t1 = [(-t106 * t100 + t107 * t101 + t98 * t114 + t99 * t115) * MDP(8) + (-t101 * t175 - t115 * t185) * MDP(14) + (t155 * t149 * t115 - t101 * t173) * MDP(15) + (-t151 * MDP(3) - t154 * MDP(4)) * qJD(2) ^ 2 + (t117 * MDP(6) - t116 * MDP(7) + (-t100 * t152 + t114 * t175) * MDP(14) + (t100 * t149 + t114 * t173) * MDP(15)) * t144; (t128 * t144 + t150 * t163 + t157) * MDP(6) + (t129 * t144 + (t163 - t165) * t153 + t181) * MDP(7) + (t183 * t106 + t182 * t107 - t98 * t160 + t99 * t180) * MDP(8) - t149 * t171 + (t102 + t149 * t159 + (-t194 - t98) * t152) * MDP(14) + (t194 * t149 + t152 * t159 + t192) * MDP(15) + t169; (t123 * t144 - t168) * MDP(6) + (t122 * t144 - t137 * t176 + t181) * MDP(7) + (t106 * t108 - t107 * t109 + (t147 * t99 - t148 * t98) * pkin(3)) * MDP(8) + t102 * MDP(14) + t192 * MDP(15) + (-MDP(6) * t167 + (-t132 * MDP(6) - MDP(7) * t184) * qJD(2)) * qJD(1) + ((-t162 - t98) * MDP(14) + MDP(15) * t161) * t152 + (MDP(14) * t161 + t162 * MDP(15) - t171) * t149 + t169; (-t149 * MDP(14) - t174) * t155; t164 * t174 + t143 * t195 + (t164 * MDP(14) - t143 * t191) * t149;];
tauc = t1;
