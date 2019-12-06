% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP5
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:44
% EndTime: 2019-12-05 15:38:47
% DurationCPUTime: 0.91s
% Computational Cost: add. (767->167), mult. (1956->222), div. (0->0), fcn. (1351->6), ass. (0->79)
t164 = sin(pkin(8));
t165 = cos(pkin(8));
t166 = sin(qJ(4));
t168 = cos(qJ(4));
t149 = t164 * t166 - t168 * t165;
t150 = t164 * t168 + t165 * t166;
t167 = sin(qJ(2));
t156 = qJD(2) * qJ(3) + qJD(1) * t167;
t199 = t164 ^ 2 + t165 ^ 2;
t181 = t199 * MDP(8);
t177 = t156 * t181;
t182 = t199 * MDP(7);
t169 = cos(qJ(2));
t171 = t149 * t169;
t208 = pkin(6) + qJ(3);
t152 = t208 * t164;
t153 = t208 * t165;
t176 = -t152 * t168 - t153 * t166;
t201 = qJD(1) * t171 - qJD(3) * t149 + qJD(4) * t176;
t125 = -t152 * t166 + t153 * t168;
t172 = t150 * t169;
t200 = -qJD(1) * t172 + qJD(3) * t150 + qJD(4) * t125;
t197 = qJD(2) * t150;
t147 = t150 * qJD(4);
t188 = MDP(14) + MDP(16);
t187 = MDP(15) - MDP(18);
t209 = t197 ^ 2;
t207 = qJD(2) * pkin(2);
t183 = pkin(6) * qJD(2) + t156;
t133 = t183 * t165;
t206 = t133 * t166;
t198 = qJD(1) * t169;
t196 = qJD(2) * t166;
t195 = qJD(2) * t167;
t194 = qJD(2) * t168;
t193 = qJD(4) * t166;
t192 = qJD(4) * t168;
t184 = t165 * t194;
t186 = t164 * t196;
t142 = -t184 + t186;
t160 = -pkin(3) * t165 - pkin(2);
t178 = qJD(3) - t198;
t148 = qJD(2) * t160 + t178;
t111 = pkin(4) * t142 - qJ(5) * t197 + t148;
t191 = t111 * MDP(19);
t190 = t142 * MDP(16);
t132 = t183 * t164;
t116 = -t132 * t168 - t206;
t189 = qJD(5) - t116;
t185 = t164 * t193;
t180 = t116 + t206;
t179 = qJD(2) * t199;
t151 = (qJD(3) + t198) * qJD(2);
t108 = -t132 * t193 + t133 * t192 + t150 * t151;
t117 = -t132 * t166 + t133 * t168;
t175 = t132 * t192 + t149 * t151;
t155 = qJD(4) * t184;
t134 = qJD(2) * t185 - t155;
t135 = qJD(2) * t147;
t161 = qJD(1) * t195;
t174 = pkin(4) * t135 + qJ(5) * t134 + t161;
t173 = -t111 * t197 - t108;
t137 = t149 * t167;
t170 = qJD(2) ^ 2;
t154 = t178 - t207;
t146 = -t165 * t192 + t185;
t138 = t142 ^ 2;
t136 = t150 * t167;
t123 = pkin(4) * t149 - qJ(5) * t150 + t160;
t122 = pkin(4) * t197 + qJ(5) * t142;
t121 = t155 + (t142 - t186) * qJD(4);
t119 = qJD(2) * t172 - qJD(4) * t137;
t118 = -qJD(2) * t171 - t167 * t147;
t113 = qJD(4) * qJ(5) + t117;
t112 = -qJD(4) * pkin(4) + t189;
t110 = pkin(4) * t147 + qJ(5) * t146 - qJD(5) * t150;
t109 = -qJD(5) * t197 + t174;
t107 = (qJD(5) - t206) * qJD(4) - t175;
t1 = [(-t118 * t142 + t119 * t197 - t134 * t136 + t135 * t137) * MDP(17) + (-t107 * t137 + t108 * t136 + t112 * t119 + t113 * t118) * MDP(19) + t188 * t142 * t195 + (-t118 * t187 - t119 * t188) * qJD(4) + (-t109 * MDP(19) - t188 * t135 + t187 * t134 + qJD(2) * t177 + (-MDP(4) + t182) * t170) * t169 + (t151 * t181 + (-t165 * MDP(5) + t164 * MDP(6) - MDP(3)) * t170 + ((t154 - t198) * MDP(8) + t191 + t187 * t197) * qJD(2)) * t167; (-t134 * t150 - t146 * t197) * MDP(9) + (t134 * t149 - t135 * t150 + t142 * t146 - t147 * t197) * MDP(10) + (t135 * t160 + t147 * t148) * MDP(14) + (-t134 * t160 - t146 * t148) * MDP(15) + (t109 * t149 + t110 * t142 + t111 * t147 + t123 * t135) * MDP(16) + (-t107 * t149 + t108 * t150 - t112 * t146 - t113 * t147 - t125 * t135 + t134 * t176 - t142 * t201 + t197 * t200) * MDP(17) + (-t109 * t150 - t110 * t197 + t111 * t146 + t123 * t134) * MDP(18) + (t107 * t125 - t108 * t176 + t109 * t123 + t110 * t111 + t200 * t112 + t201 * t113) * MDP(19) + (qJ(3) * t181 + t182) * t151 + (qJD(2) * t182 + t177) * qJD(3) + (-t146 * MDP(11) - t147 * MDP(12) - t187 * t201 - t188 * t200) * qJD(4) + ((-MDP(7) * t179 - t177) * t169 + ((-t154 - t207) * MDP(8) + (qJD(2) * t149 - t142) * MDP(14) - t190 + t197 * MDP(18) - t191) * t167) * qJD(1); (-t156 * t179 + t161) * MDP(8) + (-t138 - t209) * MDP(17) + (t113 * t142 + (-qJD(5) - t112) * t197 + t174) * MDP(19) - t170 * t182 - t187 * (-t155 + (t142 + t186) * qJD(4)) + 0.2e1 * t188 * t197 * qJD(4); (-t138 + t209) * MDP(10) + t121 * MDP(11) + (-t148 * t197 - t108) * MDP(14) + t175 * MDP(15) + t173 * MDP(16) + (pkin(4) * t134 - qJ(5) * t135 - (-t113 + t117) * t197) * MDP(17) + (t122 * t197 - t175) * MDP(18) + (-pkin(4) * t108 + qJ(5) * t107 - t111 * t122 - t112 * t117 + t113 * t189) * MDP(19) + (t197 * MDP(9) + t148 * MDP(15) - t122 * MDP(16) + (t112 - t189) * MDP(17) - t111 * MDP(18)) * t142 + ((-t164 * t194 - t165 * t196 + t197) * MDP(12) + t180 * MDP(15) + (0.2e1 * qJD(5) - t180) * MDP(18) + t188 * t117) * qJD(4); t197 * t190 + t121 * MDP(17) + (-qJD(4) ^ 2 - t209) * MDP(18) + (-qJD(4) * t113 - t173) * MDP(19);];
tauc = t1;
