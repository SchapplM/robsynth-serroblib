% Calculate Coriolis joint torque vector for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:55
% EndTime: 2021-01-15 20:08:59
% DurationCPUTime: 0.90s
% Computational Cost: add. (887->173), mult. (1723->236), div. (0->0), fcn. (927->6), ass. (0->102)
t182 = sin(qJ(4));
t184 = cos(qJ(4));
t177 = qJD(1) + qJD(2);
t185 = cos(qJ(2));
t231 = pkin(1) * qJD(1);
t201 = t185 * t231;
t161 = t177 * pkin(2) + t201;
t181 = cos(pkin(8));
t183 = sin(qJ(2));
t202 = t183 * t231;
t166 = t181 * t202;
t180 = sin(pkin(8));
t140 = t180 * t161 + t166;
t135 = t177 * pkin(7) + t140;
t197 = qJ(5) * t177 + t135;
t191 = t197 * t184;
t125 = t182 * qJD(3) + t191;
t237 = t125 * qJD(4);
t190 = t183 * MDP(5) + t185 * MDP(6);
t230 = pkin(1) * qJD(2);
t236 = t190 * t230;
t178 = t182 ^ 2;
t179 = t184 ^ 2;
t235 = (t178 - t179) * MDP(9);
t234 = pkin(4) * t178;
t233 = t181 * pkin(2);
t232 = t184 * pkin(4);
t229 = qJD(4) * pkin(4);
t228 = t125 * t184;
t165 = t180 * t202;
t139 = t181 * t161 - t165;
t134 = -t177 * pkin(3) - t139;
t227 = t134 * t177;
t176 = t177 ^ 2;
t226 = t176 * t184;
t225 = t177 * t182;
t224 = t180 * t183;
t223 = t181 * t183;
t186 = qJD(4) ^ 2;
t221 = t184 * t186;
t172 = t185 * pkin(1) + pkin(2);
t213 = pkin(1) * t223 + t180 * t172;
t152 = pkin(7) + t213;
t219 = -qJ(5) - t152;
t169 = t180 * pkin(2) + pkin(7);
t218 = -qJ(5) - t169;
t174 = t184 * qJD(3);
t124 = -t197 * t182 + t174;
t123 = t124 + t229;
t217 = t123 - t124;
t199 = -pkin(3) - t232;
t128 = t199 * t177 + qJD(5) - t139;
t154 = (t180 * t185 + t223) * t230;
t149 = qJD(1) * t154;
t206 = t182 * qJD(4);
t198 = t177 * t206;
t133 = pkin(4) * t198 + t149;
t205 = t184 * qJD(4);
t216 = t128 * t205 + t133 * t182;
t215 = t134 * t205 + t149 * t182;
t153 = t180 * t201 + t166;
t155 = t181 * t201 - t165;
t214 = t153 * t184 * t177 + t155 * t206;
t212 = -t178 - t179;
t210 = MDP(16) * t177;
t209 = MDP(17) * t177;
t208 = qJD(5) * t177;
t203 = -MDP(14) - MDP(16);
t200 = 0.2e1 * t184 * MDP(8) * t198 - 0.2e1 * t177 * qJD(4) * t235 + MDP(10) * t221;
t195 = -pkin(1) * t224 + t181 * t172;
t151 = -pkin(3) - t195;
t156 = (t181 * t185 - t224) * t230;
t196 = t151 * t177 - t156;
t194 = qJD(4) * t219;
t193 = qJD(4) * t218;
t150 = qJD(1) * t156;
t192 = qJD(4) * qJD(3) + t150;
t189 = t123 * t182 - t228;
t188 = t152 * t186 + t154 * t177;
t187 = (-qJD(5) - t128) * t177 - t192;
t175 = t184 * qJ(5);
t173 = t184 * qJD(5);
t170 = -pkin(3) - t233;
t162 = t199 - t233;
t159 = t184 * t169 + t175;
t158 = t218 * t182;
t147 = -t182 * qJD(5) + t184 * t193;
t146 = t182 * t193 + t173;
t145 = t151 - t232;
t144 = t155 * t205;
t142 = pkin(4) * t206 + t154;
t137 = t184 * t152 + t175;
t136 = t219 * t182;
t131 = t135 * t206;
t129 = t134 * t206;
t126 = t128 * t206;
t122 = (-qJD(5) - t156) * t182 + t184 * t194;
t121 = t184 * t156 + t182 * t194 + t173;
t120 = (-t150 - t208) * t182 - t237;
t119 = -qJ(5) * t198 - t131 + (t192 + t208) * t184;
t118 = t119 * t184;
t1 = [(-t139 * t154 + t140 * t156 - t149 * t195 + t150 * t213) * MDP(7) + t129 * MDP(13) + t215 * MDP(14) + t126 * MDP(15) + t216 * MDP(16) + t118 * MDP(17) + (t119 * t137 + t120 * t136 + t125 * t121 + t123 * t122 + t128 * t142 + t133 * t145) * MDP(18) + ((-t149 - t188) * MDP(13) + (-t142 * t177 - t133) * MDP(15) + t121 * t209) * t184 + (-t186 * MDP(11) + t188 * MDP(14) + t142 * t210 + (-t122 * t177 - t120) * MDP(17)) * t182 + (t122 * MDP(15) - t121 * MDP(16) + (t196 * MDP(14) + t145 * t210 + (-t136 * t177 - t123) * MDP(17)) * t184 + (t196 * MDP(13) + t145 * t177 * MDP(15) + (-t137 * t177 - t125) * MDP(17)) * t182) * qJD(4) + t200 + (-qJD(1) - t177) * t236; (t139 * t153 - t140 * t155 + (-t149 * t181 + t150 * t180) * pkin(2)) * MDP(7) + (-t149 * t184 - t169 * t221 + t129 + t214) * MDP(13) + (t144 + t215) * MDP(14) + (t147 * qJD(4) - t133 * t184 + t126 + t214) * MDP(15) + (-t146 * qJD(4) + t144 + t216) * MDP(16) + (-t123 * t205 + t118) * MDP(17) + (t119 * t159 + t120 * t158 + t123 * t147 + t125 * t146 - t128 * t153 + t133 * t162 - t155 * t228) * MDP(18) - qJD(1) * t236 + ((-t120 - t237) * MDP(17) + (t123 * t155 + t128 * t229) * MDP(18) + (t169 * MDP(14) - MDP(11)) * t186) * t182 + ((t146 * t184 - t147 * t182 + t212 * t155) * MDP(17) + t203 * t182 * t153 + t190 * t231 + (MDP(16) * t234 + (t170 * MDP(14) + t162 * MDP(16) - t158 * MDP(17)) * t184 + (t170 * MDP(13) + (t162 - t232) * MDP(15) - t159 * MDP(17)) * t182) * qJD(4)) * t177 + t200; (-t189 * qJD(4) + t119 * t182 + t120 * t184) * MDP(18) + (t203 * t184 + (-MDP(13) - MDP(15)) * t182) * t186; -t182 * MDP(8) * t226 + t176 * t235 + (-t150 - t227) * t182 * MDP(13) + (t131 + (-t182 * t135 + t174) * qJD(4) + (-t192 - t227) * t184) * MDP(14) + ((t125 - t191) * qJD(4) + (pkin(4) * t226 + t187) * t182) * MDP(15) + (-t176 * t234 + t131 + (qJ(5) * t225 + t124) * qJD(4) + t187 * t184) * MDP(16) + (t217 - t229) * t184 * t209 + (t217 * t125 + (-t128 * t225 + t120) * pkin(4)) * MDP(18); t133 * MDP(18) + t212 * MDP(17) * t176 + (t189 * MDP(18) + 0.2e1 * (t182 * MDP(15) + t184 * MDP(16)) * qJD(4)) * t177;];
tauc = t1;
