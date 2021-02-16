% Calculate Coriolis joint torque vector for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:42
% EndTime: 2021-01-15 15:22:46
% DurationCPUTime: 0.89s
% Computational Cost: add. (858->171), mult. (2149->223), div. (0->0), fcn. (1331->4), ass. (0->76)
t167 = sin(pkin(8));
t168 = cos(pkin(8));
t169 = sin(qJ(3));
t170 = cos(qJ(3));
t152 = t167 * t170 + t168 * t169;
t194 = qJD(2) * t152;
t208 = 0.2e1 * t194;
t207 = MDP(14) + MDP(17);
t206 = (t169 ^ 2 - t170 ^ 2) * MDP(6);
t189 = t170 * MDP(11);
t204 = t169 * MDP(10) + t189;
t186 = MDP(12) + MDP(16);
t144 = t194 ^ 2;
t202 = -qJ(4) - pkin(6);
t179 = qJD(3) * t202;
t141 = qJD(4) * t170 + t169 * t179;
t187 = qJD(3) * qJD(1);
t130 = qJD(2) * t141 + t170 * t187;
t175 = -qJD(4) * t169 + t170 * t179;
t173 = qJD(2) * t175 - t169 * t187;
t108 = t130 * t167 - t168 * t173;
t155 = t202 * t170;
t181 = t202 * t169;
t131 = -t155 * t167 - t168 * t181;
t201 = t108 * t131;
t197 = t168 * t170;
t151 = t167 * t169 - t197;
t200 = t108 * t151;
t182 = qJD(2) * t197;
t192 = qJD(2) * t169;
t145 = t167 * t192 - t182;
t163 = -pkin(3) * t170 - pkin(2);
t193 = qJD(2) * t163;
t154 = qJD(4) + t193;
t115 = pkin(4) * t145 - qJ(5) * t194 + t154;
t199 = t115 * t194;
t142 = qJD(1) * t169 - qJD(2) * t155;
t198 = t142 * t167;
t134 = t168 * t142;
t196 = t170 * MDP(5);
t109 = t130 * t168 + t167 * t173;
t140 = t170 * qJD(1) + qJD(2) * t181;
t137 = qJD(3) * pkin(3) + t140;
t117 = t167 * t137 + t134;
t190 = t169 * qJD(3);
t121 = t140 * t168 - t198;
t188 = qJD(5) - t121;
t185 = MDP(13) - MDP(18);
t184 = 0.2e1 * qJD(2);
t183 = pkin(3) * t192;
t180 = qJD(2) * t190;
t116 = t137 * t168 - t198;
t147 = t152 * qJD(3);
t138 = qJD(2) * t147;
t156 = t167 * t180;
t139 = qJD(3) * t182 - t156;
t159 = pkin(3) * t180;
t178 = pkin(4) * t138 - qJ(5) * t139 + t159;
t119 = t140 * t167 + t134;
t177 = qJD(3) * t119 - t108;
t120 = t141 * t167 - t168 * t175;
t122 = t168 * t141 + t167 * t175;
t132 = -t168 * t155 + t167 * t181;
t174 = t108 * t152 + t120 * t194 - t122 * t145 + t131 * t139 - t132 * t138;
t171 = qJD(3) ^ 2;
t162 = -pkin(3) * t168 - pkin(4);
t160 = pkin(3) * t167 + qJ(5);
t150 = qJD(3) * t197 - t167 * t190;
t124 = pkin(4) * t151 - qJ(5) * t152 + t163;
t118 = pkin(4) * t194 + qJ(5) * t145 + t183;
t113 = qJD(3) * qJ(5) + t117;
t112 = -qJD(3) * pkin(4) + qJD(5) - t116;
t111 = pkin(3) * t190 + pkin(4) * t147 - qJ(5) * t150 - qJD(5) * t152;
t107 = qJD(3) * qJD(5) + t109;
t106 = -qJD(5) * t194 + t178;
t1 = [(t109 * t152 - t116 * t147 + t117 * t150 + t200) * MDP(15) + (t107 * t152 + t112 * t147 + t113 * t150 + t200) * MDP(19) - t204 * t171 + (-t147 * t186 - t150 * t185) * qJD(3) + t207 * (-t138 * t152 + t139 * t151 - t145 * t150 + t147 * t194); (t138 * t163 + t147 * t154) * MDP(12) + (t139 * t163 + t150 * t154) * MDP(13) + (-t109 * t151 - t116 * t150 - t117 * t147 + t174) * MDP(14) + (t109 * t132 - t116 * t120 + t117 * t122 + t201) * MDP(15) + (t106 * t151 + t111 * t145 + t115 * t147 + t124 * t138) * MDP(16) + (-t107 * t151 + t112 * t150 - t113 * t147 + t174) * MDP(17) + (-t106 * t152 - t111 * t194 - t115 * t150 - t124 * t139) * MDP(18) + (t106 * t124 + t107 * t132 + t111 * t115 + t112 * t120 + t113 * t122 + t201) * MDP(19) + (t170 * MDP(7) - t169 * MDP(8) + (-MDP(10) * t170 + MDP(11) * t169) * pkin(6)) * t171 + (-t185 * t122 - t186 * t120 + (-pkin(2) * t189 - t206) * t184 + ((-MDP(10) * pkin(2) + t196) * t184 + ((qJD(2) * t151 + t145) * MDP(12) + MDP(13) * t208 + (t154 + t193) * MDP(15)) * pkin(3)) * t169) * qJD(3); (-t145 * t183 - t154 * t194 + t177) * MDP(12) + (qJD(3) * t121 + t145 * t154 - t183 * t194 - t109) * MDP(13) + ((t117 - t119) * t194 + (-t116 + t121) * t145 + (-t138 * t167 - t139 * t168) * pkin(3)) * MDP(14) + (t116 * t119 - t117 * t121 + (-t108 * t168 + t109 * t167 - t154 * t192) * pkin(3)) * MDP(15) + (-t118 * t145 + t177 - t199) * MDP(16) + (-t138 * t160 + t139 * t162 + (t113 - t119) * t194 + (t112 - t188) * t145) * MDP(17) + (-t115 * t145 + t118 * t194 + (0.2e1 * qJD(5) - t121) * qJD(3) + t109) * MDP(18) + (t107 * t160 + t108 * t162 - t112 * t119 + t113 * t188 - t115 * t118) * MDP(19) + (pkin(2) * t204 - t169 * t196 + t206) * qJD(2) ^ 2; (t116 * t194 + t117 * t145 + t159) * MDP(15) + (t113 * t145 + (-qJD(5) - t112) * t194 + t178) * MDP(19) + t185 * (-t156 + (-t145 + t182) * qJD(3)) + t207 * (-t145 ^ 2 - t144) + t186 * qJD(3) * t208; t194 * t145 * MDP(16) + (-t156 + (t145 + t182) * qJD(3)) * MDP(17) + (-t144 - t171) * MDP(18) + (-qJD(3) * t113 + t108 + t199) * MDP(19);];
tauc = t1;
