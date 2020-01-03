% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:15
% EndTime: 2019-12-31 18:09:17
% DurationCPUTime: 0.72s
% Computational Cost: add. (945->159), mult. (2227->217), div. (0->0), fcn. (1368->6), ass. (0->78)
t162 = sin(pkin(7)) * pkin(1) + pkin(6);
t200 = qJ(4) + t162;
t210 = MDP(12) + MDP(15);
t173 = cos(qJ(3));
t209 = t173 * MDP(5);
t172 = sin(qJ(3));
t208 = (t172 ^ 2 - t173 ^ 2) * MDP(6);
t186 = t200 * qJD(1);
t136 = t173 * qJD(2) - t186 * t172;
t137 = qJD(2) * t172 + t186 * t173;
t168 = sin(pkin(8));
t170 = cos(pkin(8));
t153 = t168 * t173 + t170 * t172;
t148 = t153 * qJD(1);
t143 = t148 ^ 2;
t192 = qJD(1) * qJD(4);
t128 = t136 * qJD(3) + t173 * t192;
t176 = -t137 * qJD(3) - t172 * t192;
t106 = t128 * t168 - t170 * t176;
t151 = t200 * t173;
t188 = t200 * t172;
t126 = t151 * t168 + t170 * t188;
t206 = t106 * t126;
t203 = t170 * t173;
t152 = t168 * t172 - t203;
t205 = t106 * t152;
t204 = t137 * t168;
t132 = t170 * t137;
t174 = qJD(3) ^ 2;
t202 = t172 * t174;
t201 = t173 * t174;
t107 = t170 * t128 + t168 * t176;
t134 = qJD(3) * pkin(3) + t136;
t113 = t168 * t134 + t132;
t164 = -cos(pkin(7)) * pkin(1) - pkin(2);
t155 = qJD(1) * t164;
t197 = t172 * qJD(1);
t196 = t172 * qJD(3);
t195 = t173 * qJD(1);
t117 = t136 * t170 - t204;
t194 = qJD(5) - t117;
t193 = qJD(1) * qJD(3);
t191 = pkin(3) * t196;
t190 = t170 * t195;
t189 = t172 * t193;
t187 = qJD(3) * t200;
t112 = t134 * t170 - t204;
t184 = -pkin(3) * t173 + t164;
t183 = 0.2e1 * qJD(3) * t155;
t147 = t153 * qJD(3);
t139 = qJD(1) * t147;
t156 = t168 * t189;
t140 = qJD(3) * t190 - t156;
t159 = pkin(3) * t189;
t182 = pkin(4) * t139 - qJ(5) * t140 + t159;
t180 = t184 * qJD(1);
t144 = qJD(4) + t180;
t145 = t168 * t197 - t190;
t120 = pkin(4) * t145 - qJ(5) * t148 + t144;
t181 = t120 * t148 + t106;
t178 = -qJD(4) * t172 - t173 * t187;
t138 = qJD(4) * t173 - t172 * t187;
t118 = t138 * t168 - t170 * t178;
t119 = t170 * t138 + t168 * t178;
t127 = t170 * t151 - t168 * t188;
t177 = t106 * t153 + t118 * t148 - t119 * t145 + t126 * t140 - t127 * t139;
t163 = -pkin(3) * t170 - pkin(4);
t160 = pkin(3) * t168 + qJ(5);
t150 = qJD(3) * t203 - t168 * t196;
t122 = pkin(4) * t152 - qJ(5) * t153 + t184;
t121 = pkin(3) * t197 + pkin(4) * t148 + qJ(5) * t145;
t116 = t136 * t168 + t132;
t114 = pkin(4) * t147 - qJ(5) * t150 - qJD(5) * t153 + t191;
t111 = qJD(3) * qJ(5) + t113;
t110 = -qJD(3) * pkin(4) + qJD(5) - t112;
t109 = -qJD(5) * t148 + t182;
t105 = qJD(3) * qJD(5) + t107;
t1 = [0.2e1 * t189 * t209 - 0.2e1 * t193 * t208 + MDP(7) * t201 - MDP(8) * t202 + (-t162 * t201 + t172 * t183) * MDP(10) + (t162 * t202 + t173 * t183) * MDP(11) + (-t107 * t152 - t112 * t150 - t113 * t147 + t177) * MDP(12) + (t206 + t107 * t127 - t112 * t118 + t113 * t119 + (t144 + t180) * t191) * MDP(13) + (-qJD(3) * t118 + t109 * t152 + t114 * t145 + t120 * t147 + t122 * t139) * MDP(14) + (-t105 * t152 + t110 * t150 - t111 * t147 + t177) * MDP(15) + (qJD(3) * t119 - t109 * t153 - t114 * t148 - t120 * t150 - t122 * t140) * MDP(16) + (t105 * t127 + t109 * t122 + t110 * t118 + t111 * t119 + t114 * t120 + t206) * MDP(17); (t107 * t153 - t112 * t147 + t113 * t150 + t205) * MDP(13) + (t105 * t153 + t110 * t147 + t111 * t150 + t205) * MDP(17) + (-t172 * MDP(10) - t173 * MDP(11)) * t174 + (-t147 * MDP(14) + t150 * MDP(16)) * qJD(3) + t210 * (-t153 * t139 + t140 * t152 - t150 * t145 + t147 * t148); ((t113 - t116) * t148 + (-t112 + t117) * t145 + (-t139 * t168 - t140 * t170) * pkin(3)) * MDP(12) + (t112 * t116 - t113 * t117 + (-t106 * t170 + t107 * t168 - t144 * t197) * pkin(3)) * MDP(13) + (qJD(3) * t116 - t121 * t145 - t181) * MDP(14) + (-t139 * t160 + t140 * t163 + (t111 - t116) * t148 + (t110 - t194) * t145) * MDP(15) + (-t120 * t145 + t121 * t148 + (0.2e1 * qJD(5) - t117) * qJD(3) + t107) * MDP(16) + (t105 * t160 + t106 * t163 - t110 * t116 + t194 * t111 - t120 * t121) * MDP(17) + (-t172 * t209 + t208) * qJD(1) ^ 2 + (-MDP(10) * t197 - MDP(11) * t195) * t155; (t112 * t148 + t113 * t145 + t159) * MDP(13) + t156 * MDP(16) + (t111 * t145 + (-qJD(5) - t110) * t148 + t182) * MDP(17) + ((t168 * t195 + t170 * t197 + t148) * MDP(14) + (t145 - t190) * MDP(16)) * qJD(3) + t210 * (-t145 ^ 2 - t143); t148 * t145 * MDP(14) + (-t156 + (t145 + t190) * qJD(3)) * MDP(15) + (-t143 - t174) * MDP(16) + (-qJD(3) * t111 + t181) * MDP(17);];
tauc = t1;
