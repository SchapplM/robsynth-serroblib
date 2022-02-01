% Calculate Coriolis joint torque vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:44
% EndTime: 2022-01-20 09:51:47
% DurationCPUTime: 0.69s
% Computational Cost: add. (591->113), mult. (1236->162), div. (0->0), fcn. (773->8), ass. (0->72)
t164 = sin(pkin(9));
t166 = cos(pkin(9));
t187 = t164 ^ 2 + t166 ^ 2;
t165 = sin(pkin(8));
t169 = sin(qJ(2));
t199 = pkin(1) * qJD(1);
t185 = t169 * t199;
t150 = t165 * t185;
t167 = cos(pkin(8));
t171 = cos(qJ(2));
t198 = pkin(1) * qJD(2);
t152 = t167 * t171 * t198;
t126 = qJD(1) * t152 - qJD(2) * t150;
t163 = qJD(1) + qJD(2);
t119 = qJD(4) * t163 + t126;
t207 = t119 * t187;
t206 = t169 * MDP(5) + t171 * MDP(6);
t168 = sin(qJ(5));
t170 = cos(qJ(5));
t143 = t164 * t170 + t166 * t168;
t129 = t143 * t163;
t205 = t163 * t187;
t184 = t171 * t199;
t136 = t167 * t184 - t150;
t202 = t136 - qJD(4);
t201 = pkin(2) * t167;
t200 = pkin(4) * t166;
t197 = t164 * t168;
t195 = t166 * MDP(8);
t193 = t166 * t170;
t192 = t167 * t169;
t144 = t163 * pkin(2) + t184;
t120 = t144 * t167 - t150;
t174 = qJD(4) - t120;
t181 = -pkin(3) - t200;
t109 = t181 * t163 + t174;
t135 = (t165 * t171 + t192) * t198;
t125 = qJD(1) * t135;
t142 = -t193 + t197;
t137 = t142 * qJD(5);
t191 = -t109 * t137 + t125 * t143;
t138 = t143 * qJD(5);
t190 = t109 * t138 + t125 * t142;
t151 = t167 * t185;
t121 = t165 * t144 + t151;
t156 = pkin(1) * t171 + pkin(2);
t188 = pkin(1) * t192 + t165 * t156;
t153 = t165 * t169 * pkin(1);
t183 = t163 * t197;
t182 = t163 * t193;
t179 = t156 * t167 - t153;
t145 = qJD(5) * t182;
t122 = -qJD(5) * t183 + t145;
t123 = t163 * t138;
t127 = -t182 + t183;
t178 = (-t122 * t142 - t123 * t143 + t127 * t137 - t129 * t138) * MDP(12) + (t122 * t143 - t129 * t137) * MDP(11) + (-t137 * MDP(13) - t138 * MDP(14)) * qJD(5);
t177 = -pkin(3) - t179;
t173 = t187 * (qJ(4) * t163 + t121);
t172 = -qJD(2) * t153 + t152;
t159 = t166 * pkin(7);
t155 = pkin(2) * t165 + qJ(4);
t146 = t181 - t201;
t140 = t155 * t166 + t159;
t139 = (-pkin(7) - t155) * t164;
t134 = t165 * t184 + t151;
t133 = qJ(4) + t188;
t132 = qJD(4) + t172;
t124 = t177 - t200;
t118 = t133 * t166 + t159;
t117 = (-pkin(7) - t133) * t164;
t113 = -pkin(3) * t163 + t174;
t1 = [(-t120 * t135 + t121 * t172 - t125 * t179 + t126 * t188) * MDP(7) + (-t135 * t163 - t125) * t195 + (t132 * t205 + t207) * MDP(9) + (t113 * t135 + t125 * t177 + t173 * t132 + t133 * t207) * MDP(10) + (t124 * t123 + t135 * t127 + ((-t117 * t168 - t118 * t170) * qJD(5) - t143 * t132) * qJD(5) + t190) * MDP(16) + (t124 * t122 + t135 * t129 + ((-t117 * t170 + t118 * t168) * qJD(5) + t142 * t132) * qJD(5) + t191) * MDP(17) + t178 + t206 * (-qJD(1) - t163) * t198; (t120 * t134 - t121 * t136 + (-t125 * t167 + t126 * t165) * pkin(2)) * MDP(7) + (t134 * t163 - t125) * t195 + (-t202 * t205 + t207) * MDP(9) + (t125 * (-pkin(3) - t201) - t113 * t134 + t155 * t207 - t202 * t173) * MDP(10) + (t146 * t123 - t134 * t127 + ((-t139 * t168 - t140 * t170) * qJD(5) + t202 * t143) * qJD(5) + t190) * MDP(16) + (t146 * t122 - t134 * t129 + ((-t139 * t170 + t140 * t168) * qJD(5) - t202 * t142) * qJD(5) + t191) * MDP(17) + t178 + t206 * (-qJD(2) + t163) * t199; (-MDP(16) * t138 + MDP(17) * t137) * qJD(5); (-t173 * t163 + t125) * MDP(10) + t145 * MDP(17) - t187 * MDP(9) * t163 ^ 2 + (0.2e1 * t129 * MDP(16) + (-t127 - t183) * MDP(17)) * qJD(5); t129 * t127 * MDP(11) + (-t127 ^ 2 + t129 ^ 2) * MDP(12) + (t145 + (t127 - t183) * qJD(5)) * MDP(13) + (-t109 * t129 - t143 * t119) * MDP(16) + (t109 * t127 + t142 * t119) * MDP(17);];
tauc = t1;
