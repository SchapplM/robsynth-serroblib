% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:24
% EndTime: 2019-12-31 18:07:29
% DurationCPUTime: 1.58s
% Computational Cost: add. (993->197), mult. (2288->276), div. (0->0), fcn. (1610->6), ass. (0->91)
t180 = cos(qJ(5));
t203 = t180 * qJD(4);
t175 = sin(pkin(8));
t179 = sin(qJ(4));
t209 = qJD(1) * t179;
t198 = t175 * t209;
t176 = cos(pkin(8));
t181 = cos(qJ(4));
t208 = qJD(1) * t181;
t200 = t176 * t208;
t150 = -t198 + t200;
t178 = sin(qJ(5));
t215 = t150 * t178;
t138 = -t203 + t215;
t154 = t175 * t181 + t176 * t179;
t186 = qJD(1) * t154;
t231 = qJD(5) + t186;
t235 = t138 * t231;
t140 = qJD(4) * t178 + t150 * t180;
t234 = t140 * t231;
t194 = t180 * t231;
t158 = qJD(4) * t198;
t206 = qJD(4) * t181;
t199 = t176 * t206;
t144 = qJD(1) * t199 - t158;
t216 = t144 * t178;
t233 = -t194 * t231 - t216;
t211 = t175 ^ 2 + t176 ^ 2;
t232 = qJD(3) * t211;
t230 = qJ(2) * MDP(6) + t175 * MDP(7) + t176 * MDP(8) + MDP(5);
t177 = -pkin(1) - qJ(3);
t228 = qJD(1) * t177;
t159 = qJD(2) + t228;
t197 = -pkin(6) * qJD(1) + t159;
t145 = t197 * t175;
t146 = t197 * t176;
t131 = t145 * t181 + t146 * t179;
t153 = t175 * t179 - t181 * t176;
t185 = t153 * qJD(3);
t120 = -qJD(1) * t185 + qJD(4) * t131;
t227 = t231 * (pkin(4) * t150 + t231 * pkin(7)) + t120;
t187 = t154 * qJD(3);
t190 = t145 * t179 - t146 * t181;
t119 = -qJD(1) * t187 - qJD(4) * t190;
t225 = -pkin(6) + t177;
t155 = t225 * t175;
t156 = t225 * t176;
t136 = t155 * t179 - t156 * t181;
t123 = -qJD(4) * t136 - t187;
t126 = -qJD(4) * pkin(4) + t190;
t174 = qJD(1) * qJ(2);
t168 = qJD(3) + t174;
t170 = t175 * pkin(3);
t157 = qJD(1) * t170 + t168;
t128 = pkin(4) * t186 - pkin(7) * t150 + t157;
t164 = qJ(2) + t170;
t134 = pkin(4) * t154 + pkin(7) * t153 + t164;
t137 = t155 * t181 + t156 * t179;
t207 = qJD(4) * t179;
t151 = -t175 * t206 - t176 * t207;
t226 = -(qJD(5) * t134 + t123) * t231 - t137 * t144 - (qJD(5) * t128 + t119) * t154 - t120 * t153 + t126 * t151;
t205 = qJD(5) * t178;
t143 = qJD(4) * t186;
t212 = qJD(5) * t203 - t180 * t143;
t121 = -t150 * t205 + t212;
t223 = t121 * t153;
t222 = t121 * t178;
t221 = t126 * t153;
t220 = t134 * t144;
t219 = t138 * t150;
t218 = t140 * t150;
t217 = t143 * t178;
t142 = t180 * t144;
t210 = qJD(1) * t153;
t204 = qJD(5) * t180;
t173 = qJD(1) * qJD(2);
t195 = qJD(1) * t211;
t127 = qJD(4) * pkin(7) + t131;
t118 = t127 * t180 + t128 * t178;
t191 = t127 * t178 - t128 * t180;
t189 = t142 + (-t178 * t186 - t205) * t231;
t188 = t151 * t180 + t153 * t205;
t183 = -pkin(7) * t144 + (t126 - t190) * t231;
t182 = qJD(1) ^ 2;
t152 = -t175 * t207 + t199;
t132 = pkin(4) * t152 - pkin(7) * t151 + qJD(2);
t129 = pkin(4) * t144 + pkin(7) * t143 + t173;
t125 = t180 * t129;
t124 = qJD(4) * t137 - t185;
t122 = qJD(5) * t140 - t217;
t1 = [0.2e1 * qJD(3) * MDP(9) * t195 + ((t168 + t174) * qJD(2) + (-t159 - t228) * t232) * MDP(10) + (t143 * t153 + t150 * t151) * MDP(11) + (t143 * t154 + t144 * t153 - t150 * t152 - t151 * t186) * MDP(12) + (0.2e1 * t186 * qJD(2) + t144 * t164 + t152 * t157) * MDP(16) + (-t143 * t164 + t151 * t157 + (t150 - t210) * qJD(2)) * MDP(17) + (t140 * t188 - t180 * t223) * MDP(18) + ((-t138 * t180 - t140 * t178) * t151 + (t222 + t122 * t180 + (-t138 * t178 + t140 * t180) * qJD(5)) * t153) * MDP(19) + (t121 * t154 + t140 * t152 - t142 * t153 + t188 * t231) * MDP(20) + (t153 * t216 - t122 * t154 - t138 * t152 + (-t151 * t178 + t153 * t204) * t231) * MDP(21) + (t144 * t154 + t152 * t231) * MDP(22) + (-t191 * t152 + t136 * t122 + t124 * t138 + t125 * t154 + (t132 * t231 + t220 + (-t127 * t154 - t137 * t231 - t221) * qJD(5)) * t180 + t226 * t178) * MDP(23) + (-t118 * t152 + t136 * t121 + t124 * t140 + (-(-qJD(5) * t137 + t132) * t231 - t220 - (-qJD(5) * t127 + t129) * t154 + qJD(5) * t221) * t178 + t226 * t180) * MDP(24) + (MDP(13) * t151 - MDP(14) * t152 - MDP(16) * t124 - MDP(17) * t123) * qJD(4) + 0.2e1 * t230 * t173; (t122 * t153 - t138 * t151 - t154 * t216) * MDP(23) + (-t140 * t151 - t154 * t142 + t223) * MDP(24) + (MDP(16) * t151 - MDP(17) * t152) * qJD(4) + ((-t168 - t232) * MDP(10) - t186 * MDP(16) - t150 * MDP(17)) * qJD(1) - t230 * t182 + ((-qJD(1) * t180 - t152 * t178 - t154 * t204) * MDP(23) + (qJD(1) * t178 - t152 * t180 + t154 * t205) * MDP(24)) * t231; (t159 * t195 + t173) * MDP(10) - t158 * MDP(16) + (t189 - t219) * MDP(23) + (-t218 + t233) * MDP(24) - t211 * MDP(9) * t182 + ((t150 + t200) * MDP(16) + (-t175 * t208 - t176 * t209 - t186) * MDP(17)) * qJD(4); t150 ^ 2 * MDP(12) + (t158 + (t150 - t200) * qJD(4)) * MDP(14) + (qJD(3) * t210 - t150 * t157) * MDP(16) + (t140 * t194 + t222) * MDP(18) + ((t121 - t235) * t180 + (-t122 - t234) * t178) * MDP(19) + (-t218 - t233) * MDP(20) + (t189 + t219) * MDP(21) - t231 * t150 * MDP(22) + (-pkin(4) * t122 - t131 * t138 + t150 * t191 + t183 * t178 - t180 * t227) * MDP(23) + (-pkin(4) * t121 + t118 * t150 - t131 * t140 + t178 * t227 + t183 * t180) * MDP(24) + (t150 * MDP(11) + (qJD(3) + t157) * MDP(17) - MDP(12) * t186) * t186; t140 * t138 * MDP(18) + (-t138 ^ 2 + t140 ^ 2) * MDP(19) + (t212 + t235) * MDP(20) + (t217 + t234) * MDP(21) + t144 * MDP(22) + (t118 * t231 - t119 * t178 - t126 * t140 + t125) * MDP(23) + (-t119 * t180 + t126 * t138 - t129 * t178 - t191 * t231) * MDP(24) + (-MDP(20) * t215 - MDP(21) * t140 - MDP(23) * t118 + MDP(24) * t191) * qJD(5);];
tauc = t1;
