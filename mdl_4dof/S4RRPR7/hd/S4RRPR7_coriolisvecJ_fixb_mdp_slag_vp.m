% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:46
% EndTime: 2019-12-31 17:06:49
% DurationCPUTime: 1.07s
% Computational Cost: add. (835->179), mult. (2242->264), div. (0->0), fcn. (1535->6), ass. (0->95)
t202 = (qJD(1) * qJD(2));
t229 = -2 * t202;
t178 = sin(qJ(2));
t228 = t178 * MDP(4);
t180 = cos(qJ(2));
t227 = (t178 ^ 2 - t180 ^ 2) * MDP(5);
t224 = -qJ(3) - pkin(5);
t196 = qJD(2) * t224;
t149 = qJD(3) * t180 + t178 * t196;
t143 = t149 * qJD(1);
t175 = sin(pkin(7));
t176 = cos(pkin(7));
t186 = -t178 * qJD(3) + t180 * t196;
t184 = qJD(1) * t186;
t117 = t143 * t175 - t176 * t184;
t206 = qJD(1) * t178;
t213 = t176 * t180;
t150 = qJD(1) * t213 - t175 * t206;
t148 = qJD(4) - t150;
t159 = t175 * t180 + t176 * t178;
t152 = t159 * qJD(1);
t169 = pkin(2) * t175 + pkin(6);
t226 = (pkin(2) * t206 + pkin(3) * t152 - pkin(6) * t150 + qJD(4) * t169) * t148 + t117;
t118 = t176 * t143 + t175 * t184;
t200 = -pkin(2) * t180 - pkin(1);
t190 = t200 * qJD(1);
t163 = qJD(3) + t190;
t122 = -pkin(3) * t150 - pkin(6) * t152 + t163;
t124 = t176 * t149 + t175 * t186;
t199 = t224 * t178;
t161 = qJD(1) * t199;
t223 = qJD(2) * pkin(2);
t157 = t161 + t223;
t164 = t224 * t180;
t162 = qJD(1) * t164;
t214 = t162 * t175;
t131 = t157 * t176 + t214;
t127 = -qJD(2) * pkin(3) - t131;
t158 = t175 * t178 - t213;
t130 = pkin(3) * t158 - pkin(6) * t159 + t200;
t154 = t158 * qJD(2);
t136 = -t176 * t164 + t175 * t199;
t151 = t159 * qJD(2);
t145 = qJD(1) * t151;
t189 = t117 * t159 - t136 * t145;
t225 = -t127 * t154 - (qJD(4) * t130 + t124) * t148 - (qJD(4) * t122 + t118) * t158 + t189;
t177 = sin(qJ(4));
t204 = qJD(4) * t177;
t197 = t180 * t202;
t198 = t178 * t202;
t146 = -t175 * t198 + t176 * t197;
t179 = cos(qJ(4));
t203 = t179 * qJD(2);
t208 = qJD(4) * t203 + t179 * t146;
t119 = -t152 * t204 + t208;
t222 = t119 * t177;
t221 = t130 * t145;
t215 = t152 * t177;
t137 = -t203 + t215;
t220 = t137 * t148;
t219 = t137 * t152;
t139 = qJD(2) * t177 + t152 * t179;
t218 = t139 * t148;
t217 = t139 * t152;
t216 = t146 * t177;
t155 = t176 * t162;
t212 = t177 * t145;
t181 = qJD(2) ^ 2;
t211 = t178 * t181;
t141 = t179 * t145;
t210 = t180 * t181;
t182 = qJD(1) ^ 2;
t209 = t180 * t182;
t132 = t175 * t157 - t155;
t205 = qJD(4) * t159;
t201 = t178 * t223;
t195 = pkin(1) * t229;
t194 = t148 * t179;
t128 = qJD(2) * pkin(6) + t132;
t114 = t122 * t179 - t128 * t177;
t115 = t122 * t177 + t128 * t179;
t188 = t141 + (t150 * t177 - t204) * t148;
t187 = -t154 * t179 - t159 * t204;
t134 = t161 * t176 + t214;
t183 = -t169 * t145 + (t127 + t134) * t148;
t170 = -pkin(2) * t176 - pkin(3);
t167 = pkin(2) * t198;
t135 = -t164 * t175 - t176 * t199;
t133 = t161 * t175 - t155;
t126 = pkin(3) * t151 + pkin(6) * t154 + t201;
t123 = t149 * t175 - t176 * t186;
t121 = pkin(3) * t145 - pkin(6) * t146 + t167;
t120 = qJD(4) * t139 + t216;
t116 = t179 * t121;
t1 = [0.2e1 * t197 * t228 + t227 * t229 + MDP(6) * t210 - MDP(7) * t211 + (-pkin(5) * t210 + t178 * t195) * MDP(9) + (pkin(5) * t211 + t180 * t195) * MDP(10) + (-t118 * t158 + t123 * t152 + t124 * t150 + t131 * t154 - t132 * t151 + t135 * t146 + t189) * MDP(11) + (t117 * t135 + t118 * t136 - t131 * t123 + t132 * t124 + (t163 + t190) * t201) * MDP(12) + (t119 * t159 * t179 + t187 * t139) * MDP(13) + (-(-t137 * t179 - t139 * t177) * t154 + (-t222 - t120 * t179 + (t137 * t177 - t139 * t179) * qJD(4)) * t159) * MDP(14) + (t119 * t158 + t139 * t151 + t159 * t141 + t187 * t148) * MDP(15) + (-t159 * t212 - t120 * t158 - t137 * t151 + (t154 * t177 - t179 * t205) * t148) * MDP(16) + (t145 * t158 + t148 * t151) * MDP(17) + (t114 * t151 + t116 * t158 + t135 * t120 + t123 * t137 + (t126 * t148 + t221 + (t127 * t159 - t128 * t158 - t136 * t148) * qJD(4)) * t179 + t225 * t177) * MDP(18) + (-t115 * t151 + t135 * t119 + t123 * t139 + (-(-qJD(4) * t136 + t126) * t148 - t221 - (-qJD(4) * t128 + t121) * t158 - t127 * t205) * t177 + t225 * t179) * MDP(19); -t209 * t228 + t182 * t227 + ((t132 - t133) * t152 + (t131 - t134) * t150 + (-t145 * t175 - t146 * t176) * pkin(2)) * MDP(11) + (t131 * t133 - t132 * t134 + (-t117 * t176 + t118 * t175 - t163 * t206) * pkin(2)) * MDP(12) + (t139 * t194 + t222) * MDP(13) + ((t119 - t220) * t179 + (-t120 - t218) * t177) * MDP(14) + (t148 * t194 + t212 - t217) * MDP(15) + (t188 + t219) * MDP(16) - t148 * t152 * MDP(17) + (-t114 * t152 + t170 * t120 - t133 * t137 + t183 * t177 - t226 * t179) * MDP(18) + (t115 * t152 + t170 * t119 - t133 * t139 + t226 * t177 + t183 * t179) * MDP(19) + (t182 * t178 * MDP(9) + MDP(10) * t209) * pkin(1); (-t150 ^ 2 - t152 ^ 2) * MDP(11) + (t131 * t152 - t132 * t150 + t167) * MDP(12) + (t188 - t219) * MDP(18) + (-t148 ^ 2 * t179 - t212 - t217) * MDP(19); t139 * t137 * MDP(13) + (-t137 ^ 2 + t139 ^ 2) * MDP(14) + (t208 + t220) * MDP(15) + (-t216 + t218) * MDP(16) + t145 * MDP(17) + (t115 * t148 - t118 * t177 - t127 * t139 + t116) * MDP(18) + (t114 * t148 - t118 * t179 - t121 * t177 + t127 * t137) * MDP(19) + (-MDP(15) * t215 - t139 * MDP(16) - t115 * MDP(18) - t114 * MDP(19)) * qJD(4);];
tauc = t1;
