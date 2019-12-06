% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:38
% EndTime: 2019-12-05 15:51:42
% DurationCPUTime: 1.52s
% Computational Cost: add. (687->196), mult. (1839->302), div. (0->0), fcn. (1378->10), ass. (0->105)
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t205 = t175 * MDP(12);
t235 = t172 * MDP(11) + t205;
t167 = sin(pkin(10));
t168 = sin(pkin(5));
t169 = cos(pkin(10));
t173 = sin(qJ(2));
t176 = cos(qJ(2));
t145 = (t167 * t173 - t169 * t176) * t168;
t165 = t172 ^ 2;
t234 = (-t175 ^ 2 + t165) * MDP(7);
t174 = cos(qJ(5));
t171 = sin(qJ(5));
t213 = qJD(4) * t171;
t216 = qJD(2) * t172;
t152 = t174 * t216 + t213;
t233 = t152 * MDP(13);
t232 = pkin(2) * t169;
t217 = qJD(1) * t168;
t201 = t176 * t217;
t155 = qJD(2) * pkin(2) + t201;
t202 = t173 * t217;
t157 = t169 * t202;
t136 = t167 * t155 + t157;
t134 = qJD(2) * pkin(7) + t136;
t170 = cos(pkin(5));
t159 = qJD(1) * t170 + qJD(3);
t127 = t134 * t175 + t159 * t172;
t143 = qJD(2) * t145;
t138 = qJD(1) * t143;
t121 = qJD(4) * t127 - t138 * t172;
t231 = t121 * t171;
t230 = t121 * t174;
t198 = t171 * t216;
t204 = qJD(2) * qJD(4);
t197 = t175 * t204;
t203 = qJD(4) * qJD(5);
t219 = (t197 + t203) * t174;
t139 = -qJD(5) * t198 + t219;
t229 = t139 * t171;
t209 = qJD(5) * t174;
t199 = t172 * t209;
t211 = qJD(4) * t175;
t140 = t171 * t203 + (t171 * t211 + t199) * qJD(2);
t228 = t140 * t175;
t206 = t174 * qJD(4);
t150 = t198 - t206;
t215 = qJD(2) * t175;
t160 = -qJD(5) + t215;
t227 = t150 * t160;
t226 = t150 * t172;
t225 = t152 * t160;
t224 = t171 * t160;
t223 = t171 * t175;
t222 = t174 * t160;
t221 = t174 * t175;
t141 = t167 * t201 + t157;
t190 = pkin(4) * t172 - pkin(8) * t175;
t154 = t190 * qJD(4);
t220 = t141 - t154;
t161 = pkin(2) * t167 + pkin(7);
t214 = qJD(4) * t161;
t212 = qJD(4) * t172;
t210 = qJD(5) * t171;
t186 = t134 * t172 - t159 * t175;
t124 = -qJD(4) * pkin(4) + t186;
t208 = t124 * qJD(5);
t200 = t172 * t210;
t196 = MDP(17) * t212;
t125 = qJD(4) * pkin(8) + t127;
t195 = t160 * t161 + t125;
t194 = -t139 * t175 + t152 * t212;
t156 = t167 * t202;
t135 = t155 * t169 - t156;
t192 = t160 * t200;
t191 = t160 * t199;
t183 = -pkin(4) * t175 - pkin(8) * t172 - pkin(3);
t129 = t183 * qJD(2) - t135;
t119 = t125 * t174 + t129 * t171;
t189 = t125 * t171 - t129 * t174;
t146 = (t167 * t176 + t169 * t173) * t168;
t132 = t146 * t175 + t170 * t172;
t188 = t132 * t174 + t145 * t171;
t187 = -t132 * t171 + t145 * t174;
t131 = t146 * t172 - t170 * t175;
t184 = qJD(2) * t165 - t160 * t175;
t142 = qJD(2) * t146;
t137 = qJD(1) * t142;
t177 = qJD(4) ^ 2;
t182 = qJD(2) * t141 - t161 * t177 - t137;
t133 = -qJD(2) * pkin(3) - t135;
t144 = t169 * t201 - t156;
t181 = qJD(4) * (qJD(2) * (-pkin(3) - t232) + t133 + t144);
t180 = t184 * t171;
t120 = -t186 * qJD(4) - t138 * t175;
t179 = qJD(4) * t124 + qJD(5) * t129 - t144 * t160 + t120;
t178 = qJD(2) ^ 2;
t153 = t190 * qJD(2);
t149 = t183 - t232;
t130 = (qJD(1) * t146 + t154) * qJD(2);
t128 = t174 * t130;
t123 = -t131 * qJD(4) - t143 * t175;
t122 = t132 * qJD(4) - t143 * t172;
t1 = [(-t135 * t142 - t136 * t143 + t137 * t145 - t138 * t146) * MDP(5) + (-(-t188 * qJD(5) - t123 * t171 + t142 * t174) * t160 + t122 * t150 + t131 * t140) * MDP(18) + ((t187 * qJD(5) + t123 * t174 + t142 * t171) * t160 + t122 * t152 + t131 * t139) * MDP(19) + (-MDP(3) * t173 - MDP(4) * t176) * t178 * t168 + (-MDP(11) * t122 - MDP(12) * t123) * qJD(4) + ((-t175 * MDP(11) + t172 * MDP(12)) * t142 + (t145 * t205 + (t145 * MDP(11) + t187 * MDP(18) - t188 * MDP(19)) * t172) * qJD(4)) * qJD(2); (t135 * t141 - t136 * t144 + (-t137 * t169 - t138 * t167) * pkin(2)) * MDP(5) - 0.2e1 * t204 * t234 - t200 * t233 + (-t150 * t174 - t152 * t171) * t211 * MDP(14) + (t184 * t206 + t192 + t194) * MDP(15) + (t191 + t228 + (-t180 - t226) * qJD(4)) * MDP(16) - t215 * t196 + (-t196 + (t149 * t210 + t220 * t174) * MDP(18) + (t149 * t209 - t220 * t171) * MDP(19)) * t160 + (t177 * MDP(8) + t182 * MDP(11) + t181 * MDP(12) + t206 * t233 + (t150 * t214 + t179 * t171 + t195 * t209 - t128) * MDP(18) + (t152 * t214 + (-t195 * qJD(5) + t130) * t171 + t179 * t174) * MDP(19)) * t175 + (0.2e1 * MDP(6) * t197 - t177 * MDP(9) + t181 * MDP(11) - t182 * MDP(12) + t139 * t174 * MDP(13) + (-t229 - t140 * t174 + (t150 * t171 - t152 * t174) * qJD(5)) * MDP(14) + (t174 * t208 + t231 + t161 * t140 - t144 * t150 + (-t161 * t224 + (t149 * t174 - t161 * t223) * qJD(2) - t189) * qJD(4)) * MDP(18) + (-t171 * t208 + t230 + t161 * t139 - t144 * t152 + (-t161 * t222 - (t149 * t171 + t161 * t221) * qJD(2) - t119) * qJD(4)) * MDP(19)) * t172; (t191 - t228) * MDP(18) + (-t192 + t194) * MDP(19) - t235 * t177 + ((-t180 + t226) * MDP(18) - t184 * MDP(19) * t174) * qJD(4); (-t152 * t222 + t229) * MDP(13) + ((t139 + t227) * t174 + (-t140 + t225) * t171) * MDP(14) + (-t160 * t209 + (t160 * t221 + (-t152 + t213) * t172) * qJD(2)) * MDP(15) + (t160 * t210 + (-t160 * t223 + (t150 + t206) * t172) * qJD(2)) * MDP(16) + t160 * MDP(17) * t216 + (-pkin(4) * t140 - t230 + (t174 * t153 + t171 * t186) * t160 - t127 * t150 + (pkin(8) * t222 + t124 * t171) * qJD(5) + (t189 * t172 + (-pkin(8) * t212 - t124 * t175) * t171) * qJD(2)) * MDP(18) + (-pkin(4) * t139 + t231 - (t171 * t153 - t174 * t186) * t160 - t127 * t152 + (-pkin(8) * t224 + t124 * t174) * qJD(5) + (-t124 * t221 + (-pkin(8) * t206 + t119) * t172) * qJD(2)) * MDP(19) + (-t172 * t175 * MDP(6) + t234) * t178 + t235 * (-qJD(2) * t133 + t138); t150 * t233 + (-t150 ^ 2 + t152 ^ 2) * MDP(14) + (t219 - t227) * MDP(15) + (-t171 * t197 - t225) * MDP(16) + qJD(2) * t196 + (-t119 * t160 - t171 * t120 - t124 * t152 + t128) * MDP(18) + (-t174 * t120 + t124 * t150 - t171 * t130 + t160 * t189) * MDP(19) + (-MDP(15) * t198 - t152 * MDP(16) - t119 * MDP(18) + t189 * MDP(19)) * qJD(5);];
tauc = t1;
