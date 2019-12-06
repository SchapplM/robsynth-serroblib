% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR2
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:25
% EndTime: 2019-12-05 18:20:28
% DurationCPUTime: 0.71s
% Computational Cost: add. (631->126), mult. (1358->194), div. (0->0), fcn. (762->8), ass. (0->82)
t158 = sin(pkin(8));
t162 = sin(qJ(2));
t211 = pkin(1) * qJD(1);
t186 = t162 * t211;
t143 = t158 * t186;
t160 = cos(pkin(8));
t164 = cos(qJ(2));
t210 = pkin(1) * qJD(2);
t145 = t160 * t164 * t210;
t129 = qJD(1) * t145 - qJD(2) * t143;
t154 = qJD(1) + qJD(2);
t123 = qJD(4) * t154 + t129;
t157 = sin(pkin(9));
t152 = t157 ^ 2;
t159 = cos(pkin(9));
t193 = t159 ^ 2 + t152;
t219 = t123 * t193;
t161 = sin(qJ(5));
t163 = cos(qJ(5));
t205 = t152 * t163;
t218 = -t161 * MDP(12) * t205 + (t161 ^ 2 - t163 ^ 2) * MDP(13) * t152;
t217 = t162 * MDP(5) + t164 * MDP(6);
t185 = t164 * t211;
t134 = t160 * t185 - t143;
t213 = qJD(4) - t134;
t212 = pkin(2) * t160;
t139 = t154 * pkin(2) + t185;
t144 = t160 * t186;
t125 = t158 * t139 + t144;
t120 = qJ(4) * t154 + t125;
t111 = -t159 * qJD(3) + t120 * t157;
t209 = t111 * t157;
t146 = t158 * t162 * pkin(1);
t167 = -qJD(2) * t146 + t145;
t130 = qJD(4) + t167;
t208 = t130 * t154;
t207 = t152 * t154;
t206 = t152 * t161;
t204 = t154 * t157;
t203 = t154 * t159;
t202 = t159 * MDP(8);
t201 = t159 * t161;
t200 = t159 * t163;
t199 = t160 * t162;
t142 = -qJD(5) + t203;
t198 = t161 * t142;
t197 = t163 * t142;
t168 = -pkin(4) * t159 - pkin(7) * t157 - pkin(3);
t124 = t139 * t160 - t143;
t172 = qJD(4) - t124;
t109 = t168 * t154 + t172;
t112 = qJD(3) * t157 + t120 * t159;
t133 = (t158 * t164 + t199) * t210;
t128 = qJD(1) * t133;
t196 = (t123 * t200 + t161 * t128 + (t109 * t163 - t112 * t161) * qJD(5)) * t159 + t123 * t205;
t149 = pkin(1) * t164 + pkin(2);
t194 = pkin(1) * t199 + t158 * t149;
t189 = qJD(5) * t161;
t188 = qJD(5) * t163;
t187 = qJD(5) + t142;
t184 = t111 * t204;
t183 = t154 * t205;
t182 = t157 * t188;
t181 = t157 * t189;
t179 = t149 * t160 - t146;
t138 = t181 * t203;
t178 = (t142 * t181 + t138) * MDP(14) + (t142 + t203) * MDP(15) * t182 + 0.2e1 * t218 * qJD(5) * t154;
t171 = -t109 * t161 - t112 * t163;
t108 = t171 * qJD(5) - t123 * t201 + t163 * t128;
t174 = -t108 * t159 + t111 * t182 + t123 * t206;
t170 = t112 * t159 + t209;
t169 = t142 * t159 + t207;
t148 = pkin(2) * t158 + qJ(4);
t166 = t148 * t188 + t213 * t161;
t151 = t154 ^ 2;
t136 = t168 - t212;
t132 = t158 * t185 + t144;
t131 = qJ(4) + t194;
t126 = t128 * t157;
t121 = t168 - t179;
t117 = -pkin(3) * t154 + t172;
t1 = [(-t124 * t133 + t125 * t167 - t128 * t179 + t129 * t194) * MDP(7) + (-t133 * t154 - t128) * t202 + (t133 * t204 + t126) * MDP(9) + (t193 * t208 + t219) * MDP(10) + (t128 * (-pkin(3) - t179) + t117 * t133 + t170 * t130 + t131 * t219) * MDP(11) + (-(-t130 * t201 + t133 * t163) * t142 + t206 * t208 + (-(-t121 * t161 - t131 * t200) * t142 + t131 * t183) * qJD(5) + t174) * MDP(17) + ((t130 * t200 + t133 * t161) * t142 + t130 * t183 + (t121 * t197 + (-t169 * t131 - t209) * t161) * qJD(5) + t196) * MDP(18) + t178 + t217 * (-qJD(1) - t154) * t210; (t124 * t132 - t125 * t134 + (-t128 * t160 + t129 * t158) * pkin(2)) * MDP(7) + (t132 * t154 - t128) * t202 + (-t132 * t204 + t126) * MDP(9) + (t213 * t154 * t193 + t219) * MDP(10) + (t128 * (-pkin(3) - t212) - t117 * t132 + t148 * t219 + t213 * t170) * MDP(11) + ((t132 * t163 + t136 * t189 + t159 * t166) * t142 + t166 * t207 + t174) * MDP(17) + (-t132 * t198 + t213 * t169 * t163 + (t136 * t197 + (-t169 * t148 - t209) * t161) * qJD(5) + t196) * MDP(18) + t178 + t217 * (-qJD(2) + t154) * t211; t138 * MDP(18) + (-MDP(18) * t198 + (t142 - t203) * MDP(17) * t163) * t157 * qJD(5); (-t170 * t154 + t128) * MDP(11) - t193 * MDP(10) * t151 + (MDP(17) * t161 + MDP(18) * t163) * (-t142 ^ 2 - t152 * t151); (t171 * t142 - t163 * t184 + t108) * MDP(17) + ((-t187 * t109 - t159 * t123) * t163 + (t187 * t112 - t128 + t184) * t161) * MDP(18) - (t161 * MDP(14) + t163 * MDP(15)) * t187 * t204 - t218 * t151;];
tauc = t1;
