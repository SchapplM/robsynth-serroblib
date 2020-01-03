% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:32
% DurationCPUTime: 1.70s
% Computational Cost: add. (4074->262), mult. (8855->291), div. (0->0), fcn. (5184->6), ass. (0->169)
t217 = pkin(2) + pkin(3);
t127 = qJDD(2) - qJDD(4);
t135 = sin(qJ(4));
t136 = sin(qJ(2));
t138 = cos(qJ(4));
t139 = cos(qJ(2));
t90 = (-t136 * t135 - t139 * t138) * qJD(1);
t175 = qJD(1) * t139;
t176 = qJD(1) * t136;
t92 = -t135 * t175 + t138 * t176;
t205 = t92 * t90;
t229 = -t127 + t205;
t236 = t229 * pkin(4);
t235 = t135 * t229;
t234 = t138 * t229;
t140 = qJD(2) ^ 2;
t132 = t136 ^ 2;
t141 = qJD(1) ^ 2;
t187 = t132 * t141;
t112 = t140 + t187;
t179 = t139 * t141;
t116 = t136 * t179;
t110 = qJDD(2) - t116;
t180 = t139 * t110;
t233 = pkin(6) * (-t136 * t112 + t180);
t128 = qJD(2) - qJD(4);
t196 = t90 * t128;
t173 = qJD(1) * qJD(2);
t170 = t139 * t173;
t172 = t136 * qJDD(1);
t103 = t170 + t172;
t123 = t139 * qJDD(1);
t171 = t136 * t173;
t104 = t123 - t171;
t59 = t90 * qJD(4) + t138 * t103 - t135 * t104;
t48 = t59 + t196;
t232 = qJ(5) * t48;
t231 = t59 - t196;
t230 = t104 - t171;
t111 = -qJD(2) * pkin(3) - pkin(7) * t176;
t133 = t139 ^ 2;
t186 = t133 * t141;
t185 = t136 * qJ(3);
t209 = t139 * pkin(2);
t158 = -t185 - t209;
t100 = t158 * qJD(1);
t218 = 2 * qJD(3);
t137 = sin(qJ(1));
t211 = cos(qJ(1));
t150 = t211 * g(1) + t137 * g(2);
t190 = qJDD(1) * pkin(6);
t94 = -t141 * pkin(1) - t150 + t190;
t82 = t139 * t94;
t164 = qJDD(2) * qJ(3) + qJD(2) * t218 + t100 * t175 + t82;
t210 = t136 * g(3);
t148 = t164 - t210;
t207 = t140 * pkin(2);
t52 = t148 - t207;
t34 = -pkin(3) * t186 - t104 * pkin(7) + qJD(2) * t111 + t52;
t169 = -qJDD(2) * pkin(2) - t140 * qJ(3) + qJDD(3);
t208 = t139 * g(3);
t151 = t169 + t208;
t165 = qJD(1) * t100 + t94;
t36 = -qJDD(2) * pkin(3) + (-t103 + t170) * pkin(7) + (-pkin(3) * t179 + t165) * t136 + t151;
t16 = t135 * t36 + t138 * t34;
t166 = t135 * t103 + t138 * t104;
t58 = -t92 * qJD(4) - t166;
t76 = -t128 * pkin(4) - t92 * qJ(5);
t149 = t58 * qJ(5) + 0.2e1 * qJD(5) * t90 + t128 * t76 + t16;
t126 = t128 ^ 2;
t89 = t92 ^ 2;
t73 = -t89 - t126;
t88 = t90 ^ 2;
t228 = -t149 + (t73 + t88) * pkin(4);
t226 = t104 * pkin(3) - pkin(7) * t186;
t225 = t165 * t136;
t224 = t180 + t136 * (-t140 + t186);
t61 = -t126 - t88;
t37 = t135 * t61 + t234;
t38 = t138 * t61 - t235;
t222 = qJ(3) * t38 - t217 * t37;
t62 = t205 + t127;
t199 = t135 * t62;
t49 = t138 * t73 + t199;
t197 = t138 * t62;
t50 = -t135 * t73 + t197;
t221 = qJ(3) * t50 - t217 * t49;
t44 = (qJD(4) + t128) * t92 + t166;
t168 = t137 * g(1) - t211 * g(2);
t93 = qJDD(1) * pkin(1) + t141 * pkin(6) + t168;
t143 = t230 * pkin(2) + t93;
t142 = t176 * t218 + t143;
t157 = t103 + t170;
t220 = t157 * qJ(3) + t142;
t219 = -t58 * pkin(4) - t88 * qJ(5) + t92 * t76 + qJDD(5);
t15 = t135 * t34 - t138 * t36;
t146 = -t15 - t232 + t236;
t191 = qJD(5) * t92;
t84 = -0.2e1 * t191;
t10 = t146 + t84;
t214 = pkin(4) * t10;
t213 = pkin(4) * t48;
t23 = -t135 * t44 - t138 * t48;
t24 = t135 * t48 - t138 * t44;
t60 = -t88 - t89;
t212 = pkin(1) * t60 + pkin(6) * (t136 * t23 + t139 * t24);
t43 = (qJD(4) - t128) * t92 + t166;
t204 = pkin(6) * (t136 * t37 + t139 * t38) + pkin(1) * t43;
t203 = pkin(6) * (t136 * t49 + t139 * t50) + pkin(1) * t231;
t105 = t123 - 0.2e1 * t171;
t114 = -t140 - t186;
t109 = qJDD(2) + t116;
t183 = t136 * t109;
t202 = pkin(6) * (t139 * t114 - t183) + pkin(1) * t105;
t195 = qJ(3) * t139;
t29 = t103 * qJ(3) + (qJD(2) * t195 + (t218 + t111) * t136) * qJD(1) + t143 + t226;
t200 = t135 * t29;
t198 = t138 * t29;
t194 = qJ(5) * t135;
t193 = qJ(5) * t138;
t189 = t128 * t135;
t188 = t128 * t138;
t184 = t136 * t105;
t177 = t132 + t133;
t107 = t177 * t141;
t178 = pkin(1) * t107 + t177 * t190;
t74 = t136 * t94 + t208;
t75 = t82 - t210;
t167 = t136 * t74 + t139 * t75;
t163 = qJ(3) * t24 - t217 * t23;
t162 = -pkin(7) * t23 + qJ(3) * t60;
t161 = -pkin(7) * t37 + qJ(3) * t43;
t160 = -pkin(7) * t49 + qJ(3) * t231;
t6 = t135 * t16 - t138 * t15;
t7 = t135 * t15 + t138 * t16;
t102 = 0.2e1 * t170 + t172;
t155 = t139 * t102 + t184;
t154 = -pkin(7) * t38 + t217 * t43;
t153 = -pkin(7) * t50 + t217 * t231;
t152 = -pkin(7) * t24 + t217 * t60;
t147 = t139 * t217 + pkin(1) + t185;
t145 = t169 + t225;
t144 = t146 + t236;
t14 = t29 + t219;
t108 = (t132 - t133) * t141;
t85 = 0.2e1 * t191;
t78 = -t89 + t126;
t77 = t88 - t126;
t72 = t183 + t139 * (t140 - t187);
t71 = t157 * t136;
t70 = t230 * t139;
t65 = t89 - t88;
t57 = t145 + t208;
t28 = (t136 * (-t135 * t92 - t138 * t90) + t139 * (t135 * t90 - t138 * t92)) * t128;
t27 = -pkin(4) * t231 + qJ(5) * t62;
t26 = t136 * (t138 * t77 + t199) + t139 * (-t135 * t77 + t197);
t25 = t136 * (-t135 * t78 + t234) + t139 * (-t138 * t78 - t235);
t19 = t136 * (t138 * t59 + t92 * t189) + t139 * (-t135 * t59 + t92 * t188);
t18 = t136 * (-t135 * t58 + t90 * t188) + t139 * (-t138 * t58 - t90 * t189);
t13 = -qJ(5) * t73 + t14;
t12 = -t88 * pkin(4) + t149;
t11 = t136 * (-t135 * t231 - t138 * t43) + t139 * (t135 * t43 - t138 * t231);
t8 = -pkin(4) * t43 + qJ(5) * t61 - t111 * t176 - t219 - t220 - t226;
t5 = -t146 + t85 + t232;
t4 = -qJ(5) * t44 + (-t60 - t88) * pkin(4) + t149;
t3 = -pkin(4) * t14 + qJ(5) * t12;
t2 = -t135 * t10 + t138 * t12;
t1 = t138 * t10 + t135 * t12;
t9 = [0, 0, 0, 0, 0, qJDD(1), t168, t150, 0, 0, t71, t155, t72, t70, t224, 0, t139 * t93 + t202, -pkin(1) * t102 - t136 * t93 - t233, t167 + t178, pkin(1) * t93 + pkin(6) * t167, t71, t72, -t155, 0, -t224, t70, t139 * (pkin(2) * t105 + t142) + (t139 * t157 + t184) * qJ(3) + t202, t139 * (pkin(2) * t107 + t164 - t207) + (qJ(3) * t107 + t145) * t136 + t178, t136 * t142 + t233 + (pkin(1) + t209) * t102 + (t102 + t157) * t185, pkin(6) * (t136 * t57 + t139 * t52) + (pkin(1) - t158) * t220, t19, t11, t25, t18, t26, t28, t136 * (t161 + t200) + t139 * (t154 + t198) + t204, t136 * (t160 + t198) + t139 * (t153 - t200) + t203, t136 * (t162 - t6) + t139 * (t152 - t7) + t212, t147 * t29 + (pkin(6) - pkin(7)) * (t136 * t6 + t139 * t7), t19, t11, t25, t18, t26, t28, t136 * (-t135 * t8 - t193 * t229 + t161) + t139 * (-t138 * t8 + t194 * t229 + t154) + t204, t136 * (t138 * t13 - t135 * t27 + t160) + t139 * (-t135 * t13 - t138 * t27 + t153) + t203, t136 * (-t135 * t4 + t138 * t5 + t162) + t139 * (-t135 * t5 - t138 * t4 + t152) + t212, t136 * (-pkin(7) * t1 - t10 * t193 - t135 * t3) + t139 * (-pkin(7) * t2 + t10 * t194 - t138 * t3) + pkin(6) * (t136 * t1 + t139 * t2) + t147 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t108, t172, t116, t123, qJDD(2), -t74, -t75, 0, 0, -t116, t172, -t108, qJDD(2), -t123, t116, pkin(2) * t109 + qJ(3) * t114 - t151 - t225, (-pkin(2) * t136 + t195) * qJDD(1), qJ(3) * t110 + (t112 - t140) * pkin(2) + t148, -pkin(2) * t57 + qJ(3) * t52, t205, -t65, -t48, -t205, t44, t127, t15 + t222, t16 + t221, t163, qJ(3) * t7 - t217 * t6, t205, -t65, -t48, -t205, t44, t127, -t144 + t85 + t222, t221 - t228, t163 + t213, qJ(3) * t2 - t217 * t1 - t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t172, -t112, t57, 0, 0, 0, 0, 0, 0, t37, t49, t23, t6, 0, 0, 0, 0, 0, 0, t37, t49, t23, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t65, t48, t205, -t44, -t127, -t15, -t16, 0, 0, -t205, t65, t48, t205, -t44, -t127, t144 + t84, t228, -t213, t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t231, t60, t14;];
tauJ_reg = t9;
