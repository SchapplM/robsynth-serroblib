% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:42:12
% EndTime: 2019-05-04 21:42:19
% DurationCPUTime: 3.29s
% Computational Cost: add. (11901->337), mult. (24670->516), div. (0->0), fcn. (18952->14), ass. (0->201)
t236 = 2 * qJD(4);
t170 = sin(qJ(6));
t163 = sin(pkin(12));
t166 = cos(pkin(12));
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t191 = t163 * t174 + t166 * t171;
t141 = t191 * qJD(2);
t173 = cos(qJ(6));
t126 = -t173 * qJD(5) + t141 * t170;
t128 = t170 * qJD(5) + t141 * t173;
t101 = t128 * t126;
t200 = t166 * qJDD(2);
t201 = t163 * qJDD(2);
t192 = t171 * t201 - t174 * t200;
t205 = t141 * qJD(5);
t116 = -t192 - t205;
t109 = qJDD(6) - t116;
t229 = -t101 + t109;
t235 = t170 * t229;
t212 = t163 * t171;
t139 = (-t166 * t174 + t212) * qJD(2);
t119 = t141 * t139;
t227 = qJDD(5) - t119;
t234 = t171 * t227;
t233 = t173 * t229;
t232 = t174 * t227;
t159 = t166 ^ 2;
t177 = t163 ^ 2;
t226 = qJD(2) ^ 2;
t147 = (t159 + t177) * t226;
t164 = sin(pkin(11));
t172 = sin(qJ(2));
t175 = cos(qJ(2));
t165 = sin(pkin(6));
t168 = cos(pkin(6));
t215 = sin(pkin(10));
t216 = cos(pkin(10));
t188 = t215 * g(1) - t216 * g(2);
t187 = t168 * t188;
t208 = -g(3) + qJDD(1);
t184 = t165 * t208 + t187;
t189 = -t216 * g(1) - t215 * g(2);
t107 = -t172 * t189 + t175 * t184;
t183 = qJDD(2) * pkin(2) + t107;
t108 = t172 * t184 + t175 * t189;
t106 = -t226 * pkin(2) + t108;
t167 = cos(pkin(11));
t211 = t167 * t106;
t82 = t164 * t183 + t211;
t182 = qJDD(2) * qJ(4) + t82;
t231 = -t226 * pkin(3) + qJD(2) * t236 + t182;
t230 = qJ(4) + pkin(8);
t135 = qJD(6) + t139;
t138 = t191 * qJDD(2);
t206 = t139 * qJD(5);
t118 = t138 - t206;
t193 = -t173 * qJDD(5) + t118 * t170;
t73 = (qJD(6) - t135) * t128 + t193;
t124 = t126 ^ 2;
t125 = t128 ^ 2;
t134 = t135 ^ 2;
t136 = t139 ^ 2;
t137 = t141 ^ 2;
t225 = t166 * pkin(4);
t151 = t168 * t208;
t199 = t151 + qJDD(3);
t186 = -t165 * t188 + t199;
t185 = t166 * t186;
t180 = -pkin(8) * t201 + t185 + (t225 * t226 - t231) * t163;
t207 = t226 * t159;
t67 = t163 * t186 + t231 * t166;
t62 = -pkin(4) * t207 + pkin(8) * t200 + t67;
t36 = t171 * t62 - t174 * t180;
t218 = t174 * t62;
t37 = t218 + t171 * t185 + (-t211 - t164 * t107 + (-t164 * pkin(2) - t230) * qJDD(2) + (-(2 * qJD(4)) + (pkin(3) + t225) * qJD(2)) * qJD(2)) * t212;
t19 = t171 * t37 - t174 * t36;
t224 = t163 * t19;
t110 = pkin(5) * t139 - pkin(9) * t141;
t176 = qJD(5) ^ 2;
t31 = -qJDD(5) * pkin(5) - t176 * pkin(9) + t110 * t141 + t36;
t223 = t170 * t31;
t85 = t101 + t109;
t222 = t170 * t85;
t162 = qJDD(2) * pkin(3);
t194 = t106 * t164 - t167 * t183;
t80 = -t226 * qJ(4) + qJDD(4) - t162 + t194;
t72 = -pkin(4) * t200 + t80 + (-t177 * t226 - t207) * pkin(8);
t221 = t171 * t72;
t220 = t173 * t31;
t219 = t173 * t85;
t217 = t174 * t72;
t214 = t135 * t170;
t213 = t135 * t173;
t113 = qJDD(5) + t119;
t210 = t171 * t113;
t209 = t174 * t113;
t203 = qJD(6) + t135;
t198 = t171 * t101;
t197 = t174 * t101;
t196 = -pkin(5) * t174 - pkin(4);
t66 = -t185 + ((-pkin(3) * qJD(2) + t236) * qJD(2) + t182) * t163;
t41 = t163 * t66 + t166 * t67;
t32 = -t176 * pkin(5) + qJDD(5) * pkin(9) - t139 * t110 + t171 * t180 + t218;
t48 = (-t118 + t206) * pkin(9) + (-t116 + t205) * pkin(5) + t72;
t22 = t170 * t32 - t173 * t48;
t23 = t170 * t48 + t173 * t32;
t11 = t170 * t22 + t173 * t23;
t20 = t171 * t36 + t174 * t37;
t195 = -t80 + t162;
t4 = t171 * t11 - t174 * t31;
t5 = t174 * t11 + t171 * t31;
t3 = -t163 * t4 + t166 * t5;
t10 = t170 * t23 - t173 * t22;
t190 = -t170 * qJDD(5) - t118 * t173;
t93 = -t126 * qJD(6) - t190;
t156 = t159 * qJDD(2);
t155 = t177 * qJDD(2);
t146 = -qJDD(2) * t164 - t226 * t167;
t145 = qJDD(2) * t167 - t226 * t164;
t144 = t156 + t155;
t143 = t166 * t147;
t142 = t163 * t147;
t131 = -t137 - t176;
t130 = -t137 + t176;
t129 = t136 - t176;
t122 = -t143 * t164 + t167 * t200;
t121 = t142 * t164 - t167 * t201;
t120 = t144 * t164 + t147 * t167;
t117 = t138 - 0.2e1 * t206;
t115 = t192 + 0.2e1 * t205;
t111 = -t176 - t136;
t105 = t135 * t126;
t104 = -t125 + t134;
t103 = t124 - t134;
t100 = -t136 - t137;
t99 = t125 - t124;
t97 = -t125 - t134;
t96 = -t171 * t131 - t209;
t95 = t174 * t131 - t210;
t94 = -t134 - t124;
t92 = -qJD(6) * t128 - t193;
t91 = t124 + t125;
t90 = t171 * t138 - t174 * t192;
t89 = -t174 * t138 - t171 * t192;
t88 = t174 * t111 - t234;
t87 = t171 * t111 + t232;
t83 = (-t126 * t173 + t128 * t170) * t135;
t78 = t203 * t126 + t190;
t77 = t105 + t93;
t76 = -t105 + t93;
t75 = -t203 * t128 - t193;
t71 = -t128 * t214 + t173 * t93;
t70 = t126 * t213 - t170 * t92;
t69 = -t163 * t95 + t166 * t96;
t68 = t163 * t96 + t166 * t95;
t65 = t103 * t173 - t222;
t64 = -t170 * t104 + t233;
t60 = -t170 * t97 - t219;
t59 = t173 * t97 - t222;
t58 = -t163 * t89 + t166 * t90;
t57 = t163 * t90 + t166 * t89;
t56 = t173 * t94 - t235;
t55 = t170 * t94 + t233;
t54 = -t163 * t87 + t166 * t88;
t53 = t163 * t88 + t166 * t87;
t52 = -t117 * t167 + t164 * t69;
t51 = t164 * t82 - t167 * t194;
t50 = -t115 * t167 + t164 * t54;
t49 = -t100 * t167 + t164 * t58;
t46 = -t170 * t76 + t173 * t75;
t45 = t170 * t77 - t173 * t73;
t44 = -t170 * t73 - t173 * t77;
t43 = -t171 * t78 + t174 * t60;
t42 = t171 * t60 + t174 * t78;
t40 = t163 * t67 - t166 * t66;
t39 = -t171 * t75 + t174 * t56;
t38 = t171 * t56 + t174 * t75;
t34 = -t171 * t91 + t174 * t45;
t33 = t171 * t45 + t174 * t91;
t30 = t164 * t41 - t167 * t80;
t29 = -pkin(9) * t59 + t220;
t28 = -pkin(9) * t55 + t223;
t27 = -t163 * t42 + t166 * t43;
t26 = t163 * t43 + t166 * t42;
t25 = -t163 * t38 + t166 * t39;
t24 = t163 * t39 + t166 * t38;
t18 = -t163 * t33 + t166 * t34;
t17 = t163 * t34 + t166 * t33;
t16 = t164 * t27 - t167 * t59;
t15 = t164 * t25 - t167 * t55;
t14 = -pkin(5) * t59 + t23;
t13 = -pkin(5) * t55 + t22;
t12 = t164 * t18 - t167 * t44;
t9 = t166 * t20 - t224;
t8 = t163 * t20 + t166 * t19;
t7 = t164 * t9 - t167 * t72;
t6 = -pkin(9) * t44 - t10;
t2 = t163 * t5 + t166 * t4;
t1 = -t10 * t167 + t164 * t3;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0, 0, 0, 0, 0, 0, (qJDD(2) * t175 - t226 * t172) * t165, (-qJDD(2) * t172 - t226 * t175) * t165, 0, t168 * t151 + (t175 * t107 + t172 * t108 - t187) * t165, 0, 0, 0, 0, 0, 0, (t145 * t175 + t146 * t172) * t165, (-t145 * t172 + t146 * t175) * t165, 0, t168 * t199 + (t172 * (t164 * t194 + t167 * t82) + t175 * t51 - t187) * t165, 0, 0, 0, 0, 0, 0, (t172 * (-t143 * t167 - t164 * t200) + t175 * t122) * t165, (t172 * (t142 * t167 + t164 * t201) + t175 * t121) * t165, (t172 * (t144 * t167 - t147 * t164) + t175 * t120) * t165, t168 * t40 + (t172 * (t164 * t80 + t167 * t41) + t175 * t30) * t165, 0, 0, 0, 0, 0, 0, t168 * t53 + (t172 * (t115 * t164 + t167 * t54) + t175 * t50) * t165, t168 * t68 + (t172 * (t117 * t164 + t167 * t69) + t175 * t52) * t165, t168 * t57 + (t172 * (t100 * t164 + t167 * t58) + t175 * t49) * t165, t168 * t8 + (t172 * (t164 * t72 + t167 * t9) + t175 * t7) * t165, 0, 0, 0, 0, 0, 0, t168 * t24 + (t172 * (t164 * t55 + t167 * t25) + t175 * t15) * t165, t168 * t26 + (t172 * (t164 * t59 + t167 * t27) + t175 * t16) * t165, t168 * t17 + (t172 * (t164 * t44 + t167 * t18) + t175 * t12) * t165, t168 * t2 + (t172 * (t10 * t164 + t167 * t3) + t175 * t1) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t107, -t108, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t145 - t194, pkin(2) * t146 - t82, 0, pkin(2) * t51, t155, 0.2e1 * t163 * t200, 0, t156, 0, 0, pkin(2) * t122 - qJ(4) * t143 + t195 * t166, pkin(2) * t121 + qJ(4) * t142 - t195 * t163, pkin(2) * t120 + pkin(3) * t147 + qJ(4) * t144 + t41, pkin(2) * t30 - pkin(3) * t80 + qJ(4) * t41, t163 * (t174 * t118 - t171 * t205) + t166 * (t171 * t118 + t174 * t205), t163 * (-t174 * t115 - t171 * t117) + t166 * (-t171 * t115 + t174 * t117), t163 * (-t171 * t130 + t232) + t166 * (t174 * t130 + t234), t163 * (-t171 * t116 + t174 * t206) + t166 * (t174 * t116 + t171 * t206), t163 * (t174 * t129 - t210) + t166 * (t171 * t129 + t209), (t163 * (-t139 * t174 + t141 * t171) + t166 * (-t139 * t171 - t141 * t174)) * qJD(5), t163 * (-pkin(8) * t87 + t221) + t166 * (-pkin(4) * t115 + pkin(8) * t88 - t217) - pkin(3) * t115 + qJ(4) * t54 + pkin(2) * t50, t163 * (-pkin(8) * t95 + t217) + t166 * (-pkin(4) * t117 + pkin(8) * t96 + t221) - pkin(3) * t117 + qJ(4) * t69 + pkin(2) * t52, t163 * (-pkin(8) * t89 - t19) + t166 * (-pkin(4) * t100 + pkin(8) * t90 + t20) - pkin(3) * t100 + qJ(4) * t58 + pkin(2) * t49, -pkin(8) * t224 + t166 * (-pkin(4) * t72 + pkin(8) * t20) - pkin(3) * t72 + qJ(4) * t9 + pkin(2) * t7, t163 * (t174 * t71 + t198) + t166 * (t171 * t71 - t197), t163 * (t171 * t99 + t174 * t46) + t166 * (t171 * t46 - t174 * t99), t163 * (t171 * t77 + t174 * t64) + t166 * (t171 * t64 - t174 * t77), t163 * (t174 * t70 - t198) + t166 * (t171 * t70 + t197), t163 * (-t171 * t73 + t174 * t65) + t166 * (t171 * t65 + t174 * t73), t163 * (t171 * t109 + t174 * t83) + t166 * (-t174 * t109 + t171 * t83), t163 * (-pkin(8) * t38 - t171 * t13 + t174 * t28) + t166 * (-pkin(4) * t55 + pkin(8) * t39 + t174 * t13 + t171 * t28) - pkin(3) * t55 + qJ(4) * t25 + pkin(2) * t15, t163 * (-pkin(8) * t42 - t171 * t14 + t174 * t29) + t166 * (-pkin(4) * t59 + pkin(8) * t43 + t174 * t14 + t171 * t29) - pkin(3) * t59 + qJ(4) * t27 + pkin(2) * t16, t163 * (-pkin(8) * t33 + t174 * t6) + t166 * (pkin(8) * t34 + t171 * t6) + qJ(4) * t18 + pkin(2) * t12 + (pkin(5) * t212 + t166 * t196 - pkin(3)) * t44, pkin(2) * t1 + (t163 * (pkin(5) * t171 - pkin(9) * t174) + t166 * (-pkin(9) * t171 + t196) - pkin(3)) * t10 + t230 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t53, t68, t57, t8, 0, 0, 0, 0, 0, 0, t24, t26, t17, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t201, -t147, t80, 0, 0, 0, 0, 0, 0, t115, t117, t100, t72, 0, 0, 0, 0, 0, 0, t55, t59, t44, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t137 - t136, t138, -t119, -t192, qJDD(5), -t36, -t37, 0, 0, t128 * t213 + t170 * t93, t170 * t75 + t173 * t76, t173 * t104 + t235, t126 * t214 + t173 * t92, t170 * t103 + t219, (-t126 * t170 - t128 * t173) * t135, pkin(5) * t75 + pkin(9) * t56 - t220, pkin(5) * t78 + pkin(9) * t60 + t223, pkin(5) * t91 + pkin(9) * t45 + t11, -pkin(5) * t31 + pkin(9) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t99, t77, -t101, -t73, t109, -t22, -t23, 0, 0;];
tauJ_reg  = t21;
