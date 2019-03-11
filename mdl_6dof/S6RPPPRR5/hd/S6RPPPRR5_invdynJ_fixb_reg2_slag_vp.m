% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:54
% EndTime: 2019-03-09 01:37:56
% DurationCPUTime: 2.17s
% Computational Cost: add. (3604->378), mult. (5501->474), div. (0->0), fcn. (3050->8), ass. (0->199)
t113 = sin(pkin(9));
t114 = cos(pkin(9));
t108 = (qJD(1) * qJD(2));
t109 = (qJ(2) * qJDD(1));
t174 = qJDD(3) + t108 + t109;
t68 = pkin(3) * qJDD(1) + t174;
t116 = -pkin(1) - qJ(3);
t246 = (qJDD(1) * t116);
t152 = qJDD(2) + t246;
t197 = (qJD(3) * qJD(1));
t69 = t152 - t197;
t31 = t113 * t68 + t114 * t69;
t29 = qJDD(1) * pkin(7) + t31;
t251 = -qJD(4) * qJD(5) - t29;
t118 = sin(qJ(5));
t206 = qJD(1) * t118;
t250 = qJD(6) * t206 - qJDD(5);
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t193 = t118 * qJDD(1);
t121 = cos(qJ(5));
t205 = qJD(1) * t121;
t23 = t117 * ((qJD(6) + t205) * qJD(5) + t193) + t250 * t120;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t63 = -t119 * t113 + t114 * t122;
t64 = t113 * t122 + t114 * t119;
t249 = t63 * pkin(4) + pkin(7) * t64;
t201 = qJD(6) * t117;
t198 = qJD(1) * qJD(5);
t177 = t118 * t198;
t96 = t121 * qJDD(1);
t62 = qJDD(6) - t96 + t177;
t84 = -qJD(6) + t205;
t145 = t120 * t62 + t201 * t84;
t248 = t145 * t118;
t242 = g(1) * t64;
t163 = -g(2) * t63 + t242;
t236 = g(3) * t121;
t91 = qJ(2) * qJD(1) + qJD(3);
t81 = pkin(3) * qJD(1) + t91;
t82 = qJD(1) * t116 + qJD(2);
t39 = t113 * t81 + t114 * t82;
t37 = qJD(1) * pkin(7) + t39;
t223 = t121 * t37;
t27 = qJD(4) * t118 + t223;
t95 = t121 * qJDD(4);
t10 = -qJD(5) * t27 - t118 * t29 + t95;
t6 = -qJDD(5) * pkin(5) - t10;
t129 = t118 * t163 - t236 - t6;
t247 = pkin(8) * qJD(6) * t84 + t129;
t225 = t118 * t37;
t26 = qJD(4) * t121 - t225;
t154 = t118 * t26 - t121 * t27;
t227 = t113 * t69 - t114 * t68;
t28 = -qJDD(1) * pkin(4) + t227;
t245 = qJD(1) * t154 + t28;
t101 = 2 * t108;
t241 = g(2) * t64;
t19 = qJD(5) * pkin(8) + t27;
t238 = pkin(5) * t121;
t153 = -pkin(8) * t118 - pkin(4) - t238;
t38 = -t113 * t82 + t114 * t81;
t21 = qJD(1) * t153 - t38;
t7 = -t117 * t19 + t120 * t21;
t240 = t7 * t84;
t8 = t117 * t21 + t120 * t19;
t239 = t8 * t84;
t237 = g(3) * t118;
t199 = t120 * qJD(5);
t65 = t117 * t206 - t199;
t235 = t65 * t84;
t204 = qJD(5) * t117;
t67 = t120 * t206 + t204;
t234 = t67 * t65;
t233 = t67 * t84;
t179 = t121 * t199;
t221 = t23 * t120;
t232 = -t118 * t221 - t65 * t179;
t180 = t118 * t199;
t212 = t117 * t121;
t52 = -t113 * t212 - t114 * t120;
t211 = t120 * t121;
t55 = t113 * t117 + t114 * t211;
t231 = t55 * qJD(1) + qJD(6) * t52 - t113 * t180;
t144 = -t113 * t120 + t114 * t212;
t54 = t113 * t211 - t114 * t117;
t230 = -t54 * qJD(1) - qJD(6) * t144 - t114 * t180;
t203 = qJD(5) * t118;
t181 = t117 * t203;
t229 = -t144 * qJD(1) - qJD(6) * t54 + t113 * t181;
t228 = -t52 * qJD(1) - qJD(6) * t55 + t114 * t181;
t115 = pkin(3) + qJ(2);
t61 = t113 * t115 + t114 * t116;
t226 = t117 * t63;
t224 = t120 * t63;
t222 = t121 * t63;
t220 = t122 * pkin(1) + t119 * qJ(2);
t36 = -qJD(1) * pkin(4) - t38;
t219 = qJD(1) * t36;
t77 = qJD(2) * t113 - qJD(3) * t114;
t218 = qJD(1) * t77;
t217 = qJD(5) * t65;
t216 = qJD(6) * t65;
t215 = qJDD(1) * pkin(1);
t124 = qJD(1) ^ 2;
t214 = t113 * t124;
t213 = t114 * t124;
t210 = g(1) * t119 - g(2) * t122;
t111 = t118 ^ 2;
t112 = t121 ^ 2;
t209 = t111 - t112;
t208 = t111 + t112;
t123 = qJD(5) ^ 2;
t207 = t123 + t124;
t202 = qJD(5) * t121;
t200 = qJD(6) * t120;
t195 = qJDD(5) * t118;
t194 = qJDD(5) * t121;
t192 = -t118 * qJDD(4) + t251 * t121;
t190 = t122 * qJ(3) + t220;
t189 = t84 * t204;
t188 = t84 * t199;
t18 = -qJD(5) * pkin(5) - t26;
t187 = t18 * t202;
t186 = t67 * t202;
t185 = t118 * t124 * t121;
t183 = t113 * t206;
t182 = t114 * t206;
t178 = t118 * t200;
t176 = t121 * t198;
t173 = -qJDD(2) + t210;
t100 = t122 * qJ(2);
t172 = -pkin(1) * t119 + t100;
t22 = -qJD(6) * t199 + (-t176 - t193) * t120 + t250 * t117;
t170 = -t22 + t216;
t60 = -t113 * t116 + t114 * t115;
t169 = t119 * pkin(3) + t190;
t168 = qJDD(1) * t208;
t167 = qJDD(2) - t215;
t165 = t67 * t178;
t164 = t118 * t176;
t162 = -g(1) * t63 - t241;
t160 = pkin(5) * t118 - pkin(8) * t121;
t143 = t160 * qJD(5);
t59 = pkin(7) + t61;
t76 = qJD(2) * t114 + qJD(3) * t113;
t161 = -qJD(6) * t121 * t59 + t143 - t76;
t159 = g(1) * t122 + g(2) * t119;
t158 = -t117 * t8 - t120 * t7;
t157 = t117 * t7 - t120 * t8;
t156 = -qJD(6) * t19 + t241;
t155 = t116 * t119 + t100;
t104 = t122 * pkin(3);
t150 = t104 + t155;
t149 = t118 * t6 + t187;
t148 = -t118 * t22 + t186;
t147 = t118 * t23 + t202 * t65;
t146 = -t117 * t62 + t200 * t84;
t9 = -t203 * t37 - t192;
t142 = t101 + (2 * t109) - t159;
t141 = 0.2e1 * t176 + t193;
t140 = 0.2e1 * t177 - t96;
t139 = t64 * pkin(4) - pkin(7) * t63 + t169;
t138 = t163 - t219;
t137 = qJD(1) * t76 + t162;
t136 = -t121 * t207 - t195;
t135 = t118 * t207 - t194;
t134 = -pkin(8) * t62 - t18 * t84;
t40 = t153 - t60;
t133 = -qJD(6) * t40 - t121 * t77 + t203 * t59;
t132 = -t121 * t163 - t237;
t5 = qJDD(5) * pkin(8) + t9;
t131 = -g(2) * t222 - qJD(6) * t21 + t237 - t5;
t58 = -pkin(4) - t60;
t130 = -qJDD(5) * t59 + (qJD(1) * t58 + t36 - t77) * qJD(5);
t12 = qJD(1) * t143 + qJDD(1) * t153 + t227;
t1 = qJD(6) * t7 + t117 * t12 + t120 * t5;
t11 = t120 * t12;
t2 = -qJD(6) * t8 - t117 * t5 + t11;
t128 = qJD(6) * t158 + t1 * t120 - t2 * t117;
t127 = -t10 * t118 + t9 * t121 + (-t118 * t27 - t121 * t26) * qJD(5);
t126 = qJDD(1) * t58 + t123 * t59 - t137 + t28;
t125 = t127 + t219;
t79 = -t118 * t123 + t194;
t78 = t121 * t123 + t195;
t75 = qJDD(1) * t114 - t214;
t74 = -qJDD(1) * t113 - t213;
t73 = t160 * qJD(1);
t47 = t67 * t203;
t25 = t211 * t64 - t226;
t24 = -t212 * t64 - t224;
t16 = t117 * t40 + t211 * t59;
t15 = t120 * t40 - t212 * t59;
t14 = t117 * t73 + t120 * t26;
t13 = -t117 * t26 + t120 * t73;
t4 = t117 * t133 + t120 * t161;
t3 = t117 * t161 - t120 * t133;
t17 = [0, 0, 0, 0, 0, qJDD(1), t210, t159, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t173 - 0.2e1 * t215, t142, -t167 * pkin(1) - g(1) * t172 - g(2) * t220 + (t101 + t109) * qJ(2), 0, 0, 0, qJDD(1), 0, 0, qJDD(3) + t142, 0, t173 + (2 * t197) - (2 * t246) -g(1) * t155 - g(2) * t190 + qJ(2) * t174 + t91 * qJD(2) - t82 * qJD(3) + t69 * t116, 0, 0, 0, 0, 0, qJDD(1), qJDD(1) * t60 + t137 - t227, -qJDD(1) * t61 + t163 - t218 - t31, 0, -g(1) * t150 - g(2) * t169 - t227 * t60 + t31 * t61 + t38 * t76 + t39 * t77, qJDD(1) * t111 + 0.2e1 * t164, 0.2e1 * t118 * t96 - 0.2e1 * t198 * t209, t78, qJDD(1) * t112 - 0.2e1 * t164, t79, 0, t118 * t130 - t121 * t126, t118 * t126 + t121 * t130, t168 * t59 + t208 * t218 + t127 - t163, t28 * t58 - t36 * t76 - g(1) * (t150 + t249) - g(2) * t139 - t154 * t77 + t127 * t59, t67 * t179 + (-t120 * t22 - t201 * t67) * t118, -t165 + (-t186 + (t22 + t216) * t118) * t117 + t232, t47 + (t22 - t188) * t121 + t248, t117 * t147 + t178 * t65 (t23 + t189) * t121 + (t146 - t217) * t118, -t121 * t62 - t203 * t84, -t117 * t242 - g(2) * t25 + t15 * t62 - t4 * t84 + (-g(1) * t224 - t2 + (t117 * t18 + t59 * t65) * qJD(5)) * t121 + (qJD(5) * t7 + t117 * t6 + t18 * t200 + t23 * t59 + t65 * t77) * t118, -t120 * t242 - g(2) * t24 - t16 * t62 + t3 * t84 + (g(1) * t226 + t1 + (t120 * t18 + t59 * t67) * qJD(5)) * t121 + (-qJD(5) * t8 + t120 * t6 - t18 * t201 - t22 * t59 + t67 * t77) * t118, t15 * t22 - t16 * t23 - t3 * t65 - t4 * t67 + t158 * t202 + (qJD(6) * t157 - t1 * t117 - t120 * t2 + t162) * t118, t1 * t16 + t8 * t3 + t2 * t15 + t7 * t4 + t59 * t187 - g(1) * (pkin(5) * t222 - qJ(3) * t119 + t104 + t172 + t249) - g(2) * (t238 * t64 + t139) + (pkin(8) * t162 + t18 * t77 + t6 * t59) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t124, -qJ(2) * t124 + t167 - t210, 0, 0, 0, 0, 0, 0, -t124, 0, -qJDD(1) (-qJD(3) - t91) * qJD(1) + t152 - t210, 0, 0, 0, 0, 0, 0, t74, -t75, 0, t113 * t227 + t114 * t31 + (-t113 * t39 - t114 * t38) * qJD(1) - t210, 0, 0, 0, 0, 0, 0, t113 * t140 + t114 * t136, t113 * t141 + t114 * t135, t114 * t168 - t208 * t214, t113 * t245 + t125 * t114 - t210, 0, 0, 0, 0, 0, 0, t114 * t147 - t144 * t62 - t183 * t65 - t228 * t84, t114 * t148 - t183 * t67 + t230 * t84 - t55 * t62, -t144 * t22 - t228 * t67 - t23 * t55 - t230 * t65, t1 * t55 + t114 * t149 - t144 * t2 - t18 * t183 + t228 * t7 + t230 * t8 - t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t124, qJD(1) * t82 - t159 + t174, 0, 0, 0, 0, 0, 0, t75, t74, 0, t113 * t31 - t114 * t227 + (-t113 * t38 + t114 * t39) * qJD(1) - t159, 0, 0, 0, 0, 0, 0, t113 * t136 - t114 * t140, t113 * t135 - t114 * t141, t113 * t168 + t208 * t213, t125 * t113 - t114 * t245 - t159, 0, 0, 0, 0, 0, 0, t113 * t147 + t182 * t65 - t229 * t84 + t52 * t62, t113 * t148 + t182 * t67 + t231 * t84 - t54 * t62, t22 * t52 - t229 * t67 - t23 * t54 - t231 * t65, t1 * t54 + t113 * t149 + t18 * t182 + t2 * t52 + t229 * t7 + t231 * t8 - t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) - g(3), 0, 0, 0, 0, 0, 0, t79, -t78, 0, -qJD(5) * t154 + t10 * t121 + t118 * t9 - g(3), 0, 0, 0, 0, 0, 0 (-t23 + t189) * t121 + (t146 + t217) * t118, t47 + (t22 + t188) * t121 - t248, t165 + (t118 * t170 + t186) * t117 + t232, -g(3) + (-qJD(5) * t157 - t6) * t121 + (qJD(5) * t18 + t128) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t209 * t124, t193, t185, t96, qJDD(5), -t236 + t95 + (t27 - t223) * qJD(5) + (t138 + t251) * t118, t237 + (t26 + t225) * qJD(5) + t138 * t121 + t192, 0, 0, -t22 * t117 - t120 * t233 (-t22 + t235) * t120 + (-t23 + t233) * t117 (-t118 * t67 + t211 * t84) * qJD(1) - t146, -t117 * t235 - t221 (t118 * t65 - t212 * t84) * qJD(1) + t145, t84 * t206, -pkin(5) * t23 + t134 * t117 + t247 * t120 + t13 * t84 - t7 * t206 - t27 * t65, pkin(5) * t22 - t247 * t117 + t134 * t120 - t14 * t84 + t8 * t206 - t27 * t67, t13 * t67 + t14 * t65 + (t1 + t240 + (qJD(6) * t67 - t23) * pkin(8)) * t120 + (pkin(8) * t170 - t2 + t239) * t117 + t132, -t7 * t13 - t8 * t14 - t18 * t27 + t129 * pkin(5) + (t128 + t132) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, -t65 ^ 2 + t67 ^ 2, -t22 - t235, -t234, -t23 - t233, t62, -g(1) * t24 + t117 * t131 + t120 * t156 - t18 * t67 + t11 - t239, g(1) * t25 + t18 * t65 - t240 + (-t12 - t156) * t117 + t131 * t120, 0, 0;];
tau_reg  = t17;
