% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [6x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:45
% EndTime: 2019-03-08 18:43:50
% DurationCPUTime: 2.55s
% Computational Cost: add. (2345->319), mult. (6243->492), div. (0->0), fcn. (6148->16), ass. (0->182)
t133 = cos(qJ(5));
t198 = qJD(3) * t133;
t106 = -qJD(6) + t198;
t238 = t106 + qJD(6);
t128 = cos(pkin(6));
t122 = sin(pkin(7));
t131 = sin(qJ(3));
t208 = t122 * t131;
t123 = sin(pkin(6));
t120 = sin(pkin(12));
t134 = cos(qJ(3));
t125 = cos(pkin(12));
t127 = cos(pkin(7));
t204 = t125 * t127;
t155 = t120 * t134 + t131 * t204;
t235 = t155 * t123;
t63 = t128 * t208 + t235;
t130 = sin(qJ(5));
t192 = qJD(6) * t130;
t237 = -qJD(3) * t192 + qJDD(5);
t129 = sin(qJ(6));
t132 = cos(qJ(6));
t187 = t130 * qJDD(3);
t59 = (qJD(5) * (qJD(6) + t198) + t187) * t129 - t237 * t132;
t105 = qJD(1) * t128 + qJD(2);
t57 = qJD(1) * t235 + t105 * t208;
t119 = sin(pkin(13));
t124 = cos(pkin(13));
t212 = t120 * t131;
t156 = t134 * t204 - t212;
t207 = t122 * t134;
t62 = t123 * t156 + t128 * t207;
t46 = t119 * t62 + t124 * t63;
t209 = t122 * t123;
t185 = t125 * t209;
t79 = t127 * t128 - t185;
t73 = t79 * t133;
t236 = -t130 * t46 + t73;
t190 = t132 * qJD(5);
t180 = t133 * t190;
t114 = t133 * qJDD(3);
t188 = qJD(3) * qJD(5);
t89 = t130 * t188 + qJDD(6) - t114;
t219 = t132 * t89;
t234 = -t106 * (-t129 * t192 + t180) + t130 * t219;
t160 = t119 * t134 + t124 * t131;
t75 = t160 * t122;
t77 = t160 * t127;
t202 = t134 * t124;
t90 = t119 * t131 - t202;
t144 = t128 * t75 + (-t120 * t90 + t125 * t77) * t123;
t126 = cos(pkin(11));
t206 = t123 * t126;
t121 = sin(pkin(11));
t203 = t126 * t128;
t80 = -t120 * t121 + t125 * t203;
t81 = t120 * t203 + t121 * t125;
t146 = -t206 * t75 + t77 * t80 - t81 * t90;
t211 = t121 * t123;
t210 = t121 * t128;
t82 = -t120 * t126 - t125 * t210;
t83 = -t120 * t210 + t125 * t126;
t147 = t211 * t75 + t77 * t82 - t83 * t90;
t205 = t123 * t127;
t64 = -t122 * t80 - t126 * t205;
t65 = t121 * t205 - t122 * t82;
t154 = g(1) * (-t130 * t147 + t133 * t65) + g(2) * (-t130 * t146 + t133 * t64) + g(3) * (-t130 * t144 + t73);
t168 = pkin(5) * t130 - pkin(10) * t133;
t213 = t120 * t123;
t177 = qJDD(1) * t213;
t102 = t128 * qJDD(1) + qJDD(2);
t176 = qJDD(1) * t204;
t226 = t134 * t123 * t176 + t102 * t207;
t41 = qJDD(3) * pkin(3) - t57 * qJD(3) - t131 * t177 + t226;
t189 = qJD(1) * qJD(3);
t42 = (qJD(3) * t105 * t134 + t102 * t131) * t122 + (qJDD(1) * t155 + t156 * t189) * t123;
t12 = t119 * t41 + t124 * t42;
t10 = qJDD(3) * pkin(9) + t12;
t54 = t124 * t57;
t200 = qJD(1) * t123;
t182 = t125 * t200;
t170 = t127 * t182;
t56 = t105 * t207 + t134 * t170 - t200 * t212;
t55 = qJD(3) * pkin(3) + t56;
t35 = t119 * t55 + t54;
t33 = qJD(3) * pkin(9) + t35;
t70 = t105 * t127 - t122 * t182 + qJD(4);
t20 = t130 * t70 + t133 * t33;
t69 = -qJDD(1) * t185 + t127 * t102 + qJDD(4);
t218 = t133 * t69;
t2 = -qJDD(5) * pkin(5) + qJD(5) * t20 + t10 * t130 - t218;
t233 = t106 * (pkin(10) * qJD(6) + t168 * qJD(3)) - t154 - t2;
t37 = t119 * t56 + t54;
t159 = -pkin(5) * t133 - pkin(10) * t130 - pkin(4);
t229 = pkin(3) * t124;
t84 = t159 - t229;
t96 = t168 * qJD(5);
t232 = t84 * t89 + (t37 - t96) * t106;
t74 = t119 * t208 - t122 * t202;
t76 = t90 * t127;
t25 = -t160 * t81 + t206 * t74 - t76 * t80;
t28 = -t160 * t83 - t211 * t74 - t76 * t82;
t48 = -t128 * t74 + (-t120 * t160 - t125 * t76) * t123;
t153 = g(1) * t28 + g(2) * t25 + g(3) * t48;
t18 = qJD(5) * pkin(10) + t20;
t108 = pkin(3) * t119 + pkin(9);
t216 = t106 * t108;
t230 = qJD(6) * (t18 + t216) - t153;
t199 = qJD(3) * t130;
t92 = t129 * t199 - t190;
t225 = t106 * t92;
t196 = qJD(5) * t129;
t94 = t132 * t199 + t196;
t224 = t106 * t94;
t53 = t119 * t57;
t179 = t133 * t188;
t58 = qJD(6) * t190 + (t179 + t187) * t132 + t237 * t129;
t223 = t129 * t58;
t221 = t130 * t69;
t220 = t130 * t79;
t217 = qJD(5) * t92;
t215 = t106 * t132;
t214 = t106 * t133;
t117 = t130 ^ 2;
t201 = -t133 ^ 2 + t117;
t197 = qJD(5) * t108;
t195 = qJD(5) * t130;
t194 = qJD(5) * t133;
t193 = qJD(6) * t129;
t191 = qJD(6) * t132;
t181 = t106 * t196;
t11 = -t119 * t42 + t124 * t41;
t34 = t124 * t55 - t53;
t174 = -t133 * t58 + t94 * t195;
t32 = -qJD(3) * pkin(4) - t34;
t173 = -qJD(3) * t32 - t10;
t169 = g(1) * t83 + g(2) * t81;
t21 = qJD(3) * t159 - t34;
t4 = t129 * t21 + t132 * t18;
t167 = t129 * t18 - t132 * t21;
t23 = t133 * t46 + t220;
t45 = t119 * t63 - t124 * t62;
t166 = t129 * t45 + t132 * t23;
t165 = -t129 * t23 + t132 * t45;
t67 = t127 * t130 + t133 * t75;
t164 = t129 * t74 + t132 * t67;
t163 = -t129 * t67 + t132 * t74;
t162 = t130 * t33 - t133 * t70;
t161 = t127 * t133 - t130 * t75;
t157 = t106 * t191 - t129 * t89;
t152 = -g(1) * t211 + g(2) * t206 - g(3) * t128;
t145 = t169 - t177;
t109 = -pkin(4) - t229;
t38 = t124 * t56 - t53;
t143 = -qJDD(5) * t108 + (qJD(3) * t109 + t32 + t38) * qJD(5);
t142 = qJD(6) * t84 * t106 - g(1) * t147 - g(2) * t146 - g(3) * t144;
t17 = -qJD(5) * pkin(5) + t162;
t141 = -pkin(10) * t89 + (-t17 + t162) * t106;
t1 = qJDD(5) * pkin(10) - qJD(5) * t162 + t10 * t133 + t221;
t140 = qJD(5) * t17 + qJD(6) * t21 - t106 * t38 - t108 * t89 + t1;
t135 = qJD(5) ^ 2;
t138 = -qJD(3) * t37 + t108 * t135 - t11 + t153 + (-pkin(4) + t109) * qJDD(3);
t137 = -g(1) * (t121 * t209 + t127 * t82) - g(2) * (-t122 * t206 + t127 * t80);
t136 = qJD(3) ^ 2;
t99 = qJDD(5) * t133 - t130 * t135;
t98 = qJDD(5) * t130 + t133 * t135;
t72 = t90 * t122 * qJD(3);
t71 = qJD(3) * t75;
t61 = t63 * qJD(3);
t60 = t62 * qJD(3);
t51 = qJD(5) * t67 - t130 * t72;
t50 = qJD(5) * t161 - t133 * t72;
t44 = -t119 * t61 + t124 * t60;
t43 = t119 * t60 + t124 * t61;
t31 = t133 * t144 + t220;
t16 = t130 * t65 + t133 * t147;
t14 = t130 * t64 + t133 * t146;
t8 = qJD(5) * t23 + t130 * t44;
t7 = t236 * qJD(5) + t133 * t44;
t6 = qJD(3) * t96 + qJDD(3) * t159 - t11;
t5 = t132 * t6;
t3 = [qJDD(1) - g(3), t102 * t128 - g(3) + (t120 ^ 2 + t125 ^ 2) * t123 ^ 2 * qJDD(1), 0, -qJD(3) * t61 + qJDD(3) * t62, -qJD(3) * t60 - qJDD(3) * t63, -t11 * t45 + t12 * t46 - t34 * t43 + t35 * t44 + t69 * t79 - g(3), 0, 0, 0, 0, 0, -t45 * t114 - qJD(5) * t8 + qJDD(5) * t236 + (-t133 * t43 + t195 * t45) * qJD(3), t45 * t187 - qJD(5) * t7 - qJDD(5) * t23 + (t130 * t43 + t194 * t45) * qJD(3), 0, 0, 0, 0, 0 -(-qJD(6) * t166 - t129 * t7 + t132 * t43) * t106 + t165 * t89 + t8 * t92 - t236 * t59 (qJD(6) * t165 + t129 * t43 + t132 * t7) * t106 - t166 * t89 + t8 * t94 - t236 * t58; 0, t152 + t102, 0 (qJDD(3) * t134 - t131 * t136) * t122 (-qJDD(3) * t131 - t134 * t136) * t122, -t11 * t74 + t12 * t75 + t127 * t69 - t34 * t71 - t35 * t72 + t152, 0, 0, 0, 0, 0, -t74 * t114 - qJD(5) * t51 + qJDD(5) * t161 + (-t133 * t71 + t195 * t74) * qJD(3), t74 * t187 - qJD(5) * t50 - qJDD(5) * t67 + (t130 * t71 + t194 * t74) * qJD(3), 0, 0, 0, 0, 0 -(-qJD(6) * t164 - t129 * t50 + t132 * t71) * t106 + t163 * t89 + t51 * t92 - t161 * t59 (qJD(6) * t163 + t129 * t71 + t132 * t50) * t106 - t164 * t89 + t51 * t94 - t161 * t58; 0, 0, qJDD(3), -g(3) * t62 + t145 * t131 + t137 * t134 + t226, g(3) * t63 + t56 * qJD(3) + ((-t105 * t122 - t170) * qJD(3) + t145) * t134 + (-t122 * t102 + (g(1) * t82 + g(2) * t80) * t127 + (t120 * t189 - t176 + (g(1) * t121 - g(2) * t126) * t122) * t123) * t131, t34 * t37 - t35 * t38 + (t11 * t124 + t12 * t119 + (g(3) * t213 + t169) * t131 + (-g(3) * (t122 * t128 + t123 * t204) + t137) * t134) * pkin(3), qJDD(3) * t117 + 0.2e1 * t130 * t179, 0.2e1 * t114 * t130 - 0.2e1 * t188 * t201, t98, t99, 0, t130 * t143 - t133 * t138, t130 * t138 + t133 * t143, t94 * t180 + (t132 * t58 - t193 * t94) * t130 (-t129 * t94 - t132 * t92) * t194 + (-t223 - t132 * t59 + (t129 * t92 - t132 * t94) * qJD(6)) * t130, t174 + t234 (t59 + t181) * t133 + (t157 - t217) * t130, -t106 * t195 - t133 * t89, t232 * t132 + t142 * t129 + (t140 * t129 + t230 * t132 + t92 * t197 - t5) * t133 + (t17 * t191 + t108 * t59 + t2 * t129 - t38 * t92 + (-t129 * t216 - t167) * qJD(5)) * t130, -t232 * t129 + t142 * t132 + (t94 * t197 + t140 * t132 + (-t230 + t6) * t129) * t133 + (-t17 * t193 + t108 * t58 + t2 * t132 - t38 * t94 + (-t108 * t215 - t4) * qJD(5)) * t130; 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t64 - g(3) * t79 + t69, 0, 0, 0, 0, 0, t99, -t98, 0, 0, 0, 0, 0 (-t59 + t181) * t133 + (t157 + t217) * t130, t174 - t234; 0, 0, 0, 0, 0, 0, -t130 * t136 * t133, t201 * t136, t187, t114, qJDD(5), t130 * t173 - t154 + t218, g(1) * t16 + g(2) * t14 + g(3) * t31 + t173 * t133 - t221, -t215 * t94 + t223 (t58 + t225) * t132 + (-t59 + t224) * t129 (-t130 * t94 + t132 * t214) * qJD(3) - t157, t106 * t193 + t219 + (-t129 * t214 + t130 * t92) * qJD(3), t106 * t199, -pkin(5) * t59 + t141 * t129 + t233 * t132 + t167 * t199 - t20 * t92, -pkin(5) * t58 - t233 * t129 + t141 * t132 + t4 * t199 - t20 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t92, -t92 ^ 2 + t94 ^ 2, t58 - t225, -t224 - t59, t89, -t129 * t1 + t5 - t17 * t94 - g(1) * (-t129 * t16 - t132 * t28) - g(2) * (-t129 * t14 - t132 * t25) - g(3) * (-t129 * t31 - t132 * t48) - t238 * t4, -t132 * t1 - t129 * t6 + t17 * t92 - g(1) * (t129 * t28 - t132 * t16) - g(2) * (t129 * t25 - t132 * t14) - g(3) * (t129 * t48 - t132 * t31) + t238 * t167;];
tau_reg  = t3;
