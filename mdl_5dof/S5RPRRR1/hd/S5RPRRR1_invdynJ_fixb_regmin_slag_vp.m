% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:06
% EndTime: 2019-07-18 13:26:12
% DurationCPUTime: 1.83s
% Computational Cost: add. (1170->287), mult. (2956->405), div. (0->0), fcn. (2308->8), ass. (0->157)
t131 = qJD(1) * qJD(3);
t71 = cos(qJ(3));
t115 = t71 * t131;
t67 = sin(qJ(3));
t56 = t67 * qJDD(1);
t193 = t115 + t56;
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t192 = g(1) * t68 - g(2) * t72;
t143 = qJD(4) * t67;
t191 = -qJD(1) * t143 + qJDD(3);
t136 = t71 * qJD(1);
t108 = qJD(4) + t136;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t12 = (t108 * qJD(3) + t56) * t66 - t191 * t70;
t137 = t66 * qJD(3);
t148 = qJD(1) * t67;
t39 = t70 * t148 + t137;
t52 = -qJD(4) + t136;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t16 = t65 * t39 + t69 * t52;
t190 = t16 * t52;
t60 = t70 * qJD(3);
t37 = t66 * t148 - t60;
t32 = qJD(5) + t37;
t142 = qJD(4) * t70;
t155 = t71 * t69;
t163 = t67 * t65;
t91 = t70 * t155 + t163;
t22 = t91 * qJD(1);
t100 = t69 * t142 - t22;
t139 = qJD(5) * t65;
t80 = -t66 * t139 + t100;
t189 = t80 * t32;
t116 = t67 * t131;
t59 = t71 * qJDD(1);
t85 = t116 - t59;
t188 = t85 * qJ(2);
t99 = -qJDD(2) + t192;
t144 = qJD(4) * t66;
t36 = qJDD(4) + t85;
t187 = t52 * t144 + t70 * t36;
t182 = g(2) * t68;
t184 = g(1) * t72;
t186 = -t184 - t182;
t95 = t108 * qJD(2);
t11 = qJD(4) * t60 + t191 * t66 + t193 * t70;
t18 = t69 * t39 - t65 * t52;
t4 = t18 * qJD(5) + t65 * t11 - t69 * t36;
t153 = t72 * t66;
t157 = t70 * t71;
t29 = t68 * t157 - t153;
t183 = g(2) * t29;
t180 = g(3) * t67;
t138 = qJD(5) * t69;
t3 = t69 * t11 - t52 * t138 - t39 * t139 + t65 * t36;
t179 = t3 * t65;
t178 = t66 * t3;
t177 = t70 * t4;
t176 = t11 * t66;
t175 = t16 * t32;
t174 = t18 * t32;
t111 = t32 * t69;
t173 = t36 * t71;
t172 = t37 * t52;
t171 = t37 * t67;
t170 = t39 * t52;
t169 = t39 * t67;
t168 = t39 * t70;
t167 = t52 * t66;
t10 = qJDD(5) + t12;
t166 = t65 * t10;
t165 = t66 * t67;
t164 = t66 * t71;
t162 = t67 * t69;
t161 = t67 * t70;
t160 = t67 * t72;
t159 = t68 * t66;
t158 = t69 * t10;
t156 = t71 * t65;
t74 = qJD(1) ^ 2;
t154 = t71 * t74;
t152 = t72 * t70;
t112 = qJ(2) * t131;
t151 = -t70 * qJDD(2) - t112 * t165;
t63 = t67 ^ 2;
t150 = -t71 ^ 2 + t63;
t149 = qJ(2) * t67;
t147 = qJD(3) * t67;
t146 = qJD(3) * t69;
t145 = qJD(3) * t71;
t141 = qJD(4) * t71;
t135 = qJ(2) * qJD(1);
t120 = t71 * t135;
t35 = t66 * qJD(2) + t70 * t120;
t140 = qJD(5) * t35;
t134 = t63 * qJDD(1);
t133 = qJ(2) * qJDD(1);
t132 = qJ(2) * qJDD(3);
t130 = qJD(2) * qJD(1);
t129 = t52 * t157;
t128 = t70 * t156;
t127 = t18 * t136;
t125 = t66 * t143;
t124 = t69 * t144;
t123 = t37 * t148;
t122 = t32 * t135;
t121 = t67 * t135;
t119 = t18 * t144 - t3 * t70;
t117 = t67 * t130;
t114 = t71 * t133;
t109 = 0.2e1 * t115;
t107 = -qJD(5) + t60;
t106 = qJD(4) * t120;
t105 = t71 * t112;
t28 = t71 * t159 + t152;
t30 = -t71 * t153 + t68 * t70;
t103 = -g(1) * t30 + g(2) * t28;
t21 = qJD(1) * t128 - t69 * t148;
t101 = t65 * t142 - t21;
t98 = qJD(1) * t63 + t52 * t71;
t97 = -t130 + t182;
t96 = -t136 * t167 + t187;
t5 = (qJDD(2) - t106) * t66 + (t95 - t188) * t70;
t94 = qJD(5) * t121 + t5;
t92 = t68 * t162 - t29 * t65;
t27 = t69 * t161 - t156;
t90 = -t128 + t162;
t26 = t65 * t161 + t155;
t89 = t101 * t32;
t88 = t52 * t142 - t66 * t36;
t87 = -t32 * t138 - t166;
t86 = t98 * t70;
t84 = -t154 + (qJD(4) + t52) * qJD(1);
t20 = t65 * t121 + t69 * t35;
t83 = t71 * t137 + t67 * t142;
t82 = t186 + 0.2e1 * t130;
t81 = t97 + t184;
t73 = qJD(3) ^ 2;
t79 = -qJ(2) * t73 + t99;
t6 = t70 * t106 + (t95 + t114) * t66 + t151;
t78 = g(3) * t165 + t103 - t6;
t77 = t180 - t95;
t76 = -t82 - t133;
t75 = t134 - t173 + (-t52 + 0.2e1 * t136) * t147;
t2 = -t20 * qJD(5) + t133 * t162 - t65 * t5 + (t105 + t117) * t69;
t34 = -t70 * qJD(2) + t66 * t120;
t31 = t71 * t152 + t159;
t19 = t69 * t121 - t65 * t35;
t15 = t65 * t160 + t31 * t69;
t14 = t69 * t160 - t31 * t65;
t8 = t107 * t155 + (-t124 + (-qJD(5) * t70 + qJD(3)) * t65) * t67;
t7 = -t65 * t125 - t71 * t139 - t67 * t146 + (t67 * t138 + t65 * t145) * t70;
t1 = t94 * t69 + (t193 * qJ(2) + t117 - t140) * t65;
t9 = [qJDD(1), t192, -t186, t99, t82 + 0.2e1 * t133, -t76 * qJ(2), t67 * t109 + t134, -0.2e1 * t150 * t131 + 0.2e1 * t67 * t59, qJDD(3) * t67 + t73 * t71, qJDD(3) * t71 - t73 * t67, 0, -t67 * t132 + t79 * t71, -t71 * t132 - t79 * t67, t11 * t161 + (t60 * t71 - t125) * t39, (-t37 * t70 - t39 * t66) * t145 + (-t176 - t12 * t70 + (t37 * t66 - t168) * qJD(4)) * t67, (-t52 * t60 - t11) * t71 + (qJD(3) * t39 + t187) * t67, (t137 * t52 + t12) * t71 + (-qJD(3) * t37 + t88) * t67, -t147 * t52 - t173, -t34 * t147 + g(1) * t29 - g(2) * t31 + t6 * t71 + (t66 * t98 + t171) * qJD(2) + (qJD(4) * t86 + t12 * t67 + t145 * t37 + t66 * t75) * qJ(2), -t35 * t147 - g(1) * t28 - g(2) * t30 + t5 * t71 + (t86 + t169) * qJD(2) + (t11 * t67 - t144 * t98 + t145 * t39 + t70 * t75) * qJ(2), t18 * t8 + t3 * t27, -t8 * t16 - t18 * t7 - t3 * t26 - t27 * t4, t27 * t10 + t3 * t165 + t18 * t83 + t8 * t32, -t26 * t10 - t16 * t83 - t4 * t165 - t7 * t32, t10 * t165 + t32 * t83, t2 * t165 + t6 * t26 + t34 * t7 - g(1) * (-t68 * t163 - t29 * t69) - g(2) * t15 + t83 * t19 + (t16 * t164 + t32 * t90) * qJD(2) + ((t107 * t32 * t65 - t137 * t16 + t158) * t67 + ((t65 * t144 + t146) * t32 + t66 * t4 + (qJD(4) * t16 + t87) * t70) * t71) * qJ(2), -t1 * t165 + t6 * t27 + t34 * t8 + g(1) * t92 - g(2) * t14 - t83 * t20 + (t18 * t164 - t32 * t91) * qJD(2) + ((t107 * t111 - t137 * t18 - t166) * t67 + (-(qJD(3) * t65 - t124) * t32 + t178 + (qJD(4) * t18 + t139 * t32 - t158) * t70) * t71) * qJ(2); 0, 0, 0, -qJDD(1), -t74, -t74 * qJ(2) - t99, 0, 0, 0, 0, 0, -t59 + 0.2e1 * t116, t56 + t109, 0, 0, 0, 0, 0, t96 - t123, (-t129 - t169) * qJD(1) + t88, 0, 0, 0, 0, 0, -t177 - t89 + (t87 - t190) * t66, (-t127 - t158) * t66 - t189 + t119; 0, 0, 0, 0, 0, 0, -t67 * t154, t150 * t74, t56, t59, qJDD(3), -g(3) * t71 + t76 * t67, t76 * t71 + t180, -t52 * t168 + t176, (t11 + t172) * t70 + (-t12 + t170) * t66, (t129 - t169) * qJD(1) - t88, t96 + t123, t52 * t148, (-g(3) * t70 + (-t37 - t60) * t135) * t71 + (qJD(1) * t34 + t81 * t70 + (-qJDD(1) * t70 + t66 * t84) * qJ(2)) * t67, (g(3) * t66 + (-t39 + t137) * t135) * t71 + (qJD(1) * t35 - t81 * t66 + (qJDD(1) * t66 + t70 * t84) * qJ(2)) * t67, t69 * t178 + t18 * t80, t22 * t16 + t18 * t21 + (-t16 * t69 - t18 * t65) * t142 + (-t179 - t4 * t69 + (t16 * t65 - t18 * t69) * qJD(5)) * t66, (-t127 + t158) * t66 + t189 + t119, t177 - t89 + (t87 + t190) * t66, -t10 * t70 - t32 * t167, -t2 * t70 - g(3) * t91 + t101 * t34 - t26 * t122 + (t34 * t138 + t19 * qJD(4) + t6 * t65 + (t149 * t16 - t19 * t71) * qJD(1)) * t66 - t186 * t27, t1 * t70 - g(3) * t90 + t100 * t34 - t27 * t122 + (-t34 * t139 - t20 * qJD(4) + t6 * t69 + (t149 * t18 + t20 * t71) * qJD(1)) * t66 + t186 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t11 - t172, -t12 - t170, t36, -t35 * t52 + (-t141 * t70 - t169) * t135 + (t77 - t114) * t66 + t103 - t151, g(1) * t31 + t183 - t66 * qJDD(2) + t34 * t52 + (t141 * t66 + t171) * t135 + (t77 + t188) * t70, t111 * t18 + t179, (t3 - t175) * t69 + (-t4 - t174) * t65, t111 * t32 - t18 * t39 + t166, -t32 ^ 2 * t65 + t16 * t39 + t158, -t32 * t39, -t35 * t16 - t19 * t39 + t69 * t78, -t35 * t18 + t20 * t39 - t65 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t16, -t16 ^ 2 + t18 ^ 2, t3 + t175, t174 - t4, t10, -g(1) * t14 - g(2) * t92 + g(3) * t26 - t34 * t18 + t20 * t32 + t2, g(1) * t15 + g(3) * t27 + t34 * t16 + t19 * t32 + (-t94 + t183) * t69 + (-t105 + t140 + (t97 - t133) * t67) * t65;];
tau_reg  = t9;
