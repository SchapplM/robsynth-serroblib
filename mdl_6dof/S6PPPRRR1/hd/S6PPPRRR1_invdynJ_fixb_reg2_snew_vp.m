% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPPRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:34:10
% EndTime: 2019-05-04 19:34:14
% DurationCPUTime: 2.82s
% Computational Cost: add. (21377->267), mult. (36860->436), div. (0->0), fcn. (33667->18), ass. (0->201)
t153 = cos(pkin(8));
t144 = sin(pkin(14));
t150 = cos(pkin(14));
t148 = sin(pkin(7));
t154 = cos(pkin(7));
t146 = sin(pkin(12));
t152 = cos(pkin(12));
t123 = -t152 * g(1) - t146 * g(2);
t145 = sin(pkin(13));
t151 = cos(pkin(13));
t122 = t146 * g(1) - t152 * g(2);
t141 = -g(3) + qJDD(1);
t149 = sin(pkin(6));
t155 = cos(pkin(6));
t171 = t122 * t155 + t141 * t149;
t81 = -t145 * t123 + t171 * t151;
t99 = -t149 * t122 + t155 * t141 + qJDD(2);
t182 = t148 * t99 + t154 * t81;
t82 = t151 * t123 + t171 * t145;
t56 = -t144 * t82 + t182 * t150;
t217 = t153 * t56;
t147 = sin(pkin(8));
t70 = -t148 * t81 + t154 * t99 + qJDD(3);
t218 = t147 * t70;
t221 = t217 + t218;
t165 = qJD(4) ^ 2;
t158 = sin(qJ(6));
t159 = sin(qJ(5));
t204 = qJD(4) * qJD(5);
t133 = t159 * t204;
t162 = cos(qJ(5));
t135 = t162 * qJDD(4);
t116 = t135 - t133;
t109 = -qJDD(6) + t116;
t161 = cos(qJ(6));
t207 = qJD(4) * t159;
t110 = -t161 * qJD(5) + t158 * t207;
t112 = t158 * qJD(5) + t161 * t207;
t212 = t112 * t110;
t167 = -t109 - t212;
t220 = t158 * t167;
t219 = t161 * t167;
t130 = t162 * qJD(4) - qJD(6);
t201 = t162 * t204;
t203 = t159 * qJDD(4);
t115 = t201 + t203;
t199 = -t161 * qJDD(5) + t158 * t115;
t73 = (qJD(6) + t130) * t112 + t199;
t107 = t110 ^ 2;
t108 = t112 ^ 2;
t128 = t130 ^ 2;
t164 = qJD(5) ^ 2;
t197 = -t162 * pkin(5) - t159 * pkin(11);
t160 = sin(qJ(4));
t163 = cos(qJ(4));
t57 = t182 * t144 + t150 * t82;
t42 = t221 * t160 + t163 * t57;
t39 = -t165 * pkin(4) + qJDD(4) * pkin(10) + t42;
t200 = t165 * t197 + t39;
t183 = -t147 * t56 + t153 * t70;
t48 = t162 * t183;
t25 = -qJDD(5) * pkin(5) - t164 * pkin(11) + t200 * t159 - t48;
t216 = t158 * t25;
t84 = t109 - t212;
t215 = t158 * t84;
t214 = t161 * t25;
t213 = t161 * t84;
t211 = t130 * t158;
t210 = t130 * t161;
t129 = t159 * t165 * t162;
t124 = qJDD(5) + t129;
t209 = t159 * t124;
t125 = qJDD(5) - t129;
t208 = t162 * t125;
t206 = qJD(6) - t130;
t202 = t162 * t212;
t166 = t159 * t183;
t26 = -t164 * pkin(5) + qJDD(5) * pkin(11) + t200 * t162 + t166;
t192 = -t116 + t133;
t193 = t115 + t201;
t198 = t160 * t57 - t221 * t163;
t38 = -qJDD(4) * pkin(4) - t165 * pkin(10) + t198;
t35 = t192 * pkin(5) - t193 * pkin(11) + t38;
t20 = t158 * t26 - t161 * t35;
t21 = t158 * t35 + t161 * t26;
t12 = t158 * t20 + t161 * t21;
t28 = t159 * t39 - t48;
t29 = t162 * t39 + t166;
t16 = t159 * t28 + t162 * t29;
t11 = t158 * t21 - t161 * t20;
t9 = t162 * t12 + t159 * t25;
t195 = -t11 * t163 + t160 * t9;
t8 = t159 * t12 - t162 * t25;
t3 = -t147 * t8 + t195 * t153;
t5 = t160 * t11 + t163 * t9;
t196 = t144 * t5 + t150 * t3;
t13 = t163 * t16 + t160 * t38;
t15 = t159 * t29 - t162 * t28;
t181 = t16 * t160 - t163 * t38;
t7 = -t147 * t15 + t181 * t153;
t194 = t13 * t144 + t150 * t7;
t180 = t160 * t42 - t163 * t198;
t23 = t147 ^ 2 * t56 + (t180 - t218) * t153;
t24 = t160 * t198 + t163 * t42;
t191 = t144 * t24 + t150 * t23;
t102 = t110 * t130;
t168 = -t158 * qJDD(5) - t161 * t115;
t88 = -t110 * qJD(6) - t168;
t77 = -t102 + t88;
t59 = t158 * t77 - t161 * t73;
t83 = t107 + t108;
t50 = -t159 * t83 + t162 * t59;
t58 = -t158 * t73 - t161 * t77;
t179 = t160 * t50 - t163 * t58;
t49 = t159 * t59 + t162 * t83;
t31 = -t147 * t49 + t179 * t153;
t43 = t160 * t58 + t163 * t50;
t190 = t144 * t43 + t150 * t31;
t89 = -t128 - t107;
t62 = t161 * t89 - t220;
t74 = -t206 * t112 - t199;
t52 = -t159 * t74 + t162 * t62;
t61 = t158 * t89 + t219;
t178 = t160 * t52 - t163 * t61;
t51 = t159 * t62 + t162 * t74;
t33 = -t147 * t51 + t178 * t153;
t44 = t160 * t61 + t163 * t52;
t189 = t144 * t44 + t150 * t33;
t92 = -t108 - t128;
t69 = -t158 * t92 + t213;
t78 = t206 * t110 + t168;
t54 = -t159 * t78 + t162 * t69;
t68 = t161 * t92 + t215;
t177 = t160 * t54 - t163 * t68;
t53 = t159 * t69 + t162 * t78;
t37 = -t147 * t53 + t177 * t153;
t45 = t160 * t68 + t163 * t54;
t188 = t144 * t45 + t150 * t37;
t187 = t144 * t57 + t150 * t56;
t117 = t135 - 0.2e1 * t133;
t140 = t162 ^ 2;
t138 = t140 * t165;
t127 = -t138 - t164;
t97 = t162 * t127 - t209;
t175 = t117 * t163 + t160 * t97;
t95 = t162 * t124 + t159 * t127;
t66 = -t147 * t95 + t175 * t153;
t79 = -t160 * t117 + t163 * t97;
t186 = t144 * t79 + t150 * t66;
t114 = 0.2e1 * t201 + t203;
t139 = t159 ^ 2;
t136 = t139 * t165;
t126 = -t136 - t164;
t98 = -t159 * t126 - t208;
t176 = -t114 * t163 + t160 * t98;
t96 = -t159 * t125 + t162 * t126;
t67 = -t147 * t96 + t176 * t153;
t80 = t160 * t114 + t163 * t98;
t185 = t144 * t80 + t150 * t67;
t118 = (t139 + t140) * qJDD(4);
t121 = t136 + t138;
t172 = t118 * t160 + t121 * t163;
t91 = t172 * t153;
t93 = t163 * t118 - t160 * t121;
t184 = t144 * t93 + t150 * t91;
t120 = -t160 * qJDD(4) - t163 * t165;
t105 = t120 * t153;
t169 = t163 * qJDD(4) - t160 * t165;
t174 = t105 * t150 - t144 * t169;
t106 = t169 * t153;
t173 = t106 * t150 + t120 * t144;
t170 = -pkin(4) + t197;
t104 = t169 * t147;
t103 = t120 * t147;
t101 = -t108 + t128;
t100 = t107 - t128;
t94 = t108 - t107;
t90 = t172 * t147;
t87 = -t112 * qJD(6) - t199;
t76 = t102 + t88;
t72 = t154 * t104 + t173 * t148;
t71 = t154 * t103 + t174 * t148;
t65 = t176 * t147 + t153 * t96;
t64 = t175 * t147 + t153 * t95;
t60 = t184 * t148 + t154 * t90;
t47 = t185 * t148 + t154 * t65;
t46 = t186 * t148 + t154 * t64;
t40 = t187 * t148 + t154 * t70;
t36 = t177 * t147 + t153 * t53;
t32 = t178 * t147 + t153 * t51;
t30 = t179 * t147 + t153 * t49;
t22 = t153 ^ 2 * t70 + (t180 - t217) * t147;
t18 = t188 * t148 + t154 * t36;
t17 = t189 * t148 + t154 * t32;
t14 = t190 * t148 + t154 * t30;
t10 = t191 * t148 + t154 * t22;
t6 = t181 * t147 + t153 * t15;
t4 = t194 * t148 + t154 * t6;
t2 = t195 * t147 + t153 * t8;
t1 = t196 * t148 + t154 * t2;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 * t99 + (t145 * t82 + t151 * t81) * t149, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 * t40 + (t145 * (-t144 * t56 + t150 * t57) + t151 * (-t148 * t70 + t187 * t154)) * t149, 0, 0, 0, 0, 0, 0, t155 * t72 + (t145 * (-t144 * t106 + t150 * t120) + t151 * (-t148 * t104 + t173 * t154)) * t149, t155 * t71 + (t145 * (-t144 * t105 - t150 * t169) + t151 * (-t148 * t103 + t174 * t154)) * t149, 0, t155 * t10 + (t145 * (-t144 * t23 + t150 * t24) + t151 * (-t148 * t22 + t191 * t154)) * t149, 0, 0, 0, 0, 0, 0, t155 * t46 + (t145 * (-t144 * t66 + t150 * t79) + t151 * (-t148 * t64 + t186 * t154)) * t149, t155 * t47 + (t145 * (-t144 * t67 + t150 * t80) + t151 * (-t148 * t65 + t185 * t154)) * t149, t155 * t60 + (t145 * (-t144 * t91 + t150 * t93) + t151 * (-t148 * t90 + t184 * t154)) * t149, t155 * t4 + (t145 * (t150 * t13 - t144 * t7) + t151 * (-t148 * t6 + t194 * t154)) * t149, 0, 0, 0, 0, 0, 0, t155 * t17 + (t145 * (-t144 * t33 + t150 * t44) + t151 * (-t148 * t32 + t189 * t154)) * t149, t155 * t18 + (t145 * (-t144 * t37 + t150 * t45) + t151 * (-t148 * t36 + t188 * t154)) * t149, t155 * t14 + (t145 * (-t144 * t31 + t150 * t43) + t151 * (-t148 * t30 + t190 * t154)) * t149, t155 * t1 + (t145 * (-t144 * t3 + t150 * t5) + t151 * (-t148 * t2 + t196 * t154)) * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t72, t71, 0, t10, 0, 0, 0, 0, 0, 0, t46, t47, t60, t4, 0, 0, 0, 0, 0, 0, t17, t18, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, t104, t103, 0, t22, 0, 0, 0, 0, 0, 0, t64, t65, t90, t6, 0, 0, 0, 0, 0, 0, t32, t36, t30, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t198, -t42, 0, 0, t193 * t159, t162 * t114 + t159 * t117, t209 + t162 * (-t136 + t164), -t192 * t162, t159 * (t138 - t164) + t208, 0, pkin(4) * t117 + pkin(10) * t97 - t162 * t38, -pkin(4) * t114 + pkin(10) * t98 + t159 * t38, pkin(4) * t121 + pkin(10) * t118 + t16, -pkin(4) * t38 + pkin(10) * t16, t159 * (t112 * t211 + t161 * t88) - t202, t159 * (-t158 * t76 + t161 * t74) - t162 * t94, t159 * (-t158 * t101 + t219) - t162 * t77, t159 * (-t110 * t210 - t158 * t87) + t202, t159 * (t161 * t100 + t215) + t162 * t73, t162 * t109 + t159 * (t110 * t161 - t112 * t158) * t130, t159 * (-pkin(11) * t61 + t216) + t162 * (-pkin(5) * t61 + t20) - pkin(4) * t61 + pkin(10) * t52, t159 * (-pkin(11) * t68 + t214) + t162 * (-pkin(5) * t68 + t21) - pkin(4) * t68 + pkin(10) * t54, pkin(10) * t50 - t159 * t11 + t170 * t58, pkin(10) * t9 + t170 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t136 - t138, t203, t129, t135, qJDD(5), -t28, -t29, 0, 0, -t112 * t210 + t158 * t88, t158 * t74 + t161 * t76, t161 * t101 + t220, -t110 * t211 + t161 * t87, t158 * t100 - t213, (t110 * t158 + t112 * t161) * t130, pkin(5) * t74 + pkin(11) * t62 - t214, pkin(5) * t78 + pkin(11) * t69 + t216, pkin(5) * t83 + pkin(11) * t59 + t12, -pkin(5) * t25 + pkin(11) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t94, t77, -t212, -t73, -t109, -t20, -t21, 0, 0;];
tauJ_reg  = t19;
