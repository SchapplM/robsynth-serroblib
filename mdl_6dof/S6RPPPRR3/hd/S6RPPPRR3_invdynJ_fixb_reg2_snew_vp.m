% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:41:17
% EndTime: 2019-05-05 13:41:29
% DurationCPUTime: 3.79s
% Computational Cost: add. (10770->334), mult. (22213->467), div. (0->0), fcn. (13524->10), ass. (0->193)
t159 = sin(qJ(5));
t155 = cos(pkin(10));
t162 = cos(qJ(5));
t153 = sin(pkin(10));
t200 = qJD(1) * t153;
t128 = -t162 * t155 * qJD(1) + t159 * t200;
t177 = t153 * t162 + t155 * t159;
t129 = t177 * qJD(1);
t205 = t128 * t129;
t227 = qJDD(5) - t205;
t229 = t159 * t227;
t228 = t162 * t227;
t158 = sin(qJ(6));
t161 = cos(qJ(6));
t110 = -t161 * qJD(5) - t129 * t158;
t112 = qJD(5) * t158 - t129 * t161;
t208 = t112 * t110;
t192 = t155 * qJDD(1);
t193 = t153 * qJDD(1);
t178 = t159 * t193 - t162 * t192;
t199 = qJD(5) * t129;
t100 = t199 + t178;
t94 = -qJDD(6) + t100;
t170 = -t94 - t208;
t226 = t158 * t170;
t225 = t161 * t170;
t148 = t153 ^ 2;
t149 = t155 ^ 2;
t201 = t148 + t149;
t222 = qJD(1) ^ 2;
t136 = t201 * t222;
t127 = t177 * qJDD(1);
t122 = -qJD(6) + t128;
t198 = t128 * qJD(5);
t103 = t198 - t127;
t184 = -t161 * qJDD(5) + t158 * t103;
t57 = (qJD(6) + t122) * t112 + t184;
t154 = sin(pkin(9));
t188 = t154 * qJ(2) + pkin(3);
t150 = qJDD(1) * qJ(2);
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t181 = t163 * g(1) + t160 * g(2);
t175 = 0.2e1 * qJD(2) * qJD(1) - t181;
t173 = t150 + t175;
t221 = pkin(1) + pkin(2);
t114 = -t221 * t222 + t173;
t156 = cos(pkin(9));
t187 = t160 * g(1) - t163 * g(2);
t179 = qJDD(2) - t187;
t172 = -t222 * qJ(2) + t179;
t169 = -t221 * qJDD(1) + t172;
t183 = t154 * t114 - t156 * t169;
t174 = -t222 * qJ(4) + qJDD(4) + t183;
t85 = qJDD(1) * pkin(3) + t174;
t224 = t188 * qJDD(1) + t85;
t108 = t110 ^ 2;
t109 = t112 ^ 2;
t121 = t122 ^ 2;
t125 = t128 ^ 2;
t126 = t129 ^ 2;
t223 = 0.2e1 * t153;
t220 = t155 * pkin(4);
t203 = t156 * t114;
t171 = -t154 * t179 - t203;
t182 = -t154 * t221 - qJ(4);
t168 = t182 * qJDD(1) - t188 * t222 - t171;
t202 = g(3) + qJDD(3);
t186 = t155 * t202;
t194 = qJD(4) * qJD(1);
t167 = pkin(7) * t193 + t194 * t223 + t186 + (t220 * t222 - t168) * t153;
t197 = t149 * t222;
t77 = t153 * t202 + (t168 - 0.2e1 * t194) * t155;
t69 = -pkin(4) * t197 - pkin(7) * t192 + t77;
t42 = t159 * t69 - t162 * t167;
t204 = t153 * t159;
t212 = t162 * t69;
t43 = t212 + t159 * t186 + ((pkin(7) - t182) * qJDD(1) + (0.2e1 * qJD(4) + (t188 + t220) * qJD(1)) * qJD(1) + t171) * t204;
t23 = t159 * t43 - t162 * t42;
t219 = t153 * t23;
t164 = qJD(5) ^ 2;
t95 = -pkin(5) * t128 + pkin(8) * t129;
t31 = -qJDD(5) * pkin(5) - pkin(8) * t164 - t129 * t95 + t42;
t218 = t158 * t31;
t66 = t94 - t208;
t217 = t158 * t66;
t75 = (pkin(3) + t220) * qJDD(1) + t174 + (-t148 * t222 - t197) * pkin(7);
t216 = t159 * t75;
t98 = qJDD(5) + t205;
t215 = t159 * t98;
t214 = t161 * t31;
t213 = t161 * t66;
t211 = t162 * t75;
t210 = t162 * t98;
t209 = qJDD(1) * pkin(1);
t207 = t122 * t158;
t206 = t122 * t161;
t196 = qJD(6) - t122;
t191 = t159 * t208;
t190 = t162 * t208;
t189 = -pkin(5) * t162 - pkin(4);
t32 = -t164 * pkin(5) + qJDD(5) * pkin(8) + t128 * t95 + t159 * t167 + t212;
t40 = (-t103 - t198) * pkin(8) + (-t100 - t199) * pkin(5) + t75;
t19 = t158 * t32 - t161 * t40;
t20 = t158 * t40 + t161 * t32;
t8 = t158 * t19 + t161 * t20;
t24 = t159 * t42 + t162 * t43;
t185 = qJ(2) * t156 - qJ(4);
t3 = t159 * t8 - t162 * t31;
t4 = t159 * t31 + t162 * t8;
t180 = t153 * t3 - t155 * t4;
t76 = t153 * (t203 + t154 * (-qJDD(1) * pkin(2) + t179 - t209) - qJDD(1) * qJ(4)) - t186 + (-t188 * qJD(1) - 0.2e1 * qJD(4)) * t200;
t49 = t153 * t76 + t155 * t77;
t7 = t158 * t20 - t161 * t19;
t176 = -t158 * qJDD(5) - t161 * t103;
t79 = -qJD(6) * t110 - t176;
t89 = t154 * t169 + t203;
t135 = qJDD(1) * t156 + t222 * t154;
t134 = -t154 * qJDD(1) + t222 * t156;
t133 = t201 * qJDD(1);
t132 = t155 * t136;
t131 = t153 * t136;
t124 = -t172 + t209;
t117 = -t126 - t164;
t116 = -t126 + t164;
t115 = t125 - t164;
t107 = -t132 * t154 - t156 * t192;
t106 = t131 * t154 + t156 * t193;
t105 = -t133 * t154 + t136 * t156;
t102 = t127 - 0.2e1 * t198;
t101 = 0.2e1 * t199 + t178;
t96 = -t164 - t125;
t93 = t110 * t122;
t92 = -t109 + t121;
t91 = t108 - t121;
t90 = -t125 - t126;
t87 = t109 - t108;
t86 = t100 - t199;
t83 = -t109 - t121;
t82 = -t117 * t159 - t210;
t81 = t117 * t162 - t215;
t80 = -t121 - t108;
t78 = -qJD(6) * t112 - t184;
t74 = t108 + t109;
t73 = -t127 * t159 + t162 * t86;
t72 = t127 * t162 + t159 * t86;
t71 = t162 * t96 - t229;
t70 = t159 * t96 + t228;
t65 = (t110 * t161 - t112 * t158) * t122;
t63 = t154 * t89 - t156 * t183;
t62 = t196 * t110 + t176;
t61 = t79 - t93;
t60 = t79 + t93;
t58 = -t196 * t112 - t184;
t56 = t112 * t207 + t161 * t79;
t55 = -t110 * t206 - t158 * t78;
t54 = -t153 * t81 + t155 * t82;
t53 = t161 * t91 + t217;
t52 = -t158 * t92 + t225;
t51 = -t158 * t83 + t213;
t50 = t161 * t83 + t217;
t48 = -t153 * t72 + t155 * t73;
t47 = t161 * t80 - t226;
t46 = t158 * t80 + t225;
t45 = -t153 * t70 + t155 * t71;
t44 = t102 * t156 + t154 * t54;
t38 = t101 * t156 + t154 * t45;
t37 = t154 * t48 - t156 * t90;
t36 = t154 * t49 - t156 * t85;
t35 = t158 * t61 - t161 * t57;
t34 = -t158 * t60 + t161 * t58;
t33 = -t158 * t57 - t161 * t61;
t30 = -t159 * t62 + t162 * t51;
t29 = t159 * t51 + t162 * t62;
t28 = -t159 * t58 + t162 * t47;
t27 = t159 * t47 + t162 * t58;
t26 = -t159 * t74 + t162 * t35;
t25 = t159 * t35 + t162 * t74;
t22 = -pkin(8) * t50 + t214;
t21 = -pkin(8) * t46 + t218;
t17 = -t153 * t29 + t155 * t30;
t16 = -t153 * t27 + t155 * t28;
t15 = -pkin(5) * t50 + t20;
t14 = -pkin(5) * t46 + t19;
t13 = -t153 * t25 + t155 * t26;
t12 = t154 * t17 - t156 * t50;
t11 = t154 * t16 - t156 * t46;
t10 = t155 * t24 - t219;
t9 = t10 * t154 - t156 * t75;
t6 = t13 * t154 - t156 * t33;
t5 = -pkin(8) * t33 - t7;
t1 = -t154 * t180 - t156 * t7;
t2 = [0, 0, 0, 0, 0, qJDD(1), t187, t181, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t179 + 0.2e1 * t209, 0, 0.2e1 * t150 + t175, pkin(1) * t124 + qJ(2) * (-t222 * pkin(1) + t173), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t134 + t221 * t135 + t183, qJ(2) * t135 + t221 * t134 + t89, 0, qJ(2) * (t154 * t183 + t156 * t89) - t221 * t63, t148 * qJDD(1), t192 * t223, 0, t149 * qJDD(1), 0, 0, -t221 * t107 - t185 * t132 + t224 * t155, -t221 * t106 + t185 * t131 - t224 * t153, qJ(2) * (-t133 * t156 - t136 * t154) - pkin(3) * t136 + qJ(4) * t133 - t221 * t105 - t49, qJ(2) * (t154 * t85 + t156 * t49) + pkin(3) * t85 - qJ(4) * t49 - t221 * t36, -t153 * (t103 * t162 + t159 * t199) - t155 * (t103 * t159 - t162 * t199), -t153 * (t101 * t162 + t102 * t159) - t155 * (t101 * t159 - t102 * t162), -t153 * (-t116 * t159 + t228) - t155 * (t116 * t162 + t229), -t153 * (-t100 * t159 - t162 * t198) - t155 * (t100 * t162 - t159 * t198), -t153 * (t115 * t162 - t215) - t155 * (t115 * t159 + t210), (-t153 * (t128 * t162 - t129 * t159) - t155 * (t128 * t159 + t129 * t162)) * qJD(5), qJ(2) * (-t101 * t154 + t156 * t45) - t153 * (-pkin(7) * t70 + t216) - t155 * (pkin(4) * t101 + pkin(7) * t71 - t211) - pkin(3) * t101 - qJ(4) * t45 - t221 * t38, qJ(2) * (-t102 * t154 + t156 * t54) - t153 * (-pkin(7) * t81 + t211) - t155 * (pkin(4) * t102 + pkin(7) * t82 + t216) - pkin(3) * t102 - qJ(4) * t54 - t221 * t44, qJ(2) * (t154 * t90 + t156 * t48) - t153 * (-pkin(7) * t72 - t23) - t155 * (-pkin(4) * t90 + pkin(7) * t73 + t24) + pkin(3) * t90 - qJ(4) * t48 - t221 * t37, qJ(2) * (t10 * t156 + t154 * t75) + pkin(7) * t219 - t155 * (-pkin(4) * t75 + pkin(7) * t24) + pkin(3) * t75 - qJ(4) * t10 - t221 * t9, -t153 * (t162 * t56 + t191) - t155 * (t159 * t56 - t190), -t153 * (t159 * t87 + t162 * t34) - t155 * (t159 * t34 - t162 * t87), -t153 * (t159 * t61 + t162 * t52) - t155 * (t159 * t52 - t162 * t61), -t153 * (t162 * t55 - t191) - t155 * (t159 * t55 + t190), -t153 * (-t159 * t57 + t162 * t53) - t155 * (t159 * t53 + t162 * t57), -t153 * (-t159 * t94 + t162 * t65) - t155 * (t159 * t65 + t162 * t94), qJ(2) * (t154 * t46 + t156 * t16) - t153 * (-pkin(7) * t27 - t14 * t159 + t162 * t21) - t155 * (-pkin(4) * t46 + pkin(7) * t28 + t14 * t162 + t159 * t21) + pkin(3) * t46 - qJ(4) * t16 - t221 * t11, qJ(2) * (t154 * t50 + t156 * t17) - t153 * (-pkin(7) * t29 - t15 * t159 + t162 * t22) - t155 * (-pkin(4) * t50 + pkin(7) * t30 + t15 * t162 + t159 * t22) + pkin(3) * t50 - qJ(4) * t17 - t221 * t12, -t153 * (-pkin(7) * t25 + t162 * t5) - t155 * (pkin(7) * t26 + t159 * t5) - t221 * t6 + t185 * t13 + (-pkin(5) * t204 - t155 * t189 + t188) * t33, -t221 * t1 + (-t153 * (pkin(5) * t159 - pkin(8) * t162) - t155 * (-pkin(8) * t159 + t189) + t188) * t7 + (-t185 + pkin(7)) * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t222, -t124, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, t63, 0, 0, 0, 0, 0, 0, t107, t106, t105, t36, 0, 0, 0, 0, 0, 0, t38, t44, t37, t9, 0, 0, 0, 0, 0, 0, t11, t12, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153 * t77 - t155 * t76, 0, 0, 0, 0, 0, 0, t153 * t71 + t155 * t70, t153 * t82 + t155 * t81, t153 * t73 + t155 * t72, t153 * t24 + t155 * t23, 0, 0, 0, 0, 0, 0, t153 * t28 + t155 * t27, t153 * t30 + t155 * t29, t153 * t26 + t155 * t25, t153 * t4 + t155 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, -t193, -t136, t85, 0, 0, 0, 0, 0, 0, -t101, -t102, t90, t75, 0, 0, 0, 0, 0, 0, t46, t50, t33, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, t126 - t125, -t127, -t205, t178, qJDD(5), -t42, -t43, 0, 0, -t112 * t206 + t158 * t79, t158 * t58 + t161 * t60, t161 * t92 + t226, -t110 * t207 + t161 * t78, t158 * t91 - t213, (t110 * t158 + t112 * t161) * t122, pkin(5) * t58 + pkin(8) * t47 - t214, pkin(5) * t62 + pkin(8) * t51 + t218, pkin(5) * t74 + pkin(8) * t35 + t8, -pkin(5) * t31 + pkin(8) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t87, t61, -t208, -t57, -t94, -t19, -t20, 0, 0;];
tauJ_reg  = t2;
