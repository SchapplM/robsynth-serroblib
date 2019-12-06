% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:36
% EndTime: 2019-12-05 16:37:45
% DurationCPUTime: 2.41s
% Computational Cost: add. (7028->308), mult. (14728->454), div. (0->0), fcn. (10581->12), ass. (0->186)
t158 = sin(qJ(5));
t162 = cos(qJ(3));
t190 = qJD(2) * qJD(3);
t183 = t162 * t190;
t159 = sin(qJ(3));
t189 = t159 * qJDD(2);
t126 = t183 + t189;
t153 = sin(pkin(10));
t155 = cos(pkin(10));
t178 = -t155 * qJDD(3) + t153 * t126;
t103 = qJDD(5) + t178;
t194 = qJD(2) * t159;
t122 = t153 * qJD(3) + t155 * t194;
t161 = cos(qJ(5));
t193 = qJD(2) * t162;
t97 = t158 * t122 + t161 * t193;
t99 = t161 * t122 - t158 * t193;
t77 = t99 * t97;
t216 = t103 - t77;
t219 = t158 * t216;
t218 = t161 * t216;
t154 = sin(pkin(5));
t156 = cos(pkin(5));
t201 = sin(pkin(9));
t202 = cos(pkin(9));
t170 = t201 * g(1) - t202 * g(2);
t169 = t156 * t170;
t195 = -g(3) + qJDD(1);
t217 = t154 * t195 + t169;
t141 = t159 * t190;
t188 = t162 * qJDD(2);
t127 = -t141 + t188;
t120 = -t155 * qJD(3) + t153 * t194;
t198 = t122 * t120;
t171 = -t127 - t198;
t215 = t153 * t171;
t214 = t155 * t171;
t110 = t122 * t193;
t85 = -t178 - t110;
t117 = qJD(5) + t120;
t105 = t153 * qJDD(3) + t155 * t126;
t179 = t158 * t105 + t161 * t127;
t47 = (qJD(5) - t117) * t99 + t179;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t114 = t117 ^ 2;
t118 = t120 ^ 2;
t119 = t122 ^ 2;
t213 = 2 * qJD(4);
t212 = pkin(4) * t153;
t168 = -t154 * t170 + t156 * t195;
t102 = t162 * t168;
t164 = qJD(3) ^ 2;
t177 = -qJDD(3) * pkin(3) - t164 * qJ(4) + qJDD(4) - t102;
t175 = -t162 * pkin(3) - t159 * qJ(4);
t124 = t175 * qJD(2);
t165 = qJD(2) ^ 2;
t131 = -t202 * g(1) - t201 * g(2);
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t92 = t163 * t131 + t217 * t160;
t79 = -t165 * pkin(2) + qJDD(2) * pkin(7) + t92;
t180 = qJD(2) * t124 + t79;
t54 = t180 * t159 + t177;
t211 = t153 * t54;
t89 = t127 - t198;
t210 = t153 * t89;
t209 = t155 * t54;
t208 = t155 * t89;
t150 = t162 ^ 2;
t147 = t150 * t165;
t174 = t126 + t183;
t176 = t160 * t131 - t217 * t163;
t78 = -qJDD(2) * pkin(2) - t165 * pkin(7) + t176;
t166 = -t174 * qJ(4) + (-t127 + t141) * pkin(3) + t78;
t167 = t159 * t168;
t55 = -t164 * pkin(3) + qJDD(3) * qJ(4) + t180 * t162 + t167;
t182 = t153 * t55 - t155 * t166;
t94 = t120 * pkin(4) - t122 * pkin(8);
t21 = t127 * pkin(4) - pkin(8) * t147 + (t213 + t94) * t122 + t182;
t207 = t158 * t21;
t68 = t103 + t77;
t206 = t158 * t68;
t205 = t159 * t79;
t204 = t161 * t21;
t203 = t161 * t68;
t200 = t117 * t158;
t199 = t117 * t161;
t139 = t159 * t165 * t162;
t132 = qJDD(3) + t139;
t197 = t159 * t132;
t133 = qJDD(3) - t139;
t196 = t162 * t133;
t191 = qJD(5) + t117;
t187 = t153 * t77;
t186 = t155 * t77;
t32 = -0.2e1 * qJD(4) * t120 + t153 * t166 + t155 * t55;
t185 = t162 * t198;
t184 = -pkin(4) * t155 - pkin(3);
t109 = t120 * t193;
t31 = t122 * t213 + t182;
t15 = t153 * t31 + t155 * t32;
t22 = -pkin(4) * t147 - t127 * pkin(8) - t120 * t94 + t32;
t34 = t178 * pkin(4) - t105 * pkin(8) + t205 + (t124 * t159 + (-pkin(4) * t122 - pkin(8) * t120) * t162) * qJD(2) + t177;
t11 = t158 * t22 - t161 * t34;
t65 = -t102 + t205;
t66 = t162 * t79 + t167;
t35 = t159 * t65 + t162 * t66;
t12 = t158 * t34 + t161 * t22;
t5 = -t161 * t11 + t158 * t12;
t6 = t158 * t11 + t161 * t12;
t14 = t153 * t32 - t155 * t31;
t173 = -t161 * t105 + t158 * t127;
t172 = -pkin(2) + t175;
t71 = -t97 * qJD(5) - t173;
t149 = t159 ^ 2;
t146 = t149 * t165;
t137 = -t147 - t164;
t136 = -t146 - t164;
t130 = t146 + t147;
t129 = (t149 + t150) * qJDD(2);
t128 = -0.2e1 * t141 + t188;
t125 = 0.2e1 * t183 + t189;
t116 = t162 * t127;
t108 = -t119 - t147;
t107 = -t119 + t147;
t106 = t118 - t147;
t101 = -t159 * t136 - t196;
t100 = t162 * t137 - t197;
t93 = -t147 - t118;
t87 = -t105 + t109;
t86 = t105 + t109;
t84 = -t110 + t178;
t83 = t117 * t97;
t82 = -t96 + t114;
t81 = t95 - t114;
t80 = -t118 - t119;
t76 = t96 - t95;
t75 = -t96 - t114;
t74 = -t114 - t95;
t73 = -t153 * t108 + t208;
t72 = t155 * t108 + t210;
t70 = -t99 * qJD(5) - t179;
t64 = t95 + t96;
t62 = t155 * t93 - t215;
t61 = t153 * t93 + t214;
t60 = (t158 * t99 - t161 * t97) * t117;
t59 = -t153 * t87 + t155 * t85;
t58 = t153 * t85 + t155 * t87;
t52 = t191 * t97 + t173;
t51 = t71 + t83;
t50 = t71 - t83;
t49 = -t191 * t99 - t179;
t46 = t159 * t86 + t162 * t73;
t45 = t161 * t71 - t99 * t200;
t44 = -t158 * t70 + t97 * t199;
t43 = t161 * t81 - t206;
t42 = -t158 * t82 + t218;
t41 = t159 * t84 + t162 * t62;
t40 = t159 * t80 + t162 * t59;
t39 = -t158 * t75 - t203;
t38 = t161 * t75 - t206;
t37 = t161 * t74 - t219;
t36 = t158 * t74 + t218;
t29 = -t158 * t50 + t161 * t49;
t28 = t158 * t51 - t161 * t47;
t27 = -t158 * t47 - t161 * t51;
t26 = -t153 * t52 + t155 * t39;
t25 = t153 * t39 + t155 * t52;
t24 = -t153 * t49 + t155 * t37;
t23 = t153 * t37 + t155 * t49;
t20 = -t153 * t64 + t155 * t28;
t19 = t153 * t28 + t155 * t64;
t18 = t159 * t38 + t162 * t26;
t17 = -pkin(8) * t38 + t204;
t16 = t159 * t36 + t162 * t24;
t13 = -pkin(8) * t36 + t207;
t10 = t162 * t15 + t159 * t54;
t9 = t159 * t27 + t162 * t20;
t8 = -pkin(4) * t38 + t12;
t7 = -pkin(4) * t36 + t11;
t4 = -pkin(8) * t27 - t5;
t3 = t153 * t21 + t155 * t6;
t2 = t153 * t6 - t155 * t21;
t1 = t159 * t5 + t162 * t3;
t30 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0, 0, 0, 0, 0, 0, (qJDD(2) * t163 - t160 * t165) * t154, (-qJDD(2) * t160 - t163 * t165) * t154, 0, t156 ^ 2 * t195 + (t160 * t92 - t163 * t176 - t169) * t154, 0, 0, 0, 0, 0, 0, t156 * (t162 * t132 + t159 * t137) + (t160 * t100 + t163 * t128) * t154, t156 * (-t159 * t133 + t162 * t136) + (t160 * t101 - t163 * t125) * t154, (t129 * t160 + t130 * t163) * t154, t156 * (t159 * t66 - t162 * t65) + (t160 * t35 - t163 * t78) * t154, 0, 0, 0, 0, 0, 0, t156 * (t159 * t62 - t162 * t84) + (t160 * t41 - t163 * t61) * t154, t156 * (t159 * t73 - t162 * t86) + (t160 * t46 - t163 * t72) * t154, t156 * (t159 * t59 - t162 * t80) + (t160 * t40 - t163 * t58) * t154, t156 * (t159 * t15 - t162 * t54) + (t160 * t10 - t163 * t14) * t154, 0, 0, 0, 0, 0, 0, t156 * (t159 * t24 - t162 * t36) + (t160 * t16 - t163 * t23) * t154, t156 * (t159 * t26 - t162 * t38) + (t160 * t18 - t163 * t25) * t154, t156 * (t159 * t20 - t162 * t27) + (t160 * t9 - t163 * t19) * t154, t156 * (t159 * t3 - t162 * t5) + (t160 * t1 - t163 * t2) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t176, -t92, 0, 0, t174 * t159, t162 * t125 + t159 * t128, t197 + t162 * (-t146 + t164), -t159 * t183 + t116, t159 * (t147 - t164) + t196, 0, pkin(2) * t128 + pkin(7) * t100 - t162 * t78, -pkin(2) * t125 + pkin(7) * t101 + t159 * t78, pkin(2) * t130 + pkin(7) * t129 + t35, -pkin(2) * t78 + pkin(7) * t35, t159 * (t155 * t105 + t153 * t110) - t185, t159 * (-t153 * t86 - t155 * t84) + t162 * (-t119 + t118), t159 * (-t153 * t107 + t214) + t162 * t87, t159 * (-t155 * t109 + t153 * t178) + t185, t159 * (t155 * t106 + t210) - t162 * t85, t116 + t159 * (t120 * t155 - t122 * t153) * t193, t159 * (-qJ(4) * t61 + t211) + t162 * (-pkin(3) * t61 + t31) - pkin(2) * t61 + pkin(7) * t41, t159 * (-qJ(4) * t72 + t209) + t162 * (-pkin(3) * t72 + t32) - pkin(2) * t72 + pkin(7) * t46, pkin(7) * t40 - t14 * t159 + t172 * t58, pkin(7) * t10 + t14 * t172, t159 * (t155 * t45 + t187) + t162 * (-t158 * t71 - t99 * t199), t159 * (t153 * t76 + t155 * t29) + t162 * (-t158 * t49 - t161 * t50), t159 * (t153 * t51 + t155 * t42) + t162 * (-t161 * t82 - t219), t159 * (t155 * t44 - t187) + t162 * (-t161 * t70 - t97 * t200), t159 * (-t153 * t47 + t155 * t43) + t162 * (-t158 * t81 - t203), t159 * (t153 * t103 + t155 * t60) + t162 * (t158 * t97 + t161 * t99) * t117, t159 * (-qJ(4) * t23 + t155 * t13 - t153 * t7) + t162 * (-pkin(3) * t23 - pkin(4) * t49 - pkin(8) * t37 + t204) - pkin(2) * t23 + pkin(7) * t16, t159 * (-qJ(4) * t25 - t153 * t8 + t155 * t17) + t162 * (-pkin(3) * t25 - pkin(4) * t52 - pkin(8) * t39 - t207) - pkin(2) * t25 + pkin(7) * t18, t159 * (-qJ(4) * t19 + t155 * t4 + t27 * t212) + t162 * (-pkin(3) * t19 - pkin(4) * t64 - pkin(8) * t28 - t6) - pkin(2) * t19 + pkin(7) * t9, t159 * (-qJ(4) * t2 + (-pkin(8) * t155 + t212) * t5) + t162 * (-pkin(3) * t2 + pkin(4) * t21 - pkin(8) * t6) - pkin(2) * t2 + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t146 - t147, t189, t139, t188, qJDD(3), -t65, -t66, 0, 0, t153 * t105 - t155 * t110, -t153 * t84 + t155 * t86, t155 * t107 + t215, -t153 * t109 - t155 * t178, t153 * t106 - t208, (t120 * t153 + t122 * t155) * t193, -pkin(3) * t84 + qJ(4) * t62 - t209, -pkin(3) * t86 + qJ(4) * t73 + t211, -pkin(3) * t80 + qJ(4) * t59 + t15, -pkin(3) * t54 + qJ(4) * t15, t153 * t45 - t186, t153 * t29 - t155 * t76, t153 * t42 - t155 * t51, t153 * t44 + t186, t153 * t43 + t155 * t47, -t103 * t155 + t153 * t60, -pkin(3) * t36 + qJ(4) * t24 + t13 * t153 + t155 * t7, -pkin(3) * t38 + qJ(4) * t26 + t153 * t17 + t155 * t8, qJ(4) * t20 + t153 * t4 + t184 * t27, qJ(4) * t3 + (-pkin(8) * t153 + t184) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t86, t80, t54, 0, 0, 0, 0, 0, 0, t36, t38, t27, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, t51, -t77, -t47, t103, -t11, -t12, 0, 0;];
tauJ_reg = t30;
