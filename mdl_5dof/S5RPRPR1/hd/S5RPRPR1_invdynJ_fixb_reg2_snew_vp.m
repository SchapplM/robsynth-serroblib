% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:44
% EndTime: 2019-12-05 17:47:54
% DurationCPUTime: 2.64s
% Computational Cost: add. (7941->278), mult. (17659->388), div. (0->0), fcn. (11764->8), ass. (0->180)
t148 = sin(pkin(8));
t149 = cos(pkin(8));
t151 = sin(qJ(3));
t154 = cos(qJ(3));
t126 = (-t148 * t154 - t149 * t151) * qJD(1);
t183 = qJD(1) * t154;
t128 = -qJD(1) * t148 * t151 + t149 * t183;
t196 = t128 * t126;
t213 = qJDD(3) + t196;
t217 = t148 * t213;
t216 = t149 * t213;
t150 = sin(qJ(5));
t143 = qJDD(3) + qJDD(5);
t153 = cos(qJ(5));
t95 = -t126 * t153 + t128 * t150;
t97 = t126 * t150 + t128 * t153;
t71 = t97 * t95;
t211 = -t71 + t143;
t215 = t150 * t211;
t214 = t153 * t211;
t176 = qJD(1) * qJD(3);
t171 = t154 * t176;
t174 = t151 * qJDD(1);
t132 = -t171 - t174;
t139 = t154 * qJDD(1);
t172 = t151 * t176;
t133 = t139 - t172;
t108 = t132 * t148 + t133 * t149;
t182 = qJD(3) * t126;
t87 = t108 - t182;
t107 = t132 * t149 - t133 * t148;
t61 = -qJD(5) * t95 + t107 * t150 + t108 * t153;
t144 = qJD(3) + qJD(5);
t91 = t144 * t95;
t212 = -t91 + t61;
t157 = qJD(1) ^ 2;
t186 = t154 * t157;
t152 = sin(qJ(1));
t155 = cos(qJ(1));
t170 = t152 * g(1) - g(2) * t155;
t165 = qJDD(2) - t170;
t185 = t157 * qJ(2);
t161 = t165 - t185;
t209 = pkin(6) + pkin(1);
t115 = -qJDD(1) * t209 + t161;
t188 = t154 * t115;
t80 = qJDD(3) * pkin(3) - t133 * qJ(4) + t188 + (-pkin(3) * t186 - qJ(4) * t176 + g(3)) * t151;
t110 = t154 * g(3) - t115 * t151;
t162 = qJD(3) * pkin(3) - qJ(4) * t183;
t146 = t151 ^ 2;
t193 = t146 * t157;
t81 = -pkin(3) * t193 + t132 * qJ(4) - qJD(3) * t162 - t110;
t167 = -0.2e1 * qJD(4) * t128 - t148 * t81 + t149 * t80;
t50 = 0.2e1 * qJD(4) * t126 + t148 * t80 + t149 * t81;
t210 = -pkin(7) * t87 + t167;
t93 = t95 ^ 2;
t94 = t97 ^ 2;
t124 = t126 ^ 2;
t125 = t128 ^ 2;
t142 = t144 ^ 2;
t158 = pkin(4) * t213 + t210;
t164 = qJD(3) * pkin(4) - pkin(7) * t128;
t34 = -t124 * pkin(4) + t107 * pkin(7) - qJD(3) * t164 + t50;
t17 = t150 * t34 - t153 * t158;
t201 = t153 * t34;
t18 = t150 * t158 + t201;
t8 = t150 * t18 - t153 * t17;
t208 = t148 * t8;
t207 = t149 * t8;
t145 = qJDD(1) * qJ(2);
t166 = t155 * g(1) + t152 * g(2);
t163 = -t145 + t166;
t175 = qJD(2) * qJD(1);
t159 = t163 - 0.2e1 * t175;
t82 = t132 * pkin(3) - qJDD(4) - t162 * t183 + (qJ(4) * t146 + t209) * t157 + t159;
t205 = t148 * t82;
t204 = t149 * t82;
t51 = t107 * pkin(4) + t124 * pkin(7) - t128 * t164 + t82;
t203 = t150 * t51;
t68 = t71 + t143;
t202 = t150 * t68;
t200 = t153 * t51;
t199 = t153 * t68;
t27 = t148 * t50 + t149 * t167;
t198 = t154 * t27;
t197 = qJDD(1) * pkin(1);
t195 = t144 * t150;
t194 = t144 * t153;
t147 = t154 ^ 2;
t192 = t147 * t157;
t101 = qJDD(3) - t196;
t191 = t148 * t101;
t190 = t149 * t101;
t173 = t151 * t186;
t189 = t151 * (qJDD(3) + t173);
t187 = t154 * (qJDD(3) - t173);
t184 = t146 + t147;
t181 = qJD(3) * t128;
t180 = qJD(3) * t148;
t179 = qJD(3) * t149;
t28 = -t148 * t167 + t149 * t50;
t9 = t150 * t17 + t153 * t18;
t168 = -t107 * t153 + t150 * t108;
t109 = t151 * g(3) + t188;
t76 = t154 * t109 - t151 * t110;
t85 = t107 + t181;
t160 = (-qJD(5) + t144) * t97 - t168;
t156 = qJD(3) ^ 2;
t140 = 0.2e1 * t175;
t135 = t184 * qJDD(1);
t134 = t139 - 0.2e1 * t172;
t131 = 0.2e1 * t171 + t174;
t123 = -t161 + t197;
t118 = -t125 - t156;
t117 = -t125 + t156;
t116 = t124 - t156;
t114 = t157 * t209 + t159;
t112 = -t189 + t154 * (-t156 - t192);
t111 = t151 * (-t156 - t193) + t187;
t99 = -t156 - t124;
t90 = -t94 + t142;
t89 = t93 - t142;
t88 = -t94 - t142;
t86 = t108 + t182;
t84 = -t107 + t181;
t83 = -t124 - t125;
t79 = -t118 * t148 - t190;
t78 = t118 * t149 - t191;
t73 = t149 * t99 - t217;
t72 = t148 * t99 + t216;
t70 = t94 - t93;
t66 = -t142 - t93;
t65 = (t150 * t97 - t153 * t95) * t144;
t64 = (-t150 * t95 - t153 * t97) * t144;
t63 = t148 * t87 + t149 * t85;
t62 = t148 * t85 - t149 * t87;
t60 = -qJD(5) * t97 - t168;
t59 = -t93 - t94;
t58 = t151 * t79 + t154 * t78;
t57 = t153 * t89 - t202;
t56 = -t150 * t90 + t214;
t55 = t150 * t89 + t199;
t54 = t153 * t90 + t215;
t53 = -t150 * t88 - t199;
t52 = t153 * t88 - t202;
t47 = t91 + t61;
t42 = (qJD(5) + t144) * t97 + t168;
t41 = t153 * t61 - t195 * t97;
t40 = t150 * t61 + t194 * t97;
t39 = -t150 * t60 + t194 * t95;
t38 = t153 * t60 + t195 * t95;
t37 = t151 * t73 + t154 * t72;
t36 = t153 * t66 - t215;
t35 = t150 * t66 + t214;
t33 = t151 * t63 + t154 * t62;
t31 = -t148 * t52 + t149 * t53;
t30 = t148 * t53 + t149 * t52;
t29 = -pkin(7) * t52 - t200;
t26 = t150 * t47 + t153 * t160;
t25 = -t150 * t212 - t153 * t42;
t24 = t150 * t160 - t153 * t47;
t23 = -t150 * t42 + t153 * t212;
t22 = -pkin(7) * t35 - t203;
t21 = -t148 * t35 + t149 * t36;
t20 = t148 * t36 + t149 * t35;
t19 = -pkin(4) * t212 + pkin(7) * t53 - t203;
t15 = -pkin(4) * t42 + pkin(7) * t36 + t200;
t14 = t151 * t31 + t154 * t30;
t13 = t151 * t28 + t198;
t12 = -t148 * t24 + t149 * t26;
t11 = t148 * t26 + t149 * t24;
t10 = t151 * t21 + t154 * t20;
t7 = pkin(4) * t51 + pkin(7) * t9;
t6 = -pkin(7) * t24 - t8;
t5 = t11 * t154 + t12 * t151;
t4 = -pkin(4) * t59 + pkin(7) * t26 + t9;
t3 = t149 * t9 - t208;
t2 = t148 * t9 + t207;
t1 = t151 * t3 + t154 * t2;
t16 = [0, 0, 0, 0, 0, qJDD(1), t170, t166, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t165 - 0.2e1 * t197, t140 + 0.2e1 * t145 - t166, pkin(1) * t123 + qJ(2) * (-pkin(1) * t157 + t140 - t163), (t133 - t172) * t154, -t131 * t154 - t134 * t151, t187 - t151 * (t156 - t192), (-t132 + t171) * t151, t154 * (-t156 + t193) - t189, 0, qJ(2) * t131 - t111 * t209 - t151 * t114, qJ(2) * t134 - t112 * t209 - t154 * t114, t135 * t209 - t184 * t185 - t76, -qJ(2) * t114 - t209 * t76, t154 * (t108 * t149 - t128 * t180) - t151 * (t108 * t148 + t128 * t179), t154 * (-t148 * t86 - t149 * t84) - t151 * (-t148 * t84 + t149 * t86), t154 * (-t117 * t148 + t216) - t151 * (t117 * t149 + t217), t154 * (-t107 * t148 - t126 * t179) - t151 * (t107 * t149 - t126 * t180), t154 * (t116 * t149 - t191) - t151 * (t116 * t148 + t190), (t154 * (t126 * t149 + t128 * t148) - t151 * (t126 * t148 - t128 * t149)) * qJD(3), t154 * (-qJ(4) * t72 - t205) - t151 * (-pkin(3) * t84 + qJ(4) * t73 + t204) + qJ(2) * t84 - t209 * t37, t154 * (-qJ(4) * t78 - t204) - t151 * (-pkin(3) * t86 + qJ(4) * t79 - t205) + qJ(2) * t86 - t209 * t58, t154 * (-qJ(4) * t62 - t27) - t151 * (-pkin(3) * t83 + qJ(4) * t63 + t28) + qJ(2) * t83 - t209 * t33, -qJ(4) * t198 - t151 * (pkin(3) * t82 + qJ(4) * t28) - qJ(2) * t82 - t209 * t13, t154 * (-t148 * t40 + t149 * t41) - t151 * (t148 * t41 + t149 * t40), t154 * (-t148 * t23 + t149 * t25) - t151 * (t148 * t25 + t149 * t23), t154 * (-t148 * t54 + t149 * t56) - t151 * (t148 * t56 + t149 * t54), t154 * (-t148 * t38 + t149 * t39) - t151 * (t148 * t39 + t149 * t38), t154 * (-t148 * t55 + t149 * t57) - t151 * (t148 * t57 + t149 * t55), t154 * (-t148 * t64 + t149 * t65) - t151 * (t148 * t65 + t149 * t64), t154 * (-qJ(4) * t20 - t148 * t15 + t149 * t22) - t151 * (-pkin(3) * t42 + qJ(4) * t21 + t148 * t22 + t149 * t15) + qJ(2) * t42 - t209 * t10, t154 * (-qJ(4) * t30 - t148 * t19 + t149 * t29) - t151 * (-pkin(3) * t212 + qJ(4) * t31 + t148 * t29 + t149 * t19) + qJ(2) * t212 - t209 * t14, t154 * (-qJ(4) * t11 - t148 * t4 + t149 * t6) - t151 * (-pkin(3) * t59 + qJ(4) * t12 + t148 * t6 + t149 * t4) + qJ(2) * t59 - t209 * t5, t154 * (-pkin(7) * t207 - qJ(4) * t2 - t148 * t7) - t151 * (pkin(3) * t51 - pkin(7) * t208 + qJ(4) * t3 + t149 * t7) - qJ(2) * t51 - t209 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t157, -t123, 0, 0, 0, 0, 0, 0, t111, t112, -t135, t76, 0, 0, 0, 0, 0, 0, t37, t58, t33, t13, 0, 0, 0, 0, 0, 0, t10, t14, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, (-t146 + t147) * t157, t139, -t173, -t174, qJDD(3), t109, t110, 0, 0, -t196, t125 - t124, t87, t196, t85, qJDD(3), pkin(3) * t72 + t167, pkin(3) * t78 - t50, pkin(3) * t62, pkin(3) * t27, t71, t70, t47, -t71, t160, t143, pkin(3) * t20 + pkin(4) * t35 - t17, -t201 - t150 * t210 + pkin(3) * t30 + (-t150 * t213 + t52) * pkin(4), pkin(3) * t11 + pkin(4) * t24, pkin(3) * t2 + pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t86, t83, -t82, 0, 0, 0, 0, 0, 0, t42, t212, t59, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, t47, -t71, t160, t143, -t17, -t18, 0, 0;];
tauJ_reg = t16;
