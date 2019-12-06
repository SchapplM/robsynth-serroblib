% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:24
% EndTime: 2019-12-05 16:23:32
% DurationCPUTime: 2.21s
% Computational Cost: add. (7268->283), mult. (16374->406), div. (0->0), fcn. (11528->10), ass. (0->175)
t215 = -2 * qJD(4);
t154 = sin(pkin(9));
t156 = cos(pkin(9));
t163 = cos(qJ(3));
t160 = sin(qJ(3));
t184 = qJD(2) * t160;
t126 = -t156 * t163 * qJD(2) + t154 * t184;
t128 = (t163 * t154 + t160 * t156) * qJD(2);
t107 = t128 * t126;
t206 = qJDD(3) - t107;
t214 = t154 * t206;
t213 = t156 * t206;
t159 = sin(qJ(5));
t150 = qJDD(3) + qJDD(5);
t162 = cos(qJ(5));
t100 = -t126 * t159 + t128 * t162;
t98 = t162 * t126 + t128 * t159;
t72 = t100 * t98;
t207 = -t72 + t150;
t212 = t159 * t207;
t211 = t162 * t207;
t178 = qJD(2) * qJD(3);
t173 = t163 * t178;
t177 = t160 * qJDD(2);
t133 = t173 + t177;
t204 = qJD(2) ^ 2;
t179 = t163 * t204;
t143 = t160 * t179;
t140 = qJDD(3) + t143;
t155 = sin(pkin(8));
t157 = cos(pkin(8));
t138 = -g(1) * t157 - g(2) * t155;
t161 = sin(qJ(2));
t164 = cos(qJ(2));
t185 = -g(3) + qJDD(1);
t116 = t164 * t138 + t161 * t185;
t111 = -t204 * pkin(2) + qJDD(2) * pkin(6) + t116;
t137 = -g(1) * t155 + g(2) * t157;
t94 = t111 * t160 - t163 * t137;
t78 = (-t133 + t173) * qJ(4) + t140 * pkin(3) - t94;
t175 = qJ(4) * t184;
t139 = qJD(3) * pkin(3) - t175;
t188 = t160 * t137;
t79 = t188 + (-t139 - t175) * qJD(3) + (-pkin(3) * t179 + qJDD(2) * qJ(4) + t111) * t163;
t169 = t128 * t215 - t154 * t79 + t156 * t78;
t174 = t160 * t178;
t176 = t163 * qJDD(2);
t168 = -t174 + t176;
t109 = t156 * t133 + t154 * t168;
t180 = t126 * qJD(3);
t170 = -t109 - t180;
t210 = pkin(7) * t170 + t169;
t108 = -t133 * t154 + t156 * t168;
t61 = -qJD(5) * t98 + t108 * t159 + t109 * t162;
t151 = qJD(3) + qJD(5);
t93 = t151 * t98;
t208 = t61 - t93;
t37 = t126 * t215 + t154 * t78 + t156 * t79;
t96 = t98 ^ 2;
t97 = t100 ^ 2;
t124 = t126 ^ 2;
t125 = t128 ^ 2;
t149 = t151 ^ 2;
t205 = t163 ^ 2;
t166 = pkin(4) * t206 + t210;
t114 = qJD(3) * pkin(4) - pkin(7) * t128;
t33 = -pkin(4) * t124 + pkin(7) * t108 - qJD(3) * t114 + t37;
t16 = t159 * t33 - t162 * t166;
t195 = t162 * t33;
t17 = t159 * t166 + t195;
t8 = t159 * t17 - t16 * t162;
t203 = t154 * t8;
t202 = t156 * t8;
t115 = -t138 * t161 + t164 * t185;
t110 = -qJDD(2) * pkin(2) - t204 * pkin(6) - t115;
t148 = t205 * t204;
t82 = -t168 * pkin(3) - qJ(4) * t148 + t139 * t184 + qJDD(4) + t110;
t200 = t154 * t82;
t199 = t156 * t82;
t51 = -t108 * pkin(4) - t124 * pkin(7) + t114 * t128 + t82;
t198 = t159 * t51;
t69 = t72 + t150;
t197 = t159 * t69;
t20 = t154 * t37 + t156 * t169;
t196 = t160 * t20;
t194 = t162 * t51;
t193 = t162 * t69;
t103 = qJDD(3) + t107;
t192 = t103 * t154;
t191 = t103 * t156;
t190 = t151 * t159;
t189 = t151 * t162;
t187 = t160 * t140;
t186 = t163 * (qJDD(3) - t143);
t183 = qJD(3) * t128;
t9 = t159 * t16 + t162 * t17;
t21 = -t154 * t169 + t156 * t37;
t95 = t163 * t111 + t188;
t64 = t160 * t94 + t163 * t95;
t171 = -t162 * t108 + t109 * t159;
t85 = t108 + t183;
t134 = -0.2e1 * t174 + t176;
t167 = (-qJD(5) + t151) * t100 - t171;
t165 = qJD(3) ^ 2;
t152 = t160 ^ 2;
t147 = t152 * t204;
t136 = t147 + t148;
t135 = (t152 + t205) * qJDD(2);
t132 = 0.2e1 * t173 + t177;
t119 = -t125 - t165;
t118 = -t125 + t165;
t117 = t124 - t165;
t113 = -t186 - t160 * (-t147 - t165);
t112 = t163 * (-t148 - t165) - t187;
t101 = -t165 - t124;
t91 = -t97 + t149;
t90 = t96 - t149;
t88 = -t97 - t149;
t86 = t109 - t180;
t84 = -t108 + t183;
t83 = -t124 - t125;
t81 = -t119 * t154 - t191;
t80 = t119 * t156 - t192;
t76 = t101 * t156 - t214;
t75 = t101 * t154 + t213;
t71 = t97 - t96;
t67 = -t149 - t96;
t66 = (t100 * t159 - t162 * t98) * t151;
t65 = (-t100 * t162 - t159 * t98) * t151;
t63 = -t154 * t170 + t156 * t85;
t62 = t154 * t85 + t156 * t170;
t60 = -qJD(5) * t100 - t171;
t59 = -t96 - t97;
t58 = -t160 * t80 + t163 * t81;
t57 = t162 * t90 - t197;
t56 = -t159 * t91 + t211;
t55 = t159 * t90 + t193;
t54 = t162 * t91 + t212;
t53 = -t159 * t88 - t193;
t52 = t162 * t88 - t197;
t50 = t61 + t93;
t45 = (qJD(5) + t151) * t100 + t171;
t44 = -t100 * t190 + t162 * t61;
t43 = t100 * t189 + t159 * t61;
t42 = -t159 * t60 + t98 * t189;
t41 = t162 * t60 + t98 * t190;
t40 = -t160 * t75 + t163 * t76;
t39 = t162 * t67 - t212;
t38 = t159 * t67 + t211;
t34 = -t160 * t62 + t163 * t63;
t31 = -t154 * t52 + t156 * t53;
t30 = t154 * t53 + t156 * t52;
t29 = -pkin(7) * t52 + t194;
t28 = t159 * t50 + t162 * t167;
t27 = -t159 * t208 - t162 * t45;
t26 = t159 * t167 - t162 * t50;
t25 = -t159 * t45 + t162 * t208;
t24 = -pkin(7) * t38 + t198;
t23 = -t154 * t38 + t156 * t39;
t22 = t154 * t39 + t156 * t38;
t19 = -pkin(4) * t208 + pkin(7) * t53 + t198;
t18 = -pkin(4) * t45 + pkin(7) * t39 - t194;
t14 = -t160 * t30 + t163 * t31;
t13 = -t154 * t26 + t156 * t28;
t12 = t154 * t28 + t156 * t26;
t11 = -t160 * t22 + t163 * t23;
t10 = t163 * t21 - t196;
t7 = -pkin(4) * t51 + pkin(7) * t9;
t6 = -t160 * t12 + t163 * t13;
t5 = -pkin(7) * t26 - t8;
t4 = -pkin(4) * t59 + pkin(7) * t28 + t9;
t3 = t156 * t9 - t203;
t2 = t154 * t9 + t202;
t1 = -t160 * t2 + t163 * t3;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t185, 0, 0, 0, 0, 0, 0, qJDD(2) * t164 - t204 * t161, -qJDD(2) * t161 - t204 * t164, 0, t115 * t164 + t116 * t161, 0, 0, 0, 0, 0, 0, t112 * t161 + t134 * t164, t113 * t161 - t132 * t164, t135 * t161 + t136 * t164, -t110 * t164 + t161 * t64, 0, 0, 0, 0, 0, 0, t161 * t40 - t164 * t84, t161 * t58 - t164 * t86, t161 * t34 - t164 * t83, t10 * t161 - t164 * t82, 0, 0, 0, 0, 0, 0, t11 * t161 - t164 * t45, t14 * t161 - t164 * t208, t161 * t6 - t164 * t59, t1 * t161 - t164 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t115, -t116, 0, 0, (t133 + t173) * t160, t132 * t163 + t134 * t160, t187 + t163 * (-t147 + t165), t134 * t163, t160 * (t148 - t165) + t186, 0, pkin(2) * t134 + pkin(6) * t112 - t110 * t163, -pkin(2) * t132 + pkin(6) * t113 + t110 * t160, pkin(2) * t136 + pkin(6) * t135 + t64, -pkin(2) * t110 + pkin(6) * t64, t160 * (t109 * t156 - t154 * t183) + t163 * (t109 * t154 + t156 * t183), t160 * (-t154 * t86 - t156 * t84) + t163 * (-t154 * t84 + t156 * t86), t160 * (-t118 * t154 + t213) + t163 * (t118 * t156 + t214), t160 * (-t108 * t154 + t156 * t180) + t163 * (t108 * t156 + t154 * t180), t160 * (t117 * t156 - t192) + t163 * (t117 * t154 + t191), (t160 * (-t126 * t156 + t128 * t154) + t163 * (-t126 * t154 - t128 * t156)) * qJD(3), t160 * (-qJ(4) * t75 + t200) + t163 * (-pkin(3) * t84 + qJ(4) * t76 - t199) - pkin(2) * t84 + pkin(6) * t40, t160 * (-qJ(4) * t80 + t199) + t163 * (-pkin(3) * t86 + qJ(4) * t81 + t200) - pkin(2) * t86 + pkin(6) * t58, t160 * (-qJ(4) * t62 - t20) + t163 * (-pkin(3) * t83 + qJ(4) * t63 + t21) - pkin(2) * t83 + pkin(6) * t34, -qJ(4) * t196 + t163 * (-pkin(3) * t82 + qJ(4) * t21) - pkin(2) * t82 + pkin(6) * t10, t160 * (-t154 * t43 + t156 * t44) + t163 * (t154 * t44 + t156 * t43), t160 * (-t154 * t25 + t156 * t27) + t163 * (t154 * t27 + t156 * t25), t160 * (-t154 * t54 + t156 * t56) + t163 * (t154 * t56 + t156 * t54), t160 * (-t154 * t41 + t156 * t42) + t163 * (t154 * t42 + t156 * t41), t160 * (-t154 * t55 + t156 * t57) + t163 * (t154 * t57 + t156 * t55), t160 * (-t154 * t65 + t156 * t66) + t163 * (t154 * t66 + t156 * t65), t160 * (-qJ(4) * t22 - t154 * t18 + t156 * t24) + t163 * (-pkin(3) * t45 + qJ(4) * t23 + t154 * t24 + t156 * t18) - pkin(2) * t45 + pkin(6) * t11, t160 * (-qJ(4) * t30 - t154 * t19 + t156 * t29) + t163 * (-pkin(3) * t208 + qJ(4) * t31 + t154 * t29 + t156 * t19) - pkin(2) * t208 + pkin(6) * t14, t160 * (-qJ(4) * t12 - t154 * t4 + t156 * t5) + t163 * (-pkin(3) * t59 + qJ(4) * t13 + t154 * t5 + t156 * t4) - pkin(2) * t59 + pkin(6) * t6, t160 * (-pkin(7) * t202 - qJ(4) * t2 - t154 * t7) + t163 * (-pkin(3) * t51 - pkin(7) * t203 + qJ(4) * t3 + t156 * t7) - pkin(2) * t51 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t147 - t148, t177, t143, t176, qJDD(3), -t94, -t95, 0, 0, t107, t125 - t124, -t170, -t107, t85, qJDD(3), pkin(3) * t75 + t169, pkin(3) * t80 - t37, pkin(3) * t62, pkin(3) * t20, t72, t71, t50, -t72, t167, t150, pkin(3) * t22 + pkin(4) * t38 - t16, -t195 - t159 * t210 + pkin(3) * t30 + (-t159 * t206 + t52) * pkin(4), pkin(3) * t12 + pkin(4) * t26, pkin(3) * t2 + pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t86, t83, t82, 0, 0, 0, 0, 0, 0, t45, t208, t59, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, t50, -t72, t167, t150, -t16, -t17, 0, 0;];
tauJ_reg = t15;
